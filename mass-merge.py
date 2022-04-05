#! /usr/bin/env python
"""
Merge signatures in bulk, based on an attribute in the spreadsheet.

mass-merge.py takes a list of databases along with spreadsheets (w/-F),
and merges signatures by values in the specified `merge-col`.

The spreadsheet must contain two columns, 'ident' and the value of `--merge-col`;
signatures are selected based on 'ident' and renamed to the value found in the merge column.
Singletons (no additional signatures to merge) will be renamed for consistency.
"""
import sys
import argparse
import csv
from collections import defaultdict

import sourmash

from sourmash.picklist import SignaturePicklist
from sourmash.logging import set_quiet, error, notify, print_results, debug
from sourmash import sourmash_args
from sourmash.cli.utils import (add_moltype_args, add_ksize_arg)
#from sourmash.sig import _check_abundance_compatibility

def _check_abundance_compatibility(sig1, sig2):
    if sig1.minhash.track_abundance != sig2.minhash.track_abundance:
        raise ValueError("incompatible signatures: track_abundance is {} in first sig, {} in second".format(sig1.minhash.track_abundance, sig2.minhash.track_abundance))


def massmerge(args):
    """
    rename one or more signatures.
    """
    set_quiet(args.quiet, args.quiet)
    moltype = sourmash_args.calculate_moltype(args)
    merge_col = args.merge_col
    #CTB _extend_signatures_with_from_file(args)

    # load spreadsheets
    merge_d = defaultdict(list)
    all_idents=set()
    for filename in args.from_spreadsheet:
        count = 0
        with open(filename, newline='') as fp:
            first_entry= True
            r = csv.DictReader(fp)
            for row in r:
                if first_entry:
                    if not merge_col in r.fieldnames:
                        error(f"ERROR on spreadsheet '{filename}'.")
                        error(f"Merge column {merge_col} is not present.")
                        sys.exit(-1)            
                merge_name = row[merge_col]
                ident = row['ident']
                
                assert ' ' not in ident, f"identifiers cannot have spaces - but '{ident}' does."
                assert ident not in all_idents, f"duplicate identifer: '{ident}'"
                all_idents.add(ident)

                merge_d[merge_name].append(ident)
                count += 1
        notify(f"loaded {count} identifiers from '{filename}'")

    num_merge_names = len(merge_d)
    notify(f"found a total of {num_merge_names} distinct values for signature merging.")

    # make new dict with the idx refs for all idents/alldbs
    merge_idx_d = defaultdict(list)
    found_idents = set()
    # go through all the database and load etc.
    idx_list = []
    for db in args.dblist:
        notify(f"loading index '{db}'")
        idx = sourmash.load_file_as_index(db)

        manifest = idx.manifest
        if manifest is None:
            error(f"ERROR on filename '{db}'.")
            error("No manifest, but a manifest is required.")
            sys.exit(-1)

        if args.check:
            # just check we have all idents
            ident_picklist = SignaturePicklist('ident')
            ident_picklist.pickset = all_idents
            idx = idx.select(ksize=args.ksize,
                                moltype=moltype,
                                picklist=ident_picklist)
            found_idents = ident_picklist.found
        else:
            # actually find sigs and keep track of locations across all dbs/manifests
            for n, (merge_name, idents) in enumerate(merge_d.items()):
                if n % 100 == 0:
                    merge_percent = float(n)/num_merge_names * 100
                    notify(f"...finding sigs for merge name {merge_name}; {merge_percent:.1f}% searched", end="\r")
                # build a new picklist for just these idents
                ident_picklist = SignaturePicklist('ident')
                ident_picklist.pickset = set(idents)

                this_idx = idx.select(ksize=args.ksize,
                                moltype=moltype,
                                picklist=ident_picklist)

                idx_list.append(this_idx)

                # store idx and found idents
                found_idents.update(ident_picklist.found)
                merge_idx_d[merge_name] += idx_list

    # make sure that we get all the things.
    if not all_idents.issubset(found_idents):
        remaining = all_idents - found_idents
        error(f"ERROR: {len(remaining)} identifiers from spreadsheet not found.")
        example_missing = "\n".join(remaining)
        error(f"Here are some examples: {example_missing}")
        sys.exit(-1)

    if args.check:
        notify("Everything looks copacetic. Exiting as requested by `--check`")
        sys.exit(0)

    notify("Everything looks copacetic. Proceeding to merge!")
        
    # go through, do merge, save.
    with sourmash_args.SaveSignaturesToLocation(args.output) as save_sigs:
        n = 0
        num_singletons = 0
        for m, (merge_name, idx_list) in enumerate(merge_idx_d.items()):
            if m % 100 == 0:
                merge_percent = float(n)/found_idents * 100
                notify(f"...at merge name {m}; {merge_percent:.1f}% processed", end="\r")

            # if only one item, just rename and save
            if len(idx_list) == 1:
                ss = idx_list[0].signatures()
                ss._name = merge_name
                if args.flatten:
                    ss.minhash = ss.minhash.flatten()
                save_sigs.add(ss)
                num_singletons+=1
                n += 1
            else: # merge sigs
                first_sig = None
                mh = None
                for idx in idx_list:
                    for ss in idx.signatures():
                        n += 1
                        # first sig? initialize some things
                        if first_sig is None:
                            first_sig = ss
                            mh = first_sig.minhash.copy_and_clear()

                            # forcibly remove abundance?
                            if args.flatten:
                                mh.track_abundance = False

                        try:
                            sigobj_mh = ss.minhash
                            if not args.flatten:
                                _check_abundance_compatibility(first_sig, ss)
                            else:
                                sigobj_mh.track_abundance = False

                            mh.merge(sigobj_mh)
                        except (TypeError, ValueError) as exc:
                            error("ERROR when merging signature '{}' ({}) from file {}",
                                ss, ss.md5sum()[:8])
                            error(str(exc))
                            sys.exit(-1)

                merged_ss = sourmash.SourmashSignature(mh, name=merge_name)
                save_sigs.add(merged_ss)

    notify(f"merged {n} signatures into {len(save_sigs)} signatures by column {merge_col}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('dblist', nargs='+')

    p.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    p.add_argument(
        '--flatten', action='store_true',
        help='remove abundances from all signatures while merging'
    )
    p.add_argument(
        '--merge-col', required=True,
        help='the column to merge signatures by (required)'
    )
    p.add_argument(
        '-o', '--output', metavar='FILE', default='-',
        help='output signature to this file (default stdout)'
    )
    p.add_argument(
        '-f', '--force', action='store_true',
        help='try to load all files as signatures'
    )
    p.add_argument(
        '--check', action='store_true',
        help='Just check for ability to merge; do not actually merge signatures.'
    )
    p.add_argument('-F', '--from-spreadsheet',
                   required=True,
                   action='append', default=[],
                   help="input spreadsheet containing 'ident' and '--merge-col` columns")

    add_ksize_arg(p, 31)
    add_moltype_args(p)

    args = p.parse_args()

    massmerge(args)


if __name__ == '__main__':
    sys.exit(main())
