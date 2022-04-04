"""
Using an existing database, aggreagate hashes from all strains 
into a single signature per species. Build zipfile database 
from these species-wide signatures. Perhaps useful for
faster search, taxonomic classification and/or decontamination.
Author: N. Tessa Pierce-Ward
"""
import pandas as pd
import sourmash

taxfile= "gtdb-rs202.taxonomy.v2.with-strain.csv"
taxDF = pd.read_csv(taxfile)
all_species = taxDF["species"].unique().tolist()

# species filenames need to have no spaces. Sig names keep spaces
speciesD = {s.replace(' ', '_'):s for s in all_species}

# sig grep --> sig merge

out_dir = "outputs.species-db"
logs_dir = os.path.join(out_dir, "logs") 
alphabet="nucleotide"
ksize=31

rule all:
    input: os.path.join(out_dir, f"gtdb-rs202.species.{alphabet}-k{ksize}.zip")


# could make the intermediate files temporary to avoid storing once we have the full zipfile
# IF do this, need sigs to be `input` for `cat_species_db` rule, otherwise they'll be deleted
#
rule make_species_sigs:
    input:
        database="gtdb-rs202.{alphabet}-k{ksize}.zip",
    output:
        sig = os.path.join(out_dir, "sigs", "{species}.{alphabet}-k{ksize}.sig.gz")
    params:
        signame = lambda w: speciesD[w.species], # get species name WITH the space
        grep = lambda w: speciesD[w.species].split("s__")[1], # i don't have "s__" in my sig names
    conda:
        "conf/env/sourmash4.yml"
    log: os.path.join(logs_dir, "make-species-sigs", "{species}.{alphabet}-k{ksize}.log")
    benchmark: os.path.join(logs_dir, "make-species-sigs", "{species}.{alphabet}-k{ksize}.benchmark")
    shell:
        """
        sourmash sig grep {params.grep:q} {input.database} | \
        sourmash sig merge --name {params.signame:q} \
        -o {output.sig} 2> {log}
        """


rule sigs_to_file:
    input: expand(os.path.join(out_dir, "sigs", "{species}.{{alphabet}}-k{{ksize}}.sig.gz"), species=speciesD.keys()), # use species names WITHOUT spaces
    output: os.path.join(out_dir, "gtdb-rs202.species.{alphabet}-k{ksize}.siglist"),
    log: os.path.join(logs_dir, "sigs_to_file", "gtdb-rs202.species.{alphabet}-k{ksize}.log")
    benchmark: os.path.join(logs_dir, "sigs_to_file", "gtdb-rs202.species.{alphabet}-k{ksize}.benchmark")
    run:
        with open(str(output)) as out:
            for inF in input:
                fn = os.path.abspath(str(inF))
                out.write(f"{out}\n")
    

rule cat_species_db:
    input:
        siglist=os.path.join(out_dir, "gtdb-rs202.species.{alphabet}-k{ksize}.siglist"),
        sigs=expand(os.path.join(out_dir, "sigs", "{species}.{{alphabet}}-k{{ksize}}.sig.gz"), species=speciesD.keys()),
    output: os.path.join(out_dir, "gtdb-rs202.species.{alphabet}-k{ksize}.zip")
    conda:
        "conf/env/sourmash4.yml"
    params:
        sigdir=os.path.join(out_dir, "sigs")
    log: os.path.join(logs_dir, "cat-species-db", "gtdb-rs202.species.{alphabet}-k{ksize}.log")
    benchmark: os.path.join(logs_dir, "cat-species-db", "gtdb-rs202.species.{alphabet}-k{ksize}.benchmark")
    shell:
        """
        sourmash sig cat --from-file {input.siglist} -o {output} 2> {log}
        """
    

