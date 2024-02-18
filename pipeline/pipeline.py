import os
import sys
import glob
import time

from tqdm import tqdm

import Bio.SearchIO.HmmerIO.hmmer3_text as hmmio
import Bio.SeqIO as seqio
from Bio.Blast import NCBIXML

sys.path.append(os.path.abspath("./CLANS-Python/"))

import clans.config as cfg
import clans.io.parser as parser
import clans.io.file_handler as fh
import clans.similarity_search.blast as blast
import clans.layouts.layout_handler as lh


def pipeline(
    Emin_hmm=1e-10,
    protein_tmdom_th=4,
):
    # BLAST Database for Cel
    cel_name = os.path.abspath(r"./cel/caenorhabditis_elegans.PRJNA13758.WBPS14.protein.fa")
    db_cmd = r"makeblastdb -dbtype prot -out cel.db -in {}".format(cel_name)
    os.system(db_cmd)

    # Build Rhodopsins
    gpcr_name = r"rhodopsins"
    base_name = os.path.abspath(r"./cel/rhodopsins_cel_alignment.fas")
    build_cmd = r"hmmbuild {}.hmm {}".format(gpcr_name, base_name)
    os.system(build_cmd)

    # Build Secretins
    gpcr_name = r"secretins"
    base_name = os.path.abspath(r"./cel/class_b_secretins.fas")
    build_cmd = r"hmmbuild {}.hmm {}".format(gpcr_name, base_name)
    os.system(build_cmd)

    # Build New_GPCRs
    gpcr_name = r"new_GPCRs"
    base_name = os.path.abspath(r"./cel/new_GPCRs.fas")
    build_cmd = r"hmmbuild {}.hmm {}".format(gpcr_name, base_name)
    os.system(build_cmd)

    # Search
    species = glob.glob("./nematodes/*.fa")
    species = [s.replace("./nematodes/", "") for s in species]
    species_names = []
    smax = len(species)

    print("### HMM Search ###")
    for s in tqdm(species[:smax]):
        species_names.append(s[: s.find(".")])
        for gpcr_name in ["rhodopsins", "secretins", "new_GPCRs"]:
            proteome_path = os.path.join("./nematodes", s)
            search_cmd = r"hmmsearch -E {0} {1}.hmm {2} > {3}_{1}.sea".format(
                Emin_hmm, gpcr_name, proteome_path, species_names[-1]
            )
            os.system(search_cmd)
            print("Search completed: {}".format(species_names[-1]))

            # Collect Sequences
            with open("{}_{}.sea".format(species_names[-1], gpcr_name), "r") as fp:
                parser = hmmio.Hmmer3TextParser(fp)
                entries = [f for f in parser]
                hits = entries[0].hit_keys

            fasta_contents = [
                f for f in seqio.parse(proteome_path, "fasta") if f.id in hits
            ]
            for f in fasta_contents:
                f.description = f.name
            hit_name = "{}_{}.hit".format(species_names[-1], gpcr_name)
            seqio.write(fasta_contents, hit_name, "fasta")

            # HMMTOP
            hmmtop_input = hit_name
            hmmtop_cmd = r"hmmtop -if={} -of={}_{}.tra".format(
                hmmtop_input, species_names[-1], gpcr_name
            )
            os.system(hmmtop_cmd)

            saves = []
            with open(r"{}_{}.tra".format(species_names[-1], gpcr_name), "r") as fp:
                while True:
                    prot = fp.readline()
                    if prot == "":
                        break
                    parsed_line = [p for p in prot.split(" ") if p != ""]
                    protein_name, protein_tmdom, tm_type = (
                        parsed_line[2],
                        int(parsed_line[4]),
                        parsed_line[3],
                    )
                    if (protein_tmdom >= protein_tmdom_th) and (tm_type == "OUT"):
                        saves.append(protein_name)

            # Prune fasta_contents by removing all proteins with less than T transmembrane domains
            fasta_contents = [x for x in fasta_contents if x.id in saves]
            matches_name = "{}_{}_tm{}.fa".format(
                species_names[-1], gpcr_name, protein_tmdom_th
            )
            seqio.write(fasta_contents, matches_name, "fasta")

    # Join all results
    all_fasta = []
    for spec in species_names:
        for gpcr_name in ["rhodopsins", "secretins", "new_GPCRs"]:
            c = [
                f
                for f in seqio.parse(
                    "{}_{}_tm{}.fa".format(spec, gpcr_name, protein_tmdom_th), "fasta"
                )
            ]
            all_fasta.extend(c)

    # Join Cel gpcrs
    for curated_name in [
        r"./curated/rhodopsins.fa",
        r"./curated/class_b_secretins.fa",
        r"./curated/new_GPCRs.fa",
    ]:
        cel_gpcrs = [f for f in seqio.parse(os.path.abspath(curated_name), "fasta")]
        all_fasta.extend(cel_gpcrs)

    seqio.write(all_fasta, "output_nematodes.fa", "fasta")

    # Move main output file
    output_name = "./output/output_nematodes_{}.fa".format(int(time.time_ns() / 1e6))
    os.makedirs("./output", exist_ok=True)
    os.rename("output_nematodes.fa", output_name)

    # Remove duplicates
    # Clean up directory
    to_delete = []
    to_delete.extend(glob.glob("./*.sea"))
    to_delete.extend(glob.glob("./*.hmm"))
    to_delete.extend(glob.glob("./*.hit"))
    to_delete.extend(glob.glob("./*.db.*"))
    to_delete.extend(glob.glob("./*.xml"))
    to_delete.extend(glob.glob("./*.tra"))
    to_delete.extend(glob.glob("./*_matches*"))
    to_delete.extend(glob.glob("./*_tm*.fa"))
    if len(to_delete) > 0:
        [os.remove(t) for t in to_delete]

    print("Correctly terminating pipeline.")
    return output_name


def clans(output_name, eval=1e-5, iters=20000):
    # Setup configuration for running CLANS
    print("### CLANS Interface ###")
    input_file = os.path.abspath(output_name)
    prefix = os.path.split(output_name)[1].replace(".fa", "")

    cfg.run_params["input_file"] = input_file
    cfg.run_params["input_format"] = "fasta"
    cfg.run_params["run_blast"] = True
    cfg.run_params["output_file"] = "./output/{}.clans".format(prefix)
    cfg.run_params["output_format"] = "clans"

    cfg.run_params["evalue_cutoff"] = eval
    cfg.run_params["scoring_matrix"] = cfg.BLAST_scoring_matrix
    cfg.run_params["num_of_rounds"] = iters
    cfg.run_params["similarity_cutoff"] = cfg.similarity_cutoff
    cfg.run_params.update(cfg.layouts["FR"]["params"])
    cfg.run_params["is_debug_mode"] = False
    cfg.run_params["dimensions_num_for_clustering"] = 2
    cfg.run_params["working_dir"] = "./CLANS-Python"

    # Read the input file (fasta/clans/delimited) and fill the relevant main data-structures
    before = time.time()
    fh.read_input_file(cfg.run_params["input_file"], cfg.run_params["input_format"])
    after = time.time()
    duration = after - before
    if cfg.run_params["is_problem"]:
        print(cfg.run_params["error"])
        exit()
    else:
        print("Reading the input file took " + str(duration) + " seconds")

    # Perform BLAST search and fill the HSP's E-values in the similarity matrix
    if cfg.run_params["run_blast"]:
        before = time.time()
        blast.find_HSPs()
        after = time.time()
        duration = after - before
        if cfg.run_params["is_problem"]:
            print(cfg.run_params["error"])
            exit()
        else:
            print("Performing the BLAST search took " + str(duration) + " seconds")

    # Run the Fruchterman-Reingold layout calculation for the defined number of rounds
    if cfg.run_params["num_of_rounds"] > 0:
        cfg.run_params["rounds_done"] = 0
        before = time.time()
        lh.calculate_layout("FR")
        after = time.time()
        duration = after - before
        print(
            "The calculation of "
            + str(cfg.run_params["rounds_done"])
            + " rounds took "
            + str(duration)
            + " seconds"
        )

    # Write the output file
    if cfg.run_params["output_file"] is not None:
        before = time.time()
        fh.write_file(cfg.run_params["output_file"], cfg.run_params["output_format"])
        after = time.time()
        duration = after - before
        print("Writing the output file took " + str(duration) + " seconds")
    return cfg


if __name__ == "__main__":
    out_seq = pipeline()
    clans_env = clans(out_seq)
    exit(0)
