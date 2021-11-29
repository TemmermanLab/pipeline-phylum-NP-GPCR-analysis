import os
import sys
import glob
import time
import Bio.SearchIO.HmmerIO.hmmer3_text as hmmio
import Bio.SeqIO as seqio
from Bio.Blast import NCBIXML

sys.path.append('../../CLANS-Python/')

def pipeline(Emin_hmm=1e-10,
             do_blast_vs_cel=False, Emin_blast=1e-10, Emin_blast_post=1e-10,
             protein_tmdom_th=6):
    # BLAST Database for Cel
    cel_name = r'../cel/caenorhabditis_elegans.PRJNA13758.WBPS14.protein.fa'
    db_cmd = r"makeblastdb -dbtype prot -out cel.db -in {}".format(cel_name)
    os.system(db_cmd)

    # Build
    gpcr_name = r'rhodopsins'
    base_name = r'../cel/rhodopsins_cel_alignment.fas'
    build_cmd = r"hmmbuild {}.hmm {}".format(gpcr_name, base_name)
    os.system(build_cmd)

    # Load Curated File
    curated_name = r'../curated/rhodopsins.fa'
    cel_gpcr_fa = [f.description for f in seqio.parse(curated_name, "fasta")]

    # Search
    species = glob.glob('../nematodes/*.fa')
    species = [s.replace('../nematodes/', '') for s in species]
    species_names = []
    for s in species:
        species_names.append(s[:s.find('.')])
        proteome_path = os.path.join('../nematodes', s)
        search_cmd = r"hmmsearch -E {} {}.hmm {} > {}.sea".format(Emin_hmm,
                                                                  gpcr_name,
                                                                  proteome_path,
                                                                  species_names[-1])
        os.system(search_cmd)
        print('Search completed: {}'.format(species_names[-1]))

        # Collect Sequences
        with open('{}.sea'.format(species_names[-1]), 'r') as fp:
            parser = hmmio.Hmmer3TextParser(fp)
            entries = [f for f in parser]
            hits = entries[0].hit_keys

        fasta_contents = [f for f in seqio.parse(
            proteome_path, "fasta") if f.id in hits]
        for f in fasta_contents:
            f.description = f.name
        hit_name = "{}.hit".format(species_names[-1])
        seqio.write(fasta_contents, hit_name, "fasta")

        if do_blast_vs_cel:
            # Blast vs. Cel
            blast_cmd = r"blastp -evalue {} -query {}.hit -db cel.db -out blast_result.xml -outfmt 5".format(
                Emin_blast, species_names[-1], cel_name)
            os.system(blast_cmd)
            records = [r for r in NCBIXML.parse(open("blast_result.xml"))]

            saves = []
            for r in records:
                # Each record is one query from the selected HMM hits of the species
                if r.alignments:
                    # Each alignment is one BLAST hit
                    for align in r.alignments:
                        hit_id = align.hit_def.split(' ')[0]
                        if hit_id in cel_gpcr_fa:
                            saves.append(r.query)
                            break

            fasta_contents = [x for x in fasta_contents if x.id in saves]
            matches_name = '{}_gpcr_matches.fa'.format(species_names[-1])
            seqio.write(fasta_contents, matches_name, 'fasta')

        hmmtop_input = matches_name if do_blast_vs_cel else hit_name

        # HMMTOP
        hmmtop_cmd = r'hmmtop -if={} -of={}.tra'.format(
            hmmtop_input, species_names[-1])
        os.system(hmmtop_cmd)

        saves = []
        with open(r"{}.tra".format(species_names[-1]), 'r') as fp:
            while True:
                prot = fp.readline()
                if (prot == ''):
                    break
                parsed_line = [p for p in prot.split(' ') if p != '']
                protein_name, protein_tmdom = parsed_line[2], int(
                    parsed_line[4])
                if protein_tmdom >= protein_tmdom_th:
                    saves.append(protein_name)

        # Prune fasta_contents by removing all proteins with less than 4 transmembrane domains
        fasta_contents = [x for x in fasta_contents if x.id in saves]
        matches_name = '{}_tm{}.fa'.format(species_names[-1], protein_tmdom_th)
        seqio.write(fasta_contents, matches_name, 'fasta')

    # Join all results
    all_fasta = []
    for s in species_names:
        c = [f for f in seqio.parse(
            '{}_tm{}.fa'.format(s, protein_tmdom_th), "fasta")]
        all_fasta.extend(c)

    # Join Cel gpcrs
    cel_gpcrs = [f for f in seqio.parse(curated_name, "fasta")]
    all_fasta.extend(cel_gpcrs)
    seqio.write(all_fasta, 'output_nematodes.fa', 'fasta')

    # Move main output file
    output_name = '../output/output_nematodes_{}.fa'.format(time.time_ns())
    os.makedirs('../output', exist_ok=True)
    os.rename('output_nematodes.fa',
              output_name)

    # Clean up directory
    to_delete = []
    to_delete.extend(glob.glob('./*.sea'))
    to_delete.extend(glob.glob('./*.hmm'))
    to_delete.extend(glob.glob('./*.hit'))
    to_delete.extend(glob.glob('./*.db.*'))
    to_delete.extend(glob.glob('./*.xml'))
    to_delete.extend(glob.glob('./*.tra'))
    to_delete.extend(glob.glob('./*_matches*'))
    to_delete.extend(glob.glob('./*_tm*.fa'))
    if len(to_delete) > 0:
        [os.remove(t) for t in to_delete]

    print('Correctly terminating pipeline.')
    return output_name


def clans(output_name, eval=1e-5):
    input_file = os.path.abspath(output_name)
    prefix = os.path.splitext(output_name)[0]
    clans_cmd = r'python clans_cmd.py -infile {} -eval {} -saveto out_{}.clans'.format(input_file, eval, prefix)
    os.system(clans_cmd)


if __name__ == "__main__":
    out = pipeline()
    clans(out)
    exit(0)
