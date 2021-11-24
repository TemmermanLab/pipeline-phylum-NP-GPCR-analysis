import os
import glob
import Bio.SearchIO.HmmerIO.hmmer3_text as hmmio
import Bio.SeqIO as seqio
from Bio.Blast import NCBIXML


def pipeline(Emin_hmm=0.001, Emin_blast=1e-10, Emin_blast_post=1e-10):
    # Clean up directory
    to_delete = []
    to_delete.extend(glob.glob('./*.sea'))
    to_delete.extend(glob.glob('./*.hmm'))
    to_delete.extend(glob.glob('./*.hit'))
    to_delete.extend(glob.glob('./*.db.*'))
    to_delete.extend(glob.glob('./*.xml'))
    if len(to_delete) > 0:
        [os.remove(t) for t in to_delete]

    # BLAST Database for Cel
    cel_name = r'../cel/caenorhabditis_elegans.PRJNA13758.WBPS14.protein.fa'
    db_cmd = r"makeblastdb -dbtype prot -out cel.db -in {}".format(cel_name)
    os.system(db_cmd)

    # Build
    gpcr_name = r'rhodopsins'
    base_name = r'../cel/rhodopsins_cel_alignment.fas'
    build_cmd = r"hmmbuild {}.hmm {}".format(gpcr_name, base_name)
    os.system(build_cmd)

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
        seqio.write(fasta_contents, "{}.hit".format(
            species_names[-1]), "fasta")

        # Blast vs. Cel
        blast_cmd = r"blastp -evalue {} -query {}.hit -db cel.db -out blast_result.xml -outfmt 5".format(
            Emin_blast, species_names[-1], cel_name)
        os.system(blast_cmd)
        records = [r for r in NCBIXML.parse(open("blast_result.xml"))]

        for r in records:
            if r.alignments:
                print("query: %s" % r.query)
                for align in r.alignments:
                    for hsp in align.hsps:
                        if hsp.expect < Emin_blast_post:
                            print("match: %s " % align.title[:100])

    return True


if __name__ == "__main__":
    if pipeline():
        print('Correctly terminating pipeline.')
    exit(0)
