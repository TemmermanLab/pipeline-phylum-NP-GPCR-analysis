import sys
import glob

import Bio.SeqIO as seqio

if True:
    sys.path.append('./CLANS-Python/')
    import clans.io.file_handler as fh

def append_to_all(fasta_dir, src_fasta_path, dst_fasta_path):
    flist = glob.glob(fasta_dir + "/*.fa")
    fasta = [f for f in seqio.parse(src_fasta_path, "fasta")]
    for f0 in flist:
        all_seqs = [f for f in seqio.parse(f0, "fasta")]
        all_seqs.extend()

        seqio.write(all_seqs, 'output_nematodes.fa', 'fasta')
    return None


if __name__ == "__main__":
    fpath = './output/panphylum_correct11052022_39groups'
    srcpath = './curated/amine_rec.fa'
    append_to_all(fpath, srcpath, fpath)
