import os
import sys
import time
import glob

import Bio.SeqIO as seqio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from tqdm import tqdm
import numpy as np

if True:
    sys.path.append('../../CLANS-Python/')
    import clans.io.file_handler as fh

def append_to_all(fasta_dir, src_fasta_path, dst_fasta_path):
    flist = glob.glob(fasta_dir + "/*.fa")
    fasta = [f for f in seqio.parse(src_fasta_path, "fasta")]
    for f0 in flist:
        all_seqs = [f for f in seqio.parse(f0, "fasta")]
        all_seqs.extend()

        seqio.write(all_seqs, 'output_nematodes.fa', 'fasta')
    # seqio.write(sequences_subset, os.path.join(
    #     clustdir, 'cluster_{}.fa'.format(cluster_count)), 'fasta')
    
    # Save groups in new file
    # 
    return None


if __name__ == "__main__":
    fpath = '../output/panphylum_correct11052022_39groups'
    srcpath = '../curated/amine_rec.fa'
    append_to_all(fpath, srcpath, fpath)
