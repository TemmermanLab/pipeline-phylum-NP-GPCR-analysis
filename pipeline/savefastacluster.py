import os
import sys
import time
import glob
import scipy.spatial.distance
import scipy.cluster.hierarchy

import Bio.SeqIO as seqio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.Applications import FastTreeCommandline

from tqdm import tqdm
import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

if True:
    sys.path.append('../../CLANS-Python/')
    import clans.io.file_handler as fh

def get_cmap(n):
    np.random.seed(111)
    return np.random.rand(1024, 3) if n < 1024 else np.random.rand(16384, 3)

def cluster(env, max_dist=1e-34, filename=None, cut_dist=1e-10, min_number=5):
    do_debug_clustering = False
    distance_matrix = env.similarity_values_mtx

     # Load curated fastas
    curated = [r'../curated/rhodopsins.fa', r'../curated/class_b_secretins.fa' , r'../curated/cel_protein_convertasis.fa']
    cel_fa = sum([[f.description for f in seqio.parse(curated_name, "fasta")] for curated_name in curated], [])

    # Save sequences subset to FASTA
        sequences_subset = [
            SeqRecord(Seq(env.sequences_array[s_id][1]), env.sequences_array[s_id][0], description="") for s_id in seqids]
        cel_subset = [s for s in sequences_subset if s.id in cel_fa]
        seqio.write(sequences_subset, os.path.join(
            clustdir, 'cluster_{}.fa'.format(cluster_count)), 'fasta')
        if cel_subset:
            seqio.write(cel_subset, os.path.join(clustdir, 'cel_{}.fa'.format(cluster_count)), 'fasta')

         # Save groups in new file
    if filename is not None:
        fh.write_file(filename, 'clans')
    return clusters

if __name__ == "__main__":
    fpath = '../output/output_nematodes_1649765541024.clans'
    fh.read_input_file(
        file_path=fpath, file_format='clans')
    clans_env = fh.cfg
    clans_env.run_params['filename'] = os.path.splitext(fpath)[0]
    cluster_at_evalues(clans_env)