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
    sys.path.append('../CLANS-Python/')
    import clans.io.file_handler as fh

def get_cmap(n):
    np.random.seed(111)
    return np.random.rand(1024, 3) if n < 1024 else np.random.rand(16384, 3)

def save_all_clusters_from_clans(env, clustdir):
    # Save sequences subset to FASTA
    for _, group in tqdm(clans_env.groups_dict.items(), "Saving cluster FASTA"):
        sequence_ids = group['seqIDs'].keys()
        cluster_name = group['name'].replace(os.path.sep, "-").replace(".", "_")
        sequences_subset = [SeqRecord(Seq(env.sequences_array[s_id][1]), env.sequences_array[s_id][0], description="") 
                            for s_id in sequence_ids]
        seqio.write(sequences_subset, os.path.join(clustdir, 'cluster_{}.fa'.format(cluster_name)), 'fasta')
    return True

if __name__ == "__main__":
    fpath = '../output/20240103.clans'
    fh.read_input_file(
        file_path=fpath, file_format='clans')
    clans_env = fh.cfg
    clans_env.run_params['filename'] = os.path.splitext(fpath)[0]
    if save_all_clusters_from_clans(clans_env, "../output"):
        exit(1)
    else:
        exit(0)
