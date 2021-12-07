import os
import sys
import scipy.spatial.distance
import scipy.cluster.hierarchy

import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

sys.path.append('../../CLANS-Python/')
import clans.io.file_handler as fh


def get_cmap(n):
    np.random.seed(111)
    return np.random.rand(1024, 3) if n < 1024 else np.random.rand(16384, 3)


CLUSTER_COLORS = get_cmap(1024)


def cluster(env, max_dist=1e-34, filename=None, cut_dist=1e-10, min_number=10):
    distance_matrix = env.similarity_values_mtx
    distance_matrix[distance_matrix > cut_dist] = 100
    for i in range(len(env.sequences_array)):
        distance_matrix[i, i] = 0

    linear_distance = scipy.spatial.distance.squareform(distance_matrix)
    clustered = scipy.cluster.hierarchy.linkage(linear_distance)
    clusters = scipy.cluster.hierarchy.fcluster(
        clustered, t=max_dist, criterion='distance')

    # Process clusters
    cluster_ids = np.unique(clusters).tolist()
    env.groups_dict = dict()
    cluster_count = 0
    for k in cluster_ids:
        seqids = np.where(clusters == k)[0]
        rgb_aslist = CLUSTER_COLORS[k].tolist() + [1.0]
        rgb_ascol = '{},{},{},{}'.format(*[int(255*c) for c in rgb_aslist])
        rgb_assep = rgb_ascol.replace(',', ';')

        if len(seqids) < min_number:
            continue

        # Create group for cluster k 
        env.groups_dict[k] = {'name': 'cluster_{}'.format(k),
                              'type': '0',
                              'order': cluster_count,
                              'size': '9',
                              'hide': '0',
                              'color': rgb_assep,
                              'color_rgb': rgb_ascol,
                              'color_array': rgb_aslist,
                              'seqIDs': dict.fromkeys(seqids, k)}
        cluster_count += 1

    # Save groups in new file
    if filename is not None:
        fh.write_file(filename, 'clans')
    return clusters


def view_graph(clans_env, clusters, run_id=None, cutoff=1e-50):
    fig = pylab.figure()
    if clans_env.run_params['dimensions_num_for_clustering'] < 3:
        nodes = np.vstack(
            (clans_env.sequences_array['x_coor'], clans_env.sequences_array['y_coor']))
        ax = fig.add_subplot(111)
        edges = [(e[0], e[1], float(e[2])) for e in clans_env.similarity_values_list
                 if float(e[2]) < cutoff]
    else:
        nodes = np.vstack(
            (clans_env.sequences_array['x_coor'], clans_env.sequences_array['y_coor'], clans_env.sequences_array['z_coor']))
        ax = fig.add_subplot(111, projection="3d")

    ax.scatter(*nodes, s=4, c='gray', alpha=1)
    ax.scatter(*nodes, s=12, facecolors='none',
               edgecolors=CLUSTER_COLORS[clusters-1])
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    if clans_env.run_params['dimensions_num_for_clustering'] < 3:
        ax.axis('equal')
        ax.axis('off')
    else:
        ax.set_zticks([])
    pylab.title('Clusters by linkage, {}'.format(run_id))


#def cluster_at_evalues(clans_env, evalues=[1e-10, 1e-20, 1e-30, 1e-40, 1e-50]):
def cluster_at_evalues(clans_env, evalues=[1e-30]):
    for n_d, max_d in enumerate(evalues):
        clusters = cluster(clans_env, max_d, '../output/clustered_{}.clans'.format(n_d))


if __name__ == "__main__":
    fh.read_input_file(
        file_path='../output/output_nematodes_1638877235761910000.clans', file_format='clans')
    clans_env = fh.cfg
    cluster_at_evalues(clans_env)
