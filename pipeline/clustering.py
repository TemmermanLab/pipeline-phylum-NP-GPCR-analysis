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
import numpy as np

if True:
    sys.path.append("./CLANS-Python/")
    import clans.io.file_handler as fh


def get_cmap(n):
    np.random.seed(111)
    return np.random.rand(1024, 3) if n < 1024 else np.random.rand(16384, 3)


CLUSTER_COLORS = get_cmap(1024)


def msalign(clusters_fasta_dir):
    all_fasta = glob.glob(clusters_fasta_dir + "/cluster*.fa")
    print("### Multiple Sequence Alignment ###")
    for f in tqdm(all_fasta, desc="Aligning Sequences"):
        alignment_path = os.path.splitext(f)[0] + ".fas"
        clustalw_cmd = ClustalwCommandline(
            infile=f, outfile=alignment_path, output="FASTA"
        )
        clustalw_cmd()


def seqtree(clusters_fasta_dir):
    all_fastas = glob.glob(clusters_fasta_dir + "/*.fas")
    print("### Tree per Cluster ###")
    for f in tqdm(all_fastas, desc="Generating Trees"):
        alignment_path = f
        tree_path = os.path.splitext(f)[0] + ".tree"
        tree_cmd = FastTreeCommandline(
            cmd="fasttree",
            input=alignment_path,
            out=tree_path,
            spr=4,
            mlacc=2,
            slownni=True,
        )
        tree_cmd()


def cluster(env, max_dist=1e-34, filename=None, cut_dist=1e-10, min_number=5):
    print("### Linkage Clustering ###")
    distance_matrix = env.similarity_values_mtx

    # All values bigger than cut_dist are set to maximum (100)
    distance_matrix[distance_matrix > cut_dist] = 100
    for i in range(len(env.sequences_array)):
        distance_matrix[i, i] = 0

    # Form dendrogram
    linear_distance = scipy.spatial.distance.squareform(distance_matrix)
    clustered = scipy.cluster.hierarchy.linkage(linear_distance)

    # Cut based on max_dist
    clusters = scipy.cluster.hierarchy.fcluster(
        clustered, t=max_dist, criterion="distance"
    )  # Stop at cophenetic distance max_dist
    cluster_ids = np.unique(clusters).tolist()
    cluster_sizes = [np.sum(clusters == k) for k in cluster_ids]

    # Process clusters
    large_clusters = [
        cluster_ids[q]
        for q in range(len(cluster_ids))
        if cluster_sizes[q] >= min_number
    ]
    n_clusters = len(large_clusters)

    ordered_id = np.argsort(-np.array(cluster_sizes))
    env.groups_dict = dict()
    cluster_count = 0

    clustdir = os.path.splitext(filename)[0]
    os.makedirs(clustdir)

    # Load curated fastas
    curated = [
        r"./curated/rhodopsins.fa",
        r"./curated/class_b_secretins.fa",
    ]
    cel_fa = sum(
        [
            [f.description for f in seqio.parse(curated_name, "fasta")]
            for curated_name in curated
        ],
        [],
    )

    marked_for_deletion = []
    for o in tqdm(ordered_id, "Processing Clusters"):
        k = cluster_ids[o]
        seqids = np.where(clusters == k)[0]
        rgb_aslist = CLUSTER_COLORS[k].tolist() + [1.0]
        rgb_ascol = "{},{},{},{}".format(*[int(255 * c) for c in rgb_aslist])
        rgb_assep = rgb_ascol.replace(",", ";")

        if len(seqids) < min_number:
            marked_for_deletion.append(seqids.tolist())
            continue

        # Save sequences subset to FASTA
        sequences_subset = [
            SeqRecord(
                Seq(env.sequences_array[s_id][1]),
                env.sequences_array[s_id][0],
                description="",
            )
            for s_id in seqids
        ]
        cel_subset = [s for s in sequences_subset if s.id in cel_fa]
        seqio.write(
            sequences_subset,
            os.path.join(clustdir, "cluster_{}.fa".format(cluster_count)),
            "fasta",
        )
        if cel_subset:
            seqio.write(
                cel_subset,
                os.path.join(clustdir, "cel_{}.fa".format(cluster_count)),
                "fasta",
            )

        # Create group for cluster k
        env.groups_dict[k] = {
            "name": "cluster_{}".format(cluster_count),
            "type": "0",
            "order": cluster_count,
            "size": "9",
            "hide": "0",
            "color": rgb_assep,
            "color_rgb": rgb_ascol,
            "color_array": rgb_aslist,
            "seqIDs": dict.fromkeys(seqids, k),
        }
        cluster_count += 1

    # Save groups in new file
    if filename is not None:
        fh.write_file(filename, "clans")
    return clusters


def cluster_at_evalues(clans_env, thresholds=[1e-30]):
    tstamp = int(time.time_ns() / 1e6)
    for n_d, max_d in enumerate(thresholds):
        run_fname = clans_env.run_params["filename"] + "_clustered_{}_{}".format(
            n_d, tstamp
        )
        cluster(clans_env, max_d, run_fname + ".clans")


def save_cluster_fastas(env):
    filename = clans_env.run_params["filename"]
    clustdir = os.path.splitext(filename)[0]
    os.makedirs(clustdir, exist_ok=True)

    for _, group in tqdm(clans_env.groups_dict.items(), "Saving cluster FASTA"):
        # Save sequences subset to FASTA
        seqids = group["seqIDs"].keys()
        sequences_subset = [
            SeqRecord(
                Seq(env.sequences_array[s_id][1]),
                env.sequences_array[s_id][0],
                description="",
            )
            for s_id in seqids
        ]
        seqio.write(
            sequences_subset,
            os.path.join(clustdir, "{}.fa".format(group["name"].replace("/", "_"))),
            "fasta",
        )
    return None


if __name__ == "__main__":
    fpath = (
        "./output/output_nematodes_1708274274759.clans"  # REPLACE HERE YOUR CLANS FILE
    )
    is_clustered = False
    fh.read_input_file(file_path=fpath, file_format="clans")
    clans_env = fh.cfg
    clans_env.run_params["filename"] = os.path.splitext(fpath)[0]
    if is_clustered:
        save_cluster_fastas(clans_env)
    else:
        cluster_at_evalues(clans_env)
