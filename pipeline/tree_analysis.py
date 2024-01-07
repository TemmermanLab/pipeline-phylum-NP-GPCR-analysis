import os
import json
import numpy as np
import pylab
import pandas as pd
import re

from Bio import Phylo, SeqIO
from Bio.Phylo.BaseTree import BranchColor
from tqdm import tqdm


def get_cmap(n):
    np.random.seed(42)
    COLOR_LIST = [[120, 120, 120],
                  [180, 120, 120],
                  [6, 230, 230],
                  [80, 50, 50],
                  [4, 200, 3],
                  [120, 120, 80],
                  [140, 140, 140],
                  [204, 5, 255],
                  [230, 230, 230],
                  [4, 250, 7],
                  [224, 5, 255],
                  [235, 255, 7],
                  [150, 5, 61],
                  [120, 120, 70],
                  [8, 255, 51],
                  [255, 6, 82],
                  [143, 255, 140],
                  [204, 255, 4],
                  [255, 51, 7],
                  [204, 70, 3],
                  [0, 102, 200],
                  [61, 230, 250],
                  [255, 6, 51],
                  [11, 102, 255],
                  [255, 7, 71],
                  [255, 9, 224],
                  [9, 7, 230],
                  [220, 220, 220],
                  [255, 9, 92],
                  [112, 9, 255],
                  [8, 255, 214],
                  [7, 255, 224],
                  [255, 184, 6],
                  [10, 255, 71],
                  [255, 41, 10],
                  [7, 255, 255],
                  [224, 255, 8],
                  [102, 8, 255],
                  [255, 61, 6],
                  [255, 194, 7],
                  [255, 122, 8],
                  [0, 255, 20],
                  [255, 8, 41],
                  [255, 5, 153],
                  [6, 51, 255],
                  [235, 12, 255],
                  [160, 150, 20],
                  [0, 163, 255],
                  [140, 140, 140],
                  [250, 10, 15],
                  [20, 255, 0],
                  [31, 255, 0],
                  [255, 31, 0],
                  [255, 224, 0],
                  [153, 255, 0],
                  [0, 0, 255],
                  [255, 71, 0],
                  [0, 235, 255],
                  [0, 173, 255],
                  [31, 0, 255],
                  [11, 200, 200],
                  [255, 82, 0],
                  [0, 255, 245],
                  [0, 61, 255],
                  [0, 255, 112],
                  [0, 255, 133],
                  [255, 0, 0],
                  [255, 163, 0],
                  [255, 102, 0],
                  [194, 255, 0],
                  [0, 143, 255],
                  [51, 255, 0],
                  [0, 82, 255],
                  [0, 255, 41],
                  [0, 255, 173],
                  [10, 0, 255],
                  [173, 255, 0],
                  [0, 255, 153],
                  [255, 92, 0],
                  [255, 0, 255],
                  [255, 0, 245],
                  [255, 0, 102],
                  [255, 173, 0],
                  [255, 0, 20],
                  [255, 184, 184],
                  [0, 31, 255],
                  [0, 255, 61],
                  [0, 71, 255],
                  [255, 0, 204],
                  [0, 255, 194],
                  [0, 255, 82],
                  [0, 10, 255],
                  [0, 112, 255],
                  [51, 0, 255],
                  [0, 194, 255],
                  [0, 122, 255],
                  [0, 255, 163],
                  [255, 153, 0],
                  [0, 255, 10],
                  [255, 112, 0],
                  [143, 255, 0],
                  [82, 0, 255],
                  [163, 255, 0],
                  [255, 235, 0],
                  [8, 184, 170],
                  [133, 0, 255],
                  [0, 255, 92],
                  [184, 0, 255],
                  [255, 0, 31],
                  [0, 184, 255],
                  [0, 214, 255],
                  [255, 0, 112],
                  [92, 255, 0],
                  [0, 224, 255],
                  [112, 224, 255],
                  [70, 184, 160],
                  [163, 0, 255],
                  [153, 0, 255],
                  [71, 255, 0],
                  [255, 0, 163],
                  [255, 204, 0],
                  [255, 0, 143],
                  [0, 255, 235],
                  [133, 255, 0],
                  [255, 0, 235],
                  [245, 0, 255],
                  [255, 0, 122],
                  [255, 245, 0],
                  [10, 190, 212],
                  [214, 255, 0],
                  [0, 204, 255],
                  [20, 0, 255],
                  [255, 255, 0],
                  [0, 153, 255],
                  [0, 41, 255],
                  [0, 255, 204],
                  [41, 0, 255],
                  [41, 255, 0],
                  [173, 0, 255],
                  [0, 245, 255],
                  [71, 0, 255],
                  [122, 0, 255],
                  [0, 255, 184],
                  [0, 92, 255],
                  [184, 255, 0],
                  [0, 133, 255],
                  [255, 214, 0],
                  [25, 194, 194],
                  [102, 255, 0],
                  [92, 0, 255]]
    np.random.shuffle(COLOR_LIST)
    reps = n // len(COLOR_LIST) + 1
    return np.tile(np.array(COLOR_LIST), (reps, 1))[:n]


CLUSTER_COLORS = get_cmap(256)


def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]


def process_tree(fname, do_draw=False):
    tree = Phylo.read(fname, "newick")
    tree.ladderize()

    # Load curated fastas
    curated = [r'../curated/rhodopsins.fa', r'../curated/class_b_secretins.fa',
               r'../curated/cel_protein_convertasis.fa',
               r'../curated/new_GPCRs.fa']
    cel_fa = sum([[f.description for f in SeqIO.parse(curated_name, "fasta")]
                  for curated_name in curated], [])

    terminals = tree.get_terminals()
    cel_terminals = [c for c in terminals if c.name in cel_fa]
    cel_terminals_names = [c.name for c in cel_terminals]
    other_terminals = [c for c in terminals if c.name not in cel_fa]
    other_terminals_names = [o.name for o in other_terminals]
    other_terminals_clustered = dict.fromkeys(other_terminals_names, False)

    # Find split points -- ancestors between Cels
    cel_ancestors = []
    for this in tqdm(cel_terminals, desc="Find Split Points"):
        for other in cel_terminals:
            if this != other:
                common = tree.common_ancestor(this, other)
                if common != tree.root and common not in cel_ancestors:
                    cel_ancestors.append(common)

    cel_clusters = dict.fromkeys(cel_terminals_names, None)
    depths = tree.depths(cel_terminals)
    ds = [depths[c] for c in cel_terminals]
    order = np.argsort(ds)
    for idx in tqdm(order, desc="Clustering Terminals"):
        c = cel_terminals[idx]
        all_clusters = []
        cluster_color = CLUSTER_COLORS[idx].astype(np.uint8).tolist()
        cluster_color = BranchColor(*cluster_color)
        queue = tree.get_path(c)[:-1]
        while queue:
            q = queue.pop()
            if q in cel_ancestors:
                break
            cluster = []
            for o in other_terminals:
                is_clustered = other_terminals_clustered[o.name]
                ancestor = tree.common_ancestor(o, c)
                if ancestor == q and not is_clustered:
                    ancestor.color = cluster_color
                    cluster.append(o.name)
                    other_terminals_clustered[o.name] = True

            # Update terminals by keeping only those not clustered yet
            other_terminals = [
                o for o in other_terminals if o.name not in cluster]
            all_clusters.append(cluster)

        cel_clusters[c.name] = {"cluster_hierarchy": all_clusters,
                                "cluster": [c for b in all_clusters for c in b],
                                "total": len([c for b in all_clusters for c in b])}

    if do_draw:
        Phylo.draw(tree)
    
    oname = os.path.splitext(fname)[0]
    Phylo.write(tree, oname + "_colored.tree", format="phyloxml")

    with open(oname + "_colored.json", "w") as fp:
        json.dump(cel_clusters, fp)
    return cel_clusters


def generate_phylo_matrix(input):
    dbdir = r"/Users/luca/dev/GPCRs for HMMBuild/nematodes"
    fas = [os.path.join(dbdir, t) for t in os.listdir(dbdir) if ".fa" in t]
    
    cels = []
    for cluster in input:
        cels.extend(cluster.keys())
        
    cels_vs_species = dict.fromkeys(cels, None)
    
    species = []
    for fname in fas:
        sname = os.path.splitext(os.path.split(fname)[-1])[0]
        species.append(sname)

    for c in cels_vs_species:
        cels_vs_species[c] = dict.fromkeys(species, 0)

    for fname in tqdm(fas, "Searching Species"):
        sname = os.path.splitext(os.path.split(fname)[-1])[0]
        entries_per_species = [f.name for f in SeqIO.parse(fname, "fasta")]

        for cluster in input:
            current_cluster_keys = cluster.keys()
            for c in current_cluster_keys:
                cels_vs_species[c][sname] = len([i for i in cluster[c]["cluster"] if i in entries_per_species])

    cel_vs_species_list = [[cels_vs_species[c][s] for c in cels] for s in species]
    df = pd.DataFrame(cels_vs_species).transpose()
    species_fixed = ['{} {}'.format(*re.split('[_.]', s)[:2]).capitalize() for s in species]
    df.columns = species_fixed
    df = df[sorted(df)]
    df.to_excel("cel_vs_sp.xls")
    # cel_mtx = np.array(cel_vs_species_list).T
    # pylab.rcParams.update({'font.size': 6})
    # pylab.matshow(cel_mtx == 0, cmap = "gray")
    # pylab.xticks(np.arange(len(species)), species, rotation="vertical")
    # pylab.yticks(np.arange(len(cels)), cels)
    # with open("clusters_vs_species.json", "w") as fp:
    #     json.dump(input, fp)
    return input


def process_dist(tree, other_terminals, cel_terminals):
    dists = [[tree.distance(cel, t) for t in other_terminals]
             for cel in tqdm(cel_terminals)]

    pylab.matshow(dists)
    pylab.xlabel('others')
    pylab.ylabel('cel')
    pylab.tight_layout()
    pylab.savefig('all_mtx.png', dpi=600)
    pylab.show(block=True)

    dists_mtx = np.array(dists)
    best_match = np.argmin(dists)

    pylab.matshow(dists_mtx < 1)
    pylab.xlabel('others')
    pylab.ylabel('cel')
    pylab.tight_layout()
    pylab.savefig('all_bin_mtx.png', dpi=600)


if __name__ == "__main__":
    # fname = "/Users/luca/dev/GPCRs for HMMBuild/output/output_nematodes_1646735506703_clustered_0_1646740529811/correct/cluster_VP.tree"
    #rootdir = "/Users/luca/dev/GPCRs for HMMBuild/output/panphylum_correct11052022/Curated"
    #rootdir = "/Users/luca/dev/GPCRs for HMMBuild/output/29062022/curated"
    rootdir = "/Users/luca/dev/GPCRs for HMMBuild/output/new_GPCRs_noTM"
    
    # fname = r"/Users/luca/dev/GPCRs for HMMBuild/output/output_nematodes_1646735506703_clustered_0_1646740529811/correct/cluster_sec.tree"
    #rootdir = "/Users/luca/dev/GPCRs for HMMBuild/output/output_nematodes_1646735506703_clustered_0_1646740529811/correct"
    
    already_done = True
    if not already_done:
        all_trees = [t for t in os.listdir(rootdir) if ".tree" in t and "colored" not in t]
        all_trees_clusters = []
        for fname in all_trees:
            out = process_tree(os.path.join(rootdir, fname))
            all_trees_clusters.append(out)
    else:
        all_trees = [t for t in os.listdir(rootdir) if ".json" in t and "colored" in t]
        all_trees_clusters = []
        for fname in all_trees:
            with open(os.path.join(rootdir, fname), "r") as fp:
                out = json.load(fp)
                all_trees_clusters.append(out)

    generate_phylo_matrix(all_trees_clusters)