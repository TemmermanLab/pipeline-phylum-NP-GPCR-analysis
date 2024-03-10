Global analysis of neuropeptide receptor conservation in nematodes
==================================================================

This repository contains the code required to reproduce the results in [1].

Setup
-----

Running the Pipeline
--------------------
- `pipeline.py` contains the entry point for the pipeline. It uses the functions hhmsearch to find hits in the provided species files and hmmtop to filter out sequences with a number of transmembrane domains lower than 4.
- `clustering.py` takes the .clans file as input and finds clusters within the hits using a linkage clustering algorithm. The output can be opened wih a cluster visualization tool (e.g., clans.jar) to visualize the clusters and manually check the correctness of the clustering. If singleton are present those can be manually added to the neighboring clusters.
- `save_clusters.py` the curated .clans file can be entered into this script to save the fasta files for each cluster.
- `append_fasta.py` to each fasta file a set of 16 _C. elegans_ amine receptors are added as root of each tree.
- `postanalysis.py` generates multiple sequence alignment .fas files for each curated fasta, and uses this output to infer maximum-likelihood phylogenetic trees. Trees can be visualized with a compatible software (e.g., figtree or iTol) to check the correctness of the output. If the rooting of the tree didn't work properly, this can be adjusted at this stage by selecting the group of _C. elegans_ amine receptors and choosing "set as root" in one of these visualization softwares.
- `tree_analysis.py` each curated tree is analysed to find subclusters containing known _C. elegans_ NP-GPCRs. These clusters are annotated and coloured. The output is used to create a table on Excel where the x axis has for each cell a different nematode, and each y cell is a different _C. elegans_ NP-GPCR. If an ortholog is present a "1" is added, otherwise the corresponding cell contains a "0". The final number will reflect the number of sequences found as potential orthologs for a given _C. elegans_ NP-GPCR. We didn't look at isoforms, thus we avereged our result into a binary output, either "0" or "1"
 

References
----------
1. 