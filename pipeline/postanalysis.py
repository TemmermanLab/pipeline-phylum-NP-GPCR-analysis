import os
import sys
import time
import glob

from clustering import msalign, seqtree

if __name__ == "__main__":
    working_dir = "../output/output_nematodes_1642085363694_clustered_0_1642087551128/Merged_files"
    print("Postanalysis: {}".format(working_dir))
    msalign(working_dir)
    seqtree(working_dir)
    print('Done.')