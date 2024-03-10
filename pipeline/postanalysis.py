import os
import sys
import time
import glob

from clustering import msalign, seqtree

if __name__ == "__main__":
    working_dir = "./output/Rod"  # REPLACE HERE WITH YOUR FOLDER OF INTEREST (IT SHOULD CONTAIN THE FASTA FILES FOR EACH CLUSTER WITH APPENDED AMINE RECEPTORS USED AS ROOT OF YOUR PHYLOGENETIC TREES)
    print("Postanalysis: {}".format(working_dir))
    msalign(working_dir)
    seqtree(working_dir)
    print('Done.') 