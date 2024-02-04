import os
import sys
import time
import glob

from clustering import msalign, seqtree

if __name__ == "__main__":
    working_dir = "../output/37"
    print("Postanalysis: {}".format(working_dir))
    msalign(working_dir)
    seqtree(working_dir)
    print('Done.') 