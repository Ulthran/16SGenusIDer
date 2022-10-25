import argparse
import os

from .db import create_db, find_similar
from .tree import build_tree, identify_genus

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--seq", help="the 16S sequence to be identified")
    p.add_argument("--id", help="the identity value to use with vsearch", default="0.9")
    p.add_argument("--keep_output", help="switch to not remove output files", action="store_true")

    args = p.parse_args(argv)

    if not os.path.exists("output/"):
        os.makedirs("output/")

    identify_genus(build_tree(find_similar(args.seq, args.id), args.seq, args.keep_output), args.seq)
