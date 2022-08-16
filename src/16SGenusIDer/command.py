import argparse

from .db import create_db, update_db, find_similar
from .tree import build_tree, identify_genus

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("seq", help="the 16S sequence to be identified")
    p.add_argument("--db", help="the path to a directory containing an existing file named '16S.db' or the filepath itself")

    args = p.parse_args(argv)

    create_db()
    identify_genus(build_tree(find_similar(args.seq)), args.seq)
