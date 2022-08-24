import argparse

from .download import 
from .db import create_db, find_similar
from .tree import build_tree, identify_genus

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("seq", help="the 16S sequence to be identified")
    p.add_argument("--api_key", help="your NCBI API key for increasing the speed of requests when creating a db (OPTIONAL)")

    args = p.parse_args(argv)


    create_db(args.api_key)
    identify_genus(build_tree(find_similar(args.seq), args.seq), args.seq)
