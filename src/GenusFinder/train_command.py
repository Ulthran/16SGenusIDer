import argparse
from ete3 import Tree

from GenusFinder.train import learn_curve

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("type_species", help="The exact name of the type species node")

    args = p.parse_args(argv)

    with open("db/tree_LTP_all_01_2022.ntree") as f:
        t_str = ''.join(([s.strip() for s in f.readlines()]))
        learn_curve(args.type_species, Tree(t_str, format=1, quoted_node_names=True))