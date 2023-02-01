import glob
import os
import shutil
import subprocess
from ete3 import Tree
from pathlib import Path

from train import learn_curve, probability_for
from Aligner import MuscleAligner


def insert_on_LTP_tree(seq: str) -> Tree:
    with open("output/query.fasta", "w") as f:
        f.write(f">UNKNOWN\n")
        f.write(f"{seq}\n")

    aligner = MuscleAligner()
    aligner.call(
        True,
        Path("db/LTP_01_2022_aligned.fasta"),
        Path("output/query.fasta"),
        Path("output/combined_aligned.fasta"),
    )

    "sed 's/U/T/g' output/combined_aligned.fasta > output"  # Replace Us with Ts for RAxML??
    "sed 's/\./-/g' output/combined_aligned.fasta > output"  # '.'s preceed alignemnts?

    full_output_path = os.path.join(pathlib.Path().resolve(), "output")

    subprocess.run(
        [
            "raxmlHPC",
            "-m",
            "GTRCAT",
            "-n",
            "combined",
            "-p",
            "10000",
            "-f",
            "y",  # "-f v" would give a more robust answer but take longer, might be worth it for just one seq to insert
            "-s",
            "output/combined_aligned.fasta",
            "-t",
            "output/tree_LTP_all_01_2022.ntree",
            "-w",
            full_output_path,
        ]
    )


def get_nearby_species() -> list:
    with open("output/RAxML_bestTree.placeQ") as f:
        t = Tree(f.readline())

    node = t.search_nodes(name="UNKNOWN")[0]
    while len(list(node.iter_leaves())) < 50 and node:
        node = node.up

    return [n.name for n in node.iter_leaves()]


def is_type_species(species: str) -> bool:
    return True


def distance_to_unknown(species: str) -> float:
    return 1


def determine_probabilities(seq: str):
    insert_on_LTP_tree(seq)
    with open("output/RAxML_bestTree.placeQ") as f:
        t = Tree(f.readline())
    nearby_species = get_nearby_species()
    nearby_type_species = [name for name in nearby_species if is_type_species(name)]
    lrs = {name: learn_curve(name, t) for name in nearby_type_species}
    distances = {name: distance_to_unknown(name) for name in nearby_type_species}
    probs = {name: probability_for(lr, distances[name]) for name, lr in lrs}
    print(probs)


# Builds tree using RAxML
# @param seqs is a list of 16S sequences for building the tree
# @param unknown is a 16S sequence that will also be included
# @param keep_true determines whether or not to delete the output files
# @return is the newick format tree built by RAxML
def build_tree(seqs: list, unknown: str, keep_output: bool) -> str:
    # with open("RAxML_bipartitions.output_bootstrap.tre") as f:
    #    return f.readline()

    with open("db/type_species.fasta") as db, open("output/nearest.fasta", "w") as f:
        f.write(f">UNKNOWN\n")
        f.write(f"{unknown}\n")

        flag = ""
        for l in db.readlines():
            if l[0] == ">" and l.split("\t")[0][1:] in seqs:
                flag = l.strip()
            elif flag != "":
                f.write(f"{flag}\n")
                f.write(f"{l.strip()}\n")
                flag = ""

    subprocess.run(
        [
            "muscle",
            "-align",
            "output/nearest.fasta",
            "-output",
            "output/nearest_aligned.fasta",
        ]
    )

    full_output_path = os.path.join(pathlib.Path().resolve(), "output")

    subprocess.run(
        [
            "raxmlHPC",
            "-b",
            "392781",
            "-N",
            "100",
            "-m",
            "GTRCAT",
            "-n",
            "genus1",
            "-p",
            "10000",
            "-s",
            "output/nearest_aligned.fasta",
            "-w",
            full_output_path,
        ]
    )

    subprocess.run(
        [
            "raxmlHPC",
            "-m",
            "GTRCAT",
            "-n",
            "genus2",
            "-p",
            "10000",
            "-s",
            "output/nearest_aligned.fasta",
            "-w",
            full_output_path,
        ]
    )

    subprocess.run(
        [
            "raxmlHPC",
            "-f",
            "b",
            "-m",
            "PROTGAMMAILG",
            "-n",
            "output_bootstrap.tre",
            "-t",
            "output/RAxML_bestTree.genus2",
            "-z",
            "output/RAxML_bootstrap.genus1",
            "-w",
            full_output_path,
        ]
    )

    with open("output/RAxML_bipartitions.output_bootstrap.tre") as f:
        tree = f.readline()

    if not keep_output:
        shutil.rmtree("output")

    return tree


# Get the genus associated with the given type species id
# @param id is the 16S id
# @param lookup is a dictionary containing already-looked-up ids
# NOTE: lookup should be updated by the caller of this function
# @return is the genus name
def get_genus(id: str, lookup: dict = None) -> str:
    if lookup:
        if id in lookup:
            return lookup[id]

    with open("db/type_species.fasta") as f:
        for l in f.readlines():
            if l[0] == ">" and id in l:
                return l.split("\t")[1].split(" ")[0]


# Display probabilities
# @param d is the probability list
def display_probs(d: list):
    for k, p in d:
        print(f"{k}: {round(p * 100, 5)}%")
    print("\n")


# Determines the genus of the unknown
# @param tree is a newick format object to use for identification
# @param unknown is the sequence to be id'ed
# @return is the best guess of the genus for the unknown in tree
def identify_genus(tree: str, unknown: str) -> str:
    t = Tree(tree)
    lookup = dict()

    # Calculate distance-based genus probabilities
    # Iterate through each node and take the inverse of the distance to the UNKNOWN
    # as the addition it makes to the probability
    dist_prob = dict()
    for n in t.traverse():
        if (
            n.name != "" and n.name != "UNKNOWN"
        ):  # Nodes with "" for a name are not leaves
            d = t.get_distance("UNKNOWN", n.name)
            g = get_genus(n.name, lookup)
            lookup[n.name] = g
            dist_prob[g] = dist_prob[g] + (1 / d) if g in dist_prob else (1 / d)

    dist_prob = {
        k: v / sum(dist_prob.values()) for k, v in dist_prob.items()
    }  # Normalize probabilities
    dist_prob = sorted(dist_prob.items(), key=lambda x: -x[1])

    # Calculate bootstrap-based probabilities
    # Iterate up through subtrees, starting with the smallest one containing the UNKNOWN, and at each level
    # the bootstrap value determines how much of the remaining total probability is claimed by the genuses
    # there
    # E.g. if the first subtree has a bootstrap of 93, we can be pretty sure the UNKNOWN should be on the same
    # branch as the others in this tree so their contributions will be 93% of the total probabilities. If the
    # next branch is 52, it will account for 52% of the remaining 7% of the total probabilities, and so on.
    boot_prob = dict()
    remaining_frac = 100

    node = t.search_nodes(name="UNKNOWN")[0]
    node = node.up
    count = 0
    while node:
        print(node)
        print(f"Dist: {node.dist}")
        print(f"Bootstrap: {node.support}")
        factor = node.support * remaining_frac
        sub_dict = dict()
        for n in node.traverse():
            if n.name != "UNKNOWN" and n.name != "":
                g = get_genus(n.name, lookup)
                sub_dict[g] = sub_dict[g] + 1 if g in sub_dict else 1

        sub_dict = {k: v / sum(sub_dict.values()) for k, v in sub_dict.items()}
        for k, v in sub_dict.items():
            boot_prob[k] = boot_prob[k] + factor * v if k in boot_prob else factor * v
        remaining_frac = remaining_frac - (node.support / 100) * remaining_frac

        count += 1
        if count > 4:
            break
        node = node.up

    boot_prob = {
        k: v / sum(boot_prob.values()) for k, v in boot_prob.items()
    }  # Normalize probabilities
    boot_prob = sorted(boot_prob.items(), key=lambda x: -x[1])

    print("Distance based probabilities:\n")
    display_probs(dist_prob)
    print("Bootstrap based probabilities:\n")
    display_probs(boot_prob)
