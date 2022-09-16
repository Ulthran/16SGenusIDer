import glob
import os
import pathlib
import shutil
import subprocess
from collections import OrderedDict
from ete3 import Tree

# Builds tree using RAxML
# @param seqs is a list of 16S sequences for building the tree
# @param unknown is a 16S sequence that will also be included
# @return is the newick format tree built by RAxML
def build_tree(seqs: list, unknown: str) -> str:
    #with open("RAxML_bipartitions.output_bootstrap.tre") as f:
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

    subprocess.run(["muscle",
    "-align", "output/nearest.fasta",
    "-output", "output/nearest_aligned.fasta"])

    full_output_path = os.path.join(pathlib.Path().resolve(), "output")

    subprocess.run(["raxmlHPC",
    "-b", "392781",
    "-N", "100",
    "-m", "GTRCAT",
    "-n", "genus1",
    "-p", "10000",
    "-s", "output/nearest_aligned.fasta",
    "-w", full_output_path])

    subprocess.run(["raxmlHPC",
    "-m", "GTRCAT",
    "-n", "genus2",
    "-p", "10000",
    "-s", "output/nearest_aligned.fasta",
    "-w", full_output_path])

    subprocess.run(["raxmlHPC",
    "-f", "b",
    "-m", "PROTGAMMAILG",
    "-n", "output_bootstrap.tre",
    "-t", "output/RAxML_bestTree.genus2",
    "-z", "output/RAxML_bootstrap.genus1",
    "-w", full_output_path])
    
    with open("output/RAxML_bipartitions.output_bootstrap.tre") as f:
        tree = f.readline()

    shutil.rmtree("output")

    return tree

# Get the genus associated with the given type species id
# @param id is the 16S id
# @return is the genus name
def get_genus(id: str) -> str:
    with open("db/type_species.fasta") as f:
        for l in f.readlines():
            if l[0] == ">" and id in l:
                return l.split("\t")[1].split(" ")[0]

# Determines the genus of the unknown
# @param tree is a newick format object to use for identification
# @param unknown is the sequence to be id'ed
# @return is the best guess of the genus for the unknown in tree
def identify_genus(tree: str, unknown: str) -> str:
    t = Tree(tree)
    #dists = OrderedDict()
    #for node in t.traverse("postorder"):
        #print(node.name)
        #node.name = get_genus(node.name)
        #print(node.name)
    #    dists[node.name] = node.get_distance("UNKNOWN")

    #dists = sorted(dists.items(), key=lambda x: x[1])
    #dists = list(filter(lambda x: x[0] != "UNKNOWN" and x[0] != "", dists))

    #for node, dist in dists:
    #    print(t.get_common_ancestor("UNKNOWN", node))

    node = t.search_nodes(name="UNKNOWN")[0]
    node = node.up
    count = 0
    while node:
        print(node)
        print(f"Dist: {node.dist}")
        print(f"Bootstrap: {node.support}")
        for n in node.traverse("postorder"):
            if(n.name != "UNKNOWN" and n.name != ""):
                print(get_genus(n.name))

        count += 1
        if count > 4:
            break
        node = node.up

    print(t.write())
    
    return ""