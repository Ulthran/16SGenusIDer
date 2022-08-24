import os
import subprocess
from collections import OrderedDict
from ete3 import Tree

# Builds tree using RAxML
# @param seqs is a list of 16S sequences for building the tree
# @param unknown is a 16S sequence that will also be included
# @return is the newick format tree built by RAxML
def build_tree(seqs: list, unknown: str) -> str:
    with open("RAxML_bestTree.genus") as f:
        return f.readline()
    
    with open("16S.db") as db, open("nearest.fasta", "w") as f:
        f.write(f">UNKNOWN\n")
        f.write(f"{unknown}\n")
        
        flag = ""
        for l in db.readlines():
            if l[0] == ">" and l.split(" ")[0][1:] in seqs:
                flag = l.strip()
            elif flag != "":
                print(flag)
                print(l)
                f.write(f"{flag}\n")
                f.write(f"{l.strip()}\n")
                flag = ""

    subprocess.run(["muscle",
    "-in", "nearest.fasta",
    "-out", "nearest_aligned.fasta"])

    subprocess.run(["raxmlHPC",
    "-m", "GTRCAT",
    "-n", "genus",
    "-p", "10000",
    "-s", "nearest_aligned.fasta"])

    
    try:
        os.remove("nearest.fasta")
        os.remove("nearest_aligned.fasta")
        os.remove("nearest_aligned.fasta.reduced")
    except OSError:
        pass

    

# Determines the genus of the unknown
# @param tree is a newick format object to use for identification
# @param unknown is the sequence to be id'ed
# @return is the best guess of the genus for the unknown in tree
def identify_genus(tree: str, unknown: str) -> str:
    print(tree)
    print(unknown)

    t = Tree(tree)
    dists = OrderedDict()
    for node in t.traverse("postorder"):
        dists[node.name] = node.get_distance("UNKNOWN")

    dists = sorted(dists.items(), key=lambda x: x)
    print(dists)

    try:
        os.remove("RAxML_bestTree.genus")
        os.remove("RAxML_info.genus")
        os.remove("RAxML_log.genus")
        os.remove("RAxML_parsimonyTree.genus")
        os.remove("RAxML_result.genus")
    except OSError:
        pass
    
    return ""