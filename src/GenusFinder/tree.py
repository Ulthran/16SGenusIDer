class Tree():
    def __init__(self, tree, unknown) -> None:
        self.tree = tree
        self.unknown = unknown

    # Get nearest sequences to unknown
    # @return is a list of nearest sequences
    def nearest() -> list:
        return list()

# Builds tree using RAxML
# @param seqs is a list of 16S sequences for building the tree
# @param unknown is a 16S sequence that will also be included
# @return is the Tree built by RAxML
def build_tree(seqs: list, unknown: str) -> Tree:
    return ""

# Determines the genus of the unknown
# @param tree is a Tree object to use for identification
# @return is the best guess of the genus for the unknown in tree
def identify_genus(tree: Tree) -> str:
    return ""