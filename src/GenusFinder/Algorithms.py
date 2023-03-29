import logging
import numpy as np
from collections import OrderedDict
from ete3 import Tree
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split


class Algorithms:
    """
    Manages the probability calculations for subtree and full tree methods
    Main methods are:
    1) distance_probs which calculates subtree distance-based probabilities
    2) bootstrap_probs which calculates subtree bootstrap-based probabilities
    3) train which calculates full tree LinReg-based probabilities
    """

    def __init__(
        self, tree_fp: Path, type_species_fp: Path, query: Path = None
    ) -> None:
        with open(tree_fp) as f:
            self.t = Tree(f.readline())

        self.type_species_fp = type_species_fp
        self.lookup = {}  # Store commonly accessed type species

        with open(query) as f:
            self.query = f.readline().strip()

    def distance_probs(self) -> OrderedDict:
        """
        Calculate distance-based genus probabilities
        Iterate through each node and take the inverse of the distance to the UNKNOWN
        as the addition it makes to the probability
        """
        logging.info("Starting distance-based probability calculations...")
        dist_prob = {}
        for n in self.t.traverse():
            logging.debug(f"Node: {n.name}")
            if (
                n.name != "" and n.name != "UNKNOWN"
            ):  # Nodes with "" for a name are not leaves
                d = self.t.get_distance("UNKNOWN", n.name)
                logging.debug(f"Distance from unknown: {d}")
                g = self.get_genus(n.name, self.lookup)
                self.lookup[n.name] = g
                dist_prob[g] = dist_prob[g] + (1 / d) if g in dist_prob else (1 / d)

        dist_prob = {
            k: v / sum(dist_prob.values()) for k, v in dist_prob.items()
        }  # Normalize probabilities
        dist_prob = OrderedDict(sorted(dist_prob.items(), key=lambda x: -x[1]))
        return dist_prob

    def bootstrap_probs(self) -> OrderedDict:
        """
        Calculate bootstrap-based probabilities
        Iterate up through subtrees, starting with the smallest one containing the UNKNOWN, and at each level
        the bootstrap value determines how much of the remaining total probability is claimed by the genuses
        there
        E.g. if the first subtree has a bootstrap of 93, we can be pretty sure the UNKNOWN should be on the same
        branch as the others in this tree so their contributions will be 93% of the total probabilities. If the
        next branch is 52, it will account for 52% of the remaining 7% of the total probabilities, and so on.
        """
        logging.info("Starting bootstrap-based probability calculations...")
        boot_prob = {}
        remaining_frac = 100

        node = self.t.search_nodes(name="UNKNOWN")[0]
        node = node.up
        count = 0
        while node:
            logging.debug(f"Node: {node.name}")
            logging.debug(f"Dist: {node.dist}")
            logging.debug(f"Bootstrap: {node.support}")
            factor = node.support * remaining_frac
            sub_dict = {}
            for n in node.traverse():
                if n.name != "UNKNOWN" and n.name != "":
                    g = self.get_genus(n.name, self.lookup)
                    sub_dict[g] = sub_dict[g] + 1 if g in sub_dict else 1

            sub_dict = {k: v / sum(sub_dict.values()) for k, v in sub_dict.items()}
            for k, v in sub_dict.items():
                boot_prob[k] = (
                    boot_prob[k] + factor * v if k in boot_prob else factor * v
                )
            remaining_frac = remaining_frac - (node.support / 100) * remaining_frac

            count += 1
            if count > 4:
                break
            node = node.up

        boot_prob = {
            k: v / sum(boot_prob.values()) for k, v in boot_prob.items()
        }  # Normalize probabilities
        boot_prob = OrderedDict(sorted(boot_prob.items(), key=lambda x: -x[1]))
        return boot_prob

    def train(self, training_tree: Tree, min_neighbors: int = 50) -> dict:
        ts = [
            s for s in self.get_nearby_species(min_neighbors) if self.is_type_species(s)
        ]
        lrs = {s: self.learn_curve(s, training_tree) for s in ts}
        probs = {
            s: self.probability_for(lr, self.distance_to_unknown)
            for s, lr in lrs.items()
        }
        logging.info(f"{probs}")
        return probs

    def get_genus(self, id: str, lookup: dict = None) -> str:
        if lookup:
            if id in lookup:
                return lookup[id]

        with open("db/type_species.fasta") as f:
            for l in f.readlines():
                if l[0] == ">" and id in l:
                    return l.split("\t")[1].split(" ")[0]

    def get_nearby_species(self, min_neighbors: int) -> list:
        node = self.t.search_nodes(name="UNKNOWN")[0]
        while len(list(node.iter_leaves())) < min_neighbors and node:
            node = node.up

        return [n.name for n in node.iter_leaves()]

    def is_type_species(self, species: str) -> bool:
        return True

    def distance_to_unknown(self, species: str) -> float:
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

    @staticmethod
    def learn_curve(type_species: str, t: Tree):
        n = t.search_nodes(name=type_species)[0]
        genus_name = type_species[0]
        n_iter = n.up

        while len(list(n_iter.iter_leaves())) < 30 and n_iter:
            n_iter = n_iter.up

        logging.info(genus_name)
        logging.debug(n_iter)

        if not n_iter:
            logging.error("Oh no")

        X = [
            n.get_distance(n_loop.name)
            for n_loop in n_iter.iter_leaves()
            if n_loop.name != ""
        ]
        logging.debug(X)
        X = np.array(X).reshape(-1, 1)
        X *= 50
        y = [
            int(n_loop.name[0] == genus_name)
            for n_loop in n_iter.iter_leaves()
            if n_loop.name != ""
        ]
        logging.debug(y)

        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )

        lr = LogisticRegression()
        lr.fit(X_train, y_train)

        y_true, y_pred = y_test, lr.predict(X_test)
        logging.info(classification_report(y_true, y_pred))

        logging.debug(lr.coef_.ravel())
        logging.debug(lr.intercept_.ravel())

        return lr

    @staticmethod
    def probability_for(lr, dist):
        return lr.predict_proba(dist)
