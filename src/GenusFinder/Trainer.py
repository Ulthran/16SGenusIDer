import logging
import numpy as np
from ete3 import Tree
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split


class Trainer:
    def __init__(self, tree_fp: Path) -> None:
        with open(tree_fp) as f:
            self.t = Tree(f.readline())

    def train(self, training_tree: Tree, min_neighbors: int = 50):
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
