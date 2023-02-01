import csv
import numpy as np
from ete3 import Tree
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split


def learn_curve(type_species: str, t: Tree):
    n = t.search_nodes(name=type_species)[0]
    genus_name = type_species[0]
    n_iter = n.up

    while len(list(n_iter.iter_leaves())) < 30 and n_iter:
        n_iter = n_iter.up

    print(n_iter)

    if not n_iter:
        print("Oh no")

    X = [
        n.get_distance(n_loop.name)
        for n_loop in n_iter.iter_leaves()
        if n_loop.name != ""
    ]
    print(X)
    X = np.array(X).reshape(-1, 1)
    X *= 50
    y = [
        int(n_loop.name[0] == genus_name)
        for n_loop in n_iter.iter_leaves()
        if n_loop.name != ""
    ]
    print(y)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    lr = LogisticRegression()
    lr.fit(X_train, y_train)

    y_true, y_pred = y_test, lr.predict(X_test)
    print(classification_report(y_true, y_pred))

    # print("\n")
    # print(lr.coef_.ravel())
    # print(lr.intercept_.ravel())

    return lr


def probability_for(lr, dist):
    return lr.predict_proba(dist)
