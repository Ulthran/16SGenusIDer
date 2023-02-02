# GenusFinder

<!-- Badges start -->
[![Tests](https://github.com/Ulthran/GenusFinder/actions/workflows/tests.yml/badge.svg)](https://github.com/Ulthran/GenusFinder/actions/workflows/tests.yml)
[![Super-Linter](https://github.com/Ulthran/GenusFinder/actions/workflows/linter.yml/badge.svg)](https://github.com/Ulthran/GenusFinder/actions/workflows/linter.yml)
[![Upload Python Package](https://github.com/Ulthran/GenusFinder/actions/workflows/python-publish.yml/badge.svg)](https://github.com/Ulthran/GenusFinder/actions/workflows/python-publish.yml)
<!-- Badges end -->

Given a bacterial 16S gene, infer the genus by placing it on a tree of similar sequences

## Installation

GenusFinder should be run on a computer with at least XX RAM and XX CPUs. To install,

```
NOT PUBLISHED YET
```

To install the dev version of GenusFinder,

```
git clone https://github.com/Ulthran/GenusFinder.git
cd GenusFinder
conda env create --file genusfinder_env.yaml
conda activate genusfinder
pip install .
```

## Running

```
idgenus ATCGATCGATCGATCG...GCTACTATACGA
```

## Steps

NEW

 - Fetch LTP tree and alignment
 - Align query sequence with LTP alignment
 - Add the query sequence to the LTP tree using the combined alignment
 - List nearby genuses to query and train curves on those genuses
 - Apply curves to tree to get probabilities of query being in those genuses

OLD

 - Fetch all type species 16S sequences from Tree of Life
 - Search this db for XX similar seqs
 - Build tree including the query
 - Examine nearest neighbors (subtrees snipped at high confidence intervals) for common genus
 - Examine nearest neighbors (non-topological distance between nodes) for common genus
