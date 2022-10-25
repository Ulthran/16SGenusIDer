# GenusFinder
Given a bacterial 16S gene, infer the genus by placing it on a tree of similar sequences

## Installation

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

 - Fetch all type species 16S sequences from Tree of Life
 - Search this db for XX similar seqs
 - Build tree including the query
 - Examine nearest neighbors (subtrees snipped at high confidence intervals) for common genus
 - Examine nearest neighbors (non-topological distance between nodes) for common genus
