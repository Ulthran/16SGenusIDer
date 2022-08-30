# GenusFinder
Given a bacterial 16S gene, infer the genus by placing it on a tree of similar sequences

Steps:
 - Fetch all type species 16S sequences from Tree of Life
 - Search this db for XX similar seqs
 - Build tree including the query
 - Examine nearest neighbors (subtrees snipped at high confidence intervals) for common genus
