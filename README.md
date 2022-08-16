# 16SGenusIDer
Given a bacterial 16S gene, infer the genus by placing it on a tree of similar sequences

Steps:
 - Fetch all RefSeq 16S seqs from NCBI
 - Search this db for XX similar seqs
 - Build tree including the query
 - Determine genuses of nearest seqs