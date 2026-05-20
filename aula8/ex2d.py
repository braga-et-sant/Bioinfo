from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo

names = ["S1", "S2", "S3", "S4"]

# Biopython quer a matriz em formato triangular inferior
matrix = [
    [0],
    [3, 0],
    [4, 6, 0],
    [5, 8, 9, 0]
]

dm = DistanceMatrix(names, matrix)

constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)

Phylo.draw_ascii(tree)
Phylo.write(tree, "upgma_tree.nwk", "newick")

print(tree.format("newick"))