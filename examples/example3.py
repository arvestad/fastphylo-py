import fastphylo as fph

msa = fph.read('tests/data/protein_aligned.fasta')
dm = fph.distance_matrix(msa)
print(dm.to_phylip())
