import fastphylo as e

seqs = fph.read('tests/data/protein_unaligned.fasta')
msa = fph.align(seqs)
dm = fph.distance_matrix(msa)
t_fnj = fph.fnj(dm)
t_nj = fph.nj(dm)
t_bionj = fph.bionj(dm)

for m, t in [('FNJ', t_fnj), ('NJ', t_nj), ('BIONJ', t_bionj)]:
    print(f'[{m}]', t.to_newick())
