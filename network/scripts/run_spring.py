"""SPRING (Python port, spring.py) on each layer: mclr of the node counts, latent correlation, MB neighbourhood selection with
StARS (20 subsamples of 80%, threshold 0.1, default lambda path 0.6 -> 0.006). Edge sign and weight = the symmetrised MB
coefficient (mean of the two directions; positive = positive partial association). Out: out/edges_spring_{layer}.tsv,
out/spring_stars_{layer}.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netcommon import *
import spring

for layer in (sys.argv[1:] or LAYERS):
    X, lib, N = load_layer(layer)
    fit = spring.spring(X.to_numpy(float), rep_num=20, seed=10010, verbose=True)
    A, B = fit['adj'], fit['beta']; ids = list(N.index)
    iu = np.triu_indices(len(ids), 1); ok = A[iu]
    i, j = iu[0][ok], iu[1][ok]
    E = pd.DataFrame({'a': [ids[k] for k in i], 'b': [ids[k] for k in j], 'weight': B[i, j].round(4), 'stability': fit['stability'][i, j].round(3),
                      'latent_r': fit['R'][i, j].round(3)})
    E['sign'] = np.sign(E.weight).astype(int)
    E[['a', 'b']] = [edge_key(a, b) for a, b in zip(E.a, E.b)]
    E.to_csv(os.path.join(OUT, f'edges_spring_{layer}.tsv'), sep='\t', index=False)
    ne = np.full(len(fit['lambdas']), np.nan); ne[:len(fit['path'])] = fit['path'][:, iu[0], iu[1]].sum(1)   # path stops after the StARS choice
    pd.DataFrame({'lambda': fit['lambdas'], 'instability': fit['summary'], 'edges_full_data': ne}).to_csv(
        os.path.join(OUT, f'spring_stars_{layer}.tsv'), sep='\t', index=False)
    print(layer, f'{len(ids)} nodes: StARS lambda {fit["lambdas"][fit["opt_index"]]:.4f} (index {fit["opt_index"] + 1}), {len(E)} edges, {(E.sign > 0).sum()} positive')
