"""t-SNE node positions per layer (the force-directed layouts are in forcelayout.py): a t-SNE of the shape of each node's profile
across libraries. Nodes are points, libraries are features; values are the fourth root of each node's share of the library (reads or
bases per library total), then each node's row is scaled to sum 1, so two genomes with the same ups and downs sit together whatever
their overall abundance (since 2026-11-09; before, the rows were not scaled and the map mostly sorted genomes by median abundance).
Nodes with an all-zero row in a layer are left out of that layer's map. The scaled rows are reduced to their top N_PC principal
components (fewer when there are fewer libraries), then sklearn TSNE, Euclidean, PCA initialisation, fixed seed, perplexity =
round(sqrt(number of nodes)) clipped to 5-50. Out: out/layout_{layer}.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netcommon import *
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
SEED = 20261108
N_PC = 50

for layer in LAYERS:
    X, lib, N = load_layer(layer)
    P = rel_abund(X, lib, layer).T ** 0.25                                 # nodes x libraries
    P = P[P.sum(axis=1) > 0]
    P = P.div(P.sum(axis=1), axis=0)                                       # each row sums to 1: profile shape, not abundance
    k = min(N_PC, P.shape[1], P.shape[0])
    Z = PCA(n_components=k, random_state=SEED).fit_transform(P.to_numpy())
    perp = int(np.clip(round(np.sqrt(len(P))), 5, 50))
    Y = TSNE(n_components=2, perplexity=perp, init='pca', random_state=SEED, metric='euclidean').fit_transform(Z)
    Y = (Y - Y.mean(0)) / np.abs(Y - Y.mean(0)).max()
    pd.DataFrame({'x': Y[:, 0].round(4), 'y': Y[:, 1].round(4)}, index=P.index).to_csv(os.path.join(OUT, f'layout_{layer}.tsv'), sep='\t')
    print(layer, len(P), 'nodes,', k, 'principal components, perplexity', perp)
