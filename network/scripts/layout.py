"""t-SNE node positions per layer (the force-directed layouts are in forcelayout.py): a plain t-SNE of the nodes' abundance profiles. Nodes are points, libraries are features; values are
the fourth root of each node's share of the library (reads or bases per library total), no further scaling. sklearn TSNE, Euclidean,
PCA initialisation, fixed seed, perplexity = round(sqrt(number of nodes)) clipped to 5-50. Out: out/layout_{layer}.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netcommon import *
from sklearn.manifold import TSNE
SEED = 20261108

for layer in LAYERS:
    X, lib, N = load_layer(layer)
    P = rel_abund(X, lib, layer).T ** 0.25                                 # nodes x libraries
    perp = int(np.clip(round(np.sqrt(len(P))), 5, 50))
    Y = TSNE(n_components=2, perplexity=perp, init='pca', random_state=SEED, metric='euclidean').fit_transform(P.to_numpy())
    Y = (Y - Y.mean(0)) / np.abs(Y - Y.mean(0)).max()
    pd.DataFrame({'x': Y[:, 0].round(4), 'y': Y[:, 1].round(4)}, index=P.index).to_csv(os.path.join(OUT, f'layout_{layer}.tsv'), sep='\t')
    print(layer, len(P), 'nodes, perplexity', perp)
