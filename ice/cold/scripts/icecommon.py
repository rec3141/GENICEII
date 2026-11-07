"""Shared helpers for the ice/cold analyses: sample design and paired ice-vs-water / TOP-vs-BOT contrasts."""
import sys, warnings; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
warnings.filterwarnings('ignore')
from common import *
C_DIR = '/data/dev/rnaseq-exploration/ice/cold/'

def design():
    m = meta().loc[samples()]
    m['layer'] = m.ice_layer.fillna('')
    return m

def blocks(m, kind):
    """list of (label, a_samples, b_samples): same season x pool x day; kind 'ice' = SI vs SW, 'top' = TOP vs BOT"""
    x = m[(m.group != 'other') & m.medium.isin(['SI', 'SW'])]
    if kind == 'top': x = x[x.layer.isin(['TOP', 'BOT'])]
    out = []
    for (se, gr, d), g in x.groupby(['season', 'group', 'day']):
        a, b = (g.index[g.medium == 'SI'], g.index[g.medium == 'SW']) if kind == 'ice' else (g.index[g.layer == 'TOP'], g.index[g.layer == 'BOT'])
        if len(a) and len(b): out.append((f'{se} {gr} d{int(d)}', list(a), list(b)))
    return out

def paired_log2(E, m, kind, pseudo=1.0):
    """E: samples x features (CPM). Returns DataFrame blocks x features of log2((mean a + p)/(mean b + p))."""
    return pd.DataFrame({lab: np.log2((E.loc[a].mean() + pseudo) / (E.loc[b].mean() + pseudo)) for lab, a, b in blocks(m, kind)}).T
