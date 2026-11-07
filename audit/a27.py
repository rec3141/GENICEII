import pandas as pd
F=pd.read_csv('/data/dev/rnaseq-exploration/hgt/hgt_candidates_flanks.tsv',sep='\t')
h=F[F.interpretation.str.startswith('HGT')]; print(len(h), h.level.value_counts().to_dict())
E=pd.read_csv('/data/dev/rnaseq-exploration/hgt/hc_expr_summary.tsv',sep='\t',index_col=0)
print(E.columns.tolist())
rc=[c for c in E.columns if 'read' in c.lower() or 'count' in c.lower() or c.startswith('n_')]
print(rc)
r=E['reads_total']
h=h.assign(rq=h.Q.map(r).fillna(0) if 'Q' in h else h.q.map(r).fillna(0), rt=h.t.map(r).fillna(0))
print(h[['q','t','category','rq','rt']].to_string())
print('pairs with either >=20', ((h.rq>=20)|(h.rt>=20)).sum(), 'both', ((h.rq>=20)&(h.rt>=20)).sum(), 'sum>=20', ((h.rq+h.rt)>=20).sum())
