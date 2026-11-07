import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
F=D['pathways']['families']
for k in ['tmoA-F (toluene MO)','todC (benzene/toluene DO)','nahAc/Ad (naphthalene DO)','bphA (biphenyl DO)','catE/xylE (catechol 2,3-DO)','nahC/bphC (extradiol DO)','dmpKLMNOP (phenol MO)','alkB','CYP153 (product)','CANT-HYD AlkB','CANT-HYD CYP153','pmoA/mmoX (methane MO)','prmA (propane MO)']:
    for g in ['2024-25 oiled','2025-26 oiled']:
        h=F[k]['hosts'].get(g)
        if not h: continue
        cy=sum(t['pct'] for t in h['top'] if 'Cycloclasticus' in t['label'])
        print(f"{k:30s} {g} n={h['pairs']:6d} cyclo(top{len(h['top'])})={cy:5.1f}  top: "+'; '.join(f"{t['label']}[{t['host']}] {t['pct']}" for t in h['top'][:4]))
