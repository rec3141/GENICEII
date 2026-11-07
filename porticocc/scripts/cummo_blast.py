"""CuMMO subunit A placement: local blastp of the comebin_2576 and Cycloclasticus metabat_9342 PmoA/AmoA-like proteins against a labelled
reference set (porticocc/out/cummo_refs_A.faa, fetched from NCBI/UniProt 2026-11-07), plus NCBI nr BLAST top hits if saved
(porticocc/out/nr_*.xml). Out: porticocc/out/cummo_blast.tsv, cummo_blast.json"""
import subprocess, json, re, os, glob, io
import pandas as pd, xml.etree.ElementTree as ET
O = '/data/dev/rnaseq-exploration/porticocc/out/'
LAB = {'Q607G3': ('Methylococcus capsulatus', 'pmoA, gammaproteobacterial methanotroph'), 'Q50541': ('Methylosinus trichosporium', 'pmoA, alphaproteobacterial methanotroph'),
       'SJM90451.1': ('Crenothrix polyspora', 'pmoA, Crenothrix methanotroph'), 'A9QPE3': ('Methylacidiphilum infernorum', 'pmoA, verrucomicrobial methanotroph'),
       'F6K8L1': ('Methylomonas sp. M5', 'pxmA (fragment)'), 'Q04507': ('Nitrosomonas europaea', 'amoA, betaproteobacterial ammonia oxidiser'),
       'AAC25091.1': ('Nitrosococcus oceani', 'amoA, gammaproteobacterial ammonia oxidiser'), 'WP_012215883.1': ('Nitrosopumilus', 'amoA, archaeal ammonia oxidiser'),
       'A0A0S4KUL9': ('Nitrospira inopinata', 'amoA, comammox'), 'A0ABY4QP70': ('Ca. Mycobacterium methanotrophicum', 'actinobacterial CuMMO'),
       'A0A1E3R5L4': ('Mycolicibacterium holsaticum', 'actinobacterial CuMMO'), 'A0A1J4N906': ('Nocardioides luteus', 'actinobacterial CuMMO'),
       'A0A9X2DBC6': ('Nocardioides bruguierae', 'actinobacterial CuMMO'),
       'ORU93661.1': ('Cycloclasticus symbiont of Bathymodiolus heckerae', 'Cycloclasticus CuMMO (propane/ethane/butane; Rubin-Blum et al. 2017)'),
       'OUR85102.1': ('Cycloclasticus sp. 44_32_T64', 'Cycloclasticus CuMMO, Deepwater Horizon plume simulation (Hu et al. 2017)'),
       'OUR81281.1': ('Cycloclasticus sp. 46_120_T64', 'Cycloclasticus CuMMO, Deepwater Horizon plume simulation (Hu et al. 2017)'),
       'PHS73401.1': ('Cycloclasticus sp. (MAG)', 'Cycloclasticus CuMMO'), 'TAL02096.1': ('Porticoccaceae bacterium (MAG)', 'Porticoccaceae CuMMO'),
       'OUS17543.1': ('Gammaproteobacteria bacterium 50_400_T64', '50-400-T64 type genome, Deepwater Horizon plume simulation (Hu et al. 2017)')}
def acc(h):
    m = re.match(r'(?:sp|tr)\|([^|]+)\|', h); return m.group(1) if m else h.split()[0]
qs = {}
k = None
for l in open(O + 'cummo_query.faa'):
    if l.startswith('>'): k = l[1:].strip()
    else: qs[k] = l.strip()
rows = []
for q in [x for x in qs if '_A|' in x]:
    open('/tmp/_q.faa' if False else O + '_q.faa', 'w').write(f'>{q}\n{qs[q]}\n')
    out = subprocess.run(['blastp', '-query', O + '_q.faa', '-subject', O + 'cummo_refs_A.faa', '-outfmt', '6 qseqid sseqid pident length qlen slen evalue bitscore', '-evalue', '1e-3'],
                         capture_output=True, text=True).stdout
    for line in out.strip().split('\n'):
        if not line: continue
        a = line.split('\t'); s = acc(a[1]) if '|' in a[1] else a[1]
        org, clade = LAB.get(s, (s, ''))
        rows.append(dict(query=q.split('|')[0].replace('_A', ''), ref=s, organism=org, clade=clade, pident=float(a[2]), alen=int(a[3]), qlen=int(a[4]), slen=int(a[5]), bits=float(a[7])))
os.remove(O + '_q.faa')
T = pd.DataFrame(rows).sort_values(['query', 'bits'], ascending=[True, False]).drop_duplicates(['query', 'ref'])
# the two queries against each other
open(O + '_a.faa', 'w').write(''.join(f'>{q}\n{qs[q]}\n' for q in qs if '_A|' in q))
qq = [q for q in qs if '_A|' in q]
open(O + '_b.faa', 'w').write(f'>{qq[1]}\n{qs[qq[1]]}\n'); open(O + '_c.faa', 'w').write(f'>{qq[0]}\n{qs[qq[0]]}\n')
pw = subprocess.run(['blastp', '-query', O + '_c.faa', '-subject', O + '_b.faa', '-outfmt', '6 pident length'], capture_output=True, text=True).stdout.split('\n')[0].split('\t')
for f in ('_a.faa', '_b.faa', '_c.faa'): os.remove(O + f)
T.to_csv(O + 'cummo_blast.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(T.to_string(index=False)); print('pairwise', qq, pw)
# nr BLAST (if retrieved)
nr = []
for f in sorted(glob.glob(O + 'nr_*.xml')):
    q = os.path.basename(f)[3:-4]; r = ET.parse(f).getroot()
    for h in list(r.iter('Hit'))[:6]:
        hsp = h.find('Hit_hsps/Hsp'); idn = int(hsp.find('Hsp_identity').text); al = int(hsp.find('Hsp_align-len').text)
        d = h.find('Hit_def').text; org = re.findall(r'\[([^\]]+)\]', d)
        nr.append({'q': q, 'desc': d.split(' [')[0][:90], 'org': org[0] if org else '', 'pid': round(100 * idn / al, 1), 'alen': al})
P = T[T['query'].str.startswith('Porticoccaceae')]; Cy = T[T['query'].str.startswith('Cycloclasticus')]
best = lambda d, rx: d[d.clade.str.contains(rx)].pident.max()
J = {'hits': nr or [{'q': ('50-400-T64' if r.query.startswith('Porticoccaceae') else 'Cycloclasticus metabat_9342') + ' PmoA-like', 'desc': r.clade, 'org': r.organism, 'pid': round(r.pident, 1), 'alen': int(r.alen)}
                    for r in T.groupby('query').head(8).itertuples()],
     'source': 'nr' if nr else 'local reference set',
     'pairwise': f'The 50-400-T64 and Cycloclasticus metabat_9342 subunit A proteins share {float(pw[0]):.0f}% identity over {pw[1]} amino acids with each other.',
     'best': {'p_type': round(best(P, '50-400-T64'), 1), 'p_cyclo': round(best(P, 'Cycloclasticus'), 1), 'p_pmo': round(best(P, 'pmoA'), 1), 'p_amo': round(best(P, 'amoA'), 1), 'p_actino': round(best(P, 'actinobacterial'), 1),
              'c_cyclo': round(best(Cy, 'Cycloclasticus'), 1), 'c_pmo': round(best(Cy, 'pmoA'), 1), 'c_amo': round(best(Cy, 'amoA'), 1)}}
json.dump(J, open(O + 'cummo_blast.json', 'w'), indent=1)
print(J['best'])
