/* Porticoccaceae 50-400-T64 (comebin_2576) section: the second big responder of the 2025-26 diesel pool.
   Data: porticocc.json (porticocc/scripts/porticocc_site.py). Uses helpers from app.js (el, fig, width, css, showTip, hideTip, fmt, fmt1, timecourse, legendPools). */
const pcSgn = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.1f')(v);
const pcPct = d3.format('.1%');
function pcSeries(J, values) {
  const out = {};
  J.samples.forEach((s, i) => { if (!['CTL', 'EXP'].includes(s.pool) || !['SW', 'SI'].includes(s.medium) || s.day == null || values[i] == null) return;
    const k = `${s.season}|${s.pool}|${s.medium}`; ((out[k] ||= {})[s.day] ||= []).push([values[i], s.name]); });
  return Object.entries(out).map(([k, byDay]) => { const [season, pool, medium] = k.split('|');
    return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a, x => x[0]), name: a.map(x => x[1]).join(', ') })) }; });
}
function pcPoints(rows, key) {
  const out = {};
  rows.forEach(x => { const v = x.v[key]; if (v == null) return; const k = `${x.season}|${x.pool}|${x.medium}`; ((out[k] ||= {})[x.day] ||= []).push(v); });
  return Object.entries(out).map(([k, byDay]) => { const [season, pool, medium] = k.split('|');
    return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a), name: `${a.length} librar${a.length > 1 ? 'ies' : 'y'}` })) }; });
}

function porticoccSec(sec) {
  const J = D.porticocc;
  if (!J) { sec.innerHTML = '<h2>Porticoccaceae 50-400-T64</h2><p class="kicker">porticocc.json is missing.</p>'; return; }
  const st = J.stats, b25 = st.baci['2025-26'];
  sec.innerHTML = `<h2>Porticoccaceae 50-400-T64</h2>
  <p class="kicker">After Cycloclasticus, the largest diesel responder of 2025-26 is an uncultured gammaproteobacterium. GTDB places it in the Porticoccaceae, genus 50-400-T64, with no species name and no cultured relative. Its genome (comebin_2576, ${d3.format('.1f')(st.size / 1e6)} Mbp, ${st.compl}% complete) has none of the classic alkane or aromatic genes. It has no alkB and no ring-hydroxylating dioxygenase, so the hydrocarbon screens on the other pages catch only one of its genes. Instead it carries a copper membrane monooxygenase (CuMMO, the enzyme family of methane and ammonia monooxygenases), three cytochrome P450s, thirteen Baeyer–Villiger monooxygenases, and a large β-oxidation and glyoxylate-shunt machinery. This page asks what it is, what it eats and how it fits beside Cycloclasticus. One pool per treatment, so contrasts are descriptive.</p>`;
  const row = el('div', { class: 'stat-row' }); sec.append(row);
  [[`${fmt(st.tpm['2025-26'].oil)} vs ${fmt(st.tpm['2025-26'].ctl)}`, 'TPM (transcripts per million), diesel vs control seawater after dosing, 2025-26 (medians)'],
   [`${pcSgn(b25.est)}`, `BACI 2025-26 (log₂ share change vs each pool’s own baseline, diesel minus control), q = ${d3.format('.2g')(b25.q)}`],
   [`${pcPct(st.tpm_max / 1e6)}`, 'of all transcripts at its peak (diesel pool, 2025-26 day 82); it is still rising when sampling stops'],
   [`${pcSgn(st.ice.delta)}`, `log₂ ice − water (CLR); lower in ice in ${st.ice.n - st.ice.n_up} of ${st.ice.n} pairs (q = ${d3.format('.2g')(st.ice.q)})`]]
    .forEach(([v, l]) => row.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  pcTiming(sec, J);
  pcDNA(sec, J);
  pcCuMMO(sec, J);
  pcToolkit(sec, J);
  pcMGE(sec, J);
  pcReading(sec, J);
}

/* ---------- timing vs Cycloclasticus ---------- */
function pcTiming(sec, J) {
  sec.append(el('h3', {}, 'A late bloomer behind Cycloclasticus'));
  sec.append(el('p', { class: 'kicker' }, 'In 2025-26, 50-400-T64 was already present in both pools before dosing: a few hundred to 1,500 TPM, slightly more in the diesel pool. It took off a few days after Cycloclasticus. Cycloclasticus sp001566825 (lorbin_184) jumped first (≈20,000 TPM by day 6–8) and levelled off at 120,000–165,000 TPM from day 17. 50-400-T64 kept climbing, from about 4,000 TPM on day 6 to 27,000 on day 20 and 70,000 on day 82. In 2024-25 it reached about 7,000 TPM at days 55–60, behind the metabat_9342 Cycloclasticus bloom. Choose a genome to compare.'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const names = { comebin_2576: 'Porticoccaceae 50-400-T64 (comebin_2576)', lorbin_184: 'Cycloclasticus sp001566825 (lorbin_184, 2025-26 bloom)', metabat_9342: 'Cycloclasticus metabat_9342 (2024-25 bloom)' };
  const sel = el('select', { id: 'pcg', 'aria-label': 'Genome' }); Object.entries(names).forEach(([k, v]) => sel.append(new Option(v, k))); c.append(el('label', { for: 'pcg' }, 'Genome'), sel);
  const v = el('div'); sec.append(v);
  function draw() { v.innerHTML = ''; const box = fig(v, `${names[sel.value]}: transcripts per million`, 'MAG share of the community transcriptome (TPM); lines join per-day means of replicate libraries.');
    legendPools(box.parentNode.insertBefore(el('div'), box)); timecourse(box, pcSeries(J, J.tpm[sel.value]), { ylab: 'TPM', unit: 'TPM', h: 200, floor: 0.5 }); }
  sel.onchange = draw; draw();
}

/* ---------- DNA ---------- */
function pcDNA(sec, J) {
  sec.append(el('h3', {}, 'Growth, not just activity'));
  sec.append(el('p', { class: 'kicker' }, 'DNA confirms the cells multiplied. In the Illumina metagenome k-mer screen, comebin_2576 rose from 0.002–0.02% of read k-mers before dosing to 0.11% by day 13 of 2025-26. In 2024-25 it rose to 0.2–0.35% at days 55–60, against under 0.02% in the control pool all season. In nanopore DNA (as a share of binned read coverage) it rose from 0.01–0.06% before dosing to 0.06–0.17% by days 4–6 of 2025-26, and reached 0.19% in the day-4 oil-snow sample. The June 2025 surface sample of the 2024-25 diesel pool, months after the experiment, still held 4× more of it than the control pool’s.'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'pcd', 'aria-label': 'DNA dataset' }); [['mg', 'Illumina metagenomes (% of read k-mers)'], ['np', 'Nanopore metagenomes (% of binned coverage)']].forEach(([k, l]) => sel.append(new Option(l, k)));
  const gsel = el('select', { id: 'pcdg', 'aria-label': 'Genome' }); [['comebin_2576', '50-400-T64'], ['lorbin_184', 'Cycloclasticus lorbin_184'], ['metabat_9342', 'Cycloclasticus metabat_9342']].forEach(([k, l]) => gsel.append(new Option(l, k)));
  c.append(el('label', { for: 'pcd' }, 'Data'), sel, el('label', { for: 'pcdg' }, 'Genome'), gsel);
  const v = el('div'); sec.append(v);
  function draw() { v.innerHTML = ''; const mg = sel.value === 'mg';
    const box = fig(v, `${gsel.options[gsel.selectedIndex].text}: ${mg ? 'Illumina k-mer share' : 'nanopore coverage share'}`, mg ? 'Percent of metagenome read k-mers matching the genome (k-mer screen, see Overview).' : 'Percent of binned nanopore read coverage.');
    legendPools(box.parentNode.insertBefore(el('div'), box)); timecourse(box, pcPoints(mg ? J.mg : J.nanopore, gsel.value), { ylab: '%', unit: '%', h: 190, floor: 0.0005 }); }
  sel.onchange = draw; gsel.onchange = draw; draw();
}

/* ---------- the CuMMO ---------- */
function pcCuMMO(sec, J) {
  sec.append(el('h3', {}, 'The copper monooxygenase'));
  const C = J.cummo;
  sec.append(el('p', { class: 'kicker' }, C ? C.summary : 'Placement of the CuMMO subunits is pending.'));
  if (C && C.hits) {
    const box = fig(sec, `Closest relatives of the CuMMO subunit A (BLASTp against a ${C.source})`, 'The query is the subunit A protein (PmoA/AmoA-like), with the Cycloclasticus metabat_9342 copy for comparison. References are type proteins of each CuMMO group, fetched from NCBI and UniProt. The top eight hits per query are ranked by bit score, with percent identity over the aligned length. Actinobacterial CuMMOs gave no significant hit.');
    box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Query</th><th>Hit</th><th>Organism</th><th class="num">Identity</th><th class="num">Aligned aa</th></tr></thead><tbody>${
      [...C.hits].sort((a, b) => (a.q.startsWith('50') ? 0 : 1) - (b.q.startsWith('50') ? 0 : 1)).map(h => `<tr><td>${h.q}</td><td style="white-space:normal">${h.desc}</td><td style="white-space:normal">${h.org}</td><td class="num">${d3.format('.0f')(h.pid)}%</td><td class="num">${h.alen}</td></tr>`).join('')}</tbody></table></div>`;
  }
  if (C && C.pairwise) sec.append(el('p', { class: 'note' }, C.pairwise));
  const ob = fig(sec, 'The CuMMO neighbourhood on contig_0004457', 'Genes in contig order with median CPM in the 2025-26 diesel seawater. The CuMMO subunits are highlighted. Upstream sit the urea-transporter operon (with urease just beyond), a CoxE/MoxR pair and an acyl-CoA dehydrogenase; downstream sit a multicopper oxidase, a ResB-like protein and an RND transporter.');
  ob.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Gene</th><th>Product</th><th>KO</th><th class="num">CPM diesel 25-26</th></tr></thead><tbody>${
    J.operon.map(g => `<tr${/onooxygenase subunit|AmoB/.test(g.p) ? ' style="background:var(--grid);font-weight:600"' : ''}><td class="mono">${g.id}</td><td style="white-space:normal">${g.p}</td><td class="mono">${g.ko || ''}</td><td class="num">${fmt1(g.oil25)}</td></tr>`).join('')}</tbody></table></div>`;
}

/* ---------- toolkit ---------- */
function pcToolkit(sec, J) {
  sec.append(el('h3', {}, 'A sub-terminal oxidation toolkit, always on'));
  sec.append(el('p', { class: 'kicker' }, 'Gene groups are ranked by their share of the genome’s own transcripts. The genome is rare outside the diesel pool, so the background (2025-26 control pool plus both pools before dosing) is pooled across samples. The striking result is how flat the profile is. Every group runs at nearly the same per-mille share before dosing, in the control pool and in diesel. The CuMMO, the Baeyer–Villiger monooxygenases, the alcohol and aldehyde dehydrogenases, β-oxidation and the glyoxylate shunt are this organism’s standing programme. Diesel changes how many cells run it, not what each cell does. Isocitrate lyase, the gateway of the glyoxylate shunt, is its single most expressed gene (about 3% of its transcripts). That is the signature of an organism growing on substrates that enter metabolism as acetyl-CoA, such as fatty acids from oxidised alkanes.'));
  const box = fig(sec, 'Functional groups: share of the genome’s transcripts', 'Per mille (‰) of comebin_2576 transcripts, pooled over samples. Diesel = 2025-26 diesel seawater after dosing; background = 2025-26 control pool and both pools before dosing; 2024-25 = diesel seawater after dosing. CPM columns are community-level medians in 2025-26 seawater, which mostly track the bloom.');
  box.innerHTML = `<div class="tablewrap" style="max-height:none"><table><thead><tr><th>Group</th><th class="num">Genes</th><th class="num">‰ diesel 25-26</th><th class="num">‰ background</th><th class="num">log₂ ratio</th><th class="num">‰ diesel 24-25</th><th class="num">CPM diesel</th><th class="num">CPM control</th></tr></thead><tbody>${
    J.toolkit.map(t => `<tr><td>${t.cls}</td><td class="num">${t.n}</td><td class="num">${fmt1(t.sh_oil)}</td><td class="num">${fmt1(t.sh_bg)}</td><td class="num">${pcSgn(t.lfc)}</td><td class="num">${fmt1(t.sh_oil24)}</td><td class="num">${fmt1(t.oil25)}</td><td class="num">${fmt1(t.ctl25)}</td></tr>`).join('')}</tbody></table></div>`;
  const c = el('div', { class: 'controls' }); sec.append(c);
  const groups = [...new Set(J.genes.map(g => g.cls))];
  const sel = el('select', { id: 'pct', 'aria-label': 'Group' }); groups.forEach(g => sel.append(new Option(g, g))); c.append(el('label', { for: 'pct' }, 'Genes in group'), sel);
  const v = el('div'); sec.append(v);
  function draw() { v.innerHTML = ''; const b = fig(v, '', 'Up to eight most-expressed genes per group, with their per-mille share of the genome’s transcripts and community CPM in 2025-26 seawater.');
    b.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Gene</th><th>Product</th><th>KO</th><th class="num">‰ diesel</th><th class="num">‰ background</th><th class="num">CPM diesel</th></tr></thead><tbody>${
      J.genes.filter(g => g.cls === sel.value).map(g => `<tr><td class="mono">${g.id}</td><td style="white-space:normal">${g.p}</td><td class="mono">${g.ko || ''}</td><td class="num">${d3.format('.2f')(g.sh_oil)}</td><td class="num">${d3.format('.2f')(g.sh_bg)}</td><td class="num">${fmt1(g.oil25)}</td></tr>`).join('')}</tbody></table></div>`; }
  sel.onchange = draw; sel.value = groups.find(g => /Baeyer/.test(g)) || groups[0]; draw();
}

/* ---------- MGEs ---------- */
function pcMGE(sec, J) {
  const M = J.mge;
  sec.append(el('h3', {}, 'Mobile elements'));
  sec.append(el('p', { class: 'kicker' }, `The CuMMO contig is ${M.cummo_contig_mobile ? 'part of an island or plasmid-like contig' : 'chromosomal'}: no IslandPath island, geNomad plasmid or provirus overlaps it. The genome has ${M.islands.length} IslandPath islands (${M.islands.map(i => `${i.id}, ${d3.format('.1f')(i.len / 1000)} kb`).join('; ')}), neither carrying hydrocarbon genes. It also has ${M.plasmid_like.length} plasmid-like contigs, ${M.plasmid_like.filter(p => p.conj).length} of them with integrative-conjugative-element transfer genes. So the oxidation toolkit looks like part of the core genome of this lineage, not recently acquired cargo.`));
}

/* ---------- reading ---------- */
function pcReading(sec, J) {
  const C = J.cummo || {};
  sec.append(el('h3', {}, 'How the pieces fit, and what is still a hypothesis'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li><b>Observed.</b> 50-400-T64 is a primed resident. It was present in both pools before dosing in 2025-26, already running the same oxidation programme. Its transcript share rose about 50-fold over its own pre-dosing level (to 36,000 TPM, against 270 in the control), and its DNA rose with it. It peaks after Cycloclasticus and is still rising at day 82.</li>
    <li><b>Observed.</b> Its oxidation genes are not induced within the cell. The per-mille shares of the CuMMO, BVMOs, P450s and β-oxidation are the same in diesel and background. Diesel selected for this organism; it did not switch the genes on.</li>
    <li><b>Hypothesis: it eats what Cycloclasticus leaves.</b> Cycloclasticus takes the aromatics first, mostly with dioxygenases. 50-400-T64 has no ring-hydroxylating dioxygenases. It has a CuMMO${C.short ? ` (${C.short})` : ''}, P450s, and many Baeyer–Villiger monooxygenases with cyclohexanone- and acetone-monooxygenase annotations, plus cyclopentanol dehydrogenase. That is the kit for oxidising alkanes and cycloalkanes in the middle of the chain (sub-terminal oxidation): alkane → secondary alcohol → ketone → ester → fatty acid. The glyoxylate shunt then assimilates the resulting acetyl-CoA. Cycloalkanes (naphthenes) make up a large share of diesel and are degraded slowly, which fits a late bloomer.</li>
    <li><b>Hypothesis: water, not ice.</b> It is depleted in ice relative to water (lower in 20 of 24 pairs) and tends to be lower in top than bottom ice (9 of 13 pairs). Oil droplets and dissolved hydrocarbons in the water column, not brine channels, are its habitat.</li>
    <li><b>To test.</b> Isolation or enrichment on cyclohexane or C₂–C₄ alkanes from stored diesel-pool water. Heterologous expression of the CuMMO. The cycloalkane fraction in GC×GC data over time.</li>`));
}

RENDER.porticocc = porticoccSec;
