/* Nitrogen section: N limitation, nitrifiers, denitrification and copper monooxygenases in the diesel pools.
   Data: nitrogen.json (nitrogen/scripts/nitrogen_site.py). Uses helpers from app.js (el, fig, width, css, showTip, hideTip, fmt, fmt1, timecourse, legendPools). */
const nSgn = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.1f')(v);
const nPct = d3.format('.0%');
const nCat = (J, frag) => J.categories.find(c => c.name.includes(frag));
function nSeries(J, values) {
  const out = {};
  J.samples.forEach((s, i) => { if (!['CTL', 'EXP'].includes(s.pool) || !['SW', 'SI'].includes(s.medium) || s.day == null || values[i] == null) return;
    const k = `${s.season}|${s.pool}|${s.medium}`; ((out[k] ||= {})[s.day] ||= []).push([values[i], s.name]); });
  return Object.entries(out).map(([k, byDay]) => { const [season, pool, medium] = k.split('|');
    return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a, x => x[0]), name: a.map(x => x[1]).join(', ') })) }; });
}

function nitrogenSec(sec) {
  const J = D.nitrogen;
  if (!J) { sec.innerHTML = '<h2>Nitrogen</h2><p class="kicker">nitrogen.json is missing.</p>'; return; }
  const nrt = nCat(J, 'transporters'), nas = nCat(J, 'nasA'), amt = nCat(J, 'amtB'), amo = nCat(J, 'nitrifier hosts'), cu = nCat(J, 'CuMMO'), nar = nCat(J, 'narGHI'), cyn = nCat(J, 'cynS');
  sec.innerHTML = `<h2>Nitrogen</h2>
  <p class="kicker">Diesel is almost pure carbon. Bacteria that grow on it need nitrogen from the water, so an oil spill should push the community towards nitrogen limitation. This page follows 30 nitrogen gene groups through both winters: nitrogen uptake and its regulators, nitrification, denitrification, DNRA, anammox and N₂ fixation. The clearest signal is a large, diesel-specific nitrogen-scavenging response in the hydrocarbon degraders and their partners. Ammonia oxidisers are pushed down. Nitrate is used for growth and partly respired, but not fully denitrified. Expression is CPM of protein-coding reads, with gene groups from the coordinate-corrected KOs and hosts from MAGs or, for unbinned contigs, Kaiju. No nutrient concentrations were measured in the pools (filtrate was collected from 2025-26 onward but not analysed), so this rests on transcripts alone. One pool per treatment, so all contrasts are descriptive.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  [[`${fmt(nrt.oil25)} vs ${fmt(nrt.ctl25)}`, `CPM of nitrate/nitrite transporters, diesel vs control seawater 2025-26 (BACI ${nSgn(nrt.baci25)} log₂)`],
   [`${nSgn(amt.lfc25)} / ${nSgn(nas.lfc25)}`, 'log₂ diesel ÷ control 2025-26: ammonium transporter / assimilatory nitrate reductase'],
   [`${nSgn(amo.lfc24)}`, `log₂ diesel ÷ control 2024-25: ammonia monooxygenase of nitrifiers (BACI ${nSgn(amo.baci24)})`],
   [`${nSgn(cu.lfc25)}`, 'log₂ diesel ÷ control 2025-26: copper monooxygenases of hydrocarbon degraders (not nitrifiers)']]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  nTable(sec, J);
  nExplorer(sec, J);
  nStress(sec, J);
  nGenomes(sec, J);
  nReading(sec, J);
}

/* ---------- inventory table ---------- */
function nTable(sec, J) {
  sec.append(el('h3', {}, 'The whole inventory at a glance'));
  const box = fig(sec, '', 'Median CPM in seawater after dosing. Ratios are log₂ (diesel ÷ control) of medians. BACI is the change from each pool’s own pre-dosing samples, diesel minus control; 2024-25 has a single pre-dosing sample per pool. Ice is paired ice minus water (24 pairs, both winters); Top is top minus bottom ice (13 pairs, 2024-25 and 2025-26). Cell colour runs from blue (higher in control, water or bottom ice) to orange (higher in diesel, ice or top ice). Hosts are the genomes carrying the most diesel-pool transcripts. Click a row to see its time course.');
  const col = v => { if (v == null) return 'transparent'; const t = Math.max(-1, Math.min(1, v / 4)); return t >= 0 ? d3.interpolateRgb(css('--panel'), css('--exp'))(t * 0.85) : d3.interpolateRgb(css('--panel'), css('--ctl'))(-t * 0.85); };
  const cell = v => `<td class="num" style="background:${col(v)}">${nSgn(v)}</td>`;
  let html = `<div class="tablewrap" style="max-height:none"><table><thead><tr><th>Gene group · main hosts in the diesel pool</th><th class="num">Genes</th><th class="num">CPM diesel 24-25</th><th class="num">CPM diesel 25-26</th><th class="num">24-25</th><th class="num">25-26</th><th class="num">BACI 24-25</th><th class="num">BACI 25-26</th><th class="num">Ice</th><th class="num">Top</th></tr></thead><tbody>`;
  J.groups.forEach(g => {
    const cs = J.categories.filter(c => c.group === g); if (!cs.length) return;
    html += `<tr><td colspan="10" style="background:var(--grid);font-weight:650">${g}</td></tr>`;
    cs.forEach(c => {
      html += `<tr data-c="${c.name}"><td style="white-space:normal;min-width:300px">${c.name}<div style="font-size:11.5px;color:var(--ink-3);margin-top:2px">${
        c.hosts.diesel.slice(0, 3).map(h => `${h.host} ${nPct(h.share)}`).join(' · ')}</div></td><td class="num">${fmt(c.n_genes)}</td><td class="num">${c.oil24 == null ? '–' : fmt1(c.oil24)}</td><td class="num">${c.oil25 == null ? '–' : fmt1(c.oil25)}</td>${cell(c.lfc24)}${cell(c.lfc25)}${cell(c.baci24)}${cell(c.baci25)}${cell(c.ice)}${cell(c.top)}</tr>`; });
  });
  box.innerHTML = html + '</tbody></table></div>';
  box.querySelectorAll('tr[data-c]').forEach(tr => { tr.onclick = () => { const s = document.getElementById('nsel'); if (s) { s.value = tr.dataset.c; s.onchange(); s.scrollIntoView({ behavior: 'smooth', block: 'center' }); } }; });
  sec.append(el('p', { class: 'note' }, 'Anammox (hydrazine synthase) has no genes in the coassembly, so it is not listed. Four ambiguous KO pairs are split by host. narGH in Nitrospira or Nitrospina counts as nitrite oxidoreductase (nxrAB), and as nitrate reductase in all other hosts. amoABC counts as ammonia monooxygenase only in nitrifier hosts. The same KOs in Cycloclasticus and a Porticoccaceae genome are listed as copper monooxygenases (CuMMO).'));
}

/* ---------- time-course explorer ---------- */
function nExplorer(sec, J) {
  sec.append(el('h3', {}, 'Time courses'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'nsel', 'aria-label': 'Gene group' });
  J.groups.forEach(g => { const cs = J.categories.filter(x => x.group === g); if (!cs.length) return; const og = document.createElement('optgroup'); og.label = g; cs.forEach(x => og.append(new Option(x.name, x.name))); sel.append(og); });
  c.append(el('label', { for: 'nsel' }, 'Gene group'), sel);
  const v = el('div'); sec.append(v);
  function draw() {
    v.innerHTML = ''; const cat = J.categories.find(x => x.name === sel.value);
    const box = fig(v, cat.name, `CPM of protein-coding reads; lines join per-day means of replicate libraries. Before dosing in 2025-26, seawater medians were ${cat.pre_exp25 ?? '–'} CPM in the diesel pool and ${cat.pre_ctl25 ?? '–'} CPM in the control.`);
    legendPools(box.parentNode.insertBefore(el('div'), box));
    timecourse(box, nSeries(J, J.timecourse[cat.name]), { ylab: 'CPM', unit: 'CPM', h: 200, floor: 0.05 });
  }
  sel.onchange = draw; sel.value = nCat(J, 'transporters').name; draw();
}

/* ---------- host-normalised N acquisition ---------- */
function nStress(sec, J) {
  sec.append(el('h3', {}, 'Nitrogen hunger inside each genome'));
  sec.append(el('p', { class: 'kicker' }, 'A bloom raises a genome’s share of every gene, so a rise in community CPM alone does not show a regulatory response. Here each genome’s nitrogen-acquisition transcripts are divided by all of its own transcripts. These are the ammonium, nitrate and urea transporters, the PII regulators, NtrBC and glutamine synthetase. A genome that leans right switches on nitrogen scavenging in the diesel pool beyond what its abundance explains. In 2025-26 nearly every major diesel responder leans right. Cycloclasticus roughly doubles its nitrogen-acquisition share, and several Sulfitobacter genomes rise 2–5×.'));
  const se = el('div', { class: 'controls' }); sec.append(se);
  const sel = el('select', { id: 'nsse', 'aria-label': 'Season' }); ['2025-26', '2024-25'].forEach(s => sel.append(new Option(s, s))); se.append(el('label', { for: 'nsse' }, 'Winter'), sel);
  const v = el('div'); sec.append(v);
  function draw() {
    v.innerHTML = '';
    const R = J.nstress.filter(x => x.season === sel.value).sort((a, b) => (b.host_oil + b.host_ctl) - (a.host_oil + a.host_ctl)).slice(0, 22).sort((a, b) => b.lfc - a.lfc);
    const box = fig(v, `Nitrogen-acquisition share of each genome’s transcripts, diesel vs control seawater (${sel.value})`, 'The 22 genomes with the most transcripts in either pool (at least 1,000 CPM). Dots show the median per-mille share in the control (blue) and diesel (orange) pools. The number is log₂ diesel ÷ control. Hover for genome size in each pool.');
    const W = width(box), lw = Math.min(300, W * 0.4), rh = 20, H = R.length * rh + 40;
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
    const mx = d3.max(R, d => Math.max(d.share_oil, d.share_ctl)) || 1;
    const x = d3.scaleLog().domain([Math.max(0.05, d3.min(R, d => Math.min(d.share_oil || 1, d.share_ctl || 1)) * 0.8), mx * 1.2]).range([lw + 10, W - 60]).clamp(true);
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${R.length * rh + 10})`).call(d3.axisBottom(x).ticks(5, '~g'));
    svg.append('text').attr('x', (lw + W) / 2).attr('y', H - 2).attr('text-anchor', 'middle').style('font-size', '11px').attr('fill', css('--ink-3')).text('‰ of the genome’s transcripts');
    R.forEach((d, i) => {
      const y = i * rh + 6;
      svg.append('text').attr('x', lw).attr('y', y + 4).attr('text-anchor', 'end').style('font-size', '12px').text(d.taxon === d.bin ? d.bin : `${d.taxon.length > 34 ? d.taxon.slice(0, 33) + '…' : d.taxon} ${d.bin}`);
      svg.append('line').attr('x1', x(d.share_ctl || 0.05)).attr('x2', x(d.share_oil || 0.05)).attr('y1', y).attr('y2', y).attr('stroke', css('--ink-3')).attr('stroke-width', 2);
      [['share_ctl', '--ctl'], ['share_oil', '--exp']].forEach(([k, c]) => svg.append('circle').attr('cx', x(d[k] || 0.05)).attr('cy', y).attr('r', 5).attr('fill', css(c)).attr('stroke', css('--panel')).attr('stroke-width', 2));
      svg.append('text').attr('x', W - 54).attr('y', y + 4).style('font-size', '12px').attr('fill', css('--ink-2')).text(nSgn(d.lfc));
      svg.append('rect').attr('x', 0).attr('y', y - rh / 2).attr('width', W).attr('height', rh).attr('fill', 'transparent')
        .on('mousemove', ev => showTip(ev, `<b>${d.taxon} ${d.bin}</b><br>N-acquisition share: ${d.share_ctl}‰ control → ${d.share_oil}‰ diesel (${nSgn(d.lfc)} log₂)<br>genome CPM: ${fmt(d.host_ctl)} control, ${fmt(d.host_oil)} diesel`)).on('mouseleave', hideTip);
    });
    box.parentNode.insertBefore(el('div', { class: 'legend' }, `<span><i style="background:var(--ctl)"></i>control pool</span><span><i style="background:var(--exp)"></i>diesel pool</span>`), box);
  }
  sel.onchange = draw; draw();
}

/* ---------- nitrifiers and CuMMO genomes ---------- */
function nGenomes(sec, J) {
  sec.append(el('h3', {}, 'Ammonia oxidisers lose ground; copper monooxygenases of hydrocarbon degraders take over'));
  sec.append(el('p', { class: 'kicker' }, 'Ammonia monooxygenase (amoABC) and methane monooxygenase (pmoABC) share KOs, because both belong to the copper membrane monooxygenase (CuMMO) family. In these pools the family splits cleanly by host. In 2024-25 the nitrifier copies belong to a Nitrosopumilus archaeon and a Nitrosomonadaceae bacterium, which together make about 200 TPM in the control pool after dosing. In the diesel pool they make about 60 TPM. That gap is larger than before dosing (190 vs 360 TPM), and the Nitrosomonadaceae falls most. Their per-genome amo share stays the same, so the nitrifiers lose abundance; the cells that remain still express amo as before. In 2025-26 nitrifiers are almost absent from both pools from the start. The CuMMO copies that rise with diesel belong to Cycloclasticus and to an uncultured Porticoccaceae genus (50-400-T64), which grows to about 3.6% of all transcripts in the 2025-26 diesel pool (BACI q = 0.02). Symbiotic Cycloclasticus use related enzymes to oxidise propane, ethane and butane, and the Deepwater Horizon Cycloclasticus carried the same genes (Rubin-Blum et al. 2017, <a href="https://doi.org/10.1038/nmicrobiol.2017.93">Nat Microbiol</a>). These are most likely short-chain alkane monooxygenases, not nitrifiers or methanotrophs.'));
  const box = fig(sec, 'Genome abundance and CuMMO expression', 'TPM = MAG share of community transcripts (seawater medians). MO ‰ = amo/CuMMO transcripts as a share of the genome’s own transcripts. BACI = MAG transcript-share change from the pool’s own pre-dosing baseline, diesel minus control (CLR, log₂). Dashes mean the genome is too rare to estimate.');
  box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Genome</th><th>Winter</th><th class="num">TPM pre<br>diesel / ctl</th><th class="num">TPM after<br>diesel / ctl</th><th class="num">MO ‰ after<br>diesel / ctl</th><th class="num">BACI</th></tr></thead><tbody>${
    J.genomes.flatMap(g => g.summary.map((s, k) => `<tr><td>${k ? '' : `<b>${g.label.split(' (')[0]}</b><br><span style="font-size:11.5px;color:var(--ink-3)">${g.label.includes(' (') ? g.label.split(' (')[1].replace(')', '') + ' · ' : ''}<span class="mono">${g.bin}</span></span>`}</td><td>${s.season}</td><td class="num">${fmt(s.pre_exp)} / ${fmt(s.pre_ctl)}</td><td class="num">${fmt(s.oil)} / ${fmt(s.ctl)}</td><td class="num">${s.oil_mo == null ? '–' : fmt(s.oil_mo)} / ${s.ctl_mo == null ? '–' : fmt(s.ctl_mo)}</td><td class="num">${nSgn(s.baci)}${s.q != null && s.q < 0.05 ? ' <span class="pill up">q ' + d3.format('.2g')(s.q) + '</span>' : ''}</td></tr>`)).join('')}</tbody></table></div>`;
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'ngen', 'aria-label': 'Genome' }); J.genomes.forEach(g => sel.append(new Option(`${g.label} (${g.bin})`, g.bin))); c.append(el('label', { for: 'ngen' }, 'Genome'), sel);
  const v = el('div'); sec.append(v);
  function draw() { v.innerHTML = ''; const g = J.genomes.find(x => x.bin === sel.value);
    const b = fig(v, `${g.label} ${g.bin}: transcripts per million`, 'MAG share of the community transcriptome (TPM).');
    legendPools(b.parentNode.insertBefore(el('div'), b)); timecourse(b, nSeries(J, g.tpm), { ylab: 'TPM', unit: 'TPM', h: 190, floor: 0.5 }); }
  sel.onchange = draw; draw();
}

/* ---------- reading ---------- */
function nReading(sec, J) {
  const C = f => nCat(J, f);
  const nrt = C('transporters'), nas = C('nasA'), nirbd = C('nirBD'), glnk = C('glnB'), glna = C('glnA'), urt = C('urtABCDE'), cyn = C('cynS'), nar = C('narGHI'), nor = C('norBC'), nos = C('nosZ'), hcp = C('hcp'), nif = C('nifHDK'), nrf = C('nrfAH');
  sec.append(el('h3', {}, 'How the pieces fit, and what is still a hypothesis'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li><b>Observed: nitrogen scavenging turns on with diesel.</b> In 2025-26 seawater the diesel pool expresses far more of every nitrogen-acquisition gene than the control. Nitrate/nitrite transporters are ${fmt(nrt.oil25)} vs ${fmt(nrt.ctl25)} CPM, assimilatory nitrate reductase ${fmt(nas.oil25)} vs ${fmt(nas.ctl25)}, and NADH nitrite reductase (nirBD, the assimilatory step) ${fmt(nirbd.oil25)} vs ${fmt(nirbd.ctl25)}. PII regulators are ${nSgn(glnk.lfc25)} log₂, glutamine synthetase ${nSgn(glna.lfc25)}, the urea transporter ${nSgn(urt.lfc25)} and cyanase ${nSgn(cyn.lfc25)}. All rise against each pool’s own baseline (BACI ${nSgn(nrt.baci25)} for the transporters). This is the classic Ntr starvation programme. It is carried by Cycloclasticus, Sulfurimonas and Sulfitobacter, and it rises within those genomes, not just with their abundance.</li>
    <li><b>Observed: in 2024-25 it is mostly pool history.</b> The diesel pool’s single pre-dosing sample already held ${fmt(nrt.pre_exp24)} CPM of nitrate/nitrite transporters, against ${fmt(nrt.pre_ctl24)} in the control’s. So the 2024-25 contrast (${nSgn(nrt.lfc24)} log₂) is not larger than the pre-existing difference, and BACI is negative (${nSgn(nrt.baci24)}). The same holds for the sulfur and hydrogen genes. That winter the pools were filled about two weeks apart after the intake lines froze, so their communities had already diverged.</li>
    <li><b>Observed: nitrate is reduced, but denitrification is not completed.</b> Membrane nitrate reductase rises in 2025-26 (narGHI ${nSgn(nar.lfc25)} log₂, BACI ${nSgn(nar.baci25)}), mostly because the Cycloclasticus genomes that carry it bloomed. Within the minor genome metabat_5432 its share of transcripts rises after dosing in both pools (about 1.5‰ before; 7.8‰ diesel, 6.8‰ control), so it is not a diesel switch; the dominant lorbin_184 barely expresses it. The nitrite/NO-stress marker hcp also rises (${nSgn(hcp.lfc25)}), on unbinned Cycloclasticus contigs. Nitric oxide reductase falls (norBC ${nSgn(nor.lfc25)}), nitrous oxide reductase stays at trace levels (nosZ ${fmt1(nos.oil25)} vs ${fmt1(nos.ctl25)} CPM), and nirK and nirS do not respond. So nitrate is reduced to nitrite mainly for assimilation; Cycloclasticus metabat_9342 keeps its assimilatory nitrite reductase (nirBD) at 2–3.5‰ of its transcripts before and after dosing. No step that would remove nitrogen as N₂O or N₂ rises above trace levels.</li>
    <li><b>Observed: nobody fixes N₂ to relieve the shortage.</b> Nitrogenase stays below 10 CPM in both pools (mostly unbinned Arcobacter-like contigs) and does not follow diesel. Anammox genes are absent. DNRA (nrfAH) is a trace activity (${fmt1(nrf.oil25)} CPM).</li>
    <li><b>Observed: nitrifiers decline in the diesel pool.</b> In 2024-25 the ammonia-oxidiser genomes make about 60 TPM in the diesel pool against 200 in the control, and ammonia monooxygenase falls with them. Hydrocarbons inhibit ammonia monooxygenase, and competition for ammonium with a nitrogen-hungry heterotroph bloom would act in the same direction. With one pool per treatment and only one pre-dosing sample, this is suggestive, not significant (MAG BACI −0.6 and −1.6, q = 1).</li>
    <li><b>Hypothesis: ice is less nitrogen-limited.</b> Nitrate transporters are higher in top ice than in bottom ice in none of 13 pairs, and urea transporters in only 1 of 13. Brine pockets concentrate nutrients, which may ease nitrogen limitation near the surface. Composition and bloom state also differ between layers, so this needs nutrient measurements in brine.</li>
    <li><b>What would test this.</b> Ammonium, nitrate, nitrite and urea measured in the filtrate already collected, plus N₂O in the headspace. Particulate C:N in the diesel pool would show whether the nitrogen hunger is real.</li>`));
}

RENDER.nitrogen = nitrogenSec;
