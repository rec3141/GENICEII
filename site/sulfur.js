/* Sulfur & hydrogen section: the chemolithotrophic side economy in the diesel pools. Data: sulfur.json (sulfur/scripts/sulfur_site.py).
   Uses helpers from app.js (el, fig, width, css, showTip, hideTip, fmt, fmt1, tickf, timecourse, legendPools, POOL). */
const sfSgn = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.1f')(v);
const sfPct = d3.format('.0%');

function sfSeries(J, values) {
  // average replicate libraries per pool x medium x day (CTL/EXP seawater and ice only)
  const out = {};
  J.samples.forEach((s, i) => { if (!['CTL', 'EXP'].includes(s.pool) || !['SW', 'SI'].includes(s.medium) || s.day == null || values[i] == null) return;
    const k = `${s.season}|${s.pool}|${s.medium}`; ((out[k] ||= {})[s.day] ||= []).push([values[i], s.name]); });
  return Object.entries(out).map(([k, byDay]) => { const [season, pool, medium] = k.split('|');
    return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a, x => x[0]), name: a.map(x => x[1]).join(', ') })) }; });
}
const sfCat = (J, frag) => J.categories.find(c => c.name.includes(frag));

function sulfurSec(sec) {
  const J = D.sulfur;
  if (!J) { sec.innerHTML = '<h2>Sulfur & hydrogen</h2><p class="kicker">sulfur.json is missing.</p>'; return; }
  const hyd = sfCat(J, 'quinone-reactive'), por = sfCat(J, 'pyruvate:ferredoxin'), acl = sfCat(J, 'ATP-citrate'), hox = sfCat(J, 'hoxFUYH'), dsr = J.direction.filter(d => d.gene === 'dsrAB');
  sec.innerHTML = `<h2>Sulfur & hydrogen</h2>
  <p class="kicker">Next to the hydrocarbon degraders, the diesel pools ran a chemolithotrophic side economy. Bacteria fixed CO₂ by burning hydrogen and reduced sulfur with oxygen. Its centre is Sulfurimonas (Campylobacterota), whose quinone-reactive NiFe hydrogenase, sulfide dehydrogenase and reverse-TCA cycle are among the most strongly diesel-associated transcripts in the dataset. This page follows 44 sulfur, hydrogen, carbon-fixation and electron-acceptor gene groups through both winters. Expression is CPM of protein-coding reads, with gene groups from the coordinate-corrected KOs, and hosts from MAGs or, for unbinned contigs, Kaiju. One pool per treatment, and the two pools carried different Sulfurimonas strains (see below), so all contrasts are descriptive.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  [[`${sfSgn(hyd.lfc24)} / ${sfSgn(hyd.lfc25)}`, 'log₂ diesel ÷ control, Sulfurimonas-type quinone-reactive hydrogenase (2024-25 / 2025-26 seawater)'],
   [`${sfPct(acl.hosts.diesel[0].share)}`, `of ATP-citrate lyase (reverse TCA) transcripts in the diesel pool are from Sulfurimonas (${sfPct(por.hosts.diesel[0].share)} of PFOR)`],
   [`${sfSgn(hox.lfc25)}`, `log₂ diesel ÷ control, bidirectional Hox hydrogenase in 2025-26, ${sfPct(hox.hosts.diesel[0].share)} from Cycloclasticus`],
   [`≤ ${sfPct(d3.max(dsr, d => d.reductive))}`, 'of dsrAB transcripts come from sulfate reducers; the rest is reverse (oxidative) Dsr']]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  sfTable(sec, J);
  sfExplorer(sec, J);
  sfHosts(sec, J);
  sfDirection(sec, J);
  sfSulfurimonas(sec, J);
  sfReading(sec, J);
}

/* ---------- overview table: every gene group, coloured log2 ratios, top diesel hosts ---------- */
function sfTable(sec, J) {
  sec.append(el('h3', {}, 'The whole inventory at a glance'));
  const box = fig(sec, '', 'Median CPM in seawater after dosing. Ratios are log₂ (diesel ÷ control) of medians. BACI = change from each pool’s own five pre-dosing samples, diesel minus control (2025-26 seawater only). Ice = paired ice minus water (24 pairs, both winters). Cell colour runs from blue (higher in control or water) to orange (higher in diesel or ice). Hosts = genomes carrying the most diesel-pool transcripts. Click a row to see its time course.');
  const col = v => { if (v == null) return 'transparent'; const t = Math.max(-1, Math.min(1, v / 4)); return t >= 0 ? d3.interpolateRgb(css('--panel'), css('--exp'))(t * 0.85) : d3.interpolateRgb(css('--panel'), css('--ctl'))(-t * 0.85); };
  const cell = v => `<td class="num" style="background:${col(v)}">${sfSgn(v)}</td>`;
  let html = `<div class="tablewrap" style="max-height:none"><table><thead><tr><th>Gene group · main hosts in the diesel pool</th><th class="num">Genes</th><th class="num">CPM diesel 24-25</th><th class="num">CPM diesel 25-26</th><th class="num">24-25</th><th class="num">25-26</th><th class="num">BACI 25-26</th><th class="num">Ice</th></tr></thead><tbody>`;
  J.groups.forEach(g => {
    html += `<tr><td colspan="8" style="background:var(--grid);font-weight:650">${g}</td></tr>`;
    J.categories.filter(c => c.group === g).forEach(c => {
      html += `<tr data-c="${c.name}"><td style="white-space:normal;min-width:300px">${c.name}<div style="font-size:11.5px;color:var(--ink-3);margin-top:2px">${
        c.hosts.diesel.slice(0, 3).map(h => `${h.host} ${sfPct(h.share)}`).join(' · ')}</div></td><td class="num">${fmt(c.n_genes)}</td><td class="num">${c.oil24 == null ? '–' : fmt1(c.oil24)}</td><td class="num">${c.oil25 == null ? '–' : fmt1(c.oil25)}</td>${cell(c.lfc24)}${cell(c.lfc25)}${cell(c.baci25)}${cell(c.ice)}</tr>`; });
  });
  box.innerHTML = html + '</tbody></table></div>';
  box.querySelectorAll('tr[data-c]').forEach(tr => { tr.onclick = () => { const s = document.getElementById('sfsel'); if (s) { s.value = tr.dataset.c; s.onchange(); s.scrollIntoView({ behavior: 'smooth', block: 'center' }); } }; });
}

/* ---------- time-course explorer ---------- */
function sfExplorer(sec, J) {
  sec.append(el('h3', {}, 'Time courses'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'sfsel', 'aria-label': 'Gene group' });
  J.groups.forEach(g => { const og = document.createElement('optgroup'); og.label = g; J.categories.filter(x => x.group === g).forEach(x => og.append(new Option(x.name, x.name))); sel.append(og); });
  c.append(el('label', { for: 'sfsel' }, 'Gene group'), sel);
  const v = el('div'); sec.append(v);
  function draw() {
    v.innerHTML = ''; const cat = J.categories.find(x => x.name === sel.value);
    const box = fig(v, cat.name, `CPM of protein-coding reads; lines join per-day means of replicate libraries. Before dosing in 2025-26, seawater medians were ${cat.pre_exp25 ?? '–'} CPM in the diesel pool and ${cat.pre_ctl25 ?? '–'} CPM in the control.`);
    legendPools(box.parentNode.insertBefore(el('div'), box));
    timecourse(box, sfSeries(J, J.timecourse[cat.name]), { ylab: 'CPM', unit: 'CPM', h: 200, floor: 0.05 });
  }
  sel.onchange = draw; sel.value = sfCat(J, 'quinone-reactive').name; draw();
}

/* ---------- who produces and who consumes hydrogen / oxidises sulfur: host shares ---------- */
function sfHosts(sec, J) {
  sec.append(el('h3', {}, 'Who does what'));
  sec.append(el('p', { class: 'kicker' }, 'Share of each gene group’s diesel-pool transcripts (both winters, after dosing) carried by its main genomes. The hydrogen-oxidising, CO₂-fixing core belongs to two Campylobacterota: Sulfurimonas MAGs and an unbinned Sulfurovum-like population. Sulfide and thiosulfate oxidation is shared with Sulfitobacter, Cycloclasticus and UBA9145. Cycloclasticus also carries the H₂-cycling Hox hydrogenase, a membrane nitrate reductase (narGH) and the methanethiol-detoxifying mddA.'));
  const pick = ['quinone-reactive', 'respiratory uptake', 'hoxFUYH', 'H2-evolving', '[FeFe]', 'nitrogenase', 'ATP-citrate', 'pyruvate:ferredoxin', 'fccAB', 'sqr', 'soxB', 'soxCD', 'tsdA', 'sdo/ETHE1', 'mddA', 'dmdA', 'narGH', 'napAB', 'nosZ', 'cbb3'];
  const C = pick.map(p => sfCat(J, p)).filter(Boolean);
  const box = fig(sec, 'Main genomes behind each gene group (diesel pool)', 'Bars show the five largest contributors; the grey remainder is everything else. Hover a segment for the genome.');
  const W = width(box), lw = Math.min(330, W * 0.44), rh = 22, H = C.length * rh + 30;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const x = d3.scaleLinear().domain([0, 1]).range([lw + 10, W - 10]);
  const fam = h => /Sulfurimonas|Sulfurovum/.test(h) ? '--s3' : /Cycloclasticus/.test(h) ? '--exp' : /Sulfitobacter|Yoonia|Roseobacter/.test(h) ? '--s7' : /Desulf/.test(h) ? '--s8' : /Colwellia|Cognaticolwellia|Glaciecola/.test(h) ? '--ctl' : '--s4';
  C.forEach((c, i) => {
    const y = i * rh + 8; let x0 = 0;
    svg.append('text').attr('x', lw).attr('y', y + 11).attr('text-anchor', 'end').text(c.name.length > 52 ? c.name.slice(0, 51) + '…' : c.name);
    c.hosts.diesel.forEach(h => { const w = h.share;
      svg.append('rect').attr('x', x(x0)).attr('y', y).attr('width', Math.max(0, x(x0 + w) - x(x0) - 2)).attr('height', 14).attr('rx', 2)
        .attr('fill', /^unbinned, unclassified/.test(h.host) ? css('--grid') : css(fam(h.host))).attr('opacity', /^unbinned/.test(h.host) ? 0.55 : 1)
        .on('mousemove', ev => showTip(ev, `<b>${c.name}</b><br>${h.host}: ${sfPct(h.share)} of diesel-pool transcripts`)).on('mouseleave', hideTip);
      x0 += w; });
    svg.append('rect').attr('x', x(x0)).attr('y', y).attr('width', Math.max(0, x(1) - x(x0))).attr('height', 14).attr('fill', css('--grid')).attr('rx', 2);
  });
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${C.length * rh + 8})`).call(d3.axisBottom(x).ticks(5).tickFormat(sfPct));
  box.parentNode.insertBefore(el('div', { class: 'legend' }, [['--s3', 'Sulfurimonas / Sulfurovum'], ['--exp', 'Cycloclasticus'], ['--s7', 'Roseobacters (Sulfitobacter, Yoonia)'], ['--ctl', 'Colwellia relatives'], ['--s8', 'sulfate reducers'], ['--s4', 'other genomes']]
    .map(([c, l]) => `<span><i style="background:var(${c});height:10px"></i>${l}</span>`).join('') + '<span><i style="background:var(--s4);opacity:.55;height:10px"></i>unbinned (Kaiju genus), faded</span>'), box);
}

/* ---------- dissimilatory sulfur: reductive vs oxidative ---------- */
function sfDirection(sec, J) {
  sec.append(el('h3', {}, 'Is anyone making sulfide?'));
  sec.append(el('p', { class: 'kicker' }, 'Dissimilatory sulfite reductase (dsrAB), APS reductase (aprAB) and the Qmo complex run either way. Sulfate reducers use them to make sulfide, and many sulfur oxidisers use the same enzymes in reverse. The direction is assigned here by host phylum (Desulfobacterota and relatives = reductive). Sulfate reduction is a trace activity in both pools. Nearly all dsrAB transcripts come from gammaproteobacterial sulfur oxidisers (Thiohalobacter-like and Thermochromatium-like contigs, SUP05/Thioglobus, UBA9214), and most aprAB from SAR11 (Pelagibacter), which carries apr for sulfur metabolism, not respiration. The one sulfate-reducer MAG (Desulforhopalus comebin_3581_sub) contributes a few percent. Whatever sulfide or H₂ feeds Sulfurimonas is not coming from detectable sulfate reduction in the water or ice samples.'));
  const box = fig(sec, '', 'Mean CPM after dosing, and the share by host direction. Unknown = unbinned contigs without a Kaiju phylum.');
  box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Genes</th><th>Season</th><th>Pool</th><th class="num">CPM</th><th class="num">Reductive (sulfate reducers)</th><th class="num">Oxidative / other</th><th class="num">Unknown</th></tr></thead><tbody>${
    J.direction.map(d => `<tr><td>${d.gene}</td><td>${d.season}</td><td>${d.pool}</td><td class="num">${fmt1(d.cpm)}</td><td class="num">${sfPct(d.reductive)}</td><td class="num">${sfPct(d.oxidative)}</td><td class="num">${sfPct(d.unknown)}</td></tr>`).join('')}</tbody></table></div>`;
}

/* ---------- Sulfurimonas genomes ---------- */
function sfSulfurimonas(sec, J) {
  const SF = J.sulfurimonas;
  sec.append(el('h3', {}, 'Sulfurimonas: three genomes, two pools, more in ice'));
  sec.append(el('p', { class: 'kicker' }, 'Three Sulfurimonas MAGs carry the response. semibin_088 (1.7 Mbp, 92% complete) is the dominant one. It has sqr and fccAB for sulfide, both NiFe uptake hydrogenases, the full reverse TCA cycle, and cbb3 plus aa3 oxidases, but no sox genes and no nitrate reductase in the bin. So in these pools it most likely oxidises H₂ and sulfide with oxygen at low concentration, not nitrate. lorbin_167 (3.1 Mbp) is the denitrifying generalist, with sox, nap, nir, nor and nos. comebin_1140_sub resembles semibin_088. All three are more abundant in ice than in water, in DNA as well as RNA. The control and diesel pools carry different semibin_088 strains (pool FST 0.35, against 0.03 between ice and water).'));
  const names = SF.inventory.map(i => i.bin);
  const cats = [...new Set(SF.inventory.flatMap(i => Object.keys(i.genes)))].filter(c => !/cysDN|maturation/.test(c));
  const order = J.categories.map(c => c.name).filter(c => cats.includes(c));
  const box = fig(sec, 'Gene inventory of the three Sulfurimonas MAGs', 'Number of genes in each group per MAG (assimilatory sulfate and hydrogenase maturation genes omitted). An empty cell means no gene in the bin; for semibin_088 at 92% completeness a missing pathway is probably absent.');
  box.innerHTML = `<div class="tablewrap" style="max-height:none"><table><thead><tr><th>Gene group</th>${SF.inventory.map(i => `<th class="num">${i.bin}<br><span style="font-weight:400">${d3.format('.1f')(i.size / 1e6)} Mbp, ${i.compl}%</span></th>`).join('')}</tr></thead><tbody>${
    order.map(c => `<tr><td>${c}</td>${SF.inventory.map(i => `<td class="num">${i.genes[c] ? `<span class="pill up">${i.genes[c]}</span>` : ''}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`;
  // stats
  const b2 = fig(sec, 'Ice, pools and the diesel effect', 'RNA = paired ice-vs-water shift in transcript share (CLR, log₂, 24 pairs). DNA = paired ice-vs-water shift in share of mapped Illumina reads (19 pairs). Pool FST = control vs diesel pool, 2024-25 ice. BACI = MAG transcript-share change from each pool’s own pre-dosing baseline, diesel minus control (CLR, log₂).');
  b2.innerHTML = `<div class="tablewrap"><table><thead><tr><th>MAG</th><th class="num">RNA ice − water</th><th class="num">DNA ice − water</th><th class="num">top − bottom (RNA)</th><th class="num">FST pools</th><th class="num">FST ice/water</th><th class="num">BACI 24-25</th><th class="num">BACI 25-26</th></tr></thead><tbody>${
    SF.stats.map(s => `<tr><td class="mono">${s.bin}</td><td class="num">${sfSgn(s.rna_ice)} (${s.rna_ice_up}/${s.rna_ice_n})</td><td class="num">${sfSgn(s.dna_ice)}${s.dna_ice_n ? ` (${s.dna_ice_up}/${s.dna_ice_n})` : ''}</td><td class="num">${sfSgn(s.rna_top)}</td><td class="num">${s.fst_pool_ice == null ? '–' : d3.format('.2f')(s.fst_pool_ice)}</td><td class="num">${s.fst_ice_water == null ? '–' : d3.format('.2f')(s.fst_ice_water)}</td><td class="num">${sfSgn(s.baci24)}</td><td class="num">${sfSgn(s.baci25)}${s.baci25_q != null && s.baci25_q < 0.05 ? ' <span class="pill up">q ' + d3.format('.2g')(s.baci25_q) + '</span>' : ''}</td></tr>`).join('')}</tbody></table></div>`;
  // RNA time course, selectable MAG
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'sfmag', 'aria-label': 'Sulfurimonas MAG' }); names.forEach(n => sel.append(new Option(n, n))); c.append(el('label', { for: 'sfmag' }, 'MAG'), sel);
  const v = el('div'); sec.append(v);
  function draw() { v.innerHTML = ''; const box = fig(v, `${sel.value}: transcripts per million`, 'Share of the community transcriptome (TPM, MAG level). In 2025-26 both pools already carried Sulfurimonas before dosing; the diesel pool held it while the control declined.');
    legendPools(box.parentNode.insertBefore(el('div'), box)); timecourse(box, sfSeries(J, SF.tpm[sel.value]), { ylab: 'TPM', unit: 'TPM', h: 190, floor: 1 }); }
  sel.onchange = draw; draw();
}

/* ---------- reading of the evidence ---------- */
function sfReading(sec, J) {
  const cbb = sfCat(J, 'cbb3'), cyd = sfCat(J, 'bd oxidase'), nar = sfCat(J, 'narGH'), mdd = sfCat(J, 'mddA'), q = sfCat(J, 'quinone-reactive');
  sec.append(el('h3', {}, 'How the pieces fit, and what is still a hypothesis'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li><b>Observed.</b> Sulfurimonas hydrogen oxidation and CO₂ fixation are far higher in the diesel pool in both winters (quinone-reactive hydrogenase ${sfSgn(q.lfc24)} and ${sfSgn(q.lfc25)} log₂). In 2025-26 this holds even against each pool’s own pre-dosing level (BACI ${sfSgn(q.baci25)} for the hydrogenase, ${sfSgn(sfCat(J, 'ATP-citrate').baci25)} for ATP-citrate lyase), so it is more than pool history. In 2024-25 the diesel pool’s single pre-dosing sample already held 903 CPM of the hydrogenase, against 15–69 CPM in the control’s, so in that winter the contrast is largely pool history (MAG-level BACI −0.3 to −1.3). The 2024-25 pools were filled about two weeks apart after the intake lines froze.</li>
    <li><b>Observed, then revised.</b> High-affinity terminal oxidases rise in the 2025-26 diesel pool (cbb3 ${sfSgn(cbb.lfc25)}, bd ${sfSgn(cyd.lfc25)} log₂), together with Cycloclasticus’s membrane nitrate reductase narGH (${sfSgn(nar.lfc25)}, BACI ${sfSgn(nar.baci25)}). At first sight this looks like falling oxygen. Measured within each genome, though, it is not: the rise comes from the Cycloclasticus bloom, whose genomes carry these genes, and Cycloclasticus does not change its own oxidase mix (see <a href="#oxygen">Oxygen</a>). The pools may or may not have lost oxygen; the transcripts cannot tell, and there are no dissolved-oxygen measurements.</li>
    <li><b>Observed.</b> Methanethiol S-methyltransferase mddA is ${sfSgn(mdd.lfc24)} and ${sfSgn(mdd.lfc25)} log₂ higher in diesel seawater, mostly from Cycloclasticus. Methanethiol and other volatile sulfur compounds may be released as diesel or DMSP is degraded.</li>
    <li><b>Hypothesis: where the H₂ and sulfide come from.</b> No sulfate-reducer or strong H₂-evolving signal is visible in the water and ice. Candidates: (1) anoxic micro-niches in oil droplets, oil snow or brine channels that the samples under-represent; (2) H₂ from fermenters and from Cycloclasticus’s bidirectional Hox hydrogenase, which can release H₂ when NADH accumulates; (3) reduced sulfur already in the diesel (organosulfur) and in DMSP breakdown. Dissolved H₂, sulfide and O₂ profiles in the pools would separate these.</li>
    <li><b>Hypothesis: why ice.</b> Sulfurimonas is consistently enriched in ice. Brine pockets are low in oxygen and concentrate solutes, which fits a microaerophilic H₂/sulfide oxidiser.</li>`));
}

RENDER.sulfur = sulfurSec;
