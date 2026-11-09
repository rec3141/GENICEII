/* Viruses section: viral DNA share (nanopore, Illumina), viral transcription and diversity, phage-host pairs, bloom vOTUs, prophages,
   ice vs water, and the virus-concentration method test. Data: viruses.json (viruses/scripts/vsite.py). Helpers from app.js. */
const vrSgn = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.1f')(v);
function vrSeries(rows, get, filt = () => true) {
  // rows: [{season, pool, medium, day, name}] -> per season/pool/medium series of per-day means
  const out = {};
  rows.filter(r => r.day != null && filt(r)).forEach(r => { const v = get(r); if (v == null) return;
    const k = `${r.season}|${r.pool}|${r.medium}`; ((out[k] ||= {})[r.day] ||= []).push([v, r.name]); });
  return Object.entries(out).map(([k, byDay]) => { const [season, pool, medium] = k.split('|');
    return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a, x => x[0]), name: a.map(x => x[1]).join(', ') })) }; });
}
function vrChips(parent, opts, init, onChange) {
  const c = el('div', { class: 'controls' }); parent.append(c); let cur = init;
  opts.forEach(([k, l]) => { const b = el('button', { class: 'chip', 'aria-pressed': k === cur }, l);
    b.onclick = () => { cur = k; c.querySelectorAll('button').forEach(x => x.setAttribute('aria-pressed', x === b)); onChange(k); }; c.append(b); });
  return c;
}

function virusesSec(sec) {
  const J = D.viruses;
  if (!J) { sec.innerHTML = '<h2>Viruses</h2><p class="kicker">viruses.json is missing.</p>'; return; }
  const st = J.stats;
  sec.innerHTML = `<h2>Viruses</h2>
  <p class="kicker">The phage-shock and other envelope-stress responses on the Metals &amp; toxicity page raised the question of a viral bloom after dosing. This page follows viruses in the pool samples through DNA (nanopore and Illumina metagenomes), transcription (the 130 RNA libraries) and diversity. The viral set is geNomad’s confident calls on the coassembly: virus score ≥ 0.95, at least two viral hallmark genes, and more viral than host genes by CheckV. Of ${fmt(st.genomad_all)} geNomad virus contigs, ${fmt(st.confident)} pass, collapsing to ${fmt(st.votus)} vOTUs (≥ 95% ANI over ≥ 85% of the shorter sequence). Most are partial genomes; CheckV rates ${fmt(st.quality['High-quality'] || 0)} high-quality and ${fmt(st.quality['Complete'] || 0)} complete. One pool per treatment, and the pools hold different strains (Ice vs water page), so diesel-vs-control contrasts are descriptive.</p>`;
  const rna = J.rna, sw = r => r.medium === 'SW';
  const med = (s, f) => d3.median(rna.filter(r => r.season === s && f(r)), r => r.viral);
  const preE = med('2025-26', r => sw(r) && r.pool === 'EXP' && r.phase === 'pre'), oilE = med('2025-26', r => sw(r) && r.pool === 'EXP' && r.phase === 'oiled');
  const preC = med('2025-26', r => sw(r) && r.pool === 'CTL' && r.phase === 'pre'), ctlC = med('2025-26', r => sw(r) && r.pool === 'CTL' && r.phase === 'control');
  const richE = d3.median(rna.filter(r => r.season === '2025-26' && sw(r) && r.pool === 'EXP' && r.phase === 'oiled'), r => r.rich);
  const richC = d3.median(rna.filter(r => r.season === '2025-26' && sw(r) && r.pool === 'CTL' && r.phase === 'control'), r => r.rich);
  const sr = el('div', { class: 'stat-row' }); sec.append(sr);
  [[fmt(st.votus), `vOTUs, ${fmt(st.active_votus)} of them transcribed (≥ 1 CPM in at least one library)`],
   [`${d3.format('.0f')(preE / 1000)}k → ${d3.format('.0f')(oilE / 1000)}k`, `viral transcripts per million CDS reads, 2025-26 diesel-pool seawater, before → after dosing (control ${d3.format('.0f')(preC / 1000)}k → ${d3.format('.0f')(ctlC / 1000)}k)`],
   [`${fmt(richE)} vs ${fmt(richC)}`, 'active vOTUs per library after dosing, 2025-26 seawater, diesel vs control pool'],
   ['~40×', 'rise of the main Cycloclasticus phage’s transcripts as its host bloomed (2025-26)']]
    .forEach(([v, l]) => sr.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  sec.append(el('h3', {}, 'Was there a viral bloom?'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li><b>Not across the community.</b> In both winters viral DNA and viral transcription went down in the diesel pool after dosing, not up; in Illumina DNA the viral share of reads in late seawater was 8.7% in the diesel pool vs 16.7% in the control (2024-25, days 25–70) and 11.7% vs 25% (2025-26, days 12–15). In 2025-26 the pools started level (about 5% of protein-coding transcripts were viral), then the diesel pool fell about fourfold while the control fell less. The number of active vOTUs and their diversity fell in the diesel pool too: fewer viruses, more dominated by a few.</li>
    <li><b>Yes for the phages of the oil-bloom bacteria.</b> The Cycloclasticus bloom carried its own phages. A 25.6 kb phage that tracks Cycloclasticus lorbin_722_sub (ρ = 0.9 across all libraries) rose from about 25 to over 1,000 transcripts per million by days 17–20, then fell to 83 by day 82 while its host stayed abundant. Its DNA rose too, in the diesel pool only (0–8 reads per library before dosing, 574 at day 11). In 2024-25, a phage of GCA-2708415 (comebin_1134), a diesel-pool resident before dosing, peaked at about 4,900 per million on day 8, before its host peaked on day 14. Its DNA was already abundant before dosing and did not rise, so that was a burst of activity, not growth.</li>
    <li><b>Giant viruses of algae collapsed in the diesel pool.</b> Megaviricetes (Phycodnaviridae and Mimiviridae-like) transcripts fell from about 5,800 to 250 per million by day 17 in the 2025-26 diesel pool, against 4,700 at day 19 in the control. That matches the decline of Phaeocystis and other algae on the Community response page. The one real viral bloom in the data is in the control pool: a virophage (Maveriviricetes, 22.8 kb, CheckV high quality), which replicates only inside cells co-infected by a giant virus, went from near zero to 7,400 reads per 4 M in DNA on day 3 and 16,000 by days 12–14, and peaked in transcripts on days 3–5. It never appears in the diesel pool, which suggests giant-virus infections of protists kept running in the control and stopped in the diesel pool.</li>
    <li><b>The phage-shock signal is not phage-driven.</b> The UBA9145 MAG behind the psp response carries no prophage and no filamentous-phage genes, and Inoviridae (Faserviricetes) transcripts are near zero in the diesel pool. Envelope stress from hydrocarbons remains the simpler explanation.</li>
    <li><b>Ice is not a viral hotspot.</b> Tailed-phage transcripts are slightly lower in ice than in the water beneath it (higher in ice in only 6 of 23 same-day pairs) and lower in the cold, salty top ice than at the bottom (2 of 13 pairs). Giant viruses and virophages show no consistent difference.</li>
    <li><b>No convincing prophage induction.</b> The provirus regions whose transcripts rise in the diesel pool express host-like genes at their edges (sugar and EPS genes, partitioning and transposition genes), not phage structural genes.</li>`));
  vrDNA(sec, J); vrActivity(sec, J); vrDiversity(sec, J); vrPairs(sec, J); vrHosts(sec, J); vrGroups(sec, J); vrProphage(sec, J); vrIce(sec, J); vrConc(sec, J);
}

function vrDNA(sec, J) {
  sec.append(el('h3', {}, 'Abundance: viral share of DNA'));
  const box = fig(sec, 'Viral contigs’ share of sequenced bases', `Nanopore metagenomes (filled markers, ${J.dna_np.length} pool libraries, to about day 11) and${J.dna_il ? ` Illumina metagenomes (hollow diamonds, ${J.dna_il.length} pool libraries, first 2 M read pairs mapped to the vOTUs at ≥ 95% identity)` : ' Illumina metagenomes (not yet mapped)'}. Filters retain cells and particles above about 0.2 µm, so this is mostly cell-associated and large viruses.`);
  legendPools(box.parentNode.insertBefore(el('div'), box));
  const ser = vrSeries(J.dna_np, r => r.viral, r => r.medium !== 'other');
  if (J.dna_il) ser.push(...vrSeries(J.dna_il, r => r.viral, r => r.medium !== 'other').map(s => ({ ...s, points_only: true })));
  timecourse(box, ser, { ylab: '% of bases', unit: '%', h: 210, log: false, fmt: d3.format('.1f') });
}

function vrActivity(sec, J) {
  sec.append(el('h3', {}, 'Activity: viral transcription'));
  const OPTS = [['all', 'All viral'], ['Caudoviricetes', 'Tailed phages'], ['Megaviricetes', 'Giant viruses'], ['Maveriviricetes', 'Virophages'],
    ['mk|major capsid', 'Capsid genes'], ['mk|terminase', 'Terminase'], ['mk|lysis (holin/endolysin/spanin)', 'Lysis genes'], ['mk|integrase/excisionase', 'Integrase/excisionase']];
  const wrap = el('div'); let k = 'all';
  vrChips(sec, OPTS, k, v => { k = v; draw(); }); sec.append(wrap);
  function draw() {
    wrap.innerHTML = '';
    const lab = OPTS.find(o => o[0] === k)[1];
    const box = fig(wrap, `${lab}: transcripts per million protein-coding read pairs`, 'Per-day means of the RNA libraries; solid = seawater, dashed = sea ice. Marker genes (capsid, terminase, lysis vs integrase) are summed over all confident viral contigs, so a rise in capsid and lysis genes signals lytic activity.');
    legendPools(box.parentNode.insertBefore(el('div'), box));
    const get = r => k === 'all' ? r.viral : k.startsWith('mk|') ? r.mk[k.slice(3)] : r.cls[k];
    timecourse(box, vrSeries(J.rna, get), { ylab: 'CPM', unit: 'CPM', h: 210, floor: 1 });
  }
  draw();
}

function vrDiversity(sec, J) {
  sec.append(el('h3', {}, 'Diversity'));
  const wrap = el('div'); let k = 'rich';
  vrChips(sec, [['rich', 'Active vOTUs (RNA)'], ['H', 'Shannon diversity (RNA)'], ...(J.dna_il ? [['ilrich', 'vOTUs present (Illumina DNA)'], ['ilH', 'Shannon (Illumina DNA)']] : [])], k, v => { k = v; draw(); });
  sec.append(wrap);
  function draw() {
    wrap.innerHTML = '';
    const il = k.startsWith('il'), key = il ? k.slice(2) : k;
    const box = fig(wrap, il ? 'Viral diversity in DNA' : 'Viral diversity in transcripts', il ? 'vOTUs with at least 2 mapped reads in 2 M read pairs; Shannon over vOTU read shares.' : 'Active vOTUs = vOTUs at ≥ 1 transcript per million protein-coding read pairs in that library; Shannon over the vOTU shares of viral transcripts. Richness depends on sequencing depth; the libraries were all subsampled to 10 M read pairs.');
    legendPools(box.parentNode.insertBefore(el('div'), box));
    timecourse(box, vrSeries(il ? J.dna_il : J.rna, r => r[key], r => r.medium !== 'other'), { ylab: key === 'rich' ? 'vOTUs' : 'Shannon H', h: 200, log: false, fmt: d3.format(key === 'rich' ? ',.0f' : '.2f') });
  }
  draw();
}

function vrPairs(sec, J) {
  sec.append(el('h3', {}, 'Phages and their hosts'));
  sec.append(el('p', { class: 'kicker' }, 'Seawater time courses of the most responsive vOTUs (solid, transcripts per million) beside their most strongly co-varying MAG (dashed, TPM; Spearman ρ across all 130 libraries). Co-variation suggests, but does not prove, a host.'));
  const sel = el('select', { id: 'vrpair', 'aria-label': 'Phage' }); J.pairs.forEach((p, i) => sel.append(new Option(p.label, i)));
  const c = el('div', { class: 'controls' }); c.append(el('label', { for: 'vrpair' }, 'vOTU'), sel); sec.append(c);
  const wrap = el('div'); sec.append(wrap);
  function draw() {
    wrap.innerHTML = '';
    const p = J.pairs[+sel.value];
    const box = fig(wrap, `${p.votu} · ${p.label}${p.len ? ` · ${d3.format('.1f')(p.len / 1000)} kb` : ''}`, p.host ? `Dashed: host candidate ${/^EE\d\d$/.test(p.host) ? 'eukaryotic MAG ' + p.host + ' (reads per million on its contigs, outlier-resistant profile)' : p.host + ' (TPM)'}. Colours: control pool blue, diesel pool orange.` : 'No host candidate shown.');
    const lg = el('div', { class: 'legend' }, `<span><i style="background:var(--ctl)"></i>Control pool</span><span><i style="background:var(--exp)"></i>Diesel pool</span><span style="color:var(--ink-2)"><i style="background:var(--ink-2)"></i>phage transcripts (CPM)</span>${p.host ? '<span style="color:var(--ink-2)"><i class="dash"></i>host MAG (TPM)</span>' : ''}`);
    box.parentNode.insertBefore(lg, box);
    const rows = p.pts.map(x => ({ ...x, medium: 'SW', name: `${x.pool} day ${x.day}` }));
    const ser = vrSeries(rows, r => r.phage);
    if (p.host) ser.push(...vrSeries(rows, r => r.host).map(s => ({ ...s, medium: 'SI', faint: true })));
    timecourse(box, ser, { ylab: 'CPM / TPM', h: 220, floor: 1 });
  }
  sel.onchange = draw; draw();
}

function vrHosts(sec, J) {
  const box = fig(sec, 'Diesel-responsive and dominant vOTUs with candidate hosts', 'Median transcripts per million in seawater after dosing, diesel vs control pool. BACI (2025-26): log₂ change in the diesel pool relative to its own pre-dosing level, minus the same for the control. Hosts: the MAG whose TPM co-varies most with the vOTU across all libraries (Spearman ρ).');
  const H = J.hosts.slice().sort((a, b) => Math.max(b['oil_2024-25'] || 0, b['oil_2025-26'] || 0) - Math.max(a['oil_2024-25'] || 0, a['oil_2025-26'] || 0));
  box.innerHTML = `<div class="tablewrap" style="max-height:420px"><table><thead><tr><th>vOTU</th><th>Class</th><th class="num">kb</th><th>CheckV</th><th class="num">24-25 diesel</th><th class="num">24-25 ctl</th><th class="num">25-26 diesel</th><th class="num">25-26 ctl</th><th class="num">BACI 25-26</th><th>Host candidate (ρ)</th></tr></thead><tbody>${
    H.map(h => `<tr><td class="mono">${h.votu}</td><td>${h.cls || ''}</td><td class="num">${h.length ? d3.format('.1f')(h.length / 1000) : ''}</td><td>${h.checkv_quality || ''}</td>${['oil_2024-25', 'ctl_2024-25', 'oil_2025-26', 'ctl_2025-26'].map(k => `<td class="num">${h[k] == null ? '–' : fmt1(h[k])}</td>`).join('')}<td class="num">${h['baci_2025-26'] == null ? '–' : `<span class="pill ${h['baci_2025-26'] > 0 ? 'up' : 'down'}">${vrSgn(h['baci_2025-26'])}</span>`}</td><td>${h.host1} (${h.rho1})</td></tr>`).join('')}</tbody></table></div>`;
}

function vrGroups(sec, J) {
  // no viral genome groups: the tailed-phage cluster is many phages (contigs stay individual vOTUs); the former V02 is plastid P01
  const p = J.phage_cluster;
  sec.append(el('p', { class: 'kicker' }, `${p ? `A cluster of tailed-phage contigs (about ${d3.format(',')(Math.round(p.contigs / 100) * 100)} contigs, ${d3.format('.0f')(p.mbp)} Mbp) sits inside the eukaryotic region of the contig map; it is many phages, not one genome. ` : ''}The group formerly listed here as virus V02 is algal plastid DNA, now plastid group P01 (<a href="#euks">Eukaryotes</a> page).`));
}

function vrProphage(sec, J) {
  sec.append(el('h3', {}, 'Prophages'));
  sec.append(el('p', { class: 'kicker' }, `geNomad found ${J.prophage_n.all} proviruses in the coassembly, ${J.prophage_n.in_mags} on MAG contigs. For each, the table gives its transcripts (CPM of protein-coding reads, median in seawater after dosing), the change in its share of the host MAG’s transcripts in 2025-26 relative to each pool’s pre-dosing level (BACI), and its activity (share of host transcripts divided by its share of host genes; > 1 = transcribed above the host’s average gene). Check which genes are transcribed before calling an induction: the strongest diesel-pool rises here come from host-like genes at provirus edges.`));
  const box = fig(sec, '', 'Sorted by the highest CPM in any library.');
  box.innerHTML = `<div class="tablewrap" style="max-height:380px"><table><thead><tr><th>Host</th><th>MAG</th><th class="num">kb</th><th class="num">Genes</th><th class="num">24-25 diesel</th><th class="num">24-25 ctl</th><th class="num">25-26 diesel</th><th class="num">25-26 ctl</th><th class="num">BACI share</th><th class="num">Activity 25-26</th></tr></thead><tbody>${
    J.prophage.map(p => `<tr><td>${p.host}</td><td class="mono">${p.bin}</td><td class="num">${p.kb}</td><td class="num">${p.n}</td><td class="num">${p.oil[0] ?? '–'}</td><td class="num">${p.ctl[0] ?? '–'}</td><td class="num">${p.oil[1] ?? '–'}</td><td class="num">${p.ctl[1] ?? '–'}</td><td class="num">${p.baci == null ? '–' : vrSgn(p.baci)}</td><td class="num">${p.act[1] == null ? '–' : d3.format('.1f')(p.act[1])}</td></tr>`).join('')}</tbody></table></div>`;
}

function vrIce(sec, J) {
  const box = fig(sec, 'Viruses in ice vs water', 'Paired log₂ ratios (same pool and day): ice minus seawater (23 pairs) and top minus bottom ice (13 pairs, 2024-25). Higher in ice = orange.');
  box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Metric</th><th class="num">ice − water</th><th class="num">pairs higher in ice</th><th class="num">top − bottom</th><th class="num">pairs higher at top</th></tr></thead><tbody>${
    J.ice.map(x => `<tr><td>${x.k === 'all viral' ? 'all viral transcripts' : x.k === 'active vOTUs' ? 'active vOTUs' : x.k}</td><td class="num"><span class="pill ${x.iw > 0 ? 'up' : 'down'}">${vrSgn(x.iw)}</span></td><td class="num">${x.iw_up}/${x.iw_n}</td><td class="num"><span class="pill ${x.tb > 0 ? 'up' : 'down'}">${vrSgn(x.tb)}</span></td><td class="num">${x.tb_up}/${x.tb_n}</td></tr>`).join('')}</tbody></table></div>`;
}

function vrConc(sec, J) {
  const C = J.conc;
  sec.append(el('h3', {}, 'Method test: virus concentration'));
  sec.append(el('p', { class: 'kicker' }, `Eight RNA libraries (VIRUS_*) test virus-concentration treatments: PEG precipitation, MgCl₂, acidification (pH) and their combinations, with VIRUS_NEG as the no-additive treatment. They are not survey samples and are left out of everything above. Their source is not recorded on the libraries; the 2025 field log suggests control-pool seawater from day 166 (5 May 2026, “Fernando’s protocol”: 8 mL of 1 M HCl per 4 L, against untreated controls). That link is inferred from the log. For reference they are compared with the latest control-pool seawater libraries (days 49–81).`));
  const box = fig(sec, 'Viral transcripts recovered per treatment', `Stacked by viral class, transcripts per million protein-coding read pairs. The treatments barely differ from each other (median Bray–Curtis ${C.bc_within} between their viral transcript profiles) and differ more from regular seawater (${C.bc_ref}). The no-additive library recovers as much as any additive.`);
  const rows = [...C.libs.map(x => ({ ...x, g: 'test' })), ...C.ref.map(x => ({ ...x, g: 'ref' }))];
  const KS = ['Caudoviricetes', 'Megaviricetes', 'Maveriviricetes'], cols = ['--s1', '--s3', '--s5'];
  const W = width(box), lw = Math.min(220, W * 0.4), rh = 22, H = rows.length * rh + 40;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const x = d3.scaleLinear().domain([0, d3.max(rows, r => r.viral) * 1.05]).range([lw + 10, W - 70]);
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${rows.length * rh + 14})`).call(d3.axisBottom(x).ticks(5).tickFormat(tickf));
  rows.forEach((r, i) => { const y = i * rh + 4 + (r.g === 'ref' ? 8 : 0); let x0 = 0;
    svg.append('text').attr('x', lw).attr('y', y + 12).attr('text-anchor', 'end').text(r.name);
    KS.forEach((k, j) => { const v = r.cls[k] || 0;
      svg.append('rect').attr('x', x(x0)).attr('y', y + 2).attr('width', Math.max(0, x(x0 + v) - x(x0) - 1)).attr('height', 14).attr('fill', css(cols[j]))
        .on('mousemove', ev => showTip(ev, `<b>${r.name}</b><br>${k}: ${fmt(v)} CPM<br>all viral ${fmt(r.viral)} · ${fmt(r.rich)} active vOTUs · H ${r.H}`)).on('mouseleave', hideTip); x0 += v; });
    svg.append('text').attr('x', x(r.viral) + 4).attr('y', y + 13).attr('class', 'mono').text(`${d3.format('.0f')(r.viral / 1000)}k`); });
  box.parentNode.insertBefore(el('div', { class: 'legend' }, KS.map((k, j) => `<span><i style="background:var(${cols[j]});height:10px"></i>${k}</span>`).join('') + '<span>lower block = control-pool seawater, days 49–81</span>'), box);
  const b2 = fig(sec, 'What concentration enriches', 'vOTUs most enriched in the concentrates (median of the eight) over late control-pool seawater (median), transcripts per million protein-coding read pairs. The virophage and algal giant viruses are recovered best.');
  b2.innerHTML = `<div class="tablewrap"><table><thead><tr><th>vOTU</th><th>Class</th><th>Family</th><th class="num">kb</th><th class="num">Concentrates</th><th class="num">Seawater</th><th class="num">log₂ ratio</th></tr></thead><tbody>${
    C.enriched.map(e => `<tr><td class="mono">${e.votu}</td><td>${e.cls}</td><td>${e.fam}</td><td class="num">${e.kb}</td><td class="num">${fmt(e.conc)}</td><td class="num">${fmt1(e.sw)}</td><td class="num">${vrSgn(e.lfc)}</td></tr>`).join('')}</tbody></table></div>`;
}

RENDER.viruses = virusesSec;
