/* Light section: oxygenic photosynthesis, aerobic anoxygenic photosynthesis (AAP), rhodopsins, light sensing; modelled clear-sky PAR.
   Data: light.json (light/scripts/light_site.py, par_model.py). Uses helpers from app.js (el, fig, width, css, showTip, hideTip, fmt, fmt1, tickf, timecourse, legendPools, OIL, POOL). */
const ltSgn = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.1f')(v);
const ltPct = d3.format('.0%');
const ltR = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.2f')(v);
const ltCat = (J, frag) => J.categories.find(c => c.name.includes(frag));
function ltSeries(J, values) {
  const out = {};
  J.samples.forEach((s, i) => { if (!['CTL', 'EXP'].includes(s.pool) || !['SW', 'SI'].includes(s.medium) || s.day == null || values[i] == null) return;
    const k = `${s.season}|${s.pool}|${s.medium}`; ((out[k] ||= {})[s.day] ||= []).push([values[i], s.name]); });
  return Object.entries(out).map(([k, byDay]) => { const [season, pool, medium] = k.split('|');
    return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a, x => x[0]), name: a.map(x => x[1]).join(', ') })) }; });
}

/* PAR strip: same geometry and x-domain as timecourse() in app.js, so it lines up with the chart drawn below it in the same box */
/* under-ice model lines (per pool) and logged roof state (vertical lines: solid open, dotted closed) */
function ltUnder(g, J, se, x, y) {
  ['CTL', 'EXP'].forEach(pool => {
    const U = J.underice.filter(d => d.season === se && d.pool === pool && d.day >= x.domain()[0] && d.day <= x.domain()[1] && d.under > 0);
    g.append('path').attr('d', d3.line().x(d => x(d.day)).y(d => y(Math.max(y.domain()[0], d.under)))(U)).attr('fill', 'none').attr('stroke', css(POOL[pool])).attr('stroke-width', 2);
  });
}
function ltRoof(g, J, se, x, ih) {
  (J.roof || []).filter(d => d.season === se && d.day >= x.domain()[0] && d.day <= x.domain()[1]).forEach(d => {
    const xx = x(d.day) + (d.pool === 'EXP' ? 1.5 : -1.5);
    g.append('line').attr('x1', xx).attr('x2', xx).attr('y1', 0).attr('y2', ih).attr('stroke', css(POOL[d.pool])).attr('stroke-width', d.open ? 1.6 : 1)
      .attr('stroke-dasharray', d.open ? null : '1 3').attr('opacity', d.open ? 0.8 : 0.6);
  });
}
function ltTip(J, se, p, d) {
  const u = pool => J.underice.find(q => q.season === se && q.pool === pool && q.day === d);
  const roof = (J.roof || []).filter(q => q.season === se && q.day === d).map(q => `${q.pool === 'EXP' ? 'diesel' : 'control'} roof ${q.open ? 'open' : 'closed'}`).join(', ');
  const row = (pool, lab) => { const q = u(pool); return q ? `<br>${lab}: ${q.under < 0.1 ? d3.format('.3f')(q.under) : fmt1(q.under)} under ${fmt(q.ice)} cm ice, ${fmt1(q.snow)} cm snow (${d3.format('.1%')(q.T)})` : ''; };
  return `<b>${p.date}</b> (day ${p.day})<br>clear-sky surface PAR ${fmt1(p.par)} mol m⁻² d⁻¹, day length ${fmt1(p.dl)} h${row('CTL', 'control')}${row('EXP', 'diesel')}${roof ? '<br>' + roof : ''}`;
}
const ltLegend = () => el('div', { class: 'legend' }, `<span><i style="background:var(--grid);height:10px"></i>clear-sky surface PAR</span><span><i style="background:var(--ctl)"></i>under ice, control pool</span><span><i style="background:var(--exp)"></i>under ice, diesel pool</span><span style="color:var(--ink-2)">│ roof open · ┆ roof closed (2025-26 log)</span>`);

function ltParStrip(box, J, series) {
  const W = width(box), gap = 28, pw = W < 620 ? W - 10 : (W - gap) / 2, cols = W < 620 ? 1 : 2, ph = 64;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', (ph + 6) * Math.ceil(2 / cols)).attr('role', 'img').attr('aria-label', 'modelled clear-sky PAR');
  const pmax = d3.max(J.par_daily, d => d.par);
  ['2024-25', '2025-26'].forEach((se, k) => {
    const gx = cols === 2 ? k * (pw + gap) : 0, gy = cols === 2 ? 0 : k * (ph + 6);
    const g = svg.append('g').attr('transform', `translate(${gx + 46},${gy + 6})`), iw = pw - 56, ih = ph - 22;
    const days = series.filter(s => s.season === se).flatMap(s => s.pts.map(p => p.day));
    const x = d3.scaleLinear().domain([Math.min(-7, d3.min(days) ?? -7), Math.max(22, d3.max(days) ?? 22)]).range([0, iw]).nice();
    const y = d3.scaleLog().domain([0.01, pmax * 1.2]).range([ih, 0]).clamp(true);
    const P = J.par_daily.filter(d => d.season === se && d.day >= x.domain()[0] && d.day <= x.domain()[1]);
    g.append('path').attr('d', d3.area().x(d => x(d.day)).y0(ih).y1(d => y(d.par))(P)).attr('fill', css('--grid')).attr('stroke', css('--ink-3')).attr('stroke-width', 1);
    ltUnder(g, J, se, x, y); ltRoof(g, J, se, x, ih);
    g.append('g').attr('class', 'axis').call(d3.axisLeft(y).tickValues([0.01, 0.1, 1, 10]).tickFormat(tickf));
    g.append('text').attr('x', iw).attr('y', 8).attr('text-anchor', 'end').style('font-size', '10.5px').attr('fill', css('--ink-3')).text('modelled PAR, mol m⁻² d⁻¹ (log): surface, under ice');
    g.append('line').attr('x1', x(OIL[se])).attr('x2', x(OIL[se])).attr('y1', 0).attr('y2', ih).attr('stroke', css('--ink-3')).attr('stroke-dasharray', '2 3');
    g.append('rect').attr('x', 0).attr('y', 0).attr('width', iw).attr('height', ih).attr('fill', 'transparent')
      .on('mousemove', ev => { const d = Math.round(x.invert(d3.pointer(ev)[0])); const p = P.find(q => q.day === d); if (p) showTip(ev, ltTip(J, se, p, d)); })
      .on('mouseleave', hideTip);
  });
}

function lightSec(sec) {
  const J = D.light;
  if (!J) { sec.innerHTML = '<h2>Light</h2><p class="kicker">light.json is missing.</p>'; return; }
  sec.innerHTML = `<h2>Light</h2><p class="kicker">${J.text.kicker}</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  J.text.stats.forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  ltSky(sec, J);
  ltTable(sec, J);
  ltExplorer(sec, J);
  ltAlgae(sec, J);
  ltGenomes(sec, J);
  ltPartition(sec, J);
  ltM1122(sec, J);
  sec.append(el('h3', {}, 'How the pieces fit, and what is still a hypothesis'));
  sec.append(el('ul', { class: 'kicker' }, J.text.reading.map(t => `<li>${t}</li>`).join('')));
}

/* ---------- modelled light across both winters ---------- */
function ltSky(sec, J) {
  sec.append(el('h3', {}, 'How much light there was to use'));
  sec.append(el('p', { class: 'kicker' }, J.text.sky));
  sec.append(el('p', { class: 'kicker' }, J.text.underice));
  const box = fig(sec, 'Modelled PAR at Churchill: clear-sky surface and under the ice', 'Daily photosynthetically active radiation above the ice and roof under a cloudless sky (Haurwitz clear-sky irradiance × 0.45 × 4.57 µmol J⁻¹, solar position from the NOAA equations at 58.76° N, 94.07° W). Grey: clear-sky surface PAR, an upper bound. Coloured lines: under-ice PAR per pool from the measured snow and ice (log scale; no roof, cloud or oil). Vertical pool-coloured lines: logged roof state (solid open, dotted closed). Ticks under the curve mark sampling days; the dashed line is dosing. Hover for date, snow, ice and transmission.');
  box.parentNode.insertBefore(ltLegend(), box);
  const W = width(box), gap = 28, cols = W < 620 ? 1 : 2, pw = cols === 2 ? (W - gap) / 2 : W - 10, ph = 190;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', (ph + 40) * Math.ceil(2 / cols)).attr('role', 'img').attr('aria-label', 'modelled PAR');
  const pmax = d3.max(J.par_daily, d => d.par);
  ['2024-25', '2025-26'].forEach((se, k) => {
    const gx = cols === 2 ? k * (pw + gap) : 0, gy = cols === 2 ? 0 : k * (ph + 40);
    const g = svg.append('g').attr('transform', `translate(${gx + 46},${gy + 18})`), iw = pw - 56, ih = ph - 34;
    const sd = J.samples.filter(s => s.season === se && ['CTL', 'EXP'].includes(s.pool) && s.day != null).map(s => s.day);
    const x = d3.scaleLinear().domain([Math.min(-10, d3.min(sd)), Math.max(30, d3.max(sd) + 3)]).range([0, iw]).nice();
    const y = d3.scaleLog().domain([0.01, pmax * 1.3]).range([ih, 0]).clamp(true);
    const P = J.par_daily.filter(d => d.season === se && d.day >= x.domain()[0] && d.day <= x.domain()[1]);
    g.append('g').attr('class', 'gridl').call(d3.axisLeft(y).tickValues([0.01, 0.1, 1, 10]).tickSize(-iw).tickFormat('')).select('.domain').remove();
    g.append('path').attr('d', d3.area().x(d => x(d.day)).y0(ih).y1(d => y(d.par))(P)).attr('fill', css('--grid'));
    g.append('path').attr('d', d3.line().x(d => x(d.day)).y(d => y(d.par))(P)).attr('fill', 'none').attr('stroke', css('--ink-3')).attr('stroke-width', 1.5);
    ltUnder(g, J, se, x, y); ltRoof(g, J, se, x, ih);
    g.append('g').attr('class', 'axis').attr('transform', `translate(0,${ih})`).call(d3.axisBottom(x).ticks(6));
    g.append('g').attr('class', 'axis').call(d3.axisLeft(y).tickValues([0.01, 0.1, 1, 10]).tickFormat(tickf));
    [...new Set(sd)].forEach(d => g.append('line').attr('x1', x(d)).attr('x2', x(d)).attr('y1', ih - 7).attr('y2', ih).attr('stroke', css('--ink'))
      .attr('stroke-width', 1.5));
    g.append('line').attr('x1', x(OIL[se])).attr('x2', x(OIL[se])).attr('y1', 0).attr('y2', ih).attr('stroke', css('--ink-3')).attr('stroke-dasharray', '2 3');
    g.append('text').attr('x', x(OIL[se]) + 4).attr('y', 10).text('diesel').attr('font-size', 10.5);
    const sol = P.reduce((a, b) => (b.par < a.par ? b : a), P[0]);
    if (sol) g.append('text').attr('x', x(sol.day)).attr('y', y(sol.par) - 8).attr('text-anchor', 'middle').style('font-size', '11px').attr('fill', css('--ink-2')).text(`${fmt1(sol.dl)} h of daylight`);
    g.append('text').attr('class', 'ttl').attr('x', 0).attr('y', -6).text(`CMO ${se}`);
    g.append('text').attr('x', iw).attr('y', ih + 28).attr('text-anchor', 'end').text('experiment day');
    if (k === 0) g.append('text').attr('transform', `translate(-34,${ih / 2}) rotate(-90)`).attr('text-anchor', 'middle').text('mol photons m⁻² d⁻¹');
    g.append('rect').attr('width', iw).attr('height', ih).attr('fill', 'transparent')
      .on('mousemove', ev => { const d = Math.round(x.invert(d3.pointer(ev)[0])); const p = P.find(q => q.day === d); if (p) showTip(ev, ltTip(J, se, p, d) + (sd.includes(d) ? '<br>sampling day' : '')); })
      .on('mouseleave', hideTip);
  });
}

/* ---------- inventory ---------- */
function ltTable(sec, J) {
  sec.append(el('h3', {}, 'The light genes at a glance'));
  const box = fig(sec, '', 'Median CPM in seawater after dosing. Ratios are log₂ (diesel ÷ control) of medians. BACI = change from each pool’s own five pre-dosing samples, diesel minus control (2025-26 seawater). Ice = paired ice minus water (24 pairs); Top = top minus bottom ice (13 pairs). r under ice = Spearman correlation with modelled under-ice PAR across 2025-26 seawater samples in each pool (diesel / control); it is largely a time trend. r | time = the same after removing a linear time trend within each pool (both pools pooled); none reaches q &lt; 0.05. Cell colour runs from blue (higher in control, water or bottom ice) to orange (higher in diesel, ice or top ice). Kaiju’s database has no eukaryotes, so it labels plastid genes on unbinned contigs as cyanobacteria (Synechococcus, Thermosynechococcus and so on): in the oxygenic rows those are mostly algal plastids. Click a row for its time course.');
  const col = v => { if (v == null) return 'transparent'; const t = Math.max(-1, Math.min(1, v / 4)); return t >= 0 ? d3.interpolateRgb(css('--panel'), css('--exp'))(t * 0.85) : d3.interpolateRgb(css('--panel'), css('--ctl'))(-t * 0.85); };
  const cell = v => `<td class="num" style="background:${col(v)}">${ltSgn(v)}</td>`;
  let html = `<div class="tablewrap" style="max-height:none"><table><thead><tr><th>Gene group · main hosts (all samples)</th><th class="num">Genes</th><th class="num">CPM diesel 24-25</th><th class="num">CPM diesel 25-26</th><th class="num">24-25</th><th class="num">25-26</th><th class="num">BACI 25-26</th><th class="num">Ice</th><th class="num">Top</th><th class="num">r under ice 25-26<br>diesel / control</th><th class="num">r | time<br>25-26</th></tr></thead><tbody>`;
  J.groups.forEach(g => {
    const cs = J.categories.filter(c => c.group === g); if (!cs.length) return;
    html += `<tr><td colspan="11" style="background:var(--grid);font-weight:650">${g}</td></tr>`;
    cs.forEach(c => { html += `<tr data-c="${c.name}"><td style="white-space:normal;min-width:230px">${c.name}<div style="font-size:11.5px;color:var(--ink-3);margin-top:2px">${
      c.hosts.all.slice(0, 3).map(h => `${h.host} ${ltPct(h.share)}`).join(' · ')}</div></td><td class="num">${fmt(c.n_genes)}</td><td class="num">${c.oil24 == null ? '–' : fmt1(c.oil24)}</td><td class="num">${c.oil25 == null ? '–' : fmt1(c.oil25)}</td>${cell(c.lfc24)}${cell(c.lfc25)}${cell(c.baci25)}${cell(c.ice)}${cell(c.top)}<td class="num">${c.runder ? ltR(c.runder.EXP) + ' / ' + ltR(c.runder.CTL) : '–'}</td><td class="num" title="q = ${c.qpart ?? '–'}">${ltR(c.rpart)}</td></tr>`; });
  });
  box.innerHTML = html + '</tbody></table></div>';
  box.querySelectorAll('tr[data-c]').forEach(tr => { tr.onclick = () => { const s = document.getElementById('ltsel'); if (s) { s.value = tr.dataset.c; s.onchange(); s.scrollIntoView({ behavior: 'smooth', block: 'center' }); } }; });
}

/* ---------- time-course explorer with PAR strip ---------- */
function ltExplorer(sec, J) {
  sec.append(el('h3', {}, 'Time courses against the light'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'ltsel', 'aria-label': 'Gene group' });
  J.groups.forEach(g => { const cs = J.categories.filter(x => x.group === g); if (!cs.length) return; const og = document.createElement('optgroup'); og.label = g; cs.forEach(x => og.append(new Option(x.name, x.name))); sel.append(og); });
  c.append(el('label', { for: 'ltsel' }, 'Gene group'), sel);
  const v = el('div'); sec.append(v);
  function draw() {
    v.innerHTML = ''; const cat = J.categories.find(x => x.name === sel.value);
    const box = fig(v, cat.name, 'Top strip: modelled clear-sky PAR on the same experiment-day axis. Below: CPM of protein-coding reads; lines join per-day means of replicate libraries.');
    legendPools(box.parentNode.insertBefore(el('div'), box));
    const ser = ltSeries(J, J.timecourse[cat.name]);
    ltParStrip(box, J, ser);
    timecourse(box, ser, { ylab: 'CPM', unit: 'CPM', h: 190, floor: 0.05 });
  }
  sel.onchange = draw; sel.value = (ltCat(J, 'reaction centre') || J.categories[0]).name; draw();
}

/* ---------- algae: photosynthesis share of plastid transcripts ---------- */
function ltAlgae(sec, J) {
  sec.append(el('h3', {}, 'Algae: fewer plastids, or dimmer ones?'));
  sec.append(el('p', { class: 'kicker' }, J.text.algae));
  const box = fig(sec, 'Photosynthesis genes as a share of plastid transcripts', `Photosystem I and II, cytochrome b6f, phycobilisome and RuBisCO transcripts on the ${J.plastid.n_bins} plastid MAGs, per mille of all transcripts on those MAGs. Samples with under 50 CPM on plastid MAGs are left out. Top strip: modelled clear-sky PAR.`);
  legendPools(box.parentNode.insertBefore(el('div'), box));
  const ser = ltSeries(J, J.plastid.share); ltParStrip(box, J, ser);
  timecourse(box, ser, { ylab: '‰ of plastid transcripts', unit: '‰', h: 180, log: false });
}

/* ---------- within-genome shares ---------- */
function ltGenomes(sec, J) {
  sec.append(el('h3', {}, 'Inside each light-using genome'));
  sec.append(el('p', { class: 'kicker' }, J.text.genomes));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const kind = el('select', { id: 'ltkind', 'aria-label': 'Light system' }); [['AAP', 'Aerobic anoxygenic photosynthesis'], ['rhodopsin', 'Proton-pumping rhodopsins']].forEach(([v_, l]) => kind.append(new Option(l, v_)));
  c.append(el('label', { for: 'ltkind' }, 'Light system'), kind);
  const v = el('div'); sec.append(v);
  function draw() {
    v.innerHTML = '';
    const R = J.genomes.filter(g => g.kind === kind.value && g['oil_2025-26'] > 0 && g['ctl_2025-26'] > 0).sort((a, b) => (b.host || 0) - (a.host || 0)).slice(0, 20)
      .sort((a, b) => Math.log2((b['oil_2025-26'] + 1e-3) / (b['ctl_2025-26'] + 1e-3)) - Math.log2((a['oil_2025-26'] + 1e-3) / (a['ctl_2025-26'] + 1e-3)));
    const box = fig(v, `${kind.value === 'AAP' ? 'AAP' : 'Rhodopsin'} share of each genome’s own transcripts, diesel vs control seawater (2025-26)`, 'Per mille of the genome’s protein-coding transcripts, median over seawater samples after dosing in which the genome had at least 200 CPM. Blue = control, orange = diesel. Up to 20 genomes, those with the most transcripts, among genomes expressing these genes in both pools. Hover for pre-dosing values, the ice-vs-water shift and the PAR correlation.');
    if (!R.length) { box.append(el('p', { class: 'note' }, 'No genome has enough reads in both pools in 2025-26.')); return; }
    const W = width(box), lw = Math.min(300, W * 0.4), rh = 20, H = R.length * rh + 40;
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
    const all = R.flatMap(d => [d['oil_2025-26'], d['ctl_2025-26']]).filter(x => x > 0);
    const x = d3.scaleLog().domain([Math.max(0.005, d3.min(all) * 0.8), d3.max(all) * 1.25]).range([lw + 10, W - 60]).clamp(true);
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${R.length * rh + 10})`).call(d3.axisBottom(x).ticks(5, '~g'));
    svg.append('text').attr('x', (lw + W) / 2).attr('y', H - 2).attr('text-anchor', 'middle').style('font-size', '11px').attr('fill', css('--ink-3')).text('‰ of the genome’s transcripts');
    R.forEach((d, i) => {
      const y = i * rh + 6, a = d['ctl_2025-26'] || 0.005, b = d['oil_2025-26'] || 0.005, l = Math.log2((b + 1e-3) / (a + 1e-3));
      svg.append('text').attr('x', lw).attr('y', y + 4).attr('text-anchor', 'end').style('font-size', '12px').text(d.taxon === d.bin ? d.bin : `${d.taxon.length > 32 ? d.taxon.slice(0, 31) + '…' : d.taxon} ${d.bin}`);
      svg.append('line').attr('x1', x(a)).attr('x2', x(b)).attr('y1', y).attr('y2', y).attr('stroke', css('--ink-3')).attr('stroke-width', 2);
      [[a, '--ctl'], [b, '--exp']].forEach(([val, cc]) => svg.append('circle').attr('cx', x(val)).attr('cy', y).attr('r', 5).attr('fill', css(cc)).attr('stroke', css('--panel')).attr('stroke-width', 2));
      svg.append('text').attr('x', W - 54).attr('y', y + 4).style('font-size', '12px').attr('fill', css('--ink-2')).text(ltSgn(l));
      svg.append('rect').attr('x', 0).attr('y', y - rh / 2).attr('width', W).attr('height', rh).attr('fill', 'transparent')
        .on('mousemove', ev => showTip(ev, `<b>${d.taxon} ${d.bin}</b> (${d.klass})<br>after dosing: ${d['ctl_2025-26']}‰ control, ${d['oil_2025-26']}‰ diesel<br>before dosing: ${d['pre_ctl_2025-26'] ?? '–'}‰ control, ${d['pre_exp_2025-26'] ?? '–'}‰ diesel<br>ice − water: ${ltSgn(d.ice_log2)} log₂ (${d.ice_up}/${d.ice_n} pairs up)<br>r with PAR, 2025-26 water: ${ltR(d['rpar_2025-26_EXP'])} diesel, ${ltR(d['rpar_2025-26_CTL'])} control`)).on('mouseleave', hideTip);
    });
    box.parentNode.insertBefore(el('div', { class: 'legend' }, '<span><i style="background:var(--ctl)"></i>control pool</span><span><i style="background:var(--exp)"></i>diesel pool</span>'), box);
    // per-genome time course
    const c2 = el('div', { class: 'controls' }); v.append(c2);
    const gsel = el('select', { id: 'ltgen', 'aria-label': 'Genome' });
    J.genomes.filter(g => g.kind === kind.value).sort((p, q) => (q.host || 0) - (p.host || 0)).forEach(g => gsel.append(new Option(`${g.taxon} (${g.bin})`, g.bin)));
    c2.append(el('label', { for: 'ltgen' }, 'Genome'), gsel);
    const v2 = el('div'); v.append(v2);
    function d2() { v2.innerHTML = ''; const g = J.genomes.find(q => q.kind === kind.value && q.bin === gsel.value);
      const b2 = fig(v2, `${g.taxon} ${g.bin}: ${kind.value === 'AAP' ? 'AAP' : 'rhodopsin'} share of its own transcripts`, 'Per mille of the genome’s protein-coding transcripts; samples where the genome had under 200 CPM are left out. Top strip: modelled clear-sky PAR.');
      legendPools(b2.parentNode.insertBefore(el('div'), b2)); const ser = ltSeries(J, g.series); ltParStrip(b2, J, ser);
      timecourse(b2, ser, { ylab: '‰ of genome', unit: '‰', h: 180, floor: 0.01 }); }
    gsel.onchange = d2; d2();
  }
  kind.onchange = draw; draw();
}

/* ---------- AAP vs rhodopsin by class ---------- */
function ltPartition(sec, J) {
  sec.append(el('h3', {}, 'Two ways to harvest light'));
  sec.append(el('p', { class: 'kicker' }, J.text.partition));
  const box = fig(sec, '', 'Share of all reaction-centre (AAP) and proton-pumping rhodopsin transcripts, summed over every library, by host class (GTDB for MAGs, Kaiju for unbinned contigs).');
  const cats = [...new Set(J.partition.map(p => p.cat))];
  box.innerHTML = `<div class="tablewrap"><table><thead><tr>${cats.map(c => `<th>${c}: host class</th><th class="num">share</th>`).join('')}</tr></thead><tbody>${
    d3.range(8).map(i => `<tr>${cats.map(c => { const p = J.partition.filter(q => q.cat === c)[i]; return p ? `<td>${p.klass}</td><td class="num">${ltPct(p.share)}</td>` : '<td></td><td></td>'; }).join('')}</tr>`).join('')}</tbody></table></div>`;
}

/* ---------- Sulfitobacter metabat_1122 ---------- */
function ltM1122(sec, J) {
  const M = J.m1122; if (!M) return;
  sec.append(el('h3', {}, 'Sulfitobacter metabat_1122: light when hungry?'));
  sec.append(el('p', { class: 'kicker' }, J.text.m1122));
  const box = fig(sec, 'metabat_1122: AAP share vs its physiological state', 'Each point is a library with enough reads on this genome. x = its PC1 score on the Physiological states page (negative = the nitrogen-hungry water state that loads PufH and the ammonium transporter); y = AAP genes (puf, bch, crt) as per mille of the genome’s transcripts. Circles seawater, triangles sea ice; open 2024-25, filled 2025-26.');
  const P = M.points.map(p => ({ ...p, s: J.samples.find(s => s.id === p.id) })).filter(p => p.s && p.aap != null);
  const W = width(box), H = Math.min(380, Math.max(260, W * 0.5)), mg = { l: 50, r: 14, t: 10, b: 38 };
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const x = d3.scaleLinear().domain(d3.extent(P, p => p.pc1)).nice().range([mg.l, W - mg.r]);
  const y = d3.scaleLog().domain([Math.max(0.01, d3.min(P, p => p.aap || 1) * 0.8), d3.max(P, p => p.aap) * 1.2]).range([H - mg.b, mg.t]).clamp(true);
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H - mg.b})`).call(d3.axisBottom(x).ticks(6));
  svg.append('g').attr('class', 'axis').attr('transform', `translate(${mg.l},0)`).call(d3.axisLeft(y).ticks(4, '~g'));
  svg.append('text').attr('x', W - mg.r).attr('y', H - 4).attr('text-anchor', 'end').style('font-size', '11px').attr('fill', css('--ink-3')).text('PC1 (Physiological states)');
  svg.append('text').attr('transform', `translate(12,${mg.t}) rotate(-90)`).attr('text-anchor', 'end').style('font-size', '11px').attr('fill', css('--ink-3')).text('AAP ‰ of genome');
  svg.append('g').selectAll('path').data(P).join('path').attr('d', p => d3.symbol(p.s.medium === 'SI' ? d3.symbolTriangle : d3.symbolCircle, 60)())
    .attr('transform', p => `translate(${x(p.pc1)},${y(Math.max(p.aap, 0.01))})`)
    .attr('fill', p => p.s.season === '2024-25' ? css('--panel') : css(POOL[p.s.pool])).attr('stroke', p => p.s.season === '2024-25' ? css(POOL[p.s.pool]) : css('--panel')).attr('stroke-width', 1.6)
    .on('mousemove', (ev, p) => showTip(ev, `<b>${p.s.name}</b><br>${p.s.season}, day ${p.s.day}, ${p.s.medium === 'SI' ? 'sea ice' : 'seawater'}<br>AAP ${p.aap}‰ · N acquisition ${p.nacq ?? '–'}‰ · PC1 ${p.pc1}<br>clear-sky PAR ${p.s.par ?? '–'} mol m⁻² d⁻¹`)).on('mouseleave', hideTip);
  box.parentNode.insertBefore(el('div', { class: 'legend' }, '<span><i style="background:var(--ctl)"></i>control pool</span><span><i style="background:var(--exp)"></i>diesel pool</span><span style="color:var(--ink-2)">● seawater</span><span style="color:var(--ink-2)">▲ sea ice</span><span style="color:var(--ink-2)">○ open = 2024-25</span>'), box);
}

RENDER.light = lightSec;
