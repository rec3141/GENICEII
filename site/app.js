/* Churchill diesel mesocosm site: one page, hash-addressable sections, D3 charts drawn from data.json (+ cyclo.json, hgt.json, ice.json, islands.json; ice and islands sections in ice.js, islands.js). */
const $ = (s, el = document) => el.querySelector(s);
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const POOL = { CTL: '--ctl', EXP: '--exp', other: '--other' };
const fmt = d3.format(',.0f'), fmt1 = d3.format(',.1f'), fmt2 = d3.format('.2f');
const tickf = v => v >= 1000 ? d3.format('~s')(v) : d3.format('~g')(v);
const SECTIONS = [
  ['overview', 'Overview'], ['community', 'Community response'], ['pathways', 'Degradation pathways'],
  ['cyclo', 'Cycloclasticus'], ['porticocc', '50-400-T64'], ['states', 'Physiological states'], ['guilds', 'Response guilds'], ['light', 'Light'], ['sulfur', 'Sulfur & hydrogen'], ['metals', 'Metals & toxicity'], ['nitrogen', 'Nitrogen'], ['oxygen', 'Oxygen'], ['euks', 'Eukaryotes'], ['viruses', 'Viruses'], ['hgt', 'Gene transfer'], ['islands', 'Genomic islands'], ['ice', 'Ice vs water'], ['qc', 'Data & caveats']];
// nav: top-level links and dropdown groups (ids from SECTIONS)
const NAV = [['overview'], ['Diesel response', ['community', 'pathways', 'cyclo', 'porticocc', 'states']], ['Biogeochemistry', ['light', 'sulfur', 'nitrogen', 'oxygen', 'metals']],
  ['Food web', ['euks', 'viruses']], ['Genomes & ice', ['ice', 'hgt', 'islands', 'guilds']], ['qc']];
const OIL = { '2024-25': 5.5, '2025-26': 1 };
let D = {}, tip;

function showTip(ev, html) { tip.innerHTML = html; tip.style.opacity = 1; const w = tip.offsetWidth, h = tip.offsetHeight;
  let x = ev.clientX + 14, y = ev.clientY + 14; if (x + w > innerWidth - 8) x = ev.clientX - w - 14; if (y + h > innerHeight - 8) y = ev.clientY - h - 14;
  tip.style.left = x + 'px'; tip.style.top = y + 'px'; }
function hideTip() { tip.style.opacity = 0; }
function el(tag, attrs = {}, html = '') { const e = document.createElement(tag); for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v); if (html) e.innerHTML = html; return e; }
// a chart title naming exactly one MAG gets that MAG's full lineage (site/lineage.json, scripts/lineage.py) as a subtitle
const MAG_RX = /\b((?:dastool-)?(?:comebin|metabat|lorbin|semibin|vamb|maxbin|concoct|binette)_\d+(?:_sub)?|E\d\d)\b/g;   // MAG ids and eukaryotic genome groups (E##)
function magLineage(title) {
  if (!title || !D.lineage) return '';
  const ids = [...new Set([...String(title).replace(/<[^>]*>/g, '').matchAll(MAG_RX)].map(m => m[1].replace('dastool-', '')))];
  return ids.length === 1 && D.lineage[ids[0]] ? D.lineage[ids[0]] : '';
}
function fig(parent, title, cap) { const f = el('div', { class: 'fig' }); if (title) f.append(el('div', { class: 'title' }, title));
  const lin = magLineage(title); if (lin) f.append(el('div', { class: 'lineage' }, lin)); const box = el('div', { class: 'chartbox' }); f.append(box);
  if (cap) f.append(el('p', { class: 'cap' }, cap)); parent.append(f); return box; }
function legendPools(parent, withMedium = true) {
  const l = el('div', { class: 'legend' });
  l.innerHTML = `<span><i style="background:var(--ctl)"></i>Control pool</span><span><i style="background:var(--exp)"></i>Diesel pool</span>` +
    (withMedium ? `<span style="color:var(--ink-2)"><i style="background:var(--ink-2)"></i>seawater</span><span style="color:var(--ink-2)"><i class="dash"></i>sea ice</span>` : '');
  parent.append(l); return l; }
function width(box, max = 1060) { return Math.max(300, Math.min(max, box.clientWidth || 700)); }

/* ---------- time-course small multiple: one panel per season ---------- */
function timecourse(box, series, opts = {}) {
  // series: [{season, pool, medium, pts:[{day, v, name}]}]
  const W = width(box), seasons = ['2024-25', '2025-26'], gap = 28, ph = opts.h || 190;
  const pw = W < 620 ? W - 10 : (W - gap) / 2, cols = W < 620 ? 1 : 2;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', (ph + 40) * Math.ceil(2 / cols) + 4).attr('role', 'img').attr('aria-label', opts.label || 'time course');
  const vals = series.flatMap(s => s.pts.map(p => p.v)).filter(v => v != null);
  const log = opts.log !== false;
  const vmax = d3.max(vals) || 1, vmin = log ? Math.max(opts.floor || 0.1, d3.min(vals.filter(v => v > 0)) || 0.1) : Math.min(0, d3.min(vals) ?? 0);
  seasons.forEach((se, k) => {
    const gx = cols === 2 ? k * (pw + gap) : 0, gy = cols === 2 ? 0 : k * (ph + 40);
    const g = svg.append('g').attr('transform', `translate(${gx + 46},${gy + 18})`);
    const iw = pw - 56, ih = ph - 34;
    const days = series.filter(s => s.season === se).flatMap(s => s.pts.map(p => p.day));
    const x = d3.scaleLinear().domain([Math.min(-7, d3.min(days) ?? -7), Math.max(22, d3.max(days) ?? 22)]).range([0, iw]).nice();
    const y = (log ? d3.scaleLog().domain([vmin, vmax * 1.3]).clamp(true) : d3.scaleLinear().domain([vmin < 0 ? vmin * 1.08 : 0, vmax * 1.08])).range([ih, 0]);
    g.append('g').attr('class', 'gridl').call(d3.axisLeft(y).ticks(4, log ? '~s' : '~s').tickSize(-iw).tickFormat('')).select('.domain').remove();
    g.append('g').attr('class', 'axis').attr('transform', `translate(0,${ih})`).call(d3.axisBottom(x).ticks(6));
    g.append('g').attr('class', 'axis').call(log ? d3.axisLeft(y).ticks(4).tickFormat(v => { const e = Math.log10(v); return Math.abs(e - Math.round(e)) < 1e-6 ? tickf(v) : ''; }) : d3.axisLeft(y).ticks(4).tickFormat(tickf));
    g.append('line').attr('x1', x(OIL[se])).attr('x2', x(OIL[se])).attr('y1', 0).attr('y2', ih).attr('stroke', css('--ink-3')).attr('stroke-dasharray', '2 3');
    g.append('text').attr('x', x(OIL[se]) + 4).attr('y', 10).text('diesel').attr('font-size', 10.5);
    g.append('text').attr('class', 'ttl').attr('x', 0).attr('y', -6).text(`CMO ${se}`);
    g.append('text').attr('x', iw).attr('y', ih + 28).attr('text-anchor', 'end').text('experiment day');
    if (k === 0 && opts.ylab) g.append('text').attr('transform', `translate(-38,${ih / 2}) rotate(-90)`).attr('text-anchor', 'middle').text(opts.ylab);
    series.filter(s => s.season === se).forEach(s => {
      const pts = s.pts.filter(p => p.v != null && (!log || p.v > 0)).sort((a, b) => a.day - b.day);
      const c = css(POOL[s.pool] || '--other');
      const line = d3.line().x(p => x(p.day)).y(p => y(Math.max(vmin, p.v)));
      if (!s.points_only) g.append('path').attr('d', line(pts)).attr('fill', 'none').attr('stroke', c).attr('stroke-width', 2).attr('stroke-dasharray', s.medium === 'SI' ? '5 3' : null).attr('opacity', s.faint ? .45 : 1);
      g.selectAll(null).data(pts).join('path').attr('d', d3.symbol(s.medium === 'SI' ? d3.symbolTriangle : (s.points_only ? d3.symbolDiamond : d3.symbolCircle), s.points_only ? 46 : 30)())
        .attr('transform', p => `translate(${x(p.day)},${y(Math.max(vmin, p.v))})`).attr('fill', s.points_only ? 'none' : c).attr('stroke', s.points_only ? c : css('--panel')).attr('stroke-width', s.points_only ? 1.6 : 1)
        .on('mousemove', (ev, p) => showTip(ev, `<b>${p.name || ''}</b><br>${s.pool} ${s.medium} · day ${p.day}<br>${opts.unit || ''} ${opts.fmt ? opts.fmt(p.v) : fmt1(p.v)}`)).on('mouseleave', hideTip);
    });
  });
}
function groupSeries(values, opts = {}) {
  // values aligned with D.samples; average replicate samples per day
  const out = {};
  D.samples.forEach((s, i) => { if (s.pool === 'other' || s.medium === 'other' || s.day == null) return; const v = values[i]; if (v == null) return;
    const k = `${s.season}|${s.pool}|${s.medium}`; (out[k] ||= {}); (out[k][s.day] ||= []).push([v, s.name]); });
  return Object.entries(out).map(([k, byDay]) => { const [season, pool, medium] = k.split('|');
    return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a, x => x[0]), name: a.map(x => x[1]).join(', ') })) }; });
}

/* ---------- sections ---------- */
function overview(sec) {
  sec.innerHTML = `<h2>What happened in the pools</h2>
  <p class="kicker">Each winter the observatory filled two outdoor pools with local seawater, let sea ice grow, and released diesel into one pool (2024-25: day 5; 2025-26: about 136 L on day 1, when the ice was still thin and patchy). We sequenced RNA from seawater and ice through the season and asked which organisms switched on, which genes they used, and whether those genes move between genomes.</p>
  <div class="stat-row" id="stats"></div>
  <div id="findings"></div>
  <div id="design"></div>`;
  const st = $('#stats', sec), n = D.samples.length;
  const stats = [[n, 'RNA libraries (all reads; median 31 M read pairs)'], [D.samples.filter(s => s.pool === 'EXP' && s.phase === 'oiled').length, 'diesel-pool samples after dosing'],
    [911, 'MAGs in the marine coassembly'], [D.dna_samples.length, 'CMO nanopore metagenomes (DNA)']];
  stats.forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${fmt(v)}</div><div class="l">${l}</div>`)));
  const F = (D.findings || []);
  const fd = $('#findings', sec);
  F.forEach((f, i) => fd.append(el('div', { class: 'finding' }, `<div class="n">${i + 1}</div><div><b>${f.h}</b><br>${f.t}</div>`)));
  // sampling design strip
  const box = fig($('#design', sec), 'Sampling design', 'Each mark is one RNA library, placed by experiment day. Circles are seawater, triangles sea ice, and grey squares are other material (frost flowers, sea foam, suspected marine oil snow, unlabelled late samples, virus-concentration tests, kelp). The dashed line marks the diesel release.');
  const W = width(box), rows = [['2024-25', 'CTL'], ['2024-25', 'EXP'], ['2025-26', 'CTL'], ['2025-26', 'EXP']], rh = 34;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', rows.length * rh + 40);
  const x = d3.scaleLinear().domain([-8, 85]).range([120, W - 10]);
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${rows.length * rh + 8})`).call(d3.axisBottom(x).ticks(10));
  rows.forEach(([se, p], k) => {
    const y = 14 + k * rh;
    svg.append('text').attr('x', 0).attr('y', y + 4).attr('class', 'ttl').text(`${se} ${p === 'CTL' ? 'control' : 'diesel'}`);
    svg.append('line').attr('x1', x(-8)).attr('x2', x(85)).attr('y1', y).attr('y2', y).attr('stroke', css('--grid'));
    if (p === 'EXP') svg.append('line').attr('x1', x(OIL[se])).attr('x2', x(OIL[se])).attr('y1', y - 12).attr('y2', y + 12).attr('stroke', css('--ink-3')).attr('stroke-dasharray', '2 2');
    const pts = D.samples.filter(s => s.season === se && (s.pool === p) && s.day != null);
    svg.selectAll(null).data(pts).join('path').attr('d', s => d3.symbol(s.medium === 'SI' ? d3.symbolTriangle : s.medium === 'SW' ? d3.symbolCircle : d3.symbolSquare, 42)())
      .attr('transform', (s, i) => `translate(${x(s.day)},${y + ((i % 3) - 1) * 6})`).attr('fill', s => s.medium === 'other' ? css('--other') : css(POOL[p])).attr('stroke', css('--panel'))
      .on('mousemove', (ev, s) => showTip(ev, `<b>${s.name}</b><br>${s.date} · day ${s.day}<br>${s.env || s.medium}`)).on('mouseleave', hideTip);
  });
}

function community(sec) {
  sec.innerHTML = `<h2>Community response</h2>
  <p class="kicker">Transcripts were assigned to MAGs through their contigs. Season separates the samples more than anything else. Within each season the diesel pool already differed from the control before dosing; in 2025-26 the gap widened only after day 20.</p>
  <div class="grid2"><div id="pca"></div><div id="perm"></div></div><div id="succ"></div><div id="explorer"></div><div id="mgk"></div><div id="euk"></div>`;
  // PCA
  const box = fig($('#pca', sec), 'Transcript composition by MAG', 'Aitchison PCA (centred log-ratio of read pairs per MAG). Hollow marks are samples taken before diesel was added. Hover a point for its sample.');
  legendPools(box.parentNode.insertBefore(el('div'), box), false);
  const W = Math.min(width(box), 560), H = 360, m = { l: 40, r: 10, t: 10, b: 34 };
  const pts = D.samples.filter(s => s.pc[0] != null);
  const x = d3.scaleLinear().domain(d3.extent(pts, s => s.pc[0])).nice().range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain(d3.extent(pts, s => s.pc[1])).nice().range([H - m.b, m.t]);
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6));
  svg.append('g').attr('class', 'axis').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6));
  svg.append('text').attr('x', W - m.r).attr('y', H - 4).attr('text-anchor', 'end').text('PC1 (47.1%, season)');
  svg.append('text').attr('transform', `translate(11,${m.t + 4}) rotate(90)`).text('PC2 (13.0%)');
  svg.selectAll(null).data(pts).join('path').attr('d', s => d3.symbol(s.medium === 'SI' ? d3.symbolTriangle : s.medium === 'SW' ? d3.symbolCircle : d3.symbolSquare, 46)())
    .attr('transform', s => `translate(${x(s.pc[0])},${y(s.pc[1])})`)
    .attr('fill', s => s.phase === 'pre' ? 'none' : css(POOL[s.pool])).attr('stroke', s => s.phase === 'pre' ? css(POOL[s.pool]) : css('--panel')).attr('stroke-width', s => s.phase === 'pre' ? 1.6 : 1)
    .on('mousemove', (ev, s) => showTip(ev, `<b>${s.name}</b><br>${s.season} · ${s.pool} ${s.medium} · day ${s.day ?? '–'}`)).on('mouseleave', hideTip);
  ['2024-25', '2025-26'].forEach(t => { const ps = pts.filter(s => s.season === t); const cx = d3.mean(ps, s => s.pc[0]);
    svg.append('text').attr('class', 'ttl').attr('x', x(cx) + (t === '2024-25' ? 60 : 0)).attr('y', t === '2024-25' ? y(0) : y(d3.max(ps, s => s.pc[1])) - 14).attr('text-anchor', 'middle').text(`CMO ${t}`); });
  // variance partition table
  const pb = fig($('#perm', sec), 'How much each factor explains', 'PERMANOVA on Aitchison distance (999 permutations; all p ≤ 0.005 except 2024-25 medium, p = 0.12). One pool per treatment means time points within a pool are pseudo-replicates, so read these as effect sizes rather than tests.');
  const rows = [['Season', 0.501], ['Pool, both seasons', 0.115], ['Phase (pre/control/oiled)', 0.199], ['Season × pool × medium', 0.724],
    ['2024-25: pool', 0.447], ['2024-25: oiled vs not', 0.427], ['2024-25: ice vs water', 0.031], ['2025-26: pool', 0.286], ['2025-26: oiled vs not', 0.318], ['2025-26: ice vs water', 0.105]];
  const bw = Math.min(width(pb), 520), bh = 22, sv = d3.select(pb).append('svg').attr('width', bw).attr('height', rows.length * bh + 24);
  const lab = bw < 460 ? 170 : 210, bx = d3.scaleLinear().domain([0, 0.75]).range([0, bw - lab - 50]);
  rows.forEach(([l, v], i) => { const g = sv.append('g').attr('transform', `translate(0,${i * bh + 6})`);
    g.append('text').attr('x', 0).attr('y', 12).text(l);
    g.append('rect').attr('x', lab).attr('y', 3).attr('height', 12).attr('width', bx(v)).attr('rx', 3).attr('fill', css('--accent'));
    g.append('text').attr('x', lab + 4 + bx(v)).attr('y', 13).text(`${v.toFixed(2)}`); });
  succession($('#succ', sec));
  magExplorer($('#explorer', sec));
  mgPanel($('#mgk', sec));
  eukaryotes($('#euk', sec));
}


function succession(parent) {
  const F = D.families; if (!F) return;
  parent.append(el('h3', {}, 'Succession by family'));
  const box = fig(parent, '', 'Share of MAG-assigned transcripts by bacterial family, seawater samples, averaged per day. Cycloclasticaceae expand in the diesel pool in both winters. Sulfurimonadaceae hold steady in the 2025-26 diesel pool while declining in the control. In 2024-25 a Sulfitobacter genome that was already more abundant in the diesel pool peaks at days 7–12, and Cycloclasticus dominates by day 55 (no samples between days 20 and 55).');
  const cols = ['--s2', '--s1', '--s4', '--s3', '--s7', '--s5', '--s6', '--rule'];
  const lg = el('div', { class: 'legend' }, F.order.map((k, i) => `<span><i style="background:var(${cols[i]});height:10px"></i>${k}</span>`).join('')); box.parentNode.insertBefore(lg, box);
  const W = width(box), pw = W < 620 ? W : (W - 24) / 2, ph = 170;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', (W < 620 ? 4 : 2) * (ph + 34));
  let k = 0;
  ['2024-25', '2025-26'].forEach(se => ['EXP', 'CTL'].forEach(pool => {
    const ss = D.samples.map((s, i) => [s, i]).filter(([s]) => s.season === se && s.pool === pool && s.medium === 'SW' && s.day != null);
    const days = [...new Set(ss.map(([s]) => s.day))].sort((a, b) => a - b);
    const rows = days.map(d => { const ii = ss.filter(([s]) => s.day === d).map(([, i]) => i); const o = { day: d }; F.order.forEach(f => o[f] = d3.mean(ii, i => F.pct[f][i] || 0)); return o; });
    const st = d3.stack().keys(F.order)(rows);
    const col = W < 620 ? 0 : (se === '2024-25' ? 0 : 1), row = W < 620 ? k : (pool === 'EXP' ? 0 : 1);
    const g = svg.append('g').attr('transform', `translate(${col * (pw + 24) + 36},${row * (ph + 34) + 16})`);
    const x = d3.scaleLinear().domain(d3.extent(days)).range([0, pw - 46]), y = d3.scaleLinear().domain([0, 100]).range([ph - 20, 0]);
    g.selectAll(null).data(st).join('path').attr('d', d3.area().x(d => x(d.data.day)).y0(d => y(d[0])).y1(d => y(d[1])).curve(d3.curveLinear)).attr('fill', (d, i) => css(cols[i])).attr('stroke', css('--panel')).attr('stroke-width', .5)
      .on('mousemove', (ev, d) => showTip(ev, `<b>${d.key}</b>`)).on('mouseleave', hideTip);
    g.append('g').attr('class', 'axis').attr('transform', `translate(0,${ph - 20})`).call(d3.axisBottom(x).ticks(6));
    g.append('g').attr('class', 'axis').call(d3.axisLeft(y).ticks(3).tickFormat(v => v + '%'));
    g.append('text').attr('class', 'ttl').attr('x', 0).attr('y', -4).text(`CMO ${se} · ${pool === 'EXP' ? 'diesel' : 'control'} pool, seawater`);
    if (pool === 'EXP') g.append('line').attr('x1', x(OIL[se])).attr('x2', x(OIL[se])).attr('y1', 0).attr('y2', ph - 20).attr('stroke', css('--ink')).attr('stroke-dasharray', '2 3');
    k++; }));
}


function eukaryotes(parent) {
  const E = (D.euk || []).slice().sort((a, b) => (a['2025-26']?.lfc ?? 0) - (b['2025-26']?.lfc ?? 0)); if (!E.length) return;
  parent.append(el('h3', {}, 'Eukaryotes'));
  const box = fig(parent, 'Change in eukaryotic rRNA share with diesel', 'log₂ ratio of each lineage’s share of eukaryotic rRNA reads (SILVA, order level) in diesel-pool versus control-pool samples after dosing. Euglenozoans decline in both winters; dinoflagellates, Phaeocystis, cryptophytes and Cercozoa decline mainly in 2025-26. Labyrinthulomycetes (thraustochytrids, known hydrocarbon-tolerant protists) rise in 2024-25. Ochrophyta, mostly diatoms, hold their share.');
  const lg = el('div', { class: 'legend' }, `<span><i style="background:var(--s7);height:8px;width:8px;border-radius:50%"></i>2024-25</span><span><i style="background:var(--s3);height:8px;width:8px;border-radius:50%"></i>2025-26</span>`); box.parentNode.insertBefore(lg, box);
  const W = Math.min(width(box), 760), lw = 280, rh = 22, svg = d3.select(box).append('svg').attr('width', W).attr('height', 44 + E.length * rh);
  const x = d3.scaleLinear().domain([-3.5, 2.5]).range([lw, W - 20]);
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${E.length * rh + 6})`).call(d3.axisBottom(x).ticks(6));
  svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 0).attr('y2', E.length * rh + 6).attr('stroke', css('--ink-3'));
  E.forEach((e, i) => { const y = i * rh + 12;
    svg.append('text').attr('x', lw - 8).attr('y', y + 4).attr('text-anchor', 'end').text(e.taxon);
    const a = e['2024-25']?.lfc, b = e['2025-26']?.lfc;
    if (a != null && b != null) svg.append('line').attr('x1', x(a)).attr('x2', x(b)).attr('y1', y).attr('y2', y).attr('stroke', css('--rule')).attr('stroke-width', 2);
    [['2024-25', '--s7'], ['2025-26', '--s3']].forEach(([se, c]) => { const v = e[se]; if (!v) return;
      svg.append('circle').attr('cx', x(Math.max(-3.5, Math.min(2.5, v.lfc)))).attr('cy', y).attr('r', 5).attr('fill', css(c)).attr('stroke', css('--panel'))
        .on('mousemove', ev => showTip(ev, `<b>${e.taxon}</b> ${se}<br>diesel ${v.oil}% vs control ${v.ctl}% of eukaryotic rRNA`)).on('mouseleave', hideTip); }); });
  svg.append('text').attr('x', W - 20).attr('y', E.length * rh + 34).attr('text-anchor', 'end').text('log₂ (diesel / control)');
}


function mgPanel(parent) {
  const M = D.mg; if (!M) return;
  parent.append(el('h3', {}, 'DNA through the whole season (Illumina metagenomes)'));
  parent.append(el('p', { class: 'kicker' }, `Nanopore DNA sampling stopped around day 11, before the 2024-25 Cycloclasticus bloom. To follow DNA through the season, ${M.samples.length} CMO Illumina metagenomes (Genohub 2599682/3694692 for 2024-25 and 1343765 for 2025-26) were screened for k-mers unique to selected genomes and plasmids. The first 400,000 reads of each library were used. Values are the percentage of sampled read k-mers that are unique to the target. Sequencing errors and k-mers shared with relatives make this an underestimate of the read share, but it is comparable across samples.`));
  const c = el('div', { class: 'controls' }); parent.append(c);
  const sel = el('select', { id: 'mgsel', 'aria-label': 'Genome or plasmid' }); M.targets.forEach(t => sel.append(new Option(t.label, t.id))); c.append(el('label', { for: 'mgsel' }, 'Target'), sel);
  const v = el('div'); parent.append(v);
  function draw() { v.innerHTML = ''; const t = sel.value, vals = M.values[t];
    const box = fig(v, M.targets.find(x => x.id === t).label, 'Each point is one metagenome; lines join the per-day means. Log scale; values below about 0.001% are near the detection floor of 400,000 reads.');
    legendPools(box.parentNode.insertBefore(el('div'), box));
    const ser = {}; M.samples.forEach((s, i) => { if (s.pool === 'other' || s.day == null || s.medium === 'other') return; const k = `${s.season}|${s.pool}|${s.medium}`; ((ser[k] ||= {})[s.day] ||= []).push([vals[i], s.name]); });
    timecourse(box, Object.entries(ser).map(([k, byDay]) => { const [season, pool, medium] = k.split('|'); return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a, x => x[0]), name: a.map(x => x[1]).join(', ') })) }; }),
      { ylab: '% of read k-mers', unit: '%', h: 190, floor: 0.0003, fmt: d3.format('.3~g') }); }
  sel.onchange = draw; sel.value = 'metabat_9342'; draw();
}

function magExplorer(parent) {
  const wrap = el('div'); parent.append(wrap);
  wrap.innerHTML = `<h3>MAG explorer</h3><p class="kicker">Every MAG that differed strongly between pools after dosing (|ΔCLR| > 2.5, sample-level q < 0.05 in any season × medium, or a before-after-control-impact interaction; q-values treat samples as replicates, see Data & caveats), plus all MAGs that carry a near-complete degradation pathway. Pick a row to see its transcript time course and its DNA abundance in the nanopore metagenomes.</p>
  <div class="controls"><input type="search" id="magq" placeholder="Filter by taxon or bin" aria-label="Filter MAGs"><span id="magn" class="mono"></span></div>
  <div class="grid2"><div class="tablewrap" id="magtab"></div><div id="magplot"></div></div>`;
  const ids = Object.keys(D.mags);
  const cols = [['bin', 'MAG', b => b], ['label', 'Taxon', b => D.mags[b].label], ['c', 'Compl.', b => D.mags[b].completeness],
    ['a', 'Δ 24-25 SW', b => D.mags[b].resp['2024-25 SW']?.dclr], ['b', 'Δ 25-26 SW', b => D.mags[b].resp['2025-26 SW']?.dclr],
    ['i', 'BACI 25-26', b => D.mags[b].baci['2025-26']?.i], ['t', 'max TPM', b => d3.max(D.mags[b].tpm)]];
  let sortK = 'b', dir = -1, sel = null;
  const tab = $('#magtab', wrap);
  function draw() {
    const q = $('#magq', wrap).value.toLowerCase();
    let rows = ids.filter(b => !q || b.toLowerCase().includes(q) || Object.values(D.mags[b].tax).join(' ').toLowerCase().includes(q));
    const f = cols.find(c => c[0] === sortK)[2];
    rows.sort((a, b) => dir * (((f(a) ?? -1e9) > (f(b) ?? -1e9)) ? 1 : -1));
    $('#magn', wrap).textContent = `${rows.length} MAGs`;
    tab.innerHTML = `<table><thead><tr>${cols.map(c => `<th class="${c[0] === 'bin' || c[0] === 'label' ? '' : 'num'}" data-k="${c[0]}" aria-sort="${sortK === c[0] ? (dir > 0 ? 'ascending' : 'descending') : 'none'}">${c[1]}</th>`).join('')}</tr></thead><tbody>` +
      rows.map(b => `<tr data-b="${b}" class="${b === sel ? 'sel' : ''}">${cols.map(c => { const v = c[2](b);
        if (c[0] === 'a' || c[0] === 'b' || c[0] === 'i') return `<td class="num">${v == null ? '–' : `<span class="pill ${v > 0 ? 'up' : 'down'}">${v > 0 ? '+' : ''}${v.toFixed(1)}</span>`}</td>`;
        if (c[0] === 'bin') return `<td class="mono">${b}</td>`; if (c[0] === 'label') return `<td>${v}</td>`;
        return `<td class="num">${v == null ? '–' : c[0] === 't' ? fmt(v) : v}</td>`; }).join('')}</tr>`).join('') + '</tbody></table>';
    tab.querySelectorAll('th').forEach(th => th.onclick = () => { const k = th.dataset.k; dir = sortK === k ? -dir : -1; sortK = k; draw(); });
    tab.querySelectorAll('tbody tr').forEach(tr => tr.onclick = () => { sel = tr.dataset.b; draw(); plot(sel); });
  }
  function plot(b) {
    const box = $('#magplot', wrap); box.innerHTML = '';
    const M = D.mags[b];
    const inner = fig(box, `${b} · ${M.label}`, `${M.tax.p} › ${M.tax.c} › ${M.tax.o} › ${M.tax.f}${M.size ? ` · ${(M.size / 1e6).toFixed(2)} Mbp, ${M.completeness}% complete, ${M.redundancy}% redundant (DAS Tool SCGs)` : ''}`);
    legendPools(inner.parentNode.insertBefore(el('div'), inner));
    timecourse(inner, groupSeries(M.tpm), { ylab: 'transcripts per million', unit: 'TPM', h: 180 });
    if (M.dna) {
      const d2 = fig(box, 'DNA abundance (nanopore metagenomes)', '% of binned coverage. Open diamonds are single DNA samples; most nanopore sampling stopped by day 11.');
      const ser = {}; D.dna_samples.forEach((s, i) => { if (s.pool === 'other' || s.day == null || s.medium === 'other') return; const k = `${s.season}|${s.pool}|${s.medium}`; (ser[k] ||= []).push({ day: s.day, v: M.dna[i], name: s.name }); });
      timecourse(d2, Object.entries(ser).map(([k, pts]) => { const [season, pool, medium] = k.split('|'); return { season, pool, medium, pts, points_only: true }; }), { ylab: '% of binned DNA', unit: '%', h: 160, floor: 0.001, fmt: d3.format('.3f') });
    }
  }
  $('#magq', wrap).oninput = draw; draw();
  sel = ids.sort((a, b) => (D.mags[b].baci['2025-26']?.i ?? -9) - (D.mags[a].baci['2025-26']?.i ?? -9))[0];
  draw(); plot(sel);
}

function qc(sec) {
  sec.innerHTML = `<h2>Data &amp; caveats</h2><div id="qcfig"></div><div id="cav"></div>`;
  const box = fig($('#qcfig', sec), 'Where each library’s reads went', 'Share of read pairs per library, ordered by season, pool, medium and day. “rRNA (escaped filter)” is rRNA that the k-mer filter missed and that then mapped to rRNA genes in the assembly; in the 2025-26 control pool most of it is one euglenid 18S gene (contig_1340804). Expression values use protein-coding genes only, so this lowers depth rather than biasing them.');
  const keys = [['rrna_filtered', 'rRNA (filtered)', '--seq2'], ['rrna_escaped', 'rRNA (escaped filter)', '--s5'], ['cds', 'protein-coding genes', '--s3'], ['other_mapped', 'mapped, no gene', '--seq1'], ['unmapped', 'unmapped or short local hit', '--grid']];
  const l = el('div', { class: 'legend' }, keys.map(k => `<span><i style="background:var(${k[2]});height:10px"></i>${k[1]}</span>`).join('')); box.parentNode.insertBefore(l, box);
  const W = Math.max(width(box), 760);
  ['2024-25', '2025-26'].forEach(se => {
    const ss = D.samples.filter(s => s.season === se).sort((a, b) => d3.ascending(a.pool === 'other', b.pool === 'other') || d3.ascending(a.pool, b.pool) || d3.ascending(a.medium, b.medium) || d3.ascending(a.day ?? 999, b.day ?? 999));
    const H = 170, svg = d3.select(box).append('svg').attr('width', W).attr('height', H + 26);
    const x = d3.scaleBand().domain(ss.map(s => s.id)).range([40, W - 4]).padding(0.12), y = d3.scaleLinear().domain([0, 100]).range([H, 16]);
    svg.append('text').attr('class', 'ttl').attr('x', 40).attr('y', 11).text(`CMO ${se}`);
    svg.append('g').attr('class', 'axis').attr('transform', 'translate(40,0)').call(d3.axisLeft(y).ticks(4).tickFormat(d => d + '%'));
    ss.forEach(s => { let acc = 0; keys.forEach(([k, , c]) => { const v = s.qc[k] || 0;
      svg.append('rect').attr('x', x(s.id)).attr('width', x.bandwidth()).attr('y', y(acc + v)).attr('height', Math.max(0, y(acc) - y(acc + v))).attr('fill', css(c))
        .on('mousemove', ev => showTip(ev, `<b>${s.name}</b><br>${s.pool} ${s.medium} day ${s.day ?? '–'}<br>${keys.map(kk => `${kk[1]}: ${s.qc[kk[0]]}%`).join('<br>')}`)).on('mouseleave', hideTip); acc += v; }); });
    svg.selectAll(null).data(ss).join('rect').attr('x', s => x(s.id)).attr('width', x.bandwidth()).attr('y', H + 4).attr('height', 5).attr('fill', s => css(POOL[s.pool]));
    svg.append('text').attr('x', W - 4).attr('y', H + 22).attr('text-anchor', 'end').text('bar under each column: pool colour');
  });
  const meth = el('div'); sec.insertBefore(meth, $('#cav', sec));
  meth.innerHTML = `<h3>Methods in brief</h3><ul>
  <li><b>RNA.</b> 130 BGI DNBSEQ libraries, used at full sequencing depth (4.37 G read pairs in total, 8.8–72 M per library; an earlier version of this site used the first 10 M pairs of each). They were quality-trimmed, rRNA was removed (bbduk k31 against SILVA, with a subsample kept for rRNA taxonomy), and human reads were removed. Reads were aligned with strobealign to the 11 Gbp marine Flye coassembly; only alignments covering ≥80% of the read were kept. Read pairs on protein-coding genes (Bakta) were counted with featureCounts, and rRNA loci found by barrnap were masked. All expression values are per million protein-coding read pairs (CPM) or length-normalised (TPM).</li>
  <li><b>Genomes.</b> MAGs come from the danaSeq <span class="mono">mag_analysis</span> pipeline: DAS Tool over seven binners, GTDB-Tk r226, and CheckM2 for the Cycloclasticus bins. Bins without a GTDB call are labelled by the majority Kaiju genus of their contigs.</li>
  <li><b>Function.</b> KEGG orthologs come from eggNOG-mapper, remapped onto the current Bakta genes by identical coordinates. KEGG module completeness is computed per MAG from those KOs, and hydrocarbon genes are confirmed with CANT-HYD HMMs at trusted cutoffs.</li>
  <li><b>DNA.</b> 72 CMO nanopore metagenomes were mapped to the same assembly; MAG abundance is the length-weighted mean depth as a share of binned coverage.</li>
  <li><b>Statistics.</b> Statistics are descriptive. Diesel-pool samples after dosing are compared with control-pool samples (Mann-Whitney, Benjamini-Hochberg q). Before-after-control-impact interactions use 2025-26 seawater with permuted pool labels. Community structure uses Aitchison PCA and PERMANOVA. Pseudoreplication applies throughout.</li>
  <li><b>Mobile DNA.</b> Mobile elements come from geNomad (strict plasmid = conjugation or hallmark genes), IslandPath and IntegronFinder, plus transposases within 10 kb. Cross-taxon sharing was found with MMseqs2 (≥95% aa, ≥80% coverage), then filtered to ≥99% nt identity, independent DNA depth profiles and flank alignment.</li></ul>`;
  const cav = $('#cav', sec);
  (D.caveats || []).forEach(c => cav.append(el('div', { class: 'note' + (c.warn ? ' warn' : '') }, `<b>${c.h}.</b> ${c.t}`)));
  if (D.identity) {
    const f = fig(cav, 'Sample identity check across platforms (2024-25)', 'Before dosing the diesel pool already carried Sulfitobacter lorbin_032, Casp-alpha2 comebin_2596 and Sulfurimonas semibin_088 (the “diesel-pool signature”), while the control pool carried Yoonia lorbin_042 and UBA5078 semibin_072. Every library of the same sample should agree. Units differ by platform: nanopore % of binned DNA coverage, Illumina % of read k-mers, RNA % of binned transcripts. Rows show all early seawater libraries plus every corrected or flagged one; the larger of the two markers is highlighted.');
    const tw = el('div', { class: 'tablewrap' }); f.append(tw);
    const now = r => r.corrected ? `<span class="pill up">relabelled</span> ${r.corrected}` : r.flag === 'low_yield' ? '<span class="pill na">excluded: failed library</span>' : r.flag === 'pool_mismatch' ? '<span class="pill na">excluded: wrong pool, no partner</span>' : '<span style="color:var(--ink-3)">as labelled</span>';
    const v = (x, hi, c) => `<td class="num" style="${hi ? `color:var(${c});font-weight:650` : ''}">${x == null ? '–' : x.toFixed(2)}%</td>`;
    tw.innerHTML = `<table><thead><tr><th>Platform</th><th>Library</th><th>Label</th><th>Now</th><th class="num">Diesel-pool signature</th><th class="num">Control markers</th></tr></thead><tbody>${D.identity.map(r =>
      `<tr><td>${r.platform}</td><td class="mono">${r.library}</td><td>${r.label}</td><td>${now(r)}</td>${v(r.diesel, r.diesel > r.control, '--exp')}${v(r.control, r.control >= r.diesel, '--ctl')}</tr>`).join('')}</tbody></table>`; }
}

/* pathway, cyclo and hgt renderers are attached below when their data exist */
const RENDER = { overview, community, qc };


function degradationMap(parent) {
  const P = D.pathways, F = P.families;
  const fc = (k, se) => { const v = F[k]?.cpm; if (!v) return null;
    const o = D.samples.map((s, i) => [s, v[i]]).filter(([s]) => s.season === se && (s.medium === 'SW' || s.medium === 'SI'));
    const a = d3.mean(o.filter(([s]) => s.phase === 'oiled'), x => x[1]), b = d3.mean(o.filter(([s]) => s.phase === 'control'), x => x[1]);
    return { r: (a + 0.5) / (b + 0.5), a, b }; };
  const host = k => { const h = F[k]?.hosts?.['2025-26 oiled']?.top?.[0], h1 = F[k]?.hosts?.['2024-25 oiled']?.top?.[0]; return [h1, h]; };
  const cols = [
    { t: 'Diesel fraction', nodes: [['n-alkanes C10–C25', 'subst'], ['short alkanes, C1–C4', 'subst'], ['BTEX: benzene, toluene, xylenes', 'subst'], ['naphthalenes, biphenyl (2–3 rings)', 'subst']] },
    { t: 'Activating oxygenase', nodes: [['alkB', 'enz', 'Alkane monooxygenase AlkB'], ['CYP153 (product)', 'enz', 'P450 alkane hydroxylase CYP153'], ['pmoA/mmoX (methane MO)', 'enz', 'Cu membrane MO (pmoA-like)'], ['prmA (propane MO)', 'enz', 'Propane MO (prmA)'],
      ['tmoA-F (toluene MO)', 'enz', 'Toluene 4-MO (tmoA–F)'], ['dmpKLMNOP (phenol MO)', 'enz', 'Phenol hydroxylase (dmpK–P)'], ['todC (benzene/toluene DO)', 'enz', 'Benzene/toluene DO (todC1C2)'], ['nahAc/Ad (naphthalene DO)', 'enz', 'Naphthalene DO (nahAc/Ad)'], ['bphA (biphenyl DO)', 'enz', 'Biphenyl DO (bphA)']] },
    { t: 'Intermediate', nodes: [['fatty acids → β-oxidation', 'mid'], ['alcohols, aldehydes', 'mid'], ['(methyl)catechols', 'mid'], ['dihydroxy-naphthalene / -biphenyl → salicylate', 'mid']] },
    { t: 'Ring cleavage', nodes: [['catE/xylE (catechol 2,3-DO)', 'enz', 'Catechol 2,3-DO, meta (catE)'], ['nahC/bphC (extradiol DO)', 'enz', 'Extradiol DO (nahC/bphC)'], ['catA (catechol 1,2-DO)', 'enz', 'Catechol 1,2-DO, ortho (catA)'], ['pcaGH (protocatechuate 3,4-DO)', 'enz', 'Protocatechuate 3,4-DO'], ['gentisate 1,2-DO', 'enz', 'Gentisate 1,2-DO']] },
    { t: 'Central metabolism', nodes: [['acetyl-CoA, pyruvate → TCA', 'end']] }];
  const edges = [[0, 0, 1, 0], [0, 0, 1, 1], [0, 1, 1, 2], [0, 1, 1, 3], [0, 2, 1, 4], [0, 2, 1, 5], [0, 2, 1, 6], [0, 3, 1, 7], [0, 3, 1, 8],
    [1, 0, 2, 1], [1, 1, 2, 1], [1, 2, 2, 1], [1, 3, 2, 1], [2, 1, 2, 0], [1, 4, 2, 2], [1, 5, 2, 2], [1, 6, 2, 2], [1, 7, 2, 3], [1, 8, 2, 3],
    [2, 2, 3, 0], [2, 2, 3, 2], [2, 3, 3, 1], [2, 3, 3, 4], [2, 2, 3, 3], [2, 0, 4, 0], [3, 0, 4, 0], [3, 1, 4, 0], [3, 2, 4, 0], [3, 3, 4, 0], [3, 4, 4, 0]];
  const box = fig(parent, 'Diesel degradation map', 'Each enzyme box shows the ratio of mean expression in diesel-pool to control-pool samples (seawater and ice) for each winter, and the genome with the largest share of that enzyme’s transcripts in the diesel pool. Ratios use transcripts per million protein-coding reads with a +0.5 pseudocount. Orange marks a rise of 4× or more. Enzymes are grouped by KEGG ortholog. The Cycloclasticus “alkB” genes are divergent alkB-like proteins that fall below the curated CANT-HYD AlkB cutoff; by CANT-HYD, canonical AlkB transcripts in the diesel pool come mainly from Pseudosulfitobacter, and canonical CYP153 transcripts from Cycloclasticus metabat_9342 and unbinned Sphingorhabdus in 2024-25 and from unbinned Alcanivorax and Sphingorhabdus in 2025-26. Ratios here are means with seawater and ice pooled, so they are smaller than the seawater medians quoted on the overview. The third line is the 2025-26 before-after-control-impact ratio in seawater: (diesel pool after ÷ before) ÷ (control after ÷ before), medians. It removes the pre-dosing difference between pools and is the more conservative number.');
  const W = Math.max(980, width(box)), colw = [150, 250, 170, 250, 120], gapx = (W - d3.sum(colw)) / 4, bh = 72, gap = 10;
  const xs = colw.map((w, i) => d3.sum(colw.slice(0, i)) + i * gapx);
  const H = 40 + 9 * (bh + gap);
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const pos = cols.map((c, i) => { const n = c.nodes.length, tot = n * bh + (n - 1) * gap, y0 = 30 + (H - 30 - tot) / 2; return c.nodes.map((_, j) => [xs[i], y0 + j * (bh + gap), colw[i]]); });
  cols.forEach((c, i) => svg.append('text').attr('class', 'ttl').attr('x', xs[i]).attr('y', 14).text(c.t));
  edges.forEach(([a, b, c2, d]) => { const p1 = pos[a][b], p2 = pos[c2][d]; const x1 = p1[0] + p1[2], y1 = p1[1] + bh / 2, x2 = p2[0], y2 = p2[1] + bh / 2;
    svg.append('path').attr('d', `M${x1},${y1} C${(x1 + x2) / 2},${y1} ${(x1 + x2) / 2},${y2} ${x2},${y2}`).attr('fill', 'none').attr('stroke', css('--rule')).attr('stroke-width', 1.5); });
  const pill = v => v == null ? '–' : v.r >= 1 ? `${v.r < 10 ? v.r.toFixed(1) : Math.round(v.r)}×` : `÷${(1 / v.r).toFixed(1)}`;
  cols.forEach((c, i) => c.nodes.forEach((n, j) => { const [x, y, w] = pos[i][j]; const g = svg.append('g').attr('transform', `translate(${x},${y})`);
    if (n[1] !== 'enz') { g.append('rect').attr('width', w).attr('height', bh).attr('rx', 8).attr('fill', css('--bg')).attr('stroke', css('--rule'));
      const words = n[0].split(' '); let line = '', ly = 0; const lines = []; words.forEach(wd => { if ((line + ' ' + wd).length > (w / 7)) { lines.push(line); line = wd; } else line = (line ? line + ' ' : '') + wd; }); lines.push(line);
      lines.forEach((l, k) => g.append('text').attr('x', 10).attr('y', bh / 2 - (lines.length - 1) * 7 + k * 14 + 4).text(l).attr('fill', css('--ink'))); return; }
    const f1 = fc(n[0], '2024-25'), f2 = fc(n[0], '2025-26'), [h1, h2] = host(n[0]);
    const strong = [f1, f2].some(f => f && f.r >= 4);
    g.append('rect').attr('width', w).attr('height', bh).attr('rx', 6).attr('fill', css('--panel')).attr('stroke', strong ? css('--exp') : css('--rule')).attr('stroke-width', strong ? 2 : 1);
    g.append('text').attr('x', 8).attr('y', 16).text(n[2]).attr('class', 'ttl').attr('font-size', 11.5);
    [[f1, '24-25', h1], [f2, '25-26', h2]].forEach(([f, lab, h], k) => { const yy = 32 + k * 15;
      g.append('text').attr('x', 8).attr('y', yy).text(`${lab}: ${pill(f)}`).attr('font-weight', 650).attr('fill', f && f.r >= 4 ? css('--exp') : f && f.r < 0.5 ? css('--ctl') : css('--ink-2'));
      if (h) g.append('text').attr('x', 82).attr('y', yy).text(`${h.label.replace(' (unbinned, Kaiju)', '*').replace(' (Kaiju)', '').slice(0, 22)} ${h.pct}%`).attr('fill', css('--ink-3')); });
    const bc = F[n[0]]?.baci26; if (bc != null) g.append('text').attr('x', 8).attr('y', 63).text(`25-26 vs own baseline: ${bc >= 1 ? (bc < 10 ? bc.toFixed(1) : Math.round(bc)) + '×' : '÷' + (1 / bc).toFixed(1)}`).attr('fill', bc >= 4 ? css('--exp') : css('--ink-3')).attr('font-size', 10.5);
    g.style('cursor', 'pointer').on('click', () => { const sel = $('#famsel'); if (sel) { sel.value = n[0]; sel.dispatchEvent(new Event('change')); sel.scrollIntoView({ behavior: 'smooth', block: 'center' }); } })
      .on('mousemove', ev => showTip(ev, `<b>${n[2]}</b><br>2024-25: ${f1 ? `${fmt1(f1.a)} vs ${fmt1(f1.b)} CPM` : '–'}<br>2025-26: ${f2 ? `${fmt1(f2.a)} vs ${fmt1(f2.b)} CPM` : '–'}<br>click to open below`)).on('mouseleave', hideTip); }));
  svg.append('text').attr('x', W - 4).attr('y', H - 4).attr('text-anchor', 'end').text('* unbinned contigs, genus from Kaiju');
}


function markers(parent) {
  const M = D.markers; if (!M) return;
  parent.append(el('h3', {}, 'Beyond hydrocarbons: a hydrogen and sulfur economy'));
  parent.append(el('p', { class: 'kicker' }, 'Diesel coincides with chemolithoautotrophy. Pyruvate:ferredoxin oxidoreductase, NiFe hydrogenase, the reverse-TCA enzyme ATP-citrate lyase and sulfide dehydrogenase, mostly from Sulfurimonas, are higher in the diesel pool. In 2025-26 they were already high in both pools before dosing; the diesel pool held that level while the control declined. Sulfate-reduction genes (dsrA) stay low, so the sulfide and hydrogen probably do not come from in-pool sulfate reducers. High-affinity cytochrome bd oxidase is 1.4–2.7× higher in diesel-pool ice and seawater, consistent with low-oxygen microsites, though no pre-dosing ice exists to confirm it.'));
  const box = fig(parent, '', 'Mean CPM per group. Filled marks: control pool after day 0 (blue) or diesel pool after dosing (orange). Hollow marks: samples before dosing. Circles are seawater, triangles sea ice. Log scale.');
  const W = Math.max(760, width(box)), lw = 290, pw = (W - lw - 30) / 2, rh = 30;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', 40 + M.length * rh);
  const all = M.flatMap(m => Object.values(m.groups)).filter(v => v > 0);
  const x = d3.scaleLog().domain([Math.max(0.05, d3.min(all)), d3.max(all) * 1.3]).range([0, pw - 20]).clamp(true);
  ['2024-25', '2025-26'].forEach((se, k) => { const g = svg.append('g').attr('transform', `translate(${lw + k * (pw + 10)},24)`);
    g.append('text').attr('class', 'ttl').attr('x', 0).attr('y', -8).text(`CMO ${se}`);
    g.append('g').attr('class', 'axis').attr('transform', `translate(0,${M.length * rh})`).call(d3.axisBottom(x).ticks(4).tickFormat(v => { const e = Math.log10(v); return Math.abs(e - Math.round(e)) < 1e-6 ? tickf(v) : ''; }));
    M.forEach((m, i) => { const y = i * rh + rh / 2;
      g.append('line').attr('x1', 0).attr('x2', pw - 20).attr('y1', y).attr('y2', y).attr('stroke', css('--grid'));
      Object.entries(m.groups).forEach(([key, v]) => { const [s2, pool, med, ph] = key.split('|'); if (s2 !== se || !(v > 0)) return;
        const pre = ph === 'pre';
        g.append('path').attr('d', d3.symbol(med === 'SI' ? d3.symbolTriangle : d3.symbolCircle, 52)()).attr('transform', `translate(${x(v)},${y})`)
          .attr('fill', pre ? 'none' : css(POOL[pool])).attr('stroke', pre ? css(POOL[pool]) : css('--panel')).attr('stroke-width', pre ? 1.6 : 1)
          .on('mousemove', ev => showTip(ev, `<b>${m.name}</b> (${m.ko})<br>${se} ${pool} ${med} ${ph}<br>${fmt1(v)} CPM`)).on('mouseleave', hideTip); }); }); });
  M.forEach((m, i) => svg.append('text').attr('x', lw - 10).attr('y', 24 + i * rh + rh / 2 + 4).attr('text-anchor', 'end').text(`${m.id} · ${m.name}`));
}


function potential(parent) {
  const R = (D.pathways.potential || []).filter(r => r.encoded > 0).sort((a, b) => b.encoded - a.encoded);
  if (!R.length) return;
  parent.append(el('h3', {}, 'Potential versus realised'));
  parent.append(el('p', { class: 'kicker' }, 'Many genomes carry degradation pathways; few use them. Grey bars count MAGs that encode at least 75% of a pathway. Orange bars count those that express it above 1 CPM on average in the diesel pool. The last column names the genome with the largest share of that expression.'));
  const box = fig(parent, '', 'Classes listed are the most common GTDB classes among the encoding MAGs; “no GTDB call” marks refined sub-bins without a GTDB-Tk classification.');
  const W = Math.max(860, width(box)), lw = 300, bw = 220, rh = 24, svg = d3.select(box).append('svg').attr('width', W).attr('height', 26 + R.length * rh);
  const x = d3.scaleSqrt().domain([0, d3.max(R, r => r.encoded)]).range([0, bw]);
  svg.append('text').attr('class', 'ttl').attr('x', lw).attr('y', 12).text('MAGs encoding / expressing (square-root scale)');
  svg.append('text').attr('class', 'ttl').attr('x', lw + bw + 70).attr('y', 12).text('largest share of expression');
  R.forEach((r, i) => { const y = 22 + i * rh;
    svg.append('text').attr('x', lw - 8).attr('y', y + 13).attr('text-anchor', 'end').text(D.pathways.defs[r.p].name.split(',')[0].slice(0, 46));
    svg.append('rect').attr('x', lw).attr('y', y + 3).attr('height', 14).attr('rx', 3).attr('width', x(r.encoded)).attr('fill', css('--rule'));
    svg.append('rect').attr('x', lw).attr('y', y + 3).attr('height', 14).attr('rx', 3).attr('width', x(r.expressed)).attr('fill', css('--exp'));
    svg.append('text').attr('x', lw + x(r.encoded) + 6).attr('y', y + 14).text(`${r.expressed} / ${r.encoded}`);
    svg.append('text').attr('x', lw + bw + 70).attr('y', y + 14).text(r.top ? `${r.top_share}% · ${r.top}` : '–')
      .on('mousemove', ev => showTip(ev, `<b>${D.pathways.defs[r.p].name}</b><br>${Object.entries(r.classes).map(([k, v]) => `${k || 'no GTDB call'}: ${v}`).join('<br>')}`)).on('mouseleave', hideTip); });
}

function pathways(sec) {
  const P = D.pathways;
  sec.innerHTML = `<h2>Degradation pathways</h2>
  <p class="kicker">Genes were assigned KEGG orthologs (eggNOG, mapped onto the current gene calls by coordinates) and grouped into KEGG modules for aromatic and alkane degradation, plus a few product-defined functions such as CYP153 and any ring-hydroxylating dioxygenase. A MAG “encodes” a pathway when its genes cover the module steps; it “expresses” it when reads land on those genes.</p>
  <div id="dmap"></div><div id="fams"></div><div id="heat"></div><div id="pot"></div><div id="mk"></div>`;
  potential($('#pot', sec));
  $('#mk', sec).append(el('p', { class: 'note' }, 'The hydrogen and sulfur economy that runs alongside the diesel response (Sulfurimonas, Cycloclasticus hydrogenases and nitrate respiration, sulfate reducers) now has its own page: <a href="#sulfur">Sulfur &amp; hydrogen</a>.'));
  degradationMap($('#dmap', sec));
  // ---- gene families: time course + hosts
  const famBox = $('#fams', sec);
  famBox.append(el('h3', {}, 'Key enzyme families'));
  const order = ['tmoA-F (toluene MO)', 'todC (benzene/toluene DO)', 'nahAc/Ad (naphthalene DO)', 'bphA (biphenyl DO)', 'nahC/bphC (extradiol DO)', 'catE/xylE (catechol 2,3-DO)',
    'dmpKLMNOP (phenol MO)', 'alkB', 'alkT/rubB', 'CYP153 (product)', 'pmoA/mmoX (methane MO)', 'prmA (propane MO)', 'catA (catechol 1,2-DO)', 'pcaGH (protocatechuate 3,4-DO)', 'gentisate 1,2-DO', 'nidAB (PAH DO)'];
  const NICE = { 'tmoA-F (toluene MO)': 'Toluene-4-monooxygenase (tmoA–F)', 'todC (benzene/toluene DO)': 'Benzene/toluene dioxygenase (todC1C2)', 'nahAc/Ad (naphthalene DO)': 'Naphthalene dioxygenase (nahAc/Ad)',
    'bphA (biphenyl DO)': 'Biphenyl dioxygenase (bphA)', 'nahC/bphC (extradiol DO)': 'Dihydroxy-naphthalene/biphenyl extradiol dioxygenase (nahC/bphC)', 'catE/xylE (catechol 2,3-DO)': 'Catechol 2,3-dioxygenase, meta cleavage (catE/xylE)',
    'dmpKLMNOP (phenol MO)': 'Multicomponent phenol hydroxylase (dmpK–P)', alkB: 'Alkane 1-monooxygenase (alkB)', 'alkT/rubB': 'Rubredoxin reductase (alkT)', 'CYP153 (product)': 'Cytochrome P450 alkane hydroxylase (CYP153)',
    'pmoA/mmoX (methane MO)': 'Copper membrane monooxygenase, pmoA-like', 'prmA (propane MO)': 'Propane monooxygenase (prmA)', 'catA (catechol 1,2-DO)': 'Catechol 1,2-dioxygenase, ortho cleavage (catA)',
    'pcaGH (protocatechuate 3,4-DO)': 'Protocatechuate 3,4-dioxygenase (pcaGH)', 'gentisate 1,2-DO': 'Gentisate 1,2-dioxygenase', 'nidAB (PAH DO)': 'PAH dioxygenase, mycobacterial type (nidAB)' };
  const ctr = el('div', { class: 'controls' }); famBox.append(ctr);
  const fsel = el('select', { id: 'famsel', 'aria-label': 'Enzyme family' });
  const og1 = el('optgroup', { label: 'KEGG ortholog / product families' }), og2 = el('optgroup', { label: 'CANT-HYD curated HMMs (trusted cutoffs)' });
  order.filter(k => P.families[k]).forEach(k => og1.append(new Option(`${NICE[k] || k}  ·  ${P.families[k].n_genes} genes`, k)));
  Object.keys(P.families).filter(k => k.startsWith('CANT-HYD') && P.families[k].n_genes >= 2).sort().forEach(k => og2.append(new Option(`${k.replace('CANT-HYD ', '')}  ·  ${P.families[k].n_genes} genes`, k)));
  fsel.append(og1, og2);
  ctr.append(el('label', { for: 'famsel' }, 'Enzyme family'), fsel);
  const fview = el('div', { class: 'grid2' }); famBox.append(fview);
  function drawFam() {
    const k = fsel.value, F = P.families[k]; fview.innerHTML = '';
    const a = el('div'), b = el('div'); fview.append(a, b);
    const box = fig(a, NICE[k] || k, 'Read pairs per million protein-coding read pairs (CPM), averaged over replicate samples per day. Log scale.');
    legendPools(box.parentNode.insertBefore(el('div'), box));
    timecourse(box, groupSeries(F.cpm), { ylab: 'CPM', unit: 'CPM', h: 170, floor: 0.1 });
    const hb = fig(b, 'Who expresses it', 'Share of the family’s reads by host MAG (or Kaiju genus for unbinned contigs), pooled over all samples in each group. “(Kaiju)” marks bins without a GTDB-Tk call, labelled by the majority genus of their contigs.');
    const groups = ['2024-25 oiled', '2024-25 control', '2025-26 oiled', '2025-26 control'];
    const W = Math.min(width(hb), 540);
    groups.forEach(gk => { const h = F.hosts[gk]; if (!h) return;
      const d = el('div', { style: 'margin:6px 0 12px' }); hb.append(d);
      d.append(el('div', { class: 'mono', style: 'color:var(--ink-2);margin-bottom:4px' }, `${gk} · ${fmt(h.pairs)} read pairs`));
      const sv = d3.select(d).append('svg').attr('width', W).attr('height', h.top.length * 19 + 2);
      const x = d3.scaleLinear().domain([0, 100]).range([0, W - 250]);
      h.top.forEach((t, i) => { const g = sv.append('g').attr('transform', `translate(0,${i * 19})`);
        g.append('rect').attr('x', 0).attr('y', 3).attr('height', 12).attr('rx', 3).attr('width', Math.max(1, x(t.pct))).attr('fill', css(gk.includes('oiled') ? '--exp' : '--ctl'));
        g.append('text').attr('x', x(t.pct) + 6).attr('y', 13).text(`${t.pct}%  ${t.label}${t.host && !t.host.startsWith('kaiju') && t.host !== 'unbinned' ? ' · ' + t.host : ''}`); }); });
  }
  fsel.onchange = drawFam; drawFam();
  // ---- heatmap
  const hbx = $('#heat', sec);
  hbx.append(el('h3', {}, 'Which genomes carry which pathways'));
  hbx.append(el('p', { class: 'kicker' }, 'The 47 MAGs with the highest degradation-pathway expression in the diesel pool, plus all unbinned contigs pooled. Circle size shows how complete the pathway is in the genome; colour shows the log₂ ratio of its expression in diesel-pool versus control-pool samples (both seasons, seawater and ice). Grey means encoded but rarely expressed. Hover for numbers; click a row to open that MAG in the community explorer.'));
  const fams = ['alkane', 'upper aromatic', 'lower aromatic', 'anaerobic', 'C1'];
  const cols = Object.keys(P.defs).filter(p => P.heat.some(h => h.p === p && h.c >= 0.5)).sort((a, b) => fams.indexOf(P.defs[a].family) - fams.indexOf(P.defs[b].family) || d3.ascending(P.defs[a].name, P.defs[b].name));
  const box = fig(hbx, '', '');
  const cell = 19, left = 250, top = 200, W = left + cols.length * cell + 150, H = top + P.rows.length * cell + 10;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const div = v => { const t = Math.max(-1, Math.min(1, v / 3)); return t < 0 ? d3.interpolateRgb(css('--other'), css('--ctl'))(-t) : d3.interpolateRgb(css('--other'), css('--exp'))(t); };
  cols.forEach((p, j) => { const x = left + j * cell + cell / 2;
    svg.append('text').attr('transform', `translate(${x + 3},${top - 8}) rotate(-58)`).text((P.defs[p].name.split(',')[0]).slice(0, 34)).attr('font-size', 10.5)
      .on('mousemove', ev => showTip(ev, `<b>${p}</b> ${P.defs[p].name}<br>${P.defs[p].family}<br><span class="mono">${P.defs[p].def}</span>`)).on('mouseleave', hideTip); });
  let fx = left; fams.forEach(f => { const n = cols.filter(p => P.defs[p].family === f).length; if (!n) return;
    svg.append('line').attr('x1', fx + 2).attr('x2', fx + n * cell - 2).attr('y1', top - 2).attr('y2', top - 2).attr('stroke', css('--ink-3')).attr('stroke-width', 2);
    svg.append('text').attr('x', fx + 2).attr('y', 12).attr('class', 'ttl').text(f).attr('font-size', 11).attr('transform', null); fx += n * cell; });
  const hm = {}; P.heat.forEach(h => hm[h.bin + '|' + h.p] = h);
  P.rows.forEach((r, i) => { const y = top + i * cell + cell / 2;
    const t = svg.append('text').attr('x', left - 8).attr('y', y + 4).attr('text-anchor', 'end').text(r.bin === 'unbinned' ? 'all unbinned contigs (pooled)' : `${r.label} · ${r.bin}`).style('cursor', 'pointer');
    t.on('click', () => { location.hash = '#community'; setTimeout(() => { const q = $('#magq'); if (q) { q.value = r.bin; q.dispatchEvent(new Event('input')); const tr = document.querySelector(`#magtab tr[data-b="${r.bin}"]`); if (tr) tr.click(); } }, 60); });
    svg.append('line').attr('x1', left).attr('x2', left + cols.length * cell).attr('y1', y).attr('y2', y).attr('stroke', css('--grid'));
    cols.forEach((p, j) => { const h = hm[r.bin + '|' + p]; if (!h || h.c <= 0) return;
      const rad = 2 + 6.5 * Math.sqrt(h.c); const expressed = (h.oil + h.ctl) > 1;
      svg.append('circle').attr('cx', left + j * cell + cell / 2).attr('cy', y).attr('r', rad).attr('fill', expressed ? div(h.lfc) : css('--rule')).attr('stroke', css('--panel'))
        .on('mousemove', ev => showTip(ev, `<b>${r.label} · ${r.bin}</b><br>${P.defs[p].name}<br>completeness ${Math.round(h.c * 100)}%<br>diesel ${fmt1(h.oil)} vs control ${fmt1(h.ctl)} CPM<br>log₂ ratio ${h.lfc}`)).on('mouseleave', hideTip); }); });
  const lg = el('div', { class: 'legend' }, `<span><i style="background:var(--ctl);height:10px"></i>higher in control</span><span><i style="background:var(--other);height:10px"></i>no change</span><span><i style="background:var(--exp);height:10px"></i>higher with diesel (full colour at log₂ ratio ≥ 3)</span><span><i style="background:var(--rule);height:10px"></i>encoded, &lt;1 CPM</span><span>circle area: module completeness</span>`);
  box.parentNode.insertBefore(lg, box);
}
RENDER.pathways = pathways;

const CATCOL = { 'RHO alpha (PAH/aromatic dioxygenase large)': '--s2', 'RHO beta (small subunit)': '--s4', 'RHO electron transfer (ferredoxin/reductase)': '--s5',
  'extradiol ring cleavage': '--s8', 'intradiol / ortho cleavage': '--s7', 'meta-cleavage lower pathway': '--s3', 'phenol/toluene monooxygenase (multicomponent)': '--s1',
  'salicylate / gentisate': '--s6', 'dihydrodiol dehydrogenase': '--seq3', 'alkane hydroxylation': '--seq4', 'Cu-membrane monooxygenase (pmo/amo-like)': '--accent', chemotaxis: '--ink-3' };
function cyclo(sec) {
  const C = D.cyclo; if (!C) { sec.innerHTML = '<h2>Cycloclasticus</h2><p class="kicker">Data not loaded.</p>'; return; }
  const short = b => b.replace('dastool-', '');
  sec.innerHTML = `<h2>Cycloclasticus</h2>
  <p class="kicker">Cycloclasticus is a marine specialist in polycyclic aromatic hydrocarbons. Five bins in the coassembly belong to it, forming three species-level groups. None is close to a cultured strain. Their nearest relatives are genomes assembled from Central Arctic sea ice, Canadian oil experiments and Canadian coastal water. A different species bloomed in each winter.</p>
  <div id="cy-stats"></div><div id="cy-ani"></div><div id="cy-inv"></div><div id="cy-tc"></div><div id="cy-op"></div><div id="cy-rd"></div><div id="cy-hom"></div>`;
  // genome table
  const st = $('#cy-stats', sec);
  const tw = el('div', { class: 'tablewrap', style: 'max-height:none' }); st.append(tw);
  tw.innerHTML = `<table><thead><tr><th>MAG</th><th class="num">Size (Mbp)</th><th class="num">Contigs</th><th class="num">N50 (kb)</th><th class="num">GC %</th><th class="num">CheckM2 compl. / contam.</th><th>Closest genome (skani ANI)</th></tr></thead><tbody>${C.stats.map(r => {
    const best = (r.skani_r232 || [])[0]; const ref = best ? (C.references.refs.find(x => best.ref.includes(x.acc.replace(/\.\d+$/, ''))) || {}) : {};
    return `<tr><td class="mono">${short(r.bin)}</td><td class="num">${r.size_mb}</td><td class="num">${r.contigs}</td><td class="num">${r.n50_kb}</td><td class="num">${r.gc}</td><td class="num">${r.checkm2_completeness ?? r.scg_completeness + ' (SCG)'} / ${r.checkm2_contamination ?? '–'}</td><td>${best ? `${best.ani}% · ${ref.name || best.name}${ref.source ? ' — ' + ref.source : ''}` : (r.gtdb_r226 || '')}</td></tr>`; }).join('')}</tbody></table>`;
  // ANI heatmap
  const A = C.references.ani, ab = fig($('#cy-ani', sec), 'Average nucleotide identity', 'skani pairwise ANI among the CMO bins and reference genomes. Blank cells are below about 80% ANI (no detectable alignment). The cultured isolates form one cluster that none of the CMO genomes joins.');
  const nm = A.names.map(n => { const r = C.references.refs.find(x => n.includes(x.acc)); return r ? r.name.replace(' (MAG)', '').replace('Cycloclasticus', 'C.') : short(n); });
  const n = nm.length, cs = 22, lw = 210, W2 = lw + n * cs + 140, sv = d3.select(ab).append('svg').attr('width', W2).attr('height', lw - 40 + n * cs + 10);
  const col = d3.scaleSequential().domain([80, 100]).interpolator(t => d3.interpolateRgb(css('--seq0'), css('--seq4'))(t));
  nm.forEach((t, i) => { sv.append('text').attr('x', lw - 6).attr('y', lw - 40 + i * cs + 12).attr('text-anchor', 'end').text(t).attr('font-size', 10.5);
    sv.append('text').attr('transform', `translate(${lw + i * cs + 11},${lw - 46}) rotate(-60)`).text(t).attr('font-size', 10); });
  A.values.forEach((row, i) => row.forEach((v, j) => { if (!v) return;
    sv.append('rect').attr('x', lw + j * cs).attr('y', lw - 40 + i * cs).attr('width', cs - 1).attr('height', cs - 1).attr('fill', col(v))
      .on('mousemove', ev => showTip(ev, `${nm[i]} × ${nm[j]}<br>ANI ${v}%`)).on('mouseleave', hideTip); }));
  // inventory
  const I = C.hc_inventory, ib = fig($('#cy-inv', sec), 'Hydrocarbon gene inventory', 'Number of genes per category in each bin, from Bakta product names (ring-hydroxylating oxygenase = RHO). The four larger genomes carry 14–19 RHO alpha subunits each (5 in the fragmentary lorbin_722_sub). That large dioxygenase repertoire is the Cycloclasticus signature.');
  const iw = 290 + I.bins.length * 90, svi = d3.select(ib).append('svg').attr('width', iw).attr('height', 40 + I.categories.length * 22);
  I.bins.forEach((b, j) => svi.append('text').attr('x', 250 + j * 90 + 30).attr('y', 14).attr('text-anchor', 'middle').text(short(b).replace('_sub', '*')).attr('font-size', 10.5));
  const mx = d3.max(I.counts.flat());
  I.categories.forEach((c, i) => { svi.append('text').attr('x', 244).attr('y', 36 + i * 22).attr('text-anchor', 'end').text(c.length > 40 ? c.slice(0, 39) + '…' : c).attr('font-size', 10.5);
    I.counts[i].forEach((v, j) => { if (!v) return; svi.append('circle').attr('cx', 250 + j * 90 + 30).attr('cy', 32 + i * 22).attr('r', 3 + 7 * Math.sqrt(v / mx)).attr('fill', css(CATCOL[c] || '--accent'));
      svi.append('text').attr('x', 250 + j * 90 + 44).attr('y', 36 + i * 22).text(v).attr('font-size', 10); }); });
  // share of transcripts
  const tc = $('#cy-tc', sec); tc.append(el('h3', {}, 'Share of all mapped transcripts'));
  const chips = el('div', { class: 'controls' }); tc.append(chips);
  const tbox = el('div'); tc.append(tbox);
  const bins = ['metabat_9342', 'lorbin_184', 'metabat_5432', 'comebin_3787_sub', 'lorbin_722_sub'];
  const byId = Object.fromEntries(C.mag_timecourse.map(r => [r.sample, r]));
  let cur = 'metabat_9342';
  function drawTc() { tbox.innerHTML = ''; chips.querySelectorAll('button').forEach(b => b.setAttribute('aria-pressed', b.dataset.b === cur));
    const box = fig(tbox, `${cur}: % of mapped read pairs`, 'Read pairs on the bin’s protein-coding genes as a percentage of all protein-coding read pairs, averaged over replicate samples per day. metabat_9342 drives the 2024-25 bloom; lorbin_184 and its sibling sub-bins drive 2025-26.');
    legendPools(box.parentNode.insertBefore(el('div'), box));
    timecourse(box, groupSeries(D.samples.map(s => byId[s.id] ? byId[s.id][cur] : null)), { ylab: '% of mapped pairs', unit: '%', floor: 0.0005, h: 190, fmt: d3.format('.3~g') }); }
  bins.forEach(b => { const btn = el('button', { class: 'chip', 'data-b': b, 'aria-pressed': 'false' }, b); btn.onclick = () => { cur = b; drawTc(); }; chips.append(btn); });
  drawTc();
  // operon browser
  const ob = $('#cy-op', sec); ob.append(el('h3', {}, 'Operon browser'));
  ob.append(el('p', { class: 'kicker' }, 'Hydrocarbon gene clusters in the Cycloclasticus bins. Coloured arrows are hydrocarbon genes; grey arrows are their neighbours. Bars under each gene show mean expression (CPM) in diesel-pool versus control-pool samples for the chosen season. The chart below follows the whole cluster through time.'));
  const oc = el('div', { class: 'controls' }); ob.append(oc);
  const osel = el('select', { id: 'opsel', 'aria-label': 'Operon' });
  const ops = C.operons.slice().sort((a, b) => b.pairs_oiled_total - a.pairs_oiled_total);
  ops.forEach(o => osel.append(new Option(`${short(o.bin)} · ${o.contig}:${fmt(o.start)}–${fmt(o.end)} · ${o.cats.map(c => c.split(' (')[0]).join(', ')}`, o.id)));
  const ssel = el('select', { id: 'opseason', 'aria-label': 'Season' }); ['2024-25', '2025-26'].forEach(x => ssel.append(new Option('CMO ' + x, x)));
  oc.append(el('label', { for: 'opsel' }, 'Cluster'), osel, el('label', { for: 'opseason' }, 'Season'), ssel);
  const obox = el('div'); ob.append(obox);
  function drawOp() {
    obox.innerHTML = ''; const o = C.operons.find(x => x.id === osel.value), se = ssel.value;
    const box = fig(obox, `${short(o.bin)} · ${o.contig} · ${fmt(o.end - o.start)} bp${o.genomad_plasmid ? ' · geNomad calls this contig a plasmid' : ''}`, 'Hover a gene for its product, KEGG ortholog and expression. “island” marks genes inside an IslandPath genomic island; “mobile” marks a transposase or integrase within 10 kb.');
    const W = Math.max(width(box), 640), x = d3.scaleLinear().domain([o.start, o.end]).range([10, W - 10]);
    const g0 = o.genes; const mxc = d3.max(g0, g => Math.max(g.expr[se]?.cpm_oiled || 0, g.expr[se]?.cpm_ctl || 0)) || 1;
    const by = d3.scaleSqrt().domain([0, mxc]).range([0, 70]);
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', 190);
    svg.append('line').attr('x1', 10).attr('x2', W - 10).attr('y1', 40).attr('y2', 40).attr('stroke', css('--rule'));
    g0.forEach(g => { const x0 = x(g.start), x1 = x(g.end), w = Math.max(3, x1 - x0), hd = Math.min(10, w * .4), y = 28, h = 24;
      const pts = g.strand === '+' ? [[x0, y], [x0 + w - hd, y], [x0 + w, y + h / 2], [x0 + w - hd, y + h], [x0, y + h]] : [[x0 + w, y], [x0 + hd, y], [x0, y + h / 2], [x0 + hd, y + h], [x0 + w, y + h]];
      const e = g.expr[se] || {};
      svg.append('polygon').attr('points', pts.map(p => p.join(',')).join(' ')).attr('fill', g.cat ? css(CATCOL[g.cat] || '--accent') : css('--rule')).attr('stroke', css('--panel'))
        .on('mousemove', ev => showTip(ev, `<b>${g.product}</b>${g.gene ? ' (' + g.gene + ')' : ''}<br>${g.id}${g.ko ? ' · ' + g.ko : ''}<br>${g.cat || 'neighbour'}${g.island ? ' · island' : ''}${g.mob_near ? ' · mobile' : ''}<br>${se}: diesel ${fmt1(e.cpm_oiled || 0)} vs control ${fmt1(e.cpm_ctl || 0)} CPM`)).on('mouseleave', hideTip);
      const cx = x0 + w / 2;
      svg.append('rect').attr('x', cx - 5).attr('width', 4.5).attr('y', 64).attr('height', by(e.cpm_oiled || 0)).attr('fill', css('--exp'));
      svg.append('rect').attr('x', cx + .5).attr('width', 4.5).attr('y', 64).attr('height', by(e.cpm_ctl || 0)).attr('fill', css('--ctl'));
      if (g.island || g.mob_near) svg.append('text').attr('x', cx).attr('y', 22).attr('text-anchor', 'middle').text(g.island ? 'island' : 'mobile').attr('font-size', 9.5); });
    svg.append('text').attr('x', 10).attr('y', 160).text(`bars: mean CPM in ${se}, diesel (orange) and control (blue); tallest = ${fmt1(mxc)} CPM`);
    svg.append('text').attr('x', 10).attr('y', 12).text(`${fmt(o.start)}`); svg.append('text').attr('x', W - 10).attr('y', 12).attr('text-anchor', 'end').text(`${fmt(o.end)}`);
    const tb = fig(obox, 'Whole cluster through time', 'Sum of the cluster’s hydrocarbon genes, CPM.');
    legendPools(tb.parentNode.insertBefore(el('div'), tb));
    timecourse(tb, groupSeries(D.samples.map(s => o.cpm_timecourse[s.id] ?? null)), { ylab: 'CPM', unit: 'CPM', h: 170, floor: 0.1 });
  }
  osel.onchange = drawOp; ssel.onchange = drawOp; ssel.value = '2025-26'; drawOp();
  // RNA:DNA
  const rd = fig($('#cy-rd', sec), 'Transcriptional activity per genome copy', 'RNA:DNA ratio is the % of mapped transcripts divided by the % of mapped nanopore DNA, for samples with both. A ratio above 1 means the genome is transcribed out of proportion to its abundance. metabat_9342 reached 8–85× in three of the four day 7–8 samples of 2024-25 while still only 0.01–0.03% of the DNA.');
  const pts = C.dna.filter(r => r.rna_dna != null && r.rna_dna > 0);
  const W = Math.min(width(rd), 760), H = 300, mL = 50, sv2 = d3.select(rd).append('svg').attr('width', W).attr('height', H);
  const xd = d3.scalePoint().domain(['2024-25', '2025-26']).range([mL + 120, W - 120]);
  const yd = d3.scaleLog().domain([d3.min(pts, r => r.rna_dna) * 0.8, d3.max(pts, r => r.rna_dna) * 1.2]).range([H - 30, 10]);
  sv2.append('g').attr('class', 'axis').attr('transform', `translate(${mL},0)`).call(d3.axisLeft(yd).ticks(5, '~g'));
  sv2.append('line').attr('x1', mL).attr('x2', W).attr('y1', yd(1)).attr('y2', yd(1)).attr('stroke', css('--ink-3')).attr('stroke-dasharray', '3 3');
  const bc = { 'dastool-metabat_9342': '--s2', 'dastool-lorbin_184': '--s1', 'dastool-metabat_5432': '--s3', 'dastool-comebin_3787_sub': '--s7', 'dastool-lorbin_722_sub': '--s5' };
  const leg = el('div', { class: 'legend' }, Object.entries(bc).map(([b, c]) => `<span><i style="background:var(${c});height:8px;width:8px;border-radius:50%"></i>${short(b)}</span>`).join('') + '<span>left of each season: control pool · right: diesel pool</span>');
  rd.parentNode.insertBefore(leg, rd);
  pts.forEach((r, i) => { const pool = /_EXP_/.test(r.sample) ? 1 : /_CTL_/.test(r.sample) ? -1 : 0;
    sv2.append('circle').attr('cx', xd(r.season) + pool * 55 + ((i * 37) % 30 - 15)).attr('cy', yd(r.rna_dna)).attr('r', 4.5).attr('fill', css(bc[r.bin] || '--other')).attr('stroke', css('--panel'))
      .on('mousemove', ev => showTip(ev, `<b>${r.sample}</b><br>${short(r.bin)}<br>RNA ${d3.format('.3~g')(r.rna_pct)}% · DNA ${d3.format('.3~g')(r.dna_pct)}%<br>ratio ${d3.format('.3~g')(r.rna_dna)}`)).on('mouseleave', hideTip); });
  ['2024-25', '2025-26'].forEach(se => sv2.append('text').attr('class', 'ttl').attr('x', xd(se)).attr('y', H - 8).attr('text-anchor', 'middle').text(`CMO ${se}`));
  sv2.append('text').attr('x', mL + 4).attr('y', yd(1) - 4).text('RNA:DNA = 1');
  // homology to references
  const hb = $('#cy-hom', sec); hb.append(el('h3', {}, 'Where the hydrocarbon genes come from'));
  hb.append(el('p', { class: 'kicker' }, 'Amino-acid identity of each CMO Cycloclasticus hydrocarbon gene to its best match in each reference genome. Genes inherited along with the genome sit near the genome-wide average amino-acid identity (72–82% to the isolates). Genes far above that average point to recent exchange. The dmp/tmo phenol–toluene cassette is 78–98% identical to C. zancles (median about 92%), yet sits on the chromosome with no mobility genes.'));
  const H2 = C.references.hc_homology, refs = Object.keys(H2[0]).filter(k => !['gene_id', 'bin', 'contig', 'product'].includes(k));
  const hw = el('div', { class: 'tablewrap' }); hb.append(hw);
  const cc = v => v == null ? '' : `background:${d3.interpolateRgb(css('--panel'), css('--seq3'))(Math.max(0, (v - 60) / 40))};color:${v > 88 ? '#fff' : 'inherit'}`;
  hw.innerHTML = `<table><thead><tr><th>Gene</th><th>Bin</th><th>Product</th>${refs.map(r => `<th class="num">${r.replace(/\(.*?\)/, '').trim()}</th>`).join('')}</tr></thead><tbody>${H2.slice().sort((a, b) => d3.ascending(a.bin, b.bin) || d3.ascending(a.product, b.product)).map(r => `<tr><td class="mono">${r.gene_id}</td><td class="mono">${short(r.bin)}</td><td>${r.product}</td>${refs.map(k => `<td class="num" style="${cc(r[k])}">${r[k] == null ? '–' : r[k].toFixed(0)}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
}
RENDER.cyclo = cyclo;

const ROLECOL = { conjugation: '--s7', 'IS/integrase': '--s4', 'mobility': '--s5', replication: '--s3' };
function courseSeries(course) { return Object.entries(course).map(([k, pts]) => { const [season, pool] = k.split(' '); return { season, pool, medium: 'SW', pts: pts.map(([d, v]) => ({ day: d, v })) }; }); }
function hgt(sec) {
  const H = D.hgt; if (!H) { sec.innerHTML = '<h2>Gene transfer</h2><p class="kicker">Data not loaded.</p>'; return; }
  const cnt = H.counts.confident, F = H.sharing_funnel;
  sec.innerHTML = `<h2>Gene transfer and diesel degradation</h2>
  <p class="kicker">Catabolic genes on plasmids and transposons can move between bacteria and spread a new capacity through a community. We asked three things: how many hydrocarbon genes sit on mobile DNA, whether the diesel response relied on them, and which genes have recently crossed between unrelated hosts. Gene sets here use Bakta products and curated CANT-HYD profiles (${fmt(cnt.genes)} “confident” hydrocarbon genes on ${fmt(cnt.contigs)} contigs).</p>
  <div class="stat-row" id="h-stats"></div><div id="h-ind"></div><div id="h-cat"></div><div id="h-pl"></div><div id="h-net"></div><div id="h-mg"></div>`;
  [[`${(100 * cnt.mge_borne / cnt.genes).toFixed(1)}%`, 'of hydrocarbon genes on a plasmid, phage, island or integron'], [`${(100 * cnt.mobile_assoc / cnt.genes).toFixed(1)}%`, 'on mobile DNA or within 10 kb of a transposase'],
   [`${H.mobile_vs_chrom.find(r => r.tier === 'confident' && r.season === '2024-25' && r.flag === 'mobile_assoc').frac_of_induction}% / ${H.mobile_vs_chrom.find(r => r.tier === 'confident' && r.season === '2025-26' && r.flag === 'mobile_assoc').frac_of_induction}%`, 'of the diesel-induced increase carried by mobile-associated genes (2024-25 / 2025-26)'],
   [F.flank_confirmed_HGT, 'flank-confirmed recent transfers between genera']].forEach(([v, l]) => $('#h-stats', sec).append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  // induction by genus
  const ib = fig($('#h-ind', sec), 'Who carries the induced transcription, and how much of it is mobile', 'Share of the total diesel-induced increase in hydrocarbon-gene transcripts (diesel minus control, CPM) by host genus. The darker segment of each bar is the part carried by mobile-associated copies. Cycloclasticus does nearly all of its work from the chromosome; the Sphingorhabdus contribution is mostly (52–73%) plasmid-borne. Pelagibacter (SAR11) appears through Rieske oxygenases that are equally high before dosing, so it reflects pool composition rather than a diesel response.');
  ['2024-25', '2025-26'].forEach(se => { const rows = H.induction_by_genus.filter(r => r.season === se && r.pct_of_induction > 0.8).sort((a, b) => b.pct_of_induction - a.pct_of_induction).slice(0, 9);
    const W = Math.min(width(ib), 760), lw = 190, rh = 22, sv = d3.select(ib).append('svg').attr('width', W).attr('height', 26 + rows.length * rh);
    sv.append('text').attr('class', 'ttl').attr('x', 0).attr('y', 14).text(`CMO ${se}`);
    const x = d3.scaleLinear().domain([0, 90]).range([0, W - lw - 120]);
    rows.forEach((r, i) => { const y = 22 + i * rh;
      sv.append('text').attr('x', lw - 8).attr('y', y + 12).attr('text-anchor', 'end').text(r.genus);
      sv.append('rect').attr('x', lw).attr('y', y + 2).attr('height', 14).attr('rx', 3).attr('width', x(r.pct_of_induction)).attr('fill', css('--exp')).attr('opacity', .45);
      sv.append('rect').attr('x', lw).attr('y', y + 2).attr('height', 14).attr('rx', 3).attr('width', x(r.pct_of_induction * r.mobile_pct_within / 100)).attr('fill', css('--exp'));
      sv.append('text').attr('x', lw + x(r.pct_of_induction) + 6).attr('y', y + 13).text(`${r.pct_of_induction}% · ${r.mobile_pct_within}% of it mobile`); }); });
  // by category
  const cb = fig($('#h-cat', sec), 'Mobile share by enzyme family', 'For each family: diesel-induced transcripts (CPM, log scale) and the percentage carried by mobile-associated copies. CYP153, alkB (mainly in 2024-25: 21% versus 3%), the PAH lower-pathway genes (81% mobile in 2025-26) and the rarer ortho-cleavage, gentisate and short-alkane genes are the most mobile. The dominant aromatic oxygenases are almost entirely chromosomal.');
  const cats = [...new Set(H.mobile_share_by_category.map(r => r.category))];
  const W2 = Math.min(Math.max(width(cb), 700), 900), lw2 = 260, rh2 = 24, sv2 = d3.select(cb).append('svg').attr('width', W2).attr('height', 40 + cats.length * rh2);
  const xs = d3.scaleLinear().domain([0, 100]).range([lw2, W2 - 30]);
  sv2.append('g').attr('class', 'axis').attr('transform', `translate(0,${cats.length * rh2 + 10})`).call(d3.axisBottom(xs).ticks(5).tickFormat(v => v + '% mobile'));
  const rr = d3.scaleSqrt().domain([0, d3.max(H.mobile_share_by_category, r => r.induced_cpm)]).range([2, 14]);
  cats.forEach((c, i) => { const y = i * rh2 + 12; sv2.append('text').attr('x', lw2 - 10).attr('y', y + 4).attr('text-anchor', 'end').text(c);
    sv2.append('line').attr('x1', lw2).attr('x2', W2 - 30).attr('y1', y).attr('y2', y).attr('stroke', css('--grid'));
    H.mobile_share_by_category.filter(r => r.category === c && r.induced_cpm > 0).forEach(r => sv2.append('circle').attr('cx', xs(r.mobile_share ?? 0)).attr('cy', y).attr('r', rr(r.induced_cpm))
      .attr('fill', css(r.season === '2024-25' ? '--s7' : '--s3')).attr('fill-opacity', .75).attr('stroke', css('--panel'))
      .on('mousemove', ev => showTip(ev, `<b>${c}</b> ${r.season}<br>induced ${fmt1(r.induced_cpm)} CPM<br>${r.mobile_share}% mobile<br>log₂ diesel/control ${r.l2fc_total}`)).on('mouseleave', hideTip)); });
  cb.parentNode.insertBefore(el('div', { class: 'legend' }, `<span><i style="background:var(--s7);height:8px;width:8px;border-radius:50%"></i>2024-25</span><span><i style="background:var(--s3);height:8px;width:8px;border-radius:50%"></i>2025-26</span><span>circle area: induced CPM</span>`), cb);
  // plasmid browser
  const pb = $('#h-pl', sec); pb.append(el('h3', {}, 'Catabolic plasmids and mobile elements'));
  pb.append(el('p', { class: 'kicker' }, 'Contigs that carry hydrocarbon genes and look mobile. The Sphingorhabdus conjugative plasmid contig_0009679 carries an alkane cassette (CYP153 and alkB) flanked by IS3 transposases. Its transcripts jump 20–60× within two days of diesel, and its DNA depth rises in most diesel-pool runs at days 7–10.'));
  const pc = el('div', { class: 'controls' }); pb.append(pc);
  const psel = el('select', { id: 'plsel', 'aria-label': 'Mobile element' });
  const mges = H.catabolic_mges.filter(m => H.catabolic_mge_courses[m.contig]).sort((a, b) => (+b.max_cpm || 0) - (+a.max_cpm || 0));
  mges.forEach(m => psel.append(new Option(`${m.contig} · ${m.host_genus || 'unknown host'} · ${(+m.length / 1000).toFixed(0)} kb · ${m.hc_categories}`, m.contig)));
  pc.append(el('label', { for: 'plsel' }, 'Element'), psel);
  const pv = el('div'); pb.append(pv);
  function drawPl() { pv.innerHTML = ''; const m = H.catabolic_mges.find(x => x.contig === psel.value);
    const box = fig(pv, `${m.contig} · ${m.kind} · ${(+m.length / 1000).toFixed(1)} kb · host ${m.host_genus || '?'} (${m.host_tax_source || ''})`, `${m.mobility ? 'Mobility: ' + m.mobility + '. ' : ''}${m.n_hc} hydrocarbon gene(s): ${m.hc_products}. Diesel versus control: ${fmt1(+m['cpm_oiled_2024-25'] || 0)} vs ${fmt1(+m['cpm_ctl_2024-25'] || 0)} CPM in 2024-25; ${fmt1(+m['cpm_oiled_2025-26'] || 0)} vs ${fmt1(+m['cpm_ctl_2025-26'] || 0)} CPM in 2025-26.`);
    const genes = H.plasmid_maps[m.contig];
    if (genes) { const W = Math.max(width(box), 700), L = +m.length, x = d3.scaleLinear().domain([0, L]).range([10, W - 10]);
      const roles = [...new Set(genes.map(g => g.role).filter(Boolean))];
      box.append(el('div', { class: 'legend' }, roles.map(r => `<span><i style="background:var(${r.startsWith('HC') ? '--exp' : ROLECOL[r] || '--accent'});height:10px"></i>${r}</span>`).join('') + '<span><i style="background:var(--rule);height:10px"></i>other</span>'));
      const svg = d3.select(box).append('svg').attr('width', W).attr('height', 70);
      svg.append('line').attr('x1', 10).attr('x2', W - 10).attr('y1', 35).attr('y2', 35).attr('stroke', css('--rule'));
      genes.forEach(g => { const x0 = x(g.start), x1 = x(g.end), w = Math.max(2, x1 - x0), hd = Math.min(7, w * .4), y = 25, h = 20;
        const pts = g.strand === '+' ? [[x0, y], [x0 + w - hd, y], [x0 + w, y + h / 2], [x0 + w - hd, y + h], [x0, y + h]] : [[x0 + w, y], [x0 + hd, y], [x0, y + h / 2], [x0 + hd, y + h], [x0 + w, y + h]];
        svg.append('polygon').attr('points', pts.map(p => p.join(',')).join(' ')).attr('fill', g.role ? css(g.role.startsWith('HC') ? '--exp' : ROLECOL[g.role] || '--accent') : css('--rule')).attr('stroke', css('--panel'))
          .on('mousemove', ev => showTip(ev, `<b>${g.product}</b><br>${g.gene_id} · ${fmt(g.start)}–${fmt(g.end)} (${g.strand})${g.role ? '<br>' + g.role : ''}`)).on('mouseleave', hideTip); });
      svg.append('text').attr('x', 10).attr('y', 64).text('0'); svg.append('text').attr('x', W - 10).attr('y', 64).attr('text-anchor', 'end').text(fmt(L) + ' bp'); }
    const C2 = H.catabolic_mge_courses[m.contig];
    const g2 = el('div', { class: 'grid2' }); pv.append(g2); const a = el('div'), b = el('div'); g2.append(a, b);
    const ra = fig(a, 'Hydrocarbon-gene transcripts on this element', 'CPM, seawater and ice averaged per day.'); legendPools(ra.parentNode.insertBefore(el('div'), ra), false);
    timecourse(ra, courseSeries(C2.rna), { ylab: 'CPM', unit: 'CPM', h: 160, floor: 0.1 });
    if (C2.dna && Object.keys(C2.dna).length) { const rb = fig(b, 'DNA depth of the element', 'Mean read depth in the nanopore metagenomes (nanopore DNA sampling ended around day 11).'); legendPools(rb.parentNode.insertBefore(el('div'), rb), false);
      timecourse(rb, courseSeries(C2.dna).map(s => ({ ...s, points_only: true })), { ylab: 'depth (×)', unit: 'depth', h: 160, floor: 0.01, fmt: d3.format('.2~f') }); }
  }
  psel.onchange = drawPl; if (mges.some(m => m.contig === 'contig_0009679')) psel.value = 'contig_0009679'; drawPl();
  // network
  const nb = $('#h-net', sec); nb.append(el('h3', {}, 'Recent transfers between genera'));
  nb.append(el('p', { class: 'kicker' }, `Near-identical hydrocarbon genes in different genera narrow quickly under scrutiny: ${fmt(F.pairs_95aa_diff_contigs)} protein pairs at ≥95% identity, ${F.nt_ge_99} at ≥99% nucleotide identity, ${F.cross_genus_nt99} between genera whose DNA abundance profiles are independent, and ${F.flank_confirmed_HGT} where an identical block sits inside otherwise unrelated surrounding DNA. Solid lines are those flank-confirmed transfers; dashed lines are unresolved (short contigs). Node colour is the host class.`));
  const box = fig(nb, '', 'Drag nodes to untangle. Hover an edge for the gene category, identity and flank evidence.');
  const E2 = H.hgt_network.edges.filter(e => !e.flank_call.startsWith('same lineage'));
  const ids = new Set(E2.flatMap(e => [e.source, e.target]));
  const N = H.hgt_network.nodes.filter(n => ids.has(n.id)).map(n => ({ ...n }));
  const classes = [...new Set(N.map(n => n.class || 'unknown'))]; const cc = ['--s1', '--s2', '--s3', '--s4', '--s5', '--s6', '--s7', '--s8'];
  const ccol = c => { const i = classes.indexOf(c || 'unknown'); return i < 8 ? css(cc[i]) : css('--other'); };
  box.parentNode.insertBefore(el('div', { class: 'legend' }, classes.slice(0, 8).map((c, i) => `<span><i style="background:var(${cc[i]});height:8px;width:8px;border-radius:50%"></i>${c}</span>`).join('') + (classes.length > 8 ? '<span><i style="background:var(--other);height:8px;width:8px;border-radius:50%"></i>other classes</span>' : '')), box);
  const W = Math.max(width(box), 640), Hh = 520, svg = d3.select(box).append('svg').attr('width', W).attr('height', Hh);
  const L = E2.map(e => ({ ...e }));
  const sim = d3.forceSimulation(N).force('link', d3.forceLink(L).id(d => d.id).distance(70)).force('charge', d3.forceManyBody().strength(-160)).force('center', d3.forceCenter(W / 2, Hh / 2)).force('x', d3.forceX(W / 2).strength(.05)).force('y', d3.forceY(Hh / 2).strength(.08));
  const conf = e => e.flank_call.startsWith('HGT');
  const ln = svg.append('g').selectAll('line').data(L).join('line').attr('stroke', e => conf(e) ? css('--exp') : css('--ink-3')).attr('stroke-width', e => conf(e) ? 2.4 : 1.2).attr('stroke-dasharray', e => conf(e) ? null : '4 3')
    .on('mousemove', (ev, e) => showTip(ev, `<b>${e.category}</b><br>${e.source.id} ↔ ${e.target.id} (${e.level})<br>${(e.aa_ident * 100).toFixed(1)}% aa · ${e.nt_ident.toFixed(1)}% nt<br>${e.flank_call}<br>${e.q_contig} / ${e.t_contig}`)).on('mouseleave', hideTip);
  const nd = svg.append('g').selectAll('g').data(N).join('g').call(d3.drag().on('start', (ev, d) => { if (!ev.active) sim.alphaTarget(.3).restart(); d.fx = d.x; d.fy = d.y; }).on('drag', (ev, d) => { d.fx = ev.x; d.fy = ev.y; }).on('end', (ev, d) => { if (!ev.active) sim.alphaTarget(0); d.fx = null; d.fy = null; }));
  nd.append('circle').attr('r', 7).attr('fill', d => ccol(d.class)).attr('stroke', css('--panel')).attr('stroke-width', 1.5);
  nd.append('text').attr('x', 10).attr('y', 4).text(d => d.id).attr('font-size', 10.5);
  nd.on('mousemove', (ev, d) => showTip(ev, `<b>${d.id}</b><br>${d.family || ''} · ${d.class || ''}<br>taxonomy: ${d.tax_source}`)).on('mouseleave', hideTip);
  sim.on('tick', () => { N.forEach(d => { d.x = Math.max(12, Math.min(W - 90, d.x)); d.y = Math.max(12, Math.min(Hh - 12, d.y)); });
    ln.attr('x1', e => e.source.x).attr('y1', e => e.source.y).attr('x2', e => e.target.x).attr('y2', e => e.target.y); nd.attr('transform', d => `translate(${d.x},${d.y})`); });
  // flank-confirmed table
  const tf = el('div', { class: 'tablewrap', style: 'margin-top:12px' }); nb.append(tf);
  const FL = H.hgt_flanks.filter(f => (f.interpretation || '').startsWith('HGT'));
  tf.innerHTML = `<table><thead><tr><th>Category</th><th>Level</th><th>Host A</th><th>Host B</th><th class="num">nt identity</th><th class="num">Shared block</th><th class="num">Genes in block</th><th>Evidence</th></tr></thead><tbody>${FL.map(f => `<tr><td>${f.category}</td><td>${f.level}</td><td>${String(f.q_host).replace(/nan/g, 'unknown')}</td><td>${String(f.t_host).replace(/nan/g, 'unknown')}</td><td class="num">${(+f.nt_ident).toFixed(1)}%</td><td class="num">${fmt(f.shared_block_bp)} bp</td><td class="num">${f.block_genes ?? '–'}</td><td>${f.interpretation}</td></tr>`).join('')}</tbody></table>`;
}
RENDER.hgt = hgt;


async function boot() {
  tip = $('#tip');
  D = await (await fetch('data.json', { cache: 'no-cache' })).json();
  try { const r = await fetch('lineage.json', { cache: 'no-cache' }); if (r.ok) D.lineage = await r.json(); } catch (e) { /* optional */ }
  for (const f of ['cyclo', 'hgt', 'ice', 'islands', 'sulfur', 'metals', 'nitrogen', 'light', 'porticocc', 'oxygen', 'euks', 'states', 'guilds', 'viruses']) { try { const r = await fetch(f + '.json', { cache: 'no-cache' }); if (r.ok) D[f] = await r.json(); } catch (e) { /* optional */ } }
  const nav = $('#nav'), main = $('#main'); main.innerHTML = '';
  const label = id => SECTIONS.find(s => s[0] === id)[1];
  NAV.forEach(([a, ids]) => {
    if (!ids) { nav.append(el('a', { href: '#' + a }, label(a))); return; }
    const d = el('details', { class: 'dd' }); d.append(el('summary', {}, a));
    const menu = el('div', { class: 'menu' }); ids.forEach(id => menu.append(el('a', { href: '#' + id }, label(id)))); d.append(menu); nav.append(d);
    d.addEventListener('toggle', () => { if (d.open) nav.querySelectorAll('details.dd').forEach(o => { if (o !== d) o.open = false; }); });
  });
  document.addEventListener('click', ev => { if (!ev.target.closest('#nav')) nav.querySelectorAll('details.dd').forEach(o => o.open = false); });
  document.addEventListener('keydown', ev => { if (ev.key === 'Escape') nav.querySelectorAll('details.dd[open]').forEach(o => { o.open = false; o.querySelector('summary').focus(); }); });
  SECTIONS.forEach(([id]) => { const s = el('section', { class: 'view', id: 'v-' + id }); s.hidden = true; main.append(s); });
  const done = {};
  function go() {
    const id = (location.hash || '#overview').slice(1); const ok = SECTIONS.some(s => s[0] === id) ? id : 'overview';
    SECTIONS.forEach(([k]) => { $('#v-' + k).hidden = k !== ok; const a = nav.querySelector(`a[href="#${k}"]`); if (a) a.setAttribute('aria-current', k === ok); });
    nav.querySelectorAll('details.dd').forEach(d => { d.open = false; d.querySelector('summary').setAttribute('data-current', !!d.querySelector(`a[href="#${ok}"]`)); });
    if (!done[ok]) { done[ok] = 1; const fn = RENDER[ok]; const sec = $('#v-' + ok);
      if (fn) fn(sec); else sec.innerHTML = `<h2>${SECTIONS.find(s => s[0] === ok)[1]}</h2><p class="kicker">This section is still being assembled.</p>`; }
    try { localStorage.setItem('cmo-sec', ok); } catch (e) { /* storage blocked */ }
  }
  // theme picker: Auto follows the browser; Light/Dark set data-theme on <html> (tokens in index.html). Charts read colours
  // when drawn, so a change clears the rendered sections and redraws the visible one, keeping the scroll position.
  const THEMES = [['auto', 'Auto'], ['light', 'Light'], ['dark', 'Dark']];
  const tp = el('div', { class: 'theme', role: 'radiogroup', 'aria-label': 'Colour theme' });
  let cur = 'auto'; try { cur = localStorage.getItem('cmo-theme') || 'auto'; } catch (e) { /* storage blocked */ }
  const redraw = () => { const y = scrollY; for (const k in done) { delete done[k]; $('#v-' + k).innerHTML = ''; } go(); scrollTo({ top: y }); };
  const setTheme = (t, draw = true) => {
    if (t === 'auto') document.documentElement.removeAttribute('data-theme'); else document.documentElement.setAttribute('data-theme', t);
    tp.querySelectorAll('button').forEach(b => b.setAttribute('aria-checked', b.dataset.t === t));
    try { localStorage.setItem('cmo-theme', t); } catch (e) { /* storage blocked */ }
    cur = t; if (draw) redraw();
  };
  THEMES.forEach(([t, l]) => { const b = el('button', { type: 'button', role: 'radio', 'data-t': t, 'aria-checked': t === cur }, l); b.onclick = () => { if (t !== cur) setTheme(t); }; tp.append(b); });
  nav.append(tp);
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { if (cur === 'auto') redraw(); });
  setTheme(cur, false);
  addEventListener('hashchange', () => { go(); scrollTo({ top: 0 }); });
  if (!location.hash) { try { const k = localStorage.getItem('cmo-sec'); if (k) history.replaceState(null, '', '#' + k); } catch (e) { /* storage blocked */ } }
  go();
  $('#foot').innerHTML = D.footer || '';
}
boot();
