/* Physiological states: PCA of each MAG's own transcriptome (within-genome profiles, so abundance drops out) and state shifts with diesel.
   Data: states.json (states/scripts/states_site.py). Uses helpers from app.js (el, fig, width, css, showTip, hideTip, fmt, fmt1, timecourse, legendPools). */
const stSgn = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.2f')(v);
const stPct = d3.format('.0%');
const stQ = q => q == null ? '' : q < 0.001 ? 'q < 0.001' : 'q = ' + d3.format('.2g')(q);
function stSeries(J, values) {
  const out = {};
  J.samples.forEach((s, i) => { if (values[i] == null) return;
    const k = `${s.season}|${s.pool}|${s.medium}`; ((out[k] ||= {})[s.day] ||= []).push([values[i], s.name]); });
  return Object.entries(out).map(([k, byDay]) => { const [season, pool, medium] = k.split('|');
    return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a, x => x[0]), name: a.map(x => x[1]).join(', ') })) }; });
}
const stPhaseCol = s => s.phase === 'pre' ? css('--other') : s.pool === 'EXP' ? css('--exp') : css('--ctl');

function statesSec(sec) {
  const J = D.states;
  if (!J) { sec.innerHTML = '<h2>Physiological states</h2><p class="kicker">states.json is missing.</p>'; return; }
  const M = J.mags, sig = M.filter(m => m.q_25 != null && m.q_25 < 0.05), tested = M.filter(m => m.q_25 != null);
  const top = M.filter(m => m.detail)[0];
  sec.innerHTML = `<h2>Physiological states</h2>
  <p class="kicker">The other pages ask who is there and which genes the community expresses. This page asks how each genome behaves. For every MAG with enough reads, the expression profile is that genome’s own transcripts as shares of its own total. A bloom makes a genome more abundant but does not change these shares, so abundance drops out. A principal component analysis (PCA) of these profiles places every sample on a map of the genome’s physiological state. Samples close together are in the same state; a shift with diesel means the cells themselves changed what they were doing. ${J.n_pass} MAGs had at least ${fmt(J.min_reads)} protein-coding reads in at least ${J.min_samples} samples. One pool per treatment, so the tests compare samples within a single pool pair and cannot separate diesel from other pool differences.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  [[`${J.n_pass}`, `MAGs with enough reads for a within-genome profile`],
   [`${sig.length} of ${tested.length}`, 'MAGs whose state differs between the diesel and control pools after dosing in 2025-26 (PERMANOVA, q < 0.05)'],
   [`${d3.format('.0%')(d3.median(tested, m => m.R2_25))}`, 'median share of within-genome variation explained by pool after dosing, 2025-26 (R²)'],
   [`${(J.axis_counts || {})['ice vs water'] || 0} of ${J.n_pass}`, 'MAGs whose main axis of within-genome variation separates ice from water (more than winter or pool)']]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  sec.append(el('div', { class: 'note' }, '<b>How to read this.</b> Each dot is one RNA library, placed by how that genome divided its transcripts among its genes. PC1 is the direction in which the genome’s expression varies most across samples, PC2 the next. The axes have no fixed meaning. The function panel shows what moves along them, for example ribosomal proteins (growth) or stress genes. Depth r is the correlation of PC1 with the genome’s read count: a value near ±1 warns that PC1 may partly reflect sampling depth rather than biology. A shift ratio above 1 means the diesel-pool water moved further from its own pre-dosing state than the control pool moved from its own.'));
  stTable(sec, J);
  stExplorer(sec, J);
  stReading(sec, J);
}

/* ---------- ranked table ---------- */
function stTable(sec, J) {
  sec.append(el('h3', {}, 'Which genomes change state with diesel'));
  const box = fig(sec, '', 'Ranked by R² in 2025-26: the share of within-genome variation (top 5 PCs) separating diesel-pool from control-pool samples after dosing, with permutations within water and ice. q values are Benjamini–Hochberg across MAGs. Shift = how far diesel-pool water moved from its own pre-dosing state, divided by how far control-pool water moved (median distances on the top 5 PCs). 2024-25 had one pre-dosing sample per pool and the pools were filled two weeks apart, so its contrast mixes diesel with pool history. Ice R² = water vs ice separation after dosing (both pools). Main axis = which of winter, pool or ice vs water explains most of the top three PCs. Depth r = correlation of PC1 with the genome’s read count; thinned r = agreement of PC1 after thinning every sample to the same depth (near 1 = robust to depth). Bold = key genome discussed on other pages. Click a row with a ▸ to open it below.');
  const fmtR = (r, q) => r == null ? '–' : `${d3.format('.2f')(r)}${q != null && q < 0.05 ? ' <span class="pill up">' + stQ(q) + '</span>' : ''}`;
  const rows = J.mags.map(m => `<tr data-b="${m.detail ? m.bin : ''}" style="${m.detail ? 'cursor:pointer' : ''}"><td>${m.detail ? '▸ ' : ''}${m.key ? '<b>' : ''}${m.taxon === m.bin ? '' : m.taxon}${m.key ? '</b>' : ''} <span class="mono" style="color:var(--ink-3)">${m.bin}</span></td>
    <td class="num">${m.n}</td><td class="num">${stPct(m.ev[0])}</td><td class="num">${stSgn(m.rd1)}</td><td class="num">${m.rr1 == null ? '–' : d3.format('.2f')(m.rr1)}</td><td class="num">${m.axis ? `${m.axis} <span style="color:var(--ink-3)">${d3.format('.2f')(m.axis_r2)}</span>` : '–'}</td>
    <td class="num">${fmtR(m.R2_25, m.q_25)}</td><td class="num">${m.shift_25 == null ? '–' : d3.format('.1f')(m.shift_25) + '×'}</td><td class="num">${fmtR(m.R2_24, m.q_24)}</td><td class="num">${fmtR(m.R2ice_25 ?? m.R2ice_24, m.R2ice_25 != null ? m.qice_25 : m.qice_24)}</td></tr>`).join('');
  box.innerHTML = `<div class="tablewrap" style="max-height:520px"><table><thead><tr><th>Genome</th><th class="num">Samples</th><th class="num">PC1</th><th class="num">Depth r</th><th class="num">Thinned r</th><th class="num">Main axis aligns with (R²)</th><th class="num">R² 25-26</th><th class="num">Shift 25-26</th><th class="num">R² 24-25</th><th class="num">Ice R²</th></tr></thead><tbody>${rows}</tbody></table></div>`;
  box.querySelectorAll('tr[data-b]').forEach(tr => { if (!tr.dataset.b) return; tr.onclick = () => { const s = document.getElementById('stsel'); if (s) { s.value = tr.dataset.b; s.onchange(); s.scrollIntoView({ behavior: 'smooth', block: 'center' }); } }; });
}

/* ---------- per-MAG explorer ---------- */
function stExplorer(sec, J) {
  sec.append(el('h3', {}, 'One genome at a time'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'stsel', 'aria-label': 'Genome' });
  const det = J.mags.filter(m => m.detail);
  const only = m => { const w = [...new Set(J.detail[m.bin].idx.map(si => J.samples[si].season))]; return w.length === 1 ? `, ${w[0]} only` : ''; };
  const og1 = document.createElement('optgroup'); og1.label = 'Key genomes'; det.filter(m => m.key).forEach(m => og1.append(new Option(`${m.taxon} (${m.bin}${only(m)})`, m.bin)));
  const og2 = document.createElement('optgroup'); og2.label = 'Strongest diesel shifts'; det.filter(m => !m.key).forEach(m => og2.append(new Option(`${m.taxon} (${m.bin}${only(m)})`, m.bin)));
  sel.append(og1, og2);
  const ss = el('select', { id: 'stsea', 'aria-label': 'Winter' }); [['both', 'Both winters'], ['2025-26', '2025-26'], ['2024-25', '2024-25']].forEach(([v, l]) => ss.append(new Option(l, v)));
  const tr = el('input', { type: 'checkbox', id: 'sttraj' }); tr.checked = true;
  c.append(el('label', { for: 'stsel' }, 'Genome'), sel, el('label', { for: 'stsea' }, 'Show'), ss, tr, el('label', { for: 'sttraj' }, 'join days within a pool'));
  const v = el('div'); sec.append(v);
  function draw() {
    v.innerHTML = '';
    const m = J.mags.find(x => x.bin === sel.value), d = J.detail[m.bin];
    v.append(el('p', { class: 'kicker' }, `<b>${m.taxon}</b> <span class="mono">${m.bin}</span>: ${m.n} samples, ${fmt(m.genes)} genes in the profile. PC1 explains ${stPct(d.ev[0])} and PC2 ${stPct(d.ev[1])} of within-genome variation; PC1 vs log read depth r = ${stSgn(m.rd1)}. 2025-26 diesel vs control after dosing: R² ${m.R2_25 == null ? '–' : d3.format('.2f')(m.R2_25)} ${stQ(m.q_25)}, shift ratio ${m.shift_25 == null ? '–' : d3.format('.1f')(m.shift_25) + '×'}. ${m.axis ? `Its PC1–3 space is best explained by ${m.axis} (R² ${d3.format('.2f')(m.axis_r2)}; winter ${m.R2_season ?? '–'}, pool ${m.R2_pool ?? '–'}, ice vs water ${m.R2_medium ?? '–'}).` : ''} ${(d.trend || []).filter(t => t.r1 != null).map(t => `${t.season} ${t.pool === 'EXP' ? 'diesel' : 'control'}: PC1 vs day r ${stSgn(t.r1)}, growth index ${t.g0 ?? '–'}→${t.g1 ?? '–'}‰`).join('; ')}.`));
    stScatter(v, J, m, d, ss.value, tr.checked);
    stTime(v, J, m, d);
    stFunc(v, J, m, d);
  }
  sel.onchange = draw; ss.onchange = draw; tr.onchange = draw;
  sel.value = (det.find(m => m.bin === 'lorbin_184') || det[0]).bin; draw();
}

const ST_LAYER = { TOP: 'top section', BOT: 'bottom section', HIGH: 'upper section (labelled “High”)', LOW: 'lower section (labelled “Low”)', MID: 'middle section', PANCAKE: 'pancake ice' };
const stLayer = l => ST_LAYER[l] || (l ? l.toLowerCase() : 'whole core (no layer recorded)');

function stScatter(v, J, m, d, season, traj) {
  const box = fig(v, `${m.bin}: PCA of the genome’s own transcriptome`, 'Grey = before dosing (either pool); blue = control pool after dosing; orange = diesel pool after dosing. Circles are seawater, triangles sea ice; open symbols are 2024-25, filled 2025-26. Small faded symbols are shallow libraries (under 10,000 reads on this genome): they are projected onto axes fitted on the deeper libraries, so their sampling noise cannot shape the axes. Lines join per-day means within one pool, season and medium, in time order (they start at the open end and finish at the arrow). Hover a point for the library.');
  const lg = el('div', { class: 'legend' }, `<span><i style="background:var(--other)"></i>before dosing</span><span><i style="background:var(--ctl)"></i>control pool</span><span><i style="background:var(--exp)"></i>diesel pool</span><span style="color:var(--ink-2)">● seawater</span><span style="color:var(--ink-2)">▲ sea ice</span><span style="color:var(--ink-2)">○ open = 2024-25</span><span style="color:var(--ink-2)">● filled = 2025-26</span>`);
  box.parentNode.insertBefore(lg, box);
  const pts = d.idx.map((si, k) => ({ s: J.samples[si], x: d.pc[k][0], y: d.pc[k][1], depth: d.depth[k], cl: d.cl[k] })).filter(p => season === 'both' || p.s.season === season);
  // no profile in the chosen winter: draw every sample greyed out, with a note, instead of an empty plot
  const empty = !pts.length;
  const have = [...new Set(d.idx.map(si => J.samples[si].season))].sort().join(' and ');
  if (empty) { pts.push(...d.idx.map((si, k) => ({ s: J.samples[si], x: d.pc[k][0], y: d.pc[k][1], depth: d.depth[k], cl: d.cl[k] }))); traj = false; }
  const W = width(box), H = Math.min(460, Math.max(300, W * 0.6)), mg = { l: 46, r: 14, t: 10, b: 38 };
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H).attr('role', 'img').attr('aria-label', 'PCA scatter');
  const ext = (a) => { const [lo, hi] = d3.extent(a); const p = (hi - lo) * 0.06 || 1; return [lo - p, hi + p]; };
  const x = d3.scaleLinear().domain(ext(pts.map(p => p.x))).range([mg.l, W - mg.r]), y = d3.scaleLinear().domain(ext(pts.map(p => p.y))).range([H - mg.b, mg.t]);
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H - mg.b})`).call(d3.axisBottom(x).ticks(6));
  svg.append('g').attr('class', 'axis').attr('transform', `translate(${mg.l},0)`).call(d3.axisLeft(y).ticks(5));
  svg.append('text').attr('x', W - mg.r).attr('y', H - 4).attr('text-anchor', 'end').style('font-size', '11px').attr('fill', css('--ink-3')).text(`PC1 (${stPct(d.ev[0])})`);
  svg.append('text').attr('transform', `translate(12,${mg.t}) rotate(-90)`).attr('text-anchor', 'end').style('font-size', '11px').attr('fill', css('--ink-3')).text(`PC2 (${stPct(d.ev[1])})`);
  if (traj) {
    const g = d3.group(pts, p => `${p.s.season}|${p.s.pool}|${p.s.medium}`);
    for (const [k, a] of g) {
      const byDay = d3.rollups(a, q => ({ x: d3.mean(q, p => p.x), y: d3.mean(q, p => p.y) }), p => p.s.day).sort((p, q) => p[0] - q[0]).map(z => z[1]);
      if (byDay.length < 2) continue;
      const pool = k.split('|')[1];
      svg.append('path').attr('d', d3.line().x(p => x(p.x)).y(p => y(p.y))(byDay)).attr('fill', 'none').attr('stroke', css(pool === 'EXP' ? '--exp' : '--ctl'))
        .attr('stroke-width', 1.5).attr('opacity', 0.45).attr('stroke-dasharray', k.endsWith('SI') ? '4 3' : null);
      const e = byDay[byDay.length - 1], f = byDay[byDay.length - 2];
      const ang = Math.atan2(y(e.y) - y(f.y), x(e.x) - x(f.x)) * 180 / Math.PI;
      svg.append('path').attr('d', d3.symbol(d3.symbolTriangle, 36)()).attr('transform', `translate(${x(e.x)},${y(e.y)}) rotate(${ang + 90})`).attr('fill', css(pool === 'EXP' ? '--exp' : '--ctl')).attr('opacity', 0.6);
    }
  }
  pts.sort((a, b) => (a.s.phase === 'pre') - (b.s.phase === 'pre')).reverse();
  svg.append('g').selectAll('path').data(pts).join('path')
    .attr('d', p => d3.symbol(p.s.medium === 'SI' ? d3.symbolTriangle : d3.symbolCircle, p.depth < (J.fit_reads || 0) ? 36 : 70)())
    .attr('opacity', p => p.depth < (J.fit_reads || 0) ? 0.45 : 1)
    .attr('transform', p => `translate(${x(p.x)},${y(p.y)})`).attr('fill', p => p.s.season === '2024-25' ? css('--panel') : stPhaseCol(p.s)).attr('stroke', p => p.s.season === '2024-25' ? stPhaseCol(p.s) : css('--panel'))
    .attr('stroke-width', p => p.s.season === '2024-25' ? 1.8 : 1.5)
    .on('mousemove', (ev, p) => showTip(ev, `<b>${p.s.name}</b><br>${p.s.season}, day ${p.s.day}, ${p.s.pool === 'EXP' ? 'diesel' : 'control'} pool, ${p.s.medium === 'SI' ? 'sea ice, ' + stLayer(p.s.layer) : 'seawater'}<br>${fmt(p.depth)} reads on this genome${p.depth < (J.fit_reads || 0) ? ' · shallow: projected onto the axes, not used to fit them' : ''}`))
    .on('mouseleave', hideTip);
  if (empty) {
    svg.selectAll('g:last-of-type path').each(function (p) { const open = p.s.season === '2024-25';
      d3.select(this).attr('fill', open ? css('--panel') : css('--ink-3')).attr('stroke', open ? css('--ink-3') : css('--panel')).attr('opacity', 0.35); });
    svg.append('text').attr('x', (mg.l + W - mg.r) / 2).attr('y', (mg.t + H - mg.b) / 2).attr('text-anchor', 'middle').style('font', '600 14px var(--f-display)').attr('fill', css('--ink-2'))
      .text(`No ${season} profile for ${m.bin} (too few reads). Shown greyed: ${have}.`);
  }
}

function stTime(v, J, m, d) {
  const pc1 = new Array(J.samples.length).fill(null), gr = new Array(J.samples.length).fill(null);
  d.idx.forEach((si, k) => { pc1[si] = d.pc[k][0]; gr[si] = d.growth[k] == null ? null : d.growth[k] * 1000; });
  const b1 = fig(v, `${m.bin}: PC1 over time`, 'PC1 score per library; lines join per-day means. Samples with too few reads on this genome are left out.' + (d.cp.length ? ' Largest single shift along PC1 (changepoint): ' + d.cp.filter(c => Math.abs(c.eff) >= 2).map(c => `${c.season} ${c.pool === 'EXP' ? 'diesel' : 'control'} ${c.medium === 'SI' ? 'ice' : 'water'} after day ${c.day} (${stSgn(c.eff)} SD)`).join('; ') + '.' : ''));
  legendPools(b1.parentNode.insertBefore(el('div'), b1)); timecourse(b1, stSeries(J, pc1), { ylab: 'PC1', unit: '', h: 180, log: false });
  if (d.growth.some(g => g != null)) {
    const b2 = fig(v, `${m.bin}: growth index`, 'Ribosomal-protein transcripts as ‰ of the genome’s own transcripts. Fast-growing bacteria devote more of their transcription to ribosomes, so this is a rough per-genome growth-rate proxy.');
    legendPools(b2.parentNode.insertBefore(el('div'), b2)); timecourse(b2, stSeries(J, gr), { ylab: '‰ ribosomal', unit: '‰', h: 180, log: false });
  }
}

function stFunc(v, J, m, d) {
  const F = d.func.filter(f => f.r1 != null).sort((a, b) => Math.abs(b.r1) - Math.abs(a.r1));
  const box = fig(v, `${m.bin}: what moves along PC1 and PC2`, 'Spearman correlation between each function’s share of the genome’s transcripts and the PC score, across samples. Purple bars (right) = the function gains share towards positive PC values; green bars (left) = towards negative values. Below, the genes with the largest PC1 loadings at each end.');
  const W = width(box), lw = Math.min(240, W * 0.4), rh = 22, H = F.length * rh + 34;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const half = (W - lw - 50) / 2;
  [[0, 'r1', 'PC1'], [1, 'r2', 'PC2']].forEach(([j, k, lab]) => {
    const x0 = lw + 10 + j * (half + 30), x = d3.scaleLinear().domain([-1, 1]).range([x0, x0 + half]);
    svg.append('text').attr('x', x(0)).attr('y', 11).attr('text-anchor', 'middle').style('font-size', '12px').style('font-weight', 600).attr('fill', css('--ink-2')).text(lab);
    svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 16).attr('y2', H - 18).attr('stroke', css('--rule'));
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H - 18})`).call(d3.axisBottom(x).tickValues([-1, -0.5, 0, 0.5, 1]).tickFormat(d3.format('~g')).tickSize(3));
    F.forEach((f, i) => { const yv = 20 + i * rh, val = f[k] ?? 0;
      svg.append('rect').attr('x', Math.min(x(0), x(val))).attr('y', yv).attr('width', Math.abs(x(val) - x(0))).attr('height', 14).attr('rx', 2)
        .attr('fill', css(val >= 0 ? '--s7' : '--s3')).on('mousemove', ev => showTip(ev, `<b>${f.cat}</b> (${f.n} genes, ${d3.format('.2%')(f.share)} of transcripts)<br>r with ${lab}: ${stSgn(val)}`)).on('mouseleave', hideTip); });
  });
  F.forEach((f, i) => svg.append('text').attr('x', lw).attr('y', 20 + i * rh + 11).attr('text-anchor', 'end').style('font-size', '12px').attr('fill', css('--ink')).text(f.cat));
  const L = d.load['0'], li = a => a.map(g => `<li>${g.p || 'hypothetical'} <span class="mono" style="color:var(--ink-3)">${g.g}</span></li>`).join('');
  const tl = el('div', { class: 'cols2', style: 'display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:16px;font-size:13px;margin-top:8px' },
    `<div><b>Towards negative PC1</b><ul style="margin:4px 0 0 18px;padding:0">${li(L.neg)}</ul></div><div><b>Towards positive PC1</b><ul style="margin:4px 0 0 18px;padding:0">${li(L.pos)}</ul></div>`);
  box.append(tl);
}

/* ---------- reading ---------- */
function stReading(sec, J) {
  const ul = el('ul', { class: 'kicker' });
  ul.innerHTML = (J.reading || []).map(t => `<li>${t}</li>`).join('');
  if (J.reading && J.reading.length) { sec.append(el('h3', {}, 'What the states show')); sec.append(ul); }
  sec.append(el('p', { class: 'note' }, `Robustness: for the key genomes, each sample’s reads were thinned down to the smallest depth used, and the PCA was repeated. ${(J.rarefy || []).map(r => `${r.bin} PC1 |r| = ${d3.format('.2f')(r.r1)}`).join(', ')} between full and thinned data. Values near 1 mean the axes are not driven by depth.`));
}

RENDER.states = statesSec;
