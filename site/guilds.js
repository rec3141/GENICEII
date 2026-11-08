/* Response guilds: all profiled genomes in one space, by what drives their within-genome state (fingerprints) and by which
   genomes change state together across samples (shared factors). Data: guilds.json (guilds/scripts/guilds.py, guilds_site.py).
   Uses helpers from app.js (el, fig, width, css, showTip, hideTip, fmt, fmt1, timecourse, legendPools, D, RENDER). */
const guP = d3.format('.0%'), guR = v => v == null ? '–' : d3.format('.2f')(v), guSgn = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.2f')(v);
const GU_DR = [['winter', 'winter'], ['ice', 'ice vs water'], ['diesel', 'diesel pool'], ['time', 'time'], ['light', 'under-ice light']];
const GU_CAT = ['--s1', '--s2', '--s3', '--s4', '--s5', '--s6', '--s7', '--s8'];
let guSel = null;
const guNm = (t, b) => t === b ? `<span class="mono">${b}</span>` : `${t} <span class="mono">${b}</span>`;

function guColour(J, mode) {
  const list = mode === 'fc' ? J.clusters : mode === 'guild' ? J.guilds : mode === 'dom' ? J.dominants : J.phyla;
  const key = g => mode === 'fc' ? g.fc : mode === 'guild' ? g.guild : mode === 'dom' ? g.dom : g.phylum_c;
  const col = v => { const i = list.indexOf(v); return v === 'other' || i < 0 || i >= GU_CAT.length ? css('--other') : css(GU_CAT[i]); };
  return { list, key, col };
}

function guilds(sec) {
  const J = D.guilds;
  if (!J) { sec.innerHTML = '<h2>Response guilds</h2><p class="kicker">guilds.json is missing.</p>'; return; }
  const M = J.model, g = J.genomes;
  sec.innerHTML = `<h2>Response guilds</h2>
  <p class="kicker">The Physiological states page looks at one genome at a time. This page puts all ${M.n_genomes} profiled genomes in one space, two ways. The <b>fingerprint map</b> places genomes by <i>what</i> drives their state: how much of each genome’s within-genome variation is explained by winter, ice versus water, the diesel pool, time within a winter and modelled under-ice light, and how much is shared among the seasonal drivers. Season is not treated as noise: microbes feel it through light, temperature, ice growth and nutrients, which change together, so the shared seasonal part is reported as a driver in its own right. The <b>guild map</b> places genomes by <i>when</i> they change: a shared-factor model finds the few sample-level patterns that many genomes follow at once, and genomes that load on the same factor change state together. Genomes with the same fingerprint but different guilds respond to the same driver on different schedules; genomes in the same guild but from different phyla are candidates for interaction, such as a degrader and its cross-feeders.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  [[`${M.n_genomes}`, `genomes with a within-genome profile, over ${M.n_samples} libraries`],
   [`${J.sig.diesel[0]} of ${J.sig.diesel[1]}`, 'genomes whose state has a unique diesel-pool component (q < 0.05, permutation)'],
   [`${guP(J.med.r2)}`, 'median share of a genome’s state variation (PC1–3) explained by the five drivers together'],
   [`${M.K}`, `shared factors that many genomes follow (split-half stability ≥ 0.7)`]]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  sec.append(el('div', { class: 'note' }, `<b>How to read this.</b> Each dot is a genome; size = mean transcript share (TPM). On the fingerprint map, nearby genomes have similar driver profiles; on the guild map, nearby genomes load on the same shared factors. Colour by cluster, guild, dominant driver or phylum. Hover for the genome’s driver bars; click to pin it below, with a link to its own view on the Physiological states page. Unique R² = variance only that driver explains; “seasonal” = variance winter, time and light explain jointly but not separately. Genomes profiled in one winter only have no winter component (shown as –) and are drawn with a dashed outline.`));
  guMaps(sec, J);
  const det = el('div', { id: 'gu-detail' }); sec.append(det);
  guFactors(sec, J);
  guPairs(sec, J);
  guTables(sec, J);
  sec.append(el('h3', {}, 'Caveats'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li>One pool per treatment: “diesel pool” is a pool effect (history plus dosing), not a replicated treatment. In 2024-25 the pools were filled about two weeks apart after the intake lines froze, so part of it is pool history.</li>
    <li>Winter, time and light change together within a winter; the commonality split shows how much each explains uniquely and how much only jointly. A small unique light fraction does not mean light is unimportant, only that it cannot be separated from the rest of the season.</li>
    <li>${g.filter(x => x.single).length} genomes were profiled in one winter only; their winter component is missing and was set to the median for the map.</li>
    <li>Within-genome states come from each MAG’s PC1–3 (Physiological states page). A strain switch inside a MAG looks like a change of state.</li>
    <li>Clusters are soft (fingerprint: ${M.fmethod}; guilds: ${M.gmethod}); the maps are better read as gradients than as discrete groups.</li>`));
}

/* ---------- the two maps ---------- */
function guMaps(sec, J) {
  sec.append(el('h3', {}, 'Two maps of the same genomes'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'gucol', 'aria-label': 'Colour by' });
  [['fc', 'fingerprint cluster'], ['guild', 'guild'], ['dom', 'dominant driver'], ['phy', 'phylum']].forEach(([v, l]) => sel.append(new Option(l, v)));
  c.append(el('label', { for: 'gucol' }, 'Colour by'), sel);
  const wrap = el('div'); sec.append(wrap);
  function draw() {
    wrap.innerHTML = '';
    const C = guColour(J, sel.value);
    const lg = el('div', { class: 'legend' }, C.list.slice(0, 8).map((v, i) => `<span><i style="background:var(${GU_CAT[i]})"></i>${v}</span>`).join('') + (C.list.length > 8 ? '<span><i style="background:var(--other)"></i>other</span>' : ''));
    wrap.append(lg);
    const row = el('div', { class: 'gu-row' }); wrap.append(row);
    const h1 = el('div', { class: 'gu-cell' }), h2 = el('div', { class: 'gu-cell' }); row.append(h1, h2);   // both cells exist before sizing
    guScatter(h1, J, 'fx', 'fy', 'Fingerprint map (what drives each genome)', 'UMAP of each genome’s driver fractions: unique R² for winter, ice vs water, diesel pool, time and light, plus the shared seasonal part.', C);
    guScatter(h2, J, 'gx', 'gy', 'Guild map (which genomes change together)', `UMAP of each genome’s loadings on the ${J.model.K} shared factors (R² of its PC1–3 on each factor).`, C);
  }
  sel.onchange = draw; draw();
}

function guScatter(holder, J, kx, ky, title, cap, C) {
  const box = fig(holder, title, cap);
  const W = width(box, 540), H = Math.min(420, Math.max(300, W * 0.8)), mg = 14;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H).attr('role', 'img').attr('aria-label', title);
  const G = J.genomes;
  const x = d3.scaleLinear().domain(d3.extent(G, d => d[kx])).nice().range([mg, W - mg]), y = d3.scaleLinear().domain(d3.extent(G, d => d[ky])).nice().range([H - mg, mg]);
  const rs = d3.scaleSqrt().domain([0, d3.max(G, d => d.tpm || 0)]).range([3, 16]);
  const pts = [...G].sort((a, b) => (b.tpm || 0) - (a.tpm || 0));
  svg.append('g').selectAll('circle').data(pts).join('circle')
    .attr('cx', d => x(d[kx])).attr('cy', d => y(d[ky])).attr('r', d => rs(d.tpm || 0))
    .attr('fill', d => C.col(C.key(d))).attr('fill-opacity', 0.85).attr('stroke', d => d.bin === guSel ? css('--ink') : css('--panel'))
    .attr('stroke-width', d => d.bin === guSel ? 2.5 : 1.5).attr('stroke-dasharray', d => d.single ? '2 2' : null).style('cursor', 'pointer')
    .on('mousemove', (ev, d) => showTip(ev, guTip(d)))
    .on('mouseleave', hideTip)
    .on('click', (ev, d) => { guSel = d.bin; document.querySelectorAll('#v-guilds svg circle').forEach(c => { const dd = d3.select(c).datum(); if (dd && dd.bin) d3.select(c).attr('stroke', dd.bin === guSel ? css('--ink') : css('--panel')).attr('stroke-width', dd.bin === guSel ? 2.5 : 1.5); }); guDetail(J, d); });
}

function guBars(d) {
  const rows = GU_DR.map(([k, l]) => [l, d.u[k]]).concat([['seasonal (shared)', d.seasonal], ['other shared', d.other]]);
  return rows.map(([l, v]) => `<div style="display:flex;align-items:center;gap:6px;font-size:12px"><span style="width:110px;color:var(--ink-2)">${l}</span><span style="display:inline-block;height:8px;width:${v == null ? 0 : Math.max(0, v) * 160}px;background:var(--accent);border-radius:2px"></span><span>${v == null ? '–' : guP(v)}</span></div>`).join('');
}
function guTip(d) {
  return `<b>${guNm(d.taxon, d.bin)}</b><br>${d.phylum} · TPM ${fmt1(d.tpm || 0)}<br>cluster: ${d.fc}<br>guild: ${d.guild}<br>drivers explain ${guP(d.r2)} of its state variation${d.single ? ' (one winter only)' : ''}<div style="margin-top:4px">${guBars(d)}</div>`;
}

function guDetail(J, d) {
  const det = document.getElementById('gu-detail'); det.innerHTML = '';
  const box = fig(det, `${d.taxon === d.bin ? d.bin : d.taxon + ' ' + d.bin}: what drives its state`, 'Unique R²: variance of the genome’s PC1–3 explained by that driver alone (q from within-strata permutations, BH across genomes). Seasonal: explained jointly by winter, time and light but by none alone. Function directions: correlation of each function’s share of the genome’s transcripts with its PC1.');
  const fl = (d.fn || []).map(([k, v]) => `<div style="font-size:12.5px">${k}: <b>${guSgn(v)}</b></div>`).join('');
  const ld = J.factors.map(f => `<div style="display:flex;align-items:center;gap:6px;font-size:12px"><span style="width:28px;color:var(--ink-2)" title="${f.top}">${f.id}</span><span style="display:inline-block;height:8px;width:${Math.max(0, d.load[f.id] || 0) * 160}px;background:var(--ink-3);border-radius:2px"></span><span>${guP(d.load[f.id] || 0)}</span></div>`).join('');
  const qs = GU_DR.map(([k, l]) => `${l} ${d.q[k] == null ? '–' : d.q[k] < 0.001 ? '&lt;0.001' : d3.format('.2g')(d.q[k])}`).join(' · ');
  box.innerHTML = `<div class="gu-row"><div class="gu-cell"><b style="font:600 13px var(--f-display)">Driver fractions</b>${guBars(d)}<div style="font-size:11.5px;color:var(--ink-3);margin-top:4px">q: ${qs}</div></div>
    <div class="gu-cell"><b style="font:600 13px var(--f-display)">Loadings on shared factors</b>${ld}<div style="font-size:11.5px;color:var(--ink-3)">${J.factors.map(f => `${f.id}: ${f.top}`).join('<br>')}</div><b style="font:600 13px var(--f-display);display:block;margin-top:8px">Function directions along PC1</b>${fl || '–'}</div></div>
    <p style="margin:8px 0 0"><a href="#states" data-bin="${d.bin}" id="gu-go">Open ${d.bin} on the Physiological states page →</a></p>`;
  box.querySelector('#gu-go').onclick = () => { try { localStorage.setItem('cmo-states-mag', d.bin); } catch (e) { /* storage blocked */ } };
  det.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}

/* ---------- shared factors ---------- */
function guSeries(J, values) {
  const out = {};
  J.samples.forEach((s, i) => { if (s.day == null || values[i] == null) return; const k = `${s.season}|${s.pool}|${s.medium}`; ((out[k] ||= {})[s.day] ||= []).push([values[i], s.name]); });
  return Object.entries(out).map(([k, byDay]) => { const [season, pool, medium] = k.split('|');
    return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a, x => x[0]), name: a.map(x => x[1]).join(', ') })) }; });
}
function guFactors(sec, J) {
  sec.append(el('h3', {}, 'The shared factors'));
  const M = J.model;
  sec.append(el('p', { class: 'kicker' }, `A factor model over all genomes’ state scores (missing where a genome had too few reads) finds ${M.K} sample-level patterns that recur across genomes; together they capture ${guP(M.ve[M.K])} of the genomes’ state variation, and each is reproduced by two random halves of the genomes (split-half match ${d3.format('.2f')(M.stab[M.K])}). Each factor is a set of samples in which many genomes shift the same way; the table shows what each tracks.`));
  const box = fig(sec, 'What each factor tracks', 'Spearman correlation between the factor’s sample scores and each driver across libraries. Orange = positive, blue = negative.');
  const drv = Object.keys(J.factors[0].r);
  const col = v => { if (v == null) return 'transparent'; const t = Math.max(-1, Math.min(1, v)); return t >= 0 ? d3.interpolateRgb(css('--panel'), css('--exp'))(t * 0.85) : d3.interpolateRgb(css('--panel'), css('--ctl'))(-t * 0.85); };
  box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Factor</th><th>Tracks</th>${drv.map(d => `<th class="num">${d}</th>`).join('')}</tr></thead><tbody>${
    J.factors.map(f => `<tr><td><b>${f.id}</b></td><td>${f.top}</td>${drv.map(d => `<td class="num" style="background:${col(f.r[d])}">${guSgn(f.r[d])}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`;
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'gufac', 'aria-label': 'Factor' }); J.factors.forEach(f => sel.append(new Option(`${f.id}: ${f.top}`, f.id))); c.append(el('label', { for: 'gufac' }, 'Factor'), sel);
  const v = el('div'); sec.append(v);
  function draw() { v.innerHTML = ''; const b = fig(v, `${sel.value} sample scores over time`, 'Factor score per library (standardised); lines join per-day means.');
    legendPools(b.parentNode.insertBefore(el('div'), b)); timecourse(b, guSeries(J, J.scores[sel.value]), { ylab: 'factor score', unit: '', h: 180, log: false }); }
  sel.onchange = draw; draw();
}

/* ---------- comparing the layers ---------- */
function guPairs(sec, J) {
  sec.append(el('h3', {}, 'Where the two maps disagree'));
  sec.append(el('p', { class: 'kicker' }, 'Same fingerprint, different guild: the same driver mix, but the genomes change state in different samples, for example at different times after dosing or in different layers. Same guild, different phylum: unrelated genomes whose states move together across samples, candidates for ecological coupling (shared substrate, cross-feeding, or a common stressor).'));
  const t = (rows, kind) => `<div class="tablewrap"><table><thead><tr><th>Genome A</th><th>Genome B</th><th>${kind === 'fp' ? 'Fingerprint cluster' : 'Guild'}</th><th>${kind === 'fp' ? 'Guilds' : 'Phyla'}</th><th class="num">${kind === 'fp' ? 'driver difference' : 'loading similarity'}</th></tr></thead><tbody>${
    rows.map(p => `<tr><td>${guNm(p.ta, p.a)}</td><td>${guNm(p.tb, p.b)}</td><td>${kind === 'fp' ? p.ca : p.ga}</td><td>${kind === 'fp' ? `${p.ga} / ${p.gb}` : `${p.pa} / ${p.pb}`}</td><td class="num">${guR(p.v)}</td></tr>`).join('') || '<tr><td colspan="5">none</td></tr>'}</tbody></table></div>`;
  const b1 = fig(sec, 'Same fingerprint, different guild', 'Abundant pairs in one fingerprint cluster whose driver fractions differ by less than 0.08 in total, but which sit in different guilds.'); b1.innerHTML = t(J.pairs_same_fp, 'fp');
  const b2 = fig(sec, 'Same guild, different phylum', 'Abundant pairs in one guild, from different phyla, with near-identical factor loadings (cosine ≥ 0.97).'); b2.innerHTML = t(J.pairs_same_guild, 'guild');
}

function guTables(sec, J) {
  sec.append(el('h3', {}, 'All genomes'));
  const xt = J.xtab, gs = Object.keys(xt), cs = [...new Set(gs.flatMap(g => Object.keys(xt[g])))];
  const b0 = fig(sec, 'Fingerprint clusters × guilds', 'Number of genomes in each combination.');
  b0.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Fingerprint cluster</th>${gs.map(g => `<th class="num">${g}</th>`).join('')}</tr></thead><tbody>${
    cs.map(c => `<tr><td>${c}</td>${gs.map(g => `<td class="num">${xt[g][c] || ''}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`;
  const b = fig(sec, '', 'Sorted by transcript share. Unique R² per driver; bold where q < 0.05. – = not estimable (single winter, or the driver does not vary in that genome’s samples).');
  const G = [...J.genomes].sort((a, c) => (c.tpm || 0) - (a.tpm || 0));
  const cell = (d, k) => { const v = d.u[k], q = d.q[k]; return `<td class="num">${v == null ? '–' : (q != null && q < 0.05 ? `<b>${guP(v)}</b>` : guP(v))}</td>`; };
  b.innerHTML = `<div class="tablewrap" style="max-height:520px"><table><thead><tr><th>Genome</th><th>Phylum</th><th>Cluster</th><th>Guild</th><th>Dominant</th>${GU_DR.map(([, l]) => `<th class="num">${l}</th>`).join('')}<th class="num">seasonal</th><th class="num">all drivers</th><th class="num">TPM</th></tr></thead><tbody>${
    G.map(d => `<tr data-b="${d.bin}" style="cursor:pointer"><td>${guNm(d.taxon, d.bin)}</td><td>${d.phylum}</td><td>${d.fc}</td><td>${d.guild}</td><td>${d.dom}</td>${GU_DR.map(([k]) => cell(d, k)).join('')}<td class="num">${guP(d.seasonal || 0)}</td><td class="num">${guP(d.r2)}</td><td class="num">${fmt1(d.tpm || 0)}</td></tr>`).join('')}</tbody></table></div>`;
  b.querySelectorAll('tr[data-b]').forEach(tr => { tr.onclick = () => { const d = J.genomes.find(x => x.bin === tr.dataset.b); guSel = d.bin; guDetail(J, d); }; });
}

RENDER.guilds = guilds;
