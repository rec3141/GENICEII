/* Network section (Food web): genome-genome association networks per layer (DNA, RNA), consensus of FlashWeave, SPRING,
   SPIEC-EASI, a design-residual Spearman baseline and a presence-only Pearson. Data: network.json (network/scripts/network_site.py). Helpers from app.js. */
const NW_DOM = { prokaryote: ['--s1', 'Bacteria and archaea'], eukaryote: ['--s3', 'Eukaryotes'], virus: ['--s8', 'Viruses'], organelle: ['--s4', 'Plastids'] };
const NW_VIEWS = [['all', 'All links'], ['grazer', 'Grazers with bacteria and algae'], ['virus', 'Viruses with their partners'], ['plastid', 'Plastids with eukaryotes']];
const NW_CONTAM = 'contamination 10–20%', NW_FRAG = 'fragment';
const nwPop = m => { let c = 0; while (m) { c += m & 1; m >>= 1; } return c; };

const NW_ROWS = {};      // per-layer data (network_{layer}.json), loaded on first use
async function nwLayer(l) {
  if (NW_ROWS[l]) return NW_ROWS[l];
  const L = await (await fetch(`network_${l}.json`, { cache: 'no-cache' })).json(), E = L.edges, n = E.a.length, M = D.network.methods, nb = M.length;
  // unpack the flags (bit 0 positive sign, then the masks of agreeing, opposing and strong-scoring methods) and expand the sparse
  // score columns (one value per edge the method reports, in edge order)
  const full = (1 << nb) - 1;
  const sameRows = E.same || []; E.same = new Uint8Array(n); sameRows.forEach(k => { E.same[k] = 1; });   // 'possible same genome'
  E.s = E.f.map(f => f & 1 ? 1 : -1); E.m = E.f.map(f => (f >> 1) & full); E.x = E.f.map(f => (f >> (1 + nb)) & full); E.st = E.f.map(f => (f >> (1 + 2 * nb)) & full);
  const dense = (bit, col) => { const v = new Array(n).fill(0); let k = 0; for (let i = 0; i < n; i++) if ((E.m[i] | E.x[i]) & bit) v[i] = col[k++]; return v; };
  M.forEach(m => { E[m.col] = dense(m.bit, E[m.col]); });
  NW_ROWS[l] = L; return L;
}

async function networkSec(sec) {
  const J = D.network;
  if (!J) { sec.innerHTML = '<h2>Network</h2><p class="kicker">network.json is missing.</p>'; return; }
  await Promise.all(J.layer_list.map(l => nwLayer(l.key)));
  const cons = l => J.layers[l].summary.find(s => s.method === 'consensus');
  sec.innerHTML = `<h2>Network</h2>
  <p class="kicker">Which genomes rise and fall together? Each point is a genome: a bacterial or archaeal MAG, a eukaryotic MAG, a virus or a plastid. There are three networks. The two DNA networks ask who is there together: one uses nanopore metagenomes, the other the ${J.layers.mgdna ? J.layers.mgdna.n_libs : ''} Illumina metagenomes of the pools. The RNA network uses the ${J.layers.rna.n_libs} Illumina RNA libraries (metatranscriptomes) and asks who is active together. A genome’s RNA level is the typical level of its contigs, so a few very busy contigs (such as rRNA genes) cannot set it. Five methods look for links, and the slider sets how much agreement a link needs: from four or five methods with strong scores to any single method. The default is the consensus: at least two methods agree on the sign, and one of them has already accounted for winter, pool, sea ice or seawater, and time. A link means two genomes co-vary more than these shared drivers explain. It is a lead for grazing, infection or a shared niche, not proof of any of them.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  const LL = J.layer_list.slice().reverse(), names = LL.map(l => l.label).join(' / ');
  [[LL.map(l => J.layers[l.key].n_nodes).join(' / '), `genomes in the ${names} network`],
   [LL.map(l => fmt(cons(l.key).edges)).join(' / '), `consensus links, ${names}`]]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  nwGraph(sec, J);
  nwPartners(sec, J);
  nwMethods(sec, J);
  nwReading(sec, J);
}

// one permissivity step (J.steps, network_site.py): the union of its clauses; a clause passes when >= k of the switched-on methods
// support the edge (strong clauses: only methods with a strong score count), no switched-on method has the opposite sign, and
// design clauses need FlashWeave or the residual Spearman among the supporters. The sign is the one the switched-on methods give
// (a method can be outvoted overall and still stand alone when it is the only one on). Returns the supporting-method mask
// (negative when the switched-on methods give the opposite of the stored sign) or 0.
function nwPass(E, i, step, bits, designBits) {
  const ag = E.m[i] & bits, op = E.x[i] & bits;
  if (ag && op) return 0;
  const on = ag || op, flip = !ag;
  if (!on) return 0;
  for (const c of step.clauses) { let sup = on; if (c.lvl === 'strong') sup &= E.st[i];
    if (nwPop(sup) < Math.min(c.k, nwPop(bits)) || (c.design && !(sup & designBits))) continue;
    return flip ? -on : on; }
  return 0;
}

function nwGraph(sec, J) {
  sec.append(el('h3', {}, 'The networks'));
  const M = J.methods, allBits = M.reduce((a, m) => a | m.bit, 0), designBits = M.filter(m => m.design_aware).reduce((a, m) => a | m.bit, 0);
  const S = { layer: J.layer_list[0].key, view: 'all', bits: allBits, step: J.default_step, sel: null, layout: 'tsne', sign: 'both', find: null, sample: null, hideSame: false, cut: {} };   // cut: score cutoff per method (single-method mode)
  const c1 = el('div', { class: 'controls' }); sec.append(c1);
  const chip = (parent, opts, key) => opts.forEach(([k, l]) => { const b = el('button', { class: 'chip', 'aria-pressed': S[key] === k }, l);
    b.onclick = () => { S[key] = k; S.sel = null; if (key === 'layer') { S.find = null; fi.value = ''; fmsg.textContent = ''; S.sample = null; stopPlay(); fillSamples(); } parent.querySelectorAll(`button[data-g="${key}"]`).forEach(x => x.setAttribute('aria-pressed', x === b)); draw(); };
    b.dataset.g = key; parent.append(b); });
  chip(c1, J.layer_list.map(l => [l.key, l.long]), 'layer');
  const vs = el('select', { id: 'nwview', 'aria-label': 'View' }); NW_VIEWS.forEach(([k, l]) => vs.append(new Option(l, k)));
  vs.onchange = () => { S.view = vs.value; S.sel = null; draw(); }; c1.append(el('label', { for: 'nwview' }, 'Show'), vs);
  const ls = el('select', { id: 'nwlayout', 'aria-label': 'Layout' }); J.layouts.forEach(o => ls.append(new Option(o.key === 'tsne' ? 't-SNE (profile shape)' : `${o.label} links`, o.key)));
  ls.onchange = () => moveTo(ls.value); c1.append(el('label', { for: 'nwlayout' }, 'Layout'), ls);
  const c0 = el('div', { class: 'controls' }); sec.append(c0);
  c0.append(el('span', { style: 'color:var(--ink-2)' }, 'Links'));
  chip(c0, [['both', 'Positive and negative'], ['pos', 'Positive only'], ['neg', 'Negative only']], 'sign');
  const hs = el('input', { type: 'checkbox', id: 'nwhidesame' }); hs.onchange = () => { S.hideSame = hs.checked; draw(); };
  c0.append(hs, el('label', { for: 'nwhidesame' }, 'Hide possible same-genome links'));
  // find a genome: id or any part of its name; selects it (its links only), rings it and scrolls the map into view
  const fi = el('input', { type: 'search', id: 'nwfind', placeholder: 'e.g. Cycloclasticus or lorbin_184', list: 'nwfind-list', style: 'min-width:240px', 'aria-label': 'Find genome' });
  const fl = el('datalist', { id: 'nwfind-list' }), fmsg = el('span', { style: 'color:var(--ink-2);font-size:13px' });
  c0.append(el('label', { for: 'nwfind', style: 'margin-left:14px' }, 'Find genome'), fi, fl, fmsg);
  const fillFind = L => { fl.innerHTML = ''; L.nodes.forEach(n => fl.append(new Option(n.label, n.id))); };
  fi.onchange = async () => {
    const q = fi.value.trim().toLowerCase(), L = await nwLayer(S.layer); if (!q) { S.find = null; fmsg.textContent = ''; draw(); return; }
    const hits = L.nodes.map((n, i) => [n, i]).filter(([n]) => n.id.toLowerCase() === q || n.label.toLowerCase().includes(q) || n.id.toLowerCase().includes(q));
    if (!hits.length) { fmsg.textContent = `Not in the ${layerLabel(S.layer)} network.`; S.find = null; draw(); return; }
    hits.sort((a, b) => (b[0].id.toLowerCase() === q) - (a[0].id.toLowerCase() === q) || (b[0].ab || 0) - (a[0].ab || 0));
    S.find = new Set(hits.map(h => h[1])); S.sel = hits.length === 1 ? hits[0][1] : null;
    fmsg.textContent = hits.length === 1 ? '' : `${hits.length} genomes match: ${hits.slice(0, 6).map(h => h[0].id).join(', ')}${hits.length > 6 ? ', …' : ''}`;
    await draw(); box.scrollIntoView({ behavior: 'smooth', block: 'center' }); };
  // sample picker: node size = that library's value for the genome (network_{layer}_samples.json, loaded on first use)
  // sticky under the section nav so the arrows stay reachable while looking at the map
  const cs = el('div', { class: 'controls', style: 'position:sticky;top:var(--nwstick,44px);z-index:4;background:var(--bg);padding:4px 0;border-bottom:1px solid var(--rule)' }); sec.append(cs);
  const ss = el('select', { id: 'nwsample', 'aria-label': 'Sample', style: 'max-width:420px' });
  const bPrev = el('button', { class: 'chip', title: 'Previous sample', 'aria-label': 'Previous sample' }, '◀'), bNext = el('button', { class: 'chip', title: 'Next sample', 'aria-label': 'Next sample' }, '▶');
  const bPlay = el('button', { class: 'chip', title: 'Play through the samples, one per second', 'aria-pressed': 'false' }, '▶▶');
  const sinfo = el('span', { style: 'color:var(--ink-2);font-size:12.5px' });
  cs.append(el('label', { for: 'nwsample' }, 'Sample'), ss, bPrev, bNext, bPlay, sinfo);
  const SMP = {};
  const smpData = async l => SMP[l] || (SMP[l] = await (await fetch(`network_${l}_samples.json`, { cache: 'no-cache' })).json());
  async function fillSamples() {
    const Sm = await smpData(S.layer); ss.innerHTML = ''; ss.append(new Option('None (size = median abundance)', ''));
    Sm.libs.forEach((b, k) => ss.append(new Option(`${b.label}  —  ${b.id}`, k)));
    Sm.max = d3.quantile(Sm.v.flat().filter(v => v > 0).sort(d3.ascending), 0.995) || 1;      // one scale for every sample of the layer
    ss.value = S.sample == null ? '' : S.sample;
  }
  const setSample = k => { S.sample = k; ss.value = k == null ? '' : k; draw(); };
  const stepSample = d => { const n = SMP[S.layer] ? SMP[S.layer].libs.length : 0; if (!n) return;
    setSample(S.sample == null ? (d > 0 ? 0 : n - 1) : Math.max(0, Math.min(n - 1, S.sample + d))); };
  let play = null;
  const stopPlay = () => { if (play) { clearInterval(play); play = null; } bPlay.setAttribute('aria-pressed', 'false'); };
  ss.onchange = () => { stopPlay(); setSample(ss.value === '' ? null : +ss.value); };
  bPrev.onclick = () => { stopPlay(); stepSample(-1); }; bNext.onclick = () => { stopPlay(); stepSample(1); };
  bPlay.onclick = () => { if (play) { stopPlay(); return; } bPlay.setAttribute('aria-pressed', 'true');
    play = setInterval(() => { const n = SMP[S.layer].libs.length; if (S.sample != null && S.sample >= n - 1) { stopPlay(); return; } stepSample(1); }, 1000); };
  document.addEventListener('keydown', ev => {
    if (location.hash !== '#network' || ev.altKey || ev.ctrlKey || ev.metaKey || !['ArrowLeft', 'ArrowRight'].includes(ev.key)) return;
    const t = ev.target; if (t && (t.tagName === 'INPUT' || t.tagName === 'SELECT' || t.tagName === 'TEXTAREA' || t.isContentEditable)) return;
    ev.preventDefault(); stopPlay(); stepSample(ev.key === 'ArrowRight' ? 1 : -1); });
  fillSamples();
  const c2 = el('div', { class: 'controls' }); sec.append(c2);
  c2.append(el('span', { style: 'color:var(--ink-2)' }, 'Methods'));
  M.forEach(m => { const id = 'nwm-' + m.key, cb = el('input', { type: 'checkbox', id }); cb.checked = true;
    cb.onchange = () => { S.bits = cb.checked ? S.bits | m.bit : S.bits & ~m.bit; draw(); }; c2.append(cb, el('label', { for: id, style: 'margin-right:6px' }, m.label)); });
  // permissivity slider
  const P = J.params, steps = J.steps;
  const pw = el('div', { class: 'controls', style: 'align-items:flex-start;gap:18px;flex-wrap:wrap' }); sec.append(pw);
  const left = el('div', { style: 'flex:1 1 340px;min-width:300px' }); pw.append(left);
  const sl = el('input', { type: 'range', id: 'nwperm', min: 1, max: steps.length, step: 1, value: S.step + 1, style: 'width:100%', 'aria-label': 'Permissivity', list: 'nwperm-ticks' });
  const dl = el('datalist', { id: 'nwperm-ticks' }); steps.forEach((_, k) => dl.append(new Option('', k + 1)));
  const lab = el('div', { style: 'font-weight:600;margin-top:2px' }), cnt = el('div', { style: 'color:var(--ink-2);font-size:13px' });
  const permBox = el('div'), cutBox = el('div', { style: 'display:none' });
  permBox.append(el('label', { for: 'nwperm', style: 'display:flex;justify-content:space-between;font-size:12.5px;color:var(--ink-2)' }, '<span>strict</span><span>Permissivity</span><span>permissive</span>'), sl, dl);
  // single-method mode: a score cutoff for that method instead of the permissivity steps (about 20 quantile steps of its exported
  // scores in the layer, from the lowest to the 99th percentile; the lowest = every link the method reports)
  const ts = el('input', { type: 'range', id: 'nwcut', min: 0, max: 1, step: 1, value: 0, style: 'width:100%', 'aria-label': 'Score cutoff' });
  cutBox.append(el('label', { for: 'nwcut', style: 'display:flex;justify-content:space-between;font-size:12.5px;color:var(--ink-2)' }, '<span>all links</span><span>Score cutoff</span><span>strongest</span>'), ts);
  left.append(permBox, cutBox, lab, cnt);
  sl.oninput = () => { S.step = +sl.value - 1; S.sel = null; draw(); };
  const THR = {};
  const thresholds = (L, m) => { const k = S.layer + '|' + m.key; if (THR[k]) return THR[k];
    const E = L.edges, v = []; for (let i = 0; i < E.a.length; i++) if ((E.m[i] | E.x[i]) & m.bit) v.push(E[m.col][i]);
    v.sort((a, b) => a - b); const t = [];
    if (v.length) for (let q = 0; q < 20; q++) { const x = v[Math.min(v.length - 1, Math.floor(q / 19 * 0.99 * (v.length - 1)))]; if (!t.length || x > t[t.length - 1]) t.push(x); }
    return (THR[k] = t.length ? t : [0]); };          // in exported units (score x 100)
  ts.oninput = () => { const m = single(); if (!m) return; S.cut[m.key] = THRcur[+ts.value]; S.sel = null; draw(); };
  let THRcur = [0];
  const lg = el('div', { class: 'legend' }); sec.append(lg);
  lg.innerHTML = Object.values(NW_DOM).map(([c, l]) => `<span><i style="background:var(${c});height:10px;width:10px;border-radius:50%"></i>${l}</span>`).join('') +
    `<span><i style="background:none;border:1.5px dashed var(--ink-2);height:9px;width:9px;border-radius:50%"></i>doubtful MAG</span>` +
    `<span><i style="background:var(--s1);outline:1.5px dotted var(--ink);outline-offset:1px;height:9px;width:9px;border-radius:50%"></i>contamination 10–20%</span>` +
    `<span><i style="background:var(--s1);opacity:.45;height:8px;width:8px;border-radius:50%"></i>fragment (under 50% complete)</span>` +
    `<span><i style="background:none;border-top:1.5px dashed var(--ink-2);height:0;width:16px;border-radius:0"></i>possible same genome</span>` +
    `<span><i style="background:none;border:3px double var(--ink-2);height:11px;width:11px;border-radius:50%"></i>several genomes together</span>` +
    `<span><i style="background:var(--ink-2)"></i>positive link</span><span><i style="background:var(--s7)"></i>negative link</span>`;
  const wl = el('div', { class: 'legend', style: 'margin-top:2px' }); sec.append(wl);     // line-width legend
  const box = fig(sec, '', `Point size follows median abundance or, with a sample picked, the genome’s level in that sample (RNA: reads per kb per million; nanopore DNA: parts per million of the aligned bases; Illumina DNA: per million mapped read pairs), on one square-root scale for all samples of the layer; genomes not detected in the sample are small hollow dots. ◀ ▶ or the arrow keys step through the samples in time order (winter, then control and diesel pool, seawater before sea ice, then day). Layouts: “t-SNE” is a t-SNE of each genome’s profile across libraries (fourth root of its share of each library), scaled so every genome’s profile sums to 1. Genomes that rise and fall in the same libraries sit close together, whether they are common or rare. The other layouts are force-directed: each method’s positive links pull their genomes together, more strongly for higher scores, and all genomes push each other apart. Genomes without a positive link in that method sit on the ring around it. With several methods switched on, line width grows with the number of supporting methods; with exactly one, width and darkness follow that method’s score, and the slider sets the lowest score shown instead of the permissivity step. Hover a point for its lineage and partners; click it to show only its links, click the background to clear. Strong scores: SPRING and SPIEC-EASI find the link in at least ${P.strong_stab * 100}% of their subsamples; the residual Spearman has |ρ| at least ${P.strong_rho} and q below ${P.strong_q}; the Presence Pearson has |r| at least ${P.strong_pp_r}, q below ${P.strong_q} and both genomes present together in at least ${P.strong_pp_n} libraries; FlashWeave counts whenever it reports a link. Links found only by the residual Spearman are included down to |ρ| ${P.sp_floor} (q below 0.05).`);
  const info = el('p', { class: 'cap', style: 'margin-top:4px' }); box.parentNode.insertBefore(info, box.nextSibling);
  const short = id => id.startsWith('contig_') ? 'c' + id.slice(7).replace(/^0+/, '') : id;
  const layerLabel = k => J.layer_list.find(l => l.key === k).label;

  // line style: by the number of supporting methods, or, with exactly one method switched on, by that method's score class
  const single = () => nwPop(S.bits) === 1 ? M.find(m => m.bit === S.bits) : null;
  const scoreClass = (m, E, i) => { const v = E[m.col][i] / 100; let k = 0; m.bins.forEach(b => { if (v >= b) k++; }); return k; };   // 0..bins.length
  function widthLegend(m, op) {
    const sw = (w, o) => `<svg width="26" height="10" aria-hidden="true"><line x1="1" y1="5" x2="25" y2="5" stroke="${css('--ink-2')}" stroke-width="${w}" stroke-opacity="${o}"/></svg>`;
    if (m) { const b = m.bins, lims = [`< ${b[0]}`, ...b.slice(0, -1).map((x, k) => `${x}–${b[k + 1]}`), `≥ ${b[b.length - 1]}`];
      wl.innerHTML = `<span style="color:var(--ink-2)">Line: ${m.score}</span>` + lims.map((t, k) => `<span>${sw(0.4 + 0.55 * k, Math.min(1, 0.35 + 0.16 * k))}${t}</span>`).join(''); }
    else { const n = nwPop(S.bits);
      wl.innerHTML = `<span style="color:var(--ink-2)">Line: methods supporting the link</span>` + d3.range(1, n + 1).map(k => `<span>${sw(0.4 + 0.55 * (k - 1), Math.max(op, 0.5))}${k}</span>`).join(''); }
  }

  let G = null;      // what is drawn now: svg, scales and selections, so a layout change can animate
  const pxy = (L, key, x, y) => { const p = L.layouts[key]; return [p.x.map(x), p.y.map(y)]; };
  function scales(L, key, W, H, pad) {
    const p = L.layouts[key];
    return [d3.scaleLinear().domain(d3.extent(p.x)).range([pad, W - pad - 60]), d3.scaleLinear().domain(d3.extent(p.y)).range([H - pad, pad])];
  }

  async function draw() {
    const L = await nwLayer(S.layer), N = L.nodes, E = L.edges, step = steps[S.step], m1 = single();
    J.layers[S.layer].data = L;
    permBox.style.display = m1 ? 'none' : ''; cutBox.style.display = m1 ? '' : 'none';
    let cut = 0;
    if (m1) { THRcur = thresholds(L, m1); let k = 0; const c = S.cut[m1.key] ?? THRcur[0]; THRcur.forEach((t, j) => { if (t <= c) k = j; });
      ts.max = THRcur.length - 1; ts.value = k; cut = THRcur[k];
      lab.textContent = `${m1.score} ≥ ${d3.format('.2f')(cut / 100)}${k === 0 ? ' (every link it reports)' : ''}`; }
    else lab.textContent = `Step ${S.step + 1} of ${steps.length}: ${step.label}`;
    box.innerHTML = '';
    const inView = (a, b) => {
      if (S.view === 'grazer') { const g = n => n.role === 'grazer', t = n => n.domain === 'prokaryote' || n.role === 'alga' || n.domain === 'organelle'; return (g(a) && t(b)) || (g(b) && t(a)); }
      if (S.view === 'virus') return (a.domain === 'virus') !== (b.domain === 'virus');
      if (S.view === 'plastid') return (a.domain === 'organelle' && b.domain === 'eukaryote') || (b.domain === 'organelle' && a.domain === 'eukaryote');
      return true; };
    const K = [];      // [a, b, sign, supporting-method mask, row, line class]
    for (let i = 0; i < E.a.length; i++) {
      // one method on: its own call and the score cutoff only (no steps); otherwise the permissivity step
      const p = m1 ? (((E.m[i] | E.x[i]) & m1.bit) && E[m1.col][i] >= cut ? (E.x[i] & m1.bit ? -m1.bit : m1.bit) : 0) : nwPass(E, i, step, S.bits, designBits), sup = Math.abs(p);
      const sgn = p < 0 ? -E.s[i] : E.s[i];
      if (sup && (S.sign === 'both' || (S.sign === 'pos') === (sgn > 0)) && !(S.hideSame && E.same[i]) && inView(N[E.a[i]], N[E.b[i]]))
        K.push([E.a[i], E.b[i], sgn, sup, i, m1 ? scoreClass(m1, E, i) : nwPop(sup) - 1, E.same[i]]); }
    const nbr = new Map(N.map((_, i) => [i, []])); K.forEach(e => { nbr.get(e[0]).push(e); nbr.get(e[1]).push(e); });
    const drv = new Map(); L.drivers.forEach(([i, d, w]) => { (drv.get(i) || drv.set(i, []).get(i)).push([d, w]); });
    const W = width(box), H = Math.max(420, Math.min(640, W * 0.62)), pad = 18;
    const [x, y] = scales(L, S.layout, W, H, pad), [X, Y] = pxy(L, S.layout, x, y);
    const rs = d3.scaleSqrt().domain([0, d3.max(N, n => n.ab || 0) || 1]).range([2.2, 13]);
    const Sm = S.sample != null ? await smpData(S.layer) : null, sv = Sm ? Sm.v[S.sample] : null;
    const rsS = Sm ? d3.scaleSqrt().domain([0, Sm.max]).range([2.2, 15]).clamp(true) : null;
    const rad = i => sv ? (sv[i] > 0 ? rsS(sv[i]) : 1.8) : rs(N[i].ab || 0), zero = i => sv && !(sv[i] > 0);
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', H).attr('role', 'img').attr('aria-label', `${layerLabel(S.layer)} genome network`);
    svg.append('rect').attr('width', W).attr('height', H).attr('fill', 'transparent').on('click', () => { S.sel = null; draw(); });
    const sel = S.sel, selSet = sel == null ? null : new Set([sel, ...nbr.get(sel).map(e => e[0] === sel ? e[1] : e[0])]);
    const Es = sel == null ? K : K.filter(e => e[0] === sel || e[1] === sel);
    const viewNodes = S.view === 'all' ? null : new Set(K.flatMap(e => [e[0], e[1]]));
    // edges as one path per sign x line class (fast for tens of thousands of links)
    const op = sel == null ? (Es.length > 8000 ? 0.06 : Es.length > 2500 ? 0.12 : Es.length > 800 ? 0.22 : 0.4) : 0.85;
    // single method: opacity rises with the score class, from a third of the base to (nearly) full; otherwise the base
    const opk = k => m1 ? Math.min(1, op * (0.35 + 0.65 * k / m1.bins.length) * (sel == null ? 2.2 : 1)) : op;
    const grp = d3.group(Es, e => e[2] > 0, e => e[5], e => e[6]);
    const lines = [];
    grp.forEach((byW, posi) => byW.forEach((bySame, k) => bySame.forEach((es, same) => lines.push({ posi, k, same, es }))));
    const segs = (es, X, Y) => es.map(e => `M${X[e[0]].toFixed(1)},${Y[e[0]].toFixed(1)}L${X[e[1]].toFixed(1)},${Y[e[1]].toFixed(1)}`).join('');
    const paths = svg.append('g').style('pointer-events', 'none').selectAll('path').data(lines).join('path').attr('d', d => segs(d.es, X, Y)).attr('fill', 'none')
      .attr('stroke', d => css(d.posi ? '--ink-2' : '--s7')).attr('stroke-width', d => 0.4 + 0.55 * d.k).attr('stroke-opacity', d => d.same ? Math.max(0.6, opk(d.k)) : opk(d.k))
      .attr('stroke-dasharray', d => d.same ? '4 3' : null);
    const order = N.map((n, i) => i).sort((a, b) => (N[b].ab || 0) - (N[a].ab || 0));
    const circles = svg.append('g').selectAll('circle').data(order).join('circle')
      .attr('cx', i => X[i]).attr('cy', i => Y[i]).attr('r', rad)
      .attr('fill', i => N[i].status.includes('doubtful') || zero(i) ? css('--panel') : css(NW_DOM[N[i].domain][0]))
      .attr('stroke', i => N[i].status.includes('doubtful') || zero(i) ? css(NW_DOM[N[i].domain][0]) : N[i].status.includes('multi-genome') || N[i].status.includes(NW_CONTAM) ? css('--ink') : css('--panel'))
      .attr('stroke-width', i => zero(i) ? 1 : N[i].status.includes('multi-genome') ? 3 : N[i].status.includes('doubtful') || N[i].status.includes(NW_CONTAM) ? 1.6 : 0.6)
      .attr('stroke-dasharray', i => zero(i) ? null : N[i].status.includes('doubtful') ? '2 1.5' : N[i].status.includes(NW_CONTAM) ? '1 1.5' : null)
      .attr('fill-opacity', i => N[i].status.includes(NW_FRAG) ? 0.45 : 1)
      .attr('opacity', i => selSet ? (selSet.has(i) ? 1 : 0.12) : viewNodes ? (viewNodes.has(i) ? 1 : 0.15) : (N[i].domain === 'prokaryote' ? 0.8 : 1))
      .style('cursor', 'pointer')
      .on('mousemove', (ev, i) => showTip(ev, nwTip(N, i, nbr.get(i), drv.get(i), E)))
      .on('mouseleave', hideTip)
      .on('click', (ev, i) => { ev.stopPropagation(); S.sel = S.sel === i ? null : i; draw(); });
    // labels: eukaryotes, viruses and plastids with links, plus featured genera (network_site.FEATURED) and found genomes always
    const must = i => N[i].feat || (S.find && S.find.has(i)) || i === sel;
    const labs = order.filter(i => must(i) || (N[i].domain !== 'prokaryote' && nbr.get(i).length && (!selSet || selSet.has(i)) && (!viewNodes || viewNodes.has(i))));
    const labSet = [...labs.filter(must), ...labs.filter(i => !must(i)).slice(0, sel == null ? 40 : 80)];
    const ltext = i => N[i].domain === 'virus' ? N[i].id.replace('contig_', 'c') : N[i].domain === 'prokaryote' ? `${N[i].genus || N[i].label} (${N[i].id})` : N[i].id;
    const texts = svg.append('g').selectAll('text').data(labSet).join('text').attr('x', i => X[i] + rad(i) + 2).attr('y', i => Y[i] + 3.5)
      .style('font-size', '10px').style('font-weight', i => must(i) && N[i].domain === 'prokaryote' ? 600 : null).style('pointer-events', 'none').attr('fill', i => must(i) ? css('--ink') : null).text(ltext);
    const rings = svg.append('g').style('pointer-events', 'none').selectAll('circle').data(S.find ? [...S.find] : []).join('circle')
      .attr('cx', i => X[i]).attr('cy', i => Y[i]).attr('r', i => rad(i) + 6).attr('fill', 'none').attr('stroke', css('--s2') || css('--ink')).attr('stroke-width', 2.2);
    G = { L, W, H, pad, X, Y, paths, circles, texts, rings, segs, rad, n_edges: Es.length };
    // size legend
    if (Sm) { const b = Sm.libs[S.sample], ref = [Sm.max / 100, Sm.max / 10, Sm.max].map(v => +d3.format('.1~r')(v));
      sinfo.innerHTML = `${S.sample + 1} of ${Sm.libs.length} · <span class="mono">${b.id}</span>`;
      const lgS = svg.append('g').attr('transform', `translate(${W - 58},${pad + 4})`).style('pointer-events', 'none');
      lgS.append('text').attr('x', 0).attr('y', 0).style('font-size', '9.5px').text(Sm.unit.split(' ').slice(0, 3).join(' '));
      let yy = 14; ref.forEach(v => { const r = rsS(v); yy += r; lgS.append('circle').attr('cx', 15).attr('cy', yy).attr('r', r).attr('fill', 'none').attr('stroke', css('--ink-2'));
        lgS.append('text').attr('x', 33).attr('y', yy + 3).style('font-size', '9.5px').text(d3.format('~g')(v)); yy += r + 4; });
      lgS.append('circle').attr('cx', 15).attr('cy', yy + 3).attr('r', 1.8).attr('fill', 'none').attr('stroke', css('--ink-2'));
      lgS.append('text').attr('x', 33).attr('y', yy + 6).style('font-size', '9.5px').text('0'); }
    else sinfo.textContent = '';
    if (fl.dataset.layer !== S.layer) { fillFind(L); fl.dataset.layer = S.layer; }
    widthLegend(m1, op);
    const npos = K.filter(e => e[2] > 0).length;
    const nsame = K.filter(e => e[6]).length;
    cnt.textContent = `${fmt(K.length)} links shown (${fmt(npos)} positive, ${fmt(K.length - npos)} negative) among ${N.length} genomes${nsame ? `; ${fmt(nsame)} possible same-genome` : ''}`;
    info.textContent = sel == null ? '' : `${N[sel].label}: ${nbr.get(sel).length} link${nbr.get(sel).length === 1 ? '' : 's'} shown.`;
  }

  // layout change: glide every point (and its links) from the old to the new positions, then redraw
  let anim = null;
  function moveTo(key) {
    if (anim) { anim.stop(); anim = null; }
    if (!G || matchMedia('(prefers-reduced-motion: reduce)').matches) { S.layout = key; draw(); return; }
    const g = G, [x, y] = scales(g.L, key, g.W, g.H, g.pad), [X1, Y1] = pxy(g.L, key, x, y), X0 = g.X, Y0 = g.Y, n = X0.length;
    const X = new Array(n), Y = new Array(n), heavy = g.n_edges > 12000, ease = d3.easeCubicInOut, dur = 750;
    if (heavy) g.paths.attr('display', 'none');     // too many links to redraw every frame; they come back at the end
    anim = d3.timer(t => {
      const u = ease(Math.min(1, t / dur));
      for (let i = 0; i < n; i++) { X[i] = X0[i] + (X1[i] - X0[i]) * u; Y[i] = Y0[i] + (Y1[i] - Y0[i]) * u; }
      g.circles.attr('cx', i => X[i]).attr('cy', i => Y[i]); g.rings.attr('cx', i => X[i]).attr('cy', i => Y[i]);
      g.texts.attr('x', i => X[i] + g.rad(i) + 2).attr('y', i => Y[i] + 3.5);
      if (!heavy) g.paths.attr('d', d => g.segs(d.es, X, Y));
      if (t >= dur) { anim.stop(); anim = null; S.layout = key; draw(); }
    });
  }
  draw();
}

function nwTip(N, i, edges, drv, E) {
  const n = N[i], lin = D.lineage && D.lineage[n.id] ? D.lineage[n.id] : '';
  const st = [n.ee_status && n.ee_status !== 'network' ? n.ee_status : '', ...n.status.filter(s => s !== 'doubtful' || !n.ee_status), n.group].filter(Boolean).join(' · ');
  const P = edges.map(e => { const j = e[0] === i ? e[1] : e[0]; return { j, s: e[2], k: nwPop(e[3]), r: e[4] }; })
    .sort((a, b) => b.k - a.k || (N[b.j].domain !== 'prokaryote') - (N[a.j].domain !== 'prokaryote')).slice(0, 7);
  const DL = { winter_2025_26: 'winter 2025-26', diesel_pool: 'diesel pool', control_pool: 'control pool', seawater: 'seawater', sea_ice: 'sea ice', top_ice: 'top ice', bottom_ice: 'bottom ice', day: 'day', under_ice_par: 'under-ice light' };
  return `<b>${n.label}</b>${lin ? `<br><span style="opacity:.8">${lin}</span>` : ''}${st ? `<br>${st}` : ''}<br>median ${n.ab == null ? '–' : fmt1(n.ab)} per million · detected in ${d3.format('.0%')(n.prev)} of libraries` +
    (drv && drv.length ? `<br>FlashWeave drivers: ${drv.map(([d, w]) => `${DL[d] || d} (${w > 0 ? '+' : '−'})`).join(', ')}` : '') +
    (P.length ? `<br><b>Partners</b> (methods supporting; residual Spearman ρ; Presence Pearson r):<br>${P.map(p => `${p.s > 0 ? '+' : '−'} ${N[p.j].label} (${p.k}${E && E.rho[p.r] ? `; ρ ${(E.rho[p.r] / 100).toFixed(2)}` : ''}${E && E.pr[p.r] ? `; r ${(E.pr[p.r] / 100).toFixed(2)}` : ''})${E && E.same[p.r] ? ' · possible same genome' : ''}`).join('<br>')}${edges.length > P.length ? `<br>… ${edges.length - P.length} more` : ''}` : '<br>No links in this view.');
}

function nwPartners(sec, J) {
  sec.append(el('h3', {}, 'Top partners of the eukaryotes and the giant viruses'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const s = el('select', { id: 'nwpart', 'aria-label': 'Genomes' }); [['eukaryote', 'Eukaryotic MAGs'], ['virus', 'Giant viruses and the virophage']].forEach(([k, l]) => s.append(new Option(l, k)));
  c.append(el('label', { for: 'nwpart' }, 'Genomes'), s);
  const box = fig(sec, '', `Consensus links (slider step ${J.default_step + 1}: at least two methods, one of them FlashWeave or the residual Spearman), in any layer; up to eight partners each, those found in the most layers first. Methods: ${J.methods.map(m => `${m.abbr} = ${m.label}`).join(', ')}. Click a row to show that genome in the network above.`);
  const ab = Object.fromEntries(J.methods.map(m => [m.key, m.abbr]));
  const ml = l => l ? l.map(m => ab[m]).join('') : '–', LL = J.layer_list.slice().reverse();
  function draw() {
    const P = J.partners.filter(p => p.domain === s.value);
    box.innerHTML = `<div class="tablewrap" style="max-height:640px"><table><thead><tr><th>Genome</th><th>Partner</th><th class="num">Sign</th>${LL.map(l => `<th>${l.label}</th>`).join('')}</tr></thead><tbody>${
      P.map(p => p.partners.map((q, k) => `<tr data-id="${p.id}"><td style="white-space:normal;min-width:200px">${k === 0 ? `<b>${p.label}</b>${p.role ? `<div style="font-size:11.5px;color:var(--ink-3)">${p.role}; ${p.n} partners</div>` : `<div style="font-size:11.5px;color:var(--ink-3)">${p.n} partners</div>`}` : ''}</td><td style="white-space:normal;min-width:240px"><span style="color:var(${NW_DOM[q.domain] ? NW_DOM[q.domain][0] : '--ink-2'})">●</span> ${q.label}</td><td class="num">${q.sign > 0 ? '+' : '−'}</td>${LL.map(l => `<td class="mono">${ml(q.layers[l.key])}</td>`).join('')}</tr>`).join('')).join('')}</tbody></table></div>`;
    box.querySelectorAll('tr[data-id]').forEach(tr => { tr.onclick = () => nwFocus(J, tr.dataset.id); });
  }
  s.onchange = draw; draw();
}

function nwFocus(J, id) {
  // select the genome in the network (in the first layer that has it, in page order) and scroll there
  const sec = document.getElementById('v-network'); if (!sec) return;
  const li = J.layer_list.findIndex(l => NW_ROWS[l.key].nodes.some(n => n.id === id)); if (li < 0) return;
  const L = J.layer_list[li].key, btn = sec.querySelectorAll('button[data-g="layer"]')[li]; if (btn) btn.click();
  const i = NW_ROWS[L].nodes.findIndex(n => n.id === id);
  setTimeout(() => { const c = [...sec.querySelectorAll('.fig svg circle')].find(e => d3.select(e).datum() === i);
    if (c) { const f = c.closest('.fig'); c.dispatchEvent(new MouseEvent('click', { bubbles: true })); if (f) f.scrollIntoView({ behavior: 'smooth', block: 'center' }); } }, 80);
}

function nwMethods(sec, J) {
  sec.append(el('h3', {}, 'How the methods compare'));
  sec.append(el('p', { class: 'kicker' }, `FlashWeave tests each pair while holding other genomes and the sample drivers (winter, pool, seawater or ice, ice layer, day, modelled under-ice light) fixed, so a link that only reflects a shared driver drops out. SPRING and SPIEC-EASI estimate direct links among the genomes from their relative abundances, but know nothing about the sampling design. The residual Spearman baseline removes the design first and then correlates what is left, pair by pair; it flags many pairs because it does not separate direct from indirect links. The Presence Pearson correlates the logarithm of two genomes’ shares, using only the libraries where both are present; absences do not count for or against a link. It reports a link at |r| ${J.params.pp_r} or more, with both genomes present together in at least ${J.params.pp_n} libraries and q below ${J.params.pp_q}. It does not know the design either. The consensus therefore needs two methods, one of them FlashWeave or the baseline.`));
  const box = fig(sec, 'Links per method and overlap', 'Overlap = shared links ÷ links found by either (Jaccard index); sign agreement = share of the shared links with the same sign.');
  const lab = Object.fromEntries(J.methods.map(m => [m.key, m.label])); lab.consensus = 'Consensus';
  box.innerHTML = J.layer_list.slice().reverse().map(({ key: l, label }) => { const L = J.layers[l];
    return `<div class="title" style="margin-top:6px">${label} network: ${L.n_nodes} genomes (${Object.entries(L.by_domain).map(([d, n]) => `${n} ${NW_DOM[d][1].toLowerCase()}`).join(', ')}), ${L.n_libs} libraries</div>
    <div class="tablewrap" style="max-height:none"><table><thead><tr><th>Method</th><th class="num">Links</th><th class="num">Positive</th><th class="num">In consensus</th></tr></thead><tbody>${
      L.summary.map(s => `<tr><td>${lab[s.method]}</td><td class="num">${fmt(s.edges)}</td><td class="num">${fmt(s.positive)}</td><td class="num">${s.method === 'consensus' ? '–' : fmt(s.in_consensus)}</td></tr>`).join('')}</tbody></table></div>
    <div class="tablewrap" style="max-height:none;margin-top:8px"><table><thead><tr><th>Pair of methods</th><th class="num">Shared links</th><th class="num">Overlap</th><th class="num">Sign agreement</th></tr></thead><tbody>${
      L.overlap.map(o => `<tr><td>${lab[o.m1]} · ${lab[o.m2]}</td><td class="num">${fmt(o.shared)}</td><td class="num">${d3.format('.3f')(o.jaccard)}</td><td class="num">${o.sign_agreement == null ? '–' : d3.format('.0%')(o.sign_agreement)}</td></tr>`).join('')}</tbody></table></div>`; }).join('');
}

function nwReading(sec, J) {
  sec.append(el('h3', {}, 'How to read it'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li><b>Genomes.</b> Bacterial and archaeal MAGs with CheckM2 completeness at least 15% and contamination under 20%. Those under 50% complete are drawn pale (“fragment”); those with 10–20% contamination have a dotted ring (most likely a second strain in the same bin, as in the Cycloclasticus bloom genomes lorbin_184 and metabat_5432). Two MAGs of the same genus whose RNA profiles agree almost perfectly (Spearman ρ at least 0.95) may be pieces of one genome; their link is dashed and can be hidden; every eukaryotic MAG (EE##), with doubtful ones drawn hollow; the viruses of the Viruses page plus every giant virus and virophage averaging at least 5 transcripts per million; plastid group P01 and the plastid MAGs. MAGs that hold several genomes (P01, and eukaryotic MAGs with two strains or species) are drawn with a double ring. A genome enters a layer when it is detected in at least a quarter of that layer’s libraries; eukaryotes, giant viruses, the virophage and P01 need only a tenth. Each contig counts for one genome only, viruses first, then plastid groups, eukaryotic MAGs and bacterial MAGs.</li>
    <li><b>DNA and RNA.</b> The nanopore DNA network has the fewest samples (${J.layers.dna.n_libs} pool samples, replicate runs summed), so it is the least certain. The Illumina DNA network uses ${J.layers.mgdna ? J.layers.mgdna.n_libs : ''} pool metagenomes (control and diesel pool, seawater and sea ice, both winters); a genome’s level is its share of the mapped reads, the same measure as for nanopore. The RNA network has ${J.layers.rna.n_libs} libraries. DNA shows who is there, RNA who is active, so a link can appear in one and not the other.</li>
    <li><b>Sign.</b> A positive link means the two genomes rise together after the drivers are accounted for; a negative link means one rises when the other falls. Grazing and infection can produce either, depending on timing.</li>
    <li><b>Limits.</b> One pool per treatment, so samples within a pool are not independent. Relative abundance is a share of all reads, so a large bloom pushes everything else down. Links between genomes that share contigs or reads (for example a virus and the genome it is integrated in) can be artefacts.</li>`));
}

RENDER.network = networkSec;
