/* Network section (Food web): genome-genome association networks per layer (DNA, RNA), consensus of FlashWeave, SPRING,
   SPIEC-EASI, a design-residual Spearman baseline and a presence-only Pearson. Data: network.json (network/scripts/network_site.py). Helpers from app.js. */
const NW_DOM = { prokaryote: ['--s1', 'Bacteria and archaea'], eukaryote: ['--s3', 'Eukaryotes'], virus: ['--s8', 'Viruses'], organelle: ['--s4', 'Plastids'] };
const NW_VIEWS = [['all', 'All links'], ['grazer', 'Grazers with bacteria and algae'], ['virus', 'Viruses with their partners'], ['plastid', 'Plastids with eukaryotes']];
const nwPop = m => { let c = 0; while (m) { c += m & 1; m >>= 1; } return c; };

const NW_ROWS = {};      // per-layer data (network_{layer}.json), loaded on first use
async function nwLayer(l) {
  if (NW_ROWS[l]) return NW_ROWS[l];
  const L = await (await fetch(`network_${l}.json`, { cache: 'no-cache' })).json(), E = L.edges, n = E.a.length, M = D.network.methods, nb = M.length;
  // unpack the flags (bit 0 positive sign, then the masks of agreeing, opposing and strong-scoring methods) and expand the sparse
  // score columns (one value per edge the method reports, in edge order)
  const full = (1 << nb) - 1;
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
  const kn = J.known, knHit = new Set(kn.filter(k => k.consensus).map(k => k.a + '|' + k.b));
  sec.innerHTML = `<h2>Network</h2>
  <p class="kicker">Which genomes rise and fall together? Each point is a genome: a bacterial or archaeal MAG, a eukaryotic MAG, a virus or a plastid. There are two networks. The DNA network uses nanopore metagenomes and asks who is there together; the RNA network uses the ${J.layers.rna.n_libs} Illumina RNA libraries (metatranscriptomes) and asks who is active together. A genome’s RNA level is the typical level of its contigs, so a few very busy contigs (such as rRNA genes) cannot set it. Five methods look for links, and the slider sets how much agreement a link needs: from four or five methods with strong scores to any single method. The default is the consensus: at least two methods agree on the sign, and one of them has already accounted for winter, pool, sea ice or seawater, and time. A link means two genomes co-vary more than these shared drivers explain. It is a lead for grazing, infection or a shared niche, not proof of any of them.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  const LL = J.layer_list.slice().reverse(), names = LL.map(l => l.label).join(' / ');
  [[LL.map(l => J.layers[l.key].n_nodes).join(' / '), `genomes in the ${names} network`],
   [LL.map(l => fmt(cons(l.key).edges)).join(' / '), `consensus links, ${names}`],
   [`${knHit.size} of ${new Set(kn.map(k => k.a + '|' + k.b)).size}`, 'expected virus–host and plastid–host pairs in the consensus of at least one layer']]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  nwGraph(sec, J);
  nwPartners(sec, J);
  nwKnown(sec, J);
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
  const S = { layer: J.layer_list[0].key, view: 'all', bits: allBits, step: J.default_step, sel: null, layout: 'tsne' };
  const c1 = el('div', { class: 'controls' }); sec.append(c1);
  const chip = (parent, opts, key) => opts.forEach(([k, l]) => { const b = el('button', { class: 'chip', 'aria-pressed': S[key] === k }, l);
    b.onclick = () => { S[key] = k; S.sel = null; parent.querySelectorAll(`button[data-g="${key}"]`).forEach(x => x.setAttribute('aria-pressed', x === b)); draw(); };
    b.dataset.g = key; parent.append(b); });
  chip(c1, J.layer_list.map(l => [l.key, l.long]), 'layer');
  const vs = el('select', { id: 'nwview', 'aria-label': 'View' }); NW_VIEWS.forEach(([k, l]) => vs.append(new Option(l, k)));
  vs.onchange = () => { S.view = vs.value; S.sel = null; draw(); }; c1.append(el('label', { for: 'nwview' }, 'Show'), vs);
  const ls = el('select', { id: 'nwlayout', 'aria-label': 'Layout' }); J.layouts.forEach(o => ls.append(new Option(o.key === 'tsne' ? 't-SNE (abundance profiles)' : `${o.label} links`, o.key)));
  ls.onchange = () => moveTo(ls.value); c1.append(el('label', { for: 'nwlayout' }, 'Layout'), ls);
  const c2 = el('div', { class: 'controls' }); sec.append(c2);
  c2.append(el('span', { style: 'color:var(--ink-2)' }, 'Methods'));
  M.forEach(m => { const id = 'nwm-' + m.key, cb = el('input', { type: 'checkbox', id }); cb.checked = true;
    cb.onchange = () => { S.bits = cb.checked ? S.bits | m.bit : S.bits & ~m.bit; draw(); }; c2.append(cb, el('label', { for: id, style: 'margin-right:6px' }, m.label)); });
  // permissivity slider + known-pair benchmark
  const P = J.params, steps = J.steps;
  const pw = el('div', { class: 'controls', style: 'align-items:flex-start;gap:18px;flex-wrap:wrap' }); sec.append(pw);
  const left = el('div', { style: 'flex:1 1 340px;min-width:300px' }); pw.append(left);
  const sl = el('input', { type: 'range', id: 'nwperm', min: 1, max: steps.length, step: 1, value: S.step + 1, style: 'width:100%', 'aria-label': 'Permissivity', list: 'nwperm-ticks' });
  const dl = el('datalist', { id: 'nwperm-ticks' }); steps.forEach((_, k) => dl.append(new Option('', k + 1)));
  const lab = el('div', { style: 'font-weight:600;margin-top:2px' }), cnt = el('div', { style: 'color:var(--ink-2);font-size:13px' });
  left.append(el('label', { for: 'nwperm', style: 'display:flex;justify-content:space-between;font-size:12.5px;color:var(--ink-2)' }, '<span>strict</span><span>Permissivity</span><span>permissive</span>'), sl, dl, lab, cnt);
  sl.oninput = () => { S.step = +sl.value - 1; S.sel = null; draw(); };
  const bench = el('div', { style: 'flex:1 1 360px;min-width:300px' }); pw.append(bench);
  const lg = el('div', { class: 'legend' }); sec.append(lg);
  lg.innerHTML = Object.values(NW_DOM).map(([c, l]) => `<span><i style="background:var(${c});height:10px;width:10px;border-radius:50%"></i>${l}</span>`).join('') +
    `<span><i style="background:none;border:1.5px dashed var(--ink-2);height:9px;width:9px;border-radius:50%"></i>doubtful MAG</span>` +
    `<span><i style="background:none;border:3px double var(--ink-2);height:11px;width:11px;border-radius:50%"></i>several genomes together</span>` +
    `<span><i style="background:var(--ink-2)"></i>positive link</span><span><i style="background:var(--s7)"></i>negative link</span>`;
  const wl = el('div', { class: 'legend', style: 'margin-top:2px' }); sec.append(wl);     // line-width legend
  const box = fig(sec, '', `Point size follows median abundance. Layouts: “t-SNE” places genomes by their abundance across the libraries (fourth root of the share of each library), so genomes with similar profiles sit close together. The other layouts are force-directed: each method’s positive links pull their genomes together, more strongly for higher scores, and all genomes push each other apart. Genomes without a positive link in that method sit on the ring around it. With several methods switched on, line width grows with the number of supporting methods; with exactly one, width and darkness follow that method’s score. Hover a point for its lineage and partners; click it to show only its links, click the background to clear. Strong scores: SPRING and SPIEC-EASI find the link in at least ${P.strong_stab * 100}% of their subsamples; the residual Spearman has |ρ| at least ${P.strong_rho} and q below ${P.strong_q}; the Presence Pearson has |r| at least ${P.strong_pp_r}, q below ${P.strong_q} and both genomes present together in at least ${P.strong_pp_n} libraries; FlashWeave counts whenever it reports a link. Links found only by the residual Spearman are included down to |ρ| ${P.sp_floor} (q below 0.05).`);
  const info = el('p', { class: 'cap', style: 'margin-top:4px' }); box.parentNode.insertBefore(info, box.nextSibling);
  const kp = [...new Map(J.known.map(k => [k.a + '|' + k.b, k])).values()];
  const short = id => id.startsWith('contig_') ? 'c' + id.slice(7).replace(/^0+/, '') : id;
  const layerLabel = k => J.layer_list.find(l => l.key === k).label;

  function benchTable(L) {
    // known pairs recovered at every step, for the current layer and method selection
    const E = L.edges, n = E.a.length, ix = new Map(L.nodes.map((d, i) => [d.id, i]));
    const pk = kp.map(p => { const a = ix.get(p.a), b = ix.get(p.b); if (a == null || b == null) return -1;
      for (let i = 0; i < n; i++) if ((E.a[i] === a && E.b[i] === b) || (E.a[i] === b && E.b[i] === a)) return i; return -1; });
    const rows = steps.map((st, k) => { let c = 0; for (let i = 0; i < n; i++) if (nwPass(E, i, st, S.bits, designBits)) c++;
      return { k, c, hit: pk.map(i => i >= 0 && nwPass(E, i, st, S.bits, designBits) ? 1 : 0) }; });
    bench.innerHTML = `<div style="font-size:12.5px;color:var(--ink-2);margin-bottom:3px">Expected pairs recovered (${layerLabel(S.layer)}): ● found, · not found</div>
    <div class="tablewrap" style="max-height:none"><table style="font-size:12px"><thead><tr><th>Step</th><th class="num">Links</th>${kp.map(p => `<th title="${p.pair}" class="num">${short(p.a)}<br>${p.b}</th>`).join('')}<th class="num">Found</th></tr></thead><tbody>${
      rows.map(r => `<tr style="${r.k === S.step ? 'background:var(--panel-2,rgba(127,127,127,.15));font-weight:600' : ''}"><td>${r.k + 1}</td><td class="num">${fmt(r.c)}</td>${r.hit.map(h => `<td class="num" style="color:${h ? 'var(--s3)' : 'var(--ink-3)'}">${h ? '●' : '·'}</td>`).join('')}<td class="num">${d3.sum(r.hit)} / ${kp.length}</td></tr>`).join('')}</tbody></table></div>`;
  }

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
    lab.textContent = `Step ${S.step + 1} of ${steps.length}: ${step.label}`;
    benchTable(L);
    box.innerHTML = '';
    const inView = (a, b) => {
      if (S.view === 'grazer') { const g = n => n.role === 'grazer', t = n => n.domain === 'prokaryote' || n.role === 'alga' || n.domain === 'organelle'; return (g(a) && t(b)) || (g(b) && t(a)); }
      if (S.view === 'virus') return (a.domain === 'virus') !== (b.domain === 'virus');
      if (S.view === 'plastid') return (a.domain === 'organelle' && b.domain === 'eukaryote') || (b.domain === 'organelle' && a.domain === 'eukaryote');
      return true; };
    const K = [];      // [a, b, sign, supporting-method mask, row, line class]
    for (let i = 0; i < E.a.length; i++) { const p = nwPass(E, i, step, S.bits, designBits), sup = Math.abs(p);
      if (sup && inView(N[E.a[i]], N[E.b[i]])) K.push([E.a[i], E.b[i], p < 0 ? -E.s[i] : E.s[i], sup, i, m1 ? scoreClass(m1, E, i) : nwPop(sup) - 1]); }
    const nbr = new Map(N.map((_, i) => [i, []])); K.forEach(e => { nbr.get(e[0]).push(e); nbr.get(e[1]).push(e); });
    const drv = new Map(); L.drivers.forEach(([i, d, w]) => { (drv.get(i) || drv.set(i, []).get(i)).push([d, w]); });
    const W = width(box), H = Math.max(420, Math.min(640, W * 0.62)), pad = 18;
    const [x, y] = scales(L, S.layout, W, H, pad), [X, Y] = pxy(L, S.layout, x, y);
    const rs = d3.scaleSqrt().domain([0, d3.max(N, n => n.ab || 0) || 1]).range([2.2, 13]);
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', H).attr('role', 'img').attr('aria-label', `${layerLabel(S.layer)} genome network`);
    svg.append('rect').attr('width', W).attr('height', H).attr('fill', 'transparent').on('click', () => { S.sel = null; draw(); });
    const sel = S.sel, selSet = sel == null ? null : new Set([sel, ...nbr.get(sel).map(e => e[0] === sel ? e[1] : e[0])]);
    const Es = sel == null ? K : K.filter(e => e[0] === sel || e[1] === sel);
    const viewNodes = S.view === 'all' ? null : new Set(K.flatMap(e => [e[0], e[1]]));
    // edges as one path per sign x line class (fast for tens of thousands of links)
    const op = sel == null ? (Es.length > 8000 ? 0.06 : Es.length > 2500 ? 0.12 : Es.length > 800 ? 0.22 : 0.4) : 0.85;
    // single method: opacity rises with the score class, from a third of the base to (nearly) full; otherwise the base
    const opk = k => m1 ? Math.min(1, op * (0.35 + 0.65 * k / m1.bins.length) * (sel == null ? 2.2 : 1)) : op;
    const grp = d3.group(Es, e => e[2] > 0, e => e[5]);
    const lines = [];
    grp.forEach((byW, posi) => byW.forEach((es, k) => lines.push({ posi, k, es })));
    const segs = (es, X, Y) => es.map(e => `M${X[e[0]].toFixed(1)},${Y[e[0]].toFixed(1)}L${X[e[1]].toFixed(1)},${Y[e[1]].toFixed(1)}`).join('');
    const paths = svg.append('g').style('pointer-events', 'none').selectAll('path').data(lines).join('path').attr('d', d => segs(d.es, X, Y)).attr('fill', 'none')
      .attr('stroke', d => css(d.posi ? '--ink-2' : '--s7')).attr('stroke-width', d => 0.4 + 0.55 * d.k).attr('stroke-opacity', d => opk(d.k));
    const order = N.map((n, i) => i).sort((a, b) => (N[b].ab || 0) - (N[a].ab || 0));
    const circles = svg.append('g').selectAll('circle').data(order).join('circle')
      .attr('cx', i => X[i]).attr('cy', i => Y[i]).attr('r', i => rs(N[i].ab || 0))
      .attr('fill', i => N[i].status.includes('doubtful') ? css('--panel') : css(NW_DOM[N[i].domain][0]))
      .attr('stroke', i => N[i].status.includes('doubtful') ? css(NW_DOM[N[i].domain][0]) : N[i].status.includes('multi-genome') ? css('--ink') : css('--panel'))
      .attr('stroke-width', i => N[i].status.includes('multi-genome') ? 3 : N[i].status.includes('doubtful') ? 1.6 : 0.6)
      .attr('stroke-dasharray', i => N[i].status.includes('doubtful') ? '2 1.5' : null)
      .attr('opacity', i => selSet ? (selSet.has(i) ? 1 : 0.12) : viewNodes ? (viewNodes.has(i) ? 1 : 0.15) : (N[i].domain === 'prokaryote' ? 0.8 : 1))
      .style('cursor', 'pointer')
      .on('mousemove', (ev, i) => showTip(ev, nwTip(N, i, nbr.get(i), drv.get(i), E)))
      .on('mouseleave', hideTip)
      .on('click', (ev, i) => { ev.stopPropagation(); S.sel = S.sel === i ? null : i; draw(); });
    const labs = order.filter(i => N[i].domain !== 'prokaryote' && nbr.get(i).length && (!selSet || selSet.has(i)) && (!viewNodes || viewNodes.has(i)));
    const texts = svg.append('g').selectAll('text').data(labs.slice(0, sel == null ? 40 : 80)).join('text').attr('x', i => X[i] + rs(N[i].ab || 0) + 2).attr('y', i => Y[i] + 3.5)
      .style('font-size', '10px').style('pointer-events', 'none').text(i => N[i].domain === 'virus' ? N[i].id.replace('contig_', 'c') : N[i].id);
    G = { L, W, H, pad, X, Y, paths, circles, texts, segs, rs, n_edges: Es.length };
    widthLegend(m1, op);
    const npos = K.filter(e => e[2] > 0).length;
    cnt.textContent = `${fmt(K.length)} links shown (${fmt(npos)} positive, ${fmt(K.length - npos)} negative) among ${N.length} genomes`;
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
      g.circles.attr('cx', i => X[i]).attr('cy', i => Y[i]);
      g.texts.attr('x', i => X[i] + g.rs(g.L.nodes[i].ab || 0) + 2).attr('y', i => Y[i] + 3.5);
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
    (P.length ? `<br><b>Partners</b> (methods supporting; residual Spearman ρ; Presence Pearson r):<br>${P.map(p => `${p.s > 0 ? '+' : '−'} ${N[p.j].label} (${p.k}${E && E.rho[p.r] ? `; ρ ${(E.rho[p.r] / 100).toFixed(2)}` : ''}${E && E.pr[p.r] ? `; r ${(E.pr[p.r] / 100).toFixed(2)}` : ''})`).join('<br>')}${edges.length > P.length ? `<br>… ${edges.length - P.length} more` : ''}` : '<br>No links in this view.');
}

function nwPartners(sec, J) {
  sec.append(el('h3', {}, 'Top partners of the eukaryotes and the giant viruses'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const s = el('select', { id: 'nwpart', 'aria-label': 'Genomes' }); [['eukaryote', 'Eukaryotic MAGs'], ['virus', 'Giant viruses and the virophage']].forEach(([k, l]) => s.append(new Option(l, k)));
  c.append(el('label', { for: 'nwpart' }, 'Genomes'), s);
  const box = fig(sec, '', `Consensus links (slider step ${J.default_step + 1}: at least two methods, one of them FlashWeave or the residual Spearman), in either layer; up to eight partners each, those found in both layers first. Methods: ${J.methods.map(m => `${m.abbr} = ${m.label}`).join(', ')}. Click a row to show that genome in the network above.`);
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

function nwKnown(sec, J) {
  sec.append(el('h3', {}, 'Check against expected pairs'));
  const pairs = [...new Map(J.known.map(k => [k.a + '|' + k.b, k])).values()];
  const name = id => { for (const { key: l } of J.layer_list) { const n = NW_ROWS[l] && NW_ROWS[l].nodes.find(x => x.id === id); if (n) return n.label; } return id; };
  const ab = Object.fromEntries(J.methods.map(m => [m.key, m.label]));
  const cell = k => !k ? '–' : !k.in_layer ? 'not in this layer' : k.methods.length ? `${k.consensus ? '<b>found</b>' : 'one method'}: ${k.methods.map(m => ab[m] + (k.signs[m] < 0 ? ' (−)' : '')).join(', ')}` : 'not found';
  const box = el('div', { class: 'note' }), LL = J.layer_list.slice().reverse(); sec.append(box);
  box.innerHTML = `<p style="margin:0 0 6px">The Viruses and Eukaryotes pages matched giant viruses, the virophage and a plastid to likely eukaryotic hosts by simple correlation. These pairs were not used to tune the methods. This table shows the consensus (slider step ${J.default_step + 1}); the table beside the slider shows every step.</p>
  <div class="tablewrap" style="max-height:none;border:0"><table><thead><tr><th>Pair</th>${LL.map(l => `<th>${l.label}</th>`).join('')}</tr></thead><tbody>${
    pairs.map(p => `<tr><td style="white-space:normal;min-width:220px">${name(p.a)} – ${name(p.b)}<div style="font-size:11.5px;color:var(--ink-3)">${p.pair}</div></td>${
      LL.map(l => `<td style="white-space:normal">${cell(J.known.find(k => k.a === p.a && k.b === p.b && k.layer === l.key))}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`;
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
    <li><b>Genomes.</b> Bacterial and archaeal MAGs of high or medium quality (CheckM2 completeness at least 50%, contamination under 10%); every eukaryotic MAG (EE##), with doubtful ones drawn hollow; the viruses of the Viruses page plus every giant virus and virophage averaging at least 5 transcripts per million; plastid group P01 and the plastid MAGs. MAGs that hold several genomes (P01, and eukaryotic MAGs with two strains or species) are drawn with a double ring. A genome enters a layer when it is detected in at least a quarter of that layer’s libraries; eukaryotes, giant viruses, the virophage and P01 need only a tenth. Each contig counts for one genome only, viruses first, then plastid groups, eukaryotic MAGs and bacterial MAGs.</li>
    <li><b>DNA and RNA.</b> The DNA network has fewer samples (${J.layers.dna.n_libs} pool samples with nanopore metagenomes, replicate runs summed) than the RNA network (${J.layers.rna.n_libs} libraries), so it is less certain.</li>
    <li><b>Sign.</b> A positive link means the two genomes rise together after the drivers are accounted for; a negative link means one rises when the other falls. Grazing and infection can produce either, depending on timing.</li>
    <li><b>Limits.</b> One pool per treatment, so samples within a pool are not independent. Relative abundance is a share of all reads, so a large bloom pushes everything else down. Links between genomes that share contigs or reads (for example a virus and the genome it is integrated in) can be artefacts.</li>`));
}

RENDER.network = networkSec;
