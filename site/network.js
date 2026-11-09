/* Network section (Food web): genome-genome association networks in DNA and RNA, consensus of FlashWeave, SPRING, SPIEC-EASI
   and a design-residual Spearman baseline. Data: network.json (network/scripts/network_site.py). Helpers from app.js. */
const NW_DOM = { prokaryote: ['--s1', 'Bacteria and archaea'], eukaryote: ['--s3', 'Eukaryotes'], virus: ['--s8', 'Viruses'], organelle: ['--s4', 'Plastids'] };
const NW_VIEWS = [['all', 'All links'], ['grazer', 'Grazers with bacteria and algae'], ['virus', 'Viruses with their partners'], ['plastid', 'Plastids with eukaryotes']];
const nwPop = m => { let c = 0; while (m) { c += m & 1; m >>= 1; } return c; };

function networkSec(sec) {
  const J = D.network;
  if (!J) { sec.innerHTML = '<h2>Network</h2><p class="kicker">network.json is missing.</p>'; return; }
  const M = J.methods, allBits = M.reduce((a, m) => a | m.bit, 0), designBits = M.filter(m => m.design_aware).reduce((a, m) => a | m.bit, 0);
  const cons = l => J.layers[l].summary.find(s => s.method === 'consensus');
  const kn = J.known, knHit = new Set(kn.filter(k => k.consensus).map(k => k.a + '|' + k.b));
  sec.innerHTML = `<h2>Network</h2>
  <p class="kicker">Which genomes rise and fall together? Each point is a genome: a bacterial or archaeal MAG, a eukaryotic MAG, a virus or a plastid. There are two networks. The DNA network uses nanopore metagenomes and asks who is there together; the RNA network uses the ${J.layers.rna.n_libs} RNA libraries and asks who is active together. Points are placed by their abundance across samples, so genomes with similar profiles sit close together. A line joins two genomes when at least two of four methods find the link with the same sign, and at least one of those methods has already accounted for winter, pool, sea ice or seawater, and time. A link means two genomes co-vary more than these shared drivers explain. It is a lead for grazing, infection or a shared niche, not proof of any of them.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  [[`${J.layers.dna.nodes.length} / ${J.layers.rna.nodes.length}`, 'genomes in the DNA / RNA network'],
   [`${fmt(cons('dna').edges)} / ${fmt(cons('rna').edges)}`, 'links found by at least two methods, DNA / RNA'],
   [`${knHit.size} of ${new Set(kn.map(k => k.a + '|' + k.b)).size}`, 'expected virus–host and plastid–host pairs recovered in at least one layer']]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  nwGraph(sec, J, allBits, designBits);
  nwPartners(sec, J);
  nwKnown(sec, J);
  nwMethods(sec, J);
  nwReading(sec, J);
}

function nwGraph(sec, J, allBits, designBits) {
  sec.append(el('h3', {}, 'The networks'));
  const S = { layer: 'rna', view: 'all', bits: allBits, min: 2, design: true, sel: null };
  const c1 = el('div', { class: 'controls' }); sec.append(c1);
  const chip = (parent, opts, key) => opts.forEach(([k, l]) => { const b = el('button', { class: 'chip', 'aria-pressed': S[key] === k }, l);
    b.onclick = () => { S[key] = k; S.sel = null; parent.querySelectorAll(`button[data-g="${key}"]`).forEach(x => x.setAttribute('aria-pressed', x === b)); draw(); };
    b.dataset.g = key; parent.append(b); });
  chip(c1, [['rna', 'RNA (activity)'], ['dna', 'DNA (presence)']], 'layer');
  const vs = el('select', { id: 'nwview', 'aria-label': 'View' }); NW_VIEWS.forEach(([k, l]) => vs.append(new Option(l, k)));
  vs.onchange = () => { S.view = vs.value; S.sel = null; draw(); }; c1.append(el('label', { for: 'nwview' }, 'Show'), vs);
  const c2 = el('div', { class: 'controls' }); sec.append(c2);
  c2.append(el('span', { style: 'color:var(--ink-2)' }, 'Methods'));
  J.methods.forEach(m => { const id = 'nwm-' + m.key, cb = el('input', { type: 'checkbox', id }); cb.checked = true;
    cb.onchange = () => { S.bits = cb.checked ? S.bits | m.bit : S.bits & ~m.bit; draw(); }; c2.append(cb, el('label', { for: id, style: 'margin-right:6px' }, m.label)); });
  const ms = el('select', { id: 'nwmin', 'aria-label': 'Minimum number of methods' }); [1, 2, 3, 4].forEach(k => ms.append(new Option(`at least ${k}`, k))); ms.value = 2;
  ms.onchange = () => { S.min = +ms.value; draw(); }; c2.append(el('label', { for: 'nwmin' }, 'agreeing'), ms);
  const dc = el('input', { type: 'checkbox', id: 'nwdes' }); dc.checked = true; dc.onchange = () => { S.design = dc.checked; draw(); };
  c2.append(dc, el('label', { for: 'nwdes' }, 'one of them FlashWeave or residual Spearman'));
  const lg = el('div', { class: 'legend' }); sec.append(lg);
  lg.innerHTML = Object.values(NW_DOM).map(([c, l]) => `<span><i style="background:var(${c});height:10px;width:10px;border-radius:50%"></i>${l}</span>`).join('') +
    `<span><i style="background:none;border:1.5px dashed var(--ink-2);height:9px;width:9px;border-radius:50%"></i>doubtful MAG</span>` +
    `<span><i style="background:none;border:3px double var(--ink-2);height:11px;width:11px;border-radius:50%"></i>several genomes together</span>` +
    `<span><i style="background:var(--ink-2)"></i>positive link</span><span><i style="background:var(--s7)"></i>negative link</span>`;
  const box = fig(sec, '', 'Point size follows median abundance. Positions: t-SNE of each genome’s abundance across the libraries (fourth root of its share of each library), one map per layer. Line width grows with the number of agreeing methods. Hover a point for its lineage and partners; click it to show only its links, click the background to clear. The method boxes and the “agreeing” menu change which links are drawn; single-method links of the residual Spearman are not included.');
  const info = el('p', { class: 'cap', style: 'margin-top:4px' }); box.parentNode.insertBefore(info, box.nextSibling);

  function draw() {
    box.innerHTML = '';
    const L = J.layers[S.layer], N = L.nodes;
    const keep = e => { const m = e[3] & S.bits; if (nwPop(m) < S.min) return false; if (e[4] & S.bits) return false; if (S.design && !(m & designBits)) return false; return true; };
    const inView = e => { const a = N[e[0]], b = N[e[1]];
      if (S.view === 'grazer') { const g = n => n.role === 'grazer', t = n => n.domain === 'prokaryote' || n.role === 'alga' || n.domain === 'organelle'; return (g(a) && t(b)) || (g(b) && t(a)); }
      if (S.view === 'virus') return (a.domain === 'virus') !== (b.domain === 'virus');
      if (S.view === 'plastid') return (a.domain === 'organelle' && b.domain === 'eukaryote') || (b.domain === 'organelle' && a.domain === 'eukaryote');
      return true; };
    const E = L.edges.filter(e => keep(e) && inView(e));
    const nbr = new Map(N.map((_, i) => [i, []])); E.forEach(e => { nbr.get(e[0]).push(e); nbr.get(e[1]).push(e); });
    const drv = new Map(); L.drivers.forEach(([i, d, w]) => { (drv.get(i) || drv.set(i, []).get(i)).push([d, w]); });
    const W = width(box), H = Math.max(420, Math.min(640, W * 0.62)), pad = 18;
    const ex = d3.extent(N, n => n.x), ey = d3.extent(N, n => n.y);
    const x = d3.scaleLinear().domain(ex).range([pad, W - pad - 60]), y = d3.scaleLinear().domain(ey).range([H - pad, pad]);
    const rs = d3.scaleSqrt().domain([0, d3.max(N, n => n.ab || 0) || 1]).range([2.2, 13]);
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', H).attr('role', 'img').attr('aria-label', `${S.layer.toUpperCase()} genome network`);
    svg.append('rect').attr('width', W).attr('height', H).attr('fill', 'transparent').on('click', () => { S.sel = null; draw(); });
    const sel = S.sel, selSet = sel == null ? null : new Set([sel, ...nbr.get(sel).map(e => e[0] === sel ? e[1] : e[0])]);
    const Es = sel == null ? E : E.filter(e => e[0] === sel || e[1] === sel);
    const viewNodes = S.view === 'all' ? null : new Set(E.flatMap(e => [e[0], e[1]]));
    svg.append('g').selectAll('line').data(Es).join('line').attr('x1', e => x(N[e[0]].x)).attr('y1', e => y(N[e[0]].y)).attr('x2', e => x(N[e[1]].x)).attr('y2', e => y(N[e[1]].y))
      .attr('stroke', e => css(e[2] > 0 ? '--ink-2' : '--s7')).attr('stroke-width', e => 0.4 + 0.55 * (nwPop(e[3] & S.bits) - 1)).attr('stroke-opacity', sel == null ? (Es.length > 2500 ? 0.12 : Es.length > 800 ? 0.22 : 0.4) : 0.85);
    const order = N.map((n, i) => i).sort((a, b) => (N[b].ab || 0) - (N[a].ab || 0));
    svg.append('g').selectAll('circle').data(order).join('circle')
      .attr('cx', i => x(N[i].x)).attr('cy', i => y(N[i].y)).attr('r', i => rs(N[i].ab || 0))
      .attr('fill', i => N[i].status.includes('doubtful') ? css('--panel') : css(NW_DOM[N[i].domain][0]))
      .attr('stroke', i => N[i].status.includes('doubtful') ? css(NW_DOM[N[i].domain][0]) : N[i].status.includes('multi-genome') ? css('--ink') : css('--panel'))
      .attr('stroke-width', i => N[i].status.includes('multi-genome') ? 3 : N[i].status.includes('doubtful') ? 1.6 : 0.6)
      .attr('stroke-dasharray', i => N[i].status.includes('doubtful') ? '2 1.5' : null)
      .attr('opacity', i => selSet ? (selSet.has(i) ? 1 : 0.12) : viewNodes ? (viewNodes.has(i) ? 1 : 0.15) : (N[i].domain === 'prokaryote' ? 0.8 : 1))
      .style('cursor', 'pointer')
      .on('mousemove', (ev, i) => showTip(ev, nwTip(N, i, nbr.get(i), drv.get(i))))
      .on('mouseleave', hideTip)
      .on('click', (ev, i) => { ev.stopPropagation(); S.sel = S.sel === i ? null : i; draw(); });
    // labels for the non-bacterial genomes that have links in view (bacteria are too many to label)
    const lab = order.filter(i => N[i].domain !== 'prokaryote' && nbr.get(i).length && (!selSet || selSet.has(i)) && (!viewNodes || viewNodes.has(i)));
    svg.append('g').selectAll('text').data(lab.slice(0, sel == null ? 40 : 80)).join('text').attr('x', i => x(N[i].x) + rs(N[i].ab || 0) + 2).attr('y', i => y(N[i].y) + 3.5)
      .style('font-size', '10px').style('pointer-events', 'none').text(i => N[i].domain === 'virus' ? N[i].id.replace('contig_', 'c') : N[i].id);
    const npos = E.filter(e => e[2] > 0).length;
    info.textContent = sel == null ? `${fmt(E.length)} links shown (${fmt(npos)} positive, ${fmt(E.length - npos)} negative) among ${N.length} genomes.`
      : `${N[sel].label}: ${nbr.get(sel).length} link${nbr.get(sel).length === 1 ? '' : 's'} shown.`;
  }
  draw();
}

function nwTip(N, i, edges, drv) {
  const n = N[i], lin = D.lineage && D.lineage[n.id] ? D.lineage[n.id] : '';
  const st = [n.ee_status && n.ee_status !== 'network' ? n.ee_status : '', ...n.status.filter(s => s !== 'doubtful' || !n.ee_status), n.group].filter(Boolean).join(' · ');
  const P = edges.map(e => { const j = e[0] === i ? e[1] : e[0]; return { j, s: e[2], k: nwPop(e[3]) }; })
    .sort((a, b) => b.k - a.k || (N[b.j].domain !== 'prokaryote') - (N[a.j].domain !== 'prokaryote')).slice(0, 7);
  const DL = { winter_2025_26: 'winter 2025-26', diesel_pool: 'diesel pool', control_pool: 'control pool', seawater: 'seawater', sea_ice: 'sea ice', top_ice: 'top ice', bottom_ice: 'bottom ice', day: 'day', under_ice_par: 'under-ice light' };
  return `<b>${n.label}</b>${lin ? `<br><span style="opacity:.8">${lin}</span>` : ''}${st ? `<br>${st}` : ''}<br>median ${n.ab == null ? '–' : fmt1(n.ab)} per million · detected in ${d3.format('.0%')(n.prev)} of libraries` +
    (drv && drv.length ? `<br>FlashWeave drivers: ${drv.map(([d, w]) => `${DL[d] || d} (${w > 0 ? '+' : '−'})`).join(', ')}` : '') +
    (P.length ? `<br><b>Partners</b> (methods agreeing):<br>${P.map(p => `${p.s > 0 ? '+' : '−'} ${N[p.j].label} (${p.k})`).join('<br>')}${edges.length > P.length ? `<br>… ${edges.length - P.length} more` : ''}` : '<br>No links in this view.');
}

function nwPartners(sec, J) {
  sec.append(el('h3', {}, 'Top partners of the eukaryotes and the giant viruses'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const s = el('select', { id: 'nwpart', 'aria-label': 'Genomes' }); [['eukaryote', 'Eukaryotic MAGs'], ['virus', 'Giant viruses and the virophage']].forEach(([k, l]) => s.append(new Option(l, k)));
  c.append(el('label', { for: 'nwpart' }, 'Genomes'), s);
  const box = fig(sec, '', 'Links found by at least two methods (with one of them FlashWeave or the residual Spearman), in either layer; up to eight partners each, those found in both layers first. Methods: F = FlashWeave, S = SPRING, E = SPIEC-EASI, R = residual Spearman. Click a row to show that genome in the network above.');
  const ab = { flashweave: 'F', spring: 'S', spieceasi: 'E', spearman: 'R' };
  const ml = l => l ? l.map(m => ab[m]).join('') : '–';
  function draw() {
    const P = J.partners.filter(p => p.domain === s.value);
    box.innerHTML = `<div class="tablewrap" style="max-height:640px"><table><thead><tr><th>Genome</th><th>Partner</th><th class="num">Sign</th><th>DNA</th><th>RNA</th></tr></thead><tbody>${
      P.map(p => p.partners.map((q, k) => `<tr data-id="${p.id}"><td style="white-space:normal;min-width:200px">${k === 0 ? `<b>${p.label}</b>${p.role ? `<div style="font-size:11.5px;color:var(--ink-3)">${p.role}; ${p.n} partners</div>` : `<div style="font-size:11.5px;color:var(--ink-3)">${p.n} partners</div>`}` : ''}</td><td style="white-space:normal;min-width:240px"><span style="color:var(${NW_DOM[q.domain] ? NW_DOM[q.domain][0] : '--ink-2'})">●</span> ${q.label}</td><td class="num">${q.sign > 0 ? '+' : '−'}</td><td class="mono">${ml(q.layers.dna)}</td><td class="mono">${ml(q.layers.rna)}</td></tr>`).join('')).join('')}</tbody></table></div>`;
    box.querySelectorAll('tr[data-id]').forEach(tr => { tr.onclick = () => nwFocus(J, tr.dataset.id); });
  }
  s.onchange = draw; draw();
}

function nwFocus(J, id) {
  // select the genome in the network (in the layer that has it, RNA first) and scroll there
  const sec = document.getElementById('v-network'); const svgs = sec.querySelectorAll('.fig svg circle'); if (!svgs.length) return;
  const L = J.layers.rna.nodes.some(n => n.id === id) ? 'rna' : 'dna';
  const btn = [...sec.querySelectorAll('button[data-g="layer"]')].find(b => b.textContent.startsWith(L.toUpperCase())); if (btn) btn.click();
  const i = J.layers[L].nodes.findIndex(n => n.id === id);
  const c = [...sec.querySelectorAll('.fig svg circle')].find(e => d3.select(e).datum() === i);
  if (c) { const f = c.closest('.fig'); c.dispatchEvent(new MouseEvent('click', { bubbles: true })); if (f) f.scrollIntoView({ behavior: 'smooth', block: 'center' }); }
}

function nwKnown(sec, J) {
  sec.append(el('h3', {}, 'Check against expected pairs'));
  const pairs = [...new Map(J.known.map(k => [k.a + '|' + k.b, k])).values()];
  const name = id => { for (const l of ['rna', 'dna']) { const n = J.layers[l].nodes.find(x => x.id === id); if (n) return n.label; } return id; };
  const ab = { flashweave: 'FlashWeave', spring: 'SPRING', spieceasi: 'SPIEC-EASI', spearman: 'residual Spearman' };
  const cell = k => !k ? '–' : !k.in_layer ? 'not in this layer' : k.methods.length ? `${k.consensus ? '<b>found</b>' : 'one method'}: ${k.methods.map(m => ab[m] + (k.signs[m] < 0 ? ' (−)' : '')).join(', ')}` : 'not found';
  const box = el('div', { class: 'note' }); sec.append(box);
  box.innerHTML = `<p style="margin:0 0 6px">The Viruses and Eukaryotes pages matched giant viruses, the virophage and a plastid to likely eukaryotic hosts by simple correlation. These pairs were not used to tune the methods.</p>
  <div class="tablewrap" style="max-height:none;border:0"><table><thead><tr><th>Pair</th><th>DNA</th><th>RNA</th></tr></thead><tbody>${
    pairs.map(p => { const d = J.known.find(k => k.a === p.a && k.b === p.b && k.layer === 'dna'), r = J.known.find(k => k.a === p.a && k.b === p.b && k.layer === 'rna');
      return `<tr><td style="white-space:normal;min-width:220px">${name(p.a)} – ${name(p.b)}<div style="font-size:11.5px;color:var(--ink-3)">${p.pair}</div></td><td style="white-space:normal">${cell(d)}</td><td style="white-space:normal">${cell(r)}</td></tr>`; }).join('')}</tbody></table></div>`;
}

function nwMethods(sec, J) {
  sec.append(el('h3', {}, 'How the methods compare'));
  sec.append(el('p', { class: 'kicker' }, 'FlashWeave tests each pair while holding other genomes and the sample drivers (winter, pool, seawater or ice, ice layer, day, modelled under-ice light) fixed, so a link that only reflects a shared driver drops out. SPRING and SPIEC-EASI estimate direct links among the genomes from their relative abundances, but know nothing about the sampling design. The residual Spearman baseline removes the design first and then correlates what is left, pair by pair; it flags many pairs because it does not separate direct from indirect links. The consensus therefore needs two methods, one of them FlashWeave or the baseline.'));
  const box = fig(sec, 'Links per method and overlap', 'Overlap = shared links ÷ links found by either (Jaccard index); sign agreement = share of the shared links with the same sign.');
  const lab = Object.fromEntries(J.methods.map(m => [m.key, m.label])); lab.consensus = 'Consensus';
  box.innerHTML = ['dna', 'rna'].map(l => { const L = J.layers[l];
    return `<div class="title" style="margin-top:6px">${l === 'dna' ? 'DNA' : 'RNA'} network: ${L.nodes.length} genomes (${Object.entries(L.by_domain).map(([d, n]) => `${n} ${NW_DOM[d][1].toLowerCase()}`).join(', ')}), ${L.n_libs} libraries</div>
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
