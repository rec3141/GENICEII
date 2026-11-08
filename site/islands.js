/* Genomic islands section: IslandPath-DIMOB islands in the marine coassembly, their mobility markers, functions, oil genes and expression.
   Data: islands.json (islands/scripts/islands_site.py). Uses helpers from app.js (el, fig, width, css, showTip, hideTip, fmt). */
const GCLS = { 'hydrocarbon degradation': '--exp', 'integrase': '--s7', 'transposase/IS': '--s7', 'phage': '--s5', 'conjugation/T4SS': '--s4',
  'restriction-modification': '--s8', 'CRISPR/defence': '--s8', 'metal resistance': '--s6', 'antibiotic resistance': '--s6',
  'EPS/cell surface': '--s3', 'secretion/adhesion': '--s3', 'transport': '--ctl' };
const GLEG = [['hydrocarbon degradation', 'oil degradation'], ['integrase', 'integrase / transposase'], ['phage', 'phage'], ['conjugation/T4SS', 'conjugation'],
  ['CRISPR/defence', 'defence / R-M / toxin-antitoxin'], ['metal resistance', 'metal / antibiotic resistance'], ['EPS/cell surface', 'cell surface / secretion'], ['transport', 'transport']];
const gcol = c => css(GCLS[c] || (c === 'hypothetical' ? '--grid' : '--other'));

function islandsSec(sec) {
  const J = D.islands;
  if (!J) { sec.innerHTML = '<h2>Genomic islands</h2><p class="kicker">islands.json is missing.</p>'; return; }
  const st = J.stats;
  sec.innerHTML = `<h2>Genomic islands</h2>
  <p class="kicker">IslandPath-DIMOB (run in mag_analysis on the whole coassembly) marks stretches of a contig whose codon usage departs from the rest and that carry mobility genes. Each island here is annotated with Bakta products, corrected KOs, CANT-HYD and the hydrocarbon gene set, plus other transfer evidence: tRNA at an edge, GC shift, geNomad plasmid or provirus overlap, and integrons. Expression is CPM of protein-coding reads in diesel-pool vs control seawater after dosing, using the corrected sample labels.</p>`;
  const sr = el('div', { class: 'stat-row' }); sec.append(sr);
  [[fmt(st.islands), 'islands in the coassembly'], [fmt(st.in_mags), `in ${st.mags} of the 911 MAGs`], [fmt(st.core_like), 'flagged as core-gene clusters (likely false islands)'],
   [`${st.hc_islands}`, `islands carry oil-degradation genes (${st.hc_genes} of ${fmt(st.hc_total_confident)} confident genes)`]]
    .forEach(([v, l]) => sr.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  sec.append(el('p', { class: 'note warn' }, '<b>Two cautions.</b> IslandPath also flags clusters of highly expressed core genes (ribosomal proteins, RNA polymerase, ATP synthase), whose codon usage is atypical. Islands with two or more such genes are marked “core-like” and left out of the induced lists. And because a blooming host raises every one of its genes, induction is shown both raw and as the island’s share of its host MAG’s transcripts.'));
  enrichChart(sec, J);
  sec.append(el('h3', {}, 'Oil genes on islands: what we found'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li><b>Oil genes avoid islands.</b> Only ${st.hc_genes} of ${fmt(st.hc_total_confident)} confident hydrocarbon genes (0.5%) sit in ${st.hc_islands} islands, 0.4× their share elsewhere on the same contigs. Islands are instead full of integrases (9.6×), transposases (7.3×), restriction-modification (4.2×), phage (4.2×) and toxin-antitoxin/defence genes (3×). The main mobile oil element in the pools is the Sphingorhabdus CYP153/alkB plasmid on the Gene transfer page, not an island.</li>
    <li><b>The most complete catabolic island is silent.</b> Colwellia semibin_042 carries a full multicomponent phenol hydroxylase (dmpK–P with its MopR regulator) beside an integrase and transposase (GI_348). It is barely transcribed in either pool.</li>
    <li><b>Oil islands that do respond mostly ride their host’s bloom.</b> A PAH dihydrodiol dehydrogenase (nahB) next to an IS66 transposase in Cycloclasticus metabat_5432 (GI_2286) is 16× higher in the 2025-26 diesel pool, and a 3-hydroxyphenylpropionate (mhp) meta-cleavage island in Parasphingorhabdus (GI_1266) is 14× higher in 2024-25, but both keep about the same share of their host’s transcripts.</li>
    <li><b>Two islands switch on within their hosts.</b> A CYP153 alkane hydroxylase with its ferredoxin reductase and an IS21 transposase (GI_2179, unclassified comebin_1103_sub) rises from 0 to about 2.4 CPM in 2025-26 and runs 24× above its host’s per-gene average. An Oleispira island with an IS1182 transposase and a BON-domain membrane-protection gene (GI_734) reaches 117 CPM in the 2025-26 diesel pool, 39× its host’s average.</li>
    <li><b>Cycloclasticus islands carry motility and osmoprotection, not oil genes.</b> Its flagellin and type IV pilus islands are among the most expressed islands in the diesel pool (560 and 470 CPM in 2025-26), and metabat_9342 carries an ectoine synthesis island (GI_935).</li>`));
  browser(sec, J);
}

/* ---------- what islands carry: enrichment vs the rest of the same contigs ---------- */
function enrichChart(sec, J) {
  sec.append(el('h3', {}, 'What islands carry'));
  const E = J.enrichment.slice().sort((a, b) => b.enrich - a.enrich);
  const box = fig(sec, 'Gene classes inside islands vs elsewhere on the same contigs', 'Ratio of each class’s share of protein-coding genes inside islands to its share outside, on log₂ scale. Classes come from Bakta products; “oil degradation” is the confident hydrocarbon gene set plus CANT-HYD hits.');
  const W = width(box), lw = Math.min(220, W * 0.38), rh = 22, H = E.length * rh + 34;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const x = d3.scaleLinear().domain([-2, 3.5]).range([lw + 10, W - 130]);
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${E.length * rh + 6})`).call(d3.axisBottom(x).ticks(6).tickFormat(v => `${d3.format('~g')(2 ** v)}×`));
  svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 0).attr('y2', E.length * rh + 6).attr('stroke', css('--ink-3'));
  E.forEach((e, i) => {
    const y = i * rh + 12, v = Math.log2(e.enrich);
    svg.append('text').attr('x', lw).attr('y', y + 4).attr('text-anchor', 'end').text(e.cls === 'hydrocarbon degradation' ? 'oil degradation' : e.cls);
    svg.append('rect').attr('x', Math.min(x(0), x(v))).attr('y', y - 6).attr('width', Math.abs(x(v) - x(0))).attr('height', 12).attr('rx', 2).attr('fill', gcol(e.cls))
      .attr('stroke', e.cls === 'hypothetical' || e.cls === 'other enzyme/function' ? css('--ink-3') : 'none')
      .on('mousemove', ev => showTip(ev, `<b>${e.cls}</b><br>${fmt(e.inside)} genes in islands, ${fmt(e.outside)} elsewhere on the same contigs<br>${d3.format('.2f')(e.enrich)}× enriched`)).on('mouseleave', hideTip);
    svg.append('text').attr('x', W - 124).attr('y', y + 4).attr('class', 'mono').text(`${d3.format('.1f')(e.enrich)}× · ${fmt(e.inside)}`);
  });
}

/* ---------- island browser: filters, table, gene map ---------- */
function browser(sec, J) {
  sec.append(el('h3', {}, 'Island browser'));
  sec.append(el('p', { class: 'kicker' }, 'Pick a view, then a row to draw the island. The arrows are genes, coloured by class; the lollipops under them give each gene’s median CPM in diesel-pool (orange) and control (blue) seawater for the chosen winter. Hover a gene for its product, KO and expression.'));
  const VIEWS = [['oil', 'Oil-degradation genes', x => x.hc > 0],
    ['ind25', 'Switched on in diesel 2025-26', x => !x.core && x.share[1] != null && x.share[1] >= 2 && x.oil[1] >= 10],
    ['ind24', 'Switched on in diesel 2024-25', x => !x.core && x.share[0] != null && x.share[0] >= 2 && x.oil[0] >= 10],
    ['phage', 'Phage / prophage', x => (x.cls.phage || 0) >= 5 || x.provirus],
    ['metal', 'Metal / antibiotic resistance', x => (x.cls['metal resistance'] || 0) + (x.cls['antibiotic resistance'] || 0) >= 2],
    ['def', 'Defence and R-M', x => (x.cls['CRISPR/defence'] || 0) + (x.cls['restriction-modification'] || 0) >= 3],
    ['all', 'All islands in MAGs', x => x.bin !== 'unbinned']];
  let view = 'oil', season = 1, sel = null, sortK = 'hc', dir = -1;
  const c = el('div', { class: 'controls' }); sec.append(c);
  const vs = el('select', { id: 'isv', 'aria-label': 'View' }); VIEWS.forEach(([k, l]) => vs.append(new Option(l, k))); c.append(el('label', { for: 'isv' }, 'View'), vs);
  const q = el('input', { type: 'search', id: 'isq', placeholder: 'Filter by taxon, MAG or gene product', 'aria-label': 'Filter islands' }); c.append(q);
  const sb = el('div', { style: 'display:flex;gap:6px' }); c.append(sb);
  ['2024-25', '2025-26'].forEach((s, i) => { const b = el('button', { class: 'chip', 'aria-pressed': i === season }, s); b.onclick = () => { season = i; sb.querySelectorAll('button').forEach((x, j) => x.setAttribute('aria-pressed', j === i)); table(); if (sel) draw(sel); }; sb.append(b); });
  const n = el('span', { class: 'mono' }); c.append(n);
  const g2 = el('div'); sec.append(g2);
  const tw = el('div', { class: 'tablewrap', style: 'max-height:380px' }); g2.append(tw);
  const mapBox = el('div'); g2.append(mapBox);
  const COLS = [['id', 'Island'], ['taxon', 'Host'], ['len', 'kb'], ['n', 'Genes'], ['hc', 'Oil genes'], ['mob', 'Mobility'], ['oilc', 'CPM diesel'], ['lfc', 'log₂ diesel/ctl'], ['share', 'log₂ share in host'], ['act', 'Activity']];
  const val = (x, k) => ({ oilc: x.oil[season], lfc: x.lfc[season], share: x.share[season], act: x.act[season] })[k] ?? x[k];
  function table() {
    const f = VIEWS.find(v => v[0] === view)[2], qq = q.value.toLowerCase();
    let R = J.islands.filter(f).filter(x => !qq || `${x.taxon} ${x.bin} ${x.family} ${x.hc_cat}`.toLowerCase().includes(qq) ||
      (J.genes[x.id] || []).some(g => (g.p || '').toLowerCase().includes(qq)));
    R.sort((a, b) => dir * ((val(a, sortK) ?? -1e9) - (val(b, sortK) ?? -1e9) || String(val(a, sortK)).localeCompare(String(val(b, sortK)))));
    n.textContent = `${R.length} islands`;
    tw.innerHTML = `<table><thead><tr>${COLS.map(([k, l]) => `<th class="${['id', 'taxon'].includes(k) ? '' : 'num'}" data-k="${k}" aria-sort="${sortK === k ? (dir > 0 ? 'ascending' : 'descending') : 'none'}">${l}</th>`).join('')}</tr></thead><tbody>${
      R.slice(0, 400).map(x => `<tr data-id="${x.id}" class="${x.id === sel ? 'sel' : ''}"><td class="mono">${x.id}${x.core ? ' <span class="pill na">core-like</span>' : ''}</td><td>${x.taxon} <span class="mono" style="color:var(--ink-3)">${x.bin}</span></td><td class="num">${d3.format('.1f')(x.len / 1000)}</td><td class="num">${x.n}</td><td class="num">${x.hc || ''}</td><td class="num">${x.mob}</td>
        <td class="num">${x.oil[season] == null ? '–' : d3.format(',.1f')(x.oil[season])}</td>${['lfc', 'share'].map(k => { const v = (x.oil[season] ?? 0) < 0.5 && (x.ctl[season] ?? 0) < 0.5 ? null : x[k][season]; return `<td class="num">${v == null ? '–' : `<span class="pill ${v > 0 ? 'up' : 'down'}">${v > 0 ? '+' : ''}${d3.format('.1f')(v)}</span>`}</td>`; }).join('')}<td class="num">${x.act[season] == null ? '–' : d3.format('.1f')(x.act[season])}</td></tr>`).join('')}</tbody></table>`;
    tw.querySelectorAll('th').forEach(th => th.onclick = () => { const k = th.dataset.k; dir = sortK === k ? -dir : -1; sortK = k; table(); });
    tw.querySelectorAll('tbody tr').forEach(tr => tr.onclick = () => { sel = tr.dataset.id; table(); draw(sel); });
    if (!sel && R.length) { sel = R[0].id; draw(sel); table(); }
  }
  function draw(id) {
    mapBox.innerHTML = '';
    const x = J.islands.find(i => i.id === id), G = J.genes[id];
    const ev = [x.tRNA ? 'tRNA at edge' : '', Math.abs(x.dGC) >= 3 ? `GC ${x.dGC > 0 ? '+' : ''}${x.dGC} points vs contig` : `GC ${x.dGC > 0 ? '+' : ''}${x.dGC}`, x.plasmid ? `geNomad plasmid score ${x.plasmid}` : '', x.provirus ? 'provirus overlap' : '', x.integron ? 'integron overlap' : ''].filter(Boolean).join(' · ');
    const box = fig(mapBox, `${x.id} · ${x.taxon} (${x.bin}) · ${x.contig}:${fmt(x.start)}–${fmt(x.end)}`, `${d3.format('.1f')(x.len / 1000)} kb, ${x.n} genes. ${ev}.${x.hc_cat ? ' Oil genes: ' + x.hc_cat + '.' : ''} Mobility score ${x.mob} counts integrase, transposase, tRNA edge, ≥3 phage genes, ≥2 conjugation genes, integron, provirus and |ΔGC| ≥ 3.`);
    if (!G) { box.append(el('p', { class: 'cap' }, 'Gene map not included for this island (only oil-gene, induced, mobile, phage, resistance and defence islands carry gene maps).')); return; }
    const lg = el('div', { class: 'legend' }, GLEG.map(([k, l]) => `<span><i style="background:var(${GCLS[k]});height:10px"></i>${l}</span>`).join('') + '<span><i style="background:var(--other);height:10px"></i>other</span>');
    box.parentNode.insertBefore(lg, box);
    const W = Math.max(width(box), 700), m = { l: 50, r: 16 }, ty = 30, th = 18, ey = 70, eh = 110, H = ey + eh + 34;
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
    const lo = d3.min(G, g => g.s), hi = d3.max(G, g => g.e);
    const xs = d3.scaleLinear().domain([lo, hi]).range([m.l, W - m.r]);
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${ty - 10})`).call(d3.axisTop(xs).ticks(6).tickFormat(v => `${d3.format('.1f')(v / 1000)} kb`));
    svg.append('line').attr('x1', xs(lo)).attr('x2', xs(hi)).attr('y1', ty + th / 2).attr('y2', ty + th / 2).attr('stroke', css('--ink-3'));
    const all = G.flatMap(g => [g.oil[season], g.ctl[season]]).filter(v => v != null);
    const ys = d3.scaleSymlog().constant(1).domain([0, Math.max(10, d3.max(all) || 10)]).range([ey + eh, ey]);
    svg.append('g').attr('class', 'axis').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(ys).ticks(4).tickFormat(tickf));
    svg.append('text').attr('transform', `translate(10,${ey}) rotate(90)`).text('CPM');
    const tipG = g => `<b>${g.g ? g.g + ' · ' : ''}${g.p}</b><br>${g.t}${g.hc ? ' · ' + g.hc : ''}${g.ko ? ' · ' + g.ko : ''}<br>${g.c}<br>${['2024-25', '2025-26'].map((s, i) => `${s}: diesel ${g.oil[i] ?? '–'}, control ${g.ctl[i] ?? '–'} CPM`).join('<br>')}`;
    G.forEach(g => {
      const a = xs(g.s), b = xs(g.e), w = Math.max(2, b - a), hd = Math.min(8, w * 0.4), y0 = ty, y1 = ty + th;
      const pts = g.st === '+' ? [[a, y0], [b - hd, y0], [b, (y0 + y1) / 2], [b - hd, y1], [a, y1]] : [[b, y0], [a + hd, y0], [a, (y0 + y1) / 2], [a + hd, y1], [b, y1]];
      svg.append('polygon').attr('points', pts.map(p => p.join(',')).join(' ')).attr('fill', g.t === 'cds' ? gcol(g.c) : 'none').attr('stroke', g.t === 'cds' ? css('--panel') : css('--ink-3'))
        .attr('stroke-width', 1).on('mousemove', ev => showTip(ev, tipG(g))).on('mouseleave', hideTip);
      if (g.t !== 'cds' || g.oil[season] == null) return;
      const mx = (a + b) / 2;
      [[g.ctl[season], '--ctl', -3], [g.oil[season], '--exp', 3]].forEach(([v, cc, dx]) => {
        if (v == null) return;
        svg.append('line').attr('x1', mx + dx).attr('x2', mx + dx).attr('y1', ys(0)).attr('y2', ys(v)).attr('stroke', css(cc)).attr('stroke-width', 2);
        svg.append('circle').attr('cx', mx + dx).attr('cy', ys(v)).attr('r', 3).attr('fill', css(cc)).on('mousemove', ev => showTip(ev, tipG(g))).on('mouseleave', hideTip);
      });
    });
    svg.append('text').attr('x', W - m.r).attr('y', H - 4).attr('text-anchor', 'end').text(`median CPM in ${['2024-25', '2025-26'][season]} seawater: diesel (orange) vs control (blue)`);
  }
  vs.onchange = () => { view = vs.value; sel = null; sortK = view === 'oil' ? 'hc' : view.startsWith('ind') ? 'share' : 'mob'; season = view === 'ind24' ? 0 : view === 'ind25' ? 1 : season;
    sb.querySelectorAll('button').forEach((x, j) => x.setAttribute('aria-pressed', j === season)); table(); };
  q.oninput = () => table();
  table();
}

RENDER.islands = islandsSec;
