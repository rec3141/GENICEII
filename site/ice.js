/* Ice vs water section: ice physics, ice-physiology panel, MAG niches (ice/water x top/bottom), rRNA taxa, cold adaptation, population genetics.
   Data: ice.json (built by ice/scripts/ice_site.py). Uses helpers from app.js (el, fig, width, css, showTip, hideTip, $). */
const LAYER = { SW: '--ctl', TOP: '--s4', BOT: '--s3' };
const sgn = v => (v > 0 ? '+' : '') + d3.format('.2f')(v);

function iceSec(sec) {
  const I = D.ice;
  if (!I) { sec.innerHTML = '<h2>Ice vs water</h2><p class="kicker">ice.json is missing.</p>'; return; }
  const top = I.cores.filter(c => c.Sbrine_top != null), bot = I.cores.filter(c => c.Sbrine_bot != null);
  const rng = (a, k, f = d3.format('.0f')) => { const lo = d3.min(a, x => x[k]), hi = d3.max(a, x => x[k]); return lo < 0 ? `${f(lo)} to ${f(hi)}` : `${f(lo)}–${f(hi)}`; };
  sec.innerHTML = `<h2>Ice vs water, top vs bottom</h2>
  <p class="kicker">Both pools froze over within days of filling. 2024-25 ice was cored at the top and bottom on the same days as seawater was sampled (13 top/bottom pairs). 2025-26 has whole-ice cores paired with seawater. Each comparison below is paired: ice minus water from the same pool on the same day, and top minus bottom from the same core. With one pool per treatment, the pairs are repeated measures in time. Read “n of N pairs” as consistency, not as replication.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  [[`${rng(I.cores.filter(c => c.T_top != null), 'T_top')} °C`, 'temperature in the top 5 cm of ice'], [`${rng(I.cores.filter(c => c.T_bot != null), 'T_bot', d3.format('.1f'))} °C`, 'temperature at the ice bottom'],
   [`${rng(top, 'Sbrine_top')} PSU`, 'brine salinity at the top (seawater ≈ 26–32)'], [`${rng(top, 'Vb_top', v => d3.format('.0f')(v * 100))}% vs ${rng(bot, 'Vb_bot', v => d3.format('.0f')(v * 100))}%`, 'brine volume, top vs bottom (≈5% = permeability threshold)']]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  physics(sec, I);
  panelChart(sec, I);
  magNiche(sec, I);
  rrnaTaxa(sec, I);
  divergence(sec, I);
  if (I.cold) coldSec(sec, I.cold); else sec.append(el('p', { class: 'note' }, 'Cold-adaptation genomics: still being computed.'));
  if (I.popgen) popgenSec(sec, I.popgen); else sec.append(el('p', { class: 'note' }, 'Population genetics: still being computed.'));
}

/* ---------- temperature profiles + brine table ---------- */
function physics(sec, I) {
  sec.append(el('h3', {}, 'The ice as a habitat'));
  sec.append(el('p', { class: 'kicker' }, 'Air temperatures reached −22 to −28 °C while the ice bottom stayed near the freezing point of seawater (about −1.8 °C). Brine pockets shrink and concentrate as the ice cools. Brine salinity comes from Assur’s freezing-point relation and brine volume from Frankenstein and Garner (1967), using each core’s top and bottom temperature and bulk salinity. Upper ice was usually below the ~5% brine volume at which sea ice stops being permeable, so top-ice cells sit in isolated, hypersaline pockets. Bottom ice stays connected to the water.'));
  const box = fig(sec, 'Ice temperature profiles', 'Each line is one core (temperature probe every 2.5–10 cm). Darker lines are later in the season. Hover for core details. 2025-26 control day 12 is omitted because its logged profile is inverted.');
  const W = width(box), seasons = ['2024-25', '2025-26'], cols = W < 620 ? 1 : 2, pw = cols === 1 ? W - 10 : (W - 28) / 2, ph = 260, m = { l: 44, r: 10, t: 22, b: 34 };
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', (ph + 12) * (cols === 1 ? 2 : 1));
  const maxDay = d3.max(I.profiles, p => p.day);
  const col = d3.scaleSequential(d3.interpolateRgb(css('--seq1'), css('--seq4'))).domain([0, Math.log1p(maxDay)]);
  seasons.forEach((s, i) => {
    const g = svg.append('g').attr('transform', `translate(${cols === 1 ? 0 : i * (pw + 28)},${cols === 1 ? i * (ph + 12) : 0})`);
    const P = I.profiles.filter(p => p.season === s);
    const x = d3.scaleLinear().domain([-17, -1]).range([m.l, pw - m.r]), y = d3.scaleLinear().domain([0, 90]).range([m.t, ph - m.b]);
    g.append('g').attr('class', 'gridl').attr('transform', `translate(0,${ph - m.b})`).call(d3.axisBottom(x).ticks(5).tickSize(-(ph - m.t - m.b)).tickFormat('')).select('.domain').remove();
    g.append('g').attr('class', 'axis').attr('transform', `translate(0,${ph - m.b})`).call(d3.axisBottom(x).ticks(5));
    g.append('g').attr('class', 'axis').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5));
    g.append('text').attr('class', 'ttl').attr('x', m.l).attr('y', 12).text(`CMO ${s}`);
    g.append('text').attr('x', pw - m.r).attr('y', ph - 4).attr('text-anchor', 'end').text('temperature, °C');
    g.append('text').attr('transform', `translate(10,${m.t}) rotate(90)`).text('depth below ice surface, cm');
    const ln = d3.line().x(d => x(d[1])).y(d => y(d[0]));
    P.sort((a, b) => a.day - b.day).forEach(p => {
      const c = col(Math.log1p(p.day));
      g.append('path').attr('d', ln(p.pts)).attr('fill', 'none').attr('stroke', c).attr('stroke-width', 2).attr('stroke-dasharray', p.pool === 'CTL' ? '5 3' : null);
      g.append('path').attr('d', ln(p.pts)).attr('fill', 'none').attr('stroke', 'transparent').attr('stroke-width', 12)
        .on('mousemove', ev => showTip(ev, `<b>${p.pool === 'CTL' ? 'Control' : 'Diesel'} pool, day ${p.day}</b><br>ice ${p.H ?? '?'} cm thick<br>${p.pts.map(q => `${q[0]} cm: ${q[1]} °C`).join('<br>')}`)).on('mouseleave', hideTip);
    });
  });
  const lg = el('div', { class: 'legend' }, `<span><i style="background:var(--seq2)"></i>early (days 2–20)</span><span><i style="background:var(--seq4)"></i>late (days 48 onward)</span><span><i class="dash" style="color:var(--ink-2)"></i>control pool</span><span><i style="background:var(--ink-2)"></i>diesel pool</span>`);
  box.parentNode.insertBefore(lg, box);
  const rows = I.cores.filter(c => c.Sbrine_top != null || c.T_top != null);
  const tw = el('div', { class: 'tablewrap', style: 'margin-top:10px;max-height:340px' });
  tw.innerHTML = `<table><thead><tr><th>Season</th><th>Pool</th><th class="num">Day</th><th class="num">Ice, cm</th><th class="num">Snow, cm</th><th class="num">Air °C</th><th class="num">Top °C</th><th class="num">Bottom °C</th><th class="num">Bulk S top</th><th class="num">Bulk S bottom</th><th class="num">Brine S top</th><th class="num">Brine S bottom</th><th class="num">Brine vol. top</th><th class="num">Brine vol. bottom</th></tr></thead><tbody>${
    rows.map(c => `<tr><td>${c.season}</td><td>${c.pool === 'CTL' ? 'control' : 'diesel'}</td>${[c.day, c.thickness_cm, c.snow_cm, c.T_air, c.T_top, c.T_bot, c.S_top, c.S_bot, c.Sbrine_top, c.Sbrine_bot].map((v, k) => `<td class="num">${v == null ? '–' : (k >= 8 ? d3.format('.0f')(v) : d3.format('~g')(v))}</td>`).join('')}${[c.Vb_top, c.Vb_bot].map(v => `<td class="num">${v == null ? '–' : d3.format('.1%')(v)}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
  box.parentNode.append(tw);
}

/* ---------- curated physiology panel: per-pair deltas as dots + mean bar ---------- */
function panelChart(sec, I) {
  if (!I.panel) return;
  sec.append(el('h3', {}, 'What ice life switches on'));
  sec.append(el('p', { class: 'kicker' }, 'Twenty curated groups of KEGG orthologs, summed per library (CPM of protein-coding reads, corrected KO assignments). Each dot is one pair; the bar is the mean log₂ ratio. Pick a comparison. Hover a row to see which genomes carry most of those transcripts in ice.'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const opts = [['iw', 'Ice minus water (24 pairs)'], ['tb', 'Top minus bottom ice (13 pairs, 2024-25)']]; let mode = 'iw';
  opts.forEach(([k, l]) => { const b = el('button', { class: 'chip', 'aria-pressed': k === mode }, l); b.onclick = () => { mode = k; c.querySelectorAll('button').forEach(x => x.setAttribute('aria-pressed', x === b)); draw(); }; c.append(b); });
  const box = fig(sec, '', 'Log₂ ratio of category CPM within each pair. Orange means higher in ice (or top ice); blue means higher in water (or bottom ice). Pairs are days within one pool, so “n of N” measures consistency, not replication.');
  function draw() {
    box.innerHTML = '';
    const P = I.panel.slice().sort((a, b) => b[mode] - a[mode]);
    const W = width(box), lw = Math.min(330, W * 0.45), rh = 24, H = P.length * rh + 40;
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
    const all = P.flatMap(p => p[mode + '_pairs']); const ext = Math.min(4, Math.max(1.5, d3.max(all, v => Math.abs(v))));
    const x = d3.scaleLinear().domain([-ext, ext]).range([lw + 10, W - 70]);
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${P.length * rh + 6})`).call(d3.axisBottom(x).ticks(7));
    svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 0).attr('y2', P.length * rh + 6).attr('stroke', css('--ink-3'));
    svg.append('text').attr('x', W - 70).attr('y', P.length * rh + 36).attr('text-anchor', 'end').text(mode === 'iw' ? '← water   log₂ ratio   ice →' : '← bottom   log₂ ratio   top →');
    P.forEach((p, i) => {
      const y = i * rh + 12, mu = p[mode], up = p[mode + '_up'], n = p[mode + '_n'];
      const g = svg.append('g');
      g.append('rect').attr('x', 0).attr('y', y - rh / 2).attr('width', W).attr('height', rh).attr('fill', 'transparent')
        .on('mousemove', ev => showTip(ev, `<b>${p.cat}</b><br>mean ${sgn(mu)} (higher in ${up} of ${n} pairs)<br>median CPM 2024-25: water ${p.cpm.SW}, top ${p.cpm.TOP}, bottom ${p.cpm.BOT}<br>${p.hosts.length ? 'main carriers in ice: ' + p.hosts.slice(0, 3).map(h => `${h.taxon} ${d3.format('.0%')(h.share)}`).join(', ') : ''}`)).on('mouseleave', hideTip);
      g.append('text').attr('x', lw).attr('y', y + 4).attr('text-anchor', 'end').text(p.cat);
      g.append('rect').attr('x', Math.min(x(0), x(mu))).attr('y', y - 5).attr('width', Math.abs(x(mu) - x(0))).attr('height', 10).attr('rx', 2)
        .attr('fill', css(mu > 0 ? '--exp' : '--ctl')).attr('opacity', 0.35).attr('pointer-events', 'none');
      p[mode + '_pairs'].forEach(v => g.append('circle').attr('cx', x(Math.max(-ext, Math.min(ext, v)))).attr('cy', y).attr('r', 3).attr('fill', css(v > 0 ? '--exp' : '--ctl')).attr('opacity', .8).attr('pointer-events', 'none'));
      g.append('text').attr('x', W - 66).attr('y', y + 4).attr('class', 'mono').text(`${up}/${n}`);
    });
  }
  draw();
}

/* ---------- MAG niches: ice-water (x) vs top-bottom (y) ---------- */
function magNiche(sec, I) {
  sec.append(el('h3', {}, 'Who lives where: niches of genomes and their siblings'));
  sec.append(el('p', { class: 'kicker' }, 'Each point is a MAG: the x-axis is its paired transcript-share shift from water to ice, the y-axis from bottom to top ice (CLR of TPM, log₂). Highlight a genus to see whether its species share a niche or split it. Several “Gloeomargaritales” MAGs are 100–280 kb algal chloroplast genomes (SILVA plastid 16S: Attheya septentrionalis, Guillardia, Costaria and others). They are labelled as plastids here.'));
  const genera = d3.rollups(I.mags.filter(m => m.genus), v => v.length, m => m.genus).filter(d => d[1] >= 3).sort((a, b) => b[1] - a[1]).map(d => d[0]);
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'nichesel', 'aria-label': 'Highlight genus' }); ['Cognaticolwellia vs Colwellia', ...genera].forEach(g => sel.append(new Option(g, g)));
  c.append(el('label', { for: 'nichesel' }, 'Highlight'), sel);
  const box = fig(sec, '', 'Point area scales with peak TPM. Filled points have q < 0.05 for ice vs water or top vs bottom (Wilcoxon over pairs, Benjamini-Hochberg).');
  const tw = el('div', { class: 'tablewrap', style: 'margin-top:10px;max-height:300px' }); box.parentNode.append(tw);
  function draw() {
    box.innerHTML = '';
    const hl = sel.value === 'Cognaticolwellia vs Colwellia' ? ['Cognaticolwellia', 'Colwellia'] : [sel.value];
    const W = width(box), H = Math.min(520, W * 0.75), m = { l: 46, r: 14, t: 14, b: 40 };
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
    const x = d3.scaleLinear().domain(d3.extent(I.mags, d => d.iw)).nice().range([m.l, W - m.r]), y = d3.scaleLinear().domain(d3.extent(I.mags.filter(d => d.tb != null), d => d.tb)).nice().range([H - m.b, m.t]);
    const r = d3.scaleSqrt().domain([0, d3.max(I.mags, d => d.tpm)]).range([2, 14]);
    svg.append('g').attr('class', 'gridl').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6).tickSize(-(H - m.t - m.b)).tickFormat('')).select('.domain').remove();
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6));
    svg.append('g').attr('class', 'axis').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6));
    svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', m.t).attr('y2', H - m.b).attr('stroke', css('--ink-3'));
    svg.append('line').attr('y1', y(0)).attr('y2', y(0)).attr('x1', m.l).attr('x2', W - m.r).attr('stroke', css('--ink-3'));
    svg.append('text').attr('x', W - m.r).attr('y', H - 6).attr('text-anchor', 'end').text('← water        ice →');
    svg.append('text').attr('transform', `translate(12,${m.t}) rotate(90)`).text('← bottom ice     top ice →');
    const colors = ['--exp', '--ctl', '--s3', '--s7'];
    const pts = I.mags.filter(d => d.tb != null).sort((a, b) => hl.includes(a.genus) - hl.includes(b.genus));
    pts.forEach(d => {
      const k = hl.indexOf(d.genus), on = k >= 0, sig = (d.iw_q < 0.05) || (d.tb_q < 0.05);
      svg.append('circle').attr('cx', x(d.iw)).attr('cy', y(d.tb)).attr('r', r(d.tpm)).attr('fill', on && sig ? css(colors[k]) : 'none')
        .attr('stroke', on ? css(colors[k]) : css('--ink-3')).attr('stroke-width', on ? 2 : 1).attr('opacity', on ? 1 : .35)
        .on('mousemove', ev => showTip(ev, `<b>${d.bin}</b> · ${d.taxon}<br>ice − water ${sgn(d.iw)} (${d.iw_up}/${d.iw_n} pairs, q ${d3.format('.2g')(d.iw_q)})<br>top − bottom ${sgn(d.tb)} (${d.tb_up}/${d.tb_n}, q ${d3.format('.2g')(d.tb_q)})<br>peak ${fmt(d.tpm)} TPM`)).on('mouseleave', hideTip);
      if (on && d.tpm > 20000) svg.append('text').attr('x', x(d.iw) + r(d.tpm) + 3).attr('y', y(d.tb) + 4).text(d.bin);
    });
    const lg = el('div', { class: 'legend' }, hl.map((g, k) => `<span><i style="background:var(${colors[k]});height:10px;width:10px;border-radius:50%"></i>${g}</span>`).join('') + '<span><i style="border:1px solid var(--ink-3);height:10px;width:10px;border-radius:50%;background:none"></i>other MAGs</span>');
    box.prepend(lg);
    const rows = I.mags.filter(d => hl.includes(d.genus)).sort((a, b) => b.iw - a.iw);
    tw.innerHTML = `<table><thead><tr><th>MAG</th><th>Taxon</th><th class="num">ice − water</th><th class="num">pairs up</th><th class="num">top − bottom</th><th class="num">pairs up</th><th class="num">peak TPM</th></tr></thead><tbody>${
      rows.map(d => `<tr><td class="mono">${d.bin}</td><td>${d.taxon}</td><td class="num"><span class="pill ${d.iw > 0 ? 'up' : 'down'}">${sgn(d.iw)}</span></td><td class="num">${d.iw_up}/${d.iw_n}</td><td class="num">${d.tb == null ? '–' : `<span class="pill ${d.tb > 0 ? 'up' : 'down'}">${sgn(d.tb)}</span>`}</td><td class="num">${d.tb_up}/${d.tb_n}</td><td class="num">${fmt(d.tpm)}</td></tr>`).join('')}</tbody></table>`;
  }
  sel.onchange = draw; draw();
}

/* ---------- rRNA taxa: diverging bars ---------- */
function rrnaTaxa(sec, I) {
  sec.append(el('h3', {}, 'Ice taxa from rRNA reads'));
  sec.append(el('p', { class: 'kicker' }, 'rRNA reads are classified against SILVA independently of the MAGs, so they also cover eukaryotes and unbinned taxa. Shares are within eukaryotes or within prokaryotes (plastids and mitochondria excluded). The bars show the strongest paired shifts.'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  let dom = 'euk', k = 'icewater';
  const mk = (arr, get, set) => arr.forEach(([v, l]) => { const b = el('button', { class: 'chip', 'aria-pressed': v === get() }, l); b.onclick = () => { set(v); b.parentNode.querySelectorAll(`button[data-g="${arr[0][0]}"]`).forEach(x => x.setAttribute('aria-pressed', x === b)); draw(); }; b.dataset.g = arr[0][0]; c.append(b); });
  mk([['euk', 'Eukaryotes'], ['prok', 'Prokaryotes']], () => dom, v => dom = v);
  c.append(el('span', { style: 'width:12px' }));
  mk([['icewater', 'Ice minus water'], ['topbot', 'Top minus bottom']], () => k, v => k = v);
  const box = fig(sec, '', 'Mean paired log₂ ratio of each genus’s share. The label gives the last three SILVA ranks. Spruce pollen, lichen photobionts (Trebouxia) and Dothideales fungi, enriched in ice and especially at its surface, point to wind-blown material from the boreal forest settling on the ice.');
  function draw() {
    box.innerHTML = '';
    const R = I.rrna[`${dom}_${k}`].slice().sort((a, b) => b.delta - a.delta);
    const W = width(box), lw = Math.min(360, W * 0.5), rh = 20, H = R.length * rh + 34;
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
    const ext = d3.max(R, d => Math.abs(d.delta));
    const x = d3.scaleLinear().domain([-ext, ext]).range([lw + 10, W - 50]);
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${R.length * rh + 4})`).call(d3.axisBottom(x).ticks(6));
    svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 0).attr('y2', R.length * rh + 4).attr('stroke', css('--ink-3'));
    R.forEach((d, i) => {
      const y = i * rh + 10;
      svg.append('text').attr('x', lw).attr('y', y + 4).attr('text-anchor', 'end').text(d.short.length > 58 ? d.short.slice(0, 57) + '…' : d.short);
      svg.append('rect').attr('x', Math.min(x(0), x(d.delta))).attr('y', y - 6).attr('width', Math.abs(x(d.delta) - x(0))).attr('height', 12).attr('rx', 2)
        .attr('fill', css(d.delta > 0 ? '--exp' : '--ctl')).attr('opacity', d.q < 0.05 ? 1 : 0.4)
        .on('mousemove', ev => showTip(ev, `<b>${d.lineage}</b><br>${sgn(d.delta)} (higher in ${d.up} of ${d.n} pairs, q ${d3.format('.2g')(d.q)})<br>mean share ${d3.format('.2%')(d.share)}`)).on('mouseleave', hideTip);
      svg.append('text').attr('x', W - 46).attr('y', y + 4).attr('class', 'mono').text(`${d.up}/${d.n}`);
    });
  }
  draw();
}

/* ---------- divergence of ice from paired water vs thickness ---------- */
function divergence(sec, I) {
  const box = fig(sec, 'Ice communities drift away from the water as the ice thickens', 'Aitchison distance (CLR of MAG transcript shares) between each ice library and the seawater sampled the same day in the same pool. For scale, replicate seawater libraries from the same day sit about 9 apart.');
  const W = width(box), H = 240, m = { l: 44, r: 14, t: 14, b: 38 };
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const P = I.divergence.filter(d => d.H != null);
  const x = d3.scaleLinear().domain([0, d3.max(P, d => d.H) + 5]).range([m.l, W - m.r]), y = d3.scaleLinear().domain([0, d3.max(P, d => d.d) + 2]).range([H - m.b, m.t]);
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6));
  svg.append('g').attr('class', 'axis').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5));
  svg.append('line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(9.2)).attr('y2', y(9.2)).attr('stroke', css('--ink-3')).attr('stroke-dasharray', '4 3');
  svg.append('text').attr('x', W - m.r).attr('y', y(9.2) - 4).attr('text-anchor', 'end').text('seawater replicates');
  svg.append('text').attr('x', W - m.r).attr('y', H - 4).attr('text-anchor', 'end').text('ice thickness, cm');
  svg.append('text').attr('transform', `translate(12,${m.t}) rotate(90)`).text('distance to paired water');
  const sh = { TOP: d3.symbolTriangle, BOT: d3.symbolCircle, SI: d3.symbolSquare, HIGH: d3.symbolDiamond, LOW: d3.symbolDiamond, PANCAKE: d3.symbolStar, MID: d3.symbolCross };
  P.forEach(d => svg.append('path').attr('d', d3.symbol(sh[d.layer] || d3.symbolSquare, 70)()).attr('transform', `translate(${x(d.H)},${y(d.d)})`)
    .attr('fill', css(d.season === '2024-25' ? '--s7' : '--s3')).attr('stroke', css('--panel'))
    .on('mousemove', ev => showTip(ev, `<b>${d.name}</b><br>${d.season}, ${d.pool === 'CTL' ? 'control' : 'diesel'} pool, day ${d.day}<br>ice ${d.H} cm · distance ${d.d}`)).on('mouseleave', hideTip));
  box.parentNode.insertBefore(el('div', { class: 'legend' }, `<span><i style="background:var(--s7);height:8px;width:8px;border-radius:50%"></i>2024-25</span><span><i style="background:var(--s3);height:8px;width:8px;border-radius:50%"></i>2025-26</span><span>▲ top ice · ● bottom ice · ■ whole core</span>`), box);
}

/* ---------- cold adaptation (fork output) ---------- */
function coldSec(sec, C) {
  sec.append(el('h3', {}, 'Cold adaptation in genomes and transcripts'));
  coldCharts(sec, C);
}
function popgenSec(sec, G) {
  sec.append(el('h3', {}, 'Population genetics: ice and water populations of the same species'));
  if (G.summary) sec.append(el('p', { class: 'kicker' }, G.summary));
  popgenCharts(sec, G);
}
/* filled in once the fork JSON structures are known */
function coldCharts(sec, C) {
  sec.append(el('p', { class: 'kicker' }, 'Fifteen cold-adaptation gene categories, built from KEGG and Pfam (ice-binding proteins are DUF3494/PF11999), measured as CPM of protein-coding reads. The response is regulation within each taxon more than a change in who is present: cold-shock proteins rise in ice within 84% of 167 MAGs tested, and DEAD-box helicases within 77%. The amino-acid make-up of the expressed proteome barely moves between ice and water (≤0.2 SD).'));
  // A. category medians with IQR, ice-water and top-bottom
  const cats = C.categories.slice().sort((a, b) => b.ice.median_log2 - a.ice.median_log2);
  const box = fig(sec, 'Cold-adaptation categories: paired log₂ ratios', 'Dot = median over pairs; line = interquartile range. Numbers give the fraction of pairs higher in ice (left) or top ice (right), and the fraction of MAGs where the category rises in ice.');
  const W = width(box), lw = Math.min(250, W * 0.36), rh = 24, H = cats.length * rh + 44, cols = W < 620 ? 1 : 2;
  const pw = (W - lw - 70) / 2;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const ext = 2.6;
  [['ice', 'ice − water', 0], ['top', 'top − bottom ice', 1]].forEach(([k, t, i]) => {
    const x0 = lw + 10 + i * (pw + 10), x = d3.scaleLinear().domain([-ext, ext]).range([x0, x0 + pw - 46]);
    svg.append('text').attr('class', 'ttl').attr('x', x(0)).attr('y', 10).attr('text-anchor', 'middle').text(t);
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${cats.length * rh + 22})`).call(d3.axisBottom(x).ticks(5));
    svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 16).attr('y2', cats.length * rh + 22).attr('stroke', css('--ink-3'));
    cats.forEach((c, j) => { const y = j * rh + 30, v = c[k]; const cl = Math.max(-ext, Math.min(ext, v.median_log2));
      svg.append('line').attr('x1', x(Math.max(-ext, v.q25))).attr('x2', x(Math.min(ext, v.q75))).attr('y1', y).attr('y2', y).attr('stroke', css(v.median_log2 > 0 ? '--exp' : '--ctl')).attr('stroke-width', 2).attr('opacity', .5);
      svg.append('circle').attr('cx', x(cl)).attr('cy', y).attr('r', 5).attr('fill', css(v.median_log2 > 0 ? '--exp' : '--ctl')).attr('stroke', css('--panel'))
        .on('mousemove', ev => showTip(ev, `<b>${c.name}</b> · ${t}<br>median ${sgn(v.median_log2)}, IQR ${sgn(v.q25)} to ${sgn(v.q75)}<br>higher in ${d3.format('.0%')(v.frac_up)} of ${v.n} pairs<br>${c.n_genes} genes${k === 'ice' && c.within_mag ? `<br>rises in ice within ${d3.format('.0%')(c.within_mag.frac_up)} of ${c.within_mag.n_mags} MAGs` : ''}`)).on('mouseleave', hideTip);
      svg.append('text').attr('x', x0 + pw - 6).attr('y', y + 4).attr('class', 'mono').attr('text-anchor', 'end').text(d3.format('.0%')(v.frac_up)); });
  });
  cats.forEach((c, j) => { const y = j * rh + 30;
    svg.append('text').attr('x', lw).attr('y', y + 4).attr('text-anchor', 'end').text(c.name);
    if (c.within_mag) svg.append('text').attr('x', W - 4).attr('y', y + 4).attr('class', 'mono').attr('text-anchor', 'end').text(d3.format('.0%')(c.within_mag.frac_up)); });
  svg.append('text').attr('x', W - 4).attr('y', 10).attr('text-anchor', 'end').attr('class', 'ttl').text('in MAGs');

  // B. genome traits vs ice preference
  const T = C.mag_tests.filter(t => t.rho_gc_family != null).sort((a, b) => b.rho_gc_family - a.rho_gc_family);
  const NAME = { arg_lys: 'Arg/(Arg+Lys)', pro: 'proline', acidic: 'acidic residues', charged: 'charged residues', gly: 'glycine', ser: 'serine', gln_thr: 'Gln+Thr', aliphatic: 'aliphatic index', gravy: 'GRAVY (hydropathy)', polar_unch: 'polar uncharged residues', genome_mb: 'genome size' };
  const b2 = fig(sec, 'Which genome traits go with an ice lifestyle?', `Spearman correlation between each MAG’s ice-vs-water transcript shift and a genome trait, after removing GC content and family effects (n = ${T[0].n_fe} MAGs with completeness ≥ 50%). Solid bars have q < 0.05. Gene counts are per Mbp.`);
  const W2 = width(b2), lw2 = Math.min(260, W2 * 0.42), rh2 = 18, H2 = T.length * rh2 + 48;
  const s2 = d3.select(b2).append('svg').attr('width', W2).attr('height', H2);
  const x2 = d3.scaleLinear().domain([-0.4, 0.4]).range([lw2 + 10, W2 - 40]);
  s2.append('g').attr('class', 'axis').attr('transform', `translate(0,${T.length * rh2 + 6})`).call(d3.axisBottom(x2).ticks(8));
  s2.append('line').attr('x1', x2(0)).attr('x2', x2(0)).attr('y1', 0).attr('y2', T.length * rh2 + 6).attr('stroke', css('--ink-3'));
  T.forEach((t, j) => { const y = j * rh2 + 10, sig = t.q_gc_family < 0.05;
    s2.append('text').attr('x', lw2).attr('y', y + 4).attr('text-anchor', 'end').text(NAME[t.variable] || t.variable);
    s2.append('rect').attr('x', Math.min(x2(0), x2(t.rho_gc_family))).attr('y', y - 5).attr('width', Math.abs(x2(t.rho_gc_family) - x2(0))).attr('height', 10).attr('rx', 2)
      .attr('fill', sig ? css(t.rho_gc_family > 0 ? '--exp' : '--ctl') : 'none').attr('stroke', css(t.rho_gc_family > 0 ? '--exp' : '--ctl'))
      .on('mousemove', ev => showTip(ev, `<b>${NAME[t.variable] || t.variable}</b><br>ρ = ${d3.format('.2f')(t.rho_gc_family)} (GC and family removed), q = ${d3.format('.2g')(t.q_gc_family)}<br>raw ρ = ${d3.format('.2f')(t.rho_raw)}`)).on('mouseleave', hideTip); });
  s2.append('text').attr('x', W2 - 40).attr('y', H2 - 2).attr('text-anchor', 'end').text('← water-leaning MAGs   ρ   ice-leaning MAGs →');

  // C. ice-binding proteins
  const B = C.ibp, pv = C.ibp_prevalence;
  sec.append(el('h3', {}, 'Ice-binding proteins: who has them, who uses them, and do they move?'));
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  const ice = C.categories.find(c => c.name === 'Ice-binding proteins');
  [[`${d3.format('.1f')(2 ** ice.ice.median_log2)}×`, `IBP transcripts in ice vs paired water (${Math.round(ice.ice.frac_up * ice.ice.n)} of ${ice.ice.n} pairs)`],
   [`${pv.ice}%`, `of ice-leaning MAGs carry an IBP gene, vs ${pv.neutral}% neutral and ${pv.water}% water-leaning (odds ratio ${pv.odds_ratio_ice_vs_rest})`],
   [`${B.mobility.ice_transcripts_mobile_pct}%`, `of IBP transcripts in ice come from mobile-associated genes (${B.mobility.ice_transcripts_mge_pct}% from genes on MGEs)`],
   [`${B.transfers.length}`, `flank-supported recent IBP transfers between genera or phyla; none expressed above 0.1 CPM`]].forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  const g3 = el('div', { class: 'grid2' }); sec.append(g3);
  const t1 = fig(g3, 'Main IBP transcript sources in ice', `Share of ice IBP transcripts. Eukaryotic IBPs (mostly centric diatoms) share a median ${d3.format('.0%')(B.euk_nearest_bacterial_aa.median_ident)} amino-acid identity with their closest bacterial IBP here: an ancient transfer, not a recent exchange.`);
  t1.innerHTML = `<div class="tablewrap" style="max-height:360px"><table><thead><tr><th>Source</th><th>Domain</th><th class="num">Genes</th><th class="num">CPM ice</th><th class="num">CPM water</th><th class="num">Share in ice</th></tr></thead><tbody>${
    B.taxa.map(t => `<tr><td>${t.domain === 'eukaryote' && /Kaiju/.test(t.genus) ? 'unassigned eukaryotic contigs' : t.genus}</td><td>${t.domain}</td><td class="num">${t.n}</td><td class="num">${d3.format('.1f')(t.cpm_SI)}</td><td class="num">${d3.format('.1f')(t.cpm_SW)}</td><td class="num">${d3.format('.1f')(t.share_SI)}%</td></tr>`).join('')}</tbody></table></div>`;
  const t2 = fig(g3, 'Candidate recent IBP transfers', `Pairs of IBP genes at ≥90% amino-acid identity on different contigs (${B.sharing_funnel.pairs_ge90aa_diff_contig} pairs), kept when a ≥95% nucleotide block is shared but the surrounding backgrounds differ.`);
  t2.innerHTML = `<div class="tablewrap" style="max-height:360px"><table><thead><tr><th>Host A</th><th>Host B</th><th class="num">aa id.</th><th class="num">Shared block</th><th>Level</th></tr></thead><tbody>${
    B.transfers.map(t => `<tr><td>${t.q_host}</td><td>${t.t_host}</td><td class="num">${d3.format('.0%')(t.aa_ident)}</td><td class="num">${d3.format(',')(t.shared95_bp)} bp</td><td>${t.cross}</td></tr>`).join('')}</tbody></table></div>`;
  const M = C.mpibp;
  sec.append(el('p', { class: 'note' }, `<b>The Marinomonas ice adhesin.</b> <i>${M.species}</i> (${M.mag.replace('dastool-', '')}) is one of the ice-enriched genomes here, and its species is known for MpIBP, a ~1.5 MDa calcium-dependent RTX adhesin that binds the bacterium to ice. The assembly carries a homolog on ${M.contig}: several short ORFs with 45–97% identity to fragments of the reference protein (UniProt A1YIY2/A1YIY3). The native gene’s long tandem repeats are collapsed and frameshifted in the long-read assembly. It is transcribed at about 0.3 CPM in both ice and water, so it is present but not induced in these samples.`));
}
function popgenCharts(sec, G) {
  sec.append(el('p', { class: 'kicker' }, '55 CMO Illumina metagenomes (2 M read pairs each; 2024-25 top/bottom ice with same-day seawater, and 2025-26 ice/water pairs) were mapped competitively to all 911 MAGs (bwa mem, MAPQ ≥ 20, ≥ 95% identity, ≥ 80% of the read aligned). Four libraries flagged in the metadata QC (three failed, one pool mismatch) were left out. Allele frequencies come from bcftools pileups of 72 well-covered MAGs.'));
  // A. species niches in DNA
  const GEN = ['Cognaticolwellia', 'Colwellia', 'Glaciecola', 'Polaribacter', 'Sulfitobacter', 'Yoonia'];
  const A = Object.entries(G.abundance).map(([b, v]) => ({ bin: b.replace('dastool-', ''), ...v })).filter(a => GEN.includes(a.genus) && a.log2_ice_water != null);
  const box = fig(sec, 'Ecotypes are species, not strains: sibling MAGs in ice and water DNA', 'Median paired log₂ ratio of each MAG’s share of mapped reads: ice vs same-day seawater (left), top vs bottom ice (right). Number = pairs in which the MAG was higher in ice (or top ice). Cognaticolwellia and Glaciecola are ice specialists, enriched toward the top. Their Colwellia relatives (83–84% ANI) stay with the water and bottom ice, matching the transcript shares above.');
  const W = width(box), lw = Math.min(280, W * 0.4), rh = 20, pw = (W - lw - 30) / 2, rows = GEN.flatMap(g => A.filter(a => a.genus === g).sort((a, b) => b.log2_ice_water - a.log2_ice_water));
  const H = rows.length * rh + 46, svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  let y0 = 28;
  [['log2_ice_water', 'ice − water', 0], ['log2_top_bot', 'top − bottom ice', 1]].forEach(([k, t, i]) => {
    const x0 = lw + 10 + i * (pw + 10), x = d3.scaleLinear().domain([-2.5, 2.5]).range([x0, x0 + pw - 40]);
    svg.append('text').attr('class', 'ttl').attr('x', x(0)).attr('y', 12).attr('text-anchor', 'middle').text(t);
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${rows.length * rh + y0})`).call(d3.axisBottom(x).ticks(5));
    svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 18).attr('y2', rows.length * rh + y0).attr('stroke', css('--ink-3'));
    rows.forEach((a, j) => { const y = j * rh + y0 + 6, v = a[k]; if (v == null) return; const cv = Math.max(-2.5, Math.min(2.5, v));
      svg.append('rect').attr('x', Math.min(x(0), x(cv))).attr('y', y - 5).attr('width', Math.abs(x(cv) - x(0))).attr('height', 10).attr('rx', 2).attr('fill', css(v > 0 ? '--exp' : '--ctl'))
        .on('mousemove', ev => showTip(ev, `<b>${a.label}</b><br>${t}: ${sgn(v)}${k === 'log2_ice_water' ? ` (higher in ice in ${a.n_pairs_up} of ${a.n_pairs} pairs)` : ''}<br>mean ${d3.format('.2f')(d3.mean(a.pct))}% of reads`)).on('mouseleave', hideTip);
      if (k === 'log2_ice_water') svg.append('text').attr('x', x0 + pw - 36).attr('y', y + 4).attr('class', 'mono').text(`${a.n_pairs_up}/${a.n_pairs}`); });
  });
  let prev = '';
  rows.forEach((a, j) => { const y = j * rh + y0 + 6;
    if (a.genus !== prev) { svg.append('line').attr('x1', 0).attr('x2', W).attr('y1', y - rh / 2).attr('y2', y - rh / 2).attr('stroke', css('--rule')); prev = a.genus; }
    svg.append('text').attr('x', lw).attr('y', y + 4).attr('text-anchor', 'end').text(a.label.length > 44 ? a.label.slice(0, 43) + '…' : a.label); });

  // B. FST: ice vs water and top vs bottom against a water-vs-water baseline, and pool differentiation
  const P = G.paired.slice().sort((a, b) => b.excess - a.excess);
  const b2 = fig(sec, 'Within a species, ice and water share one population', 'For each species with enough coverage (≥ 1,000 shared variable sites per pair, 2024-25), the median FST between ice and seawater sampled the same day (orange) and between top and bottom of the same core (green), against the median FST between seawater samples from consecutive sampling days in the same pool (grey, the sampling-noise baseline). Purple: FST between the control and diesel pools (ice, pooled). For 22 of 23 species ice-vs-water differentiation does not exceed the baseline.');
  const W2 = width(b2), lw2 = Math.min(280, W2 * 0.42), rh2 = 20, H2 = P.length * rh2 + 40;
  const s2 = d3.select(b2).append('svg').attr('width', W2).attr('height', H2);
  const x2 = d3.scaleLinear().domain([0, Math.max(0.4, d3.max(P, p => p.pool || 0))]).nice().range([lw2 + 10, W2 - 16]);
  s2.append('g').attr('class', 'axis').attr('transform', `translate(0,${P.length * rh2 + 6})`).call(d3.axisBottom(x2).ticks(6));
  s2.append('text').attr('x', W2 - 16).attr('y', H2 - 2).attr('text-anchor', 'end').text('FST');
  const MK = [['ww', '--other', 'water vs water (next day)'], ['iw', '--exp', 'ice vs water (same day)'], ['tb', '--s3', 'top vs bottom ice'], ['pool', '--s7', 'control vs diesel pool']];
  P.forEach((p, j) => { const y = j * rh2 + 12;
    s2.append('text').attr('x', lw2).attr('y', y + 4).attr('text-anchor', 'end').text(p.label);
    s2.append('line').attr('x1', x2(Math.min(p.ww, p.iw))).attr('x2', x2(Math.max(p.ww, p.iw))).attr('y1', y).attr('y2', y).attr('stroke', css('--rule')).attr('stroke-width', 3);
    MK.forEach(([k, c, l]) => { if (p[k] == null) return;
      s2.append('circle').attr('cx', x2(p[k])).attr('cy', y).attr('r', 5).attr('fill', css(c)).attr('stroke', css('--panel'))
        .on('mousemove', ev => showTip(ev, `<b>${p.label}</b><br>${MK.map(([kk, , ll]) => `${ll}: ${p[kk] == null ? '–' : d3.format('.3f')(p[kk])}`).join('<br>')}<br>π ice ${p.pi_ice ?? '–'}, water ${p.pi_sw ?? '–'}<br>${p.n_iw} ice–water pairs, ${p.n_ww} water–water pairs`)).on('mouseleave', hideTip); }); });
  b2.parentNode.insertBefore(el('div', { class: 'legend' }, MK.map(([, c, l]) => `<span><i style="background:var(${c});height:10px;width:10px;border-radius:50%"></i>${l}</span>`).join('')), b2);
  const top = P[0];
  sec.append(el('p', { class: 'note' }, `<b>One possible ice ecotype within a species.</b> ${top.label} (Flavobacteriaceae) is the only species whose same-day ice–water FST (${d3.format('.2f')(top.iw)}) clearly exceeds its water–water baseline (${d3.format('.2f')(top.ww)}). Its nucleotide diversity is lower in ice (π ${top.pi_ice} vs ${top.pi_sw}), and in the control pool on days 16–20 its ice and water populations differ at 350–580 nearly fixed sites. Ice seems to hold a subset of the strains in the water. All other well-covered species, including 20 SAR11-clade MAGs, Sulfurimonas and UBA9145, are one population across ice and water.`));
  // C. pools hold different strains
  const PO = G.pool.filter(p => p.medium === 'ice').sort((a, b) => b.fst - a.fst).slice(0, 12);
  const b3 = fig(sec, 'The two pools hold different strains', 'Pooled FST between control-pool and diesel-pool ice (2024-25), for the most differentiated MAGs with ≥ 200 variable sites, next to the same MAG’s ice-vs-water FST. Several species were founded by different strains in the two pools. Part of any diesel-vs-control contrast is therefore a strain difference that predates dosing.');
  b3.innerHTML = `<div class="tablewrap" style="max-height:340px"><table><thead><tr><th>MAG</th><th class="num">FST control vs diesel pool</th><th class="num">FST ice vs water</th><th class="num">Sites</th></tr></thead><tbody>${
    PO.map(p => `<tr><td>${p.label}</td><td class="num">${d3.format('.2f')(p.fst)}</td><td class="num">${p.iw == null ? '–' : d3.format('.2f')(p.iw)}</td><td class="num">${fmt(p.sites)}</td></tr>`).join('')}</tbody></table></div>`;
}

RENDER.ice = iceSec;
