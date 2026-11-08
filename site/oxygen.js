/* Oxygen section: is there transcriptional evidence that the diesel pools ran short of oxygen? Data: oxygen.json (oxygen/scripts/oxygen_site.py).
   Uses helpers from app.js (el, fig, width, css, showTip, hideTip, fmt, fmt1, legendPools). */
const oxSgn = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.2f')(v);
const oxSgn1 = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.1f')(v);
const oxP = p => p == null ? '' : p < 0.001 ? 'p < 0.001' : `p = ${d3.format('.2g')(p)}`;
const oxCt = (J, c, i) => J.contrasts.find(x => x.contrast === c && x.index === i) || {};
const OX_IDX = {
  comm_HA_LA: 'Community: high- ÷ low-affinity oxidases',
  within_HA_LA: 'Within genome: high- ÷ low-affinity oxidases (genome-centred)',
  composition_HA_logodds: 'Expected from genome mix alone',
  regulation_HA_logodds: 'Left over after genome mix (regulation)',
  comm_anaer_cpm: 'Community: anaerobic and microaerobic markers (log₂ CPM)',
};

function oxygenSec(sec) {
  const J = D.oxygen;
  if (!J) { sec.innerHTML = '<h2>Oxygen</h2><p class="kicker">oxygen.json is missing.</p>'; return; }
  const b25 = i => oxCt(J, 'BACI SW 2025-26', i).mean, d25 = i => oxCt(J, 'diesel minus control SW 2025-26', i);
  const HA = J.carriers.HA, cy = HA.rows.find(r => r.bin === 'lorbin_184');
  sec.innerHTML = `<h2>Oxygen</h2>
  <p class="kicker">Several pages pointed at falling oxygen in the diesel pool. In 2025-26 the high-affinity terminal oxidases rose (cbb3 and bd, which bacteria switch to when O₂ is scarce), along with Cycloclasticus’s nitrate reductase and bidirectional hydrogenase. No dissolved oxygen was logged, so this page asks the transcripts directly. It builds a low-oxygen expression index two ways. The <b>community</b> index is the ratio of high- to low-affinity oxidase transcripts (cbb3 + bd vs aa3 + bo3) in the whole sample, plus anaerobic markers. The <b>within-genome</b> index compares each genome only with itself: for the ${J.genome_count || 94} genomes that carry both oxidase types, it takes the same ratio, centres it on that genome’s own average and averages across genomes. A bloom of a cbb3-rich organism moves the first index but not the second. The answer is that the 2025-26 signal is almost entirely a change in who is there, not cells responding to low oxygen.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  [[oxSgn1(b25('comm_HA_LA')), 'log₂ BACI, community high- ÷ low-affinity oxidases, 2025-26 seawater'],
   [oxSgn1(b25('composition_HA_logodds')), 'log₂ of that shift expected from the change in genome mix alone'],
   [oxSgn1(b25('within_HA_LA')), `log₂ BACI, within-genome index: no regulatory shift (diesel − control ${oxSgn1(d25('within_HA_LA').mean)}, ${oxP(d25('within_HA_LA').p)})`],
   [d3.format('.0%')(cy.oil / HA.total_oil), `of the diesel pool’s high-affinity oxidase transcripts in 2025-26 come from one blooming genome, Cycloclasticus lorbin_184`]]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  oxDecomp(sec, J);
  oxTime(sec, J);
  oxCarriers(sec, J);
  oxIce(sec, J);
  oxOther(sec, J);
  oxReading(sec, J);
}

/* ---------- decomposition table ---------- */
function oxDecomp(sec, J) {
  sec.append(el('h3', {}, 'Composition or regulation?'));
  sec.append(el('p', { class: 'kicker' }, 'For binned genomes, the community’s high-affinity share can be split into two parts. One is what it would be if every genome kept its own average oxidase mix and only abundances changed (composition). The other is what is left (regulation). Diesel minus control is a median difference after dosing; BACI also subtracts each pool’s own pre-dosing level. In 2024-25 that is a single pre-dosing sample per pool, and the pools were filled two weeks apart, so 2024-25 differences may be pool history.'));
  const rows = ['comm_HA_LA', 'comm_cbb3_aa3', 'comm_bd_aa3', 'composition_HA_logodds', 'regulation_HA_logodds', 'within_HA_LA', 'within_HA_LA_noCyclo', 'comm_anaer_cpm', 'comm_hemN_hemF', 'comm_nrdDG_nrdAB'];
  const NAME = { ...OX_IDX, comm_cbb3_aa3: 'Community: cbb3 ÷ aa3', comm_bd_aa3: 'Community: bd ÷ aa3', within_HA_LA_noCyclo: 'Within genome, Cycloclasticus excluded',
    comm_hemN_hemF: 'Community: O₂-independent ÷ O₂-dependent haem enzyme (hemN ÷ hemF)', comm_nrdDG_nrdAB: 'Community: anaerobic ÷ aerobic ribonucleotide reductase (nrdDG ÷ nrdAB)' };
  const col = v => { if (v == null) return 'transparent'; const t = Math.max(-1, Math.min(1, v / 2)); return t >= 0 ? d3.interpolateRgb(css('--panel'), css('--exp'))(t * 0.85) : d3.interpolateRgb(css('--panel'), css('--ctl'))(-t * 0.85); };
  const cell = (c, i) => { const x = oxCt(J, c, i); return `<td class="num" style="background:${col(x.mean)}">${oxSgn(x.mean)}${x.p != null && x.p < 0.05 ? '<span style="color:var(--ink-3)">*</span>' : ''}</td>`; };
  const box = fig(sec, 'Low-oxygen indices: diesel effect, log₂', 'Positive means more high-affinity oxidase (or more anaerobic gene expression) in the diesel pool. * marks Mann-Whitney p < 0.05 for diesel vs control after dosing (BACI has no test). Colour runs blue (lower in diesel) to orange (higher).');
  box.innerHTML = `<div class="tablewrap" style="max-height:none"><table><thead><tr><th>Index</th><th class="num">Diesel − control 24-25</th><th class="num">BACI 24-25</th><th class="num">Diesel − control 25-26</th><th class="num">BACI 25-26</th></tr></thead><tbody>${
    rows.map(i => `<tr><td style="white-space:normal">${NAME[i]}</td>${cell('diesel minus control SW 2024-25', i)}${cell('BACI SW 2024-25', i)}${cell('diesel minus control SW 2025-26', i)}${cell('BACI SW 2025-26', i)}</tr>`).join('')}</tbody></table></div>`;
}

/* ---------- index time courses (signed values, own chart) ---------- */
function oxTime(sec, J) {
  sec.append(el('h3', {}, 'Through the season'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'oxsel', 'aria-label': 'Index' }); Object.entries(OX_IDX).forEach(([k, v]) => sel.append(new Option(v, k)));
  c.append(el('label', { for: 'oxsel' }, 'Index'), sel);
  const v = el('div'); sec.append(v);
  function draw() {
    v.innerHTML = ''; const k = sel.value;
    const box = fig(v, OX_IDX[k], k === 'within_HA_LA' ? 'Zero is each genome’s own average. In 2025-26 the index drifts upward through the season in both pools alike (Spearman ρ 0.53 control, 0.68 diesel water), so the drift is seasonal, not a diesel effect.' : 'Per-day means of replicate libraries. Log₂ scale.');
    legendPools(box.parentNode.insertBefore(el('div'), box));
    const W = width(box), seasons = ['2024-25', '2025-26'], cols = W < 620 ? 1 : 2, gap = 28, pw = cols === 2 ? (W - gap) / 2 : W - 10, ph = 200;
    const svg = d3.select(box).append('svg').attr('width', W).attr('height', (ph + 30) * (cols === 2 ? 1 : 2));
    const pts = J.samples.filter(s => ['CTL', 'EXP'].includes(s.pool) && ['SW', 'SI'].includes(s.medium) && s.day != null && s[k] != null);
    const y = d3.scaleLinear().domain(d3.extent(pts, s => s[k])).nice().range([ph - 30, 8]);
    seasons.forEach((se, i) => {
      const gx = cols === 2 ? i * (pw + gap) : 0, gy = cols === 2 ? 0 : i * (ph + 30);
      const g = svg.append('g').attr('transform', `translate(${gx + 46},${gy + 12})`), iw = pw - 56;
      const P = pts.filter(s => s.season === se), x = d3.scaleLinear().domain(d3.extent(P, s => s.day)).nice().range([0, iw]);
      g.append('text').attr('x', 0).attr('y', 0).style('font-weight', 650).style('font-size', '12px').text(`CMO ${se}`);
      g.append('g').attr('class', 'axis').attr('transform', `translate(0,${ph - 30})`).call(d3.axisBottom(x).ticks(6));
      g.append('g').attr('class', 'axis').call(d3.axisLeft(y).ticks(5));
      if (y.domain()[0] < 0 && y.domain()[1] > 0) g.append('line').attr('x1', 0).attr('x2', iw).attr('y1', y(0)).attr('y2', y(0)).attr('stroke', css('--ink-3')).attr('stroke-dasharray', '2 3');
      const od = se === '2024-25' ? 5.5 : 1; g.append('line').attr('x1', x(od)).attr('x2', x(od)).attr('y1', 8).attr('y2', ph - 30).attr('stroke', css('--ink-3')).attr('stroke-dasharray', '3 3');
      d3.groups(P, s => s.pool + '|' + s.medium).forEach(([key, a]) => {
        const [pool, med] = key.split('|');
        const byDay = d3.rollups(a, q => ({ v: d3.mean(q, s => s[k]), n: q.map(s => s.name).join(', ') }), s => s.day).sort((p, q) => p[0] - q[0]);
        const colr = css(pool === 'EXP' ? '--exp' : '--ctl');
        g.append('path').datum(byDay).attr('fill', 'none').attr('stroke', colr).attr('stroke-width', 2).attr('stroke-dasharray', med === 'SI' ? '5 3' : null)
          .attr('d', d3.line().x(d => x(d[0])).y(d => y(d[1].v)));
        g.selectAll(null).data(byDay).join('circle').attr('cx', d => x(d[0])).attr('cy', d => y(d[1].v)).attr('r', 3.5).attr('fill', colr)
          .on('mousemove', (ev, d) => showTip(ev, `<b>${d[1].n}</b><br>day ${d[0]}: ${oxSgn(d[1].v)}`)).on('mouseleave', hideTip);
      });
    });
  }
  sel.onchange = draw; sel.value = 'within_HA_LA'; draw();
}

/* ---------- who carries the oxidases ---------- */
function oxCarriers(sec, J) {
  sec.append(el('h3', {}, 'Why the community ratio rose in 2025-26'));
  const cy = J.cyclo, f = k => cy.find(x => x.gene === k);
  sec.append(el('p', { class: 'kicker' }, `Cycloclasticus lorbin_184 carries both cbb3 and aa3 oxidases and keeps a nearly constant mix of the two. Its cbb3 is ${d3.format('.2f')(f('cbb3').oil / f('aa3').oil)} of its aa3 in the diesel pool against ${d3.format('.2f')(f('cbb3').ctl / f('aa3').ctl)} in the control and ${d3.format('.2f')(f('cbb3').pre / f('aa3').pre)} before dosing. That mix holds more high-affinity oxidase than the aa3-heavy roseobacters (Sulfitobacter) and UBA9145 that dominate the control. When it bloomed, the community ratio moved with it. Its bidirectional Hox hydrogenase is a smaller share of its transcripts after dosing (${f('hox').oil}‰) than before (${f('hox').pre}‰), and its nitrate reductase stays below 1‰. Neither looks like a cell switching to low-oxygen metabolism. The rise in community anaerobic markers (BACI ${oxSgn1(oxCt(J, 'BACI SW 2025-26', 'comm_anaer_cpm').mean)} log₂) has the same cause: mostly Cycloclasticus genes carried along by the bloom.`));
  const box = fig(sec, 'Main carriers of oxidase transcripts, 2025-26 seawater (median CPM after dosing)', 'High-affinity = cbb3 + bd; low-affinity = aa3 + bo3.');
  const tb = nm => `<table><thead><tr><th>${nm === 'HA' ? 'High-affinity' : 'Low-affinity'} carriers</th><th class="num">Diesel</th><th class="num">Control</th></tr></thead><tbody>${
    J.carriers[nm].rows.map(r => `<tr><td>${r.taxon === r.bin ? r.bin : `${r.taxon} <span class="mono">${r.bin}</span>`}</td><td class="num">${fmt(r.oil)}</td><td class="num">${fmt(r.ctl)}</td></tr>`).join('')}<tr><td><b>All genes</b></td><td class="num"><b>${fmt(J.carriers[nm].total_oil)}</b></td><td class="num"><b>${fmt(J.carriers[nm].total_ctl)}</b></td></tr></tbody></table>`;
  box.innerHTML = `<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:16px"><div class="tablewrap">${tb('HA')}</div><div class="tablewrap">${tb('LA')}</div></div>
    <div class="tablewrap" style="margin-top:12px"><table><thead><tr><th>Cycloclasticus lorbin_184, ‰ of its own transcripts</th><th class="num">Before dosing</th><th class="num">Diesel</th><th class="num">Control</th></tr></thead><tbody>${
    cy.map(r => `<tr><td>${r.gene}</td><td class="num">${r.pre ?? '–'}</td><td class="num">${r.oil ?? '–'}</td><td class="num">${r.ctl ?? '–'}</td></tr>`).join('')}</tbody></table></div>`;
  const gs = J.genome_shift;
  if (gs.length) sec.append(el('p', { class: 'note' }, `Single genomes do not shift consistently either. The few with a detectable diesel effect on their own oxidase mix go in both directions: ${gs.map(g => `${g.taxon} ${g.bin} ${g.season} ${g.medium} ${oxSgn1(g.shift)} (q ${d3.format('.2g')(g.q)})`).join('; ')}.`));
}

/* ---------- ice ---------- */
function oxIce(sec, J) {
  const iw = i => oxCt(J, 'ice minus water', i), tb = i => oxCt(J, 'top minus bottom ice', i);
  sec.append(el('h3', {}, 'Ice and brine'));
  sec.append(el('p', { class: 'kicker' }, `Brine channels are often low in oxygen, and ice communities do carry more high-affinity oxidase. The community ratio is ${oxSgn1(iw('comm_HA_LA').mean)} log₂ higher in ice than in the paired water (${iw('comm_HA_LA').n_up} of ${iw('comm_HA_LA').n} pairs, ${oxP(iw('comm_HA_LA').p)}). Most of that is again composition (${oxSgn1(iw('composition_HA_logodds').mean)}, ${iw('composition_HA_logodds').n_up} of ${iw('composition_HA_logodds').n}): ice genera such as Colwellia carry more cbb3 and bd. Within genomes the shift is small${iw('within_HA_LA').p < 0.05 ? ', though nominally significant' : ' and not significant'} (${oxSgn(iw('within_HA_LA').mean)}, ${iw('within_HA_LA').n_up} of ${iw('within_HA_LA').n}, ${oxP(iw('within_HA_LA').p)}). Top ice, with the coldest and saltiest brine, leans slightly more towards high-affinity oxidases within genomes than bottom ice once Cycloclasticus is excluded (${oxSgn(tb('within_HA_LA_noCyclo').mean)}, ${tb('within_HA_LA_noCyclo').n_up} of ${tb('within_HA_LA_noCyclo').n} pairs, ${oxP(tb('within_HA_LA_noCyclo').p)}). Anaerobic markers are lower in ice than in water (${oxSgn1(iw('comm_anaer_cpm').mean)} log₂, ${oxP(iw('comm_anaer_cpm').p)}).`));
  const cv = J.covariates.filter(x => /ice_thickness|core_|roof/.test(x.covariate) && /within|comm_HA/.test(x.index));
  const L = { log_cyclo: 'Cycloclasticus TPM (log)', water_temp_c: 'water temperature', ice_thickness_cm: 'ice thickness (sampling log)', core_thickness_cm: 'ice thickness (cores)', core_T_top: 'top-ice temperature', core_Vb_top: 'top-ice brine volume', core_Sbrine_top: 'top-ice brine salinity', 'roof_open Yes vs No': 'roof open vs closed (difference)', day: 'experiment day' };
  if (cv.length) { const box = fig(sec, 'Physical covariates with p < 0.05', 'Spearman ρ (or median difference for roof open vs closed). Ice physics comes from the nearest ice core of the same pool (±2 days). These correlations share a seasonal trend: as winter advances the ice thickens, gets colder at the top and its brine saltier, so they cannot be separated from time.');
    box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Samples</th><th>Covariate</th><th>Index</th><th class="num">n</th><th class="num">ρ</th><th class="num">p</th></tr></thead><tbody>${
      cv.map(x => `<tr><td>${x.subset}</td><td>${L[x.covariate] || x.covariate}</td><td>${x.index.startsWith('within') ? 'within genome' + (x.index.endsWith('noCyclo') ? ', no Cycloclasticus' : '') : 'community'}</td><td class="num">${x.n}</td><td class="num">${oxSgn(x.rho)}</td><td class="num">${d3.format('.2g')(x.p)}</td></tr>`).join('')}</tbody></table></div>`; }
}

/* ---------- other markers + oil snow ---------- */
function oxOther(sec, J) {
  sec.append(el('h3', {}, 'Oxygen-independent enzymes, and the oil-snow sample'));
  const g = n => J.groups.find(x => x.group.startsWith(n)), row = n => J.mos.rows.find(x => x.group.startsWith(n));
  sec.append(el('p', { class: 'kicker' }, `Some enzymes come in an oxygen-requiring and an oxygen-independent version. Bacteria switch to the latter when O₂ runs out, so these pairs are a second check. In 2025-26 the anaerobic ribonucleotide reductase (nrdDG) shows no consistent diesel signal. It is lower in the diesel pool than in the control after dosing, but the diesel pool started even lower. The oxygen-independent haem enzyme hemN does not rise either. The one marine-oil-snow sample (CMO_EXP_DAY04_MOS) carries somewhat more fermentation genes than the water and ice taken the same day. Pyruvate formate-lyase is ${row('pflB').mos} vs ${row('pflB').SW} and ${row('pflB').SI} CPM, and adhE ${row('adhE').mos} vs ${row('adhE').SW} and ${row('adhE').SI}. Its oxidase indices are unremarkable. That fits low-oxygen pockets inside oil aggregates, but it is a single library.`));
  const box = fig(sec, 'Day 4 of 2025-26 in the diesel pool: oil snow vs water and ice (CPM)', '');
  box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Gene group</th><th class="num">Oil snow</th><th class="num">Water</th><th class="num">Ice</th></tr></thead><tbody>${
    J.mos.rows.map(r => `<tr><td>${r.group}</td><td class="num">${fmt1(r.mos)}</td><td class="num">${fmt1(r.SW)}</td><td class="num">${fmt1(r.SI)}</td></tr>`).join('')}</tbody></table></div>`;
}

/* ---------- reading ---------- */
function oxReading(sec, J) {
  const d24 = oxCt(J, 'diesel minus control SW 2024-25', 'within_HA_LA'), b24 = oxCt(J, 'BACI SW 2024-25', 'within_HA_LA');
  sec.append(el('h3', {}, 'What the index can and cannot show'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li><b>Can show:</b> whether cells change their own oxidase mix, the classic bacterial response to falling O₂. In 2025-26 they did not respond to diesel. The within-genome index is flat against each pool’s baseline (BACI ${oxSgn(oxCt(J, 'BACI SW 2025-26', 'within_HA_LA').mean)}), and the rise in high-affinity oxidases on the Sulfur & hydrogen page is the Cycloclasticus bloom. The earlier reading of “falling oxygen” is withdrawn on that evidence.</li>
    <li><b>2024-25:</b> the within-genome index is ${oxSgn(d24.mean)} log₂ higher in the diesel pool than in the control after dosing (${oxP(d24.p)}). But the diesel pool’s single pre-dosing sample was already higher (BACI ${oxSgn(b24.mean)}), and the pools were filled two weeks apart, so this is most likely pool history.</li>
    <li><b>Cannot show:</b> actual O₂ concentrations. Aerobic bacteria can keep the same oxidase mix down to quite low O₂. The index also misses anoxic micro-niches that hold few cells, such as oil droplets, aggregates and brine pockets. And it only covers genomes that carry both oxidase types (about three-quarters of oxidase transcripts are binned). No change in the index does not mean O₂ stayed at saturation.</li>
    <li><b>Recommendation for next season:</b> a logging O₂ optode under the ice in each pool, read at least daily, plus O₂ in melted core sections and in any oil snow collected. That is cheap and would settle the question. The transcripts give no reason to expect a large diesel-driven O₂ drop, but they cannot rule out a modest one. With 136 L of diesel being degraded under ice, it is worth measuring.</li>`));
}

RENDER.oxygen = oxygenSec;
