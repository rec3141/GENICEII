/* Metals & toxicity section: metal-resistance and homeostasis genes, solvent/oxidative stress for comparison, their hosts, ice/brine
   patterns and mobility. Data: metals.json (metals/scripts/mbuild.py). Uses helpers from app.js (el, fig, width, css, showTip, hideTip,
   fmt, fmt1, timecourse, groupSeries, legendPools, POOL). */
const MT_GROUPS = [['toxic metal', 'Toxic-metal resistance'], ['micronutrient', 'Metal micronutrient homeostasis'],
  ['solvent stress', 'Solvent and membrane stress (diesel toxicity)'], ['oxidative stress', 'Oxidative stress'], ['proteotoxic stress', 'Protein-folding stress']];
const mtSgn = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.1f')(v);
const mtShort = n => n.replace(/ \(.*\)$/, '').replace(' + sensors', '');

function metalsSec(sec) {
  const J = D.metals;
  if (!J) { sec.innerHTML = '<h2>Metals & toxicity</h2><p class="kicker">metals.json is missing.</p>'; return; }
  const C = J.cats, byName = Object.fromEntries(C.map(c => [c.name, c]));
  const d26 = (c, k = 'baci') => (c.diesel.find(d => d.season === '2025-26' && d.medium === 'SW') || {})[k];
  const tm = C.filter(c => c.group === 'toxic metal');
  const mob = d3.sum(tm, c => c.mob.mobile_assoc * c.n_conf) / d3.sum(tm, c => c.n_conf);
  sec.innerHTML = `<h2>Metals & toxicity</h2>
  <p class="kicker">Diesel is a refined fuel, low in the metals that make crude oil toxic (vanadium and nickel porphyrins), and no metals or other trace elements were measured in the pools: the CMO logs record only that nutrient samples were taken. This page therefore reads metal exposure indirectly, from the genes the communities transcribe. It covers resistance to toxic metals, uptake and storage of metal micronutrients, and, for comparison, the membrane, oxidative and protein-folding stress responses that hydrocarbons themselves cause. Expression is CPM of protein-coding read pairs (corrected sample labels). Gene sets use the coordinate-corrected KOs plus Bakta product names; the headline numbers use a confident tier that drops KO-only hits on bare hypothetical proteins.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  const zn = byName['Zn/Cd/Pb export (zntA, cadA, czcD) + sensors'], psp = byName['Phage-shock membrane stress (psp)'];
  [[`${d3.format('.1f')(2 ** d26(psp))}×`, 'phage-shock (psp) membrane-stress transcripts in the 2025-26 diesel pool, relative to its own pre-dosing level and the control'],
   [`${d3.format('.1f')(2 ** d26(zn))}×`, 'Zn/Cd/Pb export transcripts, same comparison; much of it carried by the Cycloclasticus bloom'],
   [(() => { const t = n => J.cats.find(c => c.name.startsWith(n)).ice; const hg = t('Mercury'), as = t('Arsenic'); return hg.tb_up === as.tb_up ? `${hg.tb_up} of ${hg.tb_n}` : `${hg.tb_up} and ${as.tb_up} of ${hg.tb_n}`; })(), 'top-vs-bottom ice pairs with more mercury and arsenic resistance transcripts at the colder, saltier top'],
   [`${d3.format('.0f')(mob)}% vs ${d3.format('.1f')(J.hc_ref.mobile_assoc)}%`, 'of toxic-metal resistance genes vs hydrocarbon genes sit on or next to mobile elements']]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  sec.append(el('h3', {}, 'What we found'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li><b>The diesel signature is membrane stress more than metal stress.</b> The strongest baseline-corrected diesel-pool responses in 2025-26 are phage-shock proteins (${mtSgn(d26(psp))} log₂) and the fatty-acid cis–trans isomerase that stiffens membranes against solvents (${mtSgn(d26(byName['Membrane cis-trans isomerase (cti)']))}, from a tiny base). Most of the psp rise comes from one water bacterium, UBA9145 vamb_4149 (Pseudohongiellaceae), and it is a response inside the cells: its psp transcripts reach 43,000 per million of its own transcripts in the 2025-26 diesel pool, against 2,200 before dosing and 2,600 in the control (19,400 vs 6,000 and 1,500 in 2024-25).</li>
    <li><b>Apparent metal responses mostly ride the oil-degrader bloom.</b> Zn/Cd/Pb export, copper, arsenic and mercury transcripts rise in the diesel pool, but in 2025-26 much of the rise is carried by Cycloclasticus lorbin_184, whose own metal-gene expression per host transcript is no higher after dosing than before. The degraders bring their constitutive metal defences with them.</li>
    <li><b>Some oil degraders switch on oxidative defences.</b> Within its own transcriptome, Oleispira expresses catalase, peroxidase and superoxide-dismutase genes at 6,600 per million host transcripts in the 2024-25 diesel pool against 720 before dosing (one library) and 740 in the control, and at 15,900 in the 2025-26 diesel pool. Cognaticolwellia lorbin_135 reaches 28,000 in the 2025-26 diesel pool. This fits reactive oxygen from oxygenase chemistry and hydrocarbon oxidation; the baselines are thin.</li>
    <li><b>Top ice looks metal-stressed.</b> Mercury and arsenic resistance and iron storage are higher in top than bottom ice in all 13 pairs, copper and tellurite in 12. Freezing concentrates solutes, metals included, into brine that reaches 61–166 PSU at the top. This fits that, but differences between the communities of the two layers could produce the same pattern.</li>
    <li><b>Metal resistance is the mobile part of the genome.</b> Mercury, Zn/Cd/Pb, Ni/Co and copper resistance genes are 2–3× as often on plasmids, islands or next to transposases as hydrocarbon genes are. ${J.plasmid_counts.both} plasmids carry both metal resistance and hydrocarbon genes, ${d3.format('.0%')(J.plasmid_counts.both / J.plasmid_counts.hc)} of the plasmids with oil genes: a route for co-selection.</li>`));
  mtDiesel(sec, J);
  mtCourse(sec, J);
  mtIce(sec, J);
  mtHosts(sec, J);
  mtDegraders(sec, J);
  mtMobility(sec, J);
  sec.append(el('p', { class: 'note warn' }, '<b>Caveats.</b> One pool per treatment, and the two pools hold different strains of several species (see Ice vs water), so diesel-vs-control contrasts mix treatment and pool. After the day-5 label correction the 2024-25 diesel pool has a single pre-dosing seawater library, so its baseline-corrected values rest on one sample. Gene transcription shows capacity and response, not metal concentrations, and KO-only assignments to broad families (TonB-dependent receptors, RND efflux, MerR regulators) also cover functions other than metal handling.'));
}

/* ---------- diesel response per category: raw and baseline-corrected, both winters ---------- */
function mtDiesel(sec, J) {
  sec.append(el('h3', {}, 'Diesel pool vs control'));
  const box = fig(sec, 'Change in the diesel pool after dosing, seawater', 'Filled dot: log₂ ratio of median CPM, diesel vs control pool after dosing. Hollow ring: the same ratio minus the pools’ ratio before dosing (before-after-control-impact), the fairer estimate where the pools already differed. 2024-25 has one pre-dosing diesel-pool library. Hover for medians.');
  const rows = MT_GROUPS.flatMap(([g, gl]) => [{ head: gl }, ...J.cats.filter(c => c.group === g)]);
  const W = Math.max(700, width(box)), lw = Math.min(310, W * 0.4), rh = 22, pw = (W - lw - 20) / 2, H = rows.length * rh + 52;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  ['2024-25', '2025-26'].forEach((s, i) => {
    const x0 = lw + 10 + i * (pw + 10), x = d3.scaleLinear().domain([-3, 3.5]).range([x0, x0 + pw - 16]);
    svg.append('text').attr('class', 'ttl').attr('x', x(0)).attr('y', 12).attr('text-anchor', 'middle').text(`CMO ${s}`);
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${rows.length * rh + 22})`).call(d3.axisBottom(x).ticks(5));
    svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 18).attr('y2', rows.length * rh + 22).attr('stroke', css('--ink-3'));
    rows.forEach((c, j) => { if (c.head) return; const y = j * rh + 30, d = c.diesel.find(z => z.season === s && z.medium === 'SW'); if (!d) return;
      const cl = v => x(Math.max(-3, Math.min(3.5, v)));
      const tip = ev => showTip(ev, `<b>${c.name}</b> · ${s} seawater<br>diesel ${fmt1(d.oil)} vs control ${fmt1(d.ctl)} CPM (log₂ ${mtSgn(d.lfc)})<br>${d.baci == null ? '' : `before dosing: diesel ${fmt1(d.pre_exp)} vs control ${fmt1(d.pre_ctl)} → baseline-corrected ${mtSgn(d.baci)}`}<br>${d.n_oil} diesel, ${d.n_ctl} control libraries; ${d.n_pre_exp} diesel pre-dosing`);
      if (d.baci != null) { svg.append('line').attr('x1', cl(d.lfc)).attr('x2', cl(d.baci)).attr('y1', y).attr('y2', y).attr('stroke', css('--rule')).attr('stroke-width', 2);
        svg.append('circle').attr('cx', cl(d.baci)).attr('cy', y).attr('r', 5).attr('fill', 'none').attr('stroke', css(d.baci > 0 ? '--exp' : '--ctl')).attr('stroke-width', 2).on('mousemove', tip).on('mouseleave', hideTip); }
      svg.append('circle').attr('cx', cl(d.lfc)).attr('cy', y).attr('r', 4.5).attr('fill', css(d.lfc > 0 ? '--exp' : '--ctl')).attr('stroke', css('--panel')).on('mousemove', tip).on('mouseleave', hideTip); });
  });
  rows.forEach((c, j) => { const y = j * rh + 30;
    if (c.head) svg.append('text').attr('x', 0).attr('y', y + 4).attr('class', 'ttl').text(c.head);
    else svg.append('text').attr('x', lw).attr('y', y + 4).attr('text-anchor', 'end').text(mtShort(c.name)); });
  svg.append('text').attr('x', W - 6).attr('y', H - 2).attr('text-anchor', 'end').text('← higher in control      log₂      higher in diesel pool →');
}

/* ---------- time course for a chosen category ---------- */
function mtCourse(sec, J) {
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'mtsel', 'aria-label': 'Gene set' });
  MT_GROUPS.forEach(([g, gl]) => { const og = document.createElement('optgroup'); og.label = gl; J.cats.filter(x => x.group === g).forEach(x => og.append(new Option(x.name, x.name))); sel.append(og); });
  c.append(el('label', { for: 'mtsel' }, 'Gene set'), sel);
  const holder = el('div'); sec.append(holder);
  const idx = Object.fromEntries(J.samples.map((s, i) => [s, i]));
  function draw() {
    holder.innerHTML = '';
    const v = J.series[sel.value], vals = D.samples.map(s => idx[s.id] == null ? null : v[idx[s.id]]);
    const box = fig(holder, sel.value, 'CPM of protein-coding read pairs per library; lines join per-day means. Dotted line: dosing.');
    legendPools(box.parentNode.insertBefore(el('div'), box));
    timecourse(box, groupSeries(vals), { ylab: 'CPM', unit: 'CPM', h: 190 });
  }
  sel.value = 'Phage-shock membrane stress (psp)'; sel.onchange = draw; draw();
}

/* ---------- ice: ice vs water and top vs bottom (paired) ---------- */
function mtIce(sec, J) {
  sec.append(el('h3', {}, 'In the ice: brine concentrates metals'));
  const cats = J.cats.slice().sort((a, b) => b.ice.tb - a.ice.tb);
  const box = fig(sec, 'Ice vs water and top vs bottom ice, paired', 'Each dot is one pair (same pool and day); the bar is the mean log₂ ratio, and the number counts pairs higher in ice (left) or in top ice (right). 24 ice–water pairs, 13 top–bottom pairs (2024-25). Top ice is colder, with saltier and smaller brine pockets (see Ice vs water).');
  const W = Math.max(700, width(box)), lw = Math.min(250, W * 0.34), rh = 22, pw = (W - lw - 20) / 2, H = cats.length * rh + 44;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  [['iw', 'ice − water'], ['tb', 'top − bottom ice']].forEach(([k, t], i) => {
    const x0 = lw + 10 + i * (pw + 10), x = d3.scaleLinear().domain([-2, 2]).range([x0, x0 + pw - 58]);
    svg.append('text').attr('class', 'ttl').attr('x', x(0)).attr('y', 12).attr('text-anchor', 'middle').text(t);
    svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${cats.length * rh + 22})`).call(d3.axisBottom(x).ticks(5));
    svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 18).attr('y2', cats.length * rh + 22).attr('stroke', css('--ink-3'));
    cats.forEach((c, j) => { const y = j * rh + 30, I = c.ice, mu = I[k];
      svg.append('rect').attr('x', Math.min(x(0), x(mu))).attr('y', y - 5).attr('width', Math.abs(x(mu) - x(0))).attr('height', 10).attr('rx', 2).attr('fill', css(mu > 0 ? '--exp' : '--ctl')).attr('opacity', .35);
      I[k + '_pairs'].forEach(v => v != null && svg.append('circle').attr('cx', x(Math.max(-2, Math.min(2, v)))).attr('cy', y).attr('r', 2.6).attr('fill', css(v > 0 ? '--exp' : '--ctl')).attr('opacity', .8));
      svg.append('rect').attr('x', x0).attr('y', y - rh / 2).attr('width', pw).attr('height', rh).attr('fill', 'transparent')
        .on('mousemove', ev => showTip(ev, `<b>${c.name}</b> · ${t}<br>mean ${mtSgn(mu)} log₂, higher in ${I[k + '_up']} of ${I[k + '_n']} pairs (q ${d3.format('.2g')(I[k + '_q'])})<br>median CPM 2024-25: water ${fmt1(I.cpm.SW)}, top ${fmt1(I.cpm.TOP)}, bottom ${fmt1(I.cpm.BOT)}`)).on('mouseleave', hideTip);
      svg.append('text').attr('x', x0 + pw - 40).attr('y', y + 4).attr('class', 'mono').text(`${I[k + '_up']}/${I[k + '_n']}`); });
  });
  cats.forEach((c, j) => svg.append('text').attr('x', lw).attr('y', j * rh + 34).attr('text-anchor', 'end').text(mtShort(c.name)));
}

/* ---------- hosts per category ---------- */
function mtHosts(sec, J) {
  sec.append(el('h3', {}, 'Who carries the transcripts'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'mthost', 'aria-label': 'Gene set' });
  J.cats.filter(x => x.hosts).forEach(x => sel.append(new Option(x.name, x.name))); c.append(el('label', { for: 'mthost' }, 'Gene set'), sel);
  const holder = el('div', { class: 'grid2' }); sec.append(holder);
  function draw() {
    holder.innerHTML = '';
    const x = J.cats.find(z => z.name === sel.value);
    ['2025-26', '2024-25'].forEach(s => {
      const box = fig(holder, `${s} seawater after dosing`, `Mean CPM of the set’s transcripts by genome (top six in each pool). Total: diesel pool ${fmt1(x.hosts['oil_' + s].total)} CPM, control ${fmt1(x.hosts['ctl_' + s].total)} CPM.`);
      const names = [...new Set([...x.hosts['oil_' + s].top.map(h => h.host), ...x.hosts['ctl_' + s].top.slice(0, 3).map(h => h.host)])];
      const get = (k, n) => (x.hosts[k + '_' + s].top.find(h => h.host === n) || {}).cpm;
      box.innerHTML = `<div class="tablewrap" style="max-height:300px"><table><thead><tr><th>Genome</th><th class="num">Diesel CPM</th><th class="num">Control CPM</th></tr></thead><tbody>${
        names.map(n => `<tr><td>${n}</td><td class="num">${get('oil', n) == null ? 'below top 6' : fmt1(get('oil', n))}</td><td class="num">${get('ctl', n) == null ? 'below top 6' : fmt1(get('ctl', n))}</td></tr>`).join('')}</tbody></table></div>`;
    });
  }
  sel.value = 'Zn/Cd/Pb export (zntA, cadA, czcD) + sensors'; sel.onchange = draw; draw();
}

/* ---------- within-host expression of oil degraders and other diesel responders ---------- */
function mtDegraders(sec, J) {
  sec.append(el('h3', {}, 'Inside the oil degraders'));
  const CAT = ['Oxidative stress (kat, sod, ahp, dps, osmC)', 'Copper (cop, cue, pco, cus sensors)', 'Zn/Cd/Pb export (zntA, cadA, czcD) + sensors', 'Mercury (mer)', 'Arsenic (ars, aio)', 'Phage-shock membrane stress (psp)'];
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'mtdeg', 'aria-label': 'Gene set' }); CAT.forEach(k => sel.append(new Option(k, k))); c.append(el('label', { for: 'mtdeg' }, 'Gene set'), sel);
  const box = fig(sec, '', 'Transcripts of the set per million protein-coding transcripts of the same genome (median over libraries with at least 5,000 host reads). This removes the effect of the host’s own growth: a rise here is a response inside the cells. A dot: the genome was too rare in those libraries.');
  const cols = [['2024-25|EXP|pre', '24-25 diesel, before'], ['2024-25|EXP|oiled', '24-25 diesel, after'], ['2024-25|CTL|control', '24-25 control'],
    ['2025-26|EXP|pre', '25-26 diesel, before'], ['2025-26|EXP|oiled', '25-26 diesel, after'], ['2025-26|CTL|pre', '25-26 control, before'], ['2025-26|CTL|control', '25-26 control']];
  function draw() {
    const R = J.degraders.filter(d => d.category === sel.value);
    const mx = d3.max(R.flatMap(d => cols.map(([k]) => d.v[k] || 0))) || 1;
    const col = d3.scaleSequential(d3.interpolateRgb(css('--seq0'), css('--seq4'))).domain([0, Math.log1p(mx)]);
    box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Genome</th><th class="num">Genes</th>${cols.map(([, l]) => `<th class="num">${l}</th>`).join('')}</tr></thead><tbody>${
      R.map(d => `<tr><td>${d.label}</td><td class="num">${d.genes}</td>${cols.map(([k]) => { const v = d.v[k]; return v == null ? '<td class="num" style="color:var(--ink-3)">·</td>'
        : `<td class="num" style="background:${col(Math.log1p(v))};color:${Math.log1p(v) > Math.log1p(mx) * 0.6 ? 'var(--panel)' : 'var(--ink)'}">${fmt(v)}</td>`; }).join('')}</tr>`).join('')}</tbody></table></div>`;
  }
  sel.onchange = draw; draw();
}

/* ---------- mobility vs hydrocarbon genes, and co-selection plasmids ---------- */
function mtMobility(sec, J) {
  sec.append(el('h3', {}, 'Mobility and co-selection'));
  const C = J.cats.filter(c => c.mob).sort((a, b) => b.mob.mobile_assoc - a.mob.mobile_assoc);
  const GC = { 'toxic metal': '--s7', 'micronutrient': '--s3', 'solvent stress': '--s4', 'oxidative stress': '--s8', 'proteotoxic stress': '--other' };
  const box = fig(sec, 'Share of genes on or next to mobile elements', 'Mobile-associated = on a strict geNomad plasmid, phage, provirus, IslandPath island or integron, or within 10 kb of a transposase, integrase or recombinase. Darker segment: the share on a strict plasmid. Dashed line: hydrocarbon genes (confident set), with the same definitions.');
  box.parentNode.insertBefore(el('div', { class: 'legend' }, Object.entries(GC).map(([g, c]) => `<span><i style="background:var(${c});height:10px"></i>${g}</span>`).join('')), box);
  const W = width(box), lw = Math.min(250, W * 0.4), rh = 20, H = C.length * rh + 44;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const x = d3.scaleLinear().domain([0, 25]).range([lw + 10, W - 50]);
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${C.length * rh + 6})`).call(d3.axisBottom(x).ticks(5).tickFormat(v => v + '%'));
  svg.append('line').attr('x1', x(J.hc_ref.mobile_assoc)).attr('x2', x(J.hc_ref.mobile_assoc)).attr('y1', 0).attr('y2', C.length * rh + 6).attr('stroke', css('--exp')).attr('stroke-dasharray', '4 3');
  svg.append('text').attr('x', x(J.hc_ref.mobile_assoc) + 4).attr('y', C.length * rh + 36).text(`oil-degradation genes: ${J.hc_ref.mobile_assoc}%`);
  C.forEach((c, j) => { const y = j * rh + 10, gc = GC[c.group];
    svg.append('text').attr('x', lw).attr('y', y + 4).attr('text-anchor', 'end').text(mtShort(c.name));
    svg.append('rect').attr('x', x(0)).attr('y', y - 6).attr('width', x(c.mob.mobile_assoc) - x(0)).attr('height', 12).attr('rx', 2).attr('fill', css(gc)).attr('opacity', .45)
      .on('mousemove', ev => showTip(ev, `<b>${c.name}</b> (${fmt(c.n_conf)} genes)<br>mobile-associated ${c.mob.mobile_assoc}%<br>strict plasmid ${c.mob.plasmid_strict}%, island ${c.mob.genomic_island}%, integron ${c.mob.integron}%<br>IS/transposase ≤10 kb ${c.mob.IS_within_10kb}%<br>within 10 kb of an oil gene ${c.mob.hc_within_10kb}%`)).on('mouseleave', hideTip);
    svg.append('rect').attr('x', x(0)).attr('y', y - 6).attr('width', x(c.mob.plasmid_strict) - x(0)).attr('height', 12).attr('rx', 2).attr('fill', css(gc)).attr('pointer-events', 'none');
    svg.append('text').attr('x', x(c.mob.mobile_assoc) + 4).attr('y', y + 4).attr('class', 'mono').text(`${d3.format('.0f')(c.mob.mobile_assoc)}%`); });
  const b2 = fig(sec, 'Plasmids that carry both metal resistance and oil-degradation genes', `${J.plasmid_counts.both} of the ${J.plasmid_counts.hc} strict geNomad plasmids with hydrocarbon genes also carry toxic-metal resistance genes (${J.plasmid_counts.metal} strict plasmids carry metal genes in all). The largest may be chromids or megaplasmids. Selection on either cargo would keep both.`);
  b2.innerHTML = `<div class="tablewrap" style="max-height:340px"><table><thead><tr><th>Contig</th><th>Host</th><th class="num">kb</th><th>Metal genes</th><th>Oil genes</th></tr></thead><tbody>${
    J.coplasmids.map(p => `<tr><td class="mono">${p.contig}</td><td>${p.host}</td><td class="num">${d3.format('.0f')(p.len / 1000)}</td><td>${p.metal.map(mtShort).join(', ')}</td><td>${p.hc.join(', ')}</td></tr>`).join('')}</tbody></table></div>`;
}

RENDER.metals = metalsSec;
