/* Eukaryotes section: algae, grazers, fungi and their viruses. Data: euks.json (euks/scripts/euks_site.py).
   Uses helpers from app.js (el, fig, width, css, showTip, hideTip, fmt, fmt1, timecourse, legendPools). */
const euSgn = v => v == null ? '–' : (v > 0 ? '+' : '') + d3.format('.1f')(v);
const euP = v => v == null ? '–' : v >= 10 ? d3.format('.0f')(v) : v >= 1 ? d3.format('.1f')(v) : d3.format('.2f')(v);
const euG = (J, n) => J.groups.find(g => g.name.startsWith(n));
const euS = (J, k) => J.stats.find(s => s.key === k);
function euSeries(J, values) {
  const out = {};
  J.samples.forEach((s, i) => { if (!['CTL', 'EXP'].includes(s.pool) || !['SW', 'SI'].includes(s.medium) || s.day == null || values[i] == null) return;
    const k = `${s.season}|${s.pool}|${s.medium}`; ((out[k] ||= {})[s.day] ||= []).push([values[i], s.name]); });
  return Object.entries(out).map(([k, byDay]) => { const [season, pool, medium] = k.split('|');
    return { season, pool, medium, pts: Object.entries(byDay).map(([d, a]) => ({ day: +d, v: d3.mean(a, x => x[0]), name: a.map(x => x[1]).join(', ') })) }; });
}

function euksSec(sec) {
  const J = D.euks;
  if (!J) { sec.innerHTML = '<h2>Eukaryotes</h2><p class="kicker">euks.json is missing.</p>'; return; }
  const pm = euS(J, 'plastid_mrna_pct'), gv = euS(J, 'giant_virus_cpm'), hap = euG(J, 'Haptophytes'), chr = euG(J, 'Chrysophytes');
  const vp = J.virus_links.find(l => l.virus === 'virophage_cpm' && /Phaeocystis/.test(l.taxon));
  sec.innerHTML = `<h2>Eukaryotes</h2>
  <p class="kicker">Microalgae, protist grazers, fungi and small animals share the pools with the bacteria. Diesel hits them differently. In the 2025-26 diesel pool the algae collapse: plastid transcription falls about tenfold, and haptophytes, dinoflagellates, cryptophytes and diatoms lose ground. Bacteria-eating flagellates take their place, above all chrysophytes and choanoflagellates. The giant viruses that infect algae disappear with their hosts. Composition comes from 18S rRNA reads (SILVA) as a share of all eukaryotic 18S, so it is not diluted by the bacterial bloom. Read counts per cell differ widely between groups, so compare each group with itself across samples, not groups with each other. The Kaiju database used for mRNA holds no eukaryotes, so eukaryotic activity is tracked through plastid-encoded transcripts (which Kaiju calls “cyanobacterial”) and plastid genomes. Mammal, fish and insect 18S (likely contamination) is removed. One pool per treatment, so all contrasts are descriptive. The 2024-25 pools were filled about two weeks apart because the intake lines froze, so pre-dosing differences that winter reflect pool history.</p>`;
  const st = el('div', { class: 'stat-row' }); sec.append(st);
  [[`${euP(pm.oil25)}% vs ${euP(pm.ctl25)}%`, `of mRNA is plastid-encoded, diesel vs control seawater 2025-26 (BACI ${euSgn(pm.baci25)} log₂)`],
   [`${euP(chr.oil25)}% vs ${euP(chr.ctl25)}%`, 'of eukaryotic 18S from chrysophytes and synurophytes, diesel vs control 2025-26; many are bacteria-eating mixotrophs'],
   [`${euSgn(hap.lfc25)}`, `log₂ diesel ÷ control, haptophytes 2025-26 (BACI ${euSgn(hap.baci25)})`],
   [`ρ ${vp ? d3.format('.2f')(vp.rho) : '–'}`, 'between the control-pool virophage bloom and the haptophyte Phaeocystis across 2025-26 samples']]
    .forEach(([v, l]) => st.append(el('div', { class: 'stat' }, `<div class="v">${v}</div><div class="l">${l}</div>`)));
  euTable(sec, J);
  euExplorer(sec, J);
  euPlastids(sec, J);
  euNuclear(sec, J);
  euBloom(sec, J);
  euViruses(sec, J);
  euIce(sec, J);
  euReading(sec, J);
}

const euCol = v => { if (v == null) return 'transparent'; const t = Math.max(-1, Math.min(1, v / 4)); return t >= 0 ? d3.interpolateRgb(css('--panel'), css('--exp'))(t * 0.85) : d3.interpolateRgb(css('--panel'), css('--ctl'))(-t * 0.85); };
const euCell = v => `<td class="num" style="background:${euCol(v)}">${euSgn(v)}</td>`;

function euTable(sec, J) {
  sec.append(el('h3', {}, 'Who is there, and how each group responds'));
  const box = fig(sec, '', 'Share of eukaryotic 18S reads (%), seawater medians after dosing. Ratios are log₂ (diesel ÷ control). BACI = change from each pool’s own pre-dosing samples, diesel minus control; 2024-25 has a single pre-dosing sample per pool. Ice = paired ice minus water (23 pairs); Top = top minus bottom ice (13 pairs). Orange = higher in diesel, ice or top ice; blue = higher in control, water or bottom ice. Click a row for its time course.');
  const roles = { phototroph: 'phototroph', mixotroph: 'mixotroph', grazer: 'grazer', parasite: 'parasite', osmotroph: 'osmotroph (absorbs food)', animal: 'animal', other: '' };
  box.innerHTML = `<div class="tablewrap" style="max-height:none"><table><thead><tr><th>Group</th><th class="num">% diesel 24-25</th><th class="num">% control 24-25</th><th class="num">% diesel 25-26</th><th class="num">% control 25-26</th><th class="num">24-25</th><th class="num">25-26</th><th class="num">BACI 24-25</th><th class="num">BACI 25-26</th><th class="num">Ice</th><th class="num">Top</th></tr></thead><tbody>${
    J.groups.map(g => `<tr data-c="g|${g.name}"><td style="white-space:normal;min-width:240px">${g.name}<div style="font-size:11.5px;color:var(--ink-3)">${roles[g.role] || ''}</div></td><td class="num">${euP(g.oil24)}</td><td class="num">${euP(g.ctl24)}</td><td class="num">${euP(g.oil25)}</td><td class="num">${euP(g.ctl25)}</td>${euCell(g.lfc24)}${euCell(g.lfc25)}${euCell(g.baci24)}${euCell(g.baci25)}${euCell(g.ice)}${euCell(g.top)}</tr>`).join('')}</tbody></table></div>`;
  box.querySelectorAll('tr[data-c]').forEach(tr => { tr.onclick = () => { const s = document.getElementById('eusel'); if (s) { s.value = tr.dataset.c; s.onchange(); s.scrollIntoView({ behavior: 'smooth', block: 'center' }); } }; });
  // genera
  const b2 = fig(sec, 'Main genera', 'Genera (or the deepest SILVA rank) averaging at least 0.3% of eukaryotic 18S. Same columns as above; click for a time course. SILVA places a few sequences oddly (for example “Prorocentrum minimum” inside Chrysophyceae); the group follows the lineage path.');
  const G = [...J.genera].sort((a, b) => b.mean - a.mean);
  b2.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Genus · group</th><th class="num">Mean %</th><th class="num">% diesel 25-26</th><th class="num">% control 25-26</th><th class="num">24-25</th><th class="num">25-26</th><th class="num">BACI 25-26</th><th class="num">Ice</th><th class="num">Top</th></tr></thead><tbody>${
    G.map(g => `<tr data-c="n|${g.name}"><td style="white-space:normal;min-width:260px">${g.name}<div style="font-size:11.5px;color:var(--ink-3)">${g.group}</div></td><td class="num">${euP(g.mean)}</td><td class="num">${euP(g.oil25)}</td><td class="num">${euP(g.ctl25)}</td>${euCell(g.lfc24)}${euCell(g.lfc25)}${euCell(g.baci25)}${euCell(g.ice)}${euCell(g.top)}</tr>`).join('')}</tbody></table></div>`;
  b2.querySelectorAll('tr[data-c]').forEach(tr => { tr.onclick = () => { const s = document.getElementById('eusel'); if (s) { s.value = tr.dataset.c; s.onchange(); s.scrollIntoView({ behavior: 'smooth', block: 'center' }); } }; });
}

function euExplorer(sec, J) {
  sec.append(el('h3', {}, 'Time courses'));
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'eusel', 'aria-label': 'Series' });
  const add = (label, items) => { const og = document.createElement('optgroup'); og.label = label; items.forEach(([v, t]) => og.append(new Option(t, v))); sel.append(og); };
  add('Activity and viruses', J.stats.map(s => [`s|${s.key}`, s.name]));
  add('Groups (% of eukaryotic 18S)', J.groups.map(g => [`g|${g.name}`, g.name]));
  add('Genera (% of eukaryotic 18S)', [...J.genera].sort((a, b) => a.name.localeCompare(b.name)).map(g => [`n|${g.name}`, g.name]));
  c.append(el('label', { for: 'eusel' }, 'Series'), sel);
  const v = el('div'); sec.append(v);
  function draw() {
    v.innerHTML = ''; const [k, n] = sel.value.split(/\|(.+)/);
    const vals = k === 's' ? J.stat_tc[n] : k === 'g' ? J.group_tc[n] : J.genus_tc[n];
    const title = k === 's' ? J.stats.find(s => s.key === n).name : `${n}: % of eukaryotic 18S`;
    const box = fig(v, title, 'Lines join per-day means of replicate libraries.');
    legendPools(box.parentNode.insertBefore(el('div'), box));
    timecourse(box, euSeries(J, vals), { ylab: k === 's' ? '' : '%', unit: k === 's' ? '' : '%', h: 200, floor: 0.01 });
  }
  sel.onchange = draw; sel.value = 's|plastid_mrna_pct'; draw();
}

function euPlastids(sec, J) {
  sec.append(el('h3', {}, 'Photosynthesis switches off in the diesel pool'));
  const pm = euS(J, 'plastid_mrna_pct'), pt = euS(J, 'plastid_mag_tpm');
  sec.append(el('p', { class: 'kicker' }, `Algal chloroplasts carry their own small genomes, and several were assembled as MAGs. Their transcripts (psbA, rbcL and the rest of the plastid genes) measure photosynthetic activity directly. In 2025-26 both pools started level: ${euP(pm.pre_exp25)}% vs ${euP(pm.pre_ctl25)}% of mRNA was plastid-encoded before dosing. After dosing it was ${euP(pm.oil25)}% in the diesel pool against ${euP(pm.ctl25)}% in the control (BACI ${euSgn(pm.baci25)} log₂); plastid-genome TPM tells the same story (BACI ${euSgn(pt.baci25)}). The Cycloclasticus bloom dilutes every other transcript, but by at most about twofold, so most of this drop is real. In 2024-25 the diesel pool fell to ${euP(pm.oil24)}% against ${euP(pm.ctl24)}% in the control, from similar starting points (${euP(pm.pre_exp24)}% and ${euP(pm.pre_ctl24)}%). Most plastid genomes fall in the diesel pool. The one that ends higher there than in the control, lorbin_664_sub, was already more than three times higher in the diesel pool before dosing, and it fell too.`));
  const box = fig(sec, 'Plastid genomes (seawater median TPM)', 'Host = SILVA best hit of the plastid 16S where one was recovered. “Unidentified alga” plastids lack a recovered 16S.');
  box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Plastid MAG</th><th class="num">pre 25-26<br>diesel / ctl</th><th class="num">after 25-26<br>diesel / ctl</th><th class="num">after 24-25<br>diesel / ctl</th></tr></thead><tbody>${
    J.plastids.map(p => `<tr><td>${p.label} <span class="mono" style="color:var(--ink-3)">${p.bin}</span></td><td class="num">${fmt(p.pre_exp25)} / ${fmt(p.pre_ctl25)}</td><td class="num"><span style="color:${p.oil25 < p.ctl25 / 2 ? css('--ctl') : 'inherit'}">${fmt(p.oil25)}</span> / ${fmt(p.ctl25)}</td><td class="num">${fmt(p.oil24)} / ${fmt(p.ctl24)}</td></tr>`).join('')}</tbody></table></div>`;
}

function euNuclear(sec, J) {
  const E = J.euk_bins || []; if (!E.length) return;
  sec.append(el('h3', {}, 'Eukaryotic nuclear genomes among the MAGs'));
  sec.append(el('p', { class: 'kicker' }, `A rule applied to every bin flags ${E.length} as eukaryotic: Tiara classifies at least half of each bin’s assembled length as eukaryotic (most of them 80–98%). GTDB-Tk, which only knows bacteria and archaea, had placed them as bacteria or Asgard archaea, so they are relabelled “Eukaryote (Tiara)” across the site and left out of the analyses that assume a bacterial genome (growth index, oxidase and nitrogen-hunger shares, AAP and rhodopsin hosts). Their taxonomy within the eukaryotes is not resolved yet.`));
  const box = fig(sec, 'Eukaryotic bins (seawater median TPM)', 'TPM = share of community transcripts mapped to the bin. Tiara = share of the bin’s length classified as eukaryotic.');
  box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Bin</th><th class="num">Size</th><th class="num">Tiara eukaryote</th><th class="num">pre 25-26<br>diesel / ctl</th><th class="num">after 25-26<br>diesel / ctl</th><th class="num">after 24-25<br>diesel / ctl</th></tr></thead><tbody>${
    E.map(p => `<tr><td class="mono">${p.bin}</td><td class="num">${d3.format('.1f')(p.mbp)} Mbp, ${fmt(p.contigs)} contigs</td><td class="num">${d3.format('.0%')(p.tiara_euk)}</td><td class="num">${fmt(p.pre_exp25)} / ${fmt(p.pre_ctl25)}</td><td class="num">${fmt(p.oil25)} / ${fmt(p.ctl25)}</td><td class="num">${fmt(p.oil24)} / ${fmt(p.ctl24)}</td></tr>`).join('')}</tbody></table></div>`;
  const c = el('div', { class: 'controls' }); sec.append(c);
  const sel = el('select', { id: 'eunuc', 'aria-label': 'Eukaryotic bin' }); E.forEach(p => sel.append(new Option(p.bin, p.bin))); c.append(el('label', { for: 'eunuc' }, 'Bin'), sel);
  const v = el('div'); sec.append(v);
  function draw() { v.innerHTML = ''; const p = E.find(x => x.bin === sel.value);
    const b = fig(v, `${p.bin}: transcripts per million`, 'Share of community transcripts; lines join per-day means of replicate libraries.');
    legendPools(b.parentNode.insertBefore(el('div'), b)); timecourse(b, euSeries(J, p.tpm), { ylab: 'TPM', unit: 'TPM', h: 180, floor: 1 }); }
  sel.onchange = draw; draw();
}

function euBloom(sec, J) {
  sec.append(el('h3', {}, 'Grazers follow the bacterial bloom'));
  sec.append(el('p', { class: 'kicker' }, 'Rank correlation between each group’s 18S share and Cycloclasticus transcripts (TPM) across the 2025-26 diesel-pool samples, water and ice. Small bacteria-eating flagellates rise with the bloom: telonemids, bicosoecids, MAST stramenopiles, cercozoans and chrysophytes. Phototrophs and the larger, algae-eating katablepharids and dinoflagellates fall. Both trends also follow time since dosing, so this is a pattern, not proof of feeding. Choanoflagellates and ciliates, the other classic bacterivores, do not track the bloom.'));
  const R = J.bloom_links.filter(d => d.rho != null);
  const box = fig(sec, 'Correlation with the Cycloclasticus bloom (2025-26 diesel pool)', 'Spearman ρ; bars coloured by trophic role. Hover for p.');
  const W = width(box), lw = Math.min(320, W * 0.42), rh = 20, H = R.length * rh + 34;
  const svg = d3.select(box).append('svg').attr('width', W).attr('height', H);
  const x = d3.scaleLinear().domain([-1, 1]).range([lw + 10, W - 16]);
  const col = { phototroph: '--s3', mixotroph: '--s7', grazer: '--exp', parasite: '--s8', osmotroph: '--s4', animal: '--s4', other: '--grid' };
  svg.append('line').attr('x1', x(0)).attr('x2', x(0)).attr('y1', 0).attr('y2', R.length * rh + 4).attr('stroke', css('--rule'));
  R.forEach((d, i) => { const y = i * rh + 4;
    svg.append('text').attr('x', lw).attr('y', y + 12).attr('text-anchor', 'end').style('font-size', '12px').text(d.group);
    svg.append('rect').attr('x', Math.min(x(0), x(d.rho))).attr('y', y + 2).attr('width', Math.abs(x(d.rho) - x(0))).attr('height', 13).attr('rx', 2).attr('fill', css(col[d.role] || '--s4'))
      .on('mousemove', ev => showTip(ev, `<b>${d.group}</b> (${d.role})<br>ρ ${d3.format('.2f')(d.rho)}, p ${d3.format('.2g')(d.p)}`)).on('mouseleave', hideTip); });
  svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${R.length * rh + 6})`).call(d3.axisBottom(x).ticks(5));
  box.parentNode.insertBefore(el('div', { class: 'legend' }, [['--s3', 'phototroph'], ['--s7', 'mixotroph'], ['--exp', 'grazer'], ['--s8', 'parasite'], ['--s4', 'osmotroph / animal']]
    .map(([c, l]) => `<span><i style="background:var(${c});height:10px"></i>${l}</span>`).join('')), box);
}

function euViruses(sec, J) {
  const gv = euS(J, 'giant_virus_cpm'), vp = euS(J, 'virophage_cpm'), hap = euG(J, 'Haptophytes');
  sec.append(el('h3', {}, 'Giant viruses and the virophage track their likely hosts'));
  sec.append(el('p', { class: 'kicker' }, `Giant viruses (Megaviricetes) infect algae and protists. In 2025-26 their transcripts fell from ${fmt(gv.pre_exp25)} CPM before dosing to ${fmt(gv.oil25)} in the diesel pool, while the control kept ${fmt(gv.ctl25)} (BACI ${euSgn(gv.baci25)}). Across samples they rise and fall with Mamiellophyceae (Micromonas, Bathycoccus), whose prasinoviruses are classic giant viruses, and with Thalassiosira and katablepharids. So the giant-virus collapse on the <a href="#viruses">Viruses</a> page most likely follows the loss of their hosts, not a direct effect on the viruses. The virophage contig_0049668, which bloomed only in the 2025-26 control pool (${fmt(vp.ctl25)} CPM after dosing, ${euP(vp.oil25)} in diesel), tracks the haptophyte Phaeocystis and the pedinellid Dictyochophyceae most closely. Haptophytes were almost absent from the diesel pool after dosing (${euP(hap.oil25)}% vs ${euP(hap.ctl25)}% of eukaryotic 18S). Phaeocystis giant viruses are known to carry virophage-like partners, which makes a Phaeocystis–giant virus–virophage chain in the control pool the leading hypothesis. Bicosoecids such as Cafeteria, the host in the classic virophage system, are negatively correlated (ρ −0.45).`));
  const vl = J.virus_links;
  const box = fig(sec, 'Strongest positive correlations (Spearman, all pool samples of one winter, q < 0.05)', 'Each eukaryotic genus or group vs giant-virus or virophage transcripts. Correlation across samples is not infection; co-occurrence is shaped by shared responses to diesel and season.');
  const lab = { giant_virus_cpm: 'Giant viruses', virophage_cpm: 'Virophage contig_0049668' };
  box.innerHTML = `<div class="tablewrap"><table><thead><tr><th>Virus</th><th>Winter</th><th>Eukaryote</th><th>Group</th><th class="num">ρ</th></tr></thead><tbody>${
    vl.map(l => `<tr><td>${lab[l.virus]}</td><td>${l.season}</td><td>${l.taxon}</td><td>${l.level === 'group' ? '(whole group)' : l.group}</td><td class="num">${d3.format('.2f')(l.rho)}</td></tr>`).join('')}</tbody></table></div>`;
}

function euIce(sec, J) {
  const f = euG(J, 'Fungi'), gr = euG(J, 'Green'), di = euG(J, 'Dinoflag'), ki = euG(J, 'Kinetoplastids');
  const cr = J.genera.find(g => /Cryothecomonas/.test(g.name)), pr = J.genera.find(g => /Protaspis/.test(g.name));
  sec.append(el('h3', {}, 'Ice and water'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li><b>Fungi</b> make up a larger share in ice than in the paired water in ${f.ice_up} of ${f.ice_n} pairs (${euSgn(f.ice)} log₂). They are mostly ascomycete yeasts and moulds (Helotiales, Capnodiales), which are known from sea ice and brine.</li>
    <li><b>Cercozoan parasitoids and grazers of diatoms</b> are enriched in ice: Cryothecomonas in ${cr ? cr.ice_up : '–'} of ${cr ? cr.ice_n : '–'} pairs and Protaspis in ${pr ? pr.ice_up : '–'} of ${pr ? pr.ice_n : '–'}. Cryothecomonas is a classic sea-ice parasitoid of diatoms.</li>
    <li><b>Kinetoplastids and euglenids</b> are lower in ice (higher in only ${ki.ice_up} of ${ki.ice_n} pairs).</li>
    <li><b>Top vs bottom:</b> green algae are higher in top ice in ${gr.top_up} of ${gr.top_n} pairs; dinoflagellates are lower in top ice in all ${di.top_n - di.top_up} of ${di.top_n}.</li>`));
}

function euReading(sec, J) {
  const ch = euG(J, 'Chrysophytes'), di = euG(J, 'Diatoms'), dn = euG(J, 'Dinoflag'), cr = euG(J, 'Cryptophytes'), ce = euG(J, 'Cercozoa'), ka = euG(J, 'Katablepharids');
  sec.append(el('h3', {}, 'How the pieces fit, and what is still a hypothesis'));
  sec.append(el('ul', { class: 'kicker' }, `
    <li><b>Observed: diesel shuts down the algae in 2025-26.</b> Plastid transcription and most plastid genomes fall about tenfold against each pool’s own baseline. Haptophytes, dinoflagellates (${euSgn(dn.baci25)}), cryptophytes (${euSgn(cr.baci25)}) and diatoms (${euSgn(di.baci25)}) lose share of the eukaryotes. Diesel is toxic to microalgae, and oil films and droplets also shade the water.</li>
    <li><b>Observed: bacteria-eating flagellates replace them.</b> Chrysophytes and synurophytes (Synura, Mallomonas, Paraphysomonas, Chrysolepidomonas) rise to ${euP(ch.oil25)}% of eukaryotic 18S in the diesel pool against ${euP(ch.ctl25)}% in the control (BACI ${euSgn(ch.baci25)}). Many of these are mixotrophs or obligate phagotrophs (Paraphysomonas). Choanoflagellates also gain (${euSgn(euG(J, 'Choanoflag').baci25)} BACI). Within the diesel pool, telonemids, bicosoecids and MAST rise in step with the Cycloclasticus bloom. Feeding on the bacterial bloom is the simplest explanation.</li>
    <li><b>Observed: 2024-25 differs.</b> There, Cercozoa (${euSgn(ce.lfc24)} log₂) and katablepharids (${euSgn(ka.lfc24)}) are higher in the diesel pool and chrysophytes are not, and the plastid decline is milder. That winter had a single pre-dosing sample per pool and pools filled two weeks apart, so the baseline is weak.</li>
    <li><b>Hypothesis: viruses follow hosts.</b> Giant viruses disappear with their algal hosts, and the control-pool virophage tracks Phaeocystis. Resolving it would need host-specific markers: Phaeocystis virus major capsid genes, or virophage integration sites in a Phaeocystis genome.</li>
    <li><b>Limits.</b> 18S read share reflects gene copy number as well as cells (dinoflagellates and ciliates carry many copies), and it shows composition, not total biomass. Chlorophyll and cell counts from the pools would turn these shifts into absolute numbers.</li>`));
}

RENDER.euks = euksSec;
