"""Cold-adaptation gene categories: KOs (corrected genemap KO), Pfam (hmmsearch --cut_ga on BAKTA_EXTRA proteins) and product regexes."""
CATS = {
 'Ice-binding proteins': dict(pf=['DUF3494', 'AFP', 'CfAFP', 'Ice_nucleation']),
 'Cold-shock proteins': dict(ko=['K03704'], pf=['CSD'], maxlen=120),
 'Cold DEAD-box helicases': dict(ko=['K05592', 'K11927', 'K05590', 'K03732'], pf=['DEADboxA']),
 'Cold-acclimation translation/RNA': dict(ko=['K00962', 'K02834', 'K02600', 'K02518', 'K02519', 'K12573', 'K03545']),
 'Chaperones': dict(ko=['K04043', 'K04077', 'K04078', 'K04079', 'K03695', 'K04080', 'K03686', 'K03687']),
 'Compatible solute synthesis': dict(ko=['K00108', 'K00130', 'K06718', 'K00836', 'K06720', 'K10674', 'K00697', 'K01087',
                                         'K05343', 'K06044', 'K01236', 'K00147', 'K00931', 'K00286'], pf=['Ectoine_synth', 'Osmo_MPGsynth']),
 'Compatible solute transport': dict(ko=['K02168', 'K05020', 'K02000', 'K02001', 'K02002', 'K05845', 'K05846', 'K05847'], pf=['BCCT', 'OpuAC']),
 'DMSP catabolism': dict(ko=['K17486', 'K16953']),
 'Fatty-acid desaturases': dict(ko=['K10255', 'K00507', 'K03921'], pf=['FA_desaturase', 'FA_desaturase_2', 'Lipid_DES']),
 'Membrane FA elongation (fabF/A/B)': dict(ko=['K09458', 'K01716', 'K00647']),
 'PUFA synthase (pfa)': dict(prod=r'polyunsaturated fatty acid synthase|omega-3 polyunsaturated|\bpfa[abcde]\b'),
 'Carotenoid synthesis': dict(ko=['K02291', 'K10027', 'K06443', 'K02294', 'K09836'], pf=['Caroten_synth']),
 'Oxidative stress': dict(ko=['K03781', 'K03782', 'K04564', 'K04565', 'K03386', 'K03387', 'K04063']),
 'EPS export': dict(ko=['K01991'], prod=r'exopolysaccharide'),
 'Osmotic channels/K+ uptake': dict(ko=['K03442', 'K03282', 'K01546', 'K01547', 'K01548']),
}
# Index order of amino-acid count columns (aacount.py)
AA = 'ACDEFGHIKLMNPQRSTVWY'
