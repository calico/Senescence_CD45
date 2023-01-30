all_immune_markers = ['Adgre1', 'Aif1', 'Arg1', 'Batf3', 'Bst2', 'Ccr2', 'Ccr7', 'Cd14',
       'Cd163', 'Cd19', 'Cd27', 'Cd274', 'Cd28', 'Cd3e', 'Cd4', 'Cd68',
       'Cd69', 'Cd80', 'Cd86', 'Cd8a', 'Clec9a', 'Csf1r', 'Csf3r',
       'Cx3cr1', 'Cxcr2', 'Cxcr4', 'Cxcr6', 'Dpp4', 'Eomes', 'Fcer1a',
       'Fcer1g', 'Fcgr3', 'Foxp3', 'Gata3', 'Havcr2', 'Ifit3', 'Ifit3b',
       'Igha', 'Ighd', 'Ighm', 'Igkc', 'Il2rb', 'Il3ra', 'Il7r', 'Itga1',
       'Itgae', 'Itgam', 'Itgax', 'Jchain', 'Kit', 'Klf2', 'Klrb1c',
       'Klrd1', 'Klrk1', 'Lag3', 'Ly6c1', 'Ly6g', 'Mafb', 'Mki67', 'Mrc1',
       'Ms4a1', 'Ncr1', 'Nkg7', 'Pdcd1', 'Ptprc', 'Rora', 'Rorc',
       'S100a8', 'Sdc1', 'Sell', 'Siglecf', 'Siglech', 'Tbx21', 'Thy1',
       'Tnfrsf17', 'Tnfrsf9', 'Tox', 'Tpsab1', 'Tpsb2']

tcell_markers = ['Cd3e', 'Cd4', 'Cd8a', 'Cd69', 'Cd44', 'Mki67', 'Ctla4', 'Ifng',
       'Lag3', 'Havcr2', 'Tox', 'Pdcd1', 'Tnfrsf9', 'Cxcl9', 'Il7r',
       'Ncr1', 'Eomes', 'Cxcr3', 'Klrb1c', 'Itga1', 'Tnfsf10', 'Il2ra',
       'Cxcr6', 'Id2', 'Rorc', 'Tcf7', 'Tigit', 'Nt5e', 'Vsir', 'Cd276',
       'Cd160', 'Cd52', 'Ncam1', 'Tnf', 'Prf1', 'Lamp1', 'Fasl', 'Gzmb',
       'Ccl5', 'Tnfrsf9', 'Tcf7', 'Cd44', 'Cd69', 'Ifng', 'Tnf', 'Prf1']

tcell_dysfunction_markers = ['Pdcd1', 'Havcr2', 'Tox', 'Lag3', 'Ctla4', 'Tigit', 'Btla',
       'Cd160', 'Ido1', 'Slamf6', 'Nt5e', 'Vsir', 'Cd276']

tcell_and_nkt_markers = ['Klf2', 'Sell', 'Klrb1c', 'Cd4', 'Cd69', 'Cd44', 'Foxp3', 'Il2ra',
       'Il7r', 'Cd8a', 'Cd27', 'Cd28', 'Ccr7', 'Lag3', 'Ifng', 'Havcr2',
       'Tox', 'Pdcd1', 'Eomes', 'Tbx21', 'Il7r', 'Ncr1', 'Klrb1c', 'Cd3e',
       'Itga1', 'Tnfsf10', 'Il2ra', 'Cxcr6', 'Id2', 'Rorc']

all_myeloid_markers = ['Adgre1', 'Aif1', 'Arg1', 'Ccl5', 'Ccr2', 'Ccr3', 'Cd163',
       'Cd200r3', 'Cd274', 'Cd38', 'Cd63', 'Cd68', 'Cd80', 'Cd86',
       'Chil3', 'Ciita', 'Csf1r', 'Cx3cr1', 'Egr2', 'Enpp3', 'Fcer1a',
       'Fcgr1', 'Fcgr3', 'Gsk3b', 'Il1b', 'Itga2', 'Itgam', 'Itgax',
       'Lilra5', 'Ly6c1', 'Ly6g', 'Mafb', 'Med16', 'Mrc1', 'Siglecf',
       'Socs3', 'Stat1', 'Trem2']

# MGI: GO:0019882
antigen_presentation = ['Azgp1', 'Ccl21a', 'Ccr7', 'Cd68', 'Cd74', 'Clec4b2', 'Ctse',
       'Ctss', 'Fcer1g', 'Fcgr2b', 'Fcgr3', 'Fcgr4', 'Fgl2', 'Flt3',
       'H2-Aa', 'H2-Ab1', 'H2-DMa', 'H2-DMb1', 'H2-DMb2', 'H2-Eb1',
       'H2-Eb2', 'H2-M2', 'H2-M5', 'H2-Ob', 'H2-Q10', 'H2-Q4', 'H2-Q6',
       'H2-Q7', 'H2-T24', 'Icam1', 'Ifi30', 'Ifng', 'Ighm', 'Nod2',
       'Psap', 'Ptpn22', 'Rab27a', 'Rab8b', 'Relb', 'Slc11a1', 'Thbs1',
       'Treml4', 'Unc93b1', 'Wdfy4']

# Source: MSigDB ifng hallmark set converted to mouse names via
# `map_human_to_mouse_orthologs` (see utils.py)
ifng_hallmark_set = ['Bank1', 'Batf2', 'C1s1', 'C1s2', 'Ccl2', 'Cd274', 'Cd38', 'Cd40',
       'Cd69', 'Cd74', 'Cd86', 'Cdkn1a', 'Ciita', 'Cmpk2', 'Csf2rb2',
       'Csf2rb', 'Cxcl10', 'Cxcl9', 'Epsti1', 'Fas', 'Fgl2', 'Fpr1',
       'Gbp7', 'Gbp8', 'Gbp10', 'Gbp4', 'Gbp9', 'Gch1', 'Gpr18', 'Gzma',
       'Herc6', 'H2-Q4', 'H2-Q7', 'H2-Q10', 'H2-M5', 'H2-Q6', 'H2-M2',
       'H2-Q4', 'H2-Q7', 'H2-Q10', 'H2-M5', 'H2-Q6', 'H2-M2', 'H2-DMa',
       'H2-Eb2', 'H2-Q4', 'H2-Q7', 'H2-Q10', 'H2-M5', 'H2-Q6', 'H2-M2',
       'Icam1', 'Ifi30', 'Ifih1', 'Ifit2', 'Ifit3', 'Ifit3b', 'Ifitm3',
       'Ifitm1', 'Ifitm3', 'Ifitm1', 'Il10ra', 'Il15', 'Il15ra', 'Il2rb',
       'Il6', 'Irf1', 'Irf4', 'Irf5', 'Irf7', 'Irf8', 'Isg15', 'Itgb7',
       'Jak2', 'Klrk1', 'Lcp2', 'Mt1', 'Mx1', 'Nampt', 'Nfkb1', 'Nfkbia',
       'Nlrc5', 'Oas2', 'P2ry14', 'Parp14', 'Pde4b', 'Peli1', 'Pfkp',
       'Pim1', 'Ptgs2', 'Ptpn6', 'Rapgef6', 'Ripk2', 'Rnf213', 'Rsad2',
       'Samd9l', 'Samhd1', 'Selp', 'Serping1', 'Slamf7', 'Slc25a28',
       'Socs1', 'Socs3', 'Sp110', 'St3gal5', 'St8sia4', 'Stat1', 'Stat2',
       'Stat4', 'Tnfaip2', 'Tnfaip3', 'Tnfaip6', 'Tnfsf10', 'Trim25',
       'Txnip', 'Upp1', 'Vcam1', 'Xcl1', 'Casp4', 'Ccl5', 'Ccl7', 'Cfh',
       'Ifit1', 'March1']

flow_markers = ["Itgam", "Itgax", "Ly6c1", "Ly6g", "Cd86", "Cd80", "Itgae", "Mrc1",
    "Adgre1", "H2-Ab1"]