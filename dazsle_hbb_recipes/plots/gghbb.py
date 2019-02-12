""" dazsle gghbb analysis plots """

from fnal_column_analysis_tools import hist
from fnal_column_analysis_tools.hist import plot

from copy import deepcopy

dataset = hist.Cat("dataset", "Primary dataset")
systematic = hist.Cat("systematic", "Systematic Variation")
gencat = hist.Bin("ak8_isHadronicV", "Matched", [-1,0,1,2,3,9,10,11])
# one can relabel intervals, although process mapping obviates this
titles = ["Data","QCD", "V(light) matched", "V(c) matched", "V(b) matched", "Top W(ud)+b", "Top W(cs)+b"]
for i,v in enumerate(gencat.identifiers()):
    setattr(v, 'label', titles[i])

jetpt = hist.Bin("ak8_pt", "Jet $p_T$", [450, 500, 550, 600, 675, 800, 1000])
jetpt_coarse = hist.Bin("ak8_pt", "Jet $p_T$", [450, 800])
jetmass = hist.Bin("ak8_msd", "Jet $m_{sd}$", 23, 40, 201)
jetmass_coarse = hist.Bin("ak8_msd", "Jet $m_{sd}$", [40, 100, 140, 200])
jetrho = hist.Bin("jetrho", r"Jet $\rho$", 13, -6, -2.1)
doubleb = hist.Bin("ak8_deepdoubleb", "Double-b", 20, 0., 1)
doublec = hist.Bin("ak8_deepdoublec", "Double-c", 20, 0., 1.)
doublecvb = hist.Bin("ak8_deepdoublecvb", "Double-cvb", 20, 0., 1.)
doubleb_coarse = [1., 0.93, 0.92, 0.89, 0.85, 0.7]
doubleb_coarse = hist.Bin("ak8_deepdoubleb", "Double-b", doubleb_coarse[::-1])
doublec_coarse = [0.87, 0.84, 0.83, 0.79, 0.69, 0.58]
doublec_coarse = hist.Bin("ak8_deepdoublec", "Double-c", doublec_coarse[::-1])
doublecvb_coarse = [0.93, 0.91, 0.86, 0.76, 0.6, 0.17, 0.12]
doublecvb_coarse = hist.Bin("ak8_deepdoublecvb", "Double-cvb", doublecvb_coarse[::-1])
n2ddt_coarse = hist.Bin("ak8_N2sdb1_ddt", "N2 DDT", [0.])


hists = {}
hists['sumw'] = hist.Hist("sumw", dataset, systematic, hist.Bin("sumw", "Weight value", [0.]))
hists['hjetpt'] = hist.Hist("Events", dataset, gencat, systematic, hist.Bin("ak8_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
hists['hjetpt_sr'] = hist.Hist("Events", dataset, gencat, systematic, hist.Bin("ak8_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
#hists['htagtensor'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, n2ddt_coarse, jetmass_coarse, doubleb, doublec, doublecvb, dtype='f')
hists['hsculpt'] = hist.Hist("Events", dataset, gencat, systematic,
                             jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')
hists['hsculpt_sr'] = hist.Hist("Events", dataset, gencat, systematic,
                                jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')

hists['pfmet_nminus1_sr'] = hist.Hist("Events", dataset, gencat, systematic,
                                      jetpt_coarse, jetmass_coarse, hist.Bin("pfmet", r"PF $p_{T}^{miss}$", 40, 0, 200))
hists['opposite_ak8_n3sdb1_sr'] = hist.Hist("Events", dataset, gencat, systematic,
                                            jetpt_coarse, jetmass_coarse,
                                            hist.Bin("opposite_ak8_n3sdb1", r"Jet $N_{3,sd}^{\beta=1}$", 40, 0.5, 3))
hists['opposite_ak8_tau32_sr'] = hist.Hist("Events", dataset, gencat, systematic,
                                           jetpt_coarse, jetmass_coarse,
                                           hist.Bin("opposite_ak8_tau32", r"Jet $\tau_{32}$", 40, 0, 1))
hists['opposite_ak8_msd_sr'] = hist.Hist("Events", dataset, gencat, systematic,
                                         jetpt_coarse, jetmass_coarse,
                                         hist.Bin("opposite_ak8_msd", r"Jet $\m_{sd}$", 40, 50, 200))
hists['opposite_ak4_leadingDeepCSV_sr'] = hist.Hist("Events", dataset, gencat, systematic,
                                                    jetpt_coarse, jetmass_coarse,
                                                    hist.Bin("opposite_ak4_leadingDeepCSV", "Max(DeepCSV) (of $\leq4$ leading)", 40, 0, 1))

#muon control region
hists['hjetpt_mucr'] = hist.Hist("Events", dataset, gencat, systematic, hist.Bin("ak8_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
hists['hsculpt_mucr'] = hist.Hist("Events", dataset, gencat, systematic,
                                jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')
hists['opposite_ak8_n3sdb1_mucr'] = hist.Hist("Events", dataset, gencat, systematic,
                                            jetpt_coarse, jetmass_coarse,
                                            hist.Bin("opposite_ak8_n3sdb1", r"Jet $N_{3,sd}^{\beta=1}$", 40, 0.5, 3))
hists['opposite_ak8_tau32_mucr'] = hist.Hist("Events", dataset, gencat, systematic,
                                           jetpt_coarse, jetmass_coarse,
                                           hist.Bin("opposite_ak8_tau32", r"Jet $\tau_{32}$", 40, 0, 1))
hists['opposite_ak8_msd_mucr'] = hist.Hist("Events", dataset, gencat, systematic,
                                         jetpt_coarse, jetmass_coarse,
                                         hist.Bin("opposite_ak8_msd", r"Jet $\m_{sd}$", 40, 50, 200))
hists['opposite_ak4_leadingDeepCSV_mucr'] = hist.Hist("Events", dataset, gencat, systematic,
                                                    jetpt_coarse, jetmass_coarse,
                                                    hist.Bin("opposite_ak4_leadingDeepCSV", "Max(DeepCSV) (of $\leq4$ leading)", 40, 0, 1))

