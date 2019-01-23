""" dazsle gghbb analysis plots """

from fnal_column_analysis_tools import hist
from fnal_column_analysis_tools.hist import plot

#common histograms
fBosonPt_fbweight = hist.Hist("Events / bin",
                              hist.Cat("dataset","Dataset"),
                              hist.Bin("pt",r"Boson $p_{T}$ (GeV)",100, 0, 1000))
fBosonPt_weight = hist.Hist("Events / bin",
                            hist.Cat("dataset","Dataset"),
                            hist.Bin("pt",r"Boson $p_{T}$ (GeV)",100, 0, 1000))
nPV = hist.Hist("Events / bin",
                hist.Cat("dataset","Dataset"),
                hist.Bin("npv",r"Number of PV",100, 0, 100),
                weight="weight")

gghbb_ak8pt_binning = [450, 500, 550, 600, 675, 800, 1000]

#main analysis region
h_rhop_v_t21_ak8 = hist.Hist("Events / bin",
                             hist.Bin("rho_ddt",r"AK8 $\rho^{DDT}$",15,-5,10),
                             hist.Bin("tau21",r"AK8 $\tau_{21}$",25, 0, 1.5)
                             )
h_rhop_v_t21_ak15 = hist.Hist("Events / bin",
                              hist.Bin("rho_ddt",r"AK8 $\rho^{DDT}$",15,-5,10),
                              hist.Bin("tau21",r"AK8 $\tau_{21}$",25, 0, 1.5)
                              )
#systs JES, JER, trigger, Pu, matched/unmatched
#regions topR6_N2,QGquark,QGgluon
h_msd_v_pt_ak8 = hist.Hist("Events / bin",
                           hist.Cat("region","Analysis Region"),
                           hist.Cat("systematic","Systematic Variations"),
                           hist.Bin("msd",r"AK8 $m_{SD}^{PUPPI}$ (GeV)",23,40,201),
                           hist.Bin("pt",r"AK8 $p_{T}$ (GeV)",gghbb_ak8pt_binning),
                           #, weight="weight",
                           )

#old analysis regions using tau21
#systs JES, JER, trigger, Pu, matched/unmatched
#regions topR6
h_msd_v_pt_ak8_old = hist.Hist("Events / bin",
                               hist.Cat("region","Analysis Region"),
                               hist.Cat("systematic","Systematic Variations"),
                               hist.Bin("msd",r"AK8 $m_{SD}^{PUPPI}$ (GeV)",23,40,201),
                               hist.Bin("pt",r"AK8 $p_{T}$ (GeV)",gghbb_ak8pt_binning),
                               #, weight="weight",
                               )

#top control region histograms
# pass / fail
# systs = jes, jer
h_msd_ak8_topR6 = hist.Hist("Events / bin",
                            hist.Cat("pass","Pass or Fail selection"),
                            hist.Cat("systematic","Systematic Variations"),
                            hist.Bin("msd",r"AK8 $m_{SD}^{PUPPI}$ (GeV)",23,40,201),
                            #, weight="weight",
                            )

#regions topR1-topR7, bbleading, muCR4, bbleading+muCR4
h_msd_v_pt_ak8_top = hist.Hist("Events / bin",
                               hist.Cat("region","Analysis Region"),
                               hist.Cat("pass","Pass or Fail Selection"),
                               hist.Bin("msd",r"AK8 $m_{SD}^{PUPPI}$ (GeV)",23,40,201),
                               hist.Bin("pt",r"AK8 $p_{T}$ (GeV)",gghbb_ak8pt_binning),
                               #, weight="weight",
                               )

#muon control region histograms
#soft drop mass
#pass / fail
#systs = jes, jer, muTrig, muId, muIso, puWeight
h_pt_mu_muCR4_N2 = hist.Hist("Events / bin",
                             hist.Bin("pt",r"leading muon $p_{T}$ (GeV)",50,30,500),
                             hist.Bin("eta",r"leading muon $\eta$ (GeV)",50,-2.5,2.5)
                             )
h_ak8_muCR4_N2 = hist.Hist("Events / bin",
                           hist.Bin("pt",r"AK8 leading $p_{T}$ (GeV)",50,300,2100),
                           hist.Bin("eta",r"AK8 leading $\eta$",50, -3, 3),
                           hist.Bin("btag",r"$p_{T}$-leading double b-tag",40,-1,1)
                           )
h_msd_ak8_muCR4_N2 = hist.Hist("Events / bin",
                               hist.Cat("pass","Pass or Fail selection"),
                               hist.Cat("systematic","Systematic Variations"),
                               hist.Bin("msd",r"AK8 m_{SD}^{PUPPI} (GeV)",23,40,201),
                               #, weight="weight",
                               )

#extras
h_Cuts = hist.Hist("Events / bin",hist.Bin("cut",r"$Cut_{i}$",8,0,8))
h_n_ak4 = hist.Hist("Events / bin",hist.Bin("n",r"AK4 $n_{jets}$, $p_{T}$ > 30 GeV",20,0,20))
h_ht = hist.Hist("Events / bin", hist.Bin("ht","HT (GeV)",50,300,2100))
h_pt_bb = hist.Hist("Events / bin",
                    hist.Bin("pt",r"AK8 leading $p_{T}$ (GeV)",50,300,2100),
                    hist.Bin("dbtag","double b-tag",40,-1,1),
                    hist.Bin("msd",r"AK8 $m_{SD}^{PUPPI}$ (GeV)",30,40,250),
                    )
h_n_ak4fwd = hist.Hist("Events / bin", hist.Bin("n",r"AK4 $n_{jets}, p_{T} > 30 GeV, 2.5<|\eta|<4.5$",20,0,20))
h_n_ak4L = hist.Hist("Events / bin", hist.Bin("n",r"AK4 $n_{L b-tags}, \DeltaR > 0.8, p_{T} > 40 GeV$",20,0,20))
h_n_ak4M = hist.Hist("Events / bin", hist.Bin("n",r"AK4 $n_{M b-tags}, \DeltaR > 0.8, p_{T} > 40 GeV$",20,0,20))
h_n_ak4T = hist.Hist("Events / bin", hist.Bin("n",r"AK4 $n_{T b-tags}, \DeltaR > 0.8, p_{T} > 40 GeV$",20,0,20))
h_n_ak4dR0p8 = hist.Hist("Events / bin", hist.Bin("n",r"AK4 $n_{jets}, \DeltaR > 0.8, p_{T} > 30 GeV$",20,0,20))
h_isolationCA15 = hist.Hist("Events / bin", hist.Bin("iso",r"AK8/CA15 $p_{T}$ ratio",50,0.5,1.5))
h_met = hist.Hist("Events / bins", hist.Bin("met",r"$E_{T}^{miss}$ (GeV)",50,0,500))
h_pteta_ak8 = hist.Hist("Events / bins",
                        hist.Bin("pt",r"AK8 leading $p_{T}$ (GeV)",50,300,2100),
                        hist.Bin("eta",r"AK8 leading \eta",50,-3,3)
                        )
#order leading, sub-leading, sub-sub-leading
h_pt_ak8 = hist.Hist("Events / bin",
                     hist.Cat("order","pt-rank"),
                     hist.Bin("pt","AK8 $p_{T}$ (GeV)",45,300,2100)
                     )
h_msdrho_ak8 = hist.Hist("Events / bin",
                         hist.Bin("msd",r"$p_{T}$-leading $m_{SD}$ (GeV)",23,40,201),
                         hist.Bin("rho",r"$p_{T}$-leading  $\rho=log(m_{SD}^{2}/p_{T}^{2})$",50,-7,-1)
                         )
#cut levels : raw,corr'd,dbtagCut, t21ddtCut, t21ddtCut_inc, N2Cut
h_msd_cuts_ak8 = hist.Hist("Events / bin",
                           hist.Cat("cut","cut level"),
                           hist.Bin("msd",r"AK8 $m_{SD}^{PUPPI}$ (GeV)",23, 40, 201))
h_dbtag_ak8 = hist.Hist("Events / bin",
                        hist.Cat("order","pt-rank"),
                        hist.Bin("dbtag","double b-tag",40,-1,1)
                        )
h_n_unMatchedAK4 = hist.Hist("Events / bin",
                             hist.Bin("n","Number of non matched AK4 jets",7,0,7))
h_Mqq = hist.Hist("Events / bin",hist.Bin("mqq","Dijet Mass (GeV)",50,0,3000))
h_QGLR = hist.Hist("Events / bin",hist.Bin("gqlr","QGLR",10,0,1))
h_Deta_qq = hist.Hist("Events / bin",hist.Bin("deta","max Deta qq",30,0,10))
h_t21_ak8 = hist.Hist("Events / bin",
                      hist.Bin("tau21",r"AK8 $\tau_{21}$",25,0,1.5),
                      hist.Bin("tau21ddt",r"AK8 $\tau_{21}^{DDT}$",25,0,1.5)
                      )
h_t32_ak8 = hist.Hist("Events / bin",hist.Bin("tau32",r"AK8 $\tau_{32}$",20,0,1.5))
h_t32_ak8_t21ddtcut = hist.Hist("Events / bin",hist.Bin("tau32",r"AK8 $\tau_{32}$",20,0,1.5))
h_n2b1sd_ak8 =  hist.Hist("Events / bin",hist.Bin("N2",r"AK8 $N_{2}^{1}$ (SD)", 25, -0.5, 0.5))
h_n2b1sdddt_ak8 =  hist.Hist("Events / bin",hist.Bin("N2",r"AK8 $N_{2}^{1,DDT}$ (SD)", 25, -0.5, 0.5))
h_n2b1sdddt_ak8_aftercut =  hist.Hist("Events / bin",hist.Bin("N2",r"AK8 $N_{2}^{1}$ (SD)", 25, -0.5, 0.5))
h_dbtag_ak8_aftercut = hist.Hist("Events / bin",hist.Bin("dbtag",r"$p_{T}$-leading double-b tagger",33,-1,1))

h_msd_ak8_raw_SR = hist.Hist("Events / bin",
                             hist.Cat("pass","pass or fail"),
                             hist.Bin("msdraw",r"AK8 $m_{SD}^{PUPPI}$ no corr (GeV)",23,40,201))
#pass or fail
#systematics JES, JER
h_msd_ak8_topR6 = hist.Hist("Events / bin",
                            hist.Cat("pass","pass or fail"),
                            hist.Cat("systematic","Systematic Variation"),
                            hist.Bin("msd",r"AK8 $m_{SD}^{PUPPI}$ (GeV)",23,40,201))
#ak4 b-tag L,M,T, pt > 100 / 150
h_n_ak4 = hist.Hist("Event / bin",
                    hist.Cat("ptcat","pt > 100 or 150"),
                    hist.Cat("bcat","L M T b-tag"),
                    hist.Bin("n",r"AK4 $n_{b-tags}$, $\DeltaR > 0.8$, $p_{T} > 100/150$ GeV",10,0,10))

#region = muCR1-6, bbleading R4
#pass fail for CR4 only
h_msd_ak8_muCR = hist.Hist("Events / bin",
                           hist.Cat("region","Analysis Region"),
                           hist.Cat("pass","pass or fail"),
                           hist.Bin("msd",r"AK8 $m_{SD}^{PUPPI}$ (GeV)",23,40,201)
                           )
# more muCR4 breakdown plots
h_pt_mu_muCR4 = hist.Hist("Events / bin",hist.Bin("pt",r"leading muon $p_{T}$ (GeV)",50,30,500))
h_eta_mu_muCR4 = hist.Hist("Events / bin",hist.Bin("eta",r"leading muon $\eta$",50,-2.5,2.5))
h_pteta_ak8_muCR4 = hist.Hist("Events / bin",
                              hist.Bin("pt",r"AK8 leading $p_{T} (GeV)",50,300,2100),
                              hist.Bin("eta",r"AK8 leading \eta",50,-3,3))
h_dbtag_t21ddt_ak8_muCR4 = hist.Hist("Events / bin",
                                     hist.Bin("dbtag",r"$p_{T}$-leading double b-tag",40,-1,1),
                                     hist.Bin("tau21ddt",r"AK8 $#tau_{21}^{DDT}$",25,0,1.5))

#region = topR1-7, R6 cut val 0p4-0p7, bbleading R6
#pass fail for R2-7
h_msd_ak8_topCR = hist.Hist("Events / bin",
                            hist.Cat("region","Analysis Region"),
                            hist.Cat("cut","Discriminator Cut Value"),
                            hist.Cat("pass","pass or fail"),
                            hist.Bin("msd",r"AK8 $m_{SD}^{PUPPI}$ (GeV)",23,40,201)
                            )
