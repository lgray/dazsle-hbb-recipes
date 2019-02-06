from fnal_column_analysis_tools import hist
from fnal_column_analysis_tools.hist import plot

jet_systs = {"central":"",
             "jes_up":"jes_up","jes_down":"jes_down",
             "jer_up":"jer_up","jer_down":"jer_down"}

muon_systs = {"central":"",
              "mutrig_up":"trig_up","mutrig_down":"trig_down",
              "muid_up":"id_up"    ,"muid_down":"id_down",
              "muiso_up":"iso_up"  ,"muiso_down":"iso_down",
              "mupu_up":"pu_up"    ,"mupu_down":"pu_down"}

#single muon control region plots
#                         cuts       mask        loose mus     ak4 jets         leading ak8
def muon_control_region(gghbbcuts,singlemuon_cr, looseMuons, ak4puppijet_pt30, leadingak8jet, plots):
    muons_singlemu_cr = looseMuons[singlemuon_cr]
    ak4_singlemu_cr = ak4puppijet_pt30[singlemuon_cr]
    ak8_singlemu_cr = leadingak8jet[singlemuon_cr]
    
    weight_single_mu = muons_singlemu_cr.weight_mu
    
    plots['h_ht'].fill(dataset="WQQ",
                       ht=ak4_singlemu_cr.pt.sum(),
                       weight=weight_single_mu.flatten())
    h_msd_ak8_muCR = plots['h_msd_ak8_muCR']
    h_msd_ak8_muCR.fill(dataset="WQQ",
                        region="muCR1",
                        passfail="none",
                        msd=ak8_singlemu_cr.msd_corr_8.max(),
                        weight=weight_single_mu.flatten())
    h_msd_ak8_muCR.fill(dataset="WQQ",
                        region="muCR2",
                        passfail="none",
                        msd=ak8_singlemu_cr.msd_corr_8[ak8_singlemu_cr.jdb_8 > gghbbcuts.DBTAGCUT].max(),
                        weight=weight_single_mu[ak8_singlemu_cr.jdb_8 > gghbbcuts.DBTAGCUT].flatten())
    h_msd_ak8_muCR.fill(dataset="WQQ",
                        region="muCR3",
                        passfail="none",
                        msd=ak8_singlemu_cr.msd_corr_8[ak8_singlemu_cr.jt21P_8 < 0.4].max(),
                        weight=weight_single_mu[ak8_singlemu_cr.jt21P_8 < 0.4].flatten())
    plots['h_t21ddt_ak8_muCR4'].fill(dataset="WQQ",
                                     tau21ddt=ak8_singlemu_cr.jt21P_8.flatten(),
                                     weight=weight_single_mu.flatten())
    
    singlemu_t21ddt_sel = (ak8_singlemu_cr.jt21P_8 < gghbbcuts.T21DDTCUT)
    print ak8_singlemu_cr.jt21P_8
    print singlemu_t21ddt_sel
    muons_singlemu_t21ddt_cr = muons_singlemu_cr[singlemu_t21ddt_sel]
    ak8_singlemu_t21ddt_cr = ak8_singlemu_cr[singlemu_t21ddt_sel]
    weight_single_mu_t21ddt = weight_single_mu[singlemu_t21ddt_sel]
    
    plots['h_dbtag_ak8_muCR4'].fill(dataset="WQQ",
                                    dbtag=ak8_singlemu_t21ddt_cr.jdb_8.flatten(),
                                    weight=weight_single_mu_t21ddt.flatten())
    h_msd_ak8_muCR.fill(dataset="WQQ",
                        region="muCR4",
                        passfail="none",
                        msd=ak8_singlemu_t21ddt_cr.msd_corr_8.max(),
                        weight=weight_single_mu_t21ddt.flatten())
    plots['h_pteta_ak8_muCR4'].fill(dataset="WQQ",
                                    pt=ak8_singlemu_t21ddt_cr.pt.flatten(),
                                    eta=ak8_singlemu_t21ddt_cr.eta.flatten(),
                                    weight=weight_single_mu_t21ddt.flatten())
    plots['h_pt_mu_muCR4'].fill(dataset="WQQ",
                                pt=muons_singlemu_t21ddt_sel.pt.flatten(),
                                weight=weight_single_mu_t21ddt.flatten())
    plots['h_eta_mu_muCR4'].fill(dataset="WQQ",
                                 eta=muons_singlemu_t21ddt_sel.eta.flatten(),
                                 weight=weight_single_mu_t21ddt.flatten())
    singlemu_t21ddt_pass_btag = (ak8_singlemu_t21ddt_cr.jdb_8 > gghbbcuts.DBTAGCUT)
    singlemu_t21ddt_fail_btag = ( (ak8_singlemu_t21ddt_cr.jdb_8 < gghbbcuts.DBTAGCUT) &
                                  (ak8_singlemu_t21ddt_cr.jdb_8 > gghbbcuts.DBTAGCUTMIN) )
    h_msd_ak8_muCR.fill(dataset="WQQ",
                        region="muCR4",
                        passfail="pass",
                        msd=ak8_singlemu_t21ddt_cr[singlemu_t21ddt_pass_btag].msd_corr_8.max(),
                        weight=weight_single_mu_t21ddt[singlemu_t21ddt_pass_btag].flatten())
    h_msd_ak8_muCR.fill(dataset="WQQ",
                        region="muCR4",
                        passfail="fail",
                        msd=ak8_singlemu_t21ddt_cr[singlemu_t21ddt_fail_btag].msd_corr_8.max(),
                        weight=weight_single_mu_t21ddt[singlemu_t21ddt_fail_btag].flatten())
    plots['h_msd_v_pt_ak8_top'].fill(dataset="WQQ",
                                     region="muCR4",
                                     passfail="pass",
                                     msd=ak8_singlemu_t21ddt_cr[singlemu_t21ddt_pass_btag].msd_corr_8.max(),
                                     pt=ak8_singlemu_t21ddt_cr[singlemu_t21ddt_pass_btag].pt.max(),
                                     weight=weight_single_mu_t21ddt[singlemu_t21ddt_pass_btag].flatten()
                                     )
    plots['h_msd_v_pt_ak8_top'].fill(dataset="WQQ",
                                     region="muCR4",
                                     passfail="fail",
                                     msd=ak8_singlemu_t21ddt_cr[singlemu_t21ddt_fail_btag].msd_corr_8.max(),
                                     pt=ak8_singlemu_t21ddt_cr[singlemu_t21ddt_fail_btag].pt.max(),
                                     weight=weight_single_mu_t21ddt[singlemu_t21ddt_fail_btag].flatten()
                                     )
                       
    muon_CR5 = (ak8_singlemu_cr.jdb_8 > 0.7) & (ak8_singlemu_cr.jt21P_8 < 0.4)
    plots['h_msd_ak8_muCR'].fill(dataset="WQQ",
                                 region="muCR5",
                                 passfail="none",
                                 msd=ak8_singlemu_cr[muon_CR5].msd_corr_8.max(),
                                 weight=weight_single_mu[muon_CR5].flatten())
    muon_CR6 = (ak8_singlemu_cr.jdb_8 > 0.7) & (ak8_singlemu_cr.jt21P_8 < gghbbcuts.T21DDTCUT)
    plots['h_msd_ak8_muCR'].fill(dataset="WQQ",
                                 region="muCR6",
                                 passfail="none",
                                 msd=ak8_singlemu_cr[muon_CR6].msd_corr_8.max(),
                                 weight=weight_single_mu[muon_CR6].flatten())
    
    mask_muCR4_N2 = (ak8_singlemu_cr.jtN2b1sdddt_8 < 0)
    muons_muCR4_N2 = muons_singlemu_cr[mask_muCR4_N2]
    ak8_muCR4_N2 = ak8_singlemu_cr[mask_muCR4_N2]
    weight_muCR4_N2 = weight_single_mu[mask_muCR4_N2]
    plots['h_mu_muCR4_N2'].fill(dataset="WQQ",
                                pt=muons_muCR4_N2.pt.max(),
                                eta=muons_muCR4_N2.eta.max(),
                                weight=weight_muCR4_N2.flatten())
    plots['h_ak8_muCR4_N2'].fill(dataset="WQQ",
                                 pt=ak8_muCR4_N2.pt.max(),
                                 eta=ak8_muCR4_N2.eta.max(),
                                 btag=ak8_muCR4_N2.jdb_8.max(),
                                 weight=weight_muCR4_N2.flatten())
    plots['h_msd_ak8_muCR4_N2'].fill(dataset="WQQ",
                                     passfail="none",
                                     systematic="central",
                                     msd=ak8_muCR4_N2.msd_corr_8.max(),
                                     weight=weight_muCR4_N2.flatten())
    mask_muCR4_N2_pass = (ak8_muCR4_N2.jdb_8 > gghbbcuts.DBTAGCUT)
    mask_muCR4_N2_fail = ( (ak8_muCR4_N2.jdb_8 < gghbbcuts.DBTAGCUT) &
                           (ak8_muCR4_N2.jdb_8 > gghbbcuts.DBTAGCUTMIN) )
    ak8_muCR4_N2_pass = ak8_muCR4_N2[mask_muCR4_N2_pass]
    muons_muCR4_N2_pass = muons_muCR4_N2[mask_muCR4_N2_pass]
    ak8_muCR4_N2_fail = ak8_muCR4_N2[mask_muCR4_N2_fail]
    muons_muCR4_N2_fail = muons_muCR4_N2[mask_muCR4_N2_fail]
    for name,vari in muon_systs.items():
        attr = "weight_mu"
        if len(vari) > 0: attr += vari
        plots['h_msd_ak8_muCR4_N2'].fill(dataset="WQQ",
                                         passfail="pass",
                                         systematic = name,
                                         msd = ak8_muCR4_N2_pass.msd_corr_8.max(),
                                         weight = getattr(muons_muCR4_N2_pass,attr).flatten()
                                        )
        plots['h_msd_vs_pt_ak8_muCR4_N2'].fill(dataset="WQQ",
                                               passfail="pass",
                                               systematic = name,
                                               msd = ak8_muCR4_N2_pass.msd_corr_8.max(),
                                               pt = ak8_muCR4_N2_pass.pt.max(),
                                               weight = getattr(muons_muCR4_N2_pass,attr).flatten()
                                              )
        plots['h_msd_ak8_muCR4_N2'].fill(dataset="WQQ",
                                         passfail="fail",
                                         systematic = name,
                                         msd = ak8_muCR4_N2_fail.msd_corr_8.max(),
                                         weight = getattr(muons_muCR4_N2_fail,attr).flatten()
                                        )
        plots['h_msd_vs_pt_ak8_muCR4_N2'].fill(dataset="WQQ",
                                               passfail="fail",
                                               systematic = name,
                                               msd = ak8_muCR4_N2_fail.msd_corr_8.max(),
                                               pt = ak8_muCR4_N2_fail.pt.max(),
                                               weight = getattr(muons_muCR4_N2_pass,attr).flatten()
                                              )
