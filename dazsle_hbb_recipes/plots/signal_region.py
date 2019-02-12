from fnal_column_analysis_tools import hist
from fnal_column_analysis_tools.hist import plot
from fnal_column_analysis_tools.util import numpy as np
from .systematics import jet_pt_systs,jet_weight_systs

def fill_plots_presel(dataset,gencat,systematic,leadingak8jet,weight,plots):
    genW = np.sign(weight)
    plots['sumw'].fill(dataset=dataset,
                       systematic=systematic,
                       sumw=genW)
    plots['hjetpt'].fill(dataset=dataset,
                         ak8_isHadronicV=gencat.sum(),
                         systematic=systematic,
                         ak8_pt=leadingak8jet.pt.sum(),
                         weight=weight)
    plots['hsculpt'].fill(dataset=dataset,
                          ak8_isHadronicV=gencat.sum(),
                          systematic=systematic,
                          ak8_pt=leadingak8jet.pt.sum(),
                          ak8_msd=leadingak8jet.msd_corr_8.sum(),
                          ak8_deepdoubleb=leadingak8jet.deepdoubleb.sum(),
                          ak8_deepdoublec=leadingak8jet.deepdoublec.sum(),
                          ak8_deepdoublecvb=leadingak8jet.deepdoublecvb.sum(),
                          weight=weight)

def fill_plots_sr(dataset,gencat,systematic,leadingak8jet,weight,plots):
    plots['hjetpt_sr'].fill(dataset=dataset,
                            ak8_isHadronicV=gencat.sum(),
                            systematic=systematic,
                            ak8_pt=leadingak8jet.pt.sum(),
                            weight=weight)
    plots['hsculpt_sr'].fill(dataset=dataset,
                             ak8_isHadronicV=gencat.sum(),
                             systematic=systematic,
                             ak8_pt=leadingak8jet.pt.sum(),
                             ak8_msd=leadingak8jet.msd_corr_8.sum(),
                             ak8_deepdoubleb=leadingak8jet.deepdoubleb.sum(),
                             ak8_deepdoublec=leadingak8jet.deepdoublec.sum(),
                             ak8_deepdoublecvb=leadingak8jet.deepdoublecvb.sum(),
                             weight=weight)

#preselection and signal region plots
def signal_region(gghbbcuts,
                  dataset,
                  gencat,
                  presel_weight, eventInfo, leadingak8jet,
                  looseMuons,looseElectrons,looseTaus,
                  hasTightVJet,
                  plots):
    #for name,vari in jet_systs:
    #    attr = "pt"
    #    if len(vari) > 0: attr += "_"+vari
    systematic = "central"
    srweight_nomet = ( #jet selection
                      ((leadingak8jet.pt > gghbbcuts.PTCUT).sum() > 0) &
                      ((leadingak8jet.msd_corr_8 > gghbbcuts.MASSCUT).sum() > 0) &
                      ((leadingak8jet.jtN2b1sdddt_8 < 0).sum() > 0) &
                      hasTightVJet &
                      #lepton vetos
                      (looseMuons.counts == 0) &
                      (looseElectrons.counts == 0) &
                      (looseTaus.counts == 0)
                     )
    #met selection (remove events with large MET, fake or real)
    pfmetweight = (eventInfo['pfmet'].sum() < gghbbcuts.METCUT)
    
    jetweight = leadingak8jet.weight.sum()
    #preselection
    weight = jetweight * presel_weight
    fill_plots_presel(dataset,gencat,systematic,leadingak8jet,weight,plots)
    
    #signal region no met cut
    weight_srnomet = weight * srweight_nomet
    plots['pfmet_nminus1_sr'].fill(dataset=dataset,
                                   ak8_isHadronicV=gencat.sum(),
                                   systematic=systematic,
                                   ak8_pt=leadingak8jet.pt.sum(),
                                   ak8_msd=leadingak8jet.msd_corr_8.sum(),
                                   pfmet=eventInfo['pfmet'].sum(),
                                   weight=weight_srnomet)
    
    #signal region variables
    weight_sr = weight_srnomet * pfmetweight
    fill_plots_sr(dataset,gencat,systematic,leadingak8jet,weight_sr,plots)




