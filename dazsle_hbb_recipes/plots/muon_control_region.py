from fnal_column_analysis_tools import hist
from fnal_column_analysis_tools.hist import plot
from fnal_column_analysis_tools.util import awkward
from fnal_column_analysis_tools.util import numpy as np
import math
from .systematics import jet_pt_systs,muon_weight_systs

def matchByDPhi(first,second,deltaPhiCut=2./3.*math.pi):
    args = first.phi._argcross(second.phi)
    argsnested = awkward.JaggedArray.fromcounts(first.phi.counts,
                                                awkward.JaggedArray.fromcounts(first.phi._broadcast(second.phi.counts).flatten(),
                                                                       args._content))
    phi0s  = first.phi.content[argsnested.content.content.i0]
    phi1s  = second.phi.content[argsnested.content.content.i1]
    offsets_outer = argsnested.offsets
    offsets_inner = argsnested.content.offsets
    dphis = (phi0s - phi1s + math.pi) % (2*math.pi) - math.pi
    passdphi = (np.abs(dphis) < deltaPhiCut)
    passdphi = awkward.JaggedArray.fromoffsets(offsets_inner,passdphi)
    return awkward.JaggedArray.fromoffsets(offsets_outer,passdphi.any())

def antiMatchByDR(first,second,deltaRCut=0.8):
    drCut2 = deltaRCut**2
    args = first.eta._argcross(second.eta)
    argsnested = awkward.JaggedArray.fromcounts(first.eta.counts,
                                                awkward.JaggedArray.fromcounts(first.eta._broadcast(second.eta.counts).flatten(),
                                                                       args._content))
    eta0s  = first.eta.content[argsnested.content.content.i0]
    eta1s  = second.eta.content[argsnested.content.content.i1]
    phi0s  = first.phi.content[argsnested.content.content.i0]
    phi1s  = second.phi.content[argsnested.content.content.i1]
    offsets_outer = argsnested.offsets
    offsets_inner = argsnested.content.offsets
    detas = np.abs(eta0s - eta1s)
    dphis = (phi0s - phi1s + math.pi) % (2*math.pi) - math.pi
    passdr = ((detas**2 + dphis**2) > drCut2)
    passdr = awkward.JaggedArray.fromoffsets(offsets_inner,passdr)
    return awkward.JaggedArray.fromoffsets(offsets_outer,passdr.any())

def fill_plots_mucr(dataset,gencat,systematic,leadingak8jet,weight,plots):    
    plots['hjetpt_mucr'].fill(dataset=dataset,
                              ak8_isHadronicV=gencat.sum(),
                              systematic=systematic,
                              ak8_pt=leadingak8jet.pt.sum(),
                              weight=weight)
    plots['hsculpt_mucr'].fill(dataset=dataset,
                               ak8_isHadronicV=gencat.sum(),
                               systematic=systematic,
                               ak8_pt=leadingak8jet.pt.sum(),
                               ak8_msd=leadingak8jet.msd_corr_8.sum(),
                               ak8_deepdoubleb=leadingak8jet.deepdoubleb.sum(),
                               ak8_deepdoublec=leadingak8jet.deepdoublec.sum(),
                               ak8_deepdoublecvb=leadingak8jet.deepdoublecvb.sum(),
                               weight=weight)

#muon control region plots
def muon_control_region(gghbbcuts,
                        dataset,
                        gencat,
                        presel_weight, eventInfo, leadingak8jet,
                        ak4jets_Mbtag,
                        looseMuons,looseElectrons,looseTaus,
                        hasTightVJet,
                        plots):
    leadingLooseMuon = looseMuons[looseMuons.pt.argmax()]
    
    AK8jet_muon_matches = leadingak8jet.fastmatch(leadingLooseMuon,matchfunc=matchByDPhi)
    AK8jet_AK4Mbjet_antimatches = leadingak8jet.fastmatch(ak4jets_Mbtag,matchfunc=antiMatchByDR)
    
    #for name,vari in jet_systs:
    #    attr = "pt"
    #    if len(vari) > 0: attr += "_"+vari
    systematic = "central"
    mucrweight = ( #jet selection
                    ((leadingak8jet.pt > gghbbcuts.PTCUTMUCR).sum() > 0) &
                    ((leadingak8jet.msd_corr_8 > gghbbcuts.MASSCUT).sum() > 0) &
                    hasTightVJet &
                    #muon selection
                    ((np.abs(leadingLooseMuon.eta) < 2.1).sum() > 0) &
                    ((leadingLooseMuon.pt > gghbbcuts.MUONPTCUT).sum() > 0) &
                    #matches and cross cleaning
                    (AK8jet_muon_matches.sum() == 0) &
                    (AK8jet_AK4Mbjet_antimatches.sum() > 0) &
                    #lepton vetos
                    (looseMuons.counts == 1) &
                    (looseElectrons.counts == 0) &
                    (looseTaus.counts == 0)
                  )
    
    #fill plots
    weight = leadingak8jet.weight.sum() * presel_weight
    weight_mucr = weight * mucrweight
    fill_plots_mucr(dataset,gencat,systematic,leadingak8jet,weight_mucr,plots)
