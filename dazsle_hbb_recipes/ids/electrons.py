"""electron IDs"""

import numpy as np
from awkward import JaggedArray
from fnal_column_analysis_tools.analysis_objects import JaggedCandidateArray

kEcalDriven    = 1 << 0
kTrackerDriven = 1 << 1

def calcElectronEffAreaIso(electron,rho):
    return electron.chHadIso + np.maximum(0.0, (electron.gammaIso + electron.neuHadIso - rho.content*electron.EffArea) )

def passEleVetoSel(electron):
    outs = (electron.isConv.content == 0)
    etaCut = (np.abs(electron.scEta.content) < 1.479)
    notEta = ~etaCut
    #barrel
    outs[etaCut] &= ( ( electron.effAreaIso.content[etaCut] <
                       0.175*electron.pt.content[etaCut]      ) &
                     ( electron.sieie.content[etaCut] < 0.01150             ) &
                     ( np.abs(electron.dEtaInSeed.content[etaCut]) < 0.00749) &
                     ( np.abs(electron.dPhiIn.content[etaCut]) < 0.22800    ) &
                     ( electron.hovere.content[etaCut] < 0.35600            ) &
                     ( np.abs(1.0 - electron.eoverp.content[etaCut]) <
                      0.29900*electron.ecalEnergy.content[etaCut]   ) &
                     ( electron.nMissingHits.content[etaCut] <= 2           )   )
    #endcap
    outs[notEta] &= ( ( electron.effAreaIso.content[notEta] <
                            0.159*electron.pt.content[notEta]      ) &
                       ( electron.sieie.content[notEta] < 0.03700             ) &
                       ( np.abs(electron.dEtaInSeed.content[notEta]) < 0.00895) &
                       ( np.abs(electron.dPhiIn.content[notEta]) < 0.21300    ) &
                       ( electron.hovere.content[notEta] < 0.21100            ) &
                       ( np.abs(1.0 - electron.eoverp.content[notEta]) <
                                0.15000*electron.ecalEnergy.content[notEta]   ) &
                       ( electron.nMissingHits.content[notEta] <= 3           )   )
    outs = JaggedArray.fromoffsets(electron.pt.offsets,outs)
    return outs

def passEleLooseSel(electron):
    outs = (electron.isConv.content == 0)
    etaCut = (np.abs(electron.scEta.content) < 1.479)
    notEta = ~etaCut
    #barrel
    outs[etaCut] &= ( ( electron.effAreaIso.content[etaCut] <
                                      0.0994*electron.pt.content[etaCut]     ) &
                      ( electron.sieie.content[etaCut] < 0.01100             ) &
                      ( np.abs(electron.dEtaInSeed.content[etaCut]) < 0.00477) &
                      ( np.abs(electron.dPhiIn.content[etaCut]) < 0.22200    ) &
                      ( electron.hovere.content[etaCut] < 0.29800            ) &
                      ( np.abs(1.0 - electron.eoverp.content[etaCut]) <
                               0.24100*electron.ecalEnergy.content[etaCut]   ) &
                      ( electron.nMissingHits.content[etaCut] <= 1           )   )
    #endcap
    outs[notEta] &= ( ( electron.effAreaIso.content[notEta] <
                                      0.107*electron.pt.content[notEta]      ) &
                      ( electron.sieie.content[notEta] < 0.03700             ) &
                      ( np.abs(electron.dEtaInSeed.content[notEta]) < 0.00895) &
                      ( np.abs(electron.dPhiIn.content[notEta]) < 0.21300    ) &
                      ( electron.hovere.content[notEta] < 0.10100            ) &
                      ( np.abs(1.0 - electron.eoverp.content[notEta]) <
                               0.14000*electron.ecalEnergy.content[notEta]   ) &
                      ( electron.nMissingHits.content[notEta] <= 1           )   )
    outs = JaggedArray.fromoffsets(electron.pt.offsets,outs)
    return outs

def passEleMediumSel(electron):
    outs = (electron.isConv.content == 0)
    etaCut = (np.abs(electron.scEta.content) < 1.479)
    notEta = ~etaCut
    #barrel
    outs[etaCut] &= ( ( electron.effAreaIso.content[etaCut] <
                       0.0695*electron.pt.content[etaCut]     ) &
                      ( electron.sieie.content[etaCut] < 0.00998             ) &
                      ( np.abs(electron.dEtaInSeed.content[etaCut]) < 0.00311) &
                      ( np.abs(electron.dPhiIn.content[etaCut]) < 0.10300    ) &
                      ( electron.hovere.content[etaCut] < 0.25300            ) &
                      ( np.abs(1.0 - electron.eoverp.content[etaCut]) <
                               0.13400*electron.ecalEnergy.content[etaCut]   ) &
                      ( electron.nMissingHits.content[etaCut] <= 1           )   )
    #endcap
    outs[notEta] &= ( ( electron.effAreaIso.content[notEta] <
                                      0.0821*electron.pt.content[notEta]     ) &
                      ( electron.sieie.content[notEta] < 0.02980             ) &
                      ( np.abs(electron.dEtaInSeed.content[notEta]) < 0.00609) &
                      ( np.abs(electron.dPhiIn.content[notEta]) < 0.04500    ) &
                      ( electron.hovere.content[notEta] < 0.08780            ) &
                      ( np.abs(1.0 - electron.eoverp.content[notEta]) <
                               0.13000*electron.ecalEnergy.content[notEta]   ) &
                      ( electron.nMissingHits.content[notEta] <= 1           )   )
    outs = JaggedArray.fromoffsets(electron.pt.offsets,outs)
    return outs

def passEleTightSel(electron):
    outs = (electron.isConv.content == 0)
    etaCut = (np.abs(electron.scEta.content) < 1.479)
    notEta = ~etaCut
    #barrel
    outs[etaCut] &= ( ( electron.effAreaIso.content[etaCut] <
                                      0.0588*electron.pt.content[etaCut]     ) &
                      ( electron.sieie.content[etaCut] < 0.00998             ) &
                      ( np.abs(electron.dEtaInSeed.content[etaCut]) < 0.00308) &
                      ( np.abs(electron.dPhiIn.content[etaCut]) < 0.08160    ) &
                      ( electron.hovere.content[etaCut] < 0.04140            ) &
                      ( np.abs(1.0 - electron.eoverp.content[etaCut]) <
                               0.12900*electron.ecalEnergy.content[etaCut]   ) &
                      ( electron.nMissingHits.content[etaCut] <= 1           )   )
    #endcap
    outs[notEta] &= ( ( electron.effAreaIso.content[notEta] <
                                      0.0571*electron.pt.content[notEta]     ) &
                      ( electron.sieie.content[notEta] < 0.02920             ) &
                      ( np.abs(electron.dEtaInSeed.content[notEta]) < 0.00605) &
                      ( np.abs(electron.dPhiIn.content[notEta]) < 0.03940    ) &
                      ( electron.hovere.content[notEta] < 0.06410            ) &
                      ( np.abs(1.0 - electron.eoverp.content[notEta]) <
                               0.12900*electron.ecalEnergy.content[notEta]   ) &
                      ( electron.nMissingHits.content[notEta] <= 1           )   )
    outs = JaggedArray.fromoffsets(electron.pt.offsets,outs)
    return outs

def passEleHEEPSel(electron,rho,met):
    outs = ( (electron.typeBits.content & kEcalDriven) > 0 )
    outs &= (met <= 35)
    absSCEta = np.abs(electron.scEta.content)
    etaBar = (absSCEta < 1.442)
    etaEnd = ( (absSCEta > 1.566) & (absSCEta < 2.5) )
    notEta = ~(etaBar | etaEnd)
    outs[notEta] = False #anything in the gap dies
    # setup all the hcalDepth1Iso cuts
    ones = electron.typeBits.ones_like()
    hcalDepth1Cuts = np.zeros_like(ones.content)
    mets = (met.content*ones).content
    rhos = (rho.content*ones).content
    hcalDepth1Cuts[etaBar] = 2.0+0.03*mets[etaBar]+0.28*rhos[etaBar]
    hcalDepth1Cuts[etaEnd & (mets < 50)] = 2.5+0.28*rhos[etaEnd & (mets < 50)]
    hcalDepth1Cuts[etaEnd & (mets >= 50)] = ( 2.5 + 0.03*(mets[etaEnd & (mets >= 50)]-50) +
                                             0.28*rhos[etaEnd & (mets >= 50)] )
        
    #barrel
    outs[etaBar] &= ( ( np.abs(electron.dEtaInSeed.content[etaBar]) < 0.00600) &
                      ( np.abs(electron.dPhiIn.content[etaBar]) < 0.06000    ) &
                      ( electron.hovere.content[etaBar] <
                           (1.0/(electron.ecalEnergy.content[etaBar])+0.05)  ) &
                      (electron.e2x5.content[etaBar]/
                            electron.e5x5.content[etaBar] > 0.94             ) &
                      (electron.e1x5.content[etaBar]/
                            electron.e5x5.content[etaBar] > 0.83             ) &
                      ( electron.hcalDepth1Iso.content[etaBar] <
                                     hcalDepth1Cuts[etaBar]                  ) &
                      ( electron.trkIso.content[etaBar] < 5                  ) &
                      ( electron.nMissingHits.content[etaBar] <= 1           ) &
                      ( electron.d0.content[etaBar] <= 0.02                  )   )
    #endcap
    outs[etaEnd] &= ( ( np.abs(electron.dEtaInSeed.content[etaEnd]) < 0.00400) &
                      ( np.abs(electron.dPhiIn.content[etaEnd]) < 0.06000    ) &
                      ( electron.hovere.content[etaEnd] <
                           (5.0/(electron.ecalEnergy.content[etaEnd])+0.05)  ) &
                      (electron.sieie.content[etaEnd] < 0.03                 ) &
                      ( electron.hcalDepth1Iso.content[etaEnd] <
                                     hcalDepth1Cuts[etaEnd]                  ) &
                      ( electron.trkIso.content[etaEnd] < 5                  ) &
                      ( electron.nMissingHits.content[etaEnd] <= 1           ) &
                      ( electron.d0.content[etaEnd] <= 0.05                  )   )
    outs = JaggedArray.fromoffsets(electron.pt.offsets,outs)
    return outs
