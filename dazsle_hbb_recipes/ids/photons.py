"""photon IDs"""

import numpy as np
from awkward import JaggedArray
from fnal_column_analysis_tools.analysis_objects import JaggedCandidateArray

def calcPhotonEffAreaIso(photon,rho):
    chHadIso  = np.maximum(photon.chHadIso  - rho.content*photon.EA_CHad,0.);
    neuHadIso = np.maximum(photon.neuHadIso - rho.content*photon.EA_NHad,0.);
    phoIso    = np.maximum(photon.gammaIso  - rho.content*photon.EA_Pho, 0.);
    
    maxNeuHadLooseIso = np.zeros_like(photon.pt.content)
    maxNeuHadMediumIso = np.zeros_like(photon.pt.content)
    maxNeuHadTightIso = np.zeros_like(photon.pt.content)
    etaCut = (photon.absSCEta.content < 1.479)
    notEta = ~etaCut
    #common bits
    ptshape = np.zeros_like(photon.pt.content)
    ptshape[etaCut] = photon.pt.content[etaCut]*(0.0148+0.000017*photon.pt.content[etaCut])
    ptshape[notEta] = photon.pt.content[notEta]*(0.0163+0.000014*photon.pt.content[notEta])
    
    #loose
    maxNeuHadLooseIso[etaCut] = 10.910 + ptshape[etaCut]
    maxNeuHadLooseIso[notEta] = 5.931  + ptshape[notEta]
    maxNeuHadLooseIso = JaggedArray.fromoffsets(photon.pt.offsets,maxNeuHadLooseIso)
    
    #medium
    maxNeuHadMediumIso[etaCut] = 2.725 + ptshape[etaCut]
    maxNeuHadMediumIso[notEta] = 1.715 + ptshape[notEta]
    maxNeuHadMediumIso = JaggedArray.fromoffsets(photon.pt.offsets,maxNeuHadMediumIso)
    
    #tight
    maxNeuHadTightIso[etaCut] = 0.264 + ptshape[etaCut]
    maxNeuHadTightIso[notEta] = 0.586 + ptshape[notEta]
    maxNeuHadTightIso = JaggedArray.fromoffsets(photon.pt.offsets,maxNeuHadTightIso)
    
    return chHadIso,neuHadIso,phoIso,maxNeuHadLooseIso,maxNeuHadMediumIso,maxNeuHadTightIso

def passPhoLooseSel(photon):
    outs = np.ones_like(photon.absSCEta.content,dtype=np.bool)
    etaCut = ( photon.absSCEta.content < 1.479 )
    notEta = ~etaCut
    #barrel
    outs[etaCut] &= ( (photon.hovere.content[etaCut] <= 0.0597      ) &
                      (photon.sieie.content[etaCut] <= 0.01031      ) &
                      (photon.chIsoCorr.content[etaCut] <= 1.295    ) &
                      (photon.nhIsoCorr.content[etaCut] <=
                           photon.neuIsoLooseMax.content[etaCut]    ) &
                      (photon.gamIsoCorr.content[etaCut] <= 3.630 +
                                   0.0047*photon.pt.content[etaCut] )   )
    #endcap
    outs[notEta] &= ( (photon.hovere.content[notEta] <= 0.0481      ) &
                      (photon.sieie.content[notEta] <= 0.03013      ) &
                      (photon.chIsoCorr.content[notEta] <= 1.011    ) &
                      (photon.nhIsoCorr.content[notEta] <=
                            photon.neuIsoLooseMax.content[notEta]   ) &
                      (photon.gamIsoCorr.content[notEta] <= 6.641 +
                                   0.0034*photon.pt.content[notEta] )   )
    outs = JaggedArray.fromoffsets(photon.pt.offsets,outs)
    return outs

def passPhoMediumSel(photon):
    outs = np.ones_like(photon.absSCEta.content,dtype=np.bool)
    etaCut = ( photon.absSCEta.content < 1.479 )
    notEta = ~etaCut
    #barrel
    outs[etaCut] &= ( (photon.hovere.content[etaCut] <= 0.0396      ) &
                      (photon.sieie.content[etaCut] <= 0.01022      ) &
                      (photon.chIsoCorr.content[etaCut] <= 0.441    ) &
                      (photon.nhIsoCorr.content[etaCut] <=
                          photon.neuIsoMediumMax.content[etaCut]    ) &
                      (photon.gamIsoCorr.content[etaCut] <= 2.571 +
                                   0.0047*photon.pt.content[etaCut] )   )
    #endcap
    outs[notEta] &= ( (photon.hovere.content[notEta] <= 0.0219      ) &
                      (photon.sieie.content[notEta] <= 0.03001      ) &
                      (photon.chIsoCorr.content[notEta] <= 0.442    ) &
                      (photon.nhIsoCorr.content[notEta] <=
                            photon.neuIsoMediumMax.content[notEta]  ) &
                      (photon.gamIsoCorr.content[notEta] <= 3.863 +
                                   0.0034*photon.pt.content[notEta] )   )
    outs = JaggedArray.fromoffsets(photon.pt.offsets,outs)
    return outs

def passPhoTightSel(photon):
    outs = np.ones_like(photon.absSCEta.content,dtype=np.bool)
    etaCut = ( photon.absSCEta.content < 1.479 )
    notEta = ~etaCut
    #barrel
    outs[etaCut] &= ( (photon.hovere.content[etaCut] <= 0.0269      ) &
                      (photon.sieie.content[etaCut] <= 0.00994      ) &
                      (photon.chIsoCorr.content[etaCut] <= 0.202    ) &
                      (photon.nhIsoCorr.content[etaCut] <=
                           photon.neuIsoTightMax.content[etaCut]    ) &
                      (photon.gamIsoCorr.content[etaCut] <= 2.362 +
                                   0.0047*photon.pt.content[etaCut] )   )
    #endcap
    outs[notEta] &= ( (photon.hovere.content[notEta] <= 0.0213      ) &
                      (photon.sieie.content[notEta] <= 0.03000      ) &
                      (photon.chIsoCorr.content[notEta] <= 0.034    ) &
                      (photon.nhIsoCorr.content[notEta] <=
                            photon.neuIsoTightMax.content[notEta]   ) &
                      (photon.gamIsoCorr.content[notEta] <= 2.617 +
                                   0.0034*photon.pt.content[notEta] )   )
    outs = JaggedArray.fromoffsets(photon.pt.offsets,outs)
    return outs
