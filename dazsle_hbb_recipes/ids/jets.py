"""jet IDs"""

import numpy as np
from awkward import JaggedArray
from fnal_column_analysis_tools.analysis_objects import JaggedCandidateArray

def passLooseJetSel(jet):
    outs = np.ones_like(jet.pt.content,dtype=np.bool)
    absEta = np.abs(jet.eta.content)
    etaVFor = (absEta <= 3.0)
    etaFor  = (absEta <= 2.7)
    etaCen  = (absEta <= 2.4)
    
    #forward jets
    outs[etaFor] &= ( (jet.neuHadFrac.content[etaFor] < 0.99) &
                      (jet.neuEmFrac.content[etaFor]  < 0.99) &
                      (jet.nParticles.content[etaFor] > 1   )  )
    #central jets
    outs[etaCen] &= ( (jet.chHadFrac.content[etaCen]  > 0.0 ) &
                      (jet.nCharged.content[etaCen]   > 0   ) &
                      (jet.chEmFrac.content[etaCen]   < 0.99 ) )
    #2.7-3.0
    etaHE = etaVFor & ~etaFor
    outs[etaHE] &= ( (jet.neuEmFrac.content[etaHE] > 0.01)  &
                     (jet.neuHadFrac.content[etaHE] < 0.98) &
                     (jet.nNeutrals.content[etaHE] > 2 )     )
    # > 3.0
    etaHF = ~etaVFor
    outs[etaHF] &= ( (jet.neuEmFrac.content[etaHF] > 0.90)  &
                     (jet.nNeutrals.content[etaHF] > 10 )     )
    
    outs = JaggedArray.fromoffsets(jet.pt.offsets,outs)
    return outs

def passJetTightLepVetoSel(jet):
    outs = np.ones_like(jet.pt.content,dtype=np.bool)
    absEta = np.abs(jet.eta.content)
    etaFor = (absEta <= 2.7)
    etaCen = (absEta <= 2.4)
    #forward jets
    outs[etaFor] &= ( (jet.neuHadFrac.content[etaFor] < 0.90) &
                      (jet.neuEmFrac.content[etaFor]  < 0.90) &
                      (jet.nParticles.content[etaFor] > 1   ) &
                      (jet.muonFrac.content[etaFor]   < 0.8 )  )
    #central jets
    outs[etaCen] &= ( (jet.chHadFrac.content[etaCen]  > 0.0 ) &
                      (jet.nCharged.content[etaCen]   > 0   ) &
                      (jet.chEmFrac.content[etaCen]   < 0.9 )  )
    outs = JaggedArray.fromoffsets(jet.pt.offsets,outs)
    return outs
