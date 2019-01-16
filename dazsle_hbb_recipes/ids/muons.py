"""muon IDs"""

import numpy as np
from awkward import JaggedArray
from fnal_column_analysis_tools.analysis_objects import JaggedCandidateArray

kPOGLooseMuon  =  1 << 0
kPOGMediumMuon =  1 << 1
kPOGTightMuon  =  1 << 2
kPOGSoftMuon   =  1 << 3
kPOGHighPtMuon =  1 << 4

def calcMuonDeltaBetaIso(muon):
    return muon.chHadIso + np.maximum(muon.neuHadIso + muon.gammaIso - 0.5*(muon.puIso), 0.0)

def passMuonLooseSel(muon):
    return (((muon.pogIDBits & kPOGLooseMuon) > 0) & (muon.deltaBetaIso < 0.25*muon.pt))

def passMuonMediumSel(muon):
    return (((muon.pogIDBits & kPOGMediumMuon) > 0) & (muon.deltaBetaIso < 0.25*muon.pt))

def passMuonTightSel(muon):
    return (((muon.pogIDBits & kPOGTightMuon) > 0) & (muon.deltaBetaIso < 0.25*muon.pt))

def passHighPtMuonSel(muon):
    return (((muon.pogIDBits & kPOGHighPtMuon) > 0) & (muon.deltaBetaIso < 0.25*muon.pt))
