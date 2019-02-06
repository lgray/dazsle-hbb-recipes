"""tau IDs"""

import numpy as np
from fnal_column_analysis_tools.util import awkward
from fnal_column_analysis_tools.analysis_objects import JaggedCandidateArray

# descriptions from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Discriminators
kByDecayModeFinding                        = 1 << 0   # You will always want to use this (see AN-10-82)

kByVLooseIsolation                         = 1 << 1   # isolation cone of 0.3, no PF Charged Candidates with pT > 1.5 GeV/c and no PF Gamma candidates with ET > 2.0 GeV
kByLooseIsolation                          = 1 << 2   # (description N/A)
kByMediumIsolation                         = 1 << 3   # (description N/A)
kByTightIsolation                          = 1 << 4   # (description N/A)
kByVLooseIsolationDBSumPtCorr              = 1 << 5   # (description N/A)
kByLooseIsolationDBSumPtCorr               = 1 << 6   # (description N/A)
kByMediumIsolationDBSumPtCorr              = 1 << 7   # (description N/A)
kByTightIsolationDBSumPtCorr               = 1 << 8   # (description N/A)
kByVLooseCombinedIsolationDBSumPtCorr      = 1 << 9   # isolation cone of 0.3, Delta Beta corrected sum pT of PF charged and PF gamma isolation candidates (pT > 0.5 GeV) less than 3 GeV
kByLooseCombinedIsolationDBSumPtCorr       = 1 << 10  # isolation cone of 0.5, Delta Beta corrected sum pT of PF charged and PF gamma isolation candidates (pT > 0.5 GeV) less than 2 GeV
kByMediumCombinedIsolationDBSumPtCorr      = 1 << 11  # isolation cone of 0.5, Delta Beta corrected sum pT of PF charged and PF gamma isolation candidates (pT > 0.5 GeV) less than 1 GeV
kByTightCombinedIsolationDBSumPtCorr       = 1 << 12  # isolation cone of 0.5, Delta Beta corrected sum pT of PF charged and PF gamma isolation candidates (pT > 0.5 GeV) less than 0.8 GeV
kByLooseCombinedIsolationDBSumPtCorr3Hits  = 1 << 13  # same as ByLooseCombinedIsolationDBSumPtCorr but requiring 3 hits (instead of 8) on track of isolation candidates
kByMediumCombinedIsolationDBSumPtCorr3Hits = 1 << 14  # same as ByMediumCombinedIsolationDBSumPtCorr but requiring 3 hits (instead of 8) on track of isolation candidates
kByTightCombinedIsolationDBSumPtCorr3Hits  = 1 << 15  # same as ByTightCombinedIsolationDBSumPtCorr but requiring 3 hits (instead of 8) on track of isolation candidates
kByLooseIsolationMVA                       = 1 << 16  # BDT based selection using isolation in rings around tau direction and shower shape variables
kByMediumIsolationMVA                      = 1 << 17  #  "
kByTightIsolationMVA                       = 1 << 18  #  "
kByLooseIsolationMVA2                      = 1 << 19  # same as ByLooseIsolationMVA with new training and improved performance
kByMediumIsolationMVA2                     = 1 << 20  # same as ByMediumIsolationMVA with new training and improved performance
kByTightIsolationMVA2                      = 1 << 21  # same as ByTightIsolationMVA with new training and improved performance

kByLooseElectronRejection                  = 1 << 22  # electron pion MVA discriminator < 0.6
kByMediumElectronRejection                 = 1 << 23  # electron pion MVA discriminator < -0.1 and not 1.4442 < |eta| < 1.566
kByTightElectronRejection                  = 1 << 24  # electron pion MVA discriminator < -0.1 and not 1.4442 < |eta| < 1.566 and Brem pattern cuts (see AN-10-387)
kByMVA3LooseElectronRejection              = 1 << 25  # anti-electron MVA discriminator with improved training
kByMVA3MediumElectronRejection             = 1 << 26  #  "
kByMVA3TightElectronRejection              = 1 << 27  #  "
kByMVA3VTightElectronRejection             = 1 << 28  #  "

kByLooseMuonRejection                      = 1 << 29  # Tau Lead Track not matched to chamber hits
kByMediumMuonRejection                     = 1 << 30  # Tau Lead Track not matched to global/tracker muon
kByTightMuonRejection                      = 1 << 31  # Tau Lead Track not matched to global/tracker muon and large enough energy deposit in ECAL+HCAL exceeding 0.2 times Lead Track momentum
kByLooseMuonRejection2                     = 1 << 32  # Same as LooseMuonRejection
kByMediumMuonRejection2                    = 1 << 33  # LooseMuonRejection2 && no DT, CSC or RPC Hits in last 2 Stations
kByTightMuonRejection2                     = 1 << 34  # MediumMuonRejection2 && large enough energy deposit in ECAL+HCAL in 1 prong + 0 strip decay mode (SUM(ECAL+HCAL)>0.2*pT)
kByLooseMuonRejection3                     = 1 << 35  # Same as LooseMuonRejection
kByTightMuonRejection3                     = 1 << 36  # MediumMuonRejection2 && large enough energy deposit in ECAL+HCAL in 1 prong + 0 strip decay mode (SUM(ECAL+HCAL)>0.2*pT)

def passTauSel(tau):
    return ( ((tau.hpsDisc & kByDecayModeFinding) > 0) & (tau.rawIso3Hits < 5) )
