#btag working points

import numpy as np
from awkward import JaggedArray
from fnal_column_analysis_tools.analysis_objects import JaggedCandidateArray

CSVv2_L_AK4 = 0.5426 # CSVv2 WP
CSVv2_M_AK4 = 0.8484
CSVv2_T_AK4 = 0.9535

CSVv2SubJet_L_AK8 = 0.5426; # CSVv2SubJet WP
CSVv2SubJet_M_AK8 = 0.8484;

def CSV_Medium(jet,ptCut=50,absEtaCut=2.5):
    return ( (jet.pt > ptCut) & (np.abs(jet.eta) < absEtaCut) & (jet.csv > CSVv2_M_AK4) )
