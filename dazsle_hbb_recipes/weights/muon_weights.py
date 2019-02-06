"""muon eff scale factors for trig/iso/id"""

import scipy.stats
from fnal_column_analysis_tools.util import awkward
from fnal_column_analysis_tools.util import numpy as np
from fnal_column_analysis_tools.analysis_objects import JaggedCandidateArray
from fnal_column_analysis_tools.lookup_tools import evaluator
from copy import deepcopy

def MuonIDSF(evaluator,mupt,mueta):
    pt = mupt
    eta = np.abs(mueta)
    if isinstance(mupt,awkward.JaggedArray):
        assert (mupt.offsets==mueta.offsets).all()
        pt = mupt.flatten()
        eta = mueta.flatten()
    else:
        assert mupt.size == mueta.size
    eff = evaluator['NUM_SoftID_DEN_genTracks/abseta_pt_value'](eta,pt)
    err = evaluator['NUM_SoftID_DEN_genTracks/abseta_pt_error'](eta,pt)
    lo = eff - err
    hi = eff + err
    if isinstance(mupt,awkward.JaggedArray):
        eff = awkward.JaggedArray.fromoffsets(mupt.offsets,eff)
        lo  = awkward.JaggedArray.fromoffsets(mupt.offsets,lo)
        hi  = awkward.JaggedArray.fromoffsets(mupt.offsets,hi)
    return eff,hi,lo

def MuonIsoSF(evaluator,mupt,mueta):
    pt = mupt
    eta = np.abs(mueta)
    if isinstance(mupt,awkward.JaggedArray):
        assert (mupt.offsets==mueta.offsets).all()
        pt = mupt.flatten()
        eta = mueta.flatten()
    else:
        assert mupt.size == mueta.size
    eff = evaluator['NUM_LooseRelIso_DEN_LooseID/abseta_pt_value'](eta,pt)
    err = evaluator['NUM_LooseRelIso_DEN_LooseID/abseta_pt_error'](eta,pt)
    lo = eff - err
    hi = eff + err
    if isinstance(mupt,awkward.JaggedArray):
        eff = awkward.JaggedArray.fromoffsets(mupt.offsets,eff)
        lo  = awkward.JaggedArray.fromoffsets(mupt.offsets,lo)
        hi  = awkward.JaggedArray.fromoffsets(mupt.offsets,hi)
    return eff,hi,lo

def MuonTrigSF(evaluator,mupt,mueta):
    pt = mupt
    eta = np.abs(mueta)
    if isinstance(mupt,awkward.JaggedArray):
        assert (mupt.offsets==mueta.offsets).all()
        pt = mupt.flatten()
        eta = mueta.flatten()
    else:
        assert mupt.size == mueta.size
    eff = evaluator['Mu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA'](pt,eta)
    err = evaluator['Mu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA_error'](pt,eta)
    lo = eff - err
    hi = eff + err
    if isinstance(mupt,awkward.JaggedArray):
        eff = awkward.JaggedArray.fromoffsets(mupt.offsets,eff)
        lo  = awkward.JaggedArray.fromoffsets(mupt.offsets,lo)
        hi  = awkward.JaggedArray.fromoffsets(mupt.offsets,hi)
    return eff,hi,lo
