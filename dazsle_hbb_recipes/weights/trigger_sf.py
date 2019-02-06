"""trigger eff scale factors (todo: )"""

from fnal_column_analysis_tools.util import awkward
from fnal_column_analysis_tools.util import numpy as np
import scipy.stats
from fnal_column_analysis_tools.analysis_objects import JaggedCandidateArray
from fnal_column_analysis_tools.lookup_tools import evaluator
from copy import deepcopy

def hackEvaluatorForVTrigSF(evaluator,alpha=0.31731):
    num = evaluator['data_obs_muCR4_numerator']
    den = evaluator['data_obs_muCR4_denominator']
    eff = deepcopy(num)
    lows = deepcopy(num)
    highs = deepcopy(num)
    
    eff._values = num._values/den._values
    #clopper-pearson interval
    lows._values = scipy.stats.beta.ppf(alpha/2, num._values, den._values-num._values+1)
    highs._values = scipy.stats.beta.ppf(1 - alpha/2, num._values+1, den._values-num._values)
    
    evaluator._functions['data_obs_muCR4_eff'] = eff
    evaluator._functions['data_obs_muCR4_eff_down'] = lows
    evaluator._functions['data_obs_muCR4_eff_up'] = highs

def VtrigSF(evaluator,ak8msd,ak8pt):
    msd = ak8msd
    pt = ak8pt
    if isinstance(ak8msd,awkward.JaggedArray):
        assert (ak8msd.offsets==ak8pt.offsets).all()
        msd = ak8msd.flatten()
        pt = ak8pt.flatten()
    else:
        assert ak8msd.size == ak8pt.size
    msd = np.clip(msd,0,300)
    pt = np.clip(pt,200,1000)
    eff = evaluator['data_obs_muCR4_eff'](msd,pt)
    lo  = evaluator['data_obs_muCR4_eff_down'](msd,pt)
    hi  = evaluator['data_obs_muCR4_eff_up'](msd,pt)
    if isinstance(ak8msd,awkward.JaggedArray):
        eff = awkward.JaggedArray.fromoffsets(ak8pt.offsets,eff)
        lo  = awkward.JaggedArray.fromoffsets(ak8pt.offsets,lo)
        hi  = awkward.JaggedArray.fromoffsets(ak8pt.offsets,hi)
    return eff,hi,lo
