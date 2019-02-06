"""PUPPI soft-drop mass correction"""

from fnal_column_analysis_tools.util import numpy as np
from fnal_column_analysis_tools.util import awkward
from fnal_column_analysis_tools.analysis_objects import JaggedCandidateArray

def corrGEN(pt):
    #self.corrGEN = ROOT.TF1("corrGEN", "[0]+[1]*pow(x*[2],-[3])", 200, 3500)
    #self.corrGEN.SetParameter(0, 1.00626)
    #self.corrGEN.SetParameter(1, -1.06161)
    #self.corrGEN.SetParameter(2, 0.0799900)
    #self.corrGEN.SetParameter(3, 1.20454)
    x = np.clip(pt,200,3500)
    return 1.00626 - 1.06161*pow(x*0.0799900,-1.20454)

def corrRECO_cen(pt):
    x = np.clip(pt,200,3500)
    return 1.09302 + x*(-0.000150068+x*(3.44866e-07+x*(-2.68100e-10+x*(8.67440e-14+ x*(-1.00114e-17)))))

def corrRECO_for(pt):
    x = np.clip(pt,200,3500)
    return 1.27212 + x*(-0.000571640+x*(8.37289e-07+x*(-5.20433e-10+x*(1.45375e-13+ x*(-1.50389e-17)))))

def PUPPIweight(ak8pt,ak8eta):
    pt = ak8pt
    eta = ak8eta
    if isinstance(ak8pt,awkward.JaggedArray):
        assert (ak8pt.offsets==ak8eta.offsets).all()
        pt = ak8pt.flatten()
        eta = ak8eta.flatten()
    else:
        assert ak8pt.size == ak8eta.size
    
    genCorr = corrGEN(pt)
    recoCorr = np.ones_like(eta)
    
    etaCut = np.abs(eta) < 1.3
    recoCorr[etaCut] = corrRECO_cen(pt[etaCut])
    recoCorr[~etaCut] = corrRECO_for(pt[~etaCut])
    
    total = genCorr*recoCorr
    
    if isinstance(ak8pt,awkward.JaggedArray):
        total = awkward.JaggedArray.fromoffsets(ak8pt.offsets,total)
    
    return total
