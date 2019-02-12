""" kfactors / reweighting for ttbar """

from fnal_column_analysis_tools.util import numpy as np
from fnal_column_analysis_tools.util import awkward

def calculateTopKFactor(top_pts,antitop_pts):
    w1 = np.exp(0.0615 - 0.0005*top_pts);
    w2 = np.exp(0.0615 - 0.0005*antitop_pts);

    return np.sqrt(w1*w2).content
