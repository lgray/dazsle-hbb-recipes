"""jet vetos (run these after decorating the objects)"""

import numpy as np
from fnal_column_analysis_tools.util import awkward
from fnal_column_analysis_tools.analysis_objects import JaggedCandidateArray

from muons import passMuonLooseSel,passMuonTightSel
from electrons import passEleLooseSel,passEleTightSel
from photons import passPhoMediumSel
from taus import passTauSel

def selectVetoMuons(muons):
    commonVeto = ( (muons.pt > 10) &
                   (np.abs(muons.eta) < 2.4)  )
    looseMus = commonVeto & passMuonLooseSel(muons)
    tightMus = commonVeto & (muons.pt > 20) & passMuonTightSel(muons)
    
    return ( tightMus |
            (looseMus & (tightMus.sum() > 0) ) )

def selectVetoElectrons(electrons):
    commonVeto = ( (electrons.pt > 10) &
                   (electrons.absEta < 2.5) &
                   ((electrons.absEta < 1.4442) | (electrons.absEta > 1.566)) )
    looseEles = commonVeto & passEleLooseSel(electrons)
    tightEles = commonVeto & (electrons.pt > 40) & passEleTightSel(electrons)
    
    return ( tightEles |
            (looseEles & (tightEles.sum() > 0) ) )

def selectVetoPhotons(photons,electrons):
    matchesVetoElectron = photons.fastmatch(electrons,deltaRCut=0.4)
    
    return ( (photons.pt > 175) &
             (photons.absSCEta < 1.4442) &
             (passPhoMediumSel(photons)) &
             (~matchesVetoElectron) )

def selectVetoTaus(taus,muons,electrons):
    matchesVetoElectron = taus.fastmatch(electrons,deltaRCut=0.4)
    matchesVetoMuon = taus.fastmatch(muons,deltaRCut=0.4)
    
    return ( (~matchesVetoElectron) &
             (~matchesVetoMuon) &
             (taus.pt > 18) &
             (np.abs(taus.eta) < 2.3) &
             passTauSel(taus) )
