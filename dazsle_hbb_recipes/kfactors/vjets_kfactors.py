"""WJetsQQ and ZJetsQQ kfactors"""

import numpy as np
from fnal_column_analysis_tools.lookup_tools import evaluator
from awkward import JaggedArray
from copy import deepcopy

#hack for 2016
def hackEvaluatorForVJetsQQ_2016(lookup):
    wscale=np.array([1.0,1.0,1.0,1.20,1.25,1.25,1.0])
    ptscale=np.array([0, 500, 600, 700, 800, 900, 1000,3000])
    zqq = deepcopy(lookup['ZJetsNLO'])
    wqq = deepcopy(lookup['WJetsNLO'])
    
    wqq._values = 1.35*wscale
    wqq._axes = ptscale
    
    zqq._values = np.array([1.45])
    zqq._axes = np.array([0,3000])
    
    lookup._functions['WJetsNLO_2016'] = wqq
    lookup._functions['ZJetsNLO_2016'] = zqq

#2016
def VJetsQQ_kFactor2016(sampleName,lookup,genVPt):
    kfactor = genVPt.ones_like()
    if 'ZJetsToQQ_' in sampleName:
        kfactor._content = lookup['ZJetsNLO_2016'](np.clip(genVPt.content,250.,1200.))
    elif 'WJetsToQQ_' in sampleName:
        kfactor._content = lookup['WJetsNLO_2016'](np.clip(genVPt.content,250.,1200.))
    return kfactor

#2017
def VJetsQQ_kFactor2017(sampleName,lookup,genVPt):
    kfactor = genVPt.ones_like()
    if 'ZJetsToQQ_' in sampleName:
        kfactor._content = lookup['ZJetsNLO'](np.clip(genVPt.content,250.,1200.))
    elif 'WJetsToQQ_' in sampleName:
        kfactor._content = lookup['WJetsNLO'](np.clip(genVPt.content,250.,1200.))
    return kfactor

def calculateBaseKFactor(sampleName,lookup,genVPt):
    kfactor = genVPt.ones_like()
    maxxed = np.maximum(genVPt,100.)
    
    if ( 'ZJets' in sampleName or  'DYJets' in sampleName or
        'ZPrime' in sampleName or 'VectorDiJet' in sampleName ):
        kfactor = kfactor * lookup["EWKcorr/Z"](maxxed)/lookup["ZJets_LO/inv_pt"](maxxed)
    if 'WJets' in sampleName:
        kfactor = kfactor * lookup["EWKcorr/W"](maxxed)/lookup["WJets_LO/inv_pt"](maxxed)
    return kfactor

def calculateNLOKFactorAndSysts(sampleName,lookup,genVPt):
    central = genVPt.ones_like()
    pdf_up,pdf_down = genVPt.ones_like(),genVPt.ones_like()
    ren_up,ren_down = genVPt.ones_like(),genVPt.ones_like()
    fac_up,fac_down = genVPt.ones_like(),genVPt.ones_like()
    
    clipped = np.clip(genVPt.content,100.,700.)
    
    if ( 'ZJets' in sampleName or  'DYJets' in sampleName or
         'ZPrime' in sampleName or 'VectorDiJet' in sampleName ):
        central._content   = lookup["ZJets_012j_NLO/nominal"](clipped)
        pdf_raw            = lookup["ZJets_012j_NLO/PDF"](clipped)
        pdf_up._content    = pdf_up._content   + pdf_raw
        pdf_down._content  = pdf_down._content - pdf_raw
        ren_up._content    = lookup["ZJets_012j_NLO/ren_up"](clipped)
        ren_down._content  = lookup["ZJets_012j_NLO/ren_down"](clipped)
        fac_up._content    = lookup["ZJets_012j_NLO/fact_up"](clipped)
        fac_down._content  = lookup["ZJets_012j_NLO/fact_down"](clipped)
    if 'WJets' in sampleName:
        central._content   = lookup["WJets_012j_NLO/nominal"](clipped)
        pdf_raw            = lookup["WJets_012j_NLO/PDF"](clipped)
        pdf_up._content    = pdf_up._content + pdf_raw
        pdf_down._content  = pdf_down._content - pdf_raw
        ren_up._content    = lookup["WJets_012j_NLO/ren_up"](clipped)
        ren_down._content  = lookup["WJets_012j_NLO/ren_down"](clipped)
        fac_up._content    = lookup["WJets_012j_NLO/fact_up"](clipped)
        fac_down._content  = lookup["WJets_012j_NLO/fact_down"](clipped)
        
    return (central,{'PDF':(pdf_up,pdf_down),'ren':(ren_up,ren_down),'fact':(fac_up,fac_down)})

def getKFactor2016(sampleName,lookup,genVPt):
    kfactor   = genVPt.ones_like()
    kfactor   = kfactor * calculateBaseKFactor(sampleName,lookup,genVPt)
    kfactor   = kfactor * VJetsQQ_kFactor2016(sampleName,lookup,genVPt)
    nlo_systs = calculateNLOKFactorAndSysts(sampleName,lookup,genVPt)[1]
    return {'central':kfactor,'systs':nlo_systs}

def getKFactor2017(sampleName,lookup,genVPt):
    kfactor   = genVPt.ones_like()
    kfactor   = kfactor * calculateBaseKFactor(sampleName,lookup,genVPt)
    kfactor   = kfactor * VJetsQQ_kFactor2017(sampleName,lookup,genVPt)
    nlo_systs = calculateNLOKFactorAndSysts(sampleName,lookup,genVPt)[1]
    return {'central':kfactor,'systs':nlo_systs}
