""" generator info manipulation for bacon_prod"""

from fnal_column_analysis_tools.util import numpy as np
from fnal_column_analysis_tools.util import awkward

#http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
activePdgIds = [1,2,3,4,5,6,21,11,13,15,22,23,24,-24,25,55,211,111]
activePdgIdLabels = list('duscbTgemtaZ') + ["W+","W-"] + list('hX')+['pi','pi0']
activePdgIds = {k:(1<<i) for i,k in enumerate(activePdgIds)}

def parseGeneratorHistory(gp_pdgId_in,gp_parent_in):
    inChain = activePdgIds
    #index manipulation
    offsets = gp_pdgId_in.offsets
    parents = gp_pdgId_in.parents
    pstarts = offsets[parents].astype('i4')
    gp_pdgId = gp_pdgId_in.content
    gp_ancestor = gp_parent_in.content + pstarts
    gp_ancestor_valid = (gp_parent_in.content >= 0)

    #create parentage bitmaps
    gp_pdgId_mapped = np.zeros(shape=gp_pdgId.shape, dtype='u4')
    for pdgId, bit in inChain.items():
        if abs(pdgId) == 24:
            gp_pdgId_mapped[gp_pdgId==pdgId] = bit
        else:
            gp_pdgId_mapped[np.abs(gp_pdgId)==pdgId] = bit
    
    gp_proc = np.zeros(shape=gp_pdgId.shape, dtype='u4')

    pdg_tmp = np.empty_like(gp_pdgId_mapped)
    parent_tmp = np.empty_like(gp_ancestor)
    niter = 0
    while np.any(gp_ancestor_valid) and niter < 50:
        np.take(gp_pdgId_mapped, gp_ancestor, out=pdg_tmp, mode='clip')
        np.take(gp_parent_in.content, gp_ancestor, out=parent_tmp, mode='clip')
        np.bitwise_or(gp_proc, pdg_tmp, where=gp_ancestor_valid, out=gp_proc)
        np.bitwise_and(gp_ancestor_valid, parent_tmp>=0, where=gp_ancestor_valid, out=gp_ancestor_valid)
        np.add(parent_tmp, pstarts, out=gp_ancestor)
        niter += 1

    #print 'Parsed ancestor tree in %d iterations'%niter

    if niter == 50 and np.any(gp_ancestor_valid):
        raise Exception('reached 50 iterations, gen particles not trustable')

    return awkward.JaggedArray.fromoffsets(offsets, gp_proc)

def tagDecayJagged(decayEndpoint,decayEndpointStatus,require,reject,parentage,status,pdgId,absPdg = True):
    pdgMatch = pdgId
    if absPdg: pdgMatch = np.abs(pdgId)
    return ((status == decayEndpointStatus) &
            (pdgMatch == decayEndpoint) &
            ((parentage & require) == require) &
            ((parentage & reject) == 0))

def tagDecay(decayEndpoint,decayEndpointStatus,require,reject,parentage,status,pdgId,absPdg = True):
    return tagDecayJagged(decayEndpoint,decayEndpointStatus,require,reject,parentage,status,pdgId,absPdg).any()

def isHadronicV(parentage,status,pdgId,statusValue=23):
    inChain = activePdgIds
    out = np.zeros_like(parentage.counts)
    
    #W decays
    W_plus = inChain[24]
    W_minus = inChain[-24]
    Top = inChain[6]
    Wlight_events = tagDecay(2,statusValue,W_plus,Top,parentage,status,pdgId)
    Wlight_events |= tagDecay(2,statusValue,W_minus,Top,parentage,status,pdgId)
    Wc_events = tagDecay(4,statusValue,W_plus,Top,parentage,status,pdgId)
    Wc_events |= tagDecay(4,statusValue,W_minus,Top,parentage,status,pdgId)
    
    #Z decays
    Z_mask = inChain[23]
    Zlight_events = tagDecay(1,statusValue,Z_mask,0,parentage,status,pdgId)
    Zlight_events |= tagDecay(2,statusValue,Z_mask,0,parentage,status,pdgId)
    Zlight_events |= tagDecay(3,statusValue,Z_mask,0,parentage,status,pdgId)
    Zc_events = tagDecay(4,statusValue,Z_mask,0,parentage,status,pdgId)
    Zb_events = tagDecay(5,statusValue,Z_mask,0,parentage,status,pdgId)
    
    #Higgs
    H_mask = inChain[25]
    Hc_events = tagDecay(4,statusValue,H_mask,0,parentage,status,pdgId)
    Hb_events = tagDecay(5,statusValue,H_mask,0,parentage,status,pdgId)
    
    #Exotica
    X_mask = inChain[55]
    Xc_events = tagDecay(4,statusValue,X_mask,0,parentage,status,pdgId)
    Xb_events = tagDecay(5,statusValue,X_mask,0,parentage,status,pdgId)
    
    #top specializations
    W_plus_top = W_plus | Top
    W_minus_top = W_minus | Top
    tW_events = tagDecay(2,statusValue,W_plus_top,0,parentage,status,pdgId)
    tW_events |= tagDecay(4,statusValue,W_plus_top,0,parentage,status,pdgId)
    tW_events |= tagDecay(2,statusValue,W_minus_top,0,parentage,status,pdgId)
    tW_events |= tagDecay(4,statusValue,W_minus_top,0,parentage,status,pdgId)
    
    out[Wlight_events | Zlight_events] = 1
    out[Wc_events | Zc_events | Hc_events | Xc_events] = 2
    out[Zb_events | Hb_events | Xb_events] = 3
    
    out[tW_events] = 9
    
    return out

def getHadronicVIndices(VPdgId,parentage,parents,status,pdgId,statusValue=23):
    inChain = activePdgIds
    
    mask = inChain[VPdgId]
    events = tagDecayJagged(1,statusValue,mask,0,parentage,status,pdgId).astype('u4')
    events = events + 2*tagDecayJagged(2,statusValue,mask,0,parentage,status,pdgId)
    events = events + 3*tagDecayJagged(3,statusValue,mask,0,parentage,status,pdgId)
    events = events + 4*tagDecayJagged(4,statusValue,mask,0,parentage,status,pdgId)
    events = events + 5*tagDecayJagged(5,statusValue,mask,0,parentage,status,pdgId)
    if abs(VPdgId) == 24:
        mask = inChain[-VPdgId]
        events = events + tagDecayJagged(1,statusValue,mask,0,parentage,status,pdgId)
        events = events + 2*tagDecayJagged(2,statusValue,mask,0,parentage,status,pdgId)
        events = events + 3*tagDecayJagged(3,statusValue,mask,0,parentage,status,pdgId)
        events = events + 4*tagDecayJagged(4,statusValue,mask,0,parentage,status,pdgId)
        events = events + 5*tagDecayJagged(5,statusValue,mask,0,parentage,status,pdgId)
    
    daughters_flat, = np.where(events.content != 0)
    dau_pdgIds = np.abs(pdgId.content[daughters_flat])
    
    #here we assume that V always decays into two things!!!!
    good_parents = parents[events != 0]
    good_dau_pdg = awkward.JaggedArray.fromcounts(np.full(daughters_flat.size//2,2),dau_pdgIds)
    good_dau_pdg = awkward.JaggedArray.fromcounts(good_parents.counts//2,good_dau_pdg)
    #good_daughters = awkward.JaggedArray.fromcounts(np.full(daughters_flat.size/2,2),daughters_flat)
    #good_daughters = awkward.JaggedArray.fromcounts(good_parents.counts/2,good_daughters)

    offset_Vs = (good_parents + parents.starts).content[::2]
    offset_Vs = awkward.JaggedArray.fromcounts(good_parents.counts//2,offset_Vs) - parents.starts
    
    return offset_Vs,good_dau_pdg.max() #return the up-type quark in each decay pair

def classifyTopDecays(argVs,Vdecays,eventTypes,parentage,status,parents,pdgId,statusValue=23):
    out = argVs.ones_like()
    out = out * eventTypes
    
    out.content[(out.content == 9) & (Vdecays.content == 4)] = 10
    
    return out

top_pdgId = 6
def getTopQuarks(pdgId,status): #this works for ttbar and single top...
    last_status_top = status[pdgId==top_pdgId].max()
    top_args = ((status == last_status_top) & (pdgId==top_pdgId)).astype(np.int).argmax()

    last_status_antitop = status[pdgId==-top_pdgId].max()
    antitop_args = ((status == last_status_antitop) & (pdgId==-top_pdgId)).astype(np.int).argmax()
    return top_args,antitop_args

#returns local index of first parent in history
# 0 if there is no parent of requested type
def getParentsOfType(parents,pdgs,theId,lastOnly=False):
    clipped_parents = np.maximum(parents,0)
    offset_parents = clipped_parents + parents.starts
    out = offset_parents.content.copy()
    getParentsOfTypeFlat(offset_parents.content,out,pdgs.content,theId,lastOnly)
    found = awkward.JaggedArray(parents.starts,parents.stops,out)
    return np.maximum(found-parents.starts,0)

def getParentsOfTypeFlat(parents_last,parents,pdgs,theId,out,lastOnly=False):
    #print 'step'
    parent_pdgs = pdgs[parents_last]
    notfound = ( (np.abs(parent_pdgs) != theId) & (parents_last != 0) )
    parents[notfound] = parents[parents[notfound]]
    if np.any(parents_last != parents):
        parents_last[notfound] = parents[notfound]
        getParentsOfTypeFlat(parents_last,parents,pdgs,theId,lastOnly)
    else:
        if lastOnly:
            parents[parent_pdgs == pdgs] = 0

def hasParentOfType(genPart,parents,pdgs,theId):
    pars = getParentsOfType(parents,pdgs,theId)
    return (pars[genPart] != 0)

def getChildrenOfType(parents,pdgs,parId,dauId):
    pars = getParentsOfType(parents,pdgs,parId)
    idxs = np.arange(parents.content.size)
    idxs = awkward.JaggedArray(parents.starts,parents.stops,idxs)
    idxs = idxs - parents.starts
    idxs = idxs[(pars != 0) & (np.abs(pdgs[idxs]) == dauId)]
    return idxs

def getHighestPtBoson(pdgs,pts,theId):
    idxs = np.arange(parents.content.size)
    idxs = awkward.JaggedArray(parents.starts,parents.stops,idxs)
    idxs = idxs - parents.starts
    Vidxs = idxs[np.abs(pdgs)==theId]
    Vpts = pts[np.abs(pdgs)==theId]
    return Vidxs[Vpts.argmax()]
