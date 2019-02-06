""" generator info manipulation for bacon_prod"""

from fnal_column_analysis_tools.util import numpy as np
from fnal_column_analysis_tools.util import awkward

activePdgIds = [6,11,13,15,22,23,24,25,55]# 1,2,3,4,5,21,211,111
activePdgIdLabels = list('TemtaZWhX')#+['pi','pi0']#udscbg

activePdgIds = {k:(1<<i) for i,k in enumerate(activePdgIds)}

def parseGeneratorHistory(gp_pdgId_in,gp_parent_in):
    #index manipulation
    offsets = gp_pdgId_in.offsets
    parents = gp_pdgId_in.parents
    pstarts = offsets[parents].astype('i4')
    gp_pdgId = gp_pdgId_in.content
    gp_ancestor = gp_parent_in.content + pstarts
    gp_ancestor_valid = (gp_parent_in.content >= 0)

    #create parentage bitmaps
    gp_pdgId_mapped = np.zeros(shape=gp_pdgId.shape, dtype='u4')
    for pdgId, bit in activePdgIds.items():
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

    print 'Parsed ancestor tree in %d iterations'%niter

    if niter == 50 and np.any(gp_ancestor_valid):
        raise Exception('reached 50 iterations, gen particles not trustable')

    return awkward.JaggedArray.fromoffsets(offsets, gp_proc)

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
