""" generator info manipulation for bacon_prod"""

import numpy as np
from awkward import JaggedArray

#returns local index of first parent in history
# 0 if there is no parent of requested type
def getParentsOfType(parents,pdgs,theId,lastOnly=False):
    offset_parents = parents + parents.starts
    found = JaggedArray(parents.starts,parents.stops,
                        getParentsOfTypeFlat(offset_parents.content,pdgs.content,theId,lastOnly))
    return np.maximum(found-parents.starts,0)

def getParentsOfTypeFlat(parents,pdgs,theId,lastOnly=False):
    clipped_parents = np.maximum(parents,0)
    parent_pdgs = pdgs[clipped_parents]
    notfound = (np.abs(parent_pdgs) != theId)
    newparents = clipped_parents.copy()
    newparents[notfound] = newparents[clipped_parents[notfound]]
    if (newparents != clipped_parents).any():
        return getParentsOfTypeFlat(newparents,pdgs,theId,lastOnly)
    else:
        if lastOnly:
            parents[parent_pdgs == pdgs] = 0
        return parents

def hasParentOfType(genPart,parents,pdgs,theId):
    pars = getParentsOfType(parents,pdgs,theId)
    return (pars[genPart] != 0)

def getChildrenOfType(parents,pdgs,parId,dauId):
    pars = getParentsOfType(parents,pdgs,parId)
    idxs = np.arange(parents.content.size)
    idxs = JaggedArray(parents.starts,parents.stops,idxs)
    idxs = idxs - parents.starts
    idxs = idxs[(pars != 0) & (np.abs(pdgs[idxs]) == dauId)]
    return idxs
