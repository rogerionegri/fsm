
#--------------------------------------------
#IMPORTS ------------------------------------
import numpy as np
import re


#==================================================================
def delta_nbr_weight(vec):
    vec_out = 1.0 - np.exp( np.abs(vec[:]) )
    return vec_out

#==================================================================
def nbr_weight(vec):
    vec_out = 1/(1+ np.exp(-1*(vec[:])))
    return vec_out

#==================================================================
def ndvi_map(vec):
    vec_out = np.exp(-1*np.abs(vec[:]))
    return vec_out

#==================================================================
def nbr_map(vec):
    vec_out = np.exp(-1*(1+vec[:]))
    return vec_out

#==================================================================
def select_attributes_from_list(list,bands):
    listSel = []
    for att in list:  #Lista completa com Atts e DeltaNBR
        exp = re.compile(str(att)+r'_\d')
        for item in bands:
            out = exp.match(item)
            if out != None: listSel.append(item)
    return listSel

#==================================================================
def get_attribute_index_from_list(att,bands):
    listSel = []
    listIndex = []
    exp = re.compile(att+r'_\d')
    for item, index in zip(bands,range(len(bands))):
        out = exp.match(item)
        if out != None:
            listSel.append(item)
            listIndex.append(index)
    return listSel, listIndex