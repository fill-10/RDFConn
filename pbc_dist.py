import scipy
from scipy.linalg import norm
from scipy import array

def pbc_dist(coor_1, coor_2, boxlo, boxhi):
    dv = array([])
    boxedge = boxhi-boxlo
    for i in range(0, boxlo.size):
        dv = scipy.hstack( (dv, \
        array( [min(   abs(coor_1[i]-coor_2[i]), boxedge[i] - abs( coor_1[i]-coor_2[i]) )  ] )   ) )
    return norm(dv)
