# This function is to load the lammps trajectory file
# This function can be used in any analysis needing to load trj file.

from scipy import *

def loadlammpstrj(fpointer):
    fpointer.seek(0)
    fpointer.readline()
    timestep = int(fpointer.readline().split()[0])
    fpointer.readline()
    Natom = int(fpointer.readline().split()[0])
    fpointer.readline()
    boxlo = array([])
    boxhi = array([])
    for i in range(0,3):
        cline = fpointer.readline().split()
        boxlo = hstack(  ( boxlo, array(  [  float(cline[0])  ]  ) )   )
        boxhi = hstack(  ( boxhi, array(  [  float(cline[1])  ]  ) )   )
    fpointer.readline()
    data_start_pt = fpointer.tell()
    out_list  = [timestep, Natom, boxhi, boxlo, data_start_pt]
    return out_list

    # Bug: this algorithm cannot check if the upboundary blows out of the box
    # some other bugs such as negative values of XXX... so, use this carefully


