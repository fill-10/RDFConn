# bead.py
# bead

from scipy import array

class bead:
	def __init__(self, location=[0.,0.,0.], beadtype ='0', atomID =0 , molID = 0, clusterID = 0):
		self.coor = array(location)
		self.type = beadtype
                self.atomid = atomID
                self.molid = molID
                self.clid = clusterID
