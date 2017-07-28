# borro the concept of bond to store
# pair info

from pbc_dist import pbc_dist

class bond :
    def __init__(self, bead1, bead2, bondtype='0'):
        self.bd1 = bead1
        self.bd2 = bead2
        self.btp = bondtype
        self.bvector = self.bd2.coor - self.bd1.coor
        
        #self.blength = pbc_dist(self.bd1.coor, self.bd2.coor, [0, 0, 0] , [0, 0, 0] )
        self.blength = None
        # bondlength is commented out!
