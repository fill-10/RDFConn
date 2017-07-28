# this function, or say, class, is to find out the connectivity(ies).
#

from scipy import *
from scipy.linalg import norm

from class_bead import bead

from pbc_dist import pbc_dist
from loadlammpstrj import loadlammpstrj

class connectivity:
    def __init__(self, atomname, clustercutoff):
        if clustercutoff>0:
            self.cutoff = clustercutoff
            self.at_1 = atomname
            self.clusters=[]
        else:
            print('cutoff Neg')
    def loadfile(self,trjfile):
        self.f = open(trjfile, 'r')
        [ self.timestep, self.Natom, self.boxhi, self.boxlo, self.data_start_pt] = loadlammpstrj(self.f)
        # The 'trjfile' in fact is a pointer.
        # So it passes the pointer into the function.
        # Every step in the function loadlammpstrj operates on the pointer 'self.f'.
    def genlist(self,atomname):
        list_1 = []
        for i in range(0, self.Natom):
            cline = self.f.readline().split()
            if cline[2] == atomname:
                coor_1 = array(  [  float(cline[3]), float(cline[4]), float(cline[5])  ]  )
                b1 = bead( coor_1, cline[2], int(cline[0]), float(cline[1]) )
                list_1.append(b1)
        return list_1
    
    def findcluster(self):
        self.alist = self.genlist(self.at_1)
        self.NWatom = len(self.alist)
        print 'Selected atoms: '
        print self.NWatom
        print '\n'
        print 'Atoms left: \n'
        while len(self.alist)> 1 :
            # create a current list to store the cluster found in a single step
            clist= [ self.alist[0] ]
            #print clist[0].coor
            del self.alist[0]
            for iatom in clist :
                j = 0
                while j<len(self.alist):
                    distance = pbc_dist(iatom.coor, self.alist[j].coor, self.boxlo, self.boxhi)
                    if distance <= self.cutoff:
                        clist.append(self.alist[j])
                        del self.alist[j]
                    else: j += 1
            if len(clist)>1:
                self.clusters.append(clist)
            print len(self.alist)



           
    def atoms2trjfile(self, outfilename):
        if self.clusters != [] :
            out = open(outfilename , 'w+')
            for cluster_1 in self.clusters:
                pass

    def caldistro(self, outfilename = ''):
        if self.clusters != []:
            # count size
            clust_size = []
            for ccluster in self.clusters:
                clust_size.append(len(ccluster))
            self.clust_distro = []
            
            Nmono = self.NWatom
            
            for i in range(2,max(clust_size)+1):
                self.clust_distro.append( [i, clust_size.count(i)])
                # add 1 atom 'cluster' first, optional.
                Nmono -= self.clust_distro[i-2][1] * i
            
            self.clust_distro.insert(0, [1,Nmono ])
                


            if outfilename :
                fout = open(outfilename, 'w')
                for j in self.clust_distro:
                    fout.write('%8d' %j[0] + 4*' ' + '%8d' %j[1] +'\n' )
                fout.close()                
