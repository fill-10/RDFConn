# rdf class
import os
from scipy import *


import scipy.linalg
from scipy.linalg import norm

from pbc_dist import pbc_dist
from loadlammpstrj import loadlammpstrj

from class_bead import bead
from class_bond import bond

class rdf:
    def __init__(self, upboundary = 10.0, increasingstep = 1.0, atom_1 = '2', atom_2 = '5'):
        if upboundary>0 and increasingstep >0 and  upboundary/increasingstep>=1:
            self.upbd = upboundary
            self.step = increasingstep
            self.bins = []
            self.boxlo = []
            self.boxhi = []
            for i in range(0, int(self.upbd / self.step) ):
                self.bins.append([  (  (i+0.5) * self.step )  ])
                self.bins[i].append(0) #RDF value
            #print(self.bins)
            self.at_1 = atom_1
            self.at_2 = atom_2

        else:
            print('rdf N/A')


    
    def loadfile(self, trjfile):
        self.f = open(trjfile, 'r+')
        [ self.timestep, self.Natom, self.boxhi, self.boxlo, self.data_start_pt] = loadlammpstrj(self.f)
        # The 'trjfile' in fact is a pointer.
        # So it passes the pointer into the function.
        # Every step in the function loadlammpstrj operates on the pointer 'self.f'.
        # That is why no need to run self.f.seek(self.data_start_pt) in the following functions.
        # 


    def genlist(self,atomname):
        list_1 = []
        spos = self.f.tell()
        for i in range(0, self.Natom):
            cline = self.f.readline().split()
            if cline[2] == atomname:
                coor_1 = array(  [  float(cline[3]), float(cline[4]), float(cline[5])  ]  )
                b1 = bead( coor_1, cline[2], int(cline[0]), float(cline[1]) )
                list_1.append(b1)
        self.f.seek(spos)
        return list_1

    def gen2list(self):
        self.list_1 = self.genlist(self.at_1)
        if self.at_1 != self.at_2:
            self.list_2 = self.genlist(self.at_2)
        else:
            self.list_2 = None
        #
        # print section
        print('Number of atoms in list1:   '+   str(len(self.list_1) )   )
        if self.list_2:
            print('Number of atoms in list2:   '+   str(len(self.list_2)  )  )
        else:
            print('No list2')
        # close file
        self.f.close()
        del self.f

    def sortpairs(self):
        self.gen2list()
        if self.list_2:
            for b1 in self.list_1:
                for b2 in self.list_2:
                    pair = bond(b2, b1)
                    pair.blength = pbc_dist(b2.coor, b1.coor , self.boxlo, self.boxhi)
                    for kbin in self.bins:
                        # The criteron might need to be reconsidered.
                        # In this code, it is a [r_in , r_out ) domain.
                        # The first bin includes the distance=0 case.
                        # This might be a problem.
                        # Fortunately, it should be impossible for 2 atoms 
                        # having exactly the same coordinates.
                        if pair.blength >= ( kbin[0] - self.step/2 ) and pair.blength < ( kbin[0] + self.step/2) :
                            kbin[1] = kbin[1] + 1
        else:
            for i in range(0, len(self.list_1)):
                # Notice that j is from i+1!!
                # I made the mistake that j went from i,
                # which made an error for the first bin!
                for j in range(i+1, len(self.list_1)):
                    pair = bond(self.list_1[i], self.list_1[j])
                    pair.blength = pbc_dist(self.list_1[i].coor, self.list_1[j].coor , self.boxlo, self.boxhi)
                    for kbin in self.bins:
                        if pair.blength >= ( kbin[0] - self.step/2 ) and pair.blength < ( kbin[0] + self.step/2) :
                            kbin[1] = kbin[1] +1

    def calrdf(self):
        self.result = []
        #
        #calculate normal factor
        if self.at_1 != self.at_2 :
            normfactor =  (self.boxhi[0] - self.boxlo[0])\
            * (self.boxhi[1] - self.boxlo[1]) \
            * (self.boxhi[2] - self.boxlo[2]) \
            / len(self.list_1) / len(self.list_2)
        else:
            
            normfactor =  (self.boxhi[0] - self.boxlo[0])\
            * (self.boxhi[1] - self.boxlo[1]) \
            * (self.boxhi[2] - self.boxlo[2]) \
            /  len(self.list_1) / (len(self.list_1)-1) *2
        print 'normal factor is :'
        print normfactor
        print '\n'


        for kbin in self.bins:
            Npairbin = kbin[1]
            #
            # Print the number of atoms in each bin.
            # This is to help us understand the data.
            print Npairbin
            #
            bindensity = Npairbin /  \
            (  4.0/3.0 * 3.14159 \
            *  ( (kbin[0] + self.step /2 )**3 - (kbin[0] - self.step /2 )**3 ) \
            )

            gofrdata = bindensity * normfactor

            self.result.append([kbin[0], gofrdata])

        del self.bins
        del self.list_1
        del self.list_2

        return self.result
    
    def savetofile(self, outfilename):
        out =  open(outfilename, 'w+')
        for i in range(0, len(self.result)):
            out.write(str(self.result[i][0])+'    '+str(self.result[i][1])+'\n' )
        out.close()
