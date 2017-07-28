# readlammpstrj.py
# Read lammps trajectory file. Save data in pandas dataframe.

import re
import pandas as pd

class readlammpstrj(object):
    def __init__(self, trjfilename):
        self.trjfn = trjfilename
        self.f = open(self.trjfn, 'r+')
    def readall(self):
        fbody  = self.f.readlines()
        indx = 0
        for line in fbody:
            fbody[indx] = line.split(' ')
            indx += 1
        print fbody

if __name__ == '__main__':
    go = readlammpstrj('./input/test.lammpstrj')
    go.readall()
    
