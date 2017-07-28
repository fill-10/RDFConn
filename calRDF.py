import os
import scipy
from scipy import array
from scipy import matrix

import scipy.linalg
from scipy.linalg import norm

from pbc_dist import pbc_dist

from class_bead import bead
from class_bond import bond
from class_rdf import rdf


gofr = rdf(10.0, 0.5, '7', '7')

#gofr.loadfile('tesxt.txt')
gofr.loadfile('step72M.lammpstrj')
# don't run this
# it's a huge file ...


gofr.sortpairs()

gofr.calrdf()

gofr.savetofile('rdf.dat')


