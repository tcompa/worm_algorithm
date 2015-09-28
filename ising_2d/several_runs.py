'''
created: 2015-09-28 -- 23 CEST
author: tc
'''

import numpy
import itertools
import subprocess
import sys
import os


def prepare_input(L, T, nsteps):
    f = open('input.dat', 'w')
    f.write('%i %.12f %i\n' % (L, T, nsteps))
    f.close()

list_T = numpy.arange(0.5, 4.5, 0.5)
L = 6
nsteps = 10 ** 8 * 5

print 'compilation -- start'
os.system('make clean')
os.system('make worm_ising_2d')
print 'compilation -- done'
prog = './worm_ising_2d.x'
if not os.path.isfile(prog):
    sys.exit('ERROR')

print 'runs -- start'
for T in list_T:
    prepare_input(L, T, nsteps)
    p = subprocess.Popen([prog])
    p.wait()
    L, T, K, Z_by_nsteps, e_av = numpy.loadtxt('output.dat')
    print '%4i %.4f %.4f %.4f %.8f' % (L, T, K, Z_by_nsteps, e_av)
os.remove('input.dat')
os.remove('output.dat')
