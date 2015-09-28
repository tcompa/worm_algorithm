#!/usr/bin/pypy

'''
program: worm_ising_2d.py
created: 2015-09-28 -- 23 CEST
author: tc
notes: worm-algorithm simulation of the Ising model on the two-dimensional
       square lattice
'''

import random

# parameters
T = 2.5
L = 4
nsteps = 10 ** 8

N = L ** 2
K = 1.0 / T
inv_K = 1.0 / K

# define neighbors table
nbr = {i: ((i // L) * L + (i + 1) % L, (i + L) % N,
       (i // L) * L + (i - 1) % L, (i - L) % N) for i in xrange(N)}

# initialize bond weights
bonds = {}
for i in xrange(N):
    for j in nbr[i]:
        if tuple(sorted([i, j])) not in bonds.keys():
            bonds[tuple(sorted([i, j]))] = 0

Z = 0.0
Nb = 0
E_tot = 0.0
I = 0
M = 0

# main Monte Carlo loop
for step in xrange(nsteps):
    if I == M:
        I = random.randrange(N)
        M = I
        Z += 1.0
        E_tot += - T * Nb
    # start shift move
    new_M = random.choice(nbr[M])
    bond = tuple(sorted([M, new_M]))
    nb = bonds[bond]
    if random.uniform(0.0, 1.0) < 0.5:
        delta_nb = 1
        P_acc = K / (nb + 1.0)
    else:
        delta_nb = - 1
        P_acc = nb * inv_K
    if random.uniform(0.0, 1.0) < P_acc:
        bonds[bond] += delta_nb
        Nb += delta_nb
        M = new_M

print 'Z / nsteps = %f' % (Z / float(nsteps))
print '<E> / N = %f' % (E_tot / Z / float(N))
