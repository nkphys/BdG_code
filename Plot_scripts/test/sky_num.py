import numpy as np
pi = np.pi


lx = 21 ; ly = 21

data = np.loadtxt('MySkyrmion.txt', usecols=(2, 3, 4))
sz, sx, sy = data[:, 0], data[:, 1], data[:, 2]

chi = 0.0;

for iy in range (ly):
  for ix in range (lx):

    #chi1 = 0.0 ; chi2 = 0.0

    i = iy*lx + ix

    jy = iy ; jx = ix + 1                                 # nearest neighbour site along positive x axis
    if jx >= lx: jx = 0
    j1 = jy*lx + jx

    jx = ix ; jy = iy + 1                                 # nearest neighbour site along positive y axis
    if jy >= ly: jy = 0
    j2 = jy*lx + jx

    jy = iy ; jx = ix - 1                                 # nearest neighbour site along left diagonal
    if jx < 0: jx = lx-1
    j3 = jy*lx + jx

    jx = ix ; jy = iy - 1                                 # nearest neighbour site along negative x axis
    if jy < 0: jy = lx - 1
    j4 = jy*lx + jx

    chi1 = (1/(8*pi)) * ( (sx[i]*(sy[j1]*sz[j2] - sz[j1]*sy[j2])) + (sy[i]*(sz[j1]*sx[j2] - sx[j1]*sz[j2])) + (sz[i]*(sx[j1]*sy[j2] -sy[j1]*sx[j2])) )

    chi2 = (1/(8*pi)) * ( (sx[i]*(sy[j3]*sz[j4] - sz[j3]*sy[j4])) + (sy[i]*(sz[j3]*sx[j4] - sx[j3]*sz[j4])) + (sz[i]*(sx[j3]*sy[j4] - sy[j3]*sx[j4])) )

    chi = chi + chi1 + chi2

print(chi)
