import numpy as np

from diff_cython import *

def diff_py(at, a, visc, dxidxi, dyidyi, dzidzi, itot, jtot, ktot): 
    ii = 1
    jj = itot
    kk = itot*jtot

    for k in range(1, ktot-1):
        for j in range(1, jtot-1):
            for i in range(1, itot-1):
                ijk = i + j*jj + k*kk
                at[ijk] += visc * (
                        + ( (a[ijk+ii] - a[ijk   ]) 
                          - (a[ijk   ] - a[ijk-ii]) ) * dxidxi 
                        + ( (a[ijk+jj] - a[ijk   ]) 
                          - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
                        + ( (a[ijk+kk] - a[ijk   ]) 
                          - (a[ijk   ] - a[ijk-kk]) ) * dzidzi
                        )

def init(a, at, ncells):
    for i in range(ncells):
        a[i]  = i*i
        at[i] = 0.

if(__name__ == "__main__"):
    nloop = 100
    itot = 256
    jtot = 256
    ktot = 256
    ncells = itot*jtot*ktot
    
    a  = np.zeros(ncells)
    at = np.zeros(ncells)
    
    init(a, at, ncells)
   
    for i in range(nloop): 
        #diff_py(at, a, 0.1, 0.1, 0.1, 0.1, itot, jtot, ktot) 
        diff_pyx(at, a, 0.1, 0.1, 0.1, 0.1, itot, jtot, ktot) 
   
    #print("at=%f"%at[itot*jtot+2*itot+2])
