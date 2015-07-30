import numpy as np
import time

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
        a[i]  = i**2./(i+1)**2.
        at[i] = 0.

if(__name__ == "__main__"):
    mode = 'cy'
    nloop = 100
    itot = 128
    jtot = 128
    ktot = 128
    ncells = itot*jtot*ktot
    
    a  = np.zeros(ncells)
    at = np.zeros(ncells)
    
    init(a, at, ncells)

    # Check results
    if(mode == 'py'):
        diff_py(at, a, 0.1, 0.1, 0.1, 0.1, itot, jtot, ktot) 
    elif(mode == 'cy'):
        diff_pyx(at, a, 0.1, 0.1, 0.1, 0.1, itot, jtot, ktot) 
    print("at=%f"%at[itot*jtot+itot+itot/2])

    # Time performance 
    start = time.time()
 
    if(mode == 'py'): 
        for i in range(nloop): 
            diff_py(at, a, 0.1, 0.1, 0.1, 0.1, itot, jtot, ktot) 
    elif(mode == 'cy'):
        for i in range(nloop): 
            diff_pyx(at, a, 0.1, 0.1, 0.1, 0.1, itot, jtot, ktot) 
 
    end = time.time()

    print("time/iter = %f s (%i iters)"%((end-start)/float(nloop),nloop))
