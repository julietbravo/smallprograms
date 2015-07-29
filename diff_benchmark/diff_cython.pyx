#cython: boundscheck=False
#cython: wraparound=False

import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t

def diff_pyx(np.ndarray[DTYPE_t, ndim=1] at, \
             np.ndarray[DTYPE_t, ndim=1] a, \
             double visc, double dxidxi, double dyidyi, double dzidzi, \
             int itot, int jtot, int ktot):

    cdef int ii = 1
    cdef int jj = itot
    cdef int kk = itot*jtot

    cdef int i,j,k,ijk

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
