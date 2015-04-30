#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <cstdio>

void init(double* const __restrict__ a, double* const __restrict__ at, const int ncells)
{
    for (int i=0; i<ncells; ++i)
    {
        a[i]  = i*i;
        at[i] = 0.;
    }
}

void diff(double* const __restrict__ at, const double* const __restrict__ a, const double visc, 
          const double dxidxi, const double dyidyi, const double dzidzi, 
          const int itot, const int jtot, const int ktot)
{
    const int ii = 1;
    const int jj = itot;
    const int kk = itot*jtot;

    for (int k=1; k<ktot-1; k++)
        for (int j=1; j<jtot-1; j++)
        #pragma ivdep
            for (int i=1; i<itot-1; i++)
            {
                const int ijk = i + j*jj + k*kk;
                at[ijk] += visc * (
                        + ( (a[ijk+ii] - a[ijk   ]) 
                          - (a[ijk   ] - a[ijk-ii]) ) * dxidxi 
                        + ( (a[ijk+jj] - a[ijk   ]) 
                          - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
                        + ( (a[ijk+kk] - a[ijk   ]) 
                          - (a[ijk   ] - a[ijk-kk]) ) * dzidzi
                        );
            }
}

int main()
{
    const int nloop = 100;
    const int itot = 256;
    const int jtot = 256;
    const int ktot = 256;
    const int ncells = itot*jtot*ktot;

    double *a  = new double[ncells];
    double *at = new double[ncells];
   
    init(a, at, ncells);
   
    for (int i=0; i<nloop; ++i)
        diff(at, a, 0.1, 0.1, 0.1, 0.1, itot, jtot, ktot); 
   
    //printf("at=%f\n",at[itot*jtot+2*itot+2]);
    
    return 0;
}
