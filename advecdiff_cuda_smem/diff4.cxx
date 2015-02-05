#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>

const double cdg0 = -1460./576.;
const double cdg1 =   783./576.;
const double cdg2 =   -54./576.;
const double cdg3 =     1./576.;
const double cgi = ( 1./24.);

inline double dg4(const double v1, const double v2, const double v3, const double v4, const double v5, const double v6, const double v7) 
{
    return (1./576.)*(v1+v7) + (-54./576.)*(v2+v6) + (783./576.)*(v3+v5) + (-1460./576.)*v4;
}

void diff_cpu(double * const __restrict__ at, const double * const __restrict__ a,
              const double dxidxi, const double dyidyi, const double dzidzi,
              const int istart, const int iend, 
              const int jstart, const int jend, 
              const int kstart, const int kend, 
              const int icells, const int ijcells)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*icells;
    const int jj2 = 2*icells;
    const int jj3 = 3*icells;
    const int kk1 = 1*ijcells;
    const int kk2 = 2*ijcells;
    const int kk3 = 3*ijcells;

    const double visc = 0.1;

    for (int k=kstart; k<kend; ++k)
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                at[ijk] += visc * dg4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3])*dxidxi;
                at[ijk] += visc * dg4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3])*dyidyi;
                at[ijk] += visc * dg4(a[ijk-kk3], a[ijk-kk2], a[ijk-kk1], a[ijk], a[ijk+kk1], a[ijk+kk2], a[ijk+kk3])*dzidzi;
            }
}


int main()
{
    const double dxi = 0.1;
    const double dyi = 0.1;
    const double dzi = 0.1;

    // Test configuration settings.
    const int itot = 256;
    const int jtot = 256;
    const int ktot = 256;
    const int gc   = 3;
    const int iter = 1;

    // Calculate the required variables.
    const int ntot    = (itot+2*gc)*(jtot+2*gc)*(ktot+2*gc);
    const int istart  = gc;
    const int jstart  = gc;
    const int kstart  = gc;
    const int iend    = itot+gc;
    const int jend    = jtot+gc;
    const int kend    = ktot+gc;
    const int icells  = itot+2*gc;
    const int jcells  = jtot+2*gc;
    const int kcells  = ktot+2*gc;
    const int ijcells = (itot+2*gc)*(jtot+2*gc);

    // Allocate the raw arrays.
    double *a  = new double[ntot];
    double *at = new double[ntot];

    // Initialize the raw arrays.
    for (int n=0; n<ntot; ++n)
    {
        a [n] = 0.001 * (std::rand() % 1000) - 0.5;
        at[n] = 0.;
    }

    for(int n=0; n<iter; ++n)
    {
       diff_cpu (at, a, dxi, dyi, dzi, istart, iend, jstart, jend, kstart, kend, icells, ijcells);
    }
    
    int ijk = icells/2 + jcells/2*icells + kcells/2*ijcells;
    printf("%i: value=%e\n",ijk,at[ijk]);

    return 0;
}
