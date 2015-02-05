#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>

const double cgi = ( 1./24.);

inline double i4(const double m2, const double m1, const double p1, const double p2) {return -1./16.*m2 + 9./16.*m1 + 9./16.*p1 - 1./16.*p2;}
inline double g4(const double m2, const double m1, const double p1, const double p2) {return m2 - 27.*m1 + 27.*p1 - p2;}

void advecu_cpu(double * const __restrict__ ut, 
                const double * const __restrict__ u, const double * const __restrict__ v, const double * const __restrict__ w, 
                const double dxi, const double dyi, const double dzi4,
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

    for (int k=kstart; k<kend; ++k)
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] -= g4( i4( u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]             ) * i4( u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ] ),
                               i4( u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]             ) * i4( u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1] ),
                               i4( u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]             ) * i4( u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2] ),
                               i4( u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]             ) * i4( u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3] )) * cgi * dxi;

                ut[ijk] -= g4( i4( v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1] ) * i4( u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ] ),
                               i4( v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ] ) * i4( u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1] ),
                               i4( v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1] ) * i4( u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2] ),
                               i4( v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2] ) * i4( u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3] )) * cgi * dyi;

                ut[ijk] -= g4( i4( w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1] ) * i4( u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ] ),
                               i4( w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ] ) * i4( u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1] ),
                               i4( w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1] ) * i4( u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2] ),
                               i4( w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2] ) * i4( u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3] )) * dzi4;
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
    double *u  = new double[ntot];
    double *v  = new double[ntot];
    double *w  = new double[ntot];
    double *ut1 = new double[ntot];
    double *ut2 = new double[ntot];

    // Initialize the raw arrays.
    for (int n=0; n<ntot; ++n)
    {
        u [n] = 0.001 * (std::rand() % 1000) - 0.5;
        v [n] = 0.001 * (std::rand() % 1000) - 0.5;
        w [n] = 0.001 * (std::rand() % 1000) - 0.5;
        ut1[n] = 0.;
        ut2[n] = 0.;
    }

    for(int n=0; n<iter; ++n)
       advecu_cpu (ut1, u, v, w, dxi, dyi, dzi, istart, iend, jstart, jend, kstart, kend, icells, ijcells);
       advecu_cpu2(ut2, u, v, w, dxi, dyi, dzi, istart, iend, jstart, jend, kstart, kend, icells, ijcells);

    
    int ijk = icells/2 + jcells/2*icells + kcells/2*ijcells;
    printf("%i: diff=%e\n",ijk,ut1[ijk]-ut2[ijk]);


    return 0;
}
