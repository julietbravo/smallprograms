#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "StencilBuilder.h"


double interp_val(const double m2, const double m1, const double p1, const double p2)
{
  return (-1./16)*m2 + (9./16)*m1 + (9./16)*p1 + (-1./16)*p2;
}

double grad_val(const double m2, const double m1, const double p1, const double p2)
{
  return (1./24.)*m2 + (-27./24.)*m1 + (27./24.)*p1 + (-1./24.)*p2;
}

void grad_interp_val(double * const __restrict__ c, const double * const __restrict__ a,
                     const double * const __restrict__ b, const int itot, const int jtot)
{
  const int ii1 = 1;
  const int ii2 = 2;
  const int ii3 = 3;
  const int jj1 = 1*(itot+6);
  const int jj2 = 2*(itot+6);
  const int jj3 = 3*(itot+6);

  for (int j=0; j<jtot; ++j)
    for (int i=0; i<itot; ++i)
    {
      int ij = i+3 + (j+3)*(itot+6);
      c[ij] += grad_val( interp_val( a[ij-ii3], a[ij-ii2], a[ij-ii1], a[ij    ] ) * interp_val( b[ij-jj3], b[ij-jj2], b[ij-jj1], b[ij    ] ) + interp_val( c[ij-ii3], c[ij-ii2], c[ij-ii1], c[ij    ] ),
                         interp_val( a[ij-ii2], a[ij-ii1], a[ij    ], a[ij+ii1] ) * interp_val( b[ij-jj2], b[ij-jj1], b[ij    ], b[ij+jj1] ) + interp_val( c[ij-ii2], c[ij-ii1], c[ij    ], c[ij+ii1] ),
                         interp_val( a[ij-ii1], a[ij    ], a[ij+ii1], a[ij+ii2] ) * interp_val( b[ij-jj1], b[ij    ], b[ij+jj1], b[ij+jj2] ) + interp_val( c[ij-ii1], c[ij    ], c[ij+ii1], c[ij+ii2] ),
                         interp_val( a[ij    ], a[ij+ii1], a[ij+ii2], a[ij+ii3] ) * interp_val( b[ij    ], b[ij+jj1], b[ij+jj2], b[ij+jj3] ) + interp_val( c[ij    ], c[ij+ii1], c[ij+ii2], c[ij+ii3] ))

             + grad_val( interp_val( a[ij-jj3], a[ij-jj2], a[ij-jj1], a[ij    ] ) * interp_val( b[ij-ii3], b[ij-ii2], b[ij-ii1], b[ij    ] ) + interp_val( c[ij-jj3], c[ij-jj2], c[ij-jj1], c[ij    ] ),
                         interp_val( a[ij-jj2], a[ij-jj1], a[ij    ], a[ij+jj1] ) * interp_val( b[ij-ii2], b[ij-ii1], b[ij    ], b[ij+ii1] ) + interp_val( c[ij-jj2], c[ij-jj1], c[ij    ], c[ij+jj1] ),
                         interp_val( a[ij-jj1], a[ij    ], a[ij+jj1], a[ij+jj2] ) * interp_val( b[ij-ii1], b[ij    ], b[ij+ii1], b[ij+ii2] ) + interp_val( c[ij-jj1], c[ij    ], c[ij+jj1], c[ij+jj2] ),
                         interp_val( a[ij    ], a[ij+jj1], a[ij+jj2], a[ij+jj3] ) * interp_val( b[ij    ], b[ij+ii1], b[ij+ii2], b[ij+ii3] ) + interp_val( c[ij    ], c[ij+jj1], c[ij+jj2], c[ij+jj3] ));
    }
}

int main()
{
  const int itot = 1024;
  const int jtot = 1024;
  const int iter = 32;

  const int ntot = (itot+6)*(jtot+6);

  double *a_data = new double[ntot];
  double *b_data = new double[ntot];
  double *c_data = new double[ntot];

  for (int i=0; i<ntot; ++i)
  {
    a_data[i] = 0.001 * (std::rand() % 1000) - 0.5;
    b_data[i] = 0.001 * (std::rand() % 1000) - 0.5;
    c_data[i] = 0.;
  }

  for (int ii=0; ii<iter; ++ii)
    grad_interp_val(c_data, a_data, b_data, itot, jtot);

  int i  = 512;
  int j  = 512;
  int ij = i+3 + (j+3)*(itot+6);
  std::cout << std::setprecision(8) << "c = " << c_data[ij] << std::endl;

  return 0;
}
