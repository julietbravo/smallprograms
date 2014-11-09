#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "StencilBuilder.h"

using namespace StencilBuilder;

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

  Slice a(a_data[0], itot, jtot);
  Slice b(b_data[0], itot, jtot);
  Slice c(c_data[0], itot, jtot);
  for (int ii=0; ii<iter; ++ii)
  {
    c += grad<0>( interp<1>(a, 1     ) * interp<1>(b, itot+6) + interp<1>(c, 1), 1     )
       + grad<0>( interp<1>(a, itot+6) * interp<1>(b, 1     ) + interp<1>(c, 1), itot+6);
  }

  int i  = 512;
  int j  = 512;
  int ij = i+3 + (j+3)*(itot+6);
  std::cout << std::setprecision(8) << "c = " << c_data[ij] << std::endl;

  return 0;
}
