#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "StencilBuilder_ref.h"

using namespace StencilBuilder;

int main()
{
  const int n    = 1024*1024;
  const int iter = 256;

  double *a_data = new double[n];
  double *b_data = new double[n];
  double *c_data = new double[n];

  for (int i=0; i<n; ++i)
  {
    a_data[i] = 0.001 * (std::rand() % 1000) - 0.5;
    b_data[i] = 0.001 * (std::rand() % 1000) - 0.5;
    c_data[i] = 0.;
  }

  Vector a(a_data[0], n);
  Vector b(b_data[0], n);
  Vector c(c_data[0], n);
  for (int ii=0; ii<iter; ++ii)
    c += grad<0>( interp<1>(a, 1) * interp<1>(b, 1), 1);

  std::cout << std::setprecision(8) << "c = " << c_data[3] << std::endl;

  return 0;
}
