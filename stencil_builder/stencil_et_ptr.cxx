#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "StencilBuilder_ptr.h"

using namespace StencilBuilder;

double interp_val(const double m2, const double m1, const double p1, const double p2)
{
  return (-1./16)*m2 + (9./16)*m1 + (9./16)*p1 + (-1./16)*p2;
}

double grad_val(const double m2, const double m1, const double p1, const double p2)
{
  return (1./24.)*m2 + (-27./24.)*m1 + (27./24.)*p1 + (-1./24.)*p2;
}

double grad_interp_val(const double * const a)
{
  return grad_val( interp_val( a[-3], a[-2], a[-1], a[0] ) * interp_val( a[-3], a[-2], a[-1], a[0] ), 
                   interp_val( a[-2], a[-1], a[ 0], a[1] ) * interp_val( a[-2], a[-1], a[ 0], a[1] ), 
                   interp_val( a[-1], a[ 0], a[ 1], a[2] ) * interp_val( a[-1], a[ 0], a[ 1], a[2] ), 
                   interp_val( a[ 0], a[ 1], a[ 2], a[3] ) * interp_val( a[ 0], a[ 1], a[ 2], a[3] )); 
}

double grad_interp_val(const double * const a, const double * const b)
{
  return grad_val( interp_val( a[-3], a[-2], a[-1], a[0] ) * interp_val( b[-3], b[-2], b[-1], b[0] ), 
                   interp_val( a[-2], a[-1], a[ 0], a[1] ) * interp_val( b[-2], b[-1], b[ 0], b[1] ), 
                   interp_val( a[-1], a[ 0], a[ 1], a[2] ) * interp_val( b[-1], b[ 0], b[ 1], b[2] ), 
                   interp_val( a[ 0], a[ 1], a[ 2], a[3] ) * interp_val( b[ 0], b[ 1], b[ 2], b[3] )); 
}

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

  double e = 0.;

  for (int ii=0; ii<iter; ++ii)
    for (int i=3; i<n-3; ++i)
    {
      Scalar a(&a_data[i]);
      Scalar b(&b_data[i]);
      Scalar c(&c_data[i]);
      c += grad<0>( interp<1>(a, 1) * interp<1>(b, 1), 1);
    }

  std::cout << std::setprecision(8) << "c = " << c_data[3] << std::endl;

  return 0;
}
