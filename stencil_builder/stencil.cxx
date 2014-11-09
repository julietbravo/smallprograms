#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "StencilBuilder.h"

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

  double a_data[7], b_data[7], c_data[7], d_data[7];

  for(int i=0; i<7; ++i)
  {
    a_data[i] = 0.01 * static_cast<double>(std::rand() % 100);
    b_data[i] = 0.01 * static_cast<double>(std::rand() % 100);
    c_data[i] = 0.;
    d_data[i] = 0.;
  }

  Scalar a(a_data[3]);
  Scalar b(b_data[3]);
  Scalar c(c_data[3]);
  Scalar d(d_data[3]);

  double e = 0.;

  c += grad<0>( interp<1>(a, 1) * interp<1>(a, 1), 1);
  d += grad<0>( interp<1>(a, 1) * interp<1>(b, 1), 1);
  e += grad<0>( interp<1>(a, 1) * interp<1>(b, 1), 1).eval() * 2.;

  double c_ref = grad_interp_val(&a_data[3]);
  double d_ref = grad_interp_val(&a_data[3], &b_data[3]);
  double e_ref = 2.*d_ref;

  std::cout << std::setprecision(8) << "c = " << c_data[3] << ", c_ref = " << c_ref << std::endl;
  std::cout << std::setprecision(8) << "d = " << d_data[3] << ", d_ref = " << d_ref << std::endl;
  std::cout << std::setprecision(8) << "e = " << e         << ", e_ref = " << e_ref << std::endl;

  // std::cout << "Validation:" << std::endl;
  // for(int i=0; i<7; ++i)
  //   std::cout << i << ": "
  //             << std::setw(5) << a_data[i] << ", "
  //             << std::setw(5) << b_data[i] << ", "
  //             << std::setw(5) << c_data[i] << ", "
  //             << std::setw(5) << d_data[i] << std::endl;

  return 0;
}
