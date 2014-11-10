#include <iostream>
#include <iomanip>
#include <cstdlib>

double interp_val(const double m2, const double m1, const double p1, const double p2)
{
  return (-1./16)*m2 + (9./16)*m1 + (9./16)*p1 + (-1./16)*p2;
}

double grad_val(const double m2, const double m1, const double p1, const double p2)
{
  return (1./24.)*m2 + (-27./24.)*m1 + (27./24.)*p1 + (-1./24.)*p2;
}

void grad_interp_val(double * const __restrict__ c, const double * const __restrict__ a,
                     const double * const __restrict__ b, const int n)
{
  for (int i=3; i<n-3; ++i)
    c[i] += grad_val( interp_val( a[i-3], a[i-2], a[i-1], a[i  ] ) * interp_val( b[i-3], b[i-2], b[i-1], b[i  ] ) + interp_val( a[i-3], a[i-2], a[i-1], a[i  ] ),
                      interp_val( a[i-2], a[i-1], a[i  ], a[i+1] ) * interp_val( b[i-2], b[i-1], b[i  ], b[i+1] ) + interp_val( a[i-2], a[i-1], a[i  ], a[i+1] ),
                      interp_val( a[i-1], a[i  ], a[i+1], a[i+2] ) * interp_val( b[i-1], b[i  ], b[i+1], b[i+2] ) + interp_val( a[i-1], a[i  ], a[i+1], a[i+2] ),
                      interp_val( a[i  ], a[i+1], a[i+2], a[i+3] ) * interp_val( b[i  ], b[i+1], b[i+2], b[i+3] ) + interp_val( a[i  ], a[i+1], a[i+2], a[i+3] ));
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

  for (int ii=0; ii<iter; ++ii)
    grad_interp_val(c_data, a_data, b_data, n);

  std::cout << std::setprecision(8) << "c = " << c_data[n-10] << std::endl;

  return 0;
}
