#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "StencilBuilder.h"

using namespace StencilBuilder;

int main()
{
  // Test configuration settings.
  const int itot = 256;
  const int jtot = 256;
  const int ktot = 256;
  const int gc   = 4;
  const int iter = 5;

  // Calculate the required variables.
  const int ntot = (itot+2*gc)*(jtot+2*gc)*(ktot+2*gc);
  const int istart = gc;
  const int jstart = gc;
  const int kstart = gc;
  const int iend = itot+gc;
  const int jend = jtot+gc;
  const int kend = ktot+gc;
  const int icells = itot+2*gc;
  const int ijcells = (itot+2*gc)*(jtot+2*gc);

  // Allocate the raw arrays.
  double *a_data  = new double[ntot];
  double *b_data  = new double[ntot];
  double *c_data  = new double[ntot];
  double *at_data = new double[ntot];

  // Initialize the raw arrays.
  for (int n=0; n<ntot; ++n)
  {
    a_data[n] = 0.001 * (std::rand() % 1000) - 0.5;
    b_data[n] = 0.001 * (std::rand() % 1000) - 0.5;
    c_data[n] = 0.001 * (std::rand() % 1000) - 0.5;
    at_data[n] = 0.;
  }

  // Map the StencilBuilder Field to the raw array.
  Field a ( a_data, istart, iend, jstart, jend, kstart, kend, icells, ijcells);
  Field b ( b_data, istart, iend, jstart, jend, kstart, kend, icells, ijcells);
  Field c ( c_data, istart, iend, jstart, jend, kstart, kend, icells, ijcells);
  Field at(at_data, istart, iend, jstart, jend, kstart, kend, icells, ijcells);

  // Define the distances in memory for the three dimensions.
  const int ii = 1;
  const int jj = icells;
  const int kk = ijcells;

  // Initialize the time step.
  const double dt = 1.e-3;

  // Execute the loop iter times.
  for (int n=0; n<iter; ++n)
  {
    // Advection operator.
    at += grad<0>( interp<1>(a, ii) * interp<1>(a, ii), ii )
        + grad<1>( interp<0>(b, ii) * interp<0>(a, jj), jj )
        + grad<1>( interp<0>(c, ii) * interp<0>(a, kk), kk );

    // Time integration.
    a += dt*at;

    // Tendency reset.
    at = 0.;
  }

  // Print a value in the middle of the field to check whether
  // both versions give the same result.
  const int ijk = itot/2 + (jtot/2)*icells + (ktot/2)*ijcells;
  std::cout << std::setprecision(8) << "a = " << a(itot/2, jtot/2, ktot/2) << std::endl;

  return 0;
}
