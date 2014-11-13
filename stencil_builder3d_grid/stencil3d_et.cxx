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

  // Initialize the grid.
  Grid grid(itot, jtot, ktot, gc);

  // Map the StencilBuilder Field to the raw array.
  Field a (grid);
  Field b (grid);
  Field c (grid);
  Field at(grid);

  // Initialize the fields.
  for (int n=0; n<grid.ntot; ++n)
  {
    a[n] = 0.001 * (std::rand() % 1000) - 0.5;
    b[n] = 0.001 * (std::rand() % 1000) - 0.5;
    c[n] = 0.001 * (std::rand() % 1000) - 0.5;

    at[n] = 0.;
  }

  // Define the distances in memory for the three dimensions.
  const int ii = 1;
  const int jj = grid.icells;
  const int kk = grid.ijcells;

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

  // Print a value in the middle of the field.
  std::cout << std::setprecision(8) << "a = " << a(itot/2, jtot/2, ktot/2) << std::endl;

  return 0;
}