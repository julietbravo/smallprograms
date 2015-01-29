#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "StencilBuilder.h"

using namespace StencilBuilder;

int main()
{
  // Test configuration settings.
  const int itot = 2048;
  const int jtot = 2048;
  const int gc = 4;
  const int niter = 100;

  // Initialize the grid.
  Grid grid(itot, jtot, gc);

  // Create fields on the grid.
  Field u (grid);
  Field v (grid);
  Field ut(grid);

  // Initialize the fields.
  for (int ij=0; ij<grid.ijcells; ++ij)
  {
    u[ij] = 0.001 * (std::rand() % 1000) - 0.5;
    v[ij] = 0.001 * (std::rand() % 1000) - 0.5;

    ut[ij] = 0.;
  }

  // Define the distances in memory for the three dimensions.
  const int ii = 1;
  const int jj = grid.icells;

  // Initialize the time step.
  const double dt = 1.e-3;

  // Execute the loop niter times.
  for (int n=0; n<niter; ++n)
  {
    // Advection operator.
    ut += grad<0>( interp<1>(u, ii) * interp<1>(u, ii), ii )
        + grad<1>( interp<0>(v, ii) * interp<0>(u, jj), jj );

    // Time integration.
    u += dt*ut;

    // Tendency reset.
    ut = 0.;
  }

  // Print a value in the middle of the field.
  std::cout << std::setprecision(8) << "a = " << u(itot/2, jtot/2) << std::endl;

  return 0;
}
