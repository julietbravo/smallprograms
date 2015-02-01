#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "StencilBuilder.h"

using namespace StencilBuilder;

int main()
{
  // Test configuration settings.
  const int itot = 1024;
  const int jtot = 1024;
  const int gc = 4;
  const int niter = 100;

  // Initialize the grid.
  Grid grid(itot, jtot, gc);

  // Create fields on the grid, default value is 0.
  Field u (grid);
  Field v (grid);
  Field ut(grid);

  // Randomize the u and v field.
  u.randomize();
  v.randomize();

  // Initialize the time step.
  const double dt = 1.e-3;

  // Execute the loop niter times.
  for (int n=0; n<niter; ++n)
  {
    // Advection operator.
    ut += Gx_h( Ix  (u) * Ix  (u) )
        + Gy  ( Ix_h(v) * Iy_h(u) );

    // Time integration.
    u += dt*ut;

    // Tendency reset.
    ut = 0.;
  }

  // Print a value in the middle of the field.
  std::cout << std::setprecision(8) << "a = " << u(itot/2, jtot/2) << std::endl;

  return 0;
}
