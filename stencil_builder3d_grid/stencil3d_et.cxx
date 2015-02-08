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
  const int gc = 4;
  const int iter = 5;

  // Initialize the grid.
  Grid grid(itot, jtot, ktot, gc);

  // Create fields on the grid.
  Field u (grid);
  Field v (grid);
  Field w (grid);
  Field ut(grid);

  // Initialize the fields.
  u.randomize();
  v.randomize();
  w.randomize();

  // Initialize the time step.
  const double dt = 1.e-3;
  const double visc = 1.5;

  // Execute the loop iter times.
  #pragma omp parallel
  {
    for (int n=0; n<iter; ++n)
    {
      // Advection and diffusion operator, split in directions.
      ut += Gx_h( Ix  (u) * Ix  (u) ) + visc * ( Gx_h( Gx  (u) ) );
      ut += Gy  ( Ix_h(v) * Iy_h(u) ) + visc * ( Gy  ( Gy_h(u) ) );
      ut += Gz  ( Ix_h(w) * Iz_h(u) ) + visc * ( Gz  ( Gz_h(u) ) );

      // Time integration.
      u += dt*ut;

      // Tendency reset.
      ut = 0.;
    }
  }

  // Print a value in the middle of the field.
  std::cout << std::setprecision(8) << "u = " << u(itot/2, jtot/2, ktot/2) << std::endl;

  // std::cout << std::endl;
  // std::cout << "Operator type: " << std::endl << getDemangledName(advection) << std::endl;

  return 0;
}
