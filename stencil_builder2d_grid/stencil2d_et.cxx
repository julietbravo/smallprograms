// StencilBuilder example
// Copyright (c) Chiel van Heerwaarden, 2015
// chielvanheerwaarden@gmail.com

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
  const int iter = 50;

  // Initialize the grid.
  Grid grid(itot, jtot, gc);

  // Create fields on the grid.
  Field u (grid);
  Field v (grid);
  Field ut(grid);

  // Initialize the fields.
  u.randomize();
  v.randomize();

  // Initialize the time step.
  const double dt = 1.e-3;

  // Execute the loop iter times.
  for (int n=0; n<iter; ++n)
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
  std::cout << std::setprecision(8) << "u = " << u(itot/2, jtot/2) << std::endl;

  // Print the nested type.
  auto advection = Gx_h( Ix  (u) * Ix  (u) )
                 + Gy  ( Ix_h(v) * Iy_h(u) );

  std::cout << std::endl;
  std::cout << "Operator type: " << std::endl << getDemangledName(advection) << std::endl;

  return 0;
}
