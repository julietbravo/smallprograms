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

  // Print the nested type.
  auto step1 = Ix  (u);
  std::cout << std::endl;
  std::cout << getDemangledName(step1) << std::endl;

  // Print the nested type.
  auto step2 = Ix  (u) * Ix  (u);
  std::cout << std::endl;
  std::cout << getDemangledName(step2) << std::endl;
             
  // Print the nested type.
  auto advection = Gx_h( Ix  (u) * Ix  (u) )
                 + Gy  ( Ix_h(v) * Iy_h(u) );
  std::cout << std::endl;
  std::cout << getDemangledName(advection) << std::endl;

  return 0;
}
