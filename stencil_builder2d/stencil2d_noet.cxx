#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "StencilBuilder.h"

using namespace StencilBuilder;

#define restrict __restrict__

// Fourth order interpolation function.
inline double interp(const double m2, const double m1, const double p1, const double p2)
{
  return (-1./16)*(m2+p2) + (9./16)*(m1+p1);
}

// Fourth order gradient function.
inline double grad(const double m2, const double m1, const double p1, const double p2)
{
  return (1./24.)*(m2-p2) + (27./24.)*(p1-m1);
}

// Test function with a similar structure as the advection operator.
void advection(double * const restrict ut, const double * const restrict u,
               const double * const restrict v,
               const int istart, const int iend,
               const int jstart, const int jend,
               const int icells)
{
  const int ii1 = 1;
  const int ii2 = 2;
  const int ii3 = 3;
  const int jj1 = 1*icells;
  const int jj2 = 2*icells;
  const int jj3 = 3*icells;

  for (int j=jstart; j<jend; ++j)
    #pragma ivdep
    for (int i=istart; i<iend; ++i)
    {
      const int ij = i + j*jj1;
      ut[ij] += grad( interp( u[ij-ii3], u[ij-ii2], u[ij-ii1], u[ij    ] ) * interp( u[ij-ii3], u[ij-ii2], u[ij-ii1], u[ij    ] ),
                      interp( u[ij-ii2], u[ij-ii1], u[ij    ], u[ij+ii1] ) * interp( u[ij-ii2], u[ij-ii1], u[ij    ], u[ij+ii1] ),
                      interp( u[ij-ii1], u[ij    ], u[ij+ii1], u[ij+ii2] ) * interp( u[ij-ii1], u[ij    ], u[ij+ii1], u[ij+ii2] ),
                      interp( u[ij    ], u[ij+ii1], u[ij+ii2], u[ij+ii3] ) * interp( u[ij    ], u[ij+ii1], u[ij+ii2], u[ij+ii3] ))

              + grad( interp( v[ij-ii2-jj1], v[ij-ii1-jj1], v[ij-jj1], v[ij+ii1-jj1] ) * interp( u[ij-jj3], u[ij-jj2], u[ij-jj1], u[ij    ] ),
                      interp( v[ij-ii2    ], v[ij-ii1    ], v[ij    ], v[ij+ii1    ] ) * interp( u[ij-jj2], u[ij-jj1], u[ij    ], u[ij+jj1] ),
                      interp( v[ij-ii2+jj1], v[ij-ii1+jj1], v[ij+jj1], v[ij+ii1+jj1] ) * interp( u[ij-jj1], u[ij    ], u[ij+jj1], u[ij+jj2] ),
                      interp( v[ij-ii2+jj2], v[ij-ii1+jj2], v[ij+jj2], v[ij+ii1+jj2] ) * interp( u[ij    ], u[ij+jj1], u[ij+jj2], u[ij+jj3] ));
    }
}

  // Test function for time integration.
void tendency(double * const restrict ut, double * const restrict u,
              const double dt,
              const int istart, const int iend,
              const int jstart, const int jend,
              const int icells)
{
  const int jj = icells;

  for (int j=jstart; j<jend; ++j)
    #pragma ivdep
    for (int i=istart; i<iend; ++i)
    {
      const int ij = i + j*jj;
      u[ij] += dt*ut[ij];
    }

  for (int j=jstart; j<jend; ++j)
    #pragma ivdep
    for (int i=istart; i<iend; ++i)
    {
      const int ij = i + j*jj;
      ut[ij] = 0.;
    }
}

int main()
{
  // Test configuration settings.
  const int itot = 256;
  const int jtot = 256;
  const int gc = 4;
  const int niter = 5;

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
    advection(ut.get_data(), u.get_data(), v.get_data(),
              grid.istart, grid.iend,
              grid.jstart, grid.jend,
              grid.icells);

    tendency(ut.get_data(), u.get_data(),
             dt,
             grid.istart, grid.iend,
             grid.jstart, grid.jend,
             grid.icells);
  }

  // Print a value in the middle of the field.
  std::cout << std::setprecision(8) << "a = " << u(itot/2, jtot/2) << std::endl;

  return 0;
}
