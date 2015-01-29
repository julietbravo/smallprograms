#include <iostream>
#include <iomanip>
#include <cstdlib>

#define restrict __restrict__

struct Grid
{
  Grid(const int itot, const int jtot, const int gc) :
    itot(itot), jtot(jtot), gc(gc),
    istart(gc),
    jstart(gc),
    iend(itot+gc),
    jend(jtot+gc),
    icells(itot+2*gc),
    ijcells((itot+2*gc)*(jtot+2*gc)) {}

  const int itot;
  const int jtot;
  const int gc;
  const int istart;
  const int jstart;
  const int iend;
  const int jend;
  const int icells;
  const int ijcells;
};

// Field class representing the field, whose operations expand compile time.
class Field
{
  public:
    Field(const Grid& grid) :
      grid_(grid),
      data_(new double[grid_.ijcells]) {}

    Field(const Field& a) :
      grid_(a.grid_),
      data_(new double[grid_.ijcells])
    {
      const int jj = grid_.icells;

      for (int j=grid_.jstart; j<grid_.jend; ++j)
        #pragma ivdep
        for (int i=grid_.istart; i<grid_.iend; ++i)
        {
          const int ij = i + j*jj;
          data_[ij] = a[ij];
        }
    }

    ~Field() { delete[] data_; }

    inline double* get_data() const { return data_; }
    inline const Grid& get_grid() const { return grid_; }

    inline double& operator[](const int i) const { return data_[i]; }

    inline double& operator()(const int i, const int j) const
    { return data_[i + j*grid_.icells]; }

    inline Field& operator= (const Field& a)
    {
      if (this != &a)
      {
        const int jj = grid_.icells;

        for (int j=grid_.jstart; j<grid_.jend; ++j)
          #pragma ivdep
          for (int i=grid_.istart; i<grid_.iend; ++i)
          {
            const int ij = i + j*jj;
            data_[ij] = a[ij];
          }
      }

      return *this;
    }

    inline Field& operator= (const double a)
    {
      const int jj = grid_.icells;

      for (int j=grid_.jstart; j<grid_.jend; ++j)
        #pragma ivdep
        for (int i=grid_.istart; i<grid_.iend; ++i)
        {
          const int ij = i + j*jj;
          data_[ij] = a;
        }

      return *this;
    }

    inline Field& operator+=(const Field& a)
    {
      const int jj = grid_.icells;

      for (int j=grid_.jstart; j<grid_.jend; ++j)
        #pragma ivdep
        for (int i=grid_.istart; i<grid_.iend; ++i)
        {
          const int ij = i + j*jj;
          data_[ij] += a[ij];
        }

      return *this;
    }

    inline Field& operator*=(const Field& a)
    {
      const int jj = grid_.icells;

      for (int j=grid_.jstart; j<grid_.jend; ++j)
        #pragma ivdep
        for (int i=grid_.istart; i<grid_.iend; ++i)
        {
          const int ij = i + j*jj;
          data_[ij] *= a[ij];
        }

      return *this;
    }

    inline Field& operator*=(const double a)
    {
      const int jj = grid_.icells;

      for (int j=grid_.jstart; j<grid_.jend; ++j)
        #pragma ivdep
        for (int i=grid_.istart; i<grid_.iend; ++i)
        {
          const int ij = i + j*jj;
          data_[ij] *= a;
        }

      return *this;
    }

  private:
    // Reference to the grid on which the field is created
    const Grid& grid_;
    // Pointer to the raw data.
    double* restrict data_;
};

Field operator+(Field a, const Field& b)
{
  a += b;
  return a;
}

Field operator*(Field a, const Field& b)
{
  a *= b;
  return a;
}

Field operator*(const double a, Field b)
{
  b *= a;
  return b;
}

struct Interp
{
  static inline double apply(const double a, const double b, const double c, const double d)
  { return (9./16.)*(b+c) - (1./16.)*(a+d); }
};

struct Grad
{
  static inline double apply(const double a, const double b, const double c, const double d)
  { return (27./24.)*(c-b) - 1./24*(d-a); }
};

template<int toCenter>
Field interp(const Field& a, const int nn)
{
  const Grid& grid = a.get_grid();
  Field ai(a.get_grid());

  const int jj = grid.icells;

  for (int j=grid.jstart; j<grid.jend; ++j)
    #pragma ivdep
    for (int i=grid.istart; i<grid.iend; ++i)
    {
      const int ij = i + j*jj;
      ai.get_data()[ij] = Interp::apply(a[i + (-2+toCenter)*nn], a[i + (-1+toCenter)*nn],
                                        a[i + (   toCenter)*nn], a[i + ( 1+toCenter)*nn]);
    }

  return ai;
}

template<int toCenter>
Field grad(const Field& a, const int nn)
{
  const Grid& grid = a.get_grid();
  Field ag(grid);

  const int jj = grid.icells;

  for (int j=grid.jstart; j<grid.jend; ++j)
    #pragma ivdep
    for (int i=grid.istart; i<grid.iend; ++i)
    {
      const int ij = i + j*jj;
      ag.get_data()[ij] = Grad::apply(a[i + (-2+toCenter)*nn], a[i + (-1+toCenter)*nn],
                                     a[i + (   toCenter)*nn], a[i + ( 1+toCenter)*nn]);
    }

  return ag;
}

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
