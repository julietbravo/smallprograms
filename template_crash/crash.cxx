struct Grid
{
  Grid(const int itot, const int gc) :
    itot(itot), gc(gc), istart(gc), iend(itot+gc), icells(itot+2*gc) {}

  const int itot;
  const int gc;
  const int istart;
  const int iend;
  const int icells;
};

template<int loc, class Inner>
struct Interp
{
  Interp(const Inner& inner) : inner_(inner) {}

  const Inner& inner_;

  inline double operator()(const int i) const
  {
    return   (-1./16)*(inner_(i + (-2+loc)) + inner_(i + ( 1+loc)))
           + ( 9./16)*(inner_(i + (-1+loc)) + inner_(i + (   loc)));
  }
};

template<class Inner>
inline Interp<1, Inner> Ix(const Inner& inner)
{ return Interp<1, Inner>(inner); }

template<class Inner>
inline Interp<0, Inner> Ix_h(const Inner& inner)
{ return Interp<0, Inner>(inner); }

class Field
{
  public:
    Field(const Grid& grid) :
      grid_(grid),
      data_(new double[grid_.icells]) {}

    inline double operator()(const int i) const
    { return data_[i]; }

    inline double& operator()(const int i)
    { return data_[i]; }

    template<class T>
    inline Field& operator+=(const T& expression)
    {
      for (int i=grid_.istart; i<grid_.iend; ++i)
        (*this)(i) += expression(i);

      return *this;
    }

  private:
    const Grid& grid_;
    double* data_;
};

int main()
{
  Grid grid(256, 4);

  Field u (grid);
  Field ut(grid);

  // Like this it crashes
  auto expression = Ix_h( Ix(u) );
  ut += expression;

  // But like this it does not
  ut += Ix_h( Ix(u) );

  return 0;
}
