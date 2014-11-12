namespace StencilBuilder
{
  // Stencil operators
  struct Interp
  {
    static double apply_m2(const double a) { return (-1./16.)*a; }
    static double apply_m1(const double a) { return ( 9./16.)*a; }
    static double apply_p1(const double a) { return ( 9./16.)*a; }
    static double apply_p2(const double a) { return (-1./16.)*a; }
  };

  struct Grad
  {
    static double apply_m2(const double a) { return (  1./24.)*a; }
    static double apply_m1(const double a) { return (-27./24.)*a; }
    static double apply_p1(const double a) { return ( 27./24.)*a; }
    static double apply_p2(const double a) { return ( -1./24.)*a; }
  };

  // Stencil aggregation class
  template<int toCenter, class Inner, class Op>
  struct Stencil
  {
    Stencil(const Inner &inner, const int nn) : inner_(inner), nn_(nn) {}

    const Inner &inner_;
    const int nn_;

    double operator[](const int i) const
    {
      const double sm2 = Op::apply_m2(inner_[i + (-2+toCenter)*nn_]);
      const double sm1 = Op::apply_m1(inner_[i + (-1+toCenter)*nn_]);
      const double sp1 = Op::apply_p1(inner_[i + (   toCenter)*nn_]);
      const double sp2 = Op::apply_p2(inner_[i + (+1+toCenter)*nn_]);
      return sm2 + sm1 + sp1 + sp2;
    }
  };

  // Template classes for the stencil operators
  template<int toCenter, class Inner>
  Stencil<toCenter, Inner, Interp> interp(const Inner &inner, const int nn)
  {
    return Stencil<toCenter, Inner, Interp>(inner, nn);
  }

  template<int toCenter, class Inner>
  Stencil<toCenter, Inner, Grad> grad(const Inner &inner, const int nn)
  {
    return Stencil<toCenter, Inner, Grad>(inner, nn);
  }

  // Scalar operators
  struct Multiply
  {
    static double apply(const double left, const double right) { return left*right; }
  };

  struct Add
  {
    static double apply(const double left, const double right) { return left+right; }
  };

  // Operator aggregation class
  template<class Left, class Op, class Right>
  struct Operator
  {
    Operator(const Left &left, const Right &right) : left_(left), right_(right) {}

    const Left &left_;
    const Right &right_;

    double operator[](const int i) const { return Op::apply(left_[i], right_[i]); }
  };

  // Template classes for the operators
  template<class Left, class Right>
  Operator<Left, Multiply, Right> operator*(const Left &left, const Right &right)
  {
    return Operator<Left, Multiply, Right>(left, right);
  }

  template<class Left, class Right>
  Operator<Left, Add, Right> operator+(const Left &left, const Right &right)
  {
    return Operator<Left, Add, Right>(left, right);
  }

  /*
  // Scalar class representing the scalar, whose operations expand compile time
  struct Scalar
  {
    Scalar(double *data) : data_(data) {}

    double operator[](const int i) const { return data_[i]; }

    template<class T> void operator= (const T &expression) { data_[0] =  expression[0]; }
    template<class T> void operator+=(const T &expression) { data_[0] += expression[0]; }

    double *data_;
  };
  */

  // Field class representing the field, whose operations expand compile time
  struct Field 
  {
    Field(double *data, 
          const int istart, const int iend,
          const int jstart, const int jend,
          const int kstart, const int kend,
          const int icells, const int ijcells)
      : data_(data),
        istart_(istart), iend_(iend),
        jstart_(jstart), jend_(jend),
        kstart_(kstart), kend_(kend),
        icells_(icells), ijcells_(ijcells) {}

    double operator[](const int i) const { return data_[i]; }

    // Assignment
    template<class T>
    void operator= (const T &expression)
    {
      const int ii = 1;
      const int jj = icells_;
      const int kk = ijcells_;

      for (int k=kstart_; k<kend_; ++k)
        for (int j=jstart_; j<jend_; ++j)
          for (int i=istart_; i<iend_; ++i)
          {
            const int ijk = i + j*jj + k*kk;
            data_[ijk] = expression[ijk];
          }
    }

    template<class T> void operator+=(const T &expression)
    {
      const int ii = 1;
      const int jj = icells_;
      const int kk = ijcells_;

      for (int k=kstart_; k<kend_; ++k)
        for (int j=jstart_; j<jend_; ++j)
          for (int i=istart_; i<iend_; ++i)
          {
            const int ijk = i + j*jj + k*kk;
            data_[ijk] += expression[ijk];
          }
    }

    double *data_;

    const int istart_;
    const int jstart_;
    const int kstart_;
    const int iend_;
    const int jend_;
    const int kend_;
    const int icells_;
    const int ijcells_;
  };

  // Specialization for assignment with a constant
  template<>
  void Field::operator= (const double &expression)
  {
    const int ii = 1;
    const int jj = icells_;
    const int kk = ijcells_;
  
    for (int k=kstart_; k<kend_; ++k)
      for (int j=jstart_; j<jend_; ++j)
        for (int i=istart_; i<iend_; ++i)
        {
          const int ijk = i + j*jj + k*kk;
          data_[ijk] = expression;
        }
  }
}

