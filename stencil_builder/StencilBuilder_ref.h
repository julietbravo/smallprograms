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

    double eval() const { return operator[](0); }
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

  // Scalar class representing the scalar, whose operations expand compile time
  struct Scalar
  {
    Scalar(double &data) : data_(data) {}

    double operator[](const int i) const { return (&data_)[i]; }

    template<class T> void operator= (T expression) { (&data_)[0] =  expression[0]; }
    template<class T> void operator+=(T expression) { (&data_)[0] += expression[0]; }

    double &data_;
  };
}
