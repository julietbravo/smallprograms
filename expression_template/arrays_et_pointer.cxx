#include <iostream>
#include <iomanip>
#include <cstdlib>

// This is a class for the add operator
struct add
{
  static double apply(const double a, const double b) { return a+b; }
};

// This is the template that parses the expression tree
template<class Left, class Op, class Right>
struct X
{
  X(Left a, Right b) : a_(a), b_(b) {}

  Left a_;
  Right b_;

  double operator[](const int i) { return Op::apply(a_[i], b_[i]); }
};

// This template describes the array that contains the data and the operators
struct Array
{
  Array(double *data, const int n) : data_(data), n_(n) {}

  template<class Left, class Op, class Right>
  void operator+=(X<Left, Op, Right> expression)
  {
    for (int i=0; i<n_; ++i)
      data_[i] += expression[i];
  }

  double operator[](const int i) { return data_[i]; }

  double *data_;
  const int n_;
};

template<class Left, class Right>
X<Left, add, Right> operator+(const Left a, const Right b)
{
  return X<Left, add, Right>(a, b);
}

int main()
{
  const int n    = 1024*1024;
  const int iter = 1024;

  double *a_data = new double[n];
  double *b_data = new double[n];
  double *c_data = new double[n];
  double *d_data = new double[n];

  for (int i=0; i<n; ++i)
  {
    a_data[i] = 0.001 * (std::rand() % 1000) - 0.5;
    b_data[i] = 0.001 * (std::rand() % 1000) - 0.5;
    c_data[i] = 0.001 * (std::rand() % 1000) - 0.5;
    d_data[i] = 0.;
  }

  Array A(a_data, n);
  Array B(b_data, n);
  Array C(c_data, n);
  Array D(d_data, n);

  for (int i=0; i<iter; ++i)
    D += A + B + C;

  // D += X<Array, add, Array>(B, C) + A;
  // D += X< X<Array, add, Array>, add, Array>( X<Array, add, Array>(B, C), A);
  // D.operator+=(X< X<Array, add, Array>, add, Array>( X<Array, add, Array>(B, C), A));

  for (int i=0; i<4; ++i)
    std::cout << std::setw(3) << i    << ": "
              << std::setw(3) << D[i] << std::endl;

  return 0;
}
