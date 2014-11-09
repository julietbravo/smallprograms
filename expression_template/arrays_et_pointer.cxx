#include <iostream>
#include <iomanip>

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
  void operator=(X<Left, Op, Right> expression)
  {
    for (int i=0; i<n_; ++i)
      data_[i] = expression[i];
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
  double a_data[] = { 0.1, 0.2, 0.3, 0.4 };
  double b_data[] = { 0.2, 0.3, 0.4, 0.5 };
  double c_data[] = { 0.3, 0.4, 0.5, 0.6 };

  double d1_data[4];
  double d2_data[4];
  double d3_data[4];
  double d4_data[4];

  Array A(a_data, 4);
  Array B(b_data, 4);
  Array C(c_data, 4);

  Array D1(d1_data, 4);
  Array D2(d2_data, 4);
  Array D3(d3_data, 4);
  Array D4(d4_data, 4);

  D1 = A + B + C;
  D2 = X<Array, add, Array>(B, C) + A;
  D3 = X< X<Array, add, Array>, add, Array>( X<Array, add, Array>(B, C), A);
  D4.operator=(X< X<Array, add, Array>, add, Array>( X<Array, add, Array>(B, C), A));

  for (int i=0; i<4; ++i)
    std::cout << std::setw(3) << i     << ": "
              << std::setw(3) << D1[i] << ", "
              << std::setw(3) << D2[i] << ", "
              << std::setw(3) << D3[i] << ", "
              << std::setw(3) << D4[i] << std::endl;

  return 0;
}
