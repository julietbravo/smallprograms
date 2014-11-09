#include <iostream>
#include <iomanip>
#include <cstdlib>

// This template describes the array that contains the data and the operators
struct Array
{
  Array(double *data, const int n) : data_(data), n_(n) {}

  void operator=(const Array &expression)
  {
    for (int i=0; i<n_; ++i)
      data_[i] = expression[i];
  }

  void operator+=(const Array &expression)
  {
    for (int i=0; i<n_; ++i)
      data_[i] += expression[i];
  }

  double operator[](const int i) const
  {
    return data_[i];
  }

  friend Array operator+(const Array &a, const Array &b);

  double *data_;
  const int n_;
};

Array operator+(const Array &a, const Array &b)
{
  double *data = new double[a.n_];
  Array c(data, a.n_);

  for (int i=0; i<c.n_; ++i)
    c.data_[i] = a.data_[i] + b.data_[i];

  return c;
}

void tripleAddition(double * const  __restrict__ d, const double * const __restrict__ a,
                    const double * const  __restrict__ b, const double * const __restrict__ c,
                    const int n)
{
  for (int i=0; i<n; ++i)
    d[i] += a[i] + b[i] + c[i];
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

  // for (int i=0; i<iter; ++i)
  //   D += A + B + C;

  for (int i=0; i<iter; ++i)
    tripleAddition(d_data, a_data, b_data, c_data, n);

  for (int i=0; i<4; ++i)
    std::cout << std::setw(3) << i    << ": "
              << std::setw(3) << D[i] << std::endl;

  return 0;
}
