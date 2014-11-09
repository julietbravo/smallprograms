#include <iostream>
#include <iomanip>

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

Array operator*(const Array &a, const Array &b)
{
  double *data = new double[a.n_];
  Array c(data, a.n_);

  for (int i=0; i<c.n_; ++i)
    c.data_[i] = a.data_[i] * b.data_[i];

  return c;
}


int main()
{
  double a_data[] = { 0.1, 0.2, 0.3, 0.4 };
  double b_data[] = { 0.2, 0.3, 0.4, 0.5 };
  double c_data[] = { 0.3, 0.4, 0.5, 0.6 };

  double d_data[4];

  Array A(a_data, 4);
  Array B(b_data, 4);
  Array C(c_data, 4);

  Array D(d_data, 4);

  D = A + B * C;

  for (int i=0; i<4; ++i)
    std::cout << std::setw(3) << i    << ": "
              << std::setw(3) << D[i] << std::endl;

  return 0;
}
