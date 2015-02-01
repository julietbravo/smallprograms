#include <iostream>

template<class Left, class Right>
class Node
{
  public:
    Node(const Left& left, const Right& right) : left_(left), right_(right) {}
    const double operator()() const { return left_() + right_(); }

  private:
    const Left& left_;
    const Right& right_;
};

template<int loc>
class Field
{
  public:
    Field(const double x) : x_(x) {}
    template<class T>
    Field& operator=(const T& a)
    {
      x_ = a();
      return *this;
    }

    const double operator()() const { return x_; }
    double& operator()() { return x_; }

  private:
    double x_;
};

template<class Left, class Right>
inline Node<Left, Right> operator+(const Left& left, const Right& right)
{
  return Node<Left, Right>(left, right);
}

int main()
{
  Field<0> a(0.123);
  Field<1> b(3.343);
  Field<0> c(0.0);

  c = a + b;
  std::cout << c() << std::endl;

  c() = 16;
  std::cout << c() << std::endl;

  return 0;
}
