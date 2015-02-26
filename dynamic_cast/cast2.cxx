#include <iostream>

struct Grid
{
  virtual void print() const = 0;
};

struct GridSecondOrder : public Grid
{
  void print() const { std::cout << "GridSecondOrder" << std::endl; }
};

struct GridFourthOrder : public Grid
{
  void print() const { std::cout << "GridFourthOrder" << std::endl; }
};

struct Advec
{
  virtual void print(const Grid& grid) const = 0;
};

struct AdvecSecondOrder : public Advec
{
  void print(const Grid& grid) const
  {
    dynamic_cast<const GridSecondOrder&>(grid).print();
    std::cout << "AdvecSecondOrder" << std::endl;
  }
};

struct AdvecFourthOrder : public Advec
{
  void print(const Grid& grid) const
  {
    dynamic_cast<const GridFourthOrder&>(grid).print();
    std::cout << "AdvecFourthOrder" << std::endl;
  }
};

int main()
{
  try
  {
    GridFourthOrder grid;
    AdvecSecondOrder advec;
    advec.print(grid);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
