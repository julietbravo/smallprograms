#include <iostream>
#include <random>

inline double sign(const double n) { return ( (n > 0) ? 1 : ( (n < 0) ? -1 : 0) ); }

int main()
{
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(-1., 1.);

  for (int i=0; i<10; ++i)
  {
    const double d_random = dist(mt);
    std::cout << d_random 
              << ", " << sign(d_random) * d_random
              << std::endl;
  }

  return 0;
}
