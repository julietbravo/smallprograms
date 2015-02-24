#include <iostream>

void printInt(const int i) { std::cout << i << std::endl; }
void printHello() { std::cout << "Hello!" << std::endl; }

template<class F, class... Args>
void doWork(F f, Args... args)
{
  f(args...);
}

int main()
{
  doWork(printInt, 10);
  doWork(printHello);

  return 0;
}
