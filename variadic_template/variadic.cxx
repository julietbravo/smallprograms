#include <iostream>
#include <thread>

void printInt(const int i) { std::cout << i << std::endl; }
void printHello() { std::cout << "Hello!" << std::endl; }
int calcSquare(const int i) { return i*i; }

template<class F, class... Args>
auto doWork(F f, Args... args) -> decltype( f(args...) )
{
  return f(args...);
  // std::thread t(f, args...);
  // t.join();
}

int main()
{
  doWork(printInt, 10);
  doWork(printHello);
  auto c = doWork(calcSquare, 4);

  return 0;
}
