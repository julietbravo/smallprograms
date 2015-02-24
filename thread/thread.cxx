#include <iostream>
#include <thread>
#include <future>

void printInt(const int n)
{ 
  std::cout << n << std::endl;
}

void printHello()
{ 
  std::cout << "Hello!" << std::endl;
}

int calcSquare(const int n)
{ 
  return n*n;
}

template<class F, class... Args>
auto doWork(F& f, Args... args) -> decltype( f(args...) )
{
  std::packaged_task<decltype( f(args...) )(Args...)> task(&f);
  auto fut = task.get_future();
  std::thread t(std::move(task), args...);
  t.join();

  return fut.get();
}

int main()
{
  doWork(printInt, 10);
  doWork(printHello);
  auto c = doWork(calcSquare, 4);

  std::cout << "The value of c = " << c << std::endl;

  return 0;
}
