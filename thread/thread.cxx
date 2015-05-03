#include <iostream>
#include <thread>
#include <future>

void printInt(const int n)
{ 
  std::cout << n << std::endl;
  //std::this_thread::sleep_for(std::chrono::seconds(2));
}

void printHello()
{ 
  std::cout << "Hello!" << std::endl;
  //std::this_thread::sleep_for(std::chrono::seconds(2));
}

int calcSquare(const int n)
{ 
  //std::this_thread::sleep_for(std::chrono::seconds(2));
  return n*n;
}

template<class F, class... Args>
auto doWork(F& f, Args... args) -> decltype( f(args...) )
{
  std::packaged_task<decltype( f(args...) )(Args...)> task(&f);
  auto fut = task.get_future();
  std::thread tloc(std::move(task), args...);
  tloc.join();
  return fut.get();
}


class ThreadPool
{
  public:
    ThreadPool() : active_(true)
    {
      t1_ = std::thread(&ThreadPool::listen, this);
      t2_ = std::thread(&ThreadPool::listen, this);
    }

    ~ThreadPool()
    {
      std::cout << "Terminating thread pool!" << std::endl;
      active_ = false;
      t1_.join();
      t2_.join();
    }

  private:
    void listen()
    {
      while (active_)
      {
        std::this_thread::sleep_for(std::chrono::seconds(3));
      }
    }

    std::atomic<bool> active_;
    std::thread t1_;
    std::thread t2_;
};

int main()
{
  ThreadPool threadPool;

  doWork(printInt, 10);
  doWork(printHello);
  auto c = doWork(calcSquare, 4);

  std::cout << "The value of c = " << c << std::endl;

  return 0;
}
