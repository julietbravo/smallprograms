#include <iostream>

template<typename T> struct IsAllowed { static const bool value = true; };
template<> struct IsAllowed<float> { static const bool value = false; };

template<typename T>
void printValue(T t)
{
  static_assert(IsAllowed<T>::value, "Illegal type");
  std::cout << t << std::endl; 
}

int main()
{
  int i = 10;
  float f = 100.;
  double d = 1000.;

  printValue(i);
  printValue(f);
  printValue(d);

  return 0;
}
