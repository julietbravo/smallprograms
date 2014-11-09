#include <iostream>

template<class LanguagePolicy>
class Printer : public LanguagePolicy
{
  public:
    Printer(const std::string name) : name_(name) {}
    void print() { LanguagePolicy::printImplementation(name_); }

  private:
    const std::string name_;
};

struct German
{
  void printImplementation(const std::string &name)
  { std::cout << "Gutentag, " << name << "!" << std::endl; }
};

struct English
{
  void printImplementation(const std::string &name)
  { std::cout << "Good afternoon, " << name << "!" << std::endl; }
};

int main()
{
  Printer<German>  printer1(std::string("Claus"));
  Printer<English> printer2(std::string("John" ));
  printer1.print();
  printer2.print();

  return 0;
}
