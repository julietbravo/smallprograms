#include <iostream>

template<class DerivedPrinter>
class Printer
{
  public:
    Printer(const std::string name) : name_(name) {}
    void print() { static_cast<DerivedPrinter &>(*this).printImplementation(); }

  protected:
    const std::string name_;
};

class GermanPrinter : public Printer<GermanPrinter>
{
  public:
    GermanPrinter(const std::string name) : Printer(name) {}
    void printImplementation() {
      std::cout << "Gutentag, " << name_ << "!" << std::endl; }
};

class EnglishPrinter : public Printer<EnglishPrinter>
{
  public:
    EnglishPrinter(const std::string name) : Printer(name) {}
    void printImplementation() {
      std::cout << "Good afternoon, " << name_ << "!" << std::endl; }
};

int main()
{
  GermanPrinter  printer1(std::string("Claus"));
  EnglishPrinter printer2(std::string("John" ));
  printer1.print();
  printer2.print();

  return 0;
}
