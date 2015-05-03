#include <exception>
#include <iostream>

struct Animal
{
  virtual void makeSound() { std::cout << "(...)" << std::endl; }
};

struct Dog : public Animal
{
  virtual void makeSound() { std::cout << "Woof!" << std::endl; }
};

struct Cat : public Animal
{
  virtual void makeSound() { std::cout << "Meaw!" << std::endl; }
};

struct ConfusedCat : public Cat
{
  virtual void makeSound() { std::cout << "Moooooh!" << std::endl; }
};

int main()
{
  try
  {
    Animal animal;
    Cat cat;
    ConfusedCat confusedCat;

    Animal& catReference = cat;

    // dynamic_cast<Animal&>(cat).makeSound(); // Meaw!
    // dynamic_cast<Cat&>(animal).makeSound(); // std::bad_cast
    // dynamic_cast<ConfusedCat&>(cat).makeSound(); // std::bad_cast
    // dynamic_cast<Cat&>(confusedCat).makeSound(); // Mooh!

    // catReference.makeSound(); // Meaw!
    // dynamic_cast<Cat&>(catReference).makeSound(); // Meaw!
    // dynamic_cast<ConfusedCat&>(catReference).makeSound(); // std::bad_cast
    dynamic_cast<Dog&>(catReference).makeSound(); // std::bad_cast
  }
  catch ( std::exception &e )
  {
    std::cout << e.what() << std::endl;
    return 1;
  }

  return 0;
}
