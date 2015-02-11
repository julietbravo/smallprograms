#include <exception>
#include <iostream>
#include <memory>

struct Animal
{
  virtual void makeSound() { std::cout << "(...)" << std::endl; }
};

struct Cat : public Animal
{
  virtual void makeSound() { std::cout << "Meaw!" << std::endl; }
};

struct Dog: public Animal
{
  virtual void makeSound() { std::cout << "Woof!" << std::endl; }
};

struct ConfusedCat : public Cat
{
  virtual void makeSound() { std::cout << "Moooooh!" << std::endl; }
};

std::unique_ptr<Animal> getAnimal(const int i)
{
  if (i == 1)
    return std::unique_ptr<Cat>(new Cat());
  else if (i == 2)
    return std::unique_ptr<ConfusedCat>(new ConfusedCat());
  else if (i == 3)
    return std::unique_ptr<Dog>(new Dog());
  else
    return std::unique_ptr<Animal>(new Animal());
}

int main()
{
  try
  {
    auto animal0 = getAnimal(0);
    auto animal1 = getAnimal(1);
    auto animal2 = getAnimal(2);
    auto animal3 = getAnimal(3);

    animal0->makeSound();
    animal1->makeSound();
    animal2->makeSound();
    animal3->makeSound();
  }
  catch ( std::exception &e )
  {
    std::cout << e.what() << std::endl;
    return 1;
  }

  return 0;
}
