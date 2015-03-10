#include <iostream>

struct Animal_non_virtual
{
    void print() { std::cout << "Animal" << std::endl; }
};

struct Animal_virtual
{
    virtual void print() { std::cout << "Animal" << std::endl; }
};

struct Animal_abstract
{
    virtual void print() = 0;
};

struct Monkey : public Animal_non_virtual
{
    void print() { std::cout << "Monkey" << std::endl; }
};

int main()
{
    Monkey monkey;
    monkey.print();

    return 0;
}

