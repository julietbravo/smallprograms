#include <vector>
#include <iostream>
#include <stdexcept>
#include "Input.h"

int main(int argc, char *argv[])
{
    try
    {
        if (argc != 2)
            throw std::runtime_error("Illegal number or arguments");

        std::string file_name = std::string(argv[1]);
        Input input(file_name);
        // input.printItemList();

        int itot = input.get_item<int>("grid", "itot");
        double xsize = input.get_item<double>("grid", "xsize");
        double zsize = input.get_item<double>("grid", "zsize");
        std::string swthermo = input.get_item<std::string>("thermo", "swthermo");
        std::vector<std::string> crosslist = input.get_list<std::string>("cross", "crosslist");
        std::vector<double> xy = input.get_list<double>("cross", "xy");
        double rndamp = input.get_item<double>("fields", "rndamp");
        double rndampb = input.get_item<double>("fields", "rndamp", "b");

        std::cout << "itot = " << itot  << std::endl;
        std::cout << "xsize = " << xsize << std::endl;
        std::cout << "zsize = " << zsize << std::endl;
        std::cout << "swthermo = " << swthermo << std::endl;
        std::cout << "crosslist = ";
        for (std::string &s : crosslist)
            std::cout << "\"" << s << "\"" << " ";
        std::cout << std::endl;

        std::cout << "xy = ";
        for (double &i : xy)
            std::cout << i << " ";
        std::cout << std::endl;

        std::cout << "rndamp = " << rndamp << std::endl;
        std::cout << "rndamp[b] = " << rndampb << std::endl;
    }
    catch (std::exception &e)
    {
        std::cout << "EXCEPTION: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
