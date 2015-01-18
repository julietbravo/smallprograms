#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "Input.h"

namespace
{
    template<typename T>
    void checkItem(const T &t) {}
    
    template<>
    void checkItem(const std::string &s)
    {
        // Check whether string is empty or whether the first character is not alpha.
        if (s.empty())
            throw std::runtime_error("Illegal string");
        else if (!isalpha(s[0]))
            throw std::runtime_error("Illegal string: " + s);
    
        // Return string if all characters are alphanumeric.
        if (find_if(s.begin(), s.end(), [](const char c) { return !std::isalnum(c); }) == s.end())
            return;
        else
            throw std::runtime_error("Illegal string: " + s);
    }
}

Input::Input(const std::string& file_name)
{
    std::string blockname;

    // Read file and throw exception on error.
    std::ifstream infile;
    // infile.exceptions(std::ifstream::failbit);
    infile.open(file_name);
    if (!infile.good())
        throw std::runtime_error("Illegal file name");

    std::string line;

    while (std::getline(infile, line))
    {
        // Strip of the comments.
        std::vector<std::string> strings;
        boost::split(strings, line, boost::is_any_of("#"));

        // Keep part that is not comment.
        if (strings.size() >= 2)
            line = strings[0];

        // Strip of all the whitespace.
        boost::trim(line);

        // Split string on = character.
        strings.clear();
        boost::split(strings, line, boost::is_any_of("="));

        // In case of one element, check whether we have a block header.
        if (strings.size() == 1)
        {
            std::string header = strings[0];
            boost::trim(header);

            // If only an empty line remains, jump to the next line.
            if (header.empty())
                continue;

            // Store the block name and check for validity of string.
            if (header.front() == '[' && header.back() == ']')
            {
                blockname = header.substr(1, header.size()-2);
                checkItem(blockname);
            }
            else
                throw std::runtime_error("Illegal line");
        }
        // Read item.
        else if (strings.size() == 2)
        {
            if (blockname.empty())
                throw std::runtime_error("No block name found");

            std::string left = strings[0];
            std::string right = strings[1];
            boost::trim(left);
            boost::trim(right);

            // Check if suboptions are defined.
            std::string leftsub;
            std::size_t openpos = left.find_first_of("[");
            std::size_t closepos = left.find_last_of("]");
            if (openpos != std::string::npos && closepos != std::string::npos)
            {
                if (openpos < closepos)
                {
                    // Save the suboption and check string.
                    leftsub = left.substr(openpos+1, closepos-openpos-1);
                    checkItem(leftsub);

                    // Strip of suboption.
                    left = left.substr(0, openpos);
                }
            }

            checkItem(left);

            // Leave the checking of the right string for later
            // when the type is known.
            itemlist[blockname][left][leftsub] = right;
        }
        // Throw an error.
        else
            throw std::runtime_error("Illegal line");
    }
}

void Input::printItemList()
{
    // Print the list as a test.
    for (auto &b : itemlist)
        for (auto &i : b.second)
            for (auto &is : i.second)
                std::cout << b.first << "," << i.first << "," << is.first << "," << is.second << ";" << std::endl;
}
    
namespace
{
    std::string getItemString(const ItemList& itemlist,
                              const std::string& blockname,
                              const std::string& itemname,
                              const std::string& subitemname)
    {
        auto itblock = itemlist.find(blockname);
        if (itblock == itemlist.end())
            throw std::runtime_error("Block does not exist");
    
        auto ititem = itblock->second.find(itemname);
        if (ititem == itblock->second.end())
            throw std::runtime_error("Item does not exist");
    
        auto itsubitem = ititem->second.find(subitemname);
        if (itsubitem == ititem->second.end())
            throw std::runtime_error("Subitem does not exist");
    
        return itsubitem->second;
    }
    
    template<typename T>
    T getItemFromStream(std::istringstream &ss)
    {
        // Read the item from the stringstream, operator >> trims automatically.
        T item;
        if (!(ss >> item))
            throw std::runtime_error("Item does not match type");
    
        // Check whether stringstream is empty, if not type is incorrect.
        std::string dummy;
        if (ss >> dummy)
            throw std::runtime_error("Partial item does not match type");
    
        return item;
    }
}

template<typename T>
T Input::getItem(const std::string& blockname,
                 const std::string& itemname,
                 const std::string& subitemname)
{
    std::string value = getItemString(itemlist, blockname, itemname, subitemname);

    std::istringstream ss(value);

    T item = getItemFromStream<T>(ss);
    checkItem<T>(item);

    return item;
}

template<typename T>
std::vector<T> Input::getList(const std::string &blockname,
                              const std::string &itemname,
                              const std::string &subitemname)
{
    std::string value = getItemString(itemlist, blockname, itemname, subitemname);

    std::vector<std::string> listitems;
    boost::split(listitems, value, boost::is_any_of(","));

    std::vector<T> list;
    for (std::string itemstring : listitems)
    {
        std::istringstream ss(itemstring);
        T item = getItemFromStream<T>(ss);
        checkItem(item);
        list.push_back(item);
    }

    return list;
}

// Explicitly instantiate templates.
template int Input::getItem<int>(const std::string&, const std::string&, const std::string&);
template double Input::getItem<double>(const std::string&, const std::string&, const std::string&);
template std::string Input::getItem<std::string>(const std::string&, const std::string&, const std::string&);

template std::vector<int> Input::getList<int>(const std::string&, const std::string&, const std::string&);
template std::vector<double> Input::getList<double>(const std::string&, const std::string&, const std::string&);
template std::vector<std::string> Input::getList<std::string>(const std::string&, const std::string&, const std::string&);

