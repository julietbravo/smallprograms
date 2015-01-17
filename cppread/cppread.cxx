#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <boost/algorithm/string.hpp>

// Item list is a global variable in this test.
std::map< std::string, std::map<std::string, std::string> > itemlist;

void checkString(const std::string &s)
{
  // Check whether string is empty or whether the first character is not alpha.
  if (s.empty())
    throw std::runtime_error("Illegal string");
  else if (!isalpha(s[0]))
    throw std::runtime_error("Illegal string: " + s);

  // Return string if all characters are alphanumeric.
  if (find_if(s.begin(), s.end(), [](const char c) { return !isalnum(c); }) == s.end())
    return;
  else
    throw std::runtime_error("Illegal string: " + s);
}

void readIniFile(char *argv[])
{
  std::string blockname;

  std::ifstream infile(argv[1]);
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
      {
        continue;
      }

      // Store the block name.
      if (header.front() == '[' && header.back() == ']')
      {
        blockname = header.substr(1, header.size()-2);
        checkString(blockname);
      }
      else
      {
        throw std::runtime_error("Illegal line");
      }
    }
    // Read item
    else if (strings.size() == 2)
    {
      if (blockname.empty())
      {
        throw std::runtime_error("No block name found");
      }
      std::string left = strings[0];
      std::string right = strings[1];
      boost::trim(left);
      boost::trim(right);
      checkString(left);

      // Leave the checking of the right string for later
      // when the type is known.
      itemlist[blockname][left] = right;
    }
    // Throw an error.
    else
    {
      throw std::runtime_error("Illegal line");
    }
  }

  // Print the list as a test
  for (auto &m : itemlist) {
    for (auto &s : m.second) {
      std::cout << m.first << "," << s.first << "," << s.second << ";" << std::endl;
    }
  }
}

std::string getItemString(const std::string &blockname, const std::string &itemname)
{
  auto itblock = itemlist.find(blockname);
  if (itblock == itemlist.end())
    throw std::runtime_error("Block does not exist");

  auto ititem = itblock->second.find(itemname);
  if (ititem == itblock->second.end())
    throw std::runtime_error("Item does not exist");

  return ititem->second;
}

template<typename T>
T getItem(const std::string &blockname, const std::string &itemname)
{
  std::string value = getItemString(blockname, itemname);

  std::istringstream ss(value);

  T item;
  if (!(ss >> item))
    throw std::runtime_error("Item does not match type");

  // Check whether stringstream is empty, if not type is incorrect.
  std::string dummy;
  if (ss >> dummy)
    throw std::runtime_error("Partial item does not match type");

  return item;
}

template<typename T>
std::vector<T> getList(const std::string &blockname, const std::string &itemname)
{
  std::string value = getItemString(blockname, itemname);

  std::vector<std::string> listitems;
  boost::split(listitems, value, boost::is_any_of(","));

  std::vector<T> list;
  for (std::string itemstring : listitems)
  {
    T item;

    std::istringstream ss(itemstring);
    ss >> item;

    list.push_back(item);
  }

  return list;
}

int main(int argc, char *argv[])
{
  try
  {
    readIniFile(argv);
    int itot = getItem<int>("grid", "itot");
    double xsize = getItem<double>("grid", "xsize");
    double zsize = getItem<double>("grid", "zsize");
    std::string swthermo = getItem<std::string>("thermo", "swthermo");
    std::vector<std::string> crosslist = getList<std::string>("cross", "crosslist");
    std::vector<double> xy = getList<double>("cross", "xy");

    std::cout << "itot = " << itot  << std::endl;
    std::cout << "xsize = " << xsize << std::endl;
    std::cout << "zsize = " << zsize << std::endl;
    std::cout << "swthermo = " << swthermo << std::endl;
    std::cout << "crosslist = ";
    for (std::string &s : crosslist)
    {
      std::cout << s << " ";
    }
    std::cout << std::endl;

    std::cout << "xy = ";
    for (double &i : xy)
    {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  }
  catch (std::exception &e)
  {
    std::cout << "EXCEPTION: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
