#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <boost/algorithm/string.hpp>

std::string checkString(const std::string &s)
{
  if (s.empty())
    throw std::runtime_error("Illegal string");
  else if (!isalpha(s[0]))
    throw std::runtime_error("Illegal string");

  if (find_if(s.begin(), s.end(), [](const char c) { return !isalnum(c); }) == s.end())
    return s;
  else
    throw std::runtime_error("Illegal string");
}

int main(int argc, char *argv[])
{
  std::map< std::string, std::map<std::string, std::string> > itemlist;
  std::string blockname;

  std::ifstream infile(argv[1]);
  std::string line;

  while (std::getline(infile, line))
  {
    // Strip of the comments
    std::vector<std::string> strings;
    boost::split(strings, line, boost::is_any_of("#"));

    // Keep part that is not comment
    if (strings.size() >= 2)
      line = strings[0];

    // Strip of all the whitespace.
    boost::trim(line);

    // Split string on = character.
    strings.clear();
    boost::split(strings, line, boost::is_any_of("="));

    // In case of no equal, check whether we have a block header.
    if (strings.size() == 1)
    {
      std::string header = strings[0];
      boost::trim(header);

      // If only an empty line remains, jump to the next line.
      if (header.empty())
      {
        continue;
      }
      if (header.front() == '[' && header.back() == ']')
      {
        blockname = checkString(header.substr(1, header.size()-2));
      }
      else
      {
        throw std::runtime_error("Illegal line");
      }
    }
    else if (strings.size() == 2)
    {
      if (blockname.empty())
      {
        throw std::runtime_error("No block name found");
      }
      std::string left  = strings[0];
      std::string right = strings[1];
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
      std::cout << m.first << ", " << s.first << ", " << s.second << std::endl;
    }
  }

  return 0;
}
