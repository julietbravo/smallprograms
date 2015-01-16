#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <boost/algorithm/string.hpp>

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
  std::map< std::string, std::map<std::string, std::string> > itemlist;
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
      std::string left  = strings[0];
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

int main(int argc, char *argv[])
{
  try
  {
    readIniFile(argv);
  }
  catch (std::runtime_error &e)
  {
    std::cout << "EXCEPTION: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
