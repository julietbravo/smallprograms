#ifndef INPUT
#include <map>

typedef std::map< std::string, std::map< std::string, std::map<std::string, std::string> > > ItemList;

class Input
{
    public:
        Input(const std::string&);
        template<typename T> T getItem(const std::string&, const std::string&, const std::string& = "");
        template<typename T> std::vector<T> getList(const std::string&, const std::string&, const std::string& = "");
        void printItemList();
    private:
        ItemList itemlist;
};
#endif
