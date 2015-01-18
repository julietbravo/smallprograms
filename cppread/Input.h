#ifndef INPUT
#include <map>

class Input
{
    public:
        Input(const std::string&);
        template<typename T> T getItem(const std::string&, const std::string&, const std::string& = "");
        template<typename T> std::vector<T> getList(const std::string&, const std::string&, const std::string& = "");
        void printItemList();

        typedef std::map< std::string, std::map< std::string, std::map<std::string, std::string> > > ItemList;

    private:
        ItemList itemlist;
};
#endif
