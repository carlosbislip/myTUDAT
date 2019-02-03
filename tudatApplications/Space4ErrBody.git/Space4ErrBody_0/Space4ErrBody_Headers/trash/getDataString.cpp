#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

namespace bislip {
std::vector< std::string > bislip::variables::getDataString ( const std::string &filename )
{
    std::ifstream inputdata;
    inputdata.open(filename.c_str());
    if(inputdata.fail())
    {
        std::cout << "File '" << filename << "' failed to open";
        std::cin.ignore(100, '\n');
        std::exit(1);
    }
    else{
        /*
         std::cout << "The following text had been read." << std::endl;
         std::cout << " " << std::endl;
         std::cout << inputdata.rdbuf();  // Prints content to screen.
         inputdata.seekg(std::ios_base::beg) ;
         std::cout << " " << std::endl;
         std::cout << "The program will now proceed to import the data." << std::endl;
        */
        std::string  var;
        std::vector< std::string > stuff;
        //while (!inputdata.fail() && !inputdata.eof())
        while(std::getline(inputdata, var))
        {
            // inputdata >> var;
            // stuff.push_back(var);

            if(var.size() > 0)
            {
                stuff.push_back(var);
            }
        }
        inputdata.close();
        return stuff;
    }
}
}
