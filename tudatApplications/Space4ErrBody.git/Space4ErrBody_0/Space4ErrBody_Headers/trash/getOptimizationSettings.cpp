#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

std::vector< double > getOptimizationSettings ( const std::string &filename )
//void readinputfile ( std::vector< double > *input_data )
{
    std::ifstream inputdata;
    inputdata.open(filename.c_str());
    if(inputdata.fail())
    {
        std::cout << "File '" << filename << "' failed to open";
        std::cin.ignore(100, '\n');
        std::exit(1);
//        return 1;
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
        double var;
        std::vector< double > OptimizationSettings;
        while (!inputdata.fail() && !inputdata.eof())
        {
              inputdata >> var;
              OptimizationSettings.push_back(var);
        }

        inputdata.close();
        return OptimizationSettings;
    }

   //
}
