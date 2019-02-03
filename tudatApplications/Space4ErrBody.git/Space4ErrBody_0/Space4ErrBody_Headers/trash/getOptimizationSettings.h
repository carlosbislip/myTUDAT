#include <iostream>
#include <ctime>
#include <sstream>

std::vector< double > getOptimizationSettings ( const std::string &filename );

std::vector< double > getConditions ( const std::string &filename );

std::string getCurrentDateTime( const bool useLocalTime = false );
