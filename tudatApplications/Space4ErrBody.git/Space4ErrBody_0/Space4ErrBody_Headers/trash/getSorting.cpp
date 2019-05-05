#include "getStuff.h"
#include <utility>

namespace bislip {

bool getSorting( const std::tuple< double,double > &a, const std::tuple< double,double > &b )
{
    return ( std::get<1>(a) < std::get<1>(b));
}

}
