#include "getStuff.h"

namespace bislip {
double getRoot_NewtRaph (
        const double &x_old,
        const double &f,
        const double &fp){


    return ( x_old - f / fp );
}

double getRoot_Halley (
        const double &x_old,
        const double &f,
        const double &fp,
        const double &f2p) {

    //!Eq. 2 from  http://sci-hub.tw/10.1002/zamm.19540340110
    return ( x_old - 2 * f * fp / ( 2 * fp * fp - f * f2p) );
}

double getRoot_Schroder (
        const double &x_old,
        const double &f,
        const double &fp,
        const double &f2p) {

    //! http://mathworld.wolfram.com/SchroedersMethod.html
    //! https://www-sciencedirect-com.tudelft.idm.oclc.org/science/article/pii/S0377042709006347#b20
    return ( x_old - ( f * fp ) / ( fp * fp - f * f2p ) );
}
}
