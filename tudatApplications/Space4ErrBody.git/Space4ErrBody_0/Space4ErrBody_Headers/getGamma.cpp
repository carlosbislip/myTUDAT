#include "getStuff.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace bislip{
double getGamma ( const double A, const double B, const double phi ) {
    //std::cout << " A: " <<A<<" B: " <<B<<" phi: " <<phi << std::endl;

    double x_new = 0.1;
    double x_old = 0;
    double f, fp, f2p, f3p;
    double dif;
    int sgn = 1;

    do {

        f   =      pow( std::sin( x_new ) , 2 ) - A * std::sin( x_new + phi ) - B;
        fp  =      std::sin( 2 * x_new )        - A * std::cos( x_new + phi );
        f2p =  2 * std::cos( 2 * x_new )        + A * std::sin( x_new + phi );
        f3p = -4 * std::sin( 2 * x_new )        - A * std::cos( x_new + phi );

        x_old = x_new;
        x_new = bislip::getRoot_Halley( x_old, f, fp, f2p );

        if ( std::abs( x_new ) > tudat::mathematical_constants::PI / 2 )
        {
            if ( x_new < 0 )
            {
                sgn = -1;
            }
            x_new = sgn * tudat::mathematical_constants::PI / 2;
        }
        //std::cout <<" old: " << x_old<< " new: " << x_new  << std::endl;
        dif = x_new - x_old;

    } while( abs( dif ) > 1e-10 );

return std::abs( x_new );

}
}
