#include "Space4ErrBody.h"
#include "updateGuidance.h"
#include "getStuff.h"

#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

namespace bislip {

Eigen::Vector3d MyGuidance::getCurrentBodyFixedThrustDirection(  )
{

    //E_hat_ = getE_hat( currentHeight_ ,  currentAirspeed_ , E_max_ );

    std::cout << "getCurrentBodyFixedThrustDirection"  << std::endl;
    // E_hat_ = bislip::getE_hat( FlightConditions_->getCurrentAltitude() ,  FlightConditions_->getCurrentAirspeed() , E_max_);
    std::cout << "E_hat_: " << E_hat_  << std::endl;

    std::cout << "eps_T: " << eps_T  << std::endl;

    if ( eps_T < parameterBounds_[ 4 ] )
    {
        eps_T = parameterBounds_[ 4 ];
    }
    if ( eps_T > parameterBounds_[ 5 ] )
    {
        eps_T = parameterBounds_[ 5 ];
    }

    currentEpsilon_ = tudat::unit_conversions::convertDegreesToRadians( eps_T );

    //! Simplified expressions becuase thrust azimuth is known to be zero. I.e. phi_T = 0
    bodyFixedThrustDirection_( 0 ) = std::cos( currentEpsilon_ );
    bodyFixedThrustDirection_( 1 ) = 0.0;
    bodyFixedThrustDirection_( 2 ) = std::sin( currentEpsilon_ );


    return bodyFixedThrustDirection_;

}
} // namespace aerodynamics
//} // namespace tudat
