#include "Space4ErrBody.h"
#include "updateGuidance.h"
#include "getStuff.h"



namespace bislip {

double MyGuidance::getCurrentThrustMagnitude(  )
{
    //E_hat_ = getE_hat( currentHeight_ , currentAirspeed_ , E_max_ );

    std::cout << "getCurrentThrustMagnitude"  << std::endl;
    FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    double currentHeight = FlightConditions_->getCurrentAltitude() ;
    double currentAirspeed = FlightConditions_->getCurrentAirspeed() ;
   // std::cout << "2 " << std::endl;
    E_hat_ =  ( 9.80665 * currentHeight + 0.5 * currentAirspeed * currentAirspeed ) / E_max_;
    std::cout << "E_hat_: " << E_hat_  << std::endl;
    // E_hat_ = getE_hat( );
    // std::cout << "E_hat_: " << E_hat_  << std::endl;


    if ( throttle < 0.0 )
    {
        throttle = 0.0;
    }
    if ( throttle > 1.0 )
    {
        throttle = 1.0;
    }
    //std::cout << "maxThrust_: " << maxThrust_ << std::endl;

    currentThrustMagnitude_ = throttle * maxThrust_;
    std::cout << "throttle: " << throttle << std::endl;
    std::cout << "currentThrustMagnitude_: " << currentThrustMagnitude_ << std::endl;


    return currentThrustMagnitude_;

}
} // namespace aerodynamics
//} // namespace tudat
