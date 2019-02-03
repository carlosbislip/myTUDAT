#include "updateGuidance.h"


namespace bislip {

double MyGuidance::getCurrentThrustMagnitude(  )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    return bislip::Variables::computeThrustMagnitude( FlightConditions_, bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
}
} // namespace bislip
