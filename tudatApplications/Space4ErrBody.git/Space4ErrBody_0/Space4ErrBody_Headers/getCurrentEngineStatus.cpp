#include "updateGuidance.h"


namespace bislip {

bool MyGuidance::getCurrentEngineStatus( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    return bislip::Variables::determineEngineStatus( bodyMap_.at( vehicleName_ )->getBodyMass(), bodyMap_.at( vehicleName_ )->getVehicleSystems()->getDryMass() );
}
} // namespace bislip
