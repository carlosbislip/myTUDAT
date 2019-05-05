#include "updateGuidance.h"

namespace bislip  {

Eigen::Vector6d MyGuidance::getCurrentCoefficients()
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    if( coefficientInterface_ == nullptr )
    {
        coefficientInterface_ = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                    bodyMap_.at( vehicleName_ )->getAerodynamicCoefficientInterface( ) );
    }

    if( bislipSystems_ == nullptr )
    {
        bislipSystems_ = std::dynamic_pointer_cast< bislip::BislipVehicleSystems >(
                    bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
    }

    return bislip::Variables::computeCurrentCoefficients( FlightConditions_, coefficientInterface_, bislipSystems_ );
}
} // namespace bislip
