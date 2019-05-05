#include "updateGuidance.h"

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

namespace bislip  {

Eigen::Vector6d MyGuidance::getCoefficients( const std::vector< double > &coefficient_input )
{
    Eigen::Vector6d newCoefficients;

    // Update and retrieve current aerodynamic coefficients
    coefficientInterface_->updateCurrentCoefficients( coefficient_input );
    //newCoefficients = coefficientInterface_->getCurrentAerodynamicCoefficients( );

    return newCoefficients = coefficientInterface_->getCurrentAerodynamicCoefficients( );
}

} // namespace bislip
