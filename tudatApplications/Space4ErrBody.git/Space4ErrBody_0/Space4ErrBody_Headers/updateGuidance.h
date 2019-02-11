//#include "getStuff.h"
#include <Tudat/Bislip/bislipVariables.h>

#include <limits>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>// std::pair, std::get
#include <boost/filesystem/operations.hpp>
#include <chrono>
#include <thread>

#include <Eigen/Core>

#include <pagmo/problem.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/island.hpp>
#include <pagmo/serialization.hpp>
#include <pagmo/types.hpp>

#include <boost/format.hpp>


#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Mathematics/RootFinders/secantRootFinder.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/nrlmsise00Atmosphere.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/Mathematics/Interpolators/createInterpolator.h>
#include <Tudat/Mathematics/Interpolators/interpolator.h>
#include <Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>
#include <Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h>
#include <Tudat/Astrodynamics/SystemModels/vehicleSystems.h>


#include <Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h>

namespace bislip {

class MyGuidance: public tudat::aerodynamics::AerodynamicGuidance
{
public:

    //! Constructor.
    MyGuidance(
            const tudat::simulation_setup::NamedBodyMap& bodyMap,
            const std::string vehicleName,
            const std::string centralBodyName ):
        bodyMap_( bodyMap ),
        vehicleName_( vehicleName ),
        centralBodyName_( centralBodyName ){ }

    //! Destructor.
    ~MyGuidance( ){ }

    void updateGuidance( const double currentTime );

    //double getE_hat( );

    //double getAoAforCmEqualToZero( const double AoA, const double Mach );

    //Eigen::Vector6d getCoefficients( const std::vector< double > &coefficient_input );

    Eigen::Vector6d getCurrentCoefficients( );

    Eigen::Vector3d getCurrentBodyFixedThrustDirection( );

    double getCurrentThrustMagnitude( );

    double getCurrentSpecificImpulse( )
    {
        return bodyMap_.at( vehicleName_ )->getBislipSystems()->getSpecificImpulse();
    }
    bool getCurrentEngineStatus( );

    double getBodyflapCmIncrement( );

    double getBodyFlapDeflection( );

    bool getReversal( const double &angularDistanceToGo_deg, const double &abs_chi_err_deg, const double &reversal_conditional );


public:
    tudat::simulation_setup::NamedBodyMap bodyMap_;
    std::string vehicleName_;
    std::string centralBodyName_;
    //std::vector< double > parameterBounds_Ascent_;
    //std::vector< double > parameterBounds_Descent_;
    //double finalMass_;
    //double E_hat_;
    //double Isp_;
    //double maxThrust_;
    // std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_alpha_deg_;
    // std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_eps_T_deg_;
    // std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator_throttle_;
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > FlightConditions_;// = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems_;
    std::shared_ptr< bislip::VehicleSystems > bislipSystems_;
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface_;

    //Eigen::Vector3d currentbodyFixedThrustDirection_;
    //double currentThrustMagnitude_;
    //double currentSpecificImpulse_;
    //bool currentEngineStatus_;

};

} // namespace bislip
