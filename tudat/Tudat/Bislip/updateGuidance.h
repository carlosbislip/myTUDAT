#ifndef UPDATEGUIDANCE_H
#define UPDATEGUIDANCE_H


//#include "getStuff.h"
#include "Tudat/Bislip/bislipVariables.h"
//#include "../../../../tudat/Tudat/Bislip/bislipVariables.h"
/*
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

*/
/*
#include <pagmo/problem.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/island.hpp>
#include <pagmo/serialization.hpp>
#include <pagmo/types.hpp>
*/
//#include "bislipUtilities.h"

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
#include <Tudat/External/SpiceInterface/spiceInterface.h>

#include <Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h>

namespace bislip {

class MyGuidance: public tudat::aerodynamics::AerodynamicGuidance
{
public:

    //! Constructor.
    MyGuidance(
            const tudat::simulation_setup::NamedBodyMap& bodyMap,
            const std::string vehicleName,
            const std::string centralBodyName,
            const double startingEpoch ):
        bodyMap_( bodyMap ),
        vehicleName_( vehicleName ),
        centralBodyName_( centralBodyName ),
        startingEpoch_( startingEpoch ){ }

    //! Destructor.
    ~MyGuidance( ){ }

    void updateGuidance( const double currentTime );

    void evaluateGuidanceFunctions(
           // std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
            std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems,
            std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
            const std::string &currentTrajectoryPhase,
            const double &timeOfFLight );

    void evaluateValidationGuidanceFunctions(
            // std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
            std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems,
            std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
            const std::string &currentTrajectoryPhase,
            const double &timeOfFLight );


   // Eigen::Vector6d getPartialCurrentCoefficients( );
    Eigen::Vector3d getCurrentBodyFixedThrustDirection( );
    double getCurrentThrustMagnitude( );
    double getCurrentSpecificImpulse( ) { return bodyMap_.at( vehicleName_ )->getBislipSystems()->getSpecificImpulse(); }
    bool getCurrentEngineStatus( );

    double getBodyFlapDeflection( );


public:
    tudat::simulation_setup::NamedBodyMap bodyMap_;
    std::string vehicleName_;
    std::string centralBodyName_;
    double startingEpoch_;
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > FlightConditions_;
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems_;
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems_;
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface_;

};

} // namespace bislip



#endif // UPDATEGUIDANCE_H

