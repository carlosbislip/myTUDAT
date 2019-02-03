#ifndef BISLIPVARIABLES_H
#define BISLIPVARIABLES_H

#include <Tudat/Bislip/bislipVehicleSystems.h>
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace bislip {

namespace Variables {

unsigned int millis_since_midnight ( );

std::string getCurrentDateTime ( const bool useLocalTime = false );

std::vector< std::string > getDataString ( const std::string &filename );

std::vector< double > getDataNumeri ( const std::string &filename );

double computeAngularDistance (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f);

double computeHeadingToTarget (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f);

double computeHeadingError (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f,
        const double &heading);

double computeSpecificEnergy (
        const double &height,
        const double &airspeed);

double computeNormalizedSpecificEnergy (
        const double &height,
        const double &airspeed,
        const double &E_max);

std::vector< double > HermiteDerivatives (
        const Eigen::VectorXd &mappedNormalizedSpecificEnergy,
        const Eigen::VectorXd &y);

std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > createOneDimensionalHermiteInterpolator (
        const Eigen::VectorXd &parameterValues,
        const Eigen::VectorXd &normalizedSpecificEnergy,
        const std::map< double, double > &mapped_data,
        const std::shared_ptr< tudat::interpolators::InterpolatorSettings > &interpolatorSettings );

bislip::Parameters::Optimization passOptimizationParameter (
        const bislip::Parameters::Optimization &parameter);

std::string passString (
      const std::string &string);

tudat::simulation_setup::NamedBodyMap& passBodyMap (
         tudat::simulation_setup::NamedBodyMap& bodyMap );

std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > chooseGuidanceInterpolator (
        const bislip::Parameters::Optimization &parameter,
        const std::shared_ptr< bislip::VehicleSystems > &bislipSystems);

std::pair < double, double > chooseGuidanceBounds (
        const bislip::Parameters::Optimization &parameter,
        const std::shared_ptr< bislip::VehicleSystems > &bislipSystems);

double evaluateGuidanceInterpolator (
        const bislip::Parameters::Optimization &parameter,
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

Eigen::Vector6d computeCurrentCoefficients (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

Eigen::Vector3d computeBodyFixedThrustDirection (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeThrustMagnitude (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName);

Eigen::Vector3d computeBodyFixedThrustVector (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

bool determineEngineStatus (
        const double &currentMass,
        const double &landingMass);

Eigen::Vector2d computeLocalGravity (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName );

double computeEquilibriumGlideLimit (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName);

double computeHeatingRate (
        const double &airdensity,
        const double &airspeed,
        const double &C,
        const double &N,
        const double &M);

double computeStagnationHeatFlux (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems);

double computeFlatPlateHeatFlux (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
        const std::shared_ptr< bislip::VehicleSystems > &bislipSystems);

double computeStagnationHeat (
        const double &airdensity,
        const double &airspeed,
        const double &C_s,
        const double &N,
        const double &M,
        const double &adiabaticWallTemperature,
        const double &WallTemperature);

double computeFlatPlateHeat (
        const double &airdensity,
        const double &airspeed,
        const double &C_FP_1,
        const double &C_FP_2,
        const double &adiabaticWallTemperature,
        const double &WallTemperature);

double computeHeatingRateTauber (
        const double &q_dot_s,
        const double &q_dot_FP,
        const double &lambda);

double computeBendingMoment (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computePenalty (
        const Eigen::VectorXd &dependentVariable_TimeHistory,
        const long &startIterator,
        const long &endIterator,
        const double &constraint, const double &fixedStepSize, const double &tof, const bool &direct );

double computeBodyFlapCmIncrement (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeBodyFlapDeflection( const double &bodyFlapCmIncrement );

//bool StopOrNot (const tudat::simulation_setup::NamedBodyMap& bodyMap,
//                const std::string &vehicleName,
//                const std::vector< double > &vehicleParameterValues,
//                const std::vector< double > &terminationConditionsValues);


}; // namespace Variables
} // namespace bislip

#endif // BISLIPVARIABLES_H
