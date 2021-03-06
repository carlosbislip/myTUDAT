/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_VARIABLE_H
#define TUDAT_JSONINTERFACE_VARIABLE_H

#include "Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace propagators
{

// VariableType

//! Map of `VariableType`s string representations.
static std::map< VariableType, std::string > variableTypes =
{
    { independentVariable, "independent" },
    { cpuTimeVariable, "cpuTime" },
    { stateVariable, "state" },
    { dependentVariable, "dependent" },
    { stateTransitionMatrix, "stateTransitionMatrix" },
    { sensitivityMatrix, "sensitivityMatrix" }
};

//! `VariableType`s not supported by `json_interface`.
static std::vector< VariableType > unsupportedVariableTypes = { };

//! Convert `VariableType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const VariableType& variableType )
{
    jsonObject = json_interface::stringFromEnum( variableType, variableTypes );
}

//! Convert `json` to `VariableType`.
inline void from_json( const nlohmann::json& jsonObject, VariableType& variableType )
{
    variableType = json_interface::enumFromString( jsonObject, variableTypes );
}


// PropagationDependentVariables

//! Map of `PropagationDependentVariables` string representations.
static std::map< PropagationDependentVariables, std::string > dependentVariableTypes =
{
    { mach_number_dependent_variable, "machNumber" },
    { altitude_dependent_variable, "altitude" },
    { airspeed_dependent_variable, "airspeed" },
    { local_density_dependent_variable, "localDensity" },
    { relative_speed_dependent_variable, "relativeSpeed" },
    { relative_position_dependent_variable, "relativePosition" },
    { relative_distance_dependent_variable, "relativeDistance" },
    { relative_velocity_dependent_variable, "relativeVelocity" },
    { radiation_pressure_dependent_variable, "radiationPressure" },
    { total_acceleration_norm_dependent_variable, "totalAccelerationNorm" },
    { single_acceleration_norm_dependent_variable, "accelerationNorm" },
    { total_acceleration_dependent_variable, "totalAcceleration" },
    { single_acceleration_dependent_variable, "acceleration" },
    { aerodynamic_force_coefficients_dependent_variable, "aerodynamicForceCoefficients" },
    { aerodynamic_moment_coefficients_dependent_variable, "aerodynamicMomentCoefficients" },
    { rotation_matrix_to_body_fixed_frame_variable, "rotationMatrixToBodyFixedFrame" },
    { intermediate_aerodynamic_rotation_matrix_variable, "intermediateAerodynamicRotationMatrix" },
    { relative_body_aerodynamic_orientation_angle_variable, "relativeBodyAerodynamicOrientationAngle" },
    { body_fixed_airspeed_based_velocity_variable, "bodyFixedAirspeedBasedVelocity" },
    { total_aerodynamic_g_load_variable, "totalAerodynamicGLoad" },
    { stagnation_point_heat_flux_dependent_variable, "stagnationPointHeatFlux" },
    { local_temperature_dependent_variable, "localTemperature" },
    { geodetic_latitude_dependent_variable, "geodeticLatitude" },
    { control_surface_deflection_dependent_variable, "controlSurfaceDeflection" },
    { total_mass_rate_dependent_variables, "totalMassRates" },
    { lvlh_to_inertial_frame_rotation_dependent_variable, "lvlhToInertialFrameRotation" },
    { periapsis_altitude_dependent_variable, "periapsisAltitude" },
    { total_torque_norm_dependent_variable, "totalTorqueNorm" },
    { single_torque_norm_dependent_variable, "torqueNorm" },
    { total_torque_dependent_variable, "totalTorque" },
    { single_torque_dependent_variable, "torque" },
    { body_fixed_groundspeed_based_velocity_variable, "bodyFixedGroundspeedBasedVelocity" },
    { keplerian_state_dependent_variable, "keplerElements" },
    { modified_equinocial_state_dependent_variable, "modifiedEquinoctialElements" },
    { spherical_harmonic_acceleration_terms_dependent_variable, "sphericalHarmonicsAccelerationTerms" },
    { body_fixed_relative_cartesian_position, "bodyFixedRelativeCartesianPosition" },
    { body_fixed_relative_spherical_position, "bodyFixedRelativeSphericalPosition" },
    { total_gravity_field_variation_acceleration, "totalGravityFieldVariationAcceleration" },
    { single_gravity_field_variation_acceleration, "singleGravityFieldVariationAcceleration" },
    { single_gravity_field_variation_acceleration_terms, "singleGravityFieldVariationAccelerationTerms" },
    { acceleration_partial_wrt_body_translational_state, "accelerationPartialWrtBodyTranslationalState" },
    { local_dynamic_pressure_dependent_variable, "localDynamicPressure" },
    { local_aerodynamic_heat_rate_dependent_variable, "localAerodynamicHeatRate" },
    { specific_energy, "specificEnergy" },
    { normalized_specific_energy, "normalizedSpecificEnergy" },
    { evaluated_throttle_setting, "evaluatedThrottleSetting" },
    { evaluated_thrust_elevation_angle, "evaluatedThrustElevationAngle" },
    { evaluated_thrust_azimuth_angle, "evaluatedThrustAzimuthAngle" },
    { evaluated_angle_of_attack, "evaluatedAngleOfAttack" },
    { evaluated_bank_angle, "evaluatedBankAngle" },
    { current_mass, "currentMass" },
    { engine_status, "engineStatus" },
    { angular_distance_traveled, "angularDistanceTraveled" },
    { angular_distance_to_go, "angularDistanceToGo" },
    { heading_to_target, "headingToTarget" },
    { heading_error, "headingError" },
    { heat_flux_tauber, "heatFluxTauber" },
    { body_fixed_thrust_vector, "bodyFixedThrustVector" },
    { bending_moment, "bendingMoment" },
    { local_gravity, "localGravity" },
    { skip_suppression_limit, "skipSuppressionLimit" },
    { bodyflap_deflection_moment_coefficient_increment, "bodyFlapDeflectionMomentCoefficientIncrement" },
    { bodyflap_deflection_moment_coefficient_increment_dif, "bodyFlapDeflectionMomentCoefficientIncrementDif" },
    { body_fixed_total_load_vector, "bodyFixedTotalLoadVector" },
    { body_fixed_total_g_load_vector, "bodyFixedTotal_g_LoadVector" },
    { body_fixed_total_g_load_magnitude, "bodyFixedTotal_g_LoadMagnitude" },
    { body_fixed_aero_load_vector, "bodyFixedAeroLoadVector" },
    { bank_reversal_trigger, "bankReversalTrigger" },
    { heat_flux_chapman, "heatFluxChapman" },
    { passenger_fixed_total_g_load_vector, "passengerFixedTotal_g_LoadVector" },
    { commanded_throttle_setting, "commandedThrottleSetting" },
    { commanded_thrust_elevation_angle, "commandedThrustElevationAngle" },
    { commanded_thrust_azimuth_angle, "commandedThrustAzimuthAngle" },
    { commanded_angle_of_attack, "commandedAngleOfAttack" },
    { commanded_bank_angle, "commandedBankAngle" },
    { current_lift_magnitude, "currentLiftForce" },
    { current_heading_error_deadband, "currentHeadingErrorDeadBand" },
    { reversal_conditional, "reversalConditional" },
    { wall_temperature_chapman, "wallTemperatureChapman" },
    { wall_temperature_tauber_stagnation, "wallTemperatureTauberStagnation" },
    { wall_temperature_tauber_flatplate, "wallTemperatureTauberFlatPlate" },
    { heat_flux_tauber_stagnation, "heatFluxTauberStagnation" },
    { heat_flux_tauber_flatplate, "heatFluxTauberFlatPlate" },
    { temp_bank_angle, "tempBankAngle" },
    { cumulative_angular_distance_travelled, "cumulativeAngularDistanceTravelled" },
    { groundtrack_difference, "groundtrackDifference" },
    { time_of_flight, "timeOfFlight" },
    { flight_path_angle_rate, "flightPathAngleRate" },
    { cumulative_cartesian_distance_travelled, "cumulativeCartesianDistanceTravelled" },
    { thrust_force_magnitude, "thrustForceMagnitude" },
    { aerodynamic_frame_aerodynamic_load_vector, "aerodynamicLoadVector" },
    { aerodynamic_frame_total_load_vector, "totalLoadAerodynamicFrameVector" },
    { aerodynamic_frame_total_acceleration_vector, "totalAccelerationAerodynamicFrameVector" },
    { passenger_frame_total_load_vector, "totalLoadPassengerFrameVector" },
    { passenger_frame_total_acceleration_vector, "totalAccelerationPassengerFrameVector" },
    { passenger_frame_jerk_vector, "passengerFrameJerkVector" },
    { trajectory_phase, "trajectoryPhase" },
    { height_dependent_variable, "height" },
    { angular_distance_covered_ratio, "angularDistanceCoveredRatio" }



};

//! `PropagationDependentVariables` not supported by `json_interface`.
static std::vector< PropagationDependentVariables > unsupportedDependentVariableTypes =
{

};

//! Convert `PropagationDependentVariables` to `json`.
inline void to_json( nlohmann::json& jsonObject, const PropagationDependentVariables& dependentVariable )
{
    jsonObject = json_interface::stringFromEnum( dependentVariable, dependentVariableTypes );
}

//! Convert `json` to `PropagationDependentVariables`.
inline void from_json( const nlohmann::json& jsonObject, PropagationDependentVariables& dependentVariable )
{
    dependentVariable = json_interface::enumFromString( jsonObject, dependentVariableTypes );
}


// VariableSettings

//! Create a `json` object from a shared pointer to a `VariableSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< VariableSettings >& variableSettings );

//! Create a shared pointer to a `VariableSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< VariableSettings >& variableSettings );


// SingleDependentVariableSaveSettings

//! Create a `json` object from a shared pointer to a `SingleDependentVariableSaveSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< SingleDependentVariableSaveSettings >& dependentVariableSettings );

//! Create a shared pointer to a `SingleDependentVariableSaveSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< SingleDependentVariableSaveSettings >& dependentVariableSettings );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_VARIABLE_H
