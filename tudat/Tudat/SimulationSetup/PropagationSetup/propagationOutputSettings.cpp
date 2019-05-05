/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h"

namespace tudat
{

namespace propagators
{

//! Function to get a string representing a 'named identification' of a variable type
std::string getVariableName( const VariableType variableType )
{
    switch ( variableType )
    {
    case independentVariable:
        return "Independent variable ";
    case cpuTimeVariable:
        return "Cumulative computation time variable ";
    case stateVariable:
        return "Integrated state ";
    case dependentVariable:
        return "Dependent variable ";
    case stateTransitionMatrix:
        return "State transition matrix ";
    case sensitivityMatrix:
        return "Sensitivity matrix ";
    default:
        throw std::runtime_error( "Error, variable " +
                                  std::to_string( variableType ) +
                                  "not found when retrieving parameter name " );
    }
}

//! Function to get a string representing a 'named identification' of a variable
std::string getVariableId( const std::shared_ptr< VariableSettings > variableSettings )
{
    std::shared_ptr< SingleDependentVariableSaveSettings > singleDependentVariableSaveSettings =
            std::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variableSettings );
    if ( singleDependentVariableSaveSettings )
    {
        return getDependentVariableId( singleDependentVariableSaveSettings );
    }
    else
    {
        return getVariableName( variableSettings->variableType_ );
    }
}


//! Function to get a string representing a 'named identification' of a dependent variable type
std::string getDependentVariableName( const PropagationDependentVariables propagationDependentVariables )
{
    std::string variableName = "";
    switch( propagationDependentVariables )
    {
    case mach_number_dependent_variable:
        variableName = "Mach number ";
        break;
    case altitude_dependent_variable:
        variableName = "Altitude ";
        break;
    case airspeed_dependent_variable:
        variableName = "Airspeed ";
        break;
    case local_density_dependent_variable:
        variableName = "Density ";
        break;
    case relative_speed_dependent_variable:
        variableName = "Relative speed ";
        break;
    case relative_position_dependent_variable:
        variableName = "Relative position ";
        break;
    case relative_distance_dependent_variable:
        variableName = "Relative distance ";
        break;
    case relative_velocity_dependent_variable:
        variableName = "Relative velocity ";
        break;
    case radiation_pressure_dependent_variable:
        variableName = "Radiation pressure ";
        break;
    case total_acceleration_norm_dependent_variable:
        variableName = "Total acceleration norm ";
        break;
    case single_acceleration_norm_dependent_variable:
        variableName = "Single acceleration norm of type ";
        break;
    case total_acceleration_dependent_variable:
        variableName = "Total acceleration in inertial frame ";
        break;
    case single_acceleration_dependent_variable:
        variableName = "Single acceleration in inertial frame of type ";
        break;
    case aerodynamic_force_coefficients_dependent_variable:
        variableName = "Aerodynamic force coefficients ";
        break;
    case aerodynamic_moment_coefficients_dependent_variable:
        variableName = "Aerodynamic moment coefficients ";
        break;
    case rotation_matrix_to_body_fixed_frame_variable:
        variableName = "Rotation matrix to body-fixed frame ";
        break;
    case intermediate_aerodynamic_rotation_matrix_variable:
        variableName = "Rotation matrix from ";
        break;
    case relative_body_aerodynamic_orientation_angle_variable:
        variableName = "Body orientation angle ";
        break;
    case body_fixed_airspeed_based_velocity_variable:
        variableName = "Airspeed-based velocity ";
        break;
    case body_fixed_groundspeed_based_velocity_variable:
        variableName = "Groundspeed-based velocity ";
        break;
    case total_aerodynamic_g_load_variable:
        variableName = "Aerodynamic g-load ";
        break;
    case stagnation_point_heat_flux_dependent_variable:
        variableName = "Stagnation-point heat flux ";
        break;
    case local_temperature_dependent_variable:
        variableName = "Local freestream temperature ";
        break;
    case local_dynamic_pressure_dependent_variable:
        variableName = "Local dynamic pressure ";
        break;
    case local_aerodynamic_heat_rate_dependent_variable:
        variableName = "Local aerodynamic heat rate ";
        break;
    case geodetic_latitude_dependent_variable:
        variableName = "Geodetic latitude ";
        break;
    case control_surface_deflection_dependent_variable:
        variableName = "Control Surface Deflection ";
        break;
    case total_mass_rate_dependent_variables:
        variableName = "Body mass rate ";
        break;
    case lvlh_to_inertial_frame_rotation_dependent_variable:
        variableName = "LVLH to inertial frame rotation matrix ";
        break;
    case periapsis_altitude_dependent_variable:
        variableName = "Periapsis altitude ";
        break;
    case single_torque_dependent_variable:
        variableName = "Single torque in body-fixed frame of type ";
        break;
    case total_torque_dependent_variable:
        variableName = "Total torque in body-fixed frame ";
        break;
    case single_torque_norm_dependent_variable:
        variableName = "Single torque norm in body-fixed frame of type ";
        break;
    case total_torque_norm_dependent_variable:
        variableName = "Total torque norm in body-fixed frame ";
        break;
    case keplerian_state_dependent_variable:
        variableName = "Kepler elements ";
        break;
    case modified_equinocial_state_dependent_variable:
        variableName = "Modified equinoctial elements ";
        break;
    case spherical_harmonic_acceleration_terms_dependent_variable:
        variableName = "Spherical harmonic acceleration terms ";
        break;
    case body_fixed_relative_cartesian_position:
        variableName = "Body-fixed relative Cartesian position ";
        break;
    case body_fixed_relative_spherical_position:
        variableName = "Body-fixed relative spherical position ";
        break;
    case euler_angles_to_body_fixed_313:
        variableName = "313 Euler angles to body-fixed frame ";
        break;
    case total_gravity_field_variation_acceleration:
        variableName = "Total time-variable gravity field acceleration correction ";
        break;
    case single_gravity_field_variation_acceleration:
        variableName = "Single-source time-variable gravity field acceleration correction ";
        break;
    case single_gravity_field_variation_acceleration_terms:
        variableName = "Single-source time-variable gravity field per-term acceleration correction ";
        break;
    case acceleration_partial_wrt_body_translational_state:
        variableName = "Acceleration partial w.r.t body state ";
        break;
    case specific_energy:
        variableName = "Specific energy ";
        break;
    case normalized_specific_energy:
        variableName = "Normalized specific energy ";
        break;
    case evaluated_throttle_setting:
        variableName = "Evaluated throttle setting ";
        break;
    case evaluated_thrust_elevation_angle:
        variableName = "Evaluated thrust elevation angle ";
        break;
    case evaluated_thrust_azimuth_angle:
        variableName = "Evaluated thrust azimuth angle ";
        break;
    case evaluated_angle_of_attack:
        variableName = "Evaluated angle of attack ";
        break;
    case evaluated_bank_angle:
        variableName = "Evaluated bank angle ";
        break;
    case commanded_throttle_setting:
        variableName = "Commanded throttle setting ";
        break;
    case commanded_thrust_elevation_angle:
        variableName = "Commanded thrust elevation angle ";
        break;
    case commanded_thrust_azimuth_angle:
        variableName = "Commanded thrust azimuth angle ";
        break;
    case commanded_angle_of_attack:
        variableName = "Commanded angle of attack ";
        break;
    case commanded_bank_angle:
        variableName = "Commanded bank angle ";
        break;
    case current_mass:
        variableName = "Current mass ";
        break;
    case engine_status:
        variableName = "Engine status ";
        break;
    case angular_distance_traveled:
        variableName = "Angular distance traveled ";
        break;
    case angular_distance_to_go:
        variableName = "Angular distance to go ";
        break;
    case heading_to_target:
        variableName = "Heading to target ";
        break;
    case heading_error:
        variableName = "Heading error ";
        break;
    case heat_flux_tauber:
        variableName = "Heat Flux - Tauber ";
        break;
    case body_fixed_thrust_vector:
        variableName = "Body-fixed thrust vector ";
        break;
    case bending_moment:
        variableName = "Bending moment ";
        break;
    case local_gravity:
        variableName = "Local gravity ";
        break;
    case skip_suppression_limit:
        variableName = "Skip suppression limit ";
        break;
    case bodyflap_deflection_moment_coefficient_increment:
        variableName = "BodyFlap deflection moment coefficient increment";
        break;
    case bodyflap_deflection_moment_coefficient_increment_dif:
        variableName = "BodyFlap deflection moment coefficient increment - differenced";
        break;
    case body_fixed_total_load_vector:
        variableName = "Body-fixed total load vector ";
        break;
    case body_fixed_total_g_load_vector:
        variableName = "Body-fixed total g-load vector ";
        break;
    case body_fixed_total_g_load_magnitude:
        variableName = "Body-fixed total g-load magnitude ";
        break;
    case body_fixed_aero_load_vector:
        variableName = "Body-fixed aerodynamic load vector ";
        break;
    case bank_reversal_trigger:
        variableName = "Bank reversal trigger ";
        break;
    case heat_flux_chapman:
        variableName = "Heat Flux - Chapman ";
        break;
    case passenger_fixed_total_g_load_vector:
        variableName = "Passenger-fixed total g-load vector ";
        break;
    case current_lift_magnitude:
        variableName = "Current Lift Force ";
        break;
    case current_heading_error_deadband:
        variableName = "Current Heading Error Deadband ";
        break;
    case reversal_conditional:
        variableName = "Reversal Conditional ";
        break;
    case wall_temperature_chapman:
        variableName = "Wall Temperature - Chapman  ";
        break;
    case wall_temperature_tauber_stagnation:
        variableName = "Wall Temperature - Tauber - Stagnation ";
        break;
    case wall_temperature_tauber_flatplate:
        variableName = "Wall Temperature - Tauber Flat Plate ";
        break;
    case heat_flux_tauber_stagnation:
        variableName = "Heat Flux - Tauber - Stagnation ";
        break;
    case heat_flux_tauber_flatplate:
        variableName = "Heat Flux - Tauber - Flat Plate ";
        break;
    case temp_bank_angle:
        variableName = "Temporary Bank Angle ";
        break;
    case cumulative_angular_distance_travelled:
        variableName = "Groundtrack Covered ";
        break;
    case groundtrack_difference:
        variableName = "Groundtrack Difference ";
        break;
    case time_of_flight:
        variableName = "Time of Flight ";
        break;
    case flight_path_angle_rate:
        variableName = "Time of Flight ";
        break;
    case cumulative_cartesian_distance_travelled:
        variableName = "Cumulative Distance Travelled ";
        break;
    case thrust_force_magnitude:
        variableName = "Thrust Force Magnitude ";
        break;
    case speed_of_sound:
        variableName = "Speed of Sound ";
        break;
    case adiabatic_wall_temperature:
        variableName = "Adiabatic Wall Temperature ";
        break;


    default:
        std::string errorMessage = "Error, dependent variable " +
                std::to_string( propagationDependentVariables ) +
                "not found when retrieving parameter name ";
        throw std::runtime_error( errorMessage );
    }
    return variableName;
}


//! Function to get a string representing a 'named identification' of a dependent variable
std::string getDependentVariableId(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{
    std::string variableId = getDependentVariableName( dependentVariableSettings->dependentVariableType_ );

    if( ( dependentVariableSettings->dependentVariableType_ == single_acceleration_dependent_variable ) ||
            ( dependentVariableSettings->dependentVariableType_ == single_acceleration_norm_dependent_variable ) )
    {
        std::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationDependentVariableSettings =
                std::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >( dependentVariableSettings );
        if( accelerationDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable ID, input is inconsistent (acceleration type )" );
        }
        else
        {
            variableId += basic_astrodynamics::getAccelerationModelName(
                        accelerationDependentVariableSettings->accelerationModelType_ );
        }
    }
    else if( ( dependentVariableSettings->dependentVariableType_ == single_torque_dependent_variable ) ||
             ( dependentVariableSettings->dependentVariableType_ == single_torque_norm_dependent_variable ) )
    {
        std::shared_ptr< SingleTorqueDependentVariableSaveSettings > torqueDependentVariableSettings =
                std::dynamic_pointer_cast< SingleTorqueDependentVariableSaveSettings >( dependentVariableSettings );
        if( torqueDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable ID, input is inconsistent (torque type )" );
        }
        else
        {
            variableId += basic_astrodynamics::getTorqueModelName(
                        torqueDependentVariableSettings->torqueModelType_ );
        }
    }
    else if( dependentVariableSettings->dependentVariableType_ == intermediate_aerodynamic_rotation_matrix_variable )
    {
        std::shared_ptr< IntermediateAerodynamicRotationVariableSaveSettings > rotationDependentVariableSettings =
                std::dynamic_pointer_cast< IntermediateAerodynamicRotationVariableSaveSettings >( dependentVariableSettings );
        if( rotationDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable ID, input is inconsistent (rotation matrix)" );
        }
        else
        {
            variableId +=
                    reference_frames::getAerodynamicFrameName( rotationDependentVariableSettings->baseFrame_ ) + "to " +
                    reference_frames::getAerodynamicFrameName( rotationDependentVariableSettings->targetFrame_ );
        }
    }

    else if( dependentVariableSettings->dependentVariableType_ == relative_body_aerodynamic_orientation_angle_variable )
    {
        std::shared_ptr< BodyAerodynamicAngleVariableSaveSettings > angleDependentVariableSettings =
                std::dynamic_pointer_cast< BodyAerodynamicAngleVariableSaveSettings >( dependentVariableSettings );
        if( angleDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable ID, input is inconsistent (angle)" );
        }
        else
        {
            variableId +=
                    reference_frames::getAerodynamicAngleName( angleDependentVariableSettings->angle_ );
        }
    }

    if( ( dependentVariableSettings->dependentVariableType_ == single_acceleration_dependent_variable ) ||
            ( dependentVariableSettings->dependentVariableType_ == single_acceleration_norm_dependent_variable ) ||
            ( dependentVariableSettings->dependentVariableType_ == spherical_harmonic_acceleration_terms_dependent_variable ) )
    {
        variableId += ", acting on " + dependentVariableSettings->associatedBody_;
        if( dependentVariableSettings->secondaryBody_ != dependentVariableSettings->associatedBody_ )
        {
            variableId += ", exerted by " + dependentVariableSettings->secondaryBody_;
        }
    }
    else
    {
        variableId += "of " + dependentVariableSettings->associatedBody_;
        if( dependentVariableSettings->secondaryBody_ != "" )
        {
            variableId += " w.r.t. " + dependentVariableSettings->secondaryBody_;
        }

    }
    return variableId;
}

}

}
