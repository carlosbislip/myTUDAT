/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationOutput.h"

namespace tudat
{

namespace propagators
{

//! Get the vector representation of a rotation matrix.
Eigen::VectorXd getVectorRepresentationForRotationMatrix(
        const Eigen::Matrix3d& currentRotationMatrix )
{
    Eigen::VectorXd vectorRepresentation = Eigen::VectorXd( 9 );
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            vectorRepresentation( i * 3 + j ) = currentRotationMatrix( i, j );
        }
    }
    return vectorRepresentation;
}

//! Get the vector representation of a rotation matrix.
Eigen::VectorXd getVectorRepresentationForRotationMatrixFunction(
        const std::function< Eigen::Matrix3d( ) > rotationFunction )
{
    return getVectorRepresentationForRotationMatrix( rotationFunction( ) );
}

//! Get the vector representation of a quaternion.
Eigen::VectorXd getVectorRepresentationForRotationQuaternion(
        const std::function< Eigen::Quaterniond( ) > rotationFunction )
{
    return getVectorRepresentationForRotationMatrix( rotationFunction( ).toRotationMatrix( ) );
}

//! Get the 3x3 matrix representation from a vector with 9 entries
Eigen::Matrix3d getMatrixFromVectorRotationRepresentation(
        const Eigen::VectorXd vectorRepresentation )
{
    if( vectorRepresentation.rows( ) != 9 )
    {
        throw std::runtime_error( "Error when putting vector in matrix representation, size is incompatible" );
    }
    Eigen::Matrix3d currentRotationMatrix;
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            currentRotationMatrix( i, j ) = vectorRepresentation( i * 3 + j );
        }
    }
    return currentRotationMatrix;
}


//! Get the quaternion formulation of an orthonormal matrix, from input of a vector with 9 entries corresponding to matrix
//! entries.
Eigen::Quaterniond getQuaternionFromVectorRotationRepresentation(
        const Eigen::VectorXd vectorRepresentation )
{
    return Eigen::Quaterniond( getMatrixFromVectorRotationRepresentation( vectorRepresentation ) );
}

//! Function to convert a matrix to the format used to save dependent variables
void getMatrixInOutputVectorRepresentation(
        const Eigen::MatrixXd& matrix, Eigen::VectorXd& vector )
{
    vector.setZero( matrix.rows( ) * matrix.cols( ) );
    for( int i = 0; i < matrix.rows( ); i++ )
    {
        vector.segment( i * matrix.cols( ), matrix.cols( ) ) =
                matrix.block( i, 0, 1, matrix.cols( ) ).transpose( );
    }
}

//! Function to convert a vector dependent variable output to its original matrix representation
void getOutputVectorInMatrixRepresentation(
        const Eigen::VectorXd& vector, Eigen::MatrixXd& matrix,
        const int rows, const int columns )
{
    if( rows * columns != vector.rows( ) )
    {
        throw std::runtime_error( "Error when getting matrix from output vector: sizes are incompatible" );
    }
    matrix.setZero( rows, columns );
    for( int i = 0; i < rows; i++ )
    {
        matrix.block( i, 0, 1, columns ) = vector.segment( i * columns, columns ).transpose( );
    }
}

//! Function to retrieve matrix block function output in vector representation
Eigen::VectorXd getVectorFunctionFromBlockFunction(
        const std::function< void( Eigen::Block< Eigen::MatrixXd > ) > blockFunction,
        const int numberOfRows, const int numberOfColumns )
{
    Eigen::MatrixXd matrixEvaluation = Eigen::MatrixXd::Zero( numberOfRows, numberOfColumns );
    blockFunction( matrixEvaluation.block( 0, 0, numberOfRows, numberOfColumns ) );

    Eigen::VectorXd vectorEvaluation;
    getMatrixInOutputVectorRepresentation( matrixEvaluation, vectorEvaluation );

    return vectorEvaluation;
}

//! Function to compute the Fay-Riddell equilibrium heat flux from body properties
double computeEquilibriumFayRiddellHeatFluxFromProperties(
        const std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions,
        const std::shared_ptr< system_models::VehicleSystems > vehicleSystems )
{
    return aerodynamics::computeEquilibriumFayRiddellHeatFlux(
                flightConditions->getCurrentDensity( ), flightConditions->getCurrentAirspeed( ),
                flightConditions->getCurrentFreestreamTemperature( ), flightConditions->getCurrentMachNumber( ),
                vehicleSystems->getNoseRadius( ), vehicleSystems->getWallEmissivity( ) );
}


//! Function to return a vector containing only one value given by doubleFunction
Eigen::VectorXd getVectorFromDoubleFunction( const std::function< double( ) >& doubleFunction )
{
    Eigen::VectorXd vector( 1 );
    vector << doubleFunction( );
    return vector;
}

//! Function to evaluate a set of vector-returning functions and concatenate the results.
Eigen::VectorXd evaluateListOfVectorFunctions(
        const std::vector< std::pair< std::function< Eigen::VectorXd( ) >, int > > vectorFunctionList,
        const int totalSize )
{
    Eigen::VectorXd variableList = Eigen::VectorXd::Zero( totalSize );
    int currentIndex = 0;

    for( std::pair< std::function< Eigen::VectorXd( ) >, int > vectorFunction: vectorFunctionList )
    {
        variableList.segment( currentIndex, vectorFunction.second ) = vectorFunction.first( );
        currentIndex += vectorFunction.second;
    }

    // Check consistency with input
    if( currentIndex != totalSize )
    {
        std::string errorMessage = "Error when evaluating lists of functions, sizes are inconsistent: " +
                std::to_string( currentIndex ) + " and " +
                std::to_string( totalSize );
        throw std::runtime_error( errorMessage );
    }

    return variableList;
}

//! Funtion to get the size of a dependent variable save settings
int getDependentVariableSaveSize(
        const std::shared_ptr< SingleDependentVariableSaveSettings >& singleDependentVariableSaveSettings )
{
    if ( singleDependentVariableSaveSettings->componentIndex_ >= 0 )
    {
        return 1;
    }
    else
    {
        return getDependentVariableSize(  singleDependentVariableSaveSettings );
    }
}

//! Funtion to get the size of a dependent variable
int getDependentVariableSize(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{
    int variableSize = -1;
    switch( dependentVariableSettings->dependentVariableType_ )
    {
    case mach_number_dependent_variable:
        variableSize = 1;
        break;
    case altitude_dependent_variable:
        variableSize = 1;
        break;
    case airspeed_dependent_variable:
        variableSize = 1;
        break;
    case local_density_dependent_variable:
        variableSize = 1;
        break;
    case relative_speed_dependent_variable:
        variableSize = 1;
        break;
    case relative_position_dependent_variable:
        variableSize = 3;
        break;
    case relative_distance_dependent_variable:
        variableSize = 1;
        break;
    case relative_velocity_dependent_variable:
        variableSize = 3;
        break;
    case radiation_pressure_dependent_variable:
        variableSize = 1;
        break;
    case total_acceleration_norm_dependent_variable:
        variableSize = 1;
        break;
    case single_acceleration_norm_dependent_variable:
        variableSize = 1;
        break;
    case total_acceleration_dependent_variable:
        variableSize = 3;
        break;
    case single_acceleration_dependent_variable:
        variableSize = 3;
        break;
    case aerodynamic_force_coefficients_dependent_variable:
        variableSize = 3;
        break;
    case aerodynamic_moment_coefficients_dependent_variable:
        variableSize = 3;
        break;
    case rotation_matrix_to_body_fixed_frame_variable:
        variableSize = 9;
        break;
    case intermediate_aerodynamic_rotation_matrix_variable:
        variableSize = 9;
        break;
    case relative_body_aerodynamic_orientation_angle_variable:
        variableSize = 1;
        break;
    case body_fixed_airspeed_based_velocity_variable:
        variableSize = 3;
        break;
    case body_fixed_groundspeed_based_velocity_variable:
        variableSize = 3;
        break;
    case total_aerodynamic_g_load_variable:
        variableSize = 1;
        break;
    case stagnation_point_heat_flux_dependent_variable:
        variableSize = 1;
        break;
    case local_temperature_dependent_variable:
        variableSize = 1;
        break;
    case local_dynamic_pressure_dependent_variable:
        variableSize = 1;
        break;
    case local_aerodynamic_heat_rate_dependent_variable:
        variableSize = 1;
        break;
    case geodetic_latitude_dependent_variable:
        variableSize = 1;
        break;
    case control_surface_deflection_dependent_variable:
        variableSize = 1;
        break;
    case total_mass_rate_dependent_variables:
        variableSize = 1;
        break;
    case lvlh_to_inertial_frame_rotation_dependent_variable:
        variableSize = 9;
        break;
    case periapsis_altitude_dependent_variable:
        variableSize = 1;
        break;
    case total_torque_dependent_variable:
        variableSize = 3;
        break;
    case single_torque_dependent_variable:
        variableSize = 3;
        break;
    case total_torque_norm_dependent_variable:
        variableSize = 1;
        break;
    case single_torque_norm_dependent_variable:
        variableSize = 3;
        break;
    case keplerian_state_dependent_variable:
        variableSize = 6;
        break;
    case spherical_harmonic_acceleration_terms_dependent_variable:
    {
        if( std::dynamic_pointer_cast< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                    dependentVariableSettings ) == nullptr )
        {
            std::string errorMessage = "Error, input for spherical_harmonic_acceleration_terms_dependent_variable inconsistent when getting parameter size ";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            variableSize = 3 * std::dynamic_pointer_cast< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                        dependentVariableSettings )->componentIndices_.size( );
        }
        break;
    }
    case modified_equinocial_state_dependent_variable:
        variableSize = 6;
        break;
    case body_fixed_relative_cartesian_position:
        variableSize = 3;
        break;
    case body_fixed_relative_spherical_position:
        variableSize = 3;
        break;
    case euler_angles_to_body_fixed_313:
        variableSize = 3;
        break;
    case total_gravity_field_variation_acceleration:
        variableSize = 3;
        break;
    case single_gravity_field_variation_acceleration:
        variableSize = 3;
        break;
    case single_gravity_field_variation_acceleration_terms:
    {
        if( std::dynamic_pointer_cast< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                    dependentVariableSettings ) == nullptr )
        {
            std::string errorMessage = "Error, input for single_gravity_field_variation_acceleration_terms inconsistent when getting parameter size ";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            variableSize = 3 * std::dynamic_pointer_cast< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                        dependentVariableSettings )->componentIndices_.size( );
        }
        break;
    }
    case acceleration_partial_wrt_body_translational_state:
        variableSize = 18;
        break;
    case current_body_mass_dependent_variable:
        variableSize = 1;
        break;
    case specific_energy:
        variableSize = 1;
        break;
    case normalized_specific_energy:
        variableSize = 1;
        break;
    case evaluated_throttle_setting:
        variableSize = 1;
        break;
    case evaluated_thrust_elevation_angle:
        variableSize = 1;
        break;
    case evaluated_thrust_azimuth_angle:
        variableSize = 1;
        break;
    case evaluated_angle_of_attack:
        variableSize = 1;
        break;
    case evaluated_bank_angle:
        variableSize = 1;
        break;
    case commanded_throttle_setting:
        variableSize = 1;
        break;
    case commanded_thrust_elevation_angle:
        variableSize = 1;
        break;
    case commanded_thrust_azimuth_angle:
        variableSize = 1;
        break;
    case commanded_angle_of_attack:
        variableSize = 1;
        break;
    case commanded_bank_angle:
        variableSize = 1;
        break;
    case current_mass:
        variableSize = 1;
        break;
    case engine_status:
        variableSize = 1;
        break;
    case angular_distance_traveled:
        variableSize = 1;
        break;
    case angular_distance_to_go:
        variableSize = 1;
        break;
    case heading_to_target:
        variableSize = 1;
        break;
    case heading_error:
        variableSize = 1;
        break;
    case heat_flux_tauber:
        variableSize = 1;
        break;
    case body_fixed_thrust_vector:
        variableSize = 3;
        break;
    case bending_moment:
        variableSize = 1;
        break;
    case local_gravity:
        variableSize = 3;
        break;
    case skip_suppression_limit:
        variableSize = 1;
        break;
    case bodyflap_deflection_moment_coefficient_increment:
        variableSize = 1;
        break;
    case bodyflap_deflection_moment_coefficient_increment_dif:
        variableSize = 1;
        break;
    case body_fixed_total_load_vector:
        variableSize = 3;
        break;
    case body_fixed_total_g_load_vector:
        variableSize = 3;
        break;
    case body_fixed_total_g_load_magnitude:
        variableSize = 1;
        break;
    case body_fixed_aero_load_vector:
        variableSize = 3;
        break;
    case bank_reversal_trigger:
        variableSize = 1;
        break;
    case heat_flux_chapman:
        variableSize = 1;
        break;
    case passenger_fixed_total_g_load_vector:
        variableSize = 3;
        break;
    case current_lift_magnitude:
        variableSize = 1;
        break;
    case current_heading_error_deadband:
        variableSize = 1;
        break;
    case reversal_conditional:
        variableSize = 1;
        break;
    case wall_temperature_chapman:
        variableSize = 1;
        break;
    case wall_temperature_tauber_stagnation:
        variableSize = 1;
        break;
    case wall_temperature_tauber_flatplate:
        variableSize = 1;
        break;
    case heat_flux_tauber_stagnation:
        variableSize = 1;
        break;
    case heat_flux_tauber_flatplate:
        variableSize = 1;
        break;
    case temp_bank_angle:
        variableSize = 1;
        break;
    case cumulative_angular_distance_travelled:
        variableSize = 1;
        break;
    case groundtrack_difference:
        variableSize = 1;
        break;
    case time_of_flight:
        variableSize = 1;
        break;
    case flight_path_angle_rate:
        variableSize = 1;
        break;
    case cumulative_cartesian_distance_travelled:
        variableSize = 1;
        break;
    case thrust_force_magnitude:
        variableSize = 1;
        break;
    case speed_of_sound:
        variableSize = 1;
        break;
    case adiabatic_wall_temperature:
        variableSize = 1;
        break;
    case freestream_temperature:
        variableSize = 1;
        break;
    case current_drag_magnitude:
        variableSize = 1;
        break;
    case estimated_flight_path_angle:
        variableSize = 1;
        break;
    case aerodynamic_frame_aerodynamic_load_vector:
        variableSize = 3;
        break;
    case aerodynamic_frame_total_load_vector:
        variableSize = 3;
        break;
    case aerodynamic_frame_total_acceleration_vector:
        variableSize = 3;
        break;
    case passenger_frame_total_load_vector:
        variableSize = 3;
        break;
    case passenger_frame_total_acceleration_vector:
        variableSize = 3;
        break;
    case passenger_frame_jerk_vector:
        variableSize = 3;
        break;
    case trajectory_phase:
        variableSize = 1;
        break;
    case height_dependent_variable:
        variableSize = 1;
        break;
    case angular_distance_covered_ratio:
        variableSize = 1;
        break;


    default:
        std::string errorMessage = "Error, did not recognize dependent variable size of type: " +
                std::to_string( dependentVariableSettings->dependentVariableType_ );
        throw std::runtime_error( errorMessage );
    }
    return variableSize;
}

template std::pair< std::function< Eigen::VectorXd( ) >, std::map< int, std::string > > createDependentVariableListFunction< double, double >(
        const std::shared_ptr< DependentVariableSaveSettings > saveSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > >& stateDerivativeModels );

//template std::pair< std::function< Eigen::VectorXd( ) >, int > getVectorDependentVariableFunction< double, double >(
//        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
//        const simulation_setup::NamedBodyMap& bodyMap,
//        const std::unordered_map< IntegratedStateType,
//        std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > >& stateDerivativeModels );

//template std::function< double( ) > getDoubleDependentVariableFunction< double, double >(
//        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
//        const simulation_setup::NamedBodyMap& bodyMap,
//        const std::unordered_map< IntegratedStateType,
//        std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > >& stateDerivativeModels );

} // namespace propagators

} // namespace tudat
