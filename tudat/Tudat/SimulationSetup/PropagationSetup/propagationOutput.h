/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONOUTPUT_H
#define TUDAT_PROPAGATIONOUTPUT_H
#include "Tudat/Bislip/bislipVariables.h"

#include <functional>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/Ephemerides/frameManager.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"
//#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/Mathematics/BasicMathematics/rotationRepresentations.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate a function with two input variables (by reference) from function pointers
/*!
 *  Function to evaluate a function with two input variables (by reference) from function pointers that return these
 *  two input variables.
 *  \param functionToEvaluate Function that is to be evaluated with input from function pointers.
 *  \param firstInput Function returning first input to functionToEvaluate.
 *  \param secondInput Function returning second input to functionToEvaluate.
 *  \return Output from functionToEvaluate, using functions firstInput and secondInput as input.
 */
template< typename OutputType, typename InputType >
OutputType evaluateBivariateReferenceFunction(
        const std::function< OutputType( const InputType&, const InputType& ) > functionToEvaluate,
        const std::function< InputType( ) > firstInput,
        const std::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

template< typename OutputType, typename InputType >
OutputType evaluateReferenceFunction(
        const std::function< OutputType( const InputType& ) > functionToEvaluate,
        const std::function< InputType( ) > firstInput )
{
    return functionToEvaluate( firstInput( ) );
}
template< typename OutputType, typename FirstInputType, typename SecondInputType >
OutputType evaluateBivariateReferenceFunction(
        const std::function< OutputType( const FirstInputType&, const SecondInputType& ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}
template< typename OutputType, typename FirstInputType, typename SecondInputType, typename ThirdInputType >
OutputType evaluateTrivariateReferenceFunction(
        const std::function< OutputType( const FirstInputType&, const SecondInputType&, const ThirdInputType& ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput,
        const std::function< ThirdInputType( ) > thirdInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ), thirdInput( ) );
}
template< typename OutputType, typename FirstInputType, typename SecondInputType, typename ThirdInputType >
OutputType evaluateTrivariateMixedFunction1(
        const std::function< OutputType( const FirstInputType, const SecondInputType&, const ThirdInputType& ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput,
        const std::function< ThirdInputType( ) > thirdInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ), thirdInput( ) );
}
template< typename OutputType, typename FirstInputType, typename SecondInputType, typename ThirdInputType, typename FourthInputType >
OutputType evaluateTetravariateReferenceFunction(
        const std::function< OutputType( const FirstInputType&, const SecondInputType&, const ThirdInputType&, const FourthInputType& ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput,
        const std::function< ThirdInputType( ) > thirdInput,
        const std::function< FourthInputType( ) > fourthInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ), thirdInput( ), fourthInput( ) );
}

template< typename OutputType, typename FirstInputType, typename SecondInputType, typename ThirdInputType, typename FourthInputType, typename FifthInputType >
OutputType evaluatePentavariateReferenceFunction(
        const std::function< OutputType( const FirstInputType&, const SecondInputType&, const ThirdInputType&, const FourthInputType&, const FifthInputType& ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput,
        const std::function< ThirdInputType( ) > thirdInput,
        const std::function< FourthInputType( ) > fourthInput,
        const std::function< FifthInputType( ) > fifthInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ), thirdInput( ), fourthInput( ), fifthInput( ) );
}

template< typename OutputType, typename FirstInputType, typename SecondInputType, typename ThirdInputType, typename FourthInputType, typename FifthInputType, typename SixthInputType >
OutputType evaluateHexavariateReferenceFunction(
        const std::function< OutputType( const FirstInputType&, const SecondInputType&, const ThirdInputType&, const FourthInputType&, const FifthInputType&, const SixthInputType& ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput,
        const std::function< ThirdInputType( ) > thirdInput,
        const std::function< FourthInputType( ) > fourthInput,
        const std::function< FifthInputType( ) > fifthInput,
        const std::function< SixthInputType( ) > sixthInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ), thirdInput( ), fourthInput( ), fifthInput( ), sixthInput( ) );
}




//! Function to evaluate a function with two input variables from function pointers
/*!
 *  Function to evaluate a function with two input variables from function pointers that return these
 *  two input variables.
 *  \param functionToEvaluate Function that is to be evaluated with input from function pointers.
 *  \param firstInput Function returning first input to functionToEvaluate.
 *  \param secondInput Function returning second input to functionToEvaluate.
 *  \return Output from functionToEvaluate, using functions firstInput and secondInput as input.
 */
/*template< typename OutputType, typename InputType >
OutputType evaluateBivariateFunction(
        const std::function< OutputType( const InputType, const InputType ) > functionToEvaluate,
        const std::function< InputType( ) > firstInput,
        const std::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}
*/
template< typename OutputType, typename FirstInputType, typename SecondInputType>
OutputType evaluateBivariateFunction(
        const std::function< OutputType( const FirstInputType, const SecondInputType ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

//! Function to evaluate a function with three input variables from function pointers
/*!
 *  Function to evaluate a function with three input variables from function pointers that return these
 *  three input variables.
 *  \param functionToEvaluate Function that is to be evaluated with input from function pointers.
 *  \param firstInput Function returning first input to functionToEvaluate.
 *  \param secondInput Function returning second input to functionToEvaluate.
 *  \param thirdInput Function returning third input to functionToEvaluate.
 *  \return Output from functionToEvaluate, using functions firstInput, secondInput and thirdInput as input.
 */
template< typename OutputType, typename FirstInputType, typename SecondInputType, typename ThirdInputType >
OutputType evaluateTrivariateFunction(
        const std::function< OutputType( const FirstInputType, const SecondInputType, const ThirdInputType ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput,
        const std::function< ThirdInputType( ) > thirdInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ), thirdInput( ) );
}
template< typename OutputType, typename FirstInputType, typename SecondInputType, typename ThirdInputType, typename FourthInputType >
OutputType evaluateTetravariateFunction(
        const std::function< OutputType( const FirstInputType, const SecondInputType, const ThirdInputType, const FourthInputType ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput,
        const std::function< ThirdInputType( ) > thirdInput,
        const std::function< FourthInputType( ) > fourthInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ), thirdInput( ), fourthInput( ) );
}

template< typename OutputType, typename FirstInputType, typename SecondInputType, typename ThirdInputType, typename FourthInputType, typename FifthInputType >
OutputType evaluatePentavariateFunction(
        const std::function< OutputType( const FirstInputType, const SecondInputType, const ThirdInputType, const FourthInputType, const FifthInputType ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput,
        const std::function< ThirdInputType( ) > thirdInput,
        const std::function< FourthInputType( ) > fourthInput,
        const std::function< FifthInputType( ) > fifthInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ), thirdInput( ), fourthInput( ), fifthInput( ) );
}

template< typename OutputType, typename FirstInputType, typename SecondInputType, typename ThirdInputType, typename FourthInputType, typename FifthInputType, typename SixthInputType >
OutputType evaluateHexavariateFunction(
        const std::function< OutputType( const FirstInputType, const SecondInputType, const ThirdInputType, const FourthInputType, const FifthInputType, const SixthInputType ) > functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput,
        const std::function< ThirdInputType( ) > thirdInput,
        const std::function< FourthInputType( ) > fourthInput,
        const std::function< FifthInputType( ) > fifthInput,
        const std::function< SixthInputType( ) > sixthInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ), thirdInput( ), fourthInput( ), fifthInput( ), sixthInput( ) );
}


//! Funtion to get the size of a dependent variable save settings
/*!
 * Funtion to get the size of a dependent variable save settings.
 * \param singleDependentVariableSaveSettings Save settings for a dependent variable.
 * \return Size of requested dependent variable to save (equal to the size of the associated dependent variable,
 * or equal to 1 if the property `component_` is set).
 */
int getDependentVariableSaveSize(
        const std::shared_ptr< SingleDependentVariableSaveSettings >& singleDependentVariableSaveSettings );

//! Funtion to get the size of a dependent variable
/*!
 * Funtion to get the size (i.e. number of values in variable: one for altitude, three for position, etc.)
 * of a dependent variable
 * \param dependentVariableSettings Dependent variable type for which size is to be determined.
 * \return Size of requested dependent variable.
 */
int getDependentVariableSize(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings );

//! Get the vector representation of a rotation matrix.
/*!
 *  Get the vector representation of a rotation matrix.
 *  \param currentRotationMatrix Rotation matrix that is to be put into vector rerpesentation
 *  \return Column vector consisting of transpose of concatenated rows of currentRotationMatrix input.
 */
Eigen::VectorXd getVectorRepresentationForRotationMatrix(
        const Eigen::Matrix3d& currentRotationMatrix );

//! Get the vector representation of a rotation matrix.
/*!
 *  Get the vector representation of a rotation matrix.
 *  \param rotationFunction Function returning the rotation matrix that is to be put into vector rerpesentation
 *  \return Column vector consisting of transpose of concatenated rows of rotationFunction input.
 */
Eigen::VectorXd getVectorRepresentationForRotationMatrixFunction(
        const std::function< Eigen::Matrix3d( ) > rotationFunction );

//! Get the vector representation of a quaternion.
/*!
 *  Get the vector representation of a quaternion. Quaternion is converted to a rotation matrix, which is then put into
 *  a vector representation.
 *  \param rotationFunction Function returning the quaternion that is to be put inot vector rerpesentation
 *  \return Column vector consisting of transpose of concatenated rows of matrix representation of rotationFunction input.
 */
Eigen::VectorXd getVectorRepresentationForRotationQuaternion(
        const std::function< Eigen::Quaterniond( ) > rotationFunction );


//! Get the 3x3 matrix representation from a vector with 9 entries
/*!
 *  Get the matrix representation from a vector with 9 entries. The vector is the transpose of the concatenated rows of
 *  the associated matrix.
 *  \param vectorRepresentation Vector representation of 3x3 matrix (transpose of the concatenated rows of
 *  the associated matrix)
 *  \return Matrix from 3x3 vector representation.
 */
Eigen::Matrix3d getMatrixFromVectorRotationRepresentation(
        const Eigen::VectorXd vectorRepresentation );

//! Get the quaternion formulation of an orthonormal matrix, from input of a vector with 9 entries corresponding to matrix
//! entries.
/*!
 *  et the quaternion formulation of an orthonormal matrix, from input of a vector with 9 entries corresponding to matrix
 *  entries. The input vector is the transpose of the concatenated rows of the associated (orthonormal) matrix.
 *  \param vectorRepresentation Vector representation of 3x3 matrix (transpose of the concatenated rows of
 *  the associated matrix)
 *  \return Quaternion representation of orthonormal matrix obtained from 3x3 vector representation.
 */
Eigen::Quaterniond getQuaternionFromVectorRotationRepresentation(
        const Eigen::VectorXd vectorRepresentation );

//! Function to convert a matrix to the format used to save dependent variables
/*!
 * Function to convert a matrix to the format used to save dependent variables
 * \param matrix Matrix that is to be converted
 * \param vector Vector storage format of matrix
 */
void getMatrixInOutputVectorRepresentation(
        const Eigen::MatrixXd& matrix, Eigen::VectorXd& vector );

//! Function to convert a vector dependent variable output to its original matrix representation
/*!
 *  Function to convert a vector dependent variable output to its original matrix representation
 *  \param vector Vector dependent variable output
 *  \param matrix Original matrix representation
 *  \param rows Number of rows in matrix output
 *  \param columns Number of columns in matrix output
 */
void getOutputVectorInMatrixRepresentation(
        const Eigen::VectorXd& vector, Eigen::MatrixXd& matrix,
        const int rows, const int columns );

//! Function to retrieve matrix block function output in vector representation
/*!
 *  Function to retrieve matrix block function output in vector representation
 * \param blockFunction Function that returns (by reference) a matrix block
 * \param numberOfRows Number of rows in matrix output
 * \param numberOfColumns Number of columns in matrix output
 * \return Block-matrix in vector representation
 */
Eigen::VectorXd getVectorFunctionFromBlockFunction(
        const std::function< void( Eigen::Block< Eigen::MatrixXd > ) > blockFunction,
        const int numberOfRows, const int numberOfColumns );

//! Function to compute the Fay-Riddell equilibrium heat flux from body properties
/*!
 * Function to compute the Fay-Riddell equilibrium heat flux from body properties
 * \param flightConditions Object describing the current atmospheric flight conditions of the vehicle
 * \param vehicleSystems Object describing the physical properties of the vehicle
 * \return Equilibrium heat flux according to Fay-Riddell model
 */
double computeEquilibriumFayRiddellHeatFluxFromProperties(
        const std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions,
        const std::shared_ptr< system_models::VehicleSystems > vehicleSystems );

//! Function to retrieve relevant spherical harmonic acceleration model for dependent variable setting
/*!
 *  Function to retrieve relevant spherical harmonic acceleration model for dependent variable setting
 *  \param dependentVariableSettings Settings for dependent variable, associatedBody_ defines body undergoing acceleration,
 *  secondaryBody_ body exerting acceleration
 *  \param stateDerivativeModels List of state derivative models from which acceleration is to be retrieved
 *  \return Relevant spherical harmonic acceleration model for dependent variable setting
 */
template< typename StateScalarType, typename TimeType >
std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel >
getSphericalHarmonicAccelerationForDependentVariables(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels )
{
    std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration;

    // Retrieve list of suitable acceleration models (size should be one to avoid ambiguities)s
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
            listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                dependentVariableSettings->associatedBody_,
                dependentVariableSettings->secondaryBody_, stateDerivativeModels,
                basic_astrodynamics::spherical_harmonic_gravity );


    // Check if third-body counterpart of acceleration is found
    if( listOfSuitableAccelerationModels.size( ) == 0 )
    {
        listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                    dependentVariableSettings->associatedBody_,
                    dependentVariableSettings->secondaryBody_,
                    stateDerivativeModels, basic_astrodynamics::getAssociatedThirdBodyAcceleration(
                        basic_astrodynamics::spherical_harmonic_gravity ) );
    }

    if( listOfSuitableAccelerationModels.size( ) != 1 )
    {
        std::string errorMessage = "Error when getting spherical harmonic acceleration components between bodies " +
                dependentVariableSettings->associatedBody_ + " and " +
                dependentVariableSettings->secondaryBody_ + " of type " +
                std::to_string(
                    basic_astrodynamics::spherical_harmonic_gravity ) +
                ", no such acceleration found";
        throw std::runtime_error( errorMessage );
    }
    else
    {
        sphericalHarmonicAcceleration =
                std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                    listOfSuitableAccelerationModels.at( 0 ) );
        if( sphericalHarmonicAcceleration == nullptr )
        {
            std::string errorMessage = "Error when getting spherical harmonic acceleration for dependent variable " +
                    dependentVariableSettings->associatedBody_ + " and " +
                    dependentVariableSettings->secondaryBody_ + " type is ionconsistent";
            throw std::runtime_error( errorMessage );
        }
    }

    return sphericalHarmonicAcceleration;
}


//! Function to create a function returning a requested dependent variable value (of type VectorXd).
/*!
 *  Function to create a function returning a requested dependent variable value (of type VectorXd), retrieved from
 *  environment and/or state derivative models.
 *  \param dependentVariableSettings Settings for dependent variable that is to be returned by function created here.
 *  \param bodyMap List of bodies to use in simulations (containing full environment).
 *  \param stateDerivativeModels List of state derivative models used in simulations (sorted by dynamics type as key).
 *  \param stateDerivativePartials List of state derivative partials used in simulations (sorted by dynamics type as key).
 *  \return Function returning requested dependent variable. NOTE: The environment and state derivative models need to
 *  be updated to current state and independent variable before computation is performed.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::function< Eigen::VectorXd( ) >, int > getVectorDependentVariableFunction(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels =
        std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ),
        const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials =
        std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >( ) )
{
    std::function< Eigen::VectorXd( ) > variableFunction;
    int parameterSize;

    // Retrieve base information on dependent variable
    PropagationDependentVariables dependentVariable = dependentVariableSettings->dependentVariableType_;
    const std::string& bodyWithProperty = dependentVariableSettings->associatedBody_;
    const std::string& secondaryBody = dependentVariableSettings->secondaryBody_;

    // Check dependent variable type and create function accordingly.
    switch( dependentVariable )
    {
    case relative_position_dependent_variable:
    {
        // Retrieve functions for positions of two bodies.
        std::function< Eigen::Vector3d( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                std::bind( &linear_algebra::computeVectorDifference< 3 >, std::placeholders::_1, std::placeholders::_2 );
        std::function< Eigen::Vector3d( ) > firstInput =
                std::bind( &simulation_setup::Body::getPosition, bodyMap.at( bodyWithProperty ) );
        std::function< Eigen::Vector3d( ) > secondInput =
                std::bind( &simulation_setup::Body::getPosition, bodyMap.at( secondaryBody ) );

        variableFunction = std::bind(
                    &evaluateBivariateReferenceFunction< Eigen::Vector3d, Eigen::Vector3d >,
                    functionToEvaluate, firstInput, secondInput );
        parameterSize = 3;
        break;
    }
    case relative_velocity_dependent_variable:
    {
        // Retrieve functions for velocities of two bodies.
        std::function< Eigen::Vector3d( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                std::bind( &linear_algebra::computeVectorDifference< 3 >, std::placeholders::_1, std::placeholders::_2 );
        std::function< Eigen::Vector3d( ) > firstInput =
                std::bind( &simulation_setup::Body::getVelocity, bodyMap.at( bodyWithProperty ) );
        std::function< Eigen::Vector3d( ) > secondInput =
                std::bind( &simulation_setup::Body::getVelocity, bodyMap.at( secondaryBody ) );

        variableFunction = std::bind(
                    &evaluateBivariateReferenceFunction< Eigen::Vector3d, Eigen::Vector3d >,
                    functionToEvaluate, firstInput, secondInput );
        parameterSize = 3;


        break;
    }
    case total_acceleration_dependent_variable:
    {
        // Retrieve model responsible for computing accelerations of requested bodies.
        std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                getTranslationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
        variableFunction =
                std::bind( &NBodyStateDerivative< StateScalarType, TimeType >::getTotalAccelerationForBody, nBodyModel,
                           bodyWithProperty );
        parameterSize = 3;


        break;
    }
    case single_acceleration_dependent_variable:
    {
        // Check input consistency.
        std::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationDependentVariableSettings =
                std::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >( dependentVariableSettings );
        if( accelerationDependentVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_acceleration_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            // Retrieve list of suitable acceleration models (size should be one to avoid ambiguities)
            std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                    listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                        accelerationDependentVariableSettings->associatedBody_,
                        accelerationDependentVariableSettings->secondaryBody_,
                        stateDerivativeModels, accelerationDependentVariableSettings->accelerationModelType_ );

            // Check if third-body counterpart of acceleration is found
            if( listOfSuitableAccelerationModels.size( ) == 0 && basic_astrodynamics::isAccelerationDirectGravitational(
                        accelerationDependentVariableSettings->accelerationModelType_ ) )
            {
                listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                            accelerationDependentVariableSettings->associatedBody_,
                            accelerationDependentVariableSettings->secondaryBody_,
                            stateDerivativeModels, basic_astrodynamics::getAssociatedThirdBodyAcceleration(
                                accelerationDependentVariableSettings->accelerationModelType_  ) );
            }

            if( listOfSuitableAccelerationModels.size( ) != 1 )
            {
                std::string errorMessage = "Error when getting acceleration between bodies " +
                        accelerationDependentVariableSettings->associatedBody_ + " and " +
                        accelerationDependentVariableSettings->secondaryBody_ + " of type " +
                        std::to_string(
                            accelerationDependentVariableSettings->accelerationModelType_ ) +
                        ", no such acceleration found VECTOR";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                //std::function< Eigen::Vector3d( ) > vectorFunction =
                variableFunction = std::bind( &basic_astrodynamics::AccelerationModel3d::getAcceleration,
                                              listOfSuitableAccelerationModels.at( 0 ) );
                parameterSize = 3;
            }
        }
        break;
    }
    case spherical_harmonic_acceleration_terms_dependent_variable:
    {
        // Check input consistency.
        std::shared_ptr< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings > accelerationComponentVariableSettings =
                std::dynamic_pointer_cast< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >( dependentVariableSettings );
        if( accelerationComponentVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_acceleration_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration =
                    getSphericalHarmonicAccelerationForDependentVariables(
                        accelerationComponentVariableSettings, stateDerivativeModels );

            sphericalHarmonicAcceleration->setSaveSphericalHarmonicTermsSeparately( true );
            variableFunction = std::bind(
                        &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getConcatenatedAccelerationComponents,
                        sphericalHarmonicAcceleration, accelerationComponentVariableSettings->componentIndices_ );
            parameterSize = 3 * accelerationComponentVariableSettings->componentIndices_.size( );
        }
        break;
    }
    case total_gravity_field_variation_acceleration:
    {

        std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration =
                getSphericalHarmonicAccelerationForDependentVariables(
                    dependentVariableSettings, stateDerivativeModels );
        std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
                std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                    bodyMap.at( dependentVariableSettings->secondaryBody_ )->getGravityFieldModel( ) );

        if( timeDependentGravityField == nullptr )
        {
            throw std::runtime_error( "Error when requesting save of gravity field variation acceleration, central body " +
                                      dependentVariableSettings->secondaryBody_ +
                                      " has no TimeDependentSphericalHarmonicsGravityField." );
        }
        else
        {
            std::function< Eigen::VectorXd( const Eigen::MatrixXd&, const Eigen::MatrixXd& ) > accelerationFunction =
                    std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getAccelerationWithAlternativeCoefficients,
                               sphericalHarmonicAcceleration, std::placeholders::_1, std::placeholders::_2 );
            std::function< Eigen::MatrixXd( ) > cosineCorrectionFunction =
                    std::bind( &gravitation::TimeDependentSphericalHarmonicsGravityField::getTotalCosineCoefficientCorrection,
                               timeDependentGravityField,
                               sphericalHarmonicAcceleration->getMaximumDegree( ),
                               sphericalHarmonicAcceleration->getMaximumOrder( ) );
            std::function< Eigen::MatrixXd( ) > sineCorrectionFunction =
                    std::bind( &gravitation::TimeDependentSphericalHarmonicsGravityField::getTotalSineCoefficientCorrection,
                               timeDependentGravityField,
                               sphericalHarmonicAcceleration->getMaximumDegree( ),
                               sphericalHarmonicAcceleration->getMaximumOrder( ) );

            variableFunction = std::bind( evaluateBivariateReferenceFunction< Eigen::VectorXd, Eigen::MatrixXd >,
                                          accelerationFunction, cosineCorrectionFunction, sineCorrectionFunction );

            parameterSize = 3;
        }

        break;
    }
    case single_gravity_field_variation_acceleration:
    {
        std::shared_ptr< SingleVariationSphericalHarmonicAccelerationSaveSettings > accelerationVariableSettings =
                std::dynamic_pointer_cast< SingleVariationSphericalHarmonicAccelerationSaveSettings >( dependentVariableSettings );
        if( accelerationVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_gravity_field_variation_acceleration";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration =
                    getSphericalHarmonicAccelerationForDependentVariables(
                        accelerationVariableSettings, stateDerivativeModels );
            std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
                    std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                        bodyMap.at( dependentVariableSettings->secondaryBody_ )->getGravityFieldModel( ) );

            if( timeDependentGravityField == nullptr )
            {
                throw std::runtime_error( "Error when requesting save of gravity field variation acceleration, central body " +
                                          dependentVariableSettings->secondaryBody_ +
                                          " has no TimeDependentSphericalHarmonicsGravityField." );
            }
            else
            {
                std::function< Eigen::VectorXd( const Eigen::MatrixXd&, const Eigen::MatrixXd& ) > accelerationFunction =
                        std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getAccelerationWithAlternativeCoefficients,
                                   sphericalHarmonicAcceleration, std::placeholders::_1, std::placeholders::_2 );

                std::shared_ptr< gravitation::GravityFieldVariations > gravityFieldVatiation =
                        timeDependentGravityField->getGravityFieldVariationsSet( )->getGravityFieldVariation(
                            accelerationVariableSettings->deformationType_,
                            accelerationVariableSettings->identifier_ ).second;

                std::function< Eigen::MatrixXd( ) > cosineCorrectionFunction =
                        std::bind( &gravitation::GravityFieldVariations::getLastCosineCorrection,
                                   gravityFieldVatiation  );
                std::function< Eigen::MatrixXd( ) > sineCorrectionFunction =
                        std::bind( &gravitation::GravityFieldVariations::getLastSineCorrection,
                                   gravityFieldVatiation );

                variableFunction = std::bind( evaluateBivariateReferenceFunction< Eigen::VectorXd, Eigen::MatrixXd >,
                                              accelerationFunction, cosineCorrectionFunction, sineCorrectionFunction );

            }
        }
        parameterSize = 3;

        break;
    }
    case single_gravity_field_variation_acceleration_terms:
    {
        std::shared_ptr< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings > accelerationVariableSettings =
                std::dynamic_pointer_cast< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >( dependentVariableSettings );
        if( accelerationVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_gravity_field_variation_acceleration_terms";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration =
                    getSphericalHarmonicAccelerationForDependentVariables(
                        accelerationVariableSettings, stateDerivativeModels );
            std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
                    std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                        bodyMap.at( dependentVariableSettings->secondaryBody_ )->getGravityFieldModel( ) );

            if( timeDependentGravityField == nullptr )
            {
                throw std::runtime_error( "Error when requesting save of gravity field variation acceleration, central body " +
                                          dependentVariableSettings->secondaryBody_ +
                                          " has no TimeDependentSphericalHarmonicsGravityField." );
            }
            else
            {
                std::function< Eigen::VectorXd( const Eigen::MatrixXd&, const Eigen::MatrixXd& ) > accelerationFunction =
                        std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::
                                   getAccelerationComponentsWithAlternativeCoefficients,
                                   sphericalHarmonicAcceleration, std::placeholders::_1, std::placeholders::_2, accelerationVariableSettings->componentIndices_ );

                std::shared_ptr< gravitation::GravityFieldVariations > gravityFieldVatiation =
                        timeDependentGravityField->getGravityFieldVariationsSet( )->getGravityFieldVariation(
                            accelerationVariableSettings->deformationType_,
                            accelerationVariableSettings->identifier_ ).second;

                std::function< Eigen::MatrixXd( ) > cosineCorrectionFunction =
                        std::bind( &gravitation::GravityFieldVariations::getLastCosineCorrection,
                                   gravityFieldVatiation );
                std::function< Eigen::MatrixXd( ) > sineCorrectionFunction =
                        std::bind( &gravitation::GravityFieldVariations::getLastSineCorrection,
                                   gravityFieldVatiation );

                variableFunction = std::bind( evaluateBivariateReferenceFunction< Eigen::VectorXd, Eigen::MatrixXd >,
                                              accelerationFunction, cosineCorrectionFunction, sineCorrectionFunction );

            }
            parameterSize = 3 * accelerationVariableSettings->componentIndices_.size( );

        }
        break;
    }
    case aerodynamic_force_coefficients_dependent_variable:
    {
        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
        {
            std::string errorMessage = "Error, no atmospheric flight conditions available when requesting density output of aerodynamic force coefficients " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        variableFunction = std::bind(
                    &aerodynamics::AerodynamicCoefficientInterface::getCurrentForceCoefficients,
                    std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )->getAerodynamicCoefficientInterface( ) );
        parameterSize = 3;

        break;
    }
    case aerodynamic_moment_coefficients_dependent_variable:
    {
        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
        {

            std::string errorMessage = "Error, no atmospheric flight conditions available when requesting density output of aerodynamic moment coefficients " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        variableFunction = std::bind(
                    &aerodynamics::AerodynamicCoefficientInterface::getCurrentMomentCoefficients,
                    std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )->getAerodynamicCoefficientInterface( ) );
        parameterSize = 3;

        break;
    }
    case rotation_matrix_to_body_fixed_frame_variable:
    {
        std::function< Eigen::Quaterniond( ) > rotationFunction =
                std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodyMap.at( bodyWithProperty ) );
        variableFunction = std::bind( &getVectorRepresentationForRotationQuaternion, rotationFunction );
        parameterSize = 9;
        break;
    }
    case intermediate_aerodynamic_rotation_matrix_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
        {
            std::string errorMessage= "Error, no flight conditions when creating dependent variable function of type intermediate_aerodynamic_rotation_matrix_variable";
            throw std::runtime_error( errorMessage );
        }

        // Check input consistency.
        std::shared_ptr< IntermediateAerodynamicRotationVariableSaveSettings >
                intermediateAerodynamicRotationVariableSaveSettings =
                std::dynamic_pointer_cast< IntermediateAerodynamicRotationVariableSaveSettings >(
                    dependentVariableSettings );
        if( intermediateAerodynamicRotationVariableSaveSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type intermediate_aerodynamic_rotation_matrix_variable";
            throw std::runtime_error( errorMessage );
        }

        std::function< Eigen::Quaterniond( ) > rotationFunction =
                std::bind( &reference_frames::AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                           bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                           intermediateAerodynamicRotationVariableSaveSettings->baseFrame_,
                           intermediateAerodynamicRotationVariableSaveSettings->targetFrame_ );

        variableFunction = std::bind( &getVectorRepresentationForRotationQuaternion, rotationFunction );
        parameterSize = 9;
        break;
    }
    case body_fixed_airspeed_based_velocity_variable:
    {
        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
        {
            std::string errorMessage= "Error, no atmospheric flight conditions when creating dependent variable function of type body_fixed_airspeed_based_velocity_variable";
            throw std::runtime_error( errorMessage );
        }

        variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeedBasedVelocity,
                                      std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                          bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
        parameterSize = 3;
        break;
    }
    case body_fixed_groundspeed_based_velocity_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
        {
            std::string errorMessage= "Error, no flight conditions when creating dependent variable function of type body_fixed_groundspeed_based_velocity_variable";
            throw std::runtime_error( errorMessage );
        }

        if(  bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ) == nullptr )
        {
            std::string errorMessage= "Error, no aerodynamic angle calculator when creating dependent variable function of type body_fixed_groundspeed_based_velocity_variable";
            throw std::runtime_error( errorMessage );
        }

        variableFunction = std::bind( &reference_frames::AerodynamicAngleCalculator::getCurrentGroundspeedBasedBodyFixedVelocity,
                                      bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ) );
        parameterSize = 3;
        break;
    }
    case lvlh_to_inertial_frame_rotation_dependent_variable:
    {
        std::function< Eigen::Vector6d( ) > vehicleStateFunction =
                std::bind( &simulation_setup::Body::getState, bodyMap.at( dependentVariableSettings->associatedBody_ ) );
        std::function< Eigen::Vector6d( ) > centralBodyStateFunction;

        if( ephemerides::isFrameInertial( dependentVariableSettings->secondaryBody_ ) )
        {
            centralBodyStateFunction =  [ ]( ){ return Eigen::Vector6d::Zero( ); };
        }
        else
        {
            centralBodyStateFunction =
                    std::bind( &simulation_setup::Body::getState, bodyMap.at( dependentVariableSettings->secondaryBody_ ) );
        }

        std::function< Eigen::Matrix3d( ) > rotationFunction =
                std::bind( &reference_frames::getVelocityBasedLvlhToInertialRotationFromFunctions,
                           vehicleStateFunction, centralBodyStateFunction, true );
        variableFunction = std::bind(
                    &getVectorRepresentationForRotationMatrixFunction, rotationFunction );

        parameterSize = 9;

        break;
    }
    case total_torque_dependent_variable:
    {

        // Retrieve model responsible for computing accelerations of requested bodies.
        std::shared_ptr< RotationalMotionStateDerivative< StateScalarType, TimeType > > rotationalDynamicsModel =
                getRotationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
        variableFunction = std::bind( &RotationalMotionStateDerivative< StateScalarType, TimeType >::getTotalTorqueForBody,
                                      rotationalDynamicsModel, bodyWithProperty );
        parameterSize = 3;


        break;
    }
    case single_torque_dependent_variable:
    {
        // Check input consistency.
        std::shared_ptr< SingleTorqueDependentVariableSaveSettings > torqueDependentVariableSettings =
                std::dynamic_pointer_cast< SingleTorqueDependentVariableSaveSettings >( dependentVariableSettings );
        if( torqueDependentVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_torque_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            // Retrieve list of suitable torque models (size should be one to avoid ambiguities)
            std::vector< std::shared_ptr< basic_astrodynamics::TorqueModel > >
                    listOfSuitableTorqueModels = getTorqueBetweenBodies(
                        torqueDependentVariableSettings->associatedBody_,
                        torqueDependentVariableSettings->secondaryBody_,
                        stateDerivativeModels, torqueDependentVariableSettings->torqueModelType_ );


            if( listOfSuitableTorqueModels.size( ) != 1 )
            {
                std::string errorMessage = "Error when getting torque between bodies " +
                        torqueDependentVariableSettings->associatedBody_ + " and " +
                        torqueDependentVariableSettings->secondaryBody_ + " of type " +
                        std::to_string(
                            torqueDependentVariableSettings->torqueModelType_ ) +
                        ", no such torque found";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                //std::function< Eigen::Vector3d( ) > vectorFunction =
                variableFunction = std::bind( &basic_astrodynamics::TorqueModel::getTorque,
                                              listOfSuitableTorqueModels.at( 0 ) );
                parameterSize = 3;
            }
        }
        break;
    }
    case keplerian_state_dependent_variable:
    {
        // Retrieve model responsible for computing accelerations of requested bodies.
        std::string orbitingBody = dependentVariableSettings->associatedBody_;
        std::string centralBody = dependentVariableSettings->secondaryBody_;

        if( bodyMap.count( orbitingBody ) == 0 )
        {
            throw std::runtime_error( "Error when requesting kepler elements as dependent variables, orbiting body: '" + orbitingBody + "' not found" );
        }
        else if( bodyMap.count( centralBody ) == 0 )
        {
            throw std::runtime_error( "Error when requesting kepler elements as dependent variables, central body: '" + centralBody + "' not found" );
        }

        std::function< double( ) > centralBodyGravitationalParameter;
        if( bodyMap.at( centralBody )->getGravityFieldModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when requesting kepler elements as dependent variables, central body: '" + centralBody + "' has no gravity field" );
        }
        else
        {
            centralBodyGravitationalParameter = std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                           bodyMap.at( centralBody )->getGravityFieldModel( ) );
        }

        std::function< double( ) > orbitingBodyGravitationalParameter;
        std::function< double( ) > effectiveGravitationalParameter;
        if( bodyMap.at( orbitingBody )->getGravityFieldModel( ) != nullptr )
        {
            orbitingBodyGravitationalParameter = std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                            bodyMap.at( orbitingBody )->getGravityFieldModel( ) );
            effectiveGravitationalParameter = std::bind( &utilities::sumFunctionReturn< double >,
                                                         orbitingBodyGravitationalParameter,
                                                         centralBodyGravitationalParameter );

        }
        else
        {
            effectiveGravitationalParameter = centralBodyGravitationalParameter;
        }

        // Retrieve functions for positions of two bodies.
        std::function< Eigen::Vector6d( const Eigen::Vector6d&, const Eigen::Vector6d& ) > functionToEvaluate =
                std::bind( &linear_algebra::computeVectorDifference< 6 >, std::placeholders::_1, std::placeholders::_2 );
        std::function< Eigen::Vector6d( ) > firstInput =
                std::bind( &simulation_setup::Body::getState, bodyMap.at( orbitingBody ) );
        std::function< Eigen::Vector6d( ) > secondInput =
                std::bind( &simulation_setup::Body::getState, bodyMap.at( centralBody ) );
        std::function< Eigen::Vector6d( ) > relativeStateFunction = std::bind(
                    &evaluateBivariateReferenceFunction< Eigen::Vector6d, Eigen::Vector6d >,
                    functionToEvaluate, firstInput, secondInput );

        variableFunction = std::bind( &orbital_element_conversions::convertCartesianToKeplerianElementsFromFunctions< double >,
                                      relativeStateFunction, effectiveGravitationalParameter );
        parameterSize = 6;


        break;
    }
    case modified_equinocial_state_dependent_variable:
    {
        // Retrieve model responsible for computing accelerations of requested bodies.
        std::string orbitingBody = dependentVariableSettings->associatedBody_;
        std::string centralBody = dependentVariableSettings->secondaryBody_;

        if( bodyMap.count( orbitingBody ) == 0 )
        {
            throw std::runtime_error( "Error when requesting modified equinoctial elements as dependent variables, orbiting body: '" + orbitingBody + "' not found" );
        }
        else if( bodyMap.count( centralBody ) == 0 )
        {
            throw std::runtime_error( "Error when requesting modified equinoctial elements as dependent variables, central body: '" + centralBody + "' not found" );
        }

        std::function< double( ) > centralBodyGravitationalParameter;
        if( bodyMap.at( centralBody )->getGravityFieldModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when requesting modified equinoctial elements as dependent variables, central body: '" + centralBody + "' has no gravity field" );
        }
        else
        {
            centralBodyGravitationalParameter = std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                           bodyMap.at( centralBody )->getGravityFieldModel( ) );
        }

        std::function< double( ) > orbitingBodyGravitationalParameter;
        std::function< double( ) > effectiveGravitationalParameter;
        if( bodyMap.at( orbitingBody )->getGravityFieldModel( ) != nullptr )
        {
            orbitingBodyGravitationalParameter = std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                            bodyMap.at( orbitingBody )->getGravityFieldModel( ) );
            effectiveGravitationalParameter = std::bind( &utilities::sumFunctionReturn< double >,
                                                         orbitingBodyGravitationalParameter,
                                                         centralBodyGravitationalParameter );

        }
        else
        {
            effectiveGravitationalParameter = centralBodyGravitationalParameter;
        }

        // Retrieve functions for positions of two bodies.
        std::function< Eigen::Vector6d( const Eigen::Vector6d&, const Eigen::Vector6d& ) > functionToEvaluate =
                std::bind( &linear_algebra::computeVectorDifference< 6 >, std::placeholders::_1, std::placeholders::_2 );
        std::function< Eigen::Vector6d( ) > firstInput =
                std::bind( &simulation_setup::Body::getState, bodyMap.at( orbitingBody ) );
        std::function< Eigen::Vector6d( ) > secondInput =
                std::bind( &simulation_setup::Body::getState, bodyMap.at( centralBody ) );
        std::function< Eigen::Vector6d( ) > relativeStateFunction = std::bind(
                    &evaluateBivariateReferenceFunction< Eigen::Vector6d, Eigen::Vector6d >,
                    functionToEvaluate, firstInput, secondInput );

        variableFunction = std::bind(
                    &orbital_element_conversions::convertCartesianToModifiedEquinoctialElementsFromStateFunction< double >,
                    relativeStateFunction, effectiveGravitationalParameter );
        parameterSize = 6;


        break;
    }
    case body_fixed_relative_cartesian_position:
    {
        std::function< Eigen::Vector3d( ) > positionFunctionOfRelativeBody =
                std::bind( &simulation_setup::Body::getPosition, bodyMap.at( bodyWithProperty ) );
        std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody =
                std::bind( &simulation_setup::Body::getPosition, bodyMap.at( secondaryBody ) );
        std::function< Eigen::Quaterniond( ) > orientationFunctionOfCentralBody =
                std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodyMap.at( secondaryBody ) );

        variableFunction = std::bind(
                    &reference_frames::getBodyFixedCartesianPosition, positionFunctionOfCentralBody,
                    positionFunctionOfRelativeBody, orientationFunctionOfCentralBody );
        parameterSize = 3;
        break;
    }
    case body_fixed_relative_spherical_position:
    {
        std::function< Eigen::Vector3d( ) > positionFunctionOfRelativeBody =
                std::bind( &simulation_setup::Body::getPosition, bodyMap.at( bodyWithProperty ) );
        std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody =
                std::bind( &simulation_setup::Body::getPosition, bodyMap.at( secondaryBody ) );
        std::function< Eigen::Quaterniond( ) > orientationFunctionOfCentralBody =
                std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodyMap.at( secondaryBody ) );

        variableFunction = std::bind(
                    &reference_frames::getBodyFixedSphericalPosition, positionFunctionOfCentralBody,
                    positionFunctionOfRelativeBody, orientationFunctionOfCentralBody );
        parameterSize = 3;
        break;
    }
    case euler_angles_to_body_fixed_313:
    {
        std::function< Eigen::Quaterniond( ) > orientationFunctionOfBody =
                std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodyMap.at( bodyWithProperty ) );

        std::function< Eigen::Vector3d( const Eigen::Quaterniond ) > eulerAngleFunction =
                std::bind( &basic_mathematics::get313EulerAnglesFromQuaternion, std::placeholders::_1 );

        variableFunction = std::bind(
                    &evaluateReferenceFunction< Eigen::Vector3d, Eigen::Quaterniond >,
                    eulerAngleFunction, orientationFunctionOfBody );
        parameterSize = 3;
        break;
    }
    case body_fixed_thrust_vector:
    {

        if( bodyMap.at( bodyWithProperty )->getFlightConditions( )   == nullptr )
        {
            std::string errorMessage = "Error no atmospheric flight conditions available when requesting body-fixed thrust vector output of" +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
        {
            std::string errorMessage = "Error, no vehicle systems available when requesting body-fixed thrust vector output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        int debugInfo = bodyMap.at( bodyWithProperty )->getBislipSystems()->getDebugInfo();
        //if( debugInfo == 1 ){ std::cout << "Computing body-fixed thrust vector as a dependent variable" << std::endl; }

        variableFunction = std::bind( &bislip::Variables::computeBodyFixedThrustVector,
                                      bodyMap, bodyWithProperty );

        //if( debugInfo == 1 ){ std::cout << "Body-fixed thrust vector as a dependent variable has been calculated" << std::endl; }

        parameterSize = 3;
        break;
    }

    case body_fixed_aero_load_vector:
    {

        if( bodyMap.at( bodyWithProperty )->getFlightConditions( )   == nullptr )
        {
            std::string errorMessage = "Error no atmospheric flight conditions available when requesting body-fixed thrust vector output of" +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
        {
            std::string errorMessage = "Error, no vehicle systems available when requesting body-fixed thrust vector output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        int debugInfo = bodyMap.at( bodyWithProperty )->getBislipSystems()->getDebugInfo();
        //if( debugInfo == 1 ){ std::cout << "Computing body-fixed aerodynamic load vector as a dependent variable" << std::endl; }


        variableFunction = std::bind( &bislip::Variables::computeBodyFixedAerodynamicLoad, bodyMap, bodyWithProperty );

        //if( debugInfo == 1 ){ std::cout << "Body-fixed aerodynamic load vector as a dependent variable has been calculated" << std::endl; }

        parameterSize = 3;
        break;
    }

    case local_gravity:
    {

        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
        {
            std::string errorMessage = "Error, no atmospheric flight conditions available when requesting local gravity output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        if( std::dynamic_pointer_cast< aerodynamics::AerodynamicCoefficientInterface >(
                    bodyMap.at( bodyWithProperty )->getAerodynamicCoefficientInterface( ) )== nullptr )
        {
            std::string errorMessage = "Error, no aerodynamic coefficient interface available when requesting local gravity output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        if( bodyMap.at( secondaryBody )->getGravityFieldModel( )->getGravitationalParameter()== 0.0 )
        {
            std::string errorMessage = "Error, no gravitational parameter available when requesting local gravity output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }
        if( bodyMap.at( secondaryBody )->getShapeModel( )->getAverageRadius()== 0.0 )
        {
            std::string errorMessage = "Error, no average radius available when requesting local gravity output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        variableFunction = std::bind( &bislip::Variables::computeLocalGravity, bodyMap, bodyWithProperty, secondaryBody );

        parameterSize = 2;
        break;
    }
    case body_fixed_total_load_vector:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( )   == nullptr )
        {
            std::string errorMessage = "Error no atmospheric flight conditions available when requesting body-fixed thrust vector output of" +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
        {
            std::string errorMessage = "Error, no vehicle systems available when requesting body-fixed thrust vector output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        variableFunction = std::bind( &bislip::Variables::computeBodyFixedTotalLoad,
                                      bodyMap, bodyWithProperty );
        parameterSize = 3;
        break;
    }
    case body_fixed_total_g_load_vector:
    {

        if( bodyMap.at( bodyWithProperty )->getFlightConditions( )   == nullptr )
        {
            std::string errorMessage = "Error no atmospheric flight conditions available when requesting body-fixed thrust vector output of" +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
        {
            std::string errorMessage = "Error, no vehicle systems available when requesting body-fixed thrust vector output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        variableFunction = std::bind( &bislip::Variables::computeBodyFixedTotal_g_Load_Vector,
                                      bodyMap, bodyWithProperty );
        parameterSize = 3;
        break;
    }
    case passenger_fixed_total_g_load_vector:
    {

        variableFunction = std::bind( &bislip::Variables::computePassengerFrameTotal_g_Load_Vector, bodyMap, bodyWithProperty );
        parameterSize = 3;
        break;
    }
    case aerodynamic_frame_aerodynamic_load_vector:
    {
        variableFunction = std::bind( &bislip::Variables::computeAerodynamicFrameAerodynamicLoad, bodyMap, bodyWithProperty );
        parameterSize = 3;
        break;
    }
    case aerodynamic_frame_total_load_vector:
    {
        variableFunction = std::bind( &bislip::Variables::computeAerodynamicFrameTotalLoad, bodyMap, bodyWithProperty );
        parameterSize = 3;
        break;
    }



#if( BUILD_WITH_ESTIMATION_TOOLS )
    case acceleration_partial_wrt_body_translational_state:
    {
        std::shared_ptr< AccelerationPartialWrtStateSaveSettings > accelerationPartialVariableSettings =
                std::dynamic_pointer_cast< AccelerationPartialWrtStateSaveSettings >( dependentVariableSettings );
        if( accelerationPartialVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type acceleration_partial_wrt_body_translational_state";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            if( stateDerivativePartials.count( translational_state ) == 0 )
            {
                throw std::runtime_error( "Error when requesting acceleration_partial_wrt_body_translational_state dependent variable, no translational state partials found." );
            }

            std::shared_ptr< acceleration_partials::AccelerationPartial > partialToUse =
                    getAccelerationPartialForBody(
                        stateDerivativePartials.at( translational_state ), accelerationPartialVariableSettings->accelerationModelType_,
                        accelerationPartialVariableSettings->associatedBody_,
                        accelerationPartialVariableSettings->secondaryBody_,
                        accelerationPartialVariableSettings->thirdBody_ );

            std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int > partialFunction =
                    partialToUse->getDerivativeFunctionWrtStateOfIntegratedBody(
                        std::make_pair( accelerationPartialVariableSettings->derivativeWrtBody_, "" ),
                        propagators::translational_state );

            if( partialFunction.second == 0 )
            {
                variableFunction = [ = ]( ){ return Eigen::VectorXd::Zero( 18 ); };
            }
            else
            {
                variableFunction = std::bind( &getVectorFunctionFromBlockFunction, partialFunction.first, 3, 6 );
            }

            parameterSize = 18;
        }
        break;
    }
#endif
    default:
        std::string errorMessage =
                "Error, did not recognize vector dependent variable type when making variable function: " +
                std::to_string( dependentVariableSettings->dependentVariableType_ );
        throw std::runtime_error( errorMessage );
    }
    return std::make_pair( variableFunction, parameterSize );
}

//! Acces element at index function
/*!
 * Acces element at index function
 * \param vectorFunction The function returning an Eigen vector
 * \param index The index to be accessed
 * \return The value of vector at index
 */
inline double elementAtIndexFunction( const std::function< Eigen::VectorXd( ) >& vectorFunction, const int index )
{
    return vectorFunction( )( index );
}

//! Function to create a function returning a requested dependent variable value (of type double).
/*!
 *  Function to create a function returning a requested dependent variable value (of type double), retrieved from
 *  environment and/or state derivative models.
 *  \param dependentVariableSettings Settings for dependent variable that is to be returned by function created here.
 *  \param bodyMap List of bodies to use in simulations (containing full environment).
 *  \param stateDerivativeModels List of state derivative models used in simulations (sorted by dynamics type as key).
 *  \param stateDerivativePartials List of state derivative partials used in simulations (sorted by dynamics type as key).
 *  \return Function returning requested dependent variable. NOTE: The environment and state derivative models need to
 *  be updated to current state and independent variable before computation is performed.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::function< double( ) > getDoubleDependentVariableFunction(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels =
        std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ),
        const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials =
        std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >( ) )
{
    const int componentIndex = dependentVariableSettings->componentIndex_;
    const int dependentVariableSize = getDependentVariableSize( dependentVariableSettings );
    if ( componentIndex >= 0 )
    {
        if ( dependentVariableSettings->componentIndex_ > dependentVariableSize - 1 )
        {
            throw std::runtime_error( "Error, cannot access component of variable because it exceeds its size" );
        }

        const std::pair< std::function< Eigen::VectorXd( ) >, int > vectorFunction =
                getVectorDependentVariableFunction< TimeType, StateScalarType >(
                    dependentVariableSettings, bodyMap, stateDerivativeModels, stateDerivativePartials );
        return std::bind( &elementAtIndexFunction, vectorFunction.first, componentIndex );
    }
    else
    {
        std::function< double( ) > variableFunction;

        // Retrieve base information on dependent variable
        PropagationDependentVariables dependentVariable = dependentVariableSettings->dependentVariableType_;
        const std::string& bodyWithProperty = dependentVariableSettings->associatedBody_;
        const std::string& secondaryBody = dependentVariableSettings->secondaryBody_;

        // Check dependent variable type and create function accordingly.
        switch( dependentVariable )
        {
        case mach_number_dependent_variable:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting Mach number output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::function< double( const double, const double ) > functionToEvaluate =
                    std::bind( &aerodynamics::computeMachNumber, std::placeholders::_1, std::placeholders::_2 );

            // Retrieve functions for airspeed and speed of sound.
            std::function< double( ) > firstInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeed,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
            std::function< double( ) > secondInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentSpeedOfSound,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );


            variableFunction = std::bind( &evaluateBivariateFunction< double, double, double >,
                                          functionToEvaluate, firstInput, secondInput );
            break;
        }
        case altitude_dependent_variable:
        {
            if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
            {
                std::string errorMessage = "Error, no flight conditions available when requesting altitude output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &aerodynamics::FlightConditions::getCurrentAltitude,
                                          bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
            break;
        }
        case airspeed_dependent_variable:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting airspeed output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeed,
                                          std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                              bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;
        }
        case local_density_dependent_variable:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting density output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentDensity,
                                          std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                              bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;
        }
        case radiation_pressure_dependent_variable:
        {
            if( bodyMap.at( bodyWithProperty )->getRadiationPressureInterfaces( ).count( secondaryBody ) == 0 )
            {
                std::string errorMessage = "Error, no radiation pressure interfaces when requesting radiation pressure output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &electro_magnetism::RadiationPressureInterface::getCurrentRadiationPressure,
                                          bodyMap.at( bodyWithProperty )->getRadiationPressureInterfaces( ).at( secondaryBody ) );
            break;
        }
        case relative_distance_dependent_variable:
        {
            // Retrieve functions for positions of two bodies.
            std::function< double( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                    std::bind( &linear_algebra::computeNormOfVectorDifference, std::placeholders::_1, std::placeholders::_2 );
            std::function< Eigen::Vector3d( ) > firstInput =
                    std::bind( &simulation_setup::Body::getPosition, bodyMap.at( bodyWithProperty ) );
            std::function< Eigen::Vector3d( ) > secondInput =
                    std::bind( &simulation_setup::Body::getPosition, bodyMap.at( secondaryBody ) );

            variableFunction = std::bind(
                        &evaluateBivariateReferenceFunction< double, Eigen::Vector3d >, functionToEvaluate, firstInput, secondInput );
            break;
        }
        case relative_speed_dependent_variable:
        {
            // Retrieve functions for velicoty of two bodies.
            std::function< double( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                    std::bind( &linear_algebra::computeNormOfVectorDifference, std::placeholders::_1, std::placeholders::_2 );
            std::function< Eigen::Vector3d( ) > firstInput =
                    std::bind( &simulation_setup::Body::getVelocity, bodyMap.at( bodyWithProperty ) );
            std::function< Eigen::Vector3d( ) > secondInput =
                    std::bind( &simulation_setup::Body::getVelocity, bodyMap.at( secondaryBody ) );

            variableFunction = std::bind(
                        &evaluateBivariateReferenceFunction< double, Eigen::Vector3d >, functionToEvaluate, firstInput, secondInput );

            break;
        }
        case single_acceleration_norm_dependent_variable:
        {
            // Check input consistency
            std::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationDependentVariableSettings =
                    std::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >(
                        dependentVariableSettings );
            if( accelerationDependentVariableSettings == nullptr )
            {
                std::string errorMessage = "Error, inconsistent inout when creating dependent variable function of type single_acceleration_norm_dependent_variable";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                // Retrieve list of suitable acceleration models (size should be one to avoid ambiguities)
                std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                        listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                            accelerationDependentVariableSettings->associatedBody_,
                            accelerationDependentVariableSettings->secondaryBody_,
                            stateDerivativeModels, accelerationDependentVariableSettings->accelerationModelType_ );

                // Check if third-body counterpart of acceleration is found
                if( listOfSuitableAccelerationModels.size( ) == 0 && basic_astrodynamics::isAccelerationDirectGravitational(
                            accelerationDependentVariableSettings->accelerationModelType_ ) )
                {
                    listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                                accelerationDependentVariableSettings->associatedBody_,
                                accelerationDependentVariableSettings->secondaryBody_,
                                stateDerivativeModels, basic_astrodynamics::getAssociatedThirdBodyAcceleration(
                                    accelerationDependentVariableSettings->accelerationModelType_  ) );
                }

                if( listOfSuitableAccelerationModels.size( ) != 1 )
                {
                    std::string errorMessage = "Error when getting acceleration between bodies " +
                            accelerationDependentVariableSettings->associatedBody_ + " and " +
                            accelerationDependentVariableSettings->secondaryBody_ + " of type " +
                            std::to_string(
                                accelerationDependentVariableSettings->accelerationModelType_ ) +
                            ", no such acceleration found NORM";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    std::function< Eigen::Vector3d( ) > vectorFunction =
                            std::bind( &basic_astrodynamics::AccelerationModel3d::getAcceleration,
                                       listOfSuitableAccelerationModels.at( 0 ) );
                    variableFunction = std::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );
                }
            }
            break;
        }
        case total_acceleration_norm_dependent_variable:
        {
            // Retrieve model responsible for computing accelerations of requested bodies.
            std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                    getTranslationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
            std::function< Eigen::Vector3d( ) > vectorFunction =
                    std::bind( &NBodyStateDerivative< StateScalarType, TimeType >::getTotalAccelerationForBody,
                               nBodyModel, bodyWithProperty );
            variableFunction = std::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );

            break;
        }
        case total_torque_norm_dependent_variable:
        {

            // Retrieve model responsible for computing accelerations of requested bodies.
            std::shared_ptr< RotationalMotionStateDerivative< StateScalarType, TimeType > > rotationalDynamicsModel =
                    getRotationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
            std::function< Eigen::Vector3d( ) > vectorFunction  =
                    std::bind( &RotationalMotionStateDerivative< StateScalarType, TimeType >::getTotalTorqueForBody,
                               rotationalDynamicsModel, bodyWithProperty );
            variableFunction = std::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );

            break;
        }
        case single_torque_norm_dependent_variable:
        {
            // Check input consistency.
            std::shared_ptr< SingleTorqueDependentVariableSaveSettings > torqueDependentVariableSettings =
                    std::dynamic_pointer_cast< SingleTorqueDependentVariableSaveSettings >( dependentVariableSettings );
            if( torqueDependentVariableSettings == nullptr )
            {
                std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_torque_norm_dependent_variable";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                // Retrieve list of suitable torque models (size should be one to avoid ambiguities)
                std::vector< std::shared_ptr< basic_astrodynamics::TorqueModel > >
                        listOfSuitableTorqueModels = getTorqueBetweenBodies(
                            torqueDependentVariableSettings->associatedBody_,
                            torqueDependentVariableSettings->secondaryBody_,
                            stateDerivativeModels, torqueDependentVariableSettings->torqueModelType_ );


                if( listOfSuitableTorqueModels.size( ) != 1 )
                {
                    std::string errorMessage = "Error when getting torque between bodies " +
                            torqueDependentVariableSettings->associatedBody_ + " and " +
                            torqueDependentVariableSettings->secondaryBody_ + " of type " +
                            std::to_string(
                                torqueDependentVariableSettings->torqueModelType_ ) +
                            ", no such torque found";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    //std::function< Eigen::Vector3d( ) > vectorFunction =
                    std::function< Eigen::Vector3d( ) > vectorFunction =
                            std::bind( &basic_astrodynamics::TorqueModel::getTorque,
                                       listOfSuitableTorqueModels.at( 0 ) );
                    variableFunction = std::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );
                }
            }
            break;
        }
        case relative_body_aerodynamic_orientation_angle_variable:
        {
            if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
            {
                std::string errorMessage = "Error when flight conditions for relative_body_aerodynamic_orientation_angle_variable output " +
                        bodyWithProperty + " has no flight conditions";
                throw std::runtime_error( errorMessage );
            }

            // Check input consistency.
            std::shared_ptr< BodyAerodynamicAngleVariableSaveSettings > bodyAerodynamicAngleVariableSaveSettings =
                    std::dynamic_pointer_cast< BodyAerodynamicAngleVariableSaveSettings >( dependentVariableSettings );
            if( bodyAerodynamicAngleVariableSaveSettings == nullptr )
            {
                std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type relative_body_aerodynamic_orientation_angle_variable";
                throw std::runtime_error( errorMessage );
            }

            std::function< double( ) > orientation_angle_rad = std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                                                                          bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                                                                          bodyAerodynamicAngleVariableSaveSettings->angle_ );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, orientation_angle_rad );

            break;
        }
        case total_aerodynamic_g_load_variable:
        {
            // Check input consistency
            std::shared_ptr< SingleAccelerationDependentVariableSaveSettings > aerodynamicAccelerationSettings =
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::aerodynamic, bodyWithProperty, secondaryBody, 0 );
            std::function< Eigen::VectorXd( ) > aerodynamicAccelerationFunction =
                    getVectorDependentVariableFunction( aerodynamicAccelerationSettings, bodyMap, stateDerivativeModels ).first;


            std::function< double( Eigen::Vector3d )> functionToEvaluate =
                    std::bind( &aerodynamics::computeAerodynamicLoadFromAcceleration, std::placeholders::_1 );
            variableFunction = std::bind(
                        &evaluateReferenceFunction< double, Eigen::Vector3d >,
                        functionToEvaluate, aerodynamicAccelerationFunction );

            break;
        }
        case stagnation_point_heat_flux_dependent_variable:
        {

            std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions =
                    std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
            if( flightConditions == nullptr )
            {
                std::string errorMessage = "Error no atmospheric flight conditions available when requesting stagnation point heating output of" +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
            {
                std::string errorMessage = "Error, no vehicle systems available when requesting stagnation point heating output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::shared_ptr< system_models::VehicleSystems > vehicleSystems =
                    bodyMap.at( bodyWithProperty )->getVehicleSystems( );


            if( !( vehicleSystems->getNoseRadius( ) == vehicleSystems->getNoseRadius( ) ) )
            {
                std::string errorMessage = "Error, no nose radius available when requesting stagnation point heating output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( !( vehicleSystems->getWallEmissivity( ) == vehicleSystems->getWallEmissivity( ) ) )
            {
                std::string errorMessage = "Error, no wall emmisivityavailable when requesting stagnation point heating output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            variableFunction = std::bind(
                        &computeEquilibriumFayRiddellHeatFluxFromProperties,
                        flightConditions, vehicleSystems );

            break;
        }
        case local_temperature_dependent_variable:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting temperature output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentFreestreamTemperature,
                                          std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                              bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;
        }
        case local_dynamic_pressure_dependent_variable:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) == nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting dynamic pressure "
                                           "output of " + bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentDynamicPressure,
                                          std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                              bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;
        }
        case local_aerodynamic_heat_rate_dependent_variable:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) == nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting heat rate "
                                           "output of " + bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAerodynamicHeatRate,
                                          std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                              bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;
        }
        case geodetic_latitude_dependent_variable:
        {
            if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
            {
                std::string errorMessage = "Error, no flight conditions available when requesting geodetic latitude output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            variableFunction = std::bind( &aerodynamics::FlightConditions::getCurrentGeodeticLatitude,
                                          bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
            break;
        }
        case control_surface_deflection_dependent_variable:
        {
            if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
            {
                std::string errorMessage = "Error, no vehicle systems available when requesting control surface deflection output of " +
                        bodyWithProperty + "with surface" + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::function< double( ) > deflectionAngle = std::bind( &system_models::VehicleSystems::getCurrentControlSurfaceDeflection,
                                                                    bodyMap.at( bodyWithProperty )->getVehicleSystems( ),
                                                                    dependentVariableSettings->secondaryBody_ );


            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, deflectionAngle );


            break;
        }
        case total_mass_rate_dependent_variables:
        {
            // Retrieve model responsible for computing mass rate of requested bodies.
            std::shared_ptr< BodyMassStateDerivative< StateScalarType, TimeType > > nBodyModel =
                    getBodyMassStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
            variableFunction =
                    std::bind( &BodyMassStateDerivative< StateScalarType, TimeType >::getTotalMassRateForBody, nBodyModel,
                               bodyWithProperty );

            break;
        }
        case periapsis_altitude_dependent_variable:
        {
            using namespace Eigen;
            std::function< double( const Vector6d&, const double, const double ) > functionToEvaluate =
                    std::bind( &basic_astrodynamics::computePeriapsisAltitudeFromCartesianState, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

            // Retrieve function for propagated body's Cartesian state in the global reference frame.
            std::function< Vector6d( ) > propagatedBodyStateFunction =
                    std::bind( &simulation_setup::Body::getState, bodyMap.at( bodyWithProperty ) );

            // Retrieve function for central body's Cartesian state in the global reference frame.
            std::function< Vector6d( ) > centralBodyStateFunction =
                    std::bind( &simulation_setup::Body::getState, bodyMap.at( secondaryBody ) );

            // Retrieve function for propagated body's Cartesian state in the propagation reference frame.
            std::function< Vector6d( ) > firstInput =
                    std::bind( &utilities::subtractFunctionReturn< Vector6d >,
                               propagatedBodyStateFunction, centralBodyStateFunction );

            // Retrieve function for central body's gravitational parameter.
            std::function< double( ) > secondInput =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodyMap.at( secondaryBody )->getGravityFieldModel( ) );

            // Retrieve function for central body's average radius.
            std::function< double( ) > thirdInput =
                    std::bind( &basic_astrodynamics::BodyShapeModel::getAverageRadius,
                               bodyMap.at( secondaryBody )->getShapeModel( ) );

            variableFunction = std::bind( &evaluateTrivariateFunction< double, Vector6d, double, double >,
                                          functionToEvaluate, firstInput, secondInput, thirdInput );
            break;
        }
        case specific_energy:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting specific energy output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::function< double( const double&, const double& ) > functionToEvaluate =
                    std::bind( &bislip::Variables::computeSpecificEnergy, std::placeholders::_1, std::placeholders::_2 );

            // Retrieve functions for airspeed and altitude.
            std::function< double( ) > firstInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAltitude,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );

            std::function< double( ) > secondInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeed,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );

            variableFunction = std::bind( &evaluateBivariateReferenceFunction< double, double, double >,
                                          functionToEvaluate, firstInput, secondInput );
            break;
        }
        case normalized_specific_energy:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting normalized specific energy output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::function< double( const double&, const double&, const double& ) > functionToEvaluate =
                    std::bind( &bislip::Variables::computeNormalizedSpecificEnergy, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

            // Retrieve functions for airspeed and altitude.
            std::function< double( ) > firstInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAltitude,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );

            std::function< double( ) > secondInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeed,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );

            // Retrieve functions for maximum Energy level.
            std::function< double( ) > thirdInput =
                    std::bind( &bislip::BislipVehicleSystems::getMaximumSpecificEnergy, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            variableFunction = std::bind( &evaluateTrivariateReferenceFunction< double, double, double, double >,
                                          functionToEvaluate, firstInput, secondInput, thirdInput );
            break;
        }
        case evaluated_throttle_setting:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting throttle setting output of " +
                        bodyWithProperty;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &bislip::Variables::evaluateGuidanceInterpolator, bislip::Parameters::Interpolators::ThrottleSetting, bodyMap, bodyWithProperty, secondaryBody );
            break;
        }
        case commanded_throttle_setting:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting throttle setting output of " +
                        bodyWithProperty;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &bislip::BislipVehicleSystems::getCurrentThrottleSetting, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            break;
        }
        case evaluated_thrust_elevation_angle:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting thrust elevation angle output of " +
                        bodyWithProperty;
                throw std::runtime_error( errorMessage );
            }
            std::function< double( ) > evaluated_thrust_elevation_angle_deg = std::bind( &bislip::Variables::evaluateGuidanceInterpolator, bislip::Parameters::Interpolators::ThrustElevationAngle, bodyMap, bodyWithProperty, secondaryBody );

            //std::function< double( const double& ) > rad2deg_Def =
            //        std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            //variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, evaluated_thrust_elevation_angle_rad );
            variableFunction = evaluated_thrust_elevation_angle_deg;

            break;
        }
        case commanded_thrust_elevation_angle:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting thrust elevation angle output of " +
                        bodyWithProperty;
                throw std::runtime_error( errorMessage );
            }
            std::function< double( ) > commanded_thrust_elevation_angle_rad = std::bind( &bislip::BislipVehicleSystems::getCurrentThrustElevationAngle, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, commanded_thrust_elevation_angle_rad );
            break;
        }
        case evaluated_thrust_azimuth_angle:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting thrust azimuth angle output of " +
                        bodyWithProperty;
                throw std::runtime_error( errorMessage );
            }
            std::function< double( ) > evaluated_thrust_azimuth_angle_deg = std::bind( &bislip::Variables::evaluateGuidanceInterpolator, bislip::Parameters::Interpolators::ThrustAzimuthAngle, bodyMap, bodyWithProperty, secondaryBody );

            //std::function< double( const double& ) > rad2deg_Def =
            //       std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            //variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, evaluated_thrust_azimuth_angle_rad );
            variableFunction = evaluated_thrust_azimuth_angle_deg;
            break;
        }
        case commanded_thrust_azimuth_angle:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting thrust azimuth angle output of " +
                        bodyWithProperty;
                throw std::runtime_error( errorMessage );
            }
            std::function< double( ) > commanded_thrust_azimuth_angle_rad = std::bind( &bislip::BislipVehicleSystems::getCurrentThrustAzimuthAngle, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, commanded_thrust_azimuth_angle_rad );
            break;
        }
        case evaluated_angle_of_attack:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting angle of attack output of " +
                        bodyWithProperty;
                throw std::runtime_error( errorMessage );
            }
            std::function< double( ) > evaluated_angle_of_attack_deg = std::bind( &bislip::Variables::evaluateGuidanceInterpolator, bislip::Parameters::Interpolators::AngleOfAttack, bodyMap, bodyWithProperty, secondaryBody );

            //std::function< double( const double& ) > rad2deg_Def =
            //        std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            //variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, evaluated_angle_of_attack_rad );
            variableFunction = evaluated_angle_of_attack_deg;
            break;
        }
        case commanded_angle_of_attack:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting angle of attack output of " +
                        bodyWithProperty;
                throw std::runtime_error( errorMessage );
            }
            std::function< double( ) > commanded_angle_of_attack_rad  = std::bind( &bislip::BislipVehicleSystems::getCurrentAngleOfAttack, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, commanded_angle_of_attack_rad );
            break;
        }
        case evaluated_bank_angle:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting bank elevation angle output of " +
                        bodyWithProperty;
                throw std::runtime_error( errorMessage );
            }
            std::function< double( ) > evaluated_bank_angle_deg = std::bind( &bislip::Variables::evaluateGuidanceInterpolator, bislip::Parameters::Interpolators::BankAngle, bodyMap, bodyWithProperty, secondaryBody );

            //std::function< double( const double& ) > rad2deg_Def =
            //      std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            //variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, evaluated_bank_angle_rad );
            variableFunction = evaluated_bank_angle_deg;

            break;
        }
        case commanded_bank_angle:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting bank angle output of " +
                        bodyWithProperty;
                throw std::runtime_error( errorMessage );
            }

            std::function< double( ) > commanded_bank_angle_rad = std::bind( &bislip::BislipVehicleSystems::getCurrentBankAngle, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, commanded_bank_angle_rad );
            break;
        }
        case current_mass:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting current mass output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &simulation_setup::Body::getBodyMass, bodyMap.at( bodyWithProperty ) );
            break;
        }
        case engine_status:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting current mass output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            //std::function< double( const double&, const double&, const double& ) > functionToEvaluate =
            //        std::bind( &bislip::Variables::determineEngineStatus, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

            std::function< double( const double&, const double& ) > functionToEvaluate =
                    std::bind( &bislip::Variables::determineEngineStatus, std::placeholders::_1, std::placeholders::_2 );

            // Retrieve functions for current mass and landing mass
            std::function< double( ) > firstInput =
                    std::bind( &simulation_setup::Body::getBodyMass, bodyMap.at( bodyWithProperty ) );
            std::function< double( ) > secondInput =
                    std::bind( &system_models::VehicleSystems::getDryMass, bodyMap.at( bodyWithProperty )->getVehicleSystems() );
            /*
            // Define function for normalized specific energy
            std::function< double( const double&, const double&, const double& ) > thirdInputDef =
                    std::bind( &bislip::Variables::computeNormalizedSpecificEnergy, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

            // Retrieve functions for airspeed and altitude.
            std::function< double( ) > firstInputForThirdInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAltitude,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );

            std::function< double( ) > secondInputForThirdInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeed,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );

            // Retrieve function for maximum Energy level.
            std::function< double( ) > thirdInputForThirdInput =
                    std::bind( &bislip::BislipVehicleSystems::getE_max, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            // Retrieve function for normalized specific energy
            std::function< double( ) > thirdInput =
                    std::bind( &evaluateTrivariateReferenceFunction< double, double, double, double >,
                               thirdInputDef, firstInputForThirdInput, secondInputForThirdInput, thirdInputForThirdInput );

            variableFunction = std::bind( &evaluateTrivariateReferenceFunction< double, double, double, double >,
                                          functionToEvaluate, firstInput, secondInput, thirdInput );
           */
            variableFunction = std::bind( &evaluateBivariateReferenceFunction< double, double, double >,
                                          functionToEvaluate, firstInput, secondInput );

            break;
        }
        case angular_distance_traveled:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting angular distance traveled output of " +
                        bodyWithProperty + "w.r.t. origin on " + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::function< double( const double&, const double&, const double&, const double& ) > functionToEvaluate =
                    std::bind( &bislip::Variables::computeAngularDistance, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4 );

            // Retrieve functions for initial coordinates.
            std::function< double( ) > firstInput =
                    std::bind( &bislip::BislipVehicleSystems::getInitialLat, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            std::function< double( ) > secondInput =
                    std::bind( &bislip::BislipVehicleSystems::getInitialLon, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            // Retrieve functions for current coordinates.
            std::function< double( ) > thirdInput =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::latitude_angle );
            std::function< double( ) > fourthInput =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::longitude_angle );

            std::function< double( ) > angular_distance_traveled_rad = std::bind( &evaluateTetravariateReferenceFunction< double, double, double, double, double >,
                                                                                  functionToEvaluate, firstInput, secondInput, thirdInput, fourthInput );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, angular_distance_traveled_rad );

            break;
        }
        case angular_distance_to_go:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting angular distance to go output of " +
                        bodyWithProperty + "w.r.t. final destination on " + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::function< double( const double&, const double&, const double&, const double& ) > functionToEvaluate =
                    std::bind( &bislip::Variables::computeAngularDistance, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4 );

            // Retrieve functions for current coordinates.
            std::function< double( ) > firstInput =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::latitude_angle );
            std::function< double( ) > secondInput =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::longitude_angle );

            // Retrieve functions for target coordinates.
            std::function< double( ) > thirdInput =
                    std::bind( &bislip::BislipVehicleSystems::getTargetLat, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            std::function< double( ) > fourthInput =
                    std::bind( &bislip::BislipVehicleSystems::getTargetLon, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            std::function< double( ) > angular_distance_to_go_rad = std::bind( &evaluateTetravariateReferenceFunction< double, double, double, double, double >,
                                                                               functionToEvaluate, firstInput, secondInput, thirdInput, fourthInput );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, angular_distance_to_go_rad );

            break;
        }
        case heading_to_target:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting heading to target output of " +
                        bodyWithProperty + "w.r.t. final destination on " + secondaryBody;
                throw std::runtime_error( errorMessage );
            }


            // Retrieve functions for current coordinates.
            std::function< double( ) > firstInput =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::latitude_angle );
            std::function< double( ) > secondInput =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::longitude_angle );

            // Retrieve functions for target coordinates.
            std::function< double( ) > thirdInput =
                    std::bind( &bislip::BislipVehicleSystems::getTargetLat, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            std::function< double( ) > fourthInput =
                    std::bind( &bislip::BislipVehicleSystems::getTargetLon, bodyMap.at( bodyWithProperty )->getBislipSystems() );


            std::function< double( const double&, const double&, const double&, const double& ) > headingToTarget_radDef =
                    std::bind( &bislip::Variables::computeHeadingToTarget, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4 );

            std::function< double( ) > headingToTarget_rad = std::bind( &evaluateTetravariateReferenceFunction< double, double, double, double, double >,
                                                                        headingToTarget_radDef, firstInput, secondInput, thirdInput, fourthInput );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, headingToTarget_rad );


            break;
        }
        case heading_error:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting heading error output of " +
                        bodyWithProperty + "w.r.t. final destination on " + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            // Retrieve functions for current coordinates.
            std::function< double( ) > firstInput =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::latitude_angle );
            std::function< double( ) > secondInput =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::longitude_angle );

            // Retrieve functions for target coordinates.
            std::function< double( ) > thirdInput =
                    std::bind( &bislip::BislipVehicleSystems::getTargetLat, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            std::function< double( ) > fourthInput =
                    std::bind( &bislip::BislipVehicleSystems::getTargetLon, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            // Retrieve function for heading angle
            std::function< double( ) > headingAngle_rad =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::heading_angle );

            std::function< double( const double& ) > neg2pos_Def =
                    std::bind( &bislip::Variables::convertNegativeAnglesToPositive, std::placeholders::_1 );

            std::function< double( ) > headingAngle_radPos =
                    std::bind( &evaluateReferenceFunction< double, double >, neg2pos_Def, headingAngle_rad );

            std::function< double( const double&, const double&, const double&, const double&, const double& ) > headingError_radDef =
                    std::bind( &bislip::Variables::computeHeadingError, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );

            std::function< double( ) > headingError_rad = std::bind( &evaluatePentavariateReferenceFunction< double, double, double, double, double, double >,
                                                                     headingError_radDef, firstInput, secondInput, thirdInput, fourthInput, headingAngle_radPos );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, headingError_rad );


            break;
        }
        case heat_flux_tauber:
        {
            if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
            {
                std::string errorMessage = "Error no atmospheric flight conditions available when requesting leading edge heating output of" +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
            {
                std::string errorMessage = "Error, no vehicle systems available when requesting leading edge heating output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions =
                    std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap.at( bodyWithProperty )->getFlightConditions( ) );

            std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( bodyWithProperty )->getBislipSystems( ) ;

            std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( bodyWithProperty )->getVehicleSystems( ) ;


            if( !( vehicleSystems->getNoseRadius( ) == vehicleSystems->getNoseRadius( ) ) )
            {
                std::string errorMessage = "Error, no nose radius available when requesting leading edge heating output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( !( vehicleSystems->getWallEmissivity( ) == vehicleSystems->getWallEmissivity( ) ) )
            {
                std::string errorMessage = "Error, no wall emmisivity available when requesting leading edge heating output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            // std::function< double( const double, const double, const double ) > functionToEvaluate =
            //       std::bind( &bislip::Variables::computeHeatingRateTauber, std::placeholders::_1, std::placeholders::_2 , std::placeholders::_3 );

            // Retrieve functions for Stagnation and Flat Plate heat fluxes.
            std::function< double( ) > firstInput =
                    std::bind( &bislip::Variables::computeStagnationHeatFlux, bodyMap, bodyWithProperty );
            std::function< double( ) > secondInput =
                    std::bind( &bislip::Variables::computeFlatPlateHeatFlux, bodyMap, bodyWithProperty );

            // Retrieve functions for wing sweep angle
            std::function< double( ) > thirdInput =
                    std::bind( &bislip::BislipVehicleSystems::getWingSweepAngle, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            // variableFunction = std::bind( &evaluateTrivariateFunction< double, double, double, double >,
            //                             functionToEvaluate, firstInput, secondInput, thirdInput );
            //variableFunction = std::bind( &bislip::Variables::computeHeatingRateTauber, bodyMap, bodyWithProperty );
            variableFunction = std::bind( &bislip::BislipVehicleSystems::getCurrentHeatFluxTauber, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            break;
        }
        case heat_flux_chapman:
        {
            //  std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( bodyWithProperty )->getBislipSystems( ) ;


            //bislipSystems->setWorkingRadius( bislipSystems->getNoseRadius( ) );


            //variableFunction = std::bind( &bislip::Variables::computeStagnationHeatFlux, bodyMap, bodyWithProperty );
            variableFunction = std::bind( &bislip::BislipVehicleSystems::getCurrentHeatFluxChapman, bodyMap.at( bodyWithProperty )->getBislipSystems() );



            break;
        }
        case bending_moment:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting bending moment output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            //std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions =
            //        std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >( bodyMap.at( bodyWithProperty )->getFlightConditions( ) );

            variableFunction = std::bind( &bislip::Variables::computeBendingMoment, bodyMap, bodyWithProperty );

            break;
        }
        case skip_suppression_limit:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting equilibrium glide limit output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( std::dynamic_pointer_cast< aerodynamics::AerodynamicCoefficientInterface >(
                        bodyMap.at( bodyWithProperty )->getAerodynamicCoefficientInterface( ) )== nullptr )
            {
                std::string errorMessage = "Error, no aerodynamic coefficient interface available when requesting equilibrium glide limit output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( bodyMap.at( bodyWithProperty )->getBislipSystems( )== nullptr )
            {
                std::string errorMessage = "Error, no Bislip systems available when requesting equilibrium glide limit output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::function< double( ) > bankAngleLimit_rad = std::bind( &bislip::Variables::computeSkipSuppressionLimit, bodyMap, bodyWithProperty, secondaryBody );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, bankAngleLimit_rad );
            break;
        }
        case bodyflap_deflection_moment_coefficient_increment:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting bodyflap deflection moment coefficient increment output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( std::dynamic_pointer_cast< aerodynamics::AerodynamicCoefficientInterface >(
                        bodyMap.at( bodyWithProperty )->getAerodynamicCoefficientInterface( ) )== nullptr )
            {
                std::string errorMessage = "Error, no aerodynamic coefficient interface available when requesting bodyflap deflection moment coefficient increment output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( bodyMap.at( bodyWithProperty )->getBislipSystems( )== nullptr )
            {
                std::string errorMessage = "Error, no Bislip systems available when requesting bodyflap deflection moment coefficient increment output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
            {
                std::string errorMessage = "Error, no vehicle systems available when requesting bodyflap deflection moment coefficient increment output of " +
                        bodyWithProperty + "with surface" + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            /*
            std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions =
                    std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap.at( bodyWithProperty )->getFlightConditions( ) );

            std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( bodyWithProperty )->getVehicleSystems( ) ;

            std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface =
                    std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >( bodyMap.at( bodyWithProperty )->getAerodynamicCoefficientInterface( ) );


            std::function< double( const double&, const double& ) > functionToEvaluate =
                    std::bind( &bislip::Variables::computeBodyFlapCmIncrementdif, std::placeholders::_1, std::placeholders::_2 );

            double bodyFlapDeflection = bodyMap.at( bodyWithProperty )->getVehicleSystems( )->getCurrentControlSurfaceDeflection( "BodyFlap" );
            double angleOfAttack = tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::AngleOfAttack, bodyMap, bodyWithProperty ) );

            //! Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
            std::vector< double > aerodynamicCoefficientInput;
            aerodynamicCoefficientInput.push_back( angleOfAttack );
            aerodynamicCoefficientInput.push_back( flightConditions->getCurrentMachNumber() );

            // Retrieve functions for Stagnation and Flat Plate heat fluxes.
            std::function< double( ) > firstInput =
                    std::bind( &bislip::Variables::computeFullPitchMomentCoefficient, coefficientInterface, bodyFlapDeflection, aerodynamicCoefficientInput );
            std::function< double( ) > secondInput =
                    std::bind( &bislip::Variables::computePitchMomentCoefficient, coefficientInterface, aerodynamicCoefficientInput );

            //variableFunction = std::bind( &evaluateBivariateReferenceFunction< double, double, double >,
                //                          functionToEvaluate, firstInput, secondInput );
            */
            //variableFunction = std::bind( &bislip::Variables::computeBendingMoment, bodyMap, bodyWithProperty );

            variableFunction = std::bind( &bislip::Variables::computeBodyFlapCmIncrement, bodyMap, bodyWithProperty );

            break;
        }
        case bodyflap_deflection_moment_coefficient_increment_dif:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                std::string errorMessage = "Error, no atmospheric flight conditions available when requesting bodyflap deflection moment coefficient increment output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( std::dynamic_pointer_cast< aerodynamics::AerodynamicCoefficientInterface >(
                        bodyMap.at( bodyWithProperty )->getAerodynamicCoefficientInterface( ) )== nullptr )
            {
                std::string errorMessage = "Error, no aerodynamic coefficient interface available when requesting bodyflap deflection moment coefficient increment output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( bodyMap.at( bodyWithProperty )->getBislipSystems( )== nullptr )
            {
                std::string errorMessage = "Error, no Bislip systems available when requesting bodyflap deflection moment coefficient increment output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
            {
                std::string errorMessage = "Error, no vehicle systems available when requesting bodyflap deflection moment coefficient increment output of " +
                        bodyWithProperty + "with surface" + secondaryBody;
                throw std::runtime_error( errorMessage );
            }


            std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions =
                    std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap.at( bodyWithProperty )->getFlightConditions( ) );

            std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( bodyWithProperty )->getBislipSystems( ) ;

            std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( bodyWithProperty )->getVehicleSystems( ) ;

            //std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface =
            //      std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >( bodyMap.at( bodyWithProperty )->getAerodynamicCoefficientInterface( ) );

            std::function< std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface >( ) > coefficientInterface =
                    std::bind( &tudat::simulation_setup::Body::getAerodynamicCoefficientInterface, bodyMap.at( bodyWithProperty ) );


            std::function< double( const double, const double ) > machNumberDef =
                    std::bind( &aerodynamics::computeMachNumber, std::placeholders::_1, std::placeholders::_2 );

            // Retrieve functions for airspeed and speed of sound.
            std::function< double( ) > airSpeed =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeed,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
            std::function< double( ) > speedOfSound =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentSpeedOfSound,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );

            std::function< double( ) > machNumber =
                    std::bind( &evaluateBivariateFunction< double, double, double >,
                               machNumberDef, airSpeed, speedOfSound );

            std::function< double( ) > angleOfAttack =
                    std::bind( &bislip::BislipVehicleSystems::getCurrentAngleOfAttack, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            std::function< double( ) > bodyFlapDeflection =
                    std::bind( &system_models::VehicleSystems::getCurrentControlSurfaceDeflection,
                               bodyMap.at( bodyWithProperty )->getVehicleSystems( ),
                               "BodyFlap" );

            std::function< double( ) > elevonDeflection =
                    std::bind( &system_models::VehicleSystems::getCurrentControlSurfaceDeflection,
                               bodyMap.at( bodyWithProperty )->getVehicleSystems( ),
                               "ElevonLeft" );

            std::function< std::vector< double >( const double&, const double& ) > aerodynamicCoefficientInputDef =
                    std::bind( &bislip::Variables::getAerodynamicCoefficientInput, std::placeholders::_1, std::placeholders::_2 );

            //  std::function< std::vector< double >( const double&, const double&, const double& ) > controlSurfaceCoefficientInputDef =
            //        std::bind( &bislip::Variables::getControlSurfaceCoefficientInput, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

            std::function< std::vector< double >( ) > aerodynamicCoefficientInput =
                    std::bind( &evaluateBivariateReferenceFunction< std::vector< double >, double, double >,
                               aerodynamicCoefficientInputDef, angleOfAttack, machNumber );

            // std::function< std::vector< double >( ) > controlSurfaceCoefficientInput =
            //        std::bind( &evaluateTrivariateReferenceFunction< std::vector< double >, double, double, double >,
            //                 controlSurfaceCoefficientInputDef, angleOfAttack, machNumber, bodyFlapDeflection );

            std::function< double( const double&, const double& ) > functionToEvaluate =
                    std::bind( &bislip::Variables::computeBodyFlapCmIncrementdif, std::placeholders::_1, std::placeholders::_2 );

            std::function< double( const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface >&, const std::vector< double >&, const double&, const double& ) > firstInputDef =
                    std::bind( &bislip::Variables::computeFullPitchMomentCoefficient, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 , std::placeholders::_4 );

            std::function< double( const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface >&, const std::vector< double >& ) > secondInputDef =
                    std::bind( &bislip::Variables::computePitchMomentCoefficient, std::placeholders::_1, std::placeholders::_2 );

            // Retrieve functions for Stagnation and Flat Plate heat fluxes.
            std::function< double( ) > firstInput =
                    std::bind( &evaluateTetravariateReferenceFunction< double, std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface >, std::vector< double >, double, double >,
                               firstInputDef, coefficientInterface, aerodynamicCoefficientInput, bodyFlapDeflection, elevonDeflection );

            //std::bind( &bislip::Variables::computeFullPitchMomentCoefficient, coefficientInterface, bodyFlapDeflection, aerodynamicCoefficientInput );
            std::function< double( ) > secondInput =
                    std::bind( &evaluateBivariateReferenceFunction< double, std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface >, std::vector< double > >,
                               secondInputDef, coefficientInterface, aerodynamicCoefficientInput );

            variableFunction = std::bind( &evaluateBivariateReferenceFunction< double, double, double >,
                                          functionToEvaluate, firstInput, secondInput );

            break;
        }
        case body_fixed_total_g_load_magnitude:
        {
            variableFunction = std::bind( &bislip::Variables::computeBodyFixedTotal_g_Load_Magnitude, bodyMap, bodyWithProperty );
            break;
        }
        case current_heading_error_deadband:
        {
            variableFunction = std::bind( &bislip::Variables::computeHeadingErrorDeadBand, bodyMap, bodyWithProperty );
            break;
        }
        case temp_bank_angle:
        {
            variableFunction = std::bind( &bislip::BislipVehicleSystems::getTempBankAngle, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            break;
        }
        case reversal_conditional:
        {
            /*
            std::function< double( ) > newBankAngle = std::bind( &bislip::BislipVehicleSystems::getTempBankAngle, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            //std::bind( &bislip::Variables::determineNewBankAngle, bodyMap, bodyWithProperty );

             //bislipSystems->setTempBankAngle( newBankAngle )
            //variableFunction = std::bind( &bislip::BislipVehicleSystems::getTempBankAngle, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            std::function< double( ) > currentLatitude =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::latitude_angle );
            std::function< double( ) > currentLongitude =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::longitude_angle );
            std::function< double( ) > currentHeading =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::heading_angle );

            std::function< double( ) > targatLat =
                    std::bind( &bislip::BislipVehicleSystems::getTargetLat, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            std::function< double( ) > targatLon =
                    std::bind( &bislip::BislipVehicleSystems::getTargetLon, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            std::function< double( const double&, const double&, const double&, const double&, const double& ) > currentHeadingError_rad_Def =
                    std::bind( &bislip::Variables::computeHeadingError, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_4 );

            std::function< double( ) > currentHeadingError_rad =
                    std::bind( &evaluatePentavariateReferenceFunction< double, double, double, double, double, double >,
                               currentHeadingError_rad_Def, currentLatitude, currentLongitude, targatLat, targatLon, currentHeading );

            std::function< double( ) > currentHeadingError_deg =
                    std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, currentHeadingError_rad );

            std::function< double( const double&, const double& ) > functionToEvaluate =
                    std::bind( &bislip::Variables::determineReversalConditional, std::placeholders::_1, std::placeholders::_2 );

            variableFunction = std::bind( &evaluateBivariateReferenceFunction< double, double, double >,
                                          functionToEvaluate, newBankAngle, currentHeadingError_deg );
*/
            variableFunction = std::bind( &bislip::BislipVehicleSystems::getReversalConditional, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            break;
        }
        case bank_reversal_trigger:
        {
            /*
            std::function< double( ) > currentLatitude =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::latitude_angle );
            std::function< double( ) > currentLongitude =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::longitude_angle );
            std::function< double( ) > currentHeading =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::heading_angle );
            std::function< double( ) > currentBankAngle =
                    std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                               bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                               tudat::reference_frames::bank_angle );
            std::function< double( ) > targatLat =
                    std::bind( &bislip::BislipVehicleSystems::getTargetLat, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            std::function< double( ) > targatLon =
                    std::bind( &bislip::BislipVehicleSystems::getTargetLon, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            std::function< double( const double& ) > absoluteValue_Def =
                    std::bind( &bislip::Variables::determineAbsoluteValue, std::placeholders::_1 );

            std::function< double( const double&, const double&, const double&, const double& ) > angularDistangeToGo_rad_Def =
                    std::bind( &bislip::Variables::computeAngularDistance, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4 );

            std::function< double( ) > angularDistangeToGo_rad =
                    std::bind( &evaluateTetravariateReferenceFunction< double, double, double, double, double >,
                               angularDistangeToGo_rad_Def, currentLatitude, currentLongitude, targatLat, targatLon );

            std::function< double( ) > angularDistangeToGo_deg =
                    std::bind( &evaluateReferenceFunction< double, double >,
                               rad2deg_Def, angularDistangeToGo_rad );

            std::function< double( const double&, const double&, const double&, const double&, const double& ) > currentHeadingError_rad_Def =
                    std::bind( &bislip::Variables::computeHeadingError, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_4 );

            std::function< double( ) > currentHeadingError_rad =
                    std::bind( &evaluatePentavariateReferenceFunction< double, double, double, double, double, double >,
                               currentHeadingError_rad_Def, currentLatitude, currentLongitude, targatLat, targatLon, currentHeading );

            std::function< double( ) > currentHeadingError_deg =
                    std::bind( &evaluateReferenceFunction< double, double >,
                               rad2deg_Def, currentHeadingError_rad );

            std::function< double( ) > absoluteValueCurrentHeadingError_deg =
                    std::bind( &evaluateReferenceFunction< double, double >,
                               absoluteValue_Def, currentHeadingError_deg );

            std::function< double( const double& ) > bankAngleSign_Def =
                    std::bind( &bislip::Variables::determineBankAngleSign, std::placeholders::_1 );

            std::function< double( ) > bankAngleSign =
                    std::bind( &evaluateReferenceFunction< double, double >,
                               bankAngleSign_Def, currentBankAngle );

            std::function< double( const double&, const double& ) > reversalConditional_Def =
                    std::bind( &bislip::Variables::determineReversalConditional, std::placeholders::_1, std::placeholders::_2 );

            std::function< double( ) > reversalConditional =
                    std::bind( &evaluateBivariateReferenceFunction< double, double, double >,
                               reversalConditional_Def, bankAngleSign, currentHeadingError_deg );



            //  This section is to extract the new bank angle function
            //! Extract Bislip Systems pointer.
            std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( bodyWithProperty )->getBislipSystems( ) ;

            //! Extract current bank angle.
            std::function< double( ) > currentBankAngle_rad = std::bind( &bislip::BislipVehicleSystems::getCurrentBankAngle, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            std::function< double( ) > currentBankAngle_deg = std::bind( &evaluateReferenceFunction< double, double >,
                                                                         rad2deg_Def, currentBankAngle_rad );

            std::function< int( const double& ) > bankAngleSign_Def = std::bind( &bislip::Variables::determineBankAngleSign, std::placeholders::_1 );

            std::function< int( ) > bankAngleSign = std::bind( &evaluateReferenceFunction< double, double >,
                                                               bankAngleSign_Def, currentBankAngle_deg );

            std::function< double( const bislip::Parameters::Optimization&, const tudat::simulation_setup::NamedBodyMap&, const std::string& ) > newBankAngle_Def =
                    std::bind( &bislip::Variables::evaluateGuidanceInterpolator, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

             std::function< bislip::Parameters::Optimization( ) > parameter = std::bind( &bislip::Variables::passOptimizationParameter, bislip::Parameters::Optimization::BankAngle );
*/

            //  std::function< double( const double&, const double& ) > newBankAngleSigned_Def =
            //        std::bind( &bislip::Variables::determineSignedBankAngle, std::placeholders::_1, std::placeholders::_2 );

            //std::function< double( ) > newBankAngleSigned = std::bind( &evaluateBivariateReferenceFunction< double, double, double >,
            //                                                         newBankAngleSigned_Def, bankAngleSign, newBankAngle );
            //*******************************************************


            /* std::function< double( const tudat::simulation_setup::NamedBodyMap&, const std::string&, const double& ) > functionToEvaluate =
                    std::bind( &bislip::Variables::determineBankAngleReversal, std::placeholders::_1, std::placeholders::_2 , std::placeholders::_3 );

            std::function< tudat::simulation_setup::NamedBodyMap&( ) > bodyMap_function = std::bind( &bislip::Variables::passBodyMap, bodyMap );

            //! Retrieve function for parameter string
            std::function< std::string( ) > bodyWithProperty_function = std::bind( &bislip::Variables::passString, bodyWithProperty );

            std::function< double( ) > newBankAngle = std::bind( &bislip::Variables::determineNewBankAngle, bodyMap, bodyWithProperty );

            variableFunction =  std::bind( &evaluateTrivariateReferenceFunction< double, tudat::simulation_setup::NamedBodyMap&, std::string, double >,
                                           functionToEvaluate, bodyMap_function, bodyWithProperty_function, newBankAngle );


*/


            variableFunction = std::bind( &bislip::BislipVehicleSystems::getBankAngleReversalTrigger, bodyMap.at( bodyWithProperty )->getBislipSystems() );


            break;
        }
        case current_lift_magnitude:
        {
            variableFunction = std::bind( &bislip::Variables::computeCurrentLiftForce, bodyMap, bodyWithProperty );
            break;
        }
        case wall_temperature_chapman:
        {
            /*

            //! Extract flight conditions pointer.
            std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions =
                    std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap.at( bodyWithProperty )->getFlightConditions( ) );

            //! Extract Bislip Systems pointer.
            std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( bodyWithProperty )->getBislipSystems( ) ;

            //! Extract Vehicle Systems pointer.
            std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( bodyWithProperty )->getVehicleSystems( ) ;


            int debugInfo = bislipSystems->getDebugInfo();

            if( debugInfo == 1 ){ std::cout << "Extract wall emissivity" << std::endl; }
            const double wallEmissivity = vehicleSystems->getWallEmissivity( );

            if( debugInfo == 1 ){ std::cout << "     currentAirTemperature = " << currentAirTemperature << std::endl; }
            if( debugInfo == 1 ){ std::cout << "     currentMachNumber     = " << currentMachNumber << std::endl; }

            if( debugInfo == 1 ){ std::cout << "     adiabaticWallTemperature = " << adiabaticWallTemperature << std::endl; }




            if( debugInfo == 1 ){ std::cout << "Find equilibrium wall temperature" << std::endl; }
            double equilibriumWallTemperature = bislip::Variables::findEquilibriumWallTemperature( heatTransferFunction, wallEmissivity, adiabaticWallTemperature);






            std::function< double( const double, const double ) > adiabaticWallTemperature_def =
                    std::bind( &bislip::Variables::computeAdiabaticWallTemperature, std::placeholders::_1, std::placeholders::_2 );

            // Retrieve functions for
            std::function< double( ) > firstInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentFreestreamTemperature,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
            std::function< double( ) > secondInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentMachNumber,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );

           std::function< double( ) >  adiabaticWallTemperature = std::bind( &evaluateBivariateReferenceFunction< double, double, double >,
                                          adiabaticWallTemperature_def, firstInput, secondInput );



           if( debugInfo == 1 ){ std::cout << "Get heat transfer function" << std::endl; }
         //  std::function< double( const tudat::simulation_setup::NamedBodyMap& , const std::string ) > heatTransferFunction_def =
         //          std::bind( &bislip::Variables::getStagnationHeatTransferFunction, std::placeholders::_1, std::placeholders::_2 );

//                   std::bind( &bislip::Variables::getStagnationHeatTransferFunction, bodyMap.at( bodyWithProperty )->getBislipSystems() );


           std::function< double( const tudat::simulation_setup::NamedBodyMap&, const std::string  ) > heatTransferFunction_def =
                   std::bind( &bislip::Variables::getStagnationHeatTransferFunction, std::placeholders::_1, std::placeholders::_2 );

           std::function< double(const tudat::simulation_setup::NamedBodyMap&, const std::string& ) > newBankAngle_Def =
                   std::bind( &bislip::Variables::getStagnationHeatTransferFunction, std::placeholders::_1, std::placeholders::_2 );



           std::function< double( const double ) > heatTransferFunction = std::bind( &bislip::Variables::getStagnationHeatTransferFunction, bodyMap, bodyWithProperty );



           std::function< double( const tudat::simulation_setup::NamedBodyMap&, const std::string& ) > functionToEvaluate =
                   std::bind( &bislip::Variables::getStagnationHeatTransferFunction, std::placeholders::_1, std::placeholders::_2 );

           std::function< tudat::simulation_setup::NamedBodyMap&( ) > bodyMap_function = std::bind( &bislip::Variables::passBodyMap, bodyMap );

           //! Retrieve function for parameter string
           std::function< std::string( ) > bodyWithProperty_function = std::bind( &bislip::Variables::passString, bodyWithProperty );


           variableFunction =  std::bind( &evaluateBivariateReferenceFunction< double, tudat::simulation_setup::NamedBodyMap&, std::string& >,
                                          functionToEvaluate, bodyMap_function, bodyWithProperty_function );









           std::function< double( const double& ) > rad2deg_Def =
                   std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

           variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, angular_distance_to_go_rad );





*/
            variableFunction = std::bind( &bislip::BislipVehicleSystems::getChapmanWallTemp, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            break;
        }
        case wall_temperature_tauber_stagnation:
        {
            variableFunction = std::bind( &bislip::BislipVehicleSystems::getTauberWallTempStagnation, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            break;
        }
        case wall_temperature_tauber_flatplate:
        {
            variableFunction = std::bind( &bislip::BislipVehicleSystems::getTauberWallTempFlatPlate, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            break;
        }
        case heat_flux_tauber_stagnation:
        {
            variableFunction = std::bind( &bislip::BislipVehicleSystems::getTauberHeatFluxStagnation, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            break;
        }
        case heat_flux_tauber_flatplate:
        {
            variableFunction = std::bind( &bislip::BislipVehicleSystems::getTauberHeatFluxFlatPlate, bodyMap.at( bodyWithProperty )->getBislipSystems() );
            break;
        }
        case cumulative_angular_distance_travelled:
        {
            variableFunction = std::bind( &bislip::Variables::computeCumulativeAngularDistanceTravelled, bodyMap, bodyWithProperty );
            break;
        }
        case groundtrack_difference:
        {
            variableFunction = std::bind( &bislip::Variables::computeCumulativeAngularDistanceTravelledDifference, bodyMap, bodyWithProperty );
            break;
        }
        case time_of_flight:
        {
            variableFunction = std::bind( &bislip::Variables::computeTimeOfFlight, bodyMap, bodyWithProperty );
            break;
        }
        case flight_path_angle_rate:
        {
            std::function< double( ) > flightPathAngle_rad  = std::bind( &bislip::Variables::computeFlightPathAngleRate, bodyMap, bodyWithProperty, secondaryBody );

            std::function< double( const double& ) > rad2deg_Def =
                  std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, flightPathAngle_rad );
            //variableFunction = std::bind( &bislip::BislipVehicleSystems::getCurrentFlightPathAngleRate, bodyMap.at( bodyWithProperty )->getBislipSystems() );

            break;
        }
        case cumulative_cartesian_distance_travelled:
        {
            variableFunction = std::bind( &bislip::Variables::computeCumulativeCartesianDistanceTravelled, bodyMap, bodyWithProperty );
            break;
        }
        case thrust_force_magnitude:
        {
            variableFunction = std::bind( &bislip::Variables::computeThrustMagnitude, bodyMap, bodyWithProperty );
            break;
        }
        case speed_of_sound:
        {
            variableFunction =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentSpeedOfSound,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;

        }
        case adiabatic_wall_temperature:
        {
            std::function< double( const double, const double ) > functionToEvaluate =
                    std::bind( &bislip::Variables::computeAdiabaticWallTemperature, std::placeholders::_1, std::placeholders::_2 );

            // Retrieve functions for
            std::function< double( ) > firstInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentFreestreamTemperature,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );
            std::function< double( ) > secondInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentMachNumber,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );

            variableFunction = std::bind( &evaluateBivariateReferenceFunction< double, double, double >,
                                          functionToEvaluate, firstInput, secondInput );


            break;

        }
        case current_drag_magnitude:
        {
            variableFunction = std::bind( &bislip::Variables::computeCurrentDragForce, bodyMap, bodyWithProperty );
            break;
        }

        case freestream_temperature:
        {

            variableFunction =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentFreestreamTemperature,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodyMap.at( bodyWithProperty )->getFlightConditions( ) ) );

            break;
        }
        case estimated_flight_path_angle:
        {

            std::function< double( const Eigen::Vector3d&, const double& ) > estimated_flight_path_angle_def_rad =
                    std::bind( &bislip::Variables::estimatedFlightPathAngle, std::placeholders::_1, std::placeholders::_2 );

            // Retrieve functions for
            std::function< Eigen::Vector3d( ) > firstInput =
                    std::bind( &bislip::Variables::computeAerodynamicFrameTotalLoad, bodyMap, bodyWithProperty );
            std::function< double( ) > secondInput =
                   std::bind( &simulation_setup::Body::getBodyMass, bodyMap.at( bodyWithProperty ) );

            std::function< double( ) > estimated_flight_path_angle_rad =
                    std::bind( &evaluateBivariateReferenceFunction< double, Eigen::Vector3d, double >,
                                          estimated_flight_path_angle_def_rad, firstInput, secondInput );

            std::function< double( const double& ) > rad2deg_Def =
                    std::bind( &bislip::Variables::convertRadiansToDegrees, std::placeholders::_1 );

            variableFunction = std::bind( &evaluateReferenceFunction< double, double >, rad2deg_Def, estimated_flight_path_angle_rad );

            break;
        }

        case trajectory_phase:
        {
            variableFunction = std::bind( &bislip::Variables::convertTrajectoryPhaseToBoolean, bodyMap, bodyWithProperty );

            break;
        }



        default:
            std::string errorMessage =
                    "Error, did not recognize double dependent variable type when making variable function: " +
                    std::to_string( dependentVariableSettings->dependentVariableType_ );
            throw std::runtime_error( errorMessage );
        }
        return variableFunction;
    }
}


//! Function to return a vector containing only one value given by doubleFunction
/*!
 *
 * \param doubleFunction Function returning the double value
 * \return The vector containing only one element
 */
Eigen::VectorXd getVectorFromDoubleFunction( const std::function< double( ) >& doubleFunction );

//! Function to evaluate a set of vector-returning functions and concatenate the results.
/*!
 * Function to evaluate a set of vector-returning functions and concatenate the results.
 * \param vectorFunctionList List of functions returning vector variables (pairs denote function and return vector size)
 * \param totalSize Total size of concatenated vector (used as input for efficiency.
 * \return Concatenated results from input functions.
 */
Eigen::VectorXd evaluateListOfVectorFunctions(
        const std::vector< std::pair< std::function< Eigen::VectorXd( ) >, int > > vectorFunctionList,
        const int totalSize );

//! Function to create a function that evaluates a list of dependent variables and concatenates the results.
/*!
 *  Function to create a function that evaluates a list of dependent variables and concatenates the results.
 *  Dependent variables functions are created inside this function from a list of settings on their required
 *  types/properties.
 *  \param saveSettings Object containing types and other properties of dependent variables.
 *  \param bodyMap List of bodies to use in simulations (containing full environment).
 *  \param stateDerivativeModels List of state derivative models used in simulations (sorted by dynamics type as key)
 *  \return Pair with function returning requested dependent variable values, and list variable names with start entries.
 *  NOTE: The environment and state derivative models need to
 *  be updated to current state and independent variable before computation is performed.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::function< Eigen::VectorXd( ) >, std::map< int, std::string > > createDependentVariableListFunction(
        const std::shared_ptr< DependentVariableSaveSettings > saveSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels =
        std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ) )
{
    // Retrieve list of save settings
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables =
            saveSettings->dependentVariables_;

    // create list of vector parameters
    std::vector< std::pair< std::function< Eigen::VectorXd( ) >, int > > vectorFunctionList;
    std::vector< std::pair< std::string, int > > vectorVariableList;

    for( std::shared_ptr< SingleDependentVariableSaveSettings > variable: dependentVariables )
    {
        std::pair< std::function< Eigen::VectorXd( ) >, int > vectorFunction;
        // Create double parameter
        if( getDependentVariableSaveSize( variable ) == 1 )
        {
#if( BUILD_WITH_ESTIMATION_TOOLS )
            std::function< double( ) > doubleFunction =
                    getDoubleDependentVariableFunction( variable, bodyMap, stateDerivativeModels,
                                                        saveSettings->stateDerivativePartials_ );
#else
            std::function< double( ) > doubleFunction =
                    getDoubleDependentVariableFunction( variable, bodyMap, stateDerivativeModels );
#endif
            vectorFunction = std::make_pair( std::bind( &getVectorFromDoubleFunction, doubleFunction ), 1 );
        }
        // Create vector parameter
        else
        {
#if( BUILD_WITH_ESTIMATION_TOOLS )
            vectorFunction = getVectorDependentVariableFunction(
                        variable, bodyMap, stateDerivativeModels, saveSettings->stateDerivativePartials_ );
#else
            vectorFunction =
                    getVectorDependentVariableFunction( variable, bodyMap, stateDerivativeModels );
#endif
        }
        vectorFunctionList.push_back( vectorFunction );
        vectorVariableList.push_back( std::make_pair( getDependentVariableId( variable ), vectorFunction.second ) );
    }

    // Set list of variable ids/indices in correc otder.
    int totalVariableSize = 0;
    std::map< int, std::string > dependentVariableIds;
    for( std::pair< std::string, int > vectorVariable: vectorVariableList )
    {
        dependentVariableIds[ totalVariableSize ] = vectorVariable.first;
        totalVariableSize += vectorVariable.second;
    }

    // Create function conatenating function results.
    return std::make_pair( std::bind( &evaluateListOfVectorFunctions, vectorFunctionList, totalVariableSize ),
                           dependentVariableIds );
}

extern template std::pair< std::function< Eigen::VectorXd( ) >, std::map< int, std::string > > createDependentVariableListFunction< double, double >(
        const std::shared_ptr< DependentVariableSaveSettings > saveSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > >& stateDerivativeModels );

//extern template std::pair< std::function< Eigen::VectorXd( ) >, int > getVectorDependentVariableFunction< double, double >(
//        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
//        const simulation_setup::NamedBodyMap& bodyMap,
//        const std::unordered_map< IntegratedStateType,
//        std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > >& stateDerivativeModels );

//extern template std::function< double( ) > getDoubleDependentVariableFunction< double, double >(
//        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
//        const simulation_setup::NamedBodyMap& bodyMap,
//        const std::unordered_map< IntegratedStateType,
//        std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > >& stateDerivativeModels );


} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONOUTPUT_H
