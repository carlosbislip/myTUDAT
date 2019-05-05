#ifndef BISLIPVARIABLES_H
#define BISLIPVARIABLES_H

#include <Tudat/Bislip/bislipVehicleSystems.h>
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/secantRootFinder.h"
#include "Tudat/Mathematics/RootFinders/bisection.h"


namespace bislip {

namespace Variables
{

std::string getOutputPath(
        const std::string& extraDirectory = "" );

Eigen::MatrixXd convertVectorOfVectorsDoubleToEigenMatrixXd(
        const std::vector< std::vector< double > > &vectorOfVectorsDouble );

void printEigenMatrixXdToFile( const Eigen::MatrixXd &matrixToPrint,
                               const std::string &fileName,
                               const std::string &outputSubFolder );

unsigned int millis_since_midnight ( );

std::string getCurrentDateTime ( const bool useLocalTime = false );

std::vector< std::string > getDataString ( const std::string &filename );

std::vector< double > getDataNumeri ( const std::string &filename );

double rootFinderBisection (
        const std::function< double( const double ) > &function,
        const double &minimum,
        const double &maximum );

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
        const Eigen::VectorXd &xValues,
        const Eigen::VectorXd &yValues);

std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > createOneDimensionalHermiteInterpolator (
        const Eigen::VectorXd &yValues,
        const Eigen::VectorXd &xValues,
        const std::map< double, double > &mapped_data,
        const std::shared_ptr< tudat::interpolators::InterpolatorSettings > &interpolatorSettings );

bislip::Parameters::Optimization passOptimizationParameter (
        const bislip::Parameters::Optimization &parameter);

std::string passString (
        const std::string &string);

tudat::simulation_setup::NamedBodyMap& passBodyMap (
        tudat::simulation_setup::NamedBodyMap& bodyMap );
/*
std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > chooseGuidanceInterpolator (
        const bislip::Parameters::Optimization &parameter,
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems);

std::pair < double, double > chooseGuidanceBounds (
        const bislip::Parameters::Optimization &parameter,
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems);
*/
double evaluateGuidanceInterpolator (
        const bislip::Parameters::Optimization &parameter,
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

double computeThrustMagnitudeOutput (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

Eigen::Vector2d computeLocalGravity (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName );

double computeCurrentLiftForce (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName);

double computeSkipSuppressionLimit (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName);

double computeFlightPathAngleRate(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName );

double computeHeatingRate (
        const double &airdensity,
        const double &airspeed,
        const double &C,
        const double &N,
        const double &M);

double computeHeatingRateTauber (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeHeatingRateChapman (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeRadiativeHeatFlux (
        const double &wallEmissivity,
        const double &bodyTemperature);

double computeStagnationHeatFlux (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

std::function< double( const double ) > getStagnationHeatTransferFunction (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );


double computeAdiabaticWallTemperature(
        const double airTemperature, const double machNumber );


double computeFlatPlateHeatFlux (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double findEquilibriumWallTemperature (
        const std::function< double( const double ) > &heatTransferFunction,
        const double &wallEmissivity,
        const double &adiabaticWallTemperature );

double computeEquilibiumWallTemperatureRootFinder (
        const std::function< double( const double ) > &heatTransferFunction,
        const double &wallEmmisivity,
        const double &wallTemperature );
//        const std::function< double( const double ) > heatTransferFunction,
//        const double wallEmmisivity,
//        const double adiabaticWallTemperature,
//        const tudat::simulation_setup::NamedBodyMap& bodyMap,
//        const std::string &vehicleName );

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

double computeCumulativeCartesianDistanceTravelled (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeCumulativeAngularDistanceTravelled (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeCumulativeAngularDistanceTravelledDifference (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeTimeOfFlight (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeBendingMoment (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computePenalty (
        const Eigen::VectorXd &dependentVariable_TimeHistory,
        const long &startIterator,
        const long &endIterator,
        const double &constraint, const double &fixedStepSize, const double &tof, const bool &direct );

double computeMonotonicPenalty (
        const Eigen::VectorXd &depVar_TimeHistory,
        const std::string &category );

double computeCompoundViolationPenalty(
        const Eigen::VectorXd &depVar_TimeHistory,
        const double &constraint,
        const double &propagationStepSize,
        const double &normalizer );

double computeConstraintViolationPenalty(
        const Eigen::VectorXd &depVar_TimeHistory,
        const double &constraint,
        const double &propagationStepSize,
        const double &normalizer );

Eigen::VectorXd computeConstraintViolations (
        const Eigen::VectorXd &dependentVariable_TimeHistory,
        const double &constraint );

double determineMaximumValueFromEigenVectorXd (
        const Eigen::VectorXd &vectorXd );

double computeBodyFlapCmIncrement (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

Eigen::Vector2d computeControlSurfaceDeflection (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeFullPitchMomentCoefficient (
        const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
        const double &bodyFlapDeflection,
        const double &elevonDeflection,
        const std::vector< double > &aerodynamicCoefficientInput );

double computePitchMomentCoefficient (
        const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
        const std::vector< double > &aerodynamicCoefficientInput );

double computeBodyFlapCmIncrementdif (
        const double &full,
        const double &partial );

Eigen::Vector6d computePartialCurrentCoefficients (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

Eigen::Vector6d computeFullCurrentCoefficients (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

std::vector< double > getAerodynamicCoefficientInput (
        const double &angleOfAttack,
        const double &machNumber );

std::map< std::string, std::vector< double > > getControlSurfaceCoefficientInput (
        const double &angleOfAttack,
        const double &machNumber,
        const double &bodyFlapDeflectionAngle,
        const double &elevonDeflectionAngle );

Eigen::Vector3d computeBodyFixedTotalLoad (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

Eigen::Vector3d computeBodyFixedAerodynamicLoad (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

Eigen::Vector3d computeBodyFixedTotal_g_Load_Vector (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeBodyFixedTotal_g_Load_Magnitude (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

Eigen::Matrix3d computeRotationMatrixONE( const double phi );

Eigen::Matrix3d computeRotationMatrixTHREE( const double psi );

Eigen::Vector3d computePassengerFrameTotal_g_Load_Vector (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeBankAngle (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double returnReversedBankAngle(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const double &newBankAngle );

double determineNewBankAngle(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double computeHeadingErrorDeadBand(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

double convertRadiansToDegrees (
        const double &angleInRadians );

double convertDegreesToRadians (
        const double &angleInDegrees );

double convertNegativeAnglesToPositive (
        const double &angleInRadians );

double determineAbsoluteValue (
        const double &value );

int determineBankAngleSign (
        const double &currentBankAngle );

double determineSignedBankAngle (
        const int &sign,
        const double &newbankAngle );

double determineReversalConditional (
        const double &currentBankAngle,
        const double &currentHeadingError );

bool determineBankAngleReversal (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const double &newBankAngle );

void createAlphaMachBoundingInterpolators (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::vector< double > &alphaMachEnvelopeUB,
        const std::vector< double > &alphaMachEnvelopeLB,
        const std::string &outputPath,
        const std::string &outputSubFolder );

void printAlphaMachBounds (
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &alphaMachEnvelopeUBInterpolator,
        const std::pair< double, double > &alphaMachEnvelopeUBInterpolatorDomainInterval,
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &alphaMachEnvelopeLBInterpolator,
        const std::pair< double, double > &alphaMachEnvelopeLBInterpolatorDomainInterval,
        const std::string &outputPath,
        const std::string &outputSubFolder );

void createHeadingErrorDeadBandInterpolator (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::vector< double > &headingErrorDeadBand,
        const std::string &outputPath,
        const std::string &outputSubFolder );

void printHeadingErrorDeadBandBounds (
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &headingErrorDeadBandInterpolator,
        const std::pair< double, double > &headingErrorDeadBandInterpolatorDomainInterval,
        const std::string &outputPath,
        const std::string &outputSubFolder );

void createKourouGuidanceInterpolators(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::vector< double > &kourouAngleofAttackHistory,
        const std::vector< double > &kourouBankAngleHistory,
        const std::string &outputPath,
        const std::string &outputSubFolder );

void printKourouGuidanceInterpolators(
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &kourouAngleOfAttackInterpolator,
        const std::pair< double, double > &kourouAngleOfAttackInterpolatorDomainInterval,
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &kourouBankAngleInterpolator,
        const std::pair< double, double > &kourouBankAngleInterpolatorDomainInterval,
        const std::string &outputPath,
        const std::string &outputSubFolder );

double computeSumOfEigenVectorXd (
        const Eigen::VectorXd &vector );

Eigen::VectorXd computeElementwiseDivisionOfEigenVectorXd (
        const Eigen::VectorXd numerator,
        const Eigen::VectorXd denominator );


/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


//! Function that is used to compute the net heat flux from a given heat input and wall temperature
class EquilibriumTemperatureFunction: public tudat::basic_mathematics::Function< double,double >
{
public:
    //! Constructor.
    /*!
     * Constructor
     * \param heatTransferFunction Function returning the feat flux as a function of wall temperature.
     * \param wallEmissivity Emmissivity of the wall to which heat transfer is taking place
     * \param adiabaticWallTemperature Adiabatic wall temperature
     */
    EquilibriumTemperatureFunction(
            const std::function< double( const double ) > heatTransferFunction,
            const double wallEmissivity,
            double adiabaticWallTemperature,
            std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems ):
        heatTransferFunction_( heatTransferFunction ), wallEmissivity_( wallEmissivity ),
        adiabaticWallTemperature_( adiabaticWallTemperature ), bislipSystems_( bislipSystems ){ }

    //! Destructor.
    ~EquilibriumTemperatureFunction(){}

    //! Compute net heat flux at given wall temperature
    /*!
     * Compute net heat flux at given wall temperature
     * \param currentWallTemperature Wall temperature to be used for computation of input and output of heat.
     * \return Net heat input to wall
     */
    double evaluate( const double currentWallTemperature )
    {

        int debugInfo = bislipSystems_->getDebugInfo();

        if( debugInfo == 1 ){ std::cout << "Wall Temp:" << currentWallTemperature << "----->  ( " << heatTransferFunction_( currentWallTemperature ) << " ) - ( " << bislip::Variables::computeRadiativeHeatFlux( wallEmissivity_, currentWallTemperature ) << " ) = " << heatTransferFunction_( currentWallTemperature )
                                           - bislip::Variables::computeRadiativeHeatFlux( wallEmissivity_, currentWallTemperature ) << std::endl; }


        return heatTransferFunction_( currentWallTemperature )
                - bislip::Variables::computeRadiativeHeatFlux( wallEmissivity_, currentWallTemperature );

    }

    //! Compute first derivative of net heat flux at given wall temperature (FUNCTION NOT IMPLEMENTED)
    double computeDerivative( const unsigned int order, const double independentVariable )
    {
        throw std::runtime_error( "Error, derivative of heat flux not defined" );
        return TUDAT_NAN;
    }

    //! Compute first derivative of net heat flux at given wall temperature (FUNCTION NOT IMPLEMENTED)
    double computeDefiniteIntegral( const unsigned int order, const double lowerBound, const double upperbound )
    {
        throw std::runtime_error( "Error, integrall of heat flux not defined" );
        return TUDAT_NAN;
    }

    //! Function to retrieve the lower bound of the wall temperature.
    double getLowerBound( ) { return 0.0; }

    //! Function to retrieve the upper bound of the wall temperature.
    double getUpperBound( ) { return adiabaticWallTemperature_; }

    //! Function to retrieve the initial guess of the wall temperature.
    double getInitialGuess( )
    {
        double initialGuess = adiabaticWallTemperature_ * 0.01;
        if ( bislipSystems_->getWallTemperature() != 0.0 ) { initialGuess = bislipSystems_->getWallTemperature(); }

        int debugInfo = bislipSystems_->getDebugInfo();

        if( debugInfo == 1 )
        {
            std::cout << "Initial guess: bislipSystems_->getWallTemperature() = " << bislipSystems_->getWallTemperature() << std::endl;
            std::cout << "Initial guess: adiabaticWallTemperature_ * 0.01 = " << adiabaticWallTemperature_ * 0.01 << std::endl;
        }
        // std::cout << "adiabaticWallTemperature_ * 0.01 = " << adiabaticWallTemperature_ * 0.01 << std::endl;

        return initialGuess;
    }

protected:

private:
    //! Function that returns the heat input as a function of wall temperature.
    std::function< double( const double ) > heatTransferFunction_;

    //! Constant wall emissivity.
    const double wallEmissivity_;

    //! Constant adiabatic wall temperature.
    double adiabaticWallTemperature_;

    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems_;
};

//! Function to compute the equilibrium wall temperature from the heat input and emmisivity
/*!
 * Function to compute the equilibrium wall temperature from the heat input and emmisivity. This function calls a root-
 * finder to determine the wall temperature where the heat input equals the radiative output.
 * \param heatTransferFunction Function that returns the heat input as a function of current temperature
 * \param wallEmmisivity Value of the material emmisivity
 * \param adiabaticWallTemperature Adiabatic wall temperature. This variables is only used as an upper bound, and to
 * generate an initial giess for the equilibrium temperature
 * \return Wall temperature at which input and output of heat are in equilibrium.
 */
double computeEquilibiumWallTemperature(
        const std::function< double( const double ) > heatTransferFunction,
        const double wallEmmisivity,
        const double adiabaticWallTemperature,
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName );

bool customTermination_FlightPathAngleCombination_Ascent(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName);

bool customTermination_FlightPathAngleCombination_Descent(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName);

double computeDeadbandViolationPenalty(
        const Eigen::VectorXd &depVar_TimeHistory,
        const Eigen::VectorXd &depVar_Constraint,
        const double &propagationStepSize,
        const double &normalizer );


//} //namespace_aerodynamics

//} //namespace_tudat

//#endif //TUDAT_EQUILIBRIUMWALLTEMPERATURE_H

//bool StopOrNot (const tudat::simulation_setup::NamedBodyMap& bodyMap,
//                const std::string &vehicleName,
//                const std::vector< double > &vehicleParameterValues,
//                const std::vector< double > &terminationConditionsValues);


} // namespace Variables
} // namespace bislip

#endif // BISLIPVARIABLES_H
