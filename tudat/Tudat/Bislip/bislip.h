/*
#ifndef BISLIPVARIABLES_H
#define BISLIPVARIABLES_H
*/
/*
#include <ctime>
#include <sstream>
#include <iomanip>
#include <utility>
#include <cmath>
#include <fstream>
//#include <vector>
#include <stdio.h>      // for sprintf()
//#include <iostream>     // for console output
#include <string>       // for std::string
#include <boost/date_time/posix_time/posix_time.hpp>

//#include <map>
//#include <memory>

//#include <Tudat/SimulationSetup/EnvironmentSetup/body.h>
//#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
//#include <Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Mathematics/Interpolators/createInterpolator.h>
#include <Tudat/Astrodynamics/SystemModels/vehicleSystems.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
*/
/*
#include <Tudat/Bislip/bislipHeaders.h>


namespace bislip {

namespace Parameters {

enum Optimization
{
    AngleOfAttack,
    BankAngle,
    ThrustElevationAngle,
    ThrustAzimuthAngle,
    ThrottleSetting,
    NodeInterval,
    InitialVelocity,
    MaximumVelocity,
    MaximumHeight,
    AdditionalMass
};
};

//! Wrapper class that contains the relevant hardware systems of a vehicle.
/*!
 *  Wrapper class that contains the relevant hardware systems of a vehicle. Not all member objects need to be set; nullptr
 *  member objects denote that the vehicle does not contain the associated hardware.
 */
/*
class VehicleSystems
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param dryMass Total dry mass of the vehicle (not defined; NaN by default).
     *//*
    VehicleSystems( const double landingMass = 0 ):
        landingMass_( landingMass ){ }

    //! Destructor
    ~VehicleSystems( ){ }


    void setInitialMass( const double initialMass )
    {
        initialMass_ = initialMass;
    }

    double getLandingMass( )
    {
        return landingMass_;
    }

    //! Function to retrieve the initial mass of the vehicle
    /*!
     * Function to retrieve the  total dry mass of the vehicle
     * \return initial mass of the vehicle
     *//*
    double getInitialMass( )
    {
        return initialMass_;
    }

    void setReferenceArea( const double referenceArea )
    {
        referenceArea_ = referenceArea;
    }

    double getReferenceArea( )
    {
        return referenceArea_;
    }

    void setMaxThrust( const double maxThrust )
    {
        maxThrust_ = maxThrust;
    }

    double getMaxThrust( )
    {
        return maxThrust_;
    }

    void setSpecificImpulse( const double specificImpulse )
    {
        specificImpulse_ = specificImpulse;
    }

    double getSpecificImpulse( )
    {
        return specificImpulse_;
    }

    void setWingSweepAngle( const double lambda )
    {
        lambda_ = lambda;
    }

    double getWingSweepAngle( )
    {
        return lambda_;
    }

    void setLocalBodyAngle( const double phi )
    {
        phi_ = phi;
    }

    double getLocalBodyAngle( )
    {
        return phi_;
    }

    void setTransitionDistance( const double x_T )
    {
        x_T_ = x_T;
    }

    double getTransitionDistance( )
    {
        return x_T_;
    }

    void setE_max( const double E_max )
    {
        E_max_ = E_max;
    }

    double getE_max( )
    {
        return E_max_;
    }

    void setAngleOfAttackInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack )
    {
        interpolatorAngleOfAttack_ = interpolatorAngleOfAttack;
    }
    void setBankAngleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngle )
    {
        interpolatorBankAngle_ = interpolatorBankAngle;
    }
    void setThrustElevationAngleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngle )
    {
        interpolatorThrustElevationAngle_ = interpolatorThrustElevationAngle;
    }
    void setThrustAzimuthAngleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle )
    {
        interpolatorThrustAzimuthAngle_ = interpolatorThrustAzimuthAngle;
    }
    void setThrottleInterpolator( const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting )
    {
        interpolatorThrottleSetting_ = interpolatorThrottleSetting;
    }

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getAngleOfAttackInterpolator( )
    {
        return interpolatorAngleOfAttack_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getBankAngleInterpolator( )
    {
        return interpolatorBankAngle_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrustElevationAngleInterpolator( )
    {
        return interpolatorThrustElevationAngle_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrustAzimuthAngleInterpolator( )
    {
        return interpolatorThrustAzimuthAngle_;
    }
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > getThrottleInterpolator( )
    {
        return interpolatorThrottleSetting_;
    }

    void setTargetLat( const double target_lat )
    {
        target_lat_ = target_lat;
    }

    void setTargetLon( const double target_lon )
    {
        target_lon_ = target_lon;
    }

    void setInitialLat( const double initial_lat )
    {
        initial_lat_ = initial_lat;
    }
    void setInitialLon( const double initial_lon )
    {
        initial_lon_ = initial_lon;
    }
    void setInitialDistanceToTarget( const double initial_d_to_target )
    {
        initial_d_to_target_ = initial_d_to_target;
    }
    void setFinalDistanceToTarget( const double final_d_to_target )
    {
        final_d_to_target_ = final_d_to_target;
    }
    double getTargetLat( )
    {
        return target_lat_;
    }
    double getTargetLon( )
    {
        return target_lon_;
    }
    double getInitialLat( )
    {
        return initial_lat_;
    }
    double getInitialLon( )
    {
        return initial_lon_;
    }
    double getInitialDistanceToTarget( )
    {
        return initial_d_to_target_;
    }
    double getFinalDistanceToTarget( )
    {
        return final_d_to_target_;
    }
    void setParameterBounds( const std::map< bislip::Parameters::Optimization, std::pair < double, double > > &Bounds )
    {

        Bounds_ = Bounds;
        //   if( direction == "ascent" ){ Bounds_Ascent_ = parameterBounds; }
        // if( direction == "descent" ){ Bounds_Descent_ = parameterBounds; }
    }

    //    std::pair < double, double > getParameterBounds( const std::string &parameter )
    std::pair < double, double > getParameterBounds( const bislip::Parameters::Optimization &parameter )
    {
        // std::map< std::string, std::pair < double, double > > Bounds;

        //        if( direction == "ascent" ){ Bounds = Bounds_Ascent_; }
        //       if( direction == "descent" ){ Bounds = Bounds_Descent_; }

        return Bounds_[ parameter ];
    }

    void setInitialCoordinates( const std::pair < double, double > &initialCoordinates )
    {
        initialCoordinates_ = initialCoordinates;
    }

    std::pair < double, double > getInitialCoordinates( )
    {
        return initialCoordinates_;
    }

    void setTargetCoordinates( const std::pair < double, double > &targetCoordinates )
    {
        targetCoordinates_ = targetCoordinates;
    }

    std::pair < double, double > getTargetCoordinates( )
    {
        return targetCoordinates_;
    }
    void setStartingEpoch( const double startingEpoch )
    {
        startingEpoch_ = startingEpoch;
    }
    double getStartingEpoch( )
    {
        return startingEpoch_;
    }

    void setReferenceValues( const Eigen::Vector3d referenceValues )
    {
        referenceValues_ = referenceValues;
    }
    Eigen::Vector3d getReferenceValues( )
    {
        return referenceValues_;
    }
    void setMomentReferenceCenter( const Eigen::Vector3d momentReferenceCenter )
    {
        momentReferenceCenter_ = momentReferenceCenter;
    }
    Eigen::Vector3d getMomentReferenceCenter( )
    {
        return momentReferenceCenter_;
    }
    void setMassReferenceCenter( const Eigen::Vector3d massReferenceCenter )
    {
        massReferenceCenter_ = massReferenceCenter;
    }
    Eigen::Vector3d getMassReferenceCenter( )
    {
        return massReferenceCenter_;
    }
    void setThrustReferenceCenter( const Eigen::Vector3d thrustReferenceCenter )
    {
        thrustReferenceCenter_ = thrustReferenceCenter;
    }
    Eigen::Vector3d getThrustReferenceCenter( )
    {
        return thrustReferenceCenter_;
    }





private:


    std::pair < double, double > initialCoordinates_;
    std::pair < double, double > targetCoordinates_;

    double lambda_;
    double phi_;
    double x_T_;
    double E_max_;
    double maxThrust_;
    double specificImpulse_;
    double startingEpoch_;
    double referenceArea_;
    // double initialMass_;

    //    std::map< std::string, std::pair < double, double > > Bounds_;
    std::map< std::string, std::pair < double, double > > Bounds_Ascent_;
    std::map< std::string, std::pair < double, double > > Bounds_Descent_;
    std::map< bislip::Parameters::Optimization, std::pair < double, double > > Bounds_;

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngle_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngle_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_;


    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngleAscent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_Ascent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngleAscent_;

    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorAngleOfAttack_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrottleSetting_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustElevationAngleDescent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorThrustAzimuthAngle_Descent_;
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolatorBankAngleDescent_;
    //! Passing target coordinates
    double target_lat_;
    double target_lon_;

    //! Passing initial coordinates
    double initial_lat_;
    double initial_lon_;

    //! Passing initial distance to target
    double initial_d_to_target_;

    //! Passing final distance to target
    double final_d_to_target_;

    //! Initial mass of the vehicle
    double initialMass_;

    //! Landing mass of the vehicle
    double landingMass_;


    double b_ref_;
    double c_ref_;
    double del_x_;
    double del_x_T_;
    double del_z_T_;
    Eigen::Vector3d referenceValues_;
    Eigen::Vector3d momentReferenceCenter_;
    Eigen::Vector3d thrustReferenceCenter_;
    Eigen::Vector3d massReferenceCenter_;



};

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

//std::string passOptimizationParameter (
//        const std::string &parameter);

bislip::Parameters::Optimization passOptimizationParameter (
        const bislip::Parameters::Optimization &parameter);

//std::string passDirection (
//      const std::string &direction);

std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > chooseGuidanceInterpolator (
        const bislip::Parameters::Optimization &parameter,
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems);

std::pair < double, double > chooseGuidanceBounds (
        const bislip::Parameters::Optimization &parameter,
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems);

double evaluateGuidanceInterpolator (
        const bislip::Parameters::Optimization &parameter,
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems,
        const double &height,
        const double &airspeed,
        const double &E_max);

Eigen::Vector6d computeCurrentCoefficients (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems);

Eigen::Vector3d computeBodyFixedThrustDirection (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems);

double computeThrustMagnitude (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems);

Eigen::Vector3d computeBodyFixedThrustVector (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems);

bool determineEngineStatus (
        const double &currentMass,
        const double &landingMass);

double getBodyFlapDeflection( const double &del_C_m_b );

Eigen::Vector2d computeLocalGravity (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions );

double computeEquilibriumGlideLimit (const tudat::simulation_setup::NamedBodyMap& bodyMap );
//        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
  //      const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
    //    const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems,
      //  const double &currentMass );

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
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems);

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
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions);

double computePenalty (
        const Eigen::VectorXd &dependentVariable_TimeHistory,
        const long &startIterator,
        const long &endIterator,
        const double &constraint, const double &fixedStepSize, const double &tof, const bool &direct );

double computeBodyflapCmIncrement (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
        const std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems );

//bool StopOrNot (const tudat::simulation_setup::NamedBodyMap& bodyMap,
//                const std::string &vehicleName,
//                const std::vector< double > &vehicleParameterValues,
//                const std::vector< double > &terminationConditionsValues);


}; // namespace Variables


} // namespace bislip

//pagmo::algorithm getPagmoAlgorithm ( const int index );

#endif // BISLIPVARIABLES_H
*/
