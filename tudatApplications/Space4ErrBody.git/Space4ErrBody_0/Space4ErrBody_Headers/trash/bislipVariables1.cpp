#include <iostream>
#include <ctime>
#include <cmath>
#include <sstream>
#include <utility>
#include <fstream>
#include <iomanip>
#include <vector>

#include <stdio.h>      // for sprintf()

#include <iostream>     // for console output
#include <string>       // for std::string

#include <boost/date_time/posix_time/posix_time.hpp>



#include <Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>
#include <Tudat/Mathematics/Interpolators/createInterpolator.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/SystemModels/vehicleSystems.h>
#include <Tudat/Astrodynamics/SystemModels/bislipSystems.h>
#include <Tudat/Astrodynamics/SystemModels/bislipOptimizationParameters.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Astrodynamics/Aerodynamics/equilibriumWallTemperature.h>
#include <Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h>

#include "updateGuidance.h"

namespace bislip {

namespace variables {

std::string getCurrentDateTime( bool useLocalTime ) {
    std::stringstream currentDateTime;

    // current date/time based on current system
    time_t ttNow = time(nullptr);
    tm * ptmNow;

    if (useLocalTime)
        ptmNow = localtime(&ttNow);
    else
        ptmNow = gmtime(&ttNow);

    currentDateTime << 1900 + ptmNow->tm_year;

    //month
    if (ptmNow->tm_mon < 9)
        //Fill in the leading 0 if less than 10
        currentDateTime << "0" << 1 + ptmNow->tm_mon;
    else
        currentDateTime << (1 + ptmNow->tm_mon);

    //day
    if (ptmNow->tm_mday < 10)
        currentDateTime << "0" << ptmNow->tm_mday << "_";
    else
        currentDateTime <<  ptmNow->tm_mday << "_";

    //hour
    if (ptmNow->tm_hour < 10)
        currentDateTime << "0" << ptmNow->tm_hour;
    else
        currentDateTime << ptmNow->tm_hour;

    //min
    if (ptmNow->tm_min < 10)
        currentDateTime << "0" << ptmNow->tm_min;
    else
        currentDateTime << ptmNow->tm_min;

    //sec
    if (ptmNow->tm_sec < 10)
        currentDateTime << "0" << ptmNow->tm_sec;
    else
        currentDateTime << ptmNow->tm_sec;
    /*
    // Get current time from the clock, using microseconds resolution
        const boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();

        // Get the time offset in current day
        const boost::posix_time::time_duration td = now.time_of_day();

        //
        // Extract hours, minutes, seconds and milliseconds.
        //
        // Since there is no direct accessor ".milliseconds()",
        // milliseconds are computed _by difference_ between total milliseconds
        // (for which there is an accessor), and the hours/minutes/seconds
        // values previously fetched.
        //
        const long hours        = td.hours();
        const long minutes      = td.minutes();
        const long seconds      = td.seconds();
        const long milliseconds = td.total_milliseconds() -
                                  ((hours * 3600 + minutes * 60 + seconds) * 1000);

        //
        // Format like this:
        //
        //      hh:mm:ss.SSS
        //
        // e.g. 02:15:40:321
        //
        //      ^          ^
        //      |          |
        //      123456789*12
        //      ---------10-     --> 12 chars + \0 --> 13 chars should suffice
        //
        //
        char buf[40];
        sprintf(buf, "%02ld:%02ld:%02ld.%03ld",
            hours, minutes, seconds, milliseconds);
*/
    //std::cout << "Got time: " << currentDateTime.str() << std::endl;


    //unsigned int derp = millis_since_midnight ( );

    //    return currentDateTime.str();

    return std::to_string( millis_since_midnight ( ) );
}


unsigned int millis_since_midnight ( )
{
    // current time
    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();

    // get midnight
    time_t tnow = std::chrono::system_clock::to_time_t(now);
    tm *date = std::localtime(&tnow);
    date->tm_hour = 0;
    date->tm_min = 0;
    date->tm_sec = 0;
    auto midnight = std::chrono::system_clock::from_time_t(std::mktime(date));

    // number of milliseconds between midnight and now, ie current time in millis
    // The same technique can be used for time since epoch
    return std::chrono::duration_cast<std::chrono::milliseconds>(now - midnight).count();
}


std::vector< std::string > getDataString ( const std::string &filename )
{
    std::ifstream inputdata;
    inputdata.open( filename.c_str( ) );
    if( inputdata.fail( ) )
    {
        std::cout << "File '" << filename << "' failed to open";
        std::cin.ignore( 100, '\n' );
        std::exit( 1 );
    }
    else{

        std::string  var;
        std::vector< std::string > stuff;

        while( std::getline( inputdata, var ) )
        {
            if( var.size( ) > 0 )
            {
                stuff.push_back( var );
            }
        }
        inputdata.close( );
        return stuff;
    }
}

std::vector< double > getDataNumeri ( const std::string &filename )
{
    std::ifstream inputdata;
    inputdata.open( filename.c_str( ) );
    if( inputdata.fail( ) )
    {
        std::cout << "File '" << filename << "' failed to open";
        std::cin.ignore( 100, '\n' );
        std::exit( 1 );
    }
    else{

        double var;
        std::vector< double > stuff;

        while (!inputdata.fail() && !inputdata.eof())
        {
            inputdata >> var;
            stuff.push_back(var);
        }

        inputdata.close( );
        return stuff;
    }
}

double computeAngularDistance (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f)
{
    //const double angular_distance = std::acos( std::sin( lat_c ) * std::sin( lat_f ) + std::cos( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );

    return std::acos( std::sin( lat_c ) * std::sin( lat_f ) + std::cos( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );
}

double computeHeadingToTarget (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f)
{
    // const double heading_to_target = std::atan2( std::sin( lon_f - lon_c ) * std::cos( lat_f ) , std::cos( lat_c ) * std::sin( lat_f ) - std::sin( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );

    return std::atan2( std::sin( lon_f - lon_c ) * std::cos( lat_f ) , std::cos( lat_c ) * std::sin( lat_f ) - std::sin( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );

}

double computeHeadingError (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f,
        const double &heading)
{
    // const double heading_to_target = computeHeadingToTarget( lat_c, lon_c, lat_f, lon_f );

    return computeHeadingToTarget( lat_c, lon_c, lat_f, lon_f ) - heading;
}

double computeSpecificEnergy (
        const double &height,
        const double &airspeed)
{
    return tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * height + 0.5 * airspeed * airspeed;
}

double computeNormalizedSpecificEnergy (
        const double &height,
        const double &airspeed,
        const double &E_max)
{
    return computeSpecificEnergy( height, airspeed ) / E_max;
}

//! https://doi.org/10.1016/j.cam.2017.09.049
std::vector< double > HermiteDerivatives( const Eigen::VectorXd &mappedNormalizedSpecificEnergy, const Eigen::VectorXd &y )
{
    long nodes = mappedNormalizedSpecificEnergy.size();
    Eigen::VectorXd h( nodes - 1 ), dely( nodes - 1 ), b( nodes ), x( nodes );
    Eigen::MatrixXd A( nodes, nodes );
    std::vector< double > x_vect;
    double mu, lambda;

    for ( long i = 0; i < nodes - 1; ++i ) { h( i ) = mappedNormalizedSpecificEnergy( i + 1 ) - mappedNormalizedSpecificEnergy( i ); dely( i ) = ( y( i + 1 ) - y( i ) ) / h( i ); }

    A = Eigen::MatrixXd::Zero( nodes, nodes );
    A( 0, 0 ) = 4.0;
    A( 0, 1 ) = -1.0;
    A( nodes - 1, nodes - 2 ) = -1.0;
    A( nodes - 1, nodes - 1 ) = 4.0;

    for ( long i = 1; i < nodes - 1; ++i )
    {
        lambda = h( i ) / ( h( i - 1 ) + h( i ) );
        mu = 1 - lambda;
        A( i , i - 1 ) = -mu;
        A( i , i ) = 4.0;
        A( i , i + 1 ) = -lambda * mu;
    }

    b( 0 ) = dely( 0 );
    b( nodes - 1 ) = dely( nodes - 2 );

    for ( long i = 1; i < nodes - 1; ++i ) { b( i ) = 3 * ( y( i + 1 ) - y( i - 1 ) ) / ( h( i - 1 ) + h( i ) ); }

    x = A.fullPivHouseholderQr().solve( b );

    for ( long i = 0; i < nodes; ++i ) { x_vect.push_back( x( i ) ); }

    return x_vect;
}

std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > createOneDimensionalHermiteInterpolator (
        const Eigen::VectorXd &parameterValues,
        const Eigen::VectorXd &mappedNormalizedSpecificEnergy,
        const std::map< double, double > &mapped_data,
        const std::shared_ptr< tudat::interpolators::InterpolatorSettings > &interpolatorSettings )
{
    return tudat::interpolators::createOneDimensionalInterpolator< double, double >(
                mapped_data,
                interpolatorSettings,
                std::make_pair( parameterValues( 0 ), parameterValues( parameterValues.size() - 1 ) ),
                bislip::variables::HermiteDerivatives( mappedNormalizedSpecificEnergy, parameterValues ) );
}

//std::string passOptimizationParameter (
//        const std::string &parameter)
//{
//    return parameter;
///}

bislip::variables::OptimizationParameter passOptimizationParameter (
        const bislip::variables::OptimizationParameter &parameter)
{
    return parameter;
}



std::string passDirection (
        const std::string &direction)
{
    return direction;
}


std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > chooseGuidanceInterpolator(
        const bislip::variables::OptimizationParameter &parameter,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems)
{
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator;

    if ( parameter == bislip::variables::OptimizationParameter::AngleOfAttack ) { interpolator = bislipSystems->getAngleOfAttackInterpolator(); }
    if ( parameter == bislip::variables::OptimizationParameter::BankAngle ) { interpolator = bislipSystems->getBankAngleInterpolator(); }
    if ( parameter == bislip::variables::OptimizationParameter::ThrustElevationAngle ) { interpolator = bislipSystems->getThrustElevationAngleInterpolator(); }
    if ( parameter == bislip::variables::OptimizationParameter::ThrustAzimuthAngle ) { interpolator = bislipSystems->getThrustAzimuthAngleInterpolator(); }
    if ( parameter == bislip::variables::OptimizationParameter::ThrottleSetting ) { interpolator = bislipSystems->getThrottleInterpolator(); }

    return interpolator;
}


std::pair < double, double > chooseGuidanceBounds (
        const bislip::variables::OptimizationParameter &parameter,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems)
{
    std::pair < double, double > bounds;

    if ( parameter == AngleOfAttack ) { bounds = bislipSystems->getParameterBounds( parameter ); }
    if ( parameter == BankAngle ) { bounds = bislipSystems->getParameterBounds( parameter ); }
    if ( parameter == ThrustElevationAngle ) { bounds = bislipSystems->getParameterBounds( parameter ); }
    if ( parameter == ThrustAzimuthAngle ) {  bounds = bislipSystems->getParameterBounds( parameter ); }
    if ( parameter == ThrottleSetting ) { bounds = bislipSystems->getParameterBounds( parameter ); }

    return bounds;
}

double evaluateGuidanceInterpolator (
        const bislip::variables::OptimizationParameter &parameter,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems,
        const double &height,
        const double &airspeed,
        const double &E_max)
{
    //! Select parameter interpolator.
    //std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator = bislip::variables::chooseGuidanceInterpolator( parameter, vehicleSystems );

    //! Evaluate interpolator.
    //double evaluation = interpolator->interpolate( bislip::variables::computeNormalizedSpecificEnergy( height, airspeed, E_max ) );
    double evaluation = ( bislip::variables::chooseGuidanceInterpolator( parameter, bislipSystems ) )->interpolate( bislip::variables::computeNormalizedSpecificEnergy( height, airspeed, E_max ) );

    //! Select parameter bounds.
    std::pair < double, double > bounds =  bislip::variables::chooseGuidanceBounds( parameter, bislipSystems );

    //! Impose bounds.
    if ( evaluation < bounds.first ){ evaluation = bounds.first; }
    if ( evaluation > bounds.second ){ evaluation = bounds.second; }

    return evaluation;
}

Eigen::Vector3d computeBodyFixedThrustDirection (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems)
{
    // std::cout << "Computing Body Fixed Thrust Direction" << std::endl;
    double eps_T = tudat::unit_conversions::convertDegreesToRadians(
                evaluateGuidanceInterpolator (
                    ThrustElevationAngle,
                    bislipSystems,
                    flightConditions->getCurrentAltitude(),
                    flightConditions->getCurrentAirspeed(),
                    bislipSystems->getE_max() ) );

    double phi_T = tudat::unit_conversions::convertDegreesToRadians(
                evaluateGuidanceInterpolator (
                    ThrustAzimuthAngle,
                    bislipSystems,
                    flightConditions->getCurrentAltitude(),
                    flightConditions->getCurrentAirspeed(),
                    bislipSystems->getE_max() ) );

    Eigen::Vector3d bodyFixedThrustDirection;
    bodyFixedThrustDirection( 0 ) = std::cos( eps_T ) * std::cos( phi_T );
    bodyFixedThrustDirection( 1 ) = std::cos( eps_T ) * std::sin( phi_T );
    bodyFixedThrustDirection( 2 ) = std::sin( eps_T );

    return bodyFixedThrustDirection;
}

double computeThrustMagnitude (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems)
{
    const double throttle = bislip::variables::evaluateGuidanceInterpolator (
                //flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle ),
                ThrottleSetting,
                bislipSystems,
                flightConditions->getCurrentAltitude(),
                flightConditions->getCurrentAirspeed(),
                bislipSystems->getE_max() );

    return throttle * bislipSystems->getMaxThrust();
}

Eigen::Vector3d computeBodyFixedThrustVector (

        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems)

//const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
// const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems)
{
    const Eigen::Vector3d BodyFixedThrustDirection = bislip::variables::computeBodyFixedThrustDirection( flightConditions, bislipSystems );
    const double thrustMagnitude = bislip::variables::computeThrustMagnitude( flightConditions, bislipSystems );

    Eigen::Vector3d BodyFixedThrustVector;

    BodyFixedThrustVector( 0 ) = thrustMagnitude * BodyFixedThrustDirection( 0 );
    BodyFixedThrustVector( 1 ) = thrustMagnitude * BodyFixedThrustDirection( 1 );
    BodyFixedThrustVector( 2 ) = thrustMagnitude * BodyFixedThrustDirection( 2 );

    return BodyFixedThrustVector;
}

bool determineEngineStatus (
        const double &currentMass,
        const double &landingMass)
{
    bool currentEngineStatus = true;
    if ( currentMass <= landingMass ){ currentEngineStatus = false; }
    return currentEngineStatus;
}

Eigen::Vector2d getGravs (
        const double &r,
        const double &latitude )
{

    const double R_E = 6.378137e6;
    const double mu = 3.986004418e14;
    const double J2 = 1082.626523e-6;
    const double J3 = 2.532153e-7;
    const double J4 = 1.6109876e-7;

    const double s_latitude = std::sin( latitude );
    const double c_latitude = std::cos( latitude );

    const double g_n_pre = -3 * ( mu / ( r * r ) ) * pow( R_E / r , 2 ) * s_latitude * c_latitude;
    const double g_n_sub_1 = J2;
    const double g_n_sub_2 = ( J3 / 2 ) * ( R_E / r ) * ( 1 / std::sin( latitude) ) * ( 5 * pow ( s_latitude , 2 ) - 1 );
    const double g_n_sub_3 = ( 5 * J4 / 6 ) * pow( R_E / r , 2) * ( 7 * pow ( s_latitude , 2 ) - 3 );
    const double g_n = g_n_pre * ( g_n_sub_1 + g_n_sub_2 + g_n_sub_3 );

    const double g_d_pre = ( mu / ( r * r ) );
    const double g_d_sub_1 = -( 3 * J2 / 2 ) * pow( R_E / r , 2 ) * ( 3 * pow  (s_latitude , 2 ) - 1 );
    const double g_d_sub_2 = -2 * J3  * pow( R_E / r , 3 ) * ( s_latitude ) * ( 5 * pow ( s_latitude , 2 ) - 3 );
    const double g_d_sub_3 = -( 5 * J4 / 8 ) * pow( R_E / r , 4 ) * ( 35 * pow ( s_latitude , 4 ) - 30 * pow ( s_latitude , 2 ) + 3 );
    const double g_d = g_d_pre * ( 1 + g_d_sub_1 + g_d_sub_2 + g_d_sub_3 );

    Eigen::Vector2d gravs ( 2 );
    gravs << g_n, g_d;

    return gravs;
}

double computeEquilibriumGlideLimit (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions )
{

    const double currentFlightPathAngle = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
    const double currentLatitude = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double currentHeading = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
    const double currentAltitude = flightConditions->getCurrentAltitude( );
    const double currentDensity = flightConditions->getCurrentDensity( );
    const double currentAirspeed = flightConditions->getCurrentAirspeed();
    const double S_ref = bodyMap_.at( vehicleName_ )->getAerodynamicCoefficientInterface( )->getReferenceArea( );
    const double currentMass = bodyMap_.at( vehicleName_ )->getBodyMass( );
    const double q_dyn = ( 0.5 ) * currentDensity * currentAirspeed * currentAirspeed;





    Eigen::Vector2d gravs ( 2 );
    gravs = bislip::variables::getGravs( currentAltitude, currentLatitude );

    const double maxBankAngle = ( currentMass / ( q_dyn * S_ref * currentCoefficients[ 2 ] ) ) * ( gravs( 0 ) * std::cos ( currentFlightPathAngle) - gravs( 0 ) * std::sin( currentFlightPathAngle) * std::cos( currentHeading ) )
            -omega * omega * currentAltitude * std::cos( currentLatitude ) * ( std::cos( currentLatitude ) * std::cos( currentFlightPathAngle )  + std::sin( currentFlightPathAngle ) * std::sin( currentLatitude ) * std::cos( currentHeading ) )
            -( currentAirspeed * currentAirspeed / currentAltitude ) * std::cos( currentFlightPathAngle ) - 2 * currentMass * omega * currentAirspeed * std::cos( currentLatitude ) * std::sin( currentHeading );


}


double getBodyFlapDeflection( const double &del_C_m_b )
{
    double bodyflap = del_C_m_b;

    return bodyflap;
}

double computeHeatingRate (
        const double &airdensity,
        const double &airspeed,
        const double &C,
        const double &N,
        const double &M)
{
    return C * std::pow( airdensity, N ) * std::pow( airspeed, M );
}

double computeStagnationHeat (const double &airdensity,
                              const double &airspeed,
                              const double &C_s,
                              const double &N,
                              const double &M,
                              const double &adiabaticWallTemperature,
                              const double &WallTemperature)
{
    double C = C_s * ( 1 - ( WallTemperature / adiabaticWallTemperature ) ) ;
    double q_dot_s = computeHeatingRate ( airdensity, airspeed, C, N, M) ;
    //std::cout << "q_dot_s: " << q_dot_s << std::endl;

    return q_dot_s;
}

double computeStagnationHeatFlux (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems)
{
    const double airDensity = flightConditions->getCurrentDensity( );
    const double airSpeed = flightConditions->getCurrentAirspeed( );
    const double airTemperature = flightConditions->getCurrentFreestreamTemperature( );
    const double machNumber = flightConditions->getCurrentMachNumber( );
    const double noseRadius = vehicleSystems->getNoseRadius( );
    const double wallEmissivity = vehicleSystems->getWallEmissivity( );

    const double nose_term = std::pow( noseRadius, -0.5 );
    const double C_s = ( 1.83 ) * nose_term * 1E-4;
    const double M = 3.0;
    const double N = 0.5;

    // Compute adiabatic wall temperature.
    double adiabaticWallTemperature
            = tudat::aerodynamics::computeAdiabaticWallTemperature( airTemperature , machNumber );

    std::function< double( const double ) > heatTransferFunction = std::bind(
                &computeStagnationHeat, airDensity, airSpeed, C_s, N, M, adiabaticWallTemperature, std::placeholders::_1 );

    return tudat::aerodynamics::computeEquilibriumHeatflux( heatTransferFunction, wallEmissivity, adiabaticWallTemperature );
}

double computeFlatPlateHeat (const double &airdensity,
                             const double &airspeed,
                             const double &C_FP_1,
                             const double &C_FP_2,
                             const double &adiabaticWallTemperature,
                             const double &WallTemperature)
{
    double C;
    double M;
    double N;

    if ( airspeed <= 3962.0 )
    {
        M = 3.37;
        N = 0.8;
        C = C_FP_1 * std::pow( 556 / WallTemperature, 1.0 / 4.0 ) * ( 1.0 - 1.11 * ( WallTemperature / adiabaticWallTemperature ) );
    }
    if ( airspeed > 3962.0 )
    {
        M = 3.7;
        N = 0.8;
        C = C_FP_2 * ( 1 - 1.11 * ( WallTemperature / adiabaticWallTemperature ) );
    }

    return computeHeatingRate ( airdensity, airspeed, C, N, M);;
}

double computeFlatPlateHeatFlux (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        const std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
        const std::shared_ptr< tudat::system_models::BislipSystems > &bislipSystems)
{
    const double airDensity = flightConditions->getCurrentDensity( );
    const double airSpeed = flightConditions->getCurrentAirspeed( );
    const double airTemperature = flightConditions->getCurrentFreestreamTemperature( );
    const double machNumber = flightConditions->getCurrentMachNumber( );
    const double wallEmissivity = vehicleSystems->getWallEmissivity( );
    const double phi = bislipSystems->getLocalBodyAngle( );
    const double x_T = bislipSystems->getTransitionDistance( );

    const double sin_phi_term = std::pow( std::sin( phi ), 1.6 ) ;
    const double cos_phi = std::cos( phi );
    const double x_T_term = std::pow( x_T, -1.0 / 5.0 ) ;
    const double C_FP_1 = ( 3.35 ) * sin_phi_term * std::pow( cos_phi, 1.78 ) * x_T_term * 1E-4;
    const double C_FP_2 = ( 2.2 ) * sin_phi_term * std::pow( cos_phi, 2.08 ) * x_T_term * 1E-5;

    // Compute adiabatic wall temperature.
    double adiabaticWallTemperature
            = tudat::aerodynamics::computeAdiabaticWallTemperature( airTemperature , machNumber );

    std::function< double( const double ) > heatTransferFunction = std::bind(
                &computeFlatPlateHeat, airDensity, airSpeed, C_FP_1, C_FP_2, adiabaticWallTemperature, std::placeholders::_1 );

    return tudat::aerodynamics::computeEquilibriumHeatflux( heatTransferFunction, wallEmissivity, adiabaticWallTemperature );
}

double computeHeatingRateTauber (
        const double &q_dot_s,
        const double &q_dot_FP,
        const double &lambda)
{
    const double sin_lambda = std::sin( lambda );
    const double cos_lambda = std::cos( lambda );

    double q_dot_LE = std::pow ( 0.5 * q_dot_s * q_dot_s * cos_lambda * cos_lambda + q_dot_FP * q_dot_FP * sin_lambda * sin_lambda, 0.5 );
    //  std::cout << "q_dot_FP: " << q_dot_FP << std::endl;
    // std::cout << "q_dot_LE: " << q_dot_LE << std::endl;

    return q_dot_LE;
}

double computeBendingMoment (
        const std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions)
{
    const double dynamicPressure = flightConditions->getCurrentDynamicPressure();
    const double angleOfAttack = flightConditions->getAerodynamicAngleCalculator()->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );

    //  double q_dot_LE = std::pow ( 0.5 * q_dot_s * q_dot_s * cos_lambda * cos_lambda + q_dot_FP * q_dot_FP * sin_lambda * sin_lambda, 0.5 );
    //  std::cout << "q_dot_FP: " << q_dot_FP << std::endl;
    // std::cout << "q_dot_LE: " << q_dot_LE << std::endl;

    return dynamicPressure * angleOfAttack;
}

double computePenalty (
        const Eigen::VectorXd &dependentVariable_TimeHistory,
        const long &startIterator,
        const long &endIterator,
        const double &constraint, const double &fixedStepSize, const double &tof, const bool &direct )
{
    Eigen::VectorXd dependentVariable_Violation( dependentVariable_TimeHistory.size() );
    dependentVariable_Violation = Eigen::VectorXd::Zero( dependentVariable_TimeHistory.size() );

    if ( direct == true )
    {
        if ( startIterator == 1 )
        {  for ( long i = startIterator; i < endIterator + 1; i++ )
            {
                if ( dependentVariable_TimeHistory( i ) < dependentVariable_TimeHistory( i - 1 ) ) { dependentVariable_Violation( i ) = dependentVariable_TimeHistory( i - 1 ) - dependentVariable_TimeHistory( i ); }
                //   std::cout << "dependentVariable_Violation( " << i << " ): " << dependentVariable_Violation( i ) << std::endl;

            }
        }
        else
        {
            for ( long i = startIterator + 1; i < endIterator; i++ )
            {
                if ( dependentVariable_TimeHistory( i ) > dependentVariable_TimeHistory( i - 1 ) ) { dependentVariable_Violation( i ) = dependentVariable_TimeHistory( i ) - dependentVariable_TimeHistory( i - 1 ); }
                //   std::cout << "dependentVariable_Violation( " << i << " ): " << dependentVariable_Violation( i ) << std::endl;

            }
        }
    }
    else
    {
        for ( long i = startIterator; i < endIterator; i++ )
        {
            if ( dependentVariable_TimeHistory ( i ) > constraint ) { dependentVariable_Violation( i ) = dependentVariable_TimeHistory( i ) - constraint; }

            // std::cout << "dependentVariable_Violation( " << i << " ): " << dependentVariable_Violation( i ) << std::endl;
        }
    }
    double penalty = 0;
    std::ptrdiff_t index_MaximumViolation;
    double maximum_Violation = dependentVariable_Violation.maxCoeff( &index_MaximumViolation );

    if ( direct == true ) { penalty = dependentVariable_Violation.sum(); }
    else { penalty =  ( maximum_Violation / constraint ) + ( fixedStepSize * dependentVariable_Violation.sum() ) / ( tof * constraint ); }

    //std::cout << "penalty = " << penalty << std::endl;


    return penalty;
}



/*

Eigen::MatrixXd getDependentVariableMatrix( const tudat::propagators::SingleArcDynamicsSimulator< double > simulatedDynamics, const double simulationStartEpoch, const double fixedStepSize )
{

    const std::shared_ptr< tudat::propagators::SingleArcDynamicsSimulator< double, double > >& singleArcDynamicsSimulator,


    //! Extract map of dependent variables.
    const std::map< double, Eigen::VectorXd > dependentVariableMap = simulatedDynamics.SingleArcDynamicsSimulator::getDependentVariableHistory( );

    //! Declare number of rows (time) / columns (variables).
    unsigned long rows = dependentVariableMap.size();
    unsigned long columns = ( ( simulatedDynamics.getDependentVariableHistory( ).begin() )->second ).size();

    //! Declare and initialize dependent variable matrix.
    Eigen::MatrixXd dependentVariableMatrix( rows, columns );
    dependentVariableMatrix = Eigen::MatrixXd::Zero( rows, columns );

    //! Loop to populate the matrix with the extraction of the map.
    for ( unsigned long i = 0; i < rows; i++ )
    {
        dependentVariableMatrix.row( i ) = dependentVariableMap.at( simulationStartEpoch + i * fixedStepSize );

        std::cout << "dependentVariableMatrix.row( " << i << " ) : " << dependentVariableMatrix.row( i ) << std::endl;

    }

    return dependentVariableMatrix;
}

*/


} // namespace variables
} // namespace bislip
