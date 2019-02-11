#include <Tudat/Bislip/bislipVariables.h>

namespace bislip { namespace Variables {

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
    return std::chrono::duration_cast<std::chrono::milliseconds>( now - midnight ).count();
}

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
    return std::acos( std::sin( lat_c ) * std::sin( lat_f ) + std::cos( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );
}

double computeHeadingToTarget (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f)
{
    return std::atan2( std::sin( lon_f - lon_c ) * std::cos( lat_f ) , std::cos( lat_c ) * std::sin( lat_f ) - std::sin( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );
}

double computeHeadingError (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f,
        const double &heading)
{
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
std::vector< double > HermiteDerivatives(
        const Eigen::VectorXd &mappedNormalizedSpecificEnergy,
        const Eigen::VectorXd &y )
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
                HermiteDerivatives( mappedNormalizedSpecificEnergy, parameterValues ) );
}


bislip::Parameters::Optimization passOptimizationParameter (
        const bislip::Parameters::Optimization &parameter )
{
    return parameter;
}

std::string passString (
        const std::string &string )
{
    return string;
}

tudat::simulation_setup::NamedBodyMap& passBodyMap (
        tudat::simulation_setup::NamedBodyMap& bodyMap )
{
    return bodyMap;
}


std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > chooseGuidanceInterpolator(
        const bislip::Parameters::Optimization &parameter,
        const std::shared_ptr< bislip::VehicleSystems > &bislipSystems)
{
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > interpolator;

    if ( parameter == bislip::Parameters::Optimization::AngleOfAttack ) { interpolator = bislipSystems->getAngleOfAttackInterpolator(); }
    if ( parameter == bislip::Parameters::Optimization::BankAngle ) { interpolator = bislipSystems->getBankAngleInterpolator(); }
    if ( parameter == bislip::Parameters::Optimization::ThrustElevationAngle ) { interpolator = bislipSystems->getThrustElevationAngleInterpolator(); }
    if ( parameter == bislip::Parameters::Optimization::ThrustAzimuthAngle ) { interpolator = bislipSystems->getThrustAzimuthAngleInterpolator(); }
    if ( parameter == bislip::Parameters::Optimization::ThrottleSetting ) { interpolator = bislipSystems->getThrottleInterpolator(); }

    return interpolator;
}

std::pair < double, double > chooseGuidanceBounds (
        const bislip::Parameters::Optimization &parameter,
        const std::shared_ptr< bislip::VehicleSystems > &bislipSystems)
{
    std::pair < double, double > bounds;

    if ( parameter == bislip::Parameters::Optimization::AngleOfAttack ) { bounds = bislipSystems->getParameterBounds( parameter ); }
    if ( parameter == bislip::Parameters::Optimization::BankAngle ) { bounds = bislipSystems->getParameterBounds( parameter ); }
    if ( parameter == bislip::Parameters::Optimization::ThrustElevationAngle ) { bounds = bislipSystems->getParameterBounds( parameter ); }
    if ( parameter == bislip::Parameters::Optimization::ThrustAzimuthAngle ) {  bounds = bislipSystems->getParameterBounds( parameter ); }
    if ( parameter == bislip::Parameters::Optimization::ThrottleSetting ) { bounds = bislipSystems->getParameterBounds( parameter ); }

    return bounds;
}

double evaluateGuidanceInterpolator (
        const bislip::Parameters::Optimization &parameter,
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions =
            std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap.at( vehicleName )->getFlightConditions( ) );
    std::shared_ptr< bislip::VehicleSystems > bislipSystems =
            std::dynamic_pointer_cast< bislip::VehicleSystems >( bodyMap.at( vehicleName )->getBislipSystems( ) );

    //! Evaluate interpolator.
    double evaluation = ( bislip::Variables::chooseGuidanceInterpolator( parameter, bislipSystems ) )->interpolate( bislip::Variables::computeNormalizedSpecificEnergy( flightConditions->getCurrentAltitude(), flightConditions->getCurrentAirspeed(), bislipSystems->getE_max() ) );

    //! Select parameter bounds.
    std::pair < double, double > bounds = bislip::Variables::chooseGuidanceBounds( parameter, bislipSystems );

    //! Impose bounds.
    if ( evaluation < bounds.first ){ evaluation = bounds.first; }
    if ( evaluation > bounds.second ){ evaluation = bounds.second; }

    return evaluation;
}

Eigen::Vector6d computeCurrentCoefficients (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions =
            std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap.at( vehicleName )->getFlightConditions( ) );

    std::shared_ptr< bislip::VehicleSystems > bislipSystems =
            std::dynamic_pointer_cast< bislip::VehicleSystems >( bodyMap.at( vehicleName )->getBislipSystems( ) );

    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface =
            std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >( bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    double angleOfAttack = tudat::unit_conversions::convertDegreesToRadians(
                evaluateGuidanceInterpolator (
                    bislip::Parameters::Optimization::AngleOfAttack,
                    bodyMap,
                    vehicleName ) );

    //! Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
    std::vector< double > aerodynamicCoefficientInput;
    aerodynamicCoefficientInput.push_back( angleOfAttack );
    aerodynamicCoefficientInput.push_back( flightConditions->getCurrentMachNumber() );

    // Define values of independent variables of control surface aerodynamics
    std::map< std::string, std::vector< double > > controlSurfaceCoefficientInput;
    controlSurfaceCoefficientInput[ "BodyFlap" ] = aerodynamicCoefficientInput;
    controlSurfaceCoefficientInput[ "BodyFlap" ].push_back( 0.0 );



    // Update and retrieve current aerodynamic coefficients
    coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );

    return coefficientInterface->getCurrentAerodynamicCoefficients( );
}

Eigen::Vector3d computeBodyFixedThrustDirection (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap.at( vehicleName )->getFlightConditions( ) );
    std::shared_ptr< bislip::VehicleSystems > bislipSystems = std::dynamic_pointer_cast< bislip::VehicleSystems >( bodyMap.at( vehicleName )->getBislipSystems( ) );

    // std::cout << "Computing Body Fixed Thrust Direction" << std::endl;
    double eps_T = tudat::unit_conversions::convertDegreesToRadians(
                evaluateGuidanceInterpolator (
                    bislip::Parameters::Optimization::ThrustElevationAngle,
                    bodyMap,
                    vehicleName ) );

    double phi_T = tudat::unit_conversions::convertDegreesToRadians(
                evaluateGuidanceInterpolator (
                    bislip::Parameters::Optimization::ThrustAzimuthAngle,
                    bodyMap,
                    vehicleName ) );

    Eigen::Vector3d bodyFixedThrustDirection;
    bodyFixedThrustDirection( 0 ) = std::cos( eps_T ) * std::cos( phi_T );
    bodyFixedThrustDirection( 1 ) = std::cos( eps_T ) * std::sin( phi_T );
    bodyFixedThrustDirection( 2 ) = std::sin( eps_T );

    return bodyFixedThrustDirection;
}

double computeThrustMagnitude (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap.at( vehicleName )->getFlightConditions( ) );
    std::shared_ptr< bislip::VehicleSystems > bislipSystems = std::dynamic_pointer_cast< bislip::VehicleSystems >( bodyMap.at( vehicleName )->getBislipSystems( ) );

    const double throttle = evaluateGuidanceInterpolator (
                bislip::Parameters::Optimization::ThrottleSetting,
                bodyMap,
                vehicleName );

    return throttle * bislipSystems->getMaxThrust();
}

Eigen::Vector3d computeBodyFixedThrustVector (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >( bodyMap.at( vehicleName )->getFlightConditions( ) );
    std::shared_ptr< bislip::VehicleSystems > bislipSystems = std::dynamic_pointer_cast< bislip::VehicleSystems >( bodyMap.at( vehicleName )->getBislipSystems( ) );

    const Eigen::Vector3d BodyFixedThrustDirection = computeBodyFixedThrustDirection( bodyMap, vehicleName );
    const double thrustMagnitude = computeThrustMagnitude( bodyMap, vehicleName );

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

Eigen::Vector2d computeLocalGravity (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName )
{   
    const double currentAltitude = bodyMap.at( vehicleName )->getFlightConditions( )->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );
    const double currentLatitude = bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double mu = bodyMap.at( centralBodyName )->getGravityFieldModel( )->getGravitationalParameter();
    const double R_E = bodyMap.at( centralBodyName )->getShapeModel( )->getAverageRadius();
    const double J2 = 1082.626523e-6;
    const double J3 = 2.532153e-7;
    const double J4 = 1.6109876e-7;


    const double s_latitude = std::sin( currentLatitude );
    const double s_latitude_sq = s_latitude * s_latitude;
    const double c_latitude = std::cos( currentLatitude );

    const double g_n_pre = -3 * ( mu / ( currentAltitude * currentAltitude ) ) * pow( R_E / currentAltitude , 2 ) * s_latitude * c_latitude;
    const double g_n_sub_1 = J2;
    const double g_n_sub_2 = ( J3 / 2 ) * ( R_E / currentAltitude ) * ( 1 / s_latitude ) * ( 5 * s_latitude_sq - 1 );
    const double g_n_sub_3 = ( 5 * J4 / 6 ) * pow( R_E / currentAltitude , 2) * ( 7 * s_latitude_sq - 3 );
    const double g_n = g_n_pre * ( g_n_sub_1 + g_n_sub_2 + g_n_sub_3 );

    const double g_d_pre = ( mu / ( currentAltitude * currentAltitude ) );
    const double g_d_sub_1 = -( 3 * J2 / 2 ) * pow( R_E / currentAltitude , 2 ) * ( 3 * s_latitude_sq - 1 );
    const double g_d_sub_2 = -2 * J3  * pow( R_E / currentAltitude , 3 ) * ( s_latitude ) * ( 5 * s_latitude_sq - 3 );
    const double g_d_sub_3 = -( 5 * J4 / 8 ) * pow( R_E / currentAltitude , 4 ) * ( 35 * s_latitude_sq * s_latitude_sq - 30 * s_latitude_sq + 3 );
    const double g_d = g_d_pre * ( 1 + g_d_sub_1 + g_d_sub_2 + g_d_sub_3 );

    Eigen::Vector2d localGravity ( 2 );
    localGravity << g_n, g_d;

    return localGravity;
}

double computeEquilibriumGlideLimit (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName )
{
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    std::shared_ptr< bislip::VehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    const double currentFlightPathAngle = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
    const double currentLatitude = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double currentHeading = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
    const double currentAltitude = flightConditions->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );
    const double currentAirspeed = flightConditions->getCurrentAirspeed( );
    const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    const double currentMass = bodyMap.at( vehicleName )->getBodyMass();

    const double omega = 7.292115*1E-5;

    Eigen::Vector6d currentCoefficients = computeCurrentCoefficients ( bodyMap, vehicleName );
    const double currentLift = currentDynamicPressure * ( bislipSystems->getReferenceArea( ) ) * currentCoefficients[ 2 ];

    Eigen::Vector2d gravs = bislip::Variables::computeLocalGravity( bodyMap, vehicleName, centralBodyName );

    const double currentThrustMagnitude = bislip::Variables::computeThrustMagnitude( bodyMap, vehicleName );
    const double thrustElevationAngle = tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustElevationAngle, bodyMap, vehicleName ) );
    const double thrustAzimuthAngle = tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustAzimuthAngle, bodyMap, vehicleName ) );

    const double s_delta = std::sin( currentLatitude );
    const double c_delta = std::cos( currentLatitude );
    const double s_gamma = std::cos ( currentFlightPathAngle );
    const double c_gamma = std::cos ( currentFlightPathAngle );
    const double s_chi = std::sin( currentHeading );
    const double c_chi = std::cos( currentHeading );


    const double a = currentLift + currentThrustMagnitude * ( s_gamma * std::cos( thrustAzimuthAngle ) * std::cos( thrustElevationAngle ) + c_gamma * std::sin( thrustElevationAngle ) );
    const double b = -currentThrustMagnitude * std::sin( thrustAzimuthAngle ) * std::cos( thrustElevationAngle );
    const double c = std::sqrt( a * a + b * b );
    const double phi = std::atan2( b, a );
    const double term1 = omega * omega * currentAltitude * c_delta * ( c_delta * c_gamma + s_gamma * s_delta * c_chi );
    const double term2 = ( currentAirspeed * currentAirspeed / currentAltitude ) * c_gamma;
    const double term3 = 2.0 * omega * currentAirspeed * c_delta * s_chi;
    const double term4 = gravs( 1 ) * c_gamma;
    const double term5 = gravs( 0 ) * s_gamma * c_chi;
    //const double argument =  phi - ( currentMass / c ) * ( term1 + term2 + term3 - term4 + term5 );
    const double limit = tudat::unit_conversions::convertRadiansToDegrees( std::asin( ( phi - ( currentMass / c ) * ( term1 + term2 + term3 - term4 + term5 ) ) ) );
    /*
    std::cout << "CL: " << currentCoefficients[ 2 ] << std::endl;
    std::cout << "currentMass: " << currentMass << std::endl;
    std::cout << "currentAirspeed: " << currentAirspeed << std::endl;
    std::cout << "currentLift: " << currentLift << std::endl;
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "c: " << c << std::endl;
    std::cout << "phi: " << phi << std::endl;
    std::cout << "term1: " << term1 << std::endl;
    std::cout << "term2: " << term2 << std::endl;
    std::cout << "term3: " << term3 << std::endl;
    std::cout << "term4: " << term4 << std::endl;
    std::cout << "term5: " << term5 << std::endl;
    std::cout << "argument: " << argument << std::endl;
    std::cout << "limit: " << limit << std::endl;
*/

    return limit;
}

double computeBendingMoment (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    std::shared_ptr< bislip::VehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    const double dynamicPressure = flightConditions->getCurrentDynamicPressure();
    const double angleOfAttack = flightConditions->getAerodynamicAngleCalculator()->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );

    //  double q_dot_LE = std::pow ( 0.5 * q_dot_s * q_dot_s * cos_lambda * cos_lambda + q_dot_FP * q_dot_FP * sin_lambda * sin_lambda, 0.5 );
    //  std::cout << "q_dot_FP: " << q_dot_FP << std::endl;
    // std::cout << "q_dot_LE: " << q_dot_LE << std::endl;

    return dynamicPressure * angleOfAttack;
}

double computeBodyFlapDeflection(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    std::shared_ptr< bislip::VehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    //std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > controlSurfaceCoefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
    //          bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );


    double angleOfAttack = tudat::unit_conversions::convertDegreesToRadians(
                evaluateGuidanceInterpolator (
                    bislip::Parameters::Optimization::AngleOfAttack,
                    bodyMap,
                    vehicleName ) );

    // double bodyFlapCmIncrement = 0;
    // double CmIncrement = 0;

    //! Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
    std::vector< double > aerodynamicCoefficientInput;
    aerodynamicCoefficientInput.push_back( angleOfAttack );
    aerodynamicCoefficientInput.push_back( flightConditions->getCurrentMachNumber() );
    //aerodynamicCoefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput );

    // Define values of independent variables of control surface aerodynamics
    std::map< std::string, std::vector< double > > controlSurfaceCoefficientInput;
    controlSurfaceCoefficientInput[ "BodyFlap" ] = aerodynamicCoefficientInput;
    controlSurfaceCoefficientInput[ "BodyFlap" ].push_back( tudat::unit_conversions::convertDegreesToRadians( 0.0 ) );

    // Eigen::Vector3d momentWithIncrement, momentWithoutIncrement;
    //Eigen::Vector6d currentCoefficients;
    // momentWithIncrement = aerodynamicCoefficientInterface->getCurrentMomentCoefficients( );



    // Update and retrieve current aerodynamic coefficients
    coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
    // momentWithoutIncrement = controlSurfaceCoefficientInterface->getCurrentMomentCoefficients( );
    // bodyFlapCmIncrement = bislip::Variables::computeBodyFlapCmIncrement( bodyMap, vehicleName, controlSurfaceCoefficientInterface->getCurrentAerodynamicCoefficients() );
    // CmIncrement = momentWithoutIncrement( 1 ) - momentWithIncrement( 1 );
    //controlSurfaceCoefficientInput[ "BodyFlap" ][ 2 ] += 0.001;

    // Determine function for which the root is to be determined.
    std::function< double( const double ) > coefficientFunction =
            std::bind( &bislip::Variables::computeFullPitchMomentCoefficient,
                       coefficientInterface, std::placeholders::_1, aerodynamicCoefficientInput, controlSurfaceCoefficientInput );

    double bodyFlapDeflection = TUDAT_NAN;

    //! Object to iteratively find the root of the equations C_m(alpha)=0, i.e. to determine the
    //!  angle of attack for which the pitch moment is zero.
    std::shared_ptr< tudat::root_finders::RootFinderCore< double > > rootfinder = std::make_shared< tudat::root_finders::SecantRootFinderCore< double > >(
                std::bind(
                    &tudat::root_finders::termination_conditions::RootAbsoluteToleranceTerminationCondition< double >::checkTerminationCondition,
                    std::make_shared< tudat::root_finders::termination_conditions::RootAbsoluteToleranceTerminationCondition
                    < double > >( 1.0E-15, 1000 ), std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 ), 0.5 );


    // Find root of pitch moment function
    try
    {
        bodyFlapDeflection = rootfinder->execute(
                    std::make_shared< tudat::basic_mathematics::FunctionProxy< double, double > >( coefficientFunction ),
                    controlSurfaceCoefficientInput[ "BodyFlap" ][ 2 ] );
    }
    // Throw error if not converged
    catch( std::runtime_error )
    {
        throw std::runtime_error( "Error when bodyflap deflection, root finder did not converge." );

    }

    return bodyFlapDeflection;

}


double computeFullPitchMomentCoefficient(
        const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
        const double &bodyFlapDeflection,
        const std::vector< double > &aerodynamicCoefficientInput,
        const std::map< std::string, std::vector< double > > &controlSurfaceCoefficientInput )
{
    // Update coefficients to perturbed independent variables
    std::map< std::string, std::vector< double > > perturbedControlSurfaceConditions =
            controlSurfaceCoefficientInput;


    perturbedControlSurfaceConditions[ "BodyFlap" ][ 2 ] = bodyFlapDeflection;

    coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );

    return coefficientInterface->getCurrentMomentCoefficients( )( 1 );


}









double computeBodyFlapCmIncrement (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const Eigen::Vector6d &currentCoefficients )
{
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
    //            bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    std::shared_ptr< bislip::VehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    const Eigen::Vector3d referenceValues = bislipSystems->getReferenceValues();
    const Eigen::Vector3d aerodynamicReferenceCenter = bislipSystems->getMassReferenceCenter();
    const Eigen::Vector3d thrustReferenceCenter = bislipSystems->getThrustReferenceCenter();

    //const Eigen::Vector6d currentCoefficients = bislip::Variables::computeCurrentCoefficients( bodyMap, vehicleName );
    const double currentThrustMagnitude = bislip::Variables::computeThrustMagnitude( bodyMap, vehicleName );
    const double currentDrag = currentDynamicPressure * referenceValues[ 0 ] * currentCoefficients[ 0 ];
    const double currentLift = currentDynamicPressure * referenceValues[ 0 ] * currentCoefficients[ 2 ];
    const double surfaceArea = referenceValues[ 0 ];
    // const double b_ref = referenceValues[ 1 ];
    const double c_ref = referenceValues[ 2 ];
    const double del_x = aerodynamicReferenceCenter[ 0 ];
    const double del_x_T = thrustReferenceCenter[ 0 ];
    const double del_z_T = thrustReferenceCenter[ 2 ];

    double angleOfAttack = tudat::unit_conversions::convertDegreesToRadians(
                bislip::Variables::evaluateGuidanceInterpolator (
                    bislip::Parameters::Optimization::AngleOfAttack,
                    bodyMap,
                    vehicleName ) );

    double thrustElevationAngle = tudat::unit_conversions::convertDegreesToRadians(
                bislip::Variables::evaluateGuidanceInterpolator (
                    bislip::Parameters::Optimization::ThrustElevationAngle,
                    bodyMap,
                    vehicleName ) );

    double thrustAzimuthAngle = tudat::unit_conversions::convertDegreesToRadians(
                bislip::Variables::evaluateGuidanceInterpolator (
                    bislip::Parameters::Optimization::ThrustAzimuthAngle,
                    bodyMap,
                    vehicleName ) );

    double term1 = ( del_x / currentDynamicPressure * surfaceArea * c_ref ) * ( currentDrag * std::sin( angleOfAttack ) - currentLift * std::cos( angleOfAttack ) );
    double term2 = ( currentThrustMagnitude / ( currentDynamicPressure * surfaceArea * c_ref ) ) * ( del_x_T * std::sin( thrustElevationAngle ) - del_z_T * std::cos( thrustElevationAngle ) * std::cos( thrustAzimuthAngle ) );

    return ( term1 - term2 - currentCoefficients[ 4 ] );
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

double computeStagnationHeat (
        const double &airdensity,
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

double computeFlatPlateHeat (
        const double &airdensity,
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
        const std::shared_ptr< bislip::VehicleSystems > &bislipSystems)
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

double computePenalty (
        const Eigen::VectorXd &dependentVariable_TimeHistory,
        const long &startIterator,
        const long &endIterator,
        const double &constraint,
        const double &fixedStepSize,
        const double &tof,
        const bool &direct )
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
bool StopOrNot( const tudat::simulation_setup::NamedBodyMap& bodyMap,
                const std::string &vehicleName,
                const std::vector< double > &vehicleParameterValues,
                const std::vector< double > &terminationConditionsValues) {

    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > FlightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract current latitude
    const double lat_c_rad =  FlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );

    //! Extract current longitude
    const double lon_c_rad =  FlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );

    //! Extract initial coordinates
    const double lat_i_rad = bodyMap.at( vehicleName )->getBislipSystems()->getInitialCoordinates().first;
    const double lon_i_rad = bodyMap.at( vehicleName )->getBislipSystems()->getInitialCoordinates().second;

    //! Extract Final distance to target
    const double final_d_to_target_rad = bodyMap.at( vehicleName )->getBislipSystems()->getFinalDistanceToTarget();
    const double h_UP = terminationConditionsValues[ 3 ];
    const double h_DN = terminationConditionsValues[ 4 ];
    const double V_max = terminationConditionsValues[ 5 ];
    const double n_max = terminationConditionsValues[ 6 ];
    const double q_dot_max = terminationConditionsValues[ 7 ];
    const double q_dyn_max = terminationConditionsValues[ 8 ];

    //! Extract and convert target coordinates
    const double lat_f_rad = tudat::unit_conversions::convertDegreesToRadians( bodyMap.at( vehicleName )->getBislipSystems()->getTargetCoordinates().first );
    const double lon_f_rad = tudat::unit_conversions::convertDegreesToRadians( bodyMap.at( vehicleName )->getBislipSystems()->getTargetCoordinates().second );

    //! Calculate current Distance to target.
    const double current_d_to_target_rad = bislip::Variables::computeAngularDistance( lat_c_rad , lon_c_rad , lat_f_rad , lon_f_rad );

    //! Calculate total distance traveled
    const double total_d_traveled_rad = bislip::Variables::computeAngularDistance( lat_i_rad , lon_i_rad , lat_c_rad , lon_c_rad );

    //! Extract Initial distance to target
    const double initial_d_to_target_rad = bodyMap.at( vehicleName )->getBislipSystems()->getInitialDistanceToTarget();

    //! Extract current altitude
    const double current_height = FlightConditions->getCurrentAltitude();

    //! Extract current airspeed
    const double current_V = FlightConditions->getCurrentAirspeed();

    //! Extract current angle of attack
    const double current_AoA = FlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );

    //! Extract current flight-path angle
    const double current_gamma = FlightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );

    //! Exract current mass
    const double currentMass = bodyMap.at( vehicleName )->getBodyMass( );

    //! Calculate current normalized specific energy.
    const double E_hat = bislip::Variables::computeNormalizedSpecificEnergy( current_height, current_V, bodyMap.at( vehicleName )->getBislipSystems()->getE_max() );

    //! Evaluate current throttle setting and thrust elevation angle.
    double throttle = bislip::Variables::evaluateGuidanceInterpolator(
                    //current_gamma,
                    bislip::Parameters::Optimization::ThrottleSetting,
                    bodyMap.at( vehicleName )->getBislipSystems(),
                    current_height,
                    current_V,
                    bodyMap.at( vehicleName )->getBislipSystems()->getE_max() );

    const double finalMass =  bodyMap.at( vehicleName )->getVehicleSystems()->getDryMass();
    const double current_rho = FlightConditions->getCurrentDensity( );
    const double current_q_dot = FlightConditions->getCurrentAerodynamicHeatRate( );
    const double current_q_dyn = FlightConditions->getCurrentDynamicPressure( );

    //! Declare and initialize boolean that is returned to terminate the simulation.
    bool done = false;

    //! Currently 3 conditions could terminate a simulation
    //!     When distance to target is less than a specified number
    //!     When altitude is more than a specified number
    //!     When altitude is less than a specified number
    //!     When angular distance traveled is larger than angular distance among endpoints

*/
/* if ( ( current_d_to_target_rad <= final_d_to_target_rad  ) || ( current_height <= h_DN ) || ( current_height >= h_UP ) || ( total_d_traveled_rad >= initial_d_to_target_rad ) )
    {
        done = true;
        std::cout << " first " << std::endl;
        std::cout << "current_d_to_target_rad:  " << current_d_to_target_rad << std::endl;
        std::cout << "final_d_to_target_rad :  " << final_d_to_target_rad  << std::endl;
        std::cout << "current_height:  " << current_height << std::endl;
        std::cout << "h_DN :  " << h_DN  << std::endl;
        std::cout << "h_UP :  " << h_UP  << std::endl;
        std::cout << "total_d_traveled_rad:  " << total_d_traveled_rad << std::endl;
        std::cout << "initial_d_to_target_rad:  " << initial_d_to_target_rad << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
   else if ( current_q_dot > q_dot_max*1 )
    {
        done = true;
        std::cout << " second " << std::endl;
        std::cout << "elapsed time:  " << FlightConditions->getCurrentTime() - bodyMap.at( vehicleName )->getVehicleSystems()->getStartingEpoch() << std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "Expected thrust:  " << throttle * bodyMap.at( vehicleName )->getVehicleSystems()->getMaxThrust()<< std::endl;
        std::cout << "E_hat:  " << E_hat << std::endl;
        std::cout << "current_V:  " << current_V << std::endl;
        std::cout << "current_rho:  " << current_rho << std::endl;
        std::cout << "current_q_dot:  " << current_q_dot << std::endl;
        std::cout << "q_dot_max:  " << q_dot_max << std::endl;
        std::cout << "currentMass:  " << currentMass << std::endl;
        std::cout << "finalMass:  " << finalMass << std::endl;
        std::cout << "Stopping propagation" << std::endl;

    }
    else if ( current_q_dyn > q_dyn_max*1 )
    {
        done = true;
        std::cout << " third " << std::endl;
        std::cout << "elapsed time:  " << FlightConditions->getCurrentTime() - bodyMap.at( vehicleName )->getVehicleSystems()->getStartingEpoch() << std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "E_hat:  " << E_hat << std::endl;
        std::cout << "Expected thrust:  " << throttle * bodyMap.at( vehicleName )->getVehicleSystems()->getMaxThrust() << std::endl;
        std::cout << "current_V:  " << current_V << std::endl;
        std::cout << "current_rho:  " << current_rho << std::endl;
        std::cout << "current_q_dyn:  " << current_q_dyn << std::endl;
        std::cout << "q_dyn_max:  " << q_dyn_max << std::endl;
        std::cout << "currentMass:  " << currentMass << std::endl;
        std::cout << "finalMass:  " << finalMass << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
    else if ( currentMass < finalMass )
    {
        done = true;
        std::cout << " fourth " << std::endl;
        std::cout << "elapsed time:  " << FlightConditions->getCurrentTime() - bodyMap.at( vehicleName )->getVehicleSystems()->getStartingEpoch()<< std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "Expected thrust:  " << throttle * bodyMap.at( vehicleName )->getVehicleSystems()->getMaxThrust() << std::endl;
        std::cout << "E_hat:  " << E_hat << std::endl;
        std::cout << "current_q_dot:  " << current_q_dot << std::endl;
        std::cout << "q_dot_max:  " << q_dot_max << std::endl;
        std::cout << "currentMass:  " << currentMass << std::endl;
        std::cout << "finalMass:  " << finalMass << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
    else if ( current_gamma < 0 )
    {
        done = true;
        std::cout << " fifth " << std::endl;
        std::cout << "elapsed time:  " << FlightConditions->getCurrentTime() - additional_data[ 2 ]<< std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "Expected thrust:  " << throttle * additional_data[ 3 ]<< std::endl;
        std::cout << "E_hat:  " << E_hat << std::endl;
        std::cout << "current_gamma:  " << current_gamma << std::endl;
        std::cout << "current_q_dot:  " << current_q_dot << std::endl;
        std::cout << "q_dot_max:  " << q_dot_max << std::endl;
        std::cout << "currentMass:  " << currentMass << std::endl;
        std::cout << "finalMass:  " << finalMass << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
*/
/*    if ( ( ( throttle >= 0 ) && ( throttle <= 1.0 ) )
              &&  ( ( current_gamma > tudat::mathematical_constants::PI / 2.0 ) ) )
    {
        //! This condition terminates the simulation if flight-path angle is negative or beyond 90 deg while the engine is on.
        done = true;
        std::cout << " sixth " << std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "current_gamma:  " << tudat::unit_conversions::convertRadiansToDegrees( current_gamma ) << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
    else if ( E_hat > 1.0 )
    {
        done = true;
        std::cout << " seventh " << std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "E_hat:  " << E_hat << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }

    // else if ( ( current_AoA < tudat::unit_conversions::convertDegreesToRadians( parameterBounds_Ascent[ 2 ] ) ) || ( current_AoA > tudat::unit_conversions::convertDegreesToRadians( parameterBounds_Ascent[ 3 ] ) ) )
       // {
            //! This condition terminates the simulation if angle of attack is beyond the specified bounds while the engine is on.
       //     done = true;
       //     std::cout << "throttle:  " << throttle << std::endl;
            std::cout << "current_AoA:  " << tudat::unit_conversions::convertRadiansToDegrees( current_AoA ) << std::endl;
        //    std::cout << "Stopping propagation" << std::endl;
        //}
        //else if ( ( eps_T_deg < parameterBounds_Ascent[ 4 ] ) || ( eps_T_deg > parameterBounds_Ascent[ 5 ] ) )
       // {
            //! This condition terminates the simulation if thrust elevation angle is beyond the specified bounds while the engine is on.
        //    done = true;
        //    std::cout << "throttle:  " << throttle << std::endl;
            std::cout << "eps_T_deg:  " << eps_T_deg << std::endl;
            std::cout << "Stopping propagation" << std::endl;
        }

    //}
*/
/*    else
    {
        done = false;
    }




    return done;
}

*/
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
/*
pagmo::algorithm getPagmoAlgorithm( const int index )
{
    switch( index )
    {
    case 0:
    {
        pagmo::algorithm algo{ pagmo::nsga2( ) };
        return algo;
        break;
    }
    case 1:
    {
        pagmo::algorithm algo{ pagmo::moead( ) };
        return algo;
        break;
    }
    case 2:
    {
        pagmo::algorithm algo{ pagmo::ihs( ) };
        return algo;
        break;
    }
    default:
    {
        throw std::runtime_error( "Error, multi-objective pagmo algorithm " + std::to_string( index ) + " was not found." );
    }
    }
}

*/
} //namespace Variables
                 } // namespace bislip
