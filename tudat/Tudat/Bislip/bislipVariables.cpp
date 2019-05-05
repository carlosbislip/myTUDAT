#include <Tudat/Bislip/bislipVariables.h>

namespace bislip { namespace Variables {


//! Get path for output directory.
std::string getOutputPath(
        const std::string& extraDirectory )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    //    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
    //                                                std::string( "tudatBundle/tudatApplications/Space4ErrBody_Executables_testing/Space4ErrBody_Executables_and_Headers_testing/applicationOutput_tudat.h" ).length( ) );
    //    std::string outputPath = reducedPath + "Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/code/SimulationOutput/";
    std::string outputPath = "/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/tudatBundle.git/tudatApplications/Space4ErrBody.git/Space4ErrBody_0/matlab/SimulationOutput/";
    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

Eigen::MatrixXd convertVectorOfVectorsDoubleToEigenMatrixXd( const std::vector< std::vector< double > > &vectorOfVectorsDouble )
{
    Eigen::MatrixXd matrix( vectorOfVectorsDouble.size( ), vectorOfVectorsDouble.at( 0 ).size( ) );
    for( unsigned int i = 0; i < vectorOfVectorsDouble.size( ); i++ )
    {
        for( unsigned int j = 0; j < vectorOfVectorsDouble.at( 0 ).size( ); j++ )
        {
            matrix( i, j ) = vectorOfVectorsDouble.at( i ).at( j );
        }
    }

    return matrix;
}


void printEigenMatrixXdToFile( const Eigen::MatrixXd & matrixToPrint,
                               const std::string &filename,
                               const std::string &outputSubFolder )
{

    /* Eigen::MatrixXd matrixToPrint( population.size( ), population.at( 0 ).size( ) );
    for( unsigned int i = 0; i < population.size( ); i++ )
    {
        for( unsigned int j = 0; j < population.at( 0 ).size( ); j++ )
        {
            matrixToPrint( i, j ) = population.at( i ).at( j );
        }
    }

    */
    //if( !isFitness )
    //{
    tudat::input_output::writeMatrixToFile( matrixToPrint,
                                            filename,
                                            16,
                                            bislip::Variables::getOutputPath( ) + outputSubFolder);
    //}
    //else
    /* {
        tudat::input_output::writeMatrixToFile( matrixToPrint,
                                                "fitness_" + fileSuffix + ".dat",
                                                16,
                                                bislip::Variables::getOutputPath( ) + outputSubFolder);
    } */
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

double rootFinderBisection(
        const std::function< double( const double ) > &function,
        const double &minimum,
        const double &maximum )
{

    double a = minimum;
    double b = maximum;
    double root = 0.0;
    double f_a, f_root;
    if ( function( a ) * function( b ) < 0.0 )
    {
        root = ( a + b ) / 2.0;

        // std::cout << "     Starting Bisection search:" << std::endl;
        while ( std::abs( function( root ) ) > 1e-7 )
        {

            f_a       = function( a );
            f_root    = function( root );

            if ( f_a * f_root < 0.0 ) { b = root; }
            else { a = root; }
            root = ( a + b ) / 2.0;
            //   std::cout << "          root = " << root << "     |     " << "function( root )  = " << std::abs( function( root ) ) << std::endl;
        }
    }

    return root;
}

double computeAngularDistance (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f)
{
    //! Compute the Angular Distance between two points. Based on spherical great circles.
    double angularDistance = std::acos( std::sin( lat_c ) * std::sin( lat_f ) + std::cos( lat_c ) * std::cos( lat_f ) * std::cos( lon_c - lon_f ) );

    //! Determine if calculated distance is not a number.
    if( std::isnan( angularDistance) == true ) { angularDistance = 0; }

    return angularDistance;
}

double computeHeadingToTarget (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f)
{
    //! Declare and calculate the Current Heading Angle to Target.
    double headingToTarget = std::atan2( std::sin( lon_f - lon_c ) * std::cos( lat_f ) , std::sin( lat_f ) * std::cos( lat_c ) - std::sin( lat_c ) * std::cos( lat_f ) * std::cos( lon_f - lon_c ) );

    //! Convert to a positive angle if required.
    if ( headingToTarget < 0.0 ) { headingToTarget += 2 * tudat::mathematical_constants::PI; }

    return headingToTarget;
}
double computeHeadingError (
        const double &lat_c,
        const double &lon_c,
        const double &lat_f,
        const double &lon_f,
        const double &currentHeading)
{
    //!  Convert the Current Heading Angle to a positive angle if required.
    double positiveCurrentHeading = currentHeading;
    if ( currentHeading < 0.0 ) { positiveCurrentHeading += 2 * tudat::mathematical_constants::PI; }

    //! Declare and determine the Current Heading Angle to Target. Convert to a positive angle if required.
    double positiveHeadingToTarget = bislip::Variables::computeHeadingToTarget( lat_c, lon_c, lat_f, lon_f );
    if ( positiveHeadingToTarget < 0.0 ) { positiveHeadingToTarget += 2 * tudat::mathematical_constants::PI; }

    return positiveCurrentHeading - positiveHeadingToTarget;
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
        const Eigen::VectorXd &xValues,
        const Eigen::VectorXd &yValues )
{
    long dataPoints = xValues.size();
    Eigen::VectorXd h( dataPoints - 1 ), dely( dataPoints - 1 ), b( dataPoints ), x( dataPoints );
    Eigen::MatrixXd A( dataPoints, dataPoints );
    std::vector< double > x_vect;
    double mu, lambda;

    for ( long i = 0; i < dataPoints - 1; ++i ) { h( i ) = xValues( i + 1 ) - xValues( i ); dely( i ) = ( yValues( i + 1 ) - yValues( i ) ) / h( i ); }

    A = Eigen::MatrixXd::Zero( dataPoints, dataPoints );
    A( 0, 0 ) = 4.0;
    A( 0, 1 ) = -1.0;
    A( dataPoints - 1, dataPoints - 2 ) = -1.0;
    A( dataPoints - 1, dataPoints - 1 ) = 4.0;

    for ( long i = 1; i < dataPoints - 1; ++i )
    {
        lambda = h( i ) / ( h( i - 1 ) + h( i ) );
        mu = 1 - lambda;
        A( i , i - 1 ) = -mu;
        A( i , i ) = 4.0;
        A( i , i + 1 ) = -lambda * mu;
    }

    b( 0 ) = dely( 0 );
    b( dataPoints - 1 ) = dely( dataPoints - 2 );

    for ( long i = 1; i < dataPoints - 1; ++i ) { b( i ) = 3 * ( yValues( i + 1 ) - yValues( i - 1 ) ) / ( h( i - 1 ) + h( i ) ); }

    x = A.fullPivHouseholderQr().solve( b );

    for ( long i = 0; i < dataPoints; ++i ) { x_vect.push_back( x( i ) ); }

    //std::cout << "dataPoints: " << dataPoints << std::endl;
    //std::cout << "x_vect size: " << x_vect.size() << std::endl;
    //std::cout << "x size: " << x.size() << std::endl;
    return x_vect;
}

std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > createOneDimensionalHermiteInterpolator (
        const Eigen::VectorXd &yValues,
        const Eigen::VectorXd &xValues,
        const std::map< double, double > &mapped_data,
        const std::shared_ptr< tudat::interpolators::InterpolatorSettings > &interpolatorSettings )
{

    // std::cout << "yValues size: " << yValues.size()  << std::endl;
    // std::cout << "xValues size: " << xValues.size()  << std::endl;
    // std::cout << "derivatives size: " <<( bislip::Variables::HermiteDerivatives( xValues, yValues ) ).size()  << std::endl;

    return tudat::interpolators::createOneDimensionalInterpolator< double, double >(
                mapped_data,
                interpolatorSettings,
                std::make_pair( yValues( 0 ), yValues( yValues.size() - 1 ) ),
                bislip::Variables::HermiteDerivatives( xValues, yValues ) );
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

double evaluateGuidanceInterpolator (
        const bislip::Parameters::Optimization &parameter,
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    double currentAltitude = 0.0;
    double currentAirspeed = 0.0;

    if( debugInfo == 1 ){ std::cout << "Selecting source of airspeed/altitude" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Current Time  = " << flightConditions->getCurrentTime() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Starting Time = " << bislipSystems->getStartingEpoch() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Guidance Step = " << bislipSystems->getGuidanceStepSize() << std::endl; }


    if ( std::isnan( flightConditions->getCurrentTime() ) == true )
    {
        if( debugInfo == 1 ){ std::cout << "     Selecting initial airspeed/altitude" << std::endl; }

        currentAltitude = bislipSystems->getInitialAltitude();
        currentAirspeed = bislipSystems->getInitialAirspeed();
    }
    else if( std::isnan( flightConditions->getCurrentTime() ) == false )
    {
        if( debugInfo == 1 ){ std::cout << "     Selecting propagated airspeed/altitude" << std::endl; }

        currentAltitude = flightConditions->getCurrentAltitude();
        currentAirspeed = flightConditions->getCurrentAirspeed();
    }

    if( debugInfo == 1 ){ std::cout << "Current Altitude = " << currentAltitude << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Airspeed = " << currentAirspeed << std::endl; }
    if( debugInfo == 1 ){ std::cout << "E_max = " <<  bislipSystems->getE_max() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "E_hat = " <<  bislip::Variables::computeNormalizedSpecificEnergy( currentAltitude, currentAirspeed, bislipSystems->getE_max() ) << std::endl; }

    if( debugInfo == 1 ){ std::cout << "Evaluating Interpolator for: " << parameter << std::endl; }
    //! Evaluate interpolator.
    double evaluation = ( bislipSystems->getParameterInterpolator( parameter ) )->interpolate( bislip::Variables::computeNormalizedSpecificEnergy( currentAltitude, currentAirspeed, bislipSystems->getE_max() ) );

    if( debugInfo == 1 ){ std::cout << "Select parameter bounds" << std::endl; }
    //! Declare and assign parameter bounds.
    std::pair < double, double > bounds = bislipSystems->getParameterBounds( parameter );

    if( debugInfo == 1 ){ std::cout << "Impose parameter bounds" << std::endl; }
    //! Impose bounds.
    if ( evaluation < bounds.first ){ evaluation = bounds.first; }
    if ( evaluation > bounds.second ){ evaluation = bounds.second; }

    return evaluation;
}

Eigen::Vector3d computeBodyFixedThrustDirection (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "|||||||||||||||||||| computeBodyFixedThrustDirection" << std::endl; }


    //! Evaluate the Thrust Elevation Angle Interpolator and convert to radians.
    double thrustElevationAngle = bislipSystems->getCurrentThrustElevationAngle();

    if( debugInfo == 1 ){ std::cout << "Thrust Elevation Angle for BodyFixedThrustDirection= " << thrustElevationAngle << std::endl; }

    //! Evaluate the Thrust Azimuth Angle Interpolator and convert to radians.
    double thrustAzimuthAngle = bislipSystems->getCurrentThrustAzimuthAngle();

    if( debugInfo == 1 ){ std::cout << "Thrust Azimuth Angle for BodyFixedThrustDirection= " << thrustAzimuthAngle << std::endl; }

    //! Declare and calculate the body-fixed thrust direction.
    Eigen::Vector3d bodyFixedThrustDirection;
    bodyFixedThrustDirection( 0 ) = std::cos( thrustElevationAngle ) * std::cos( thrustAzimuthAngle );
    bodyFixedThrustDirection( 1 ) = std::cos( thrustElevationAngle ) * std::sin( thrustAzimuthAngle );
    bodyFixedThrustDirection( 2 ) = std::sin( thrustElevationAngle );

    return bodyFixedThrustDirection;
}

double computeThrustMagnitude (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    //! Evaluate the Throttle Setting Interpolator.
    const double throttle = bislipSystems->getCurrentThrottleSetting();
    if( debugInfo == 1 ){ std::cout << "Throttle Setting for computeThrustMagnitude = " << throttle << std::endl; }


    double currentThrustMagnitude = throttle * bislipSystems->getMaxThrust();
    if( debugInfo == 1 ){ std::cout << "Current Thrust Magnitude for computeThrustMagnitude = " << currentThrustMagnitude << std::endl; }


    return currentThrustMagnitude;
}

Eigen::Vector3d computeBodyFixedThrustVector (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();
    if( debugInfo == 1 ){ std::cout << "|||||||||||||||||||| computeBodyFixedThrustVector" << std::endl; }


    //! Determine the Body-Fixed Thrust Direction Vector.
    const Eigen::Vector3d BodyFixedThrustDirection = bislip::Variables::computeBodyFixedThrustDirection( bodyMap, vehicleName );
    if( debugInfo == 1 ){ std::cout << "Body-fixed Thrust Direction Vector = [ " << BodyFixedThrustDirection(0) << " , " << BodyFixedThrustDirection(1) << " , " << BodyFixedThrustDirection(2) << " ]" << std::endl; }


    //! Determine the Thrust Magnitude.
    const double thrustMagnitudeOutput = bislip::Variables::computeThrustMagnitudeOutput( bodyMap, vehicleName );
    if( debugInfo == 1 ){ std::cout << "Thrust Magnitude Output for computeBodyFixedThrustVector = " << thrustMagnitudeOutput << std::endl; }


    //! Declare and calculate the body-fixed thrust vector.
    Eigen::Vector3d BodyFixedThrustVector;
    BodyFixedThrustVector( 0 ) = thrustMagnitudeOutput * BodyFixedThrustDirection( 0 );
    BodyFixedThrustVector( 1 ) = thrustMagnitudeOutput * BodyFixedThrustDirection( 1 );
    BodyFixedThrustVector( 2 ) = thrustMagnitudeOutput * BodyFixedThrustDirection( 2 );

    if( debugInfo == 1 ){ std::cout << "Body-fixed Thrust Force Vector = [ " << BodyFixedThrustVector(0) << " , " << BodyFixedThrustVector(1) << " , " << BodyFixedThrustVector(2) << " ]" << std::endl; }

    return BodyFixedThrustVector;
}

double computeThrustMagnitudeOutput (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    const double engineStatus = double( bislip::Variables::determineEngineStatus( bodyMap.at( vehicleName )->getBodyMass(), bislipSystems->getLandingMass() ) );
    const double thrustMagnitude = bislip::Variables::computeThrustMagnitude( bodyMap, vehicleName );

    return thrustMagnitude * engineStatus;
}

bool determineEngineStatus (
        const double &currentMass,
        const double &landingMass )
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

double computeCurrentLiftForce(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName)
{

    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract coefficient interface pointer.
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Extract Vehicle Systems pointer.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( ) ;

    const Eigen::Vector6d currentCoefficients = bislip::Variables::computeFullCurrentCoefficients( bodyMap, vehicleName );

    return ( flightConditions->getCurrentDynamicPressure() ) * ( bislipSystems->getReferenceArea( ) ) * currentCoefficients[ 2 ];;
}

double computeSkipSuppressionLimit (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract coefficient interface pointer.
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    const double currentFlightPathAngle = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
    const double currentLatitude        = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double currentHeading         = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
    const double currentAltitude        = flightConditions->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );
    const double currentAirspeed        = flightConditions->getCurrentAirspeed( );
    //const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    const double currentMass            = bodyMap.at( vehicleName )->getBodyMass();

    const double rotationRateEarth = 7.292115*1E-5;

    //const Eigen::Vector6d currentCoefficients = bislip::Variables::computePartialCurrentCoefficients( bodyMap, vehicleName );
    const double currentLift = bislip::Variables::computeCurrentLiftForce( bodyMap, vehicleName );

    const Eigen::Vector2d gravs = bislip::Variables::computeLocalGravity( bodyMap, vehicleName, centralBodyName );

    const double currentAngleOfAttack        = bislipSystems->getCurrentAngleOfAttack();
    const double currentThrustMagnitude      = bislipSystems->getCurrentThrustMagnitude();
    const double currentThrustElevationAngle = bislipSystems->getCurrentThrustElevationAngle();
    const double currentThrustAzimuthAngle   = bislipSystems->getCurrentThrustAzimuthAngle();
    const double currentThrottleSetting      = bislipSystems->getCurrentThrottleSetting();


    const double s_alpha = std::sin( currentAngleOfAttack );
    const double c_alpha = std::cos( currentAngleOfAttack );
    const double s_delta = std::sin( currentLatitude );
    const double c_delta = std::cos( currentLatitude );
    const double s_gamma = std::sin( currentFlightPathAngle );
    const double c_gamma = std::cos( currentFlightPathAngle );
    const double s_chi = std::sin( currentHeading );
    const double c_chi = std::cos( currentHeading );
    const double s_thrustAzimuthAngle = std::sin( currentThrustAzimuthAngle );
    const double c_thrustAzimuthAngle = std::cos( currentThrustAzimuthAngle );
    const double s_thrustElevationAngle = std::sin( currentThrustElevationAngle );
    const double c_thrustElevationAngle = std::cos( currentThrustElevationAngle );

    //! Linear sum of sine and cosine waves with same frequency and phase shift.
    //!     From the EoM, the amplitudes correspond to the coefficients as follows:
    //!     a * cos( theta ) + b * sin( theta ) = c * cos( theta - phi )
    //!     c = sqrt( a * a + b * b )
    //!     phi = arctan( b, a )
    //!
    //!     This implementation, considering the additional terms from the EoM,
    //!     resulted in the following structure:
    //!     B * sin( sigma ) + C * cos( sigma ) = E * cos( theta - phi )
    //!                                         = -( D + A * m )
    //!     E = sqrt( B * B + C * C )
    //!     phi = arctan( B, C )

    const double A1 = 2.0 * rotationRateEarth * currentAirspeed * c_delta * s_chi;
    const double A2 = ( currentAirspeed * currentAirspeed / currentAltitude ) * c_gamma;
    const double A3 = rotationRateEarth * rotationRateEarth * currentAltitude * c_delta * ( c_delta * c_gamma + s_gamma * s_delta * c_chi );

    const double A = A1 + A2 + A3;
    const double B = currentThrustMagnitude * s_thrustAzimuthAngle * c_thrustElevationAngle;
    const double C = currentLift + currentThrustMagnitude * ( s_alpha * c_thrustAzimuthAngle * c_thrustElevationAngle + c_alpha * s_thrustElevationAngle );
    const double D = currentMass * ( gravs( 0 ) * s_gamma * c_chi - gravs( 1 ) * c_gamma );
    const double E = std::sqrt( ( B * B ) + ( C * C ) );
    const double arccosArgument = -( D + A * currentMass ) / E;

    //const double argument = -( currentMass / c ) * ( term1 + term2 + term3 - term4 + term5 );
    double limit = std::acos( arccosArgument ) + std::atan2( B, C );

    if( limit < 0.0 ) { limit = 0.0; }
    if( limit > tudat::mathematical_constants::PI / 2 ) { limit = 0.0; }
    if ( std::isnan( limit ) == true ) { limit = 0.0; }

    /*
    std::cout << "Current Phase:          " << bislipSystems->getCurrentTrajectoryPhase() << std::endl;
    std::cout << "CL:                     " << currentCoefficients[ 2 ] << std::endl;
    std::cout << "currentMass:            " << currentMass << std::endl;
    std::cout << "currentThrottleSetting: " << currentThrottleSetting << std::endl;
    std::cout << "currentThrustMagnitude: " << currentThrustMagnitude << std::endl;
    std::cout << "currentMass:            " << currentMass << std::endl;
    std::cout << "currentAirspeed:        " << currentAirspeed << std::endl;
    std::cout << "currentLift:            " << currentLift << std::endl;
    std::cout << "a:                      " << a << std::endl;
    std::cout << "b:                      " << b << std::endl;
    std::cout << "c:                      " << c << std::endl;
    std::cout << "phi:                    " << phi << std::endl;
    std::cout << "term1:                  " << term1 << std::endl;
    std::cout << "term2:                  " << term2 << std::endl;
    std::cout << "term3:                  " << term3 << std::endl;
    std::cout << "term4:                  " << term4 << std::endl;
    std::cout << "term5:                  " << term5 << std::endl;
    std::cout << "argument:               " << argument << std::endl;
   std::cout << "limit:                  " << limit << std::endl;
*/

    return limit;
}


double computeFlightPathAngleRate(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName )
{

    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract coefficient interface pointer.
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    const double currentFlightPathAngle = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
    const double currentLatitude        = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double currentHeading         = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
    const double currentAltitude        = flightConditions->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );
    const double currentAirspeed        = flightConditions->getCurrentAirspeed( );
    //const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    const double currentMass            = bodyMap.at( vehicleName )->getBodyMass();

    const double rotationRateEarth = 7.292115*1E-5;

    //const Eigen::Vector6d currentCoefficients = bislip::Variables::computePartialCurrentCoefficients( bodyMap, vehicleName );
    const double currentLift = bislip::Variables::computeCurrentLiftForce( bodyMap, vehicleName );

    const Eigen::Vector2d gravs = bislip::Variables::computeLocalGravity( bodyMap, vehicleName, centralBodyName );

    const double currentAngleOfAttack        = bislipSystems->getCurrentAngleOfAttack();
    const double currentBankAngle            = bislipSystems->getCurrentBankAngle();
    const double currentThrustMagnitude      = bislipSystems->getCurrentThrustMagnitude();
    const double currentThrustElevationAngle = bislipSystems->getCurrentThrustElevationAngle();
    const double currentThrustAzimuthAngle   = bislipSystems->getCurrentThrustAzimuthAngle();

    const double s_alpha = std::sin( currentAngleOfAttack );
    const double c_alpha = std::cos( currentAngleOfAttack );
    const double s_sigma = std::sin( currentBankAngle );
    const double c_sigma = std::cos( currentBankAngle );
    const double s_delta = std::sin( currentLatitude );
    const double c_delta = std::cos( currentLatitude );
    const double s_gamma = std::sin( currentFlightPathAngle );
    const double c_gamma = std::cos( currentFlightPathAngle );
    const double s_chi = std::sin( currentHeading );
    const double c_chi = std::cos( currentHeading );
    const double s_thrustAzimuthAngle = std::sin( currentThrustAzimuthAngle );
    const double c_thrustAzimuthAngle = std::cos( currentThrustAzimuthAngle );
    const double s_thrustElevationAngle = std::sin( currentThrustElevationAngle );
    const double c_thrustElevationAngle = std::cos( currentThrustElevationAngle );

    const double F_gamma_1 = currentThrustMagnitude * s_thrustAzimuthAngle * c_thrustElevationAngle * s_sigma;
    const double F_gamma_2 = ( currentLift + currentThrustMagnitude * ( s_alpha * c_thrustAzimuthAngle * c_thrustElevationAngle  + c_alpha * s_thrustElevationAngle ) ) * c_sigma;
    const double F_gamma_3 = currentMass * ( gravs( 0 ) * s_gamma * c_chi - gravs( 1 ) * c_gamma );

    const double F_gamma = F_gamma_1 + F_gamma_2 + F_gamma_3;

    const double A1 = 2.0 * rotationRateEarth * currentAirspeed * c_delta * s_chi;
    const double A2 = ( currentAirspeed * currentAirspeed / currentAltitude ) * c_gamma;
    const double A3 = rotationRateEarth * rotationRateEarth * currentAltitude * c_delta * ( c_delta * c_gamma + s_gamma * s_delta * c_chi );

    const double A = A1 + A2 + A3;

    const double flightPathAngleRate = ( ( F_gamma / currentMass ) + A ) / currentAirspeed;

    bislipSystems->setCurrentFlightPathAngleRate( flightPathAngleRate );


    return flightPathAngleRate;
}



double computeCumulativeCartesianDistanceTravelled (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    Eigen::Vector3d previousCoordinates = bislipSystems->getPreviousCartesianCoordinates();
    Eigen::Vector3d currentCoordinates  = ( flightConditions->getCurrentBodyCenteredBodyFixedState() ).segment( 0, 3 );
    bislipSystems->setPreviousCartesianCoordinates( currentCoordinates );

    double cartesianDistanceTravelled = std::sqrt( std::pow( currentCoordinates( 0 ) - previousCoordinates( 0 ), 2 ) +
                                                   std::pow( currentCoordinates( 1 ) - previousCoordinates( 1 ), 2 ) +
                                                   std::pow( currentCoordinates( 2 ) - previousCoordinates( 2 ), 2 ) );

    double cumulativeDistanceTravelled = bislipSystems->getCumulativeDistanceTravelled( ) + cartesianDistanceTravelled;

    if( debugInfo == 1 ){ std::cout << "Previous Coordinates = [ " << previousCoordinates(0) << ", "<<previousCoordinates(1) << ", "<<previousCoordinates(2) << " ]" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Coordinates  = [ " << currentCoordinates(0)  << ", "<<currentCoordinates(1)  << ", "<<currentCoordinates(2)  << " ]" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Cartesian Distance   = " << cartesianDistanceTravelled << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Cumulative Distance travelled = " << cumulativeDistanceTravelled << std::endl; }

    bislipSystems->setCumulativeDistanceTravelled( cumulativeDistanceTravelled );

    return cumulativeDistanceTravelled;
}



double computeCumulativeAngularDistanceTravelled (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    Eigen::Vector2d previousCoordinates = bislipSystems->getPreviousCoordinates();
    Eigen::Vector2d currentCoordinates( 2 );
    const double currentLatitude  = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double currentLongitude = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );
    currentCoordinates << currentLatitude, currentLongitude;

    bislipSystems->setPreviousCoordinates( currentCoordinates );

    double distanceTravelled = tudat::unit_conversions::convertRadiansToDegrees(
                bislip::Variables::computeAngularDistance( previousCoordinates( 0 ), previousCoordinates( 1 ),
                                                           currentCoordinates( 0 ), currentCoordinates( 1 ) ) );


    double cumulativeAngularDistanceTravelled = bislipSystems->getCumulativeAngularDistanceTravelled( ) + distanceTravelled;

    if( debugInfo == 1 ){ std::cout << "Previous Coordinates = [ " << previousCoordinates(0) << ", "<<previousCoordinates(1) << " ]" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Coordinates  = [ " << currentCoordinates(0)  << ", "<<currentCoordinates(1)  << " ]" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Angular Distance Travelled = " << distanceTravelled << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Cumulative Angular Distance Travelled = " << cumulativeAngularDistanceTravelled << std::endl; }

    bislipSystems->setCumulativeAngularDistanceTravelled( cumulativeAngularDistanceTravelled );

    return cumulativeAngularDistanceTravelled;
}


double computeCumulativeAngularDistanceTravelledDifference (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    double cumulativeAngularDisplacement = bislipSystems->getCumulativeAngularDistanceTravelled();

    const double currentLatitude  = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double currentLongitude = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );

    double netAngularDisplacement = tudat::unit_conversions::convertRadiansToDegrees(
                bislip::Variables::computeAngularDistance( bislipSystems->getInitialLat(), bislipSystems->getInitialLon(),
                                                           currentLatitude, currentLongitude ) );

    double groundtrackDifference = cumulativeAngularDisplacement - netAngularDisplacement;

    if( debugInfo == 1 ){ std::cout << "Cumulative Angular Displacement = " << cumulativeAngularDisplacement << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Net Angular Displacement        = " << netAngularDisplacement << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Groundtrack Difference          = " << groundtrackDifference << std::endl; }

    return groundtrackDifference;
}

double computeTimeOfFlight(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    return flightConditions->getCurrentTime() - bislipSystems->getStartingEpoch();;
}

double computeBendingMoment (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Extract current conditions.
    const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    const double currentAngleOfAttack   = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );

    return currentDynamicPressure * currentAngleOfAttack;
}


Eigen::Vector2d computeControlSurfaceDeflection (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract coefficient interface pointer.
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( );

    //! Extract Vehicle Systems pointer.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( );

    int debugInfo = bislipSystems->getDebugInfo();

    //! Extract current conditions.
    const double currentAngleOfAttack = bislipSystems->getCurrentAngleOfAttack();
    const double currentMachNumber    = flightConditions->getCurrentMachNumber();

    if( debugInfo == 1 ){ std::cout << "Create primary input for aerodynamic coefficient interface" << std::endl; }
    //! Create primary input for aerodynamic coefficient interface.
    //!     Take care of order of input (this depends on how the coefficients are created)!
    const std::vector < double > aerodynamicCoefficientInput = bislip::Variables::getAerodynamicCoefficientInput( currentAngleOfAttack, currentMachNumber );


    if( debugInfo == 1 ){ std::cout << "Create control surface input for aerodynamic coefficient interface" << std::endl; }
    //! Create control surface input for aerodynamic coefficient interface.
    //!     Take care of order of input (this depends on how the coefficients are created)!
    std::map< std::string, std::vector< double > > controlSurfaceCoefficientInput = bislip::Variables::getControlSurfaceCoefficientInput( currentAngleOfAttack, currentMachNumber, 0.0, 0.0 );


    Eigen::VectorXd CmBound ( 2 );
    //Eigen::VectorXd absCmBound ( 2 );

    controlSurfaceCoefficientInput[ "BodyFlap" ][ 2 ] = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getBodyFlapDeflectionLimits() ).first );
    controlSurfaceCoefficientInput[ "ElevonLeft" ][ 2 ] = tudat::unit_conversions::convertDegreesToRadians( 0.0 );
    controlSurfaceCoefficientInput[ "ElevonRight" ][ 2 ] = tudat::unit_conversions::convertDegreesToRadians( 0.0 );
    coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
    CmBound( 0 ) = coefficientInterface->getCurrentMomentCoefficients( )( 1 );

    if( debugInfo == 1 ){ std::cout << "    C_m( delta_bf_min, 0.0 ) = " << CmBound( 0 ) << std::endl; }

    controlSurfaceCoefficientInput[ "BodyFlap" ][ 2 ] = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getBodyFlapDeflectionLimits() ).second );
    //controlSurfaceCoefficientInput[ "ElevonLeft" ][ 3 ] = tudat::unit_conversions::convertDegreesToRadians( 0.0 );
    //controlSurfaceCoefficientInput[ "ElevonRight" ][ 4 ] = tudat::unit_conversions::convertDegreesToRadians( 0.0 );
    coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
    CmBound( 1 ) = coefficientInterface->getCurrentMomentCoefficients( )( 1 );

    if( debugInfo == 1 ){ std::cout << "    C_m( delta_bf_max, 0.0 ) = " << CmBound( 1 ) << std::endl; }

    //absCmBound( 0 ) = std::abs( CmBound ( 0 ) );
    //absCmBound( 1 ) = std::abs( CmBound ( 1 ) );

    double a = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getBodyFlapDeflectionLimits() ).first );
    double b = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getBodyFlapDeflectionLimits() ).second );
    double delta_bf = 0.0;
    double f_a = 0.0;
    double f_delta_bf = 0.0;

    Eigen::Vector2d controlSurfaceDeflections( 2 );


    if ( CmBound( 0 ) * CmBound( 1 ) < 0 )
    {
        delta_bf = ( a + b ) / 2.0;

        controlSurfaceCoefficientInput[ "BodyFlap" ][ 2 ] = delta_bf;
        coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );

        if( debugInfo == 1 ){ std::cout << "     Starting bodyflap search:" << std::endl; }
        while ( std::abs( coefficientInterface->getCurrentMomentCoefficients( )( 1 ) ) > 1e-7 )
        {

            controlSurfaceCoefficientInput[ "BodyFlap" ][ 2 ] = a;
            coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
            f_a       = coefficientInterface->getCurrentMomentCoefficients( )( 1 );


            controlSurfaceCoefficientInput[ "BodyFlap" ][ 2 ] = delta_bf;
            coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
            f_delta_bf = coefficientInterface->getCurrentMomentCoefficients( )( 1 );

            if (f_a * f_delta_bf < 0 ) { b = delta_bf; }
            else { a = delta_bf; }
            delta_bf = ( a + b ) / 2.0;
            controlSurfaceCoefficientInput[ "BodyFlap" ][ 2 ] = delta_bf;
            coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
            // err = std::abs( coefficientInterface->getCurrentMomentCoefficients( )( 1 ) );

        }
        if( debugInfo == 1 ){ std::cout << "     delta_bf = " << delta_bf << "     |     " << "C_m( delta_bf, 0.0 ) = " << std::abs( coefficientInterface->getCurrentMomentCoefficients( )( 1 ) ) << std::endl; }

        controlSurfaceDeflections << delta_bf, 0.0 ;

    }
    else
    {
        const double currentBodyFlapAngleDeflection = vehicleSystems->getCurrentControlSurfaceDeflection( "BodyFlap" );
        //const double currentElevonAngle = vehicleSystems->getCurrentControlSurfaceDeflection( "ElevonLeft" );
        const double elevonAngleLowerLimit = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getElevonDeflectionLimits() ).first );
        const double elevonAngleUpperLimit = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getElevonDeflectionLimits() ).second );

        controlSurfaceDeflections( 0 ) = currentBodyFlapAngleDeflection;

        controlSurfaceCoefficientInput[ "BodyFlap" ][ 2 ] = currentBodyFlapAngleDeflection;
        controlSurfaceCoefficientInput[ "ElevonLeft" ][ 2 ] = elevonAngleLowerLimit;
        controlSurfaceCoefficientInput[ "ElevonRight" ][ 2 ] = elevonAngleLowerLimit;
        coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
        CmBound( 0 ) = coefficientInterface->getCurrentMomentCoefficients( )( 1 );

        //controlSurfaceCoefficientInput[ "BodyFlap" ][ 2 ] = tudat::unit_conversions::convertDegreesToRadians( 0.0 );
        controlSurfaceCoefficientInput[ "ElevonLeft" ][ 2 ] = elevonAngleUpperLimit;
        controlSurfaceCoefficientInput[ "ElevonRight" ][ 2 ] = elevonAngleUpperLimit;
        coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
        CmBound( 1 ) = coefficientInterface->getCurrentMomentCoefficients( )( 1 );

        if( debugInfo == 1 ){ std::cout << "   C_m( delta_bf, delta_el_min ) = " << CmBound( 0 ) << std::endl; }
        if( debugInfo == 1 ){ std::cout << "   C_m( delta_bf, delta_el_max ) = " << CmBound( 1 ) << std::endl; }

        if ( CmBound( 0 ) * CmBound( 1 ) < 0 )
        {

            a = elevonAngleLowerLimit;
            b = elevonAngleUpperLimit;
            double delta_el = ( a + b ) / 2;
            double f_delta_el;
            controlSurfaceCoefficientInput[ "ElevonLeft" ][ 2 ] = delta_el;
            controlSurfaceCoefficientInput[ "ElevonRight" ][ 2 ] = delta_el;
            coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );

            if( debugInfo == 1 ){ std::cout << "     Starting elevon search:" << std::endl; }
            while ( std::abs( coefficientInterface->getCurrentMomentCoefficients( )( 1 ) ) > 1e-7  )
            {

                controlSurfaceCoefficientInput[ "ElevonLeft" ][ 2 ] = a;
                controlSurfaceCoefficientInput[ "ElevonRight" ][ 2 ] = a;
                coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
                f_a       = coefficientInterface->getCurrentMomentCoefficients( )( 1 );

                controlSurfaceCoefficientInput[ "ElevonLeft" ][ 2 ] = delta_el;
                controlSurfaceCoefficientInput[ "ElevonRight" ][ 2 ] = delta_el;
                coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
                f_delta_el = coefficientInterface->getCurrentMomentCoefficients( )( 1 );

                if ( f_a * f_delta_el < 0 ) { b = delta_el; }
                else { a = delta_el; }

                delta_el = ( a + b ) / 2;
                controlSurfaceCoefficientInput[ "ElevonLeft" ][ 2 ] = delta_el;
                controlSurfaceCoefficientInput[ "ElevonRight" ][ 2 ] = delta_el;
                coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
                // err = std::abs( coefficientInterface->getCurrentMomentCoefficients( )( 1 ) );
                if( debugInfo == 1 ){ std::cout << "     delta_bf = " << controlSurfaceDeflections( 0 ) << "     |     " << "     delta_el = " << delta_el << "     |     " << "C_m( delta_bf, delta_el ) = " << std::abs( coefficientInterface->getCurrentMomentCoefficients( )( 1 ) ) << std::endl; }

            }
            controlSurfaceDeflections( 1 ) = delta_el;
        }
    }

    if( debugInfo == 1 ){ std::cout << "PA FUERA!" << std::endl; }

    return controlSurfaceDeflections;
}


double computeFullPitchMomentCoefficient(
        const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
        const double &bodyFlapDeflection,
        const double &elevonDeflection,
        const std::vector< double > &aerodynamicCoefficientInput )
{
    //! Define values of independent variables of control surface aerodynamics.
    std::map< std::string, std::vector< double > > controlSurfaceCoefficientInput = bislip::Variables::getControlSurfaceCoefficientInput( aerodynamicCoefficientInput[ 0 ], aerodynamicCoefficientInput[ 1 ], bodyFlapDeflection, elevonDeflection );

    //! Update full coefficient interface.
    coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );

    return coefficientInterface->getCurrentMomentCoefficients( )( 1 );
}

double computePitchMomentCoefficient(
        const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
        const std::vector< double > &aerodynamicCoefficientInput )
{
    //! Update coefficient interface.
    coefficientInterface->updateCurrentCoefficients( aerodynamicCoefficientInput );

    return coefficientInterface->getCurrentMomentCoefficients( )( 1 );
}

double computeBodyFlapCmIncrementdif (
        const double &full,
        const double &partial)
{
    return full - partial;
}


Eigen::Vector6d computePartialCurrentCoefficients(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract coefficient interface pointer.
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Define input to aerodynamic coefficients: take care of order of input (this depends on how the coefficients are created)!
    std::vector< double > aerodynamicCoefficientInput = bislip::Variables::getAerodynamicCoefficientInput( bislipSystems->getCurrentAngleOfAttack(), flightConditions->getCurrentMachNumber() );

    //! Update and retrieve current aerodynamic coefficients
    coefficientInterface->updateCurrentCoefficients( aerodynamicCoefficientInput );

    return coefficientInterface->getCurrentAerodynamicCoefficients( );
}

Eigen::Vector6d computeFullCurrentCoefficients(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract coefficient interface pointer.
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Extract Vehicle Systems pointer.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "    computeFullCurrentCoefficients -----> Extracting current data" << std::endl; }

    const double currentAngleOfAttack           = bislipSystems->getCurrentAngleOfAttack();
    const double currentMachNumber              = flightConditions->getCurrentMachNumber();
    const double currentBodyFlapAngleDeflection = vehicleSystems->getCurrentControlSurfaceDeflection( "BodyFlap" );
    const double currentElevonAngleDeflection   = vehicleSystems->getCurrentControlSurfaceDeflection( "ElevonLeft" );
    if( debugInfo == 1 ){ std::cout << "                                          currentAngleOfAttack           = " << currentAngleOfAttack << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                                          currentMachNumber              = " << currentMachNumber << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                                          currentBodyFlapAngleDeflection = " << currentBodyFlapAngleDeflection << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                                          currentElevonAngleDeflection   = " << currentElevonAngleDeflection << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    computeFullCurrentCoefficients -----> Create primary input for aerodynamic coefficient interface." << std::endl; }
    //! Create primary input for aerodynamic coefficient interface.
    //!     Take care of order of input (this depends on how the coefficients are created)!
    std::vector< double > aerodynamicCoefficientInput = bislip::Variables::getAerodynamicCoefficientInput( currentAngleOfAttack, currentMachNumber );


    if( debugInfo == 1 ){ std::cout << "    computeFullCurrentCoefficients -----> Create control surface input for aerodynamic coefficient interface." << std::endl; }
    //! Create control surface input for aerodynamic coefficient interface.
    //!     Take care of order of input (this depends on how the coefficients are created)!
    std::map< std::string, std::vector< double > > controlSurfaceCoefficientInput = bislip::Variables::getControlSurfaceCoefficientInput( currentAngleOfAttack, currentMachNumber, currentBodyFlapAngleDeflection, currentElevonAngleDeflection );

    if( debugInfo == 1 ){ std::cout << "    computeFullCurrentCoefficients -----> Update and retrieve current aerodynamic coefficients." << std::endl; }
    //! Update and retrieve current aerodynamic coefficients.
    coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );

    if( debugInfo == 1 ){ std::cout << "    computeFullCurrentCoefficients -----> returning result" << std::endl; }

    return coefficientInterface->getCurrentAerodynamicCoefficients( );
}


double computeBodyFlapCmIncrement (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
    //            bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();


    if( debugInfo == 1 ){ std::cout << "    computeBodyFlapCmIncrement -----> Extracting current dynamic pressure" << std::endl; }

    const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    //const Eigen::Vector3d referenceValues = bislipSystems->getReferenceValues();
    //const Eigen::Vector3d aerodynamicReferenceCenter = bislipSystems->getMassReferenceCenter();
    //const Eigen::Vector3d thrustReferenceCenter = bislipSystems->getThrustReferenceCenter();
    // const double surfaceArea = referenceValues[ 0 ];
    const double surfaceArea = bislipSystems->getReferenceValues( )( 0 );

    // const double b_ref = referenceValues[ 1 ];
    //    const double c_ref = referenceValues[ 2 ];
    const double c_ref = bislipSystems->getReferenceValues( )( 2 );
    //const double del_x = aerodynamicReferenceCenter[ 0 ];
    const double del_x = bislipSystems->getMassReferenceCenter( )( 0 );
    // const double del_x_T = thrustReferenceCenter[ 0 ];
    //const double del_z_T = thrustReferenceCenter[ 2 ];
    const double del_x_T = bislipSystems->getThrustReferenceCenter( )( 0 );
    const double del_z_T = bislipSystems->getThrustReferenceCenter( )( 2 );

    //Eigen::Vector6d currentCoefficients = bislip::Variables::computePartialCurrentCoefficients( bodyMap, vehicleName );
    const Eigen::Vector6d currentCoefficients = bislip::Variables::computePartialCurrentCoefficients( bodyMap, vehicleName );
    const double currentThrustMagnitude = bislipSystems->getCurrentThrustMagnitude();
    const double currentDragCoefficient = currentCoefficients( 0 );
    const double currentLiftCoefficient = currentCoefficients( 2 );
    const double currentPitchMomentCoefficient = currentCoefficients( 4 );



    if( debugInfo == 1 ){ std::cout << "    computeBodyFlapCmIncrement -----> Extracting current angle of attack" << std::endl; }

    double angleOfAttack = bislipSystems->getCurrentAngleOfAttack();

    if( debugInfo == 1 ){ std::cout << "    computeBodyFlapCmIncrement -----> Extracting current thrust elevation angle" << std::endl; }

    double thrustElevationAngle = bislipSystems->getCurrentThrustElevationAngle();

    /*


            tudat::unit_conversions::convertDegreesToRadians(
                bislip::Variables::evaluateGuidanceInterpolator (
                    bislip::Parameters::Optimization::ThrustElevationAngle,
                    bodyMap,
                    vehicleName ) );*/
    if( debugInfo == 1 ){ std::cout << "    computeBodyFlapCmIncrement -----> Extracting current thrust azimuth angle" << std::endl; }

    double thrustAzimuthAngle = bislipSystems->getCurrentThrustAzimuthAngle();


    /*tudat::unit_conversions::convertDegreesToRadians(
                bislip::Variables::evaluateGuidanceInterpolator (
                    bislip::Parameters::Optimization::ThrustAzimuthAngle,
                    bodyMap,
                    vehicleName ) );
                    */

    if( debugInfo == 1 ){ std::cout << "    computeBodyFlapCmIncrement -----> Calculating current drag" << std::endl; }

    const double currentDrag = currentDynamicPressure * surfaceArea * currentDragCoefficient;

    if( debugInfo == 1 ){ std::cout << "    computeBodyFlapCmIncrement -----> Calculating current lift" << std::endl; }

    const double currentLift = currentDynamicPressure * surfaceArea * currentLiftCoefficient;


    double term1 = ( del_x / currentDynamicPressure * surfaceArea * c_ref ) * ( currentDrag * std::sin( angleOfAttack ) - currentLift * std::cos( angleOfAttack ) );
    double term2 = ( currentThrustMagnitude / ( currentDynamicPressure * surfaceArea * c_ref ) ) * ( del_x_T * std::sin( thrustElevationAngle ) + del_z_T * std::cos( thrustElevationAngle ) * std::cos( thrustAzimuthAngle ) );
    /*
    std::cout << "del_x:                  " << del_x << std::endl;
    std::cout << "del_x_T:                " << del_x_T << std::endl;
    std::cout << "surfaceArea:            " << surfaceArea << std::endl;
    std::cout << "c_ref:                  " << c_ref << std::endl;
    std::cout << "currentDrag:            " << currentDrag << std::endl;
    std::cout << "currentLift:            " << currentLift << std::endl;
    std::cout << "currentThrustMagnitude: " << currentThrustMagnitude << std::endl;
    std::cout << "angleOfAttack:          " << angleOfAttack << std::endl;
    std::cout << "thrustElevationAngle:   " << thrustElevationAngle << std::endl;
    std::cout << "thrustAzimuthAngle:     " << thrustAzimuthAngle << std::endl;
    std::cout << "term1: " << term1 << std::endl;
    std::cout << "term2: " << term2 << std::endl;
    std::cout << "currentPitchMomentCoefficient: " << currentPitchMomentCoefficient  << std::endl;
    */
    if( debugInfo == 1 ){ std::cout << "    computeBodyFlapCmIncrement = " << ( term1 - term2 - currentPitchMomentCoefficient ) << std::endl; }


    //std::cout << "computeBodyFlapCmIncrement: " << term1 - term2 - currentPitchMomentCoefficient  << std::endl;

    return ( term1 - term2 - currentPitchMomentCoefficient );
}

//! Create input for aerodynamic coefficient interface.
std::vector< double > getAerodynamicCoefficientInput(
        const double &angleOfAttack,
        const double &machNumber )
{
    std::vector< double > aerodynamicCoefficientInput( 2 );
    aerodynamicCoefficientInput[ 0 ] = angleOfAttack;
    aerodynamicCoefficientInput[ 1 ] = machNumber;

    return aerodynamicCoefficientInput;
}

//! Create input for aerodynamic coefficient interface with a control surface (BodyFlap).
std::map< std::string, std::vector< double > > getControlSurfaceCoefficientInput(
        const double &angleOfAttack,
        const double &machNumber,
        const double &bodyFlapDeflectionAngle,
        const double &elevonDeflectionAngle )
{
    std::map< std::string, std::vector< double > > controlSurfaceCoefficientInput;
    controlSurfaceCoefficientInput[ "BodyFlap" ] = bislip::Variables::getAerodynamicCoefficientInput( angleOfAttack, machNumber );
    controlSurfaceCoefficientInput[ "BodyFlap" ].push_back( bodyFlapDeflectionAngle );
    controlSurfaceCoefficientInput[ "ElevonLeft" ] = bislip::Variables::getAerodynamicCoefficientInput( angleOfAttack, machNumber );
    controlSurfaceCoefficientInput[ "ElevonLeft" ].push_back( elevonDeflectionAngle );
    controlSurfaceCoefficientInput[ "ElevonRight" ] = bislip::Variables::getAerodynamicCoefficientInput( angleOfAttack, machNumber );
    controlSurfaceCoefficientInput[ "ElevonRight" ].push_back( elevonDeflectionAngle );

    return controlSurfaceCoefficientInput;
}

//! Compute the total load experienced by a vehicle.
Eigen::Vector3d computeBodyFixedTotalLoad(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "    computeBodyFixedTotalLoad -----> Calculating the Body-Fixed Aerodynamic Load" << std::endl; }

    Eigen::Vector3d bodyfixedAerodynamicLoad = bislip::Variables::computeBodyFixedAerodynamicLoad( bodyMap, vehicleName );

    if( debugInfo == 1 ){ std::cout << "    computeBodyFixedTotalLoad -----> Calculating the Body-Fixed Thrust Load" << std::endl; }

    Eigen::Vector3d bodyfixedThrustLoad = bislip::Variables::computeBodyFixedThrustVector( bodyMap, vehicleName );


    if( debugInfo == 1 ){ std::cout << "    computeBodyFixedTotalLoad -----> Calculating the Total Body-Fixed Load" << std::endl; }
    Eigen::Vector3d totalLoadBodyFrame = bodyfixedAerodynamicLoad + bodyfixedThrustLoad;
    if( debugInfo == 1 ){ std::cout << "Total Body-fixed Load = [ " << totalLoadBodyFrame(0) << " , " << totalLoadBodyFrame(1) << " , " << totalLoadBodyFrame(2) << " ]" << std::endl; }



    return totalLoadBodyFrame ;
}

//! Compute the aerodynamic load experienced by a vehicle.
Eigen::Vector3d computeBodyFixedAerodynamicLoad(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract coefficient interface pointer.
    std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > coefficientInterface = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                bodyMap.at( vehicleName )->getAerodynamicCoefficientInterface( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    //! Extract Vehicle Systems pointer.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( ) ;

    if( debugInfo == 1 ){ std::cout << "    computeBodyFixedAerodynamicLoad -----> Extracting current data" << std::endl; }

    //! Extract current conditions.
    const double currentAngleOfAttack   = bislipSystems->getCurrentAngleOfAttack();
    const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    //const double airDensity    = flightConditions->getCurrentDensity();
    //const double airSpeed      = flightConditions->getCurrentAirspeed();
    //const double machNumber    = flightConditions->getCurrentMachNumber();
    if( debugInfo == 1 ){ std::cout << "                                           currentAngleOfAttack           = " << currentAngleOfAttack << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                                           currentDynamicPressure         = " << currentDynamicPressure << std::endl; }

    //! Extract information from Bislip/Vehicle Systems.
    const double referenceArea = bislipSystems->getReferenceValues( )( 0 );
    // const double bodyFlapAngleDeflection = vehicleSystems->getCurrentControlSurfaceDeflection( "BodyFlap" );
    if( debugInfo == 1 ){ std::cout << "                                           reference area                 = " << currentAngleOfAttack << std::endl; }

    //! Calculate aerodynamic load in the Aerodynamic Frame.
    Eigen::Vector3d aerodynamicLoad = currentDynamicPressure * referenceArea * ( bislip::Variables::computeFullCurrentCoefficients( bodyMap, vehicleName ) ).segment( 0, 3 ) ;
    if( debugInfo == 1 ){ std::cout << "Aerodynamic Load = [ " << aerodynamicLoad(0) << " , " << aerodynamicLoad(1) << " , " << aerodynamicLoad(2) << " ]" << std::endl; }



    //! Transform aerodynamic load from the Aerodynamic Frame to the Body Frame.
    Eigen::Vector3d aerodynamicLoadBodyFrame = tudat::reference_frames::getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix( currentAngleOfAttack, 0.0 ) * aerodynamicLoad ;
    if( debugInfo == 1 ){ std::cout << "Body-fixed Aerodynamic Load = [ " << aerodynamicLoadBodyFrame(0) << " , " << aerodynamicLoadBodyFrame(1) << " , " << aerodynamicLoadBodyFrame(2) << " ]" << std::endl; }



    return aerodynamicLoadBodyFrame ;
}

Eigen::Vector3d computeBodyFixedTotal_g_Load_Vector (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Determine the Body-Fixed Total Load g-load vector.
    Eigen::Vector3d totalLoadBodyFrame_g_Load_Vector = ( bislip::Variables::computeBodyFixedTotalLoad( bodyMap, vehicleName ) ) / ( tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * bodyMap.at( vehicleName )->getBodyMass() );

    return totalLoadBodyFrame_g_Load_Vector;
}

double computeBodyFixedTotal_g_Load_Magnitude (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{

    //! Determine the Body-Fixed Total Load g-load magnitude.
    double totalLoadBodyFrame_g_Load_Magnitude = ( bislip::Variables::computeBodyFixedTotal_g_Load_Vector( bodyMap, vehicleName ) ).norm( );

    return totalLoadBodyFrame_g_Load_Magnitude;
}

Eigen::Matrix3d computeRotationMatrixONE( const double phi )
{
    Eigen::Matrix3d rotationMatrixONE;
    rotationMatrixONE << 1, 0, 0, 0, std::cos( phi ), std::sin( phi ), 0, -std::sin( phi ), std::cos( phi );

    return  rotationMatrixONE;
}

Eigen::Matrix3d computeRotationMatrixTHREE( const double psi )
{
    Eigen::Matrix3d rotationMatrixTHREE;
    rotationMatrixTHREE << std::cos( psi ), std::sin( psi ), 0, -std::sin( psi ), std::cos( psi ), 0, 0, 0, 1;

    return  rotationMatrixTHREE;
}

Eigen::Vector3d computePassengerFrameTotal_g_Load_Vector(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Determine the Passenger-Fixed Total Load g-load vector.
    Eigen::Vector3d totalLoadPassengerFrame_g_Load_Vector = ( bislipSystems->getBodyFrameToPassengerFrameTransformationMatrix() ) * ( bislip::Variables::computeBodyFixedTotal_g_Load_Vector( bodyMap, vehicleName ) );

    return totalLoadPassengerFrame_g_Load_Vector;
}

double computeBankAngle(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    // std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
    //           bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    //! Declare and determine the New Bank Angle variable.
    double newBankAngle = bislip::Variables::determineNewBankAngle( bodyMap, vehicleName );





    //! Used entirely for dependent variable calculation.
    bislipSystems->setTempBankAngle( newBankAngle );


    /*
    //! Determine if bank angle reversal is required.
    const bool bankReversal = bislip::Variables::determineBankAngleReversal( bodyMap, vehicleName, newBankAngle );

    double returnedBankAngle = newBankAngle;
    //! Impose bank angle reversal on new bank angle.
    if ( bankReversal == true ) { returnedBankAngle = -newBankAngle; }

    if( debugInfo == 10 ){ std::cout << "newBankAngle = " << newBankAngle << "   |   returnedBankAngle = " << returnedBankAngle << "   |   bislipSystems->getBankAngleReversalTrigger( ) = " << bislipSystems->getBankAngleReversalTrigger(  ) << std::endl; }
*/

    double returnedBankAngle = bislip::Variables::returnReversedBankAngle( bodyMap, vehicleName, newBankAngle );

    //if( debugInfo == 10 ){ std::cout << "newBankAngle = " << newBankAngle << "   |   returnedBankAngle = " << returnedBankAngle << "   |   bislipSystems->getBankAngleReversalTrigger( ) = " << bislipSystems->getBankAngleReversalTrigger(  ) << std::endl; }

    return returnedBankAngle;
}



double returnReversedBankAngle(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const double &newBankAngle )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    const bool bankReversal = bislip::Variables::determineBankAngleReversal( bodyMap, vehicleName, newBankAngle );

    double returnedReversedBankAngle = newBankAngle;

    //! Impose bank angle reversal on new bank angle.
    if ( bankReversal == true ) { returnedReversedBankAngle = -newBankAngle; }

    if( debugInfo == 2 ){ std::cout << "newBankAngle = " << newBankAngle << "   |   returnedBankAngle = " << returnedReversedBankAngle << "   |   bislipSystems->getBankAngleReversalTrigger( ) = " << bislipSystems->getBankAngleReversalTrigger(  ) << std::endl; }

    return returnedReversedBankAngle;
}

double determineNewBankAngle(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Declare and initialize the New Bank Angle variable.
    double newBankAngle = 0.0;

    //! Check trajectory phase to fix Ascent bank angle. Otherwise, interpolate a new Bank Angle.
    //  if ( bislipSystems->getCurrentTrajectoryPhase() == "Descent" )
    //{
    //! Interpolate new bank angle and assign previous SIGN for consistency.
    newBankAngle = bislip::Variables::determineSignedBankAngle(
                bislip::Variables::determineBankAngleSign( bislipSystems->getCurrentBankAngle() ),
                std::abs( bislipSystems->getEvaluatedBankAngle( ) ) );
    //}

    return newBankAngle;
}

int determineBankAngleSign ( const double &bankAngle )
{
    //! This ensures that the new bank angle has the same sign as in current bank angle.
    int sign = 1;
    if ( bankAngle < 0.0 ){ sign = -1; }
    return sign;
}


double determineSignedBankAngle ( const int &sign, const double &newbankAngle )
{ return double( sign ) * newbankAngle; }


double determineReversalConditional (
        const double &bankAngle,
        const double &headingError )
{
    //! Determine bank angle Reversal conditional. If result is lower than zero,
    //! bank reversal could happen. Should be less than zero if the vehicle is
    //! in fact aiming away from the target.

    //! This will happen if the bank angle is positive and the heading error is
    //! negative or if the bank angle is negative and the heading error is positive.

    return  bankAngle * headingError;

}

bool determineBankAngleReversal(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const double &newBankAngle )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( );

    //! Extract current conditions.
    // const double currentLatitude  = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    // const double currentLongitude = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );
    const double currentHeading   = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );

    //! Calculate angular distance to go.
    //const double angularDistanceToGo_deg = tudat::unit_conversions::convertRadiansToDegrees(
    //          bislip::Variables::computeAngularDistance (
    //            currentLatitude,
    //          currentLongitude,
    //        bodyMap.at( vehicleName )->getBislipSystems()->getTargetLat(),
    //      bodyMap.at( vehicleName )->getBislipSystems()->getTargetLon() ) );

    //! Calculate current Heading Angle Error
    const double currentHeadingAngleError_deg = tudat::unit_conversions::convertRadiansToDegrees(
                bislip::Variables::computeHeadingError (
                    flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle ),
                    flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle ),
                    bodyMap.at( vehicleName )->getBislipSystems()->getTargetLat(),
                    bodyMap.at( vehicleName )->getBislipSystems()->getTargetLon(),
                    currentHeading ) );

    //! Determine if vehicle is aiming away from target.
    const double reversalConditional = bislip::Variables::determineReversalConditional( newBankAngle, currentHeadingAngleError_deg );

    bislipSystems->setReversalConditional( reversalConditional );

    //! Determine if bank reversal should occur.
    bool reversal = false;

    if ( ( reversalConditional < 0.0 ) && ( std::abs( currentHeadingAngleError_deg ) >= bislip::Variables::computeHeadingErrorDeadBand( bodyMap, vehicleName ) ) )
    { reversal = true; }

    //bislipSystems->setBankAngleReversalTrigger( reversal );



    return reversal;
}

double computeHeadingErrorDeadBand(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( );

    //! Extract current conditions.
    const double currentLatitude  = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double currentLongitude = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );

    //! Calculate angular distance to go.
    const double angularDistanceToGo_deg = tudat::unit_conversions::convertRadiansToDegrees(
                bislip::Variables::computeAngularDistance (
                    currentLatitude,
                    currentLongitude,
                    bodyMap.at( vehicleName )->getBislipSystems()->getTargetLat(),
                    bodyMap.at( vehicleName )->getBislipSystems()->getTargetLon() ) );

    return bislipSystems->getHeadingErrorDeadBandInterpolator()->interpolate( angularDistanceToGo_deg );
}

double convertRadiansToDegrees( const double &angleInRadians )
{ return angleInRadians * 180.0 / tudat::mathematical_constants::PI; }

double convertDegreesToRadians( const double &angleInDegrees )
{ return angleInDegrees * tudat::mathematical_constants::PI / 180.0; }

double convertNegativeAnglesToPositive( const double &angleInRadians )
{
    double positiveAngle = angleInRadians;
    if( angleInRadians < 0.0 ) { positiveAngle = angleInRadians + 2 * tudat::mathematical_constants::PI; }
    return positiveAngle;
}

double determineAbsoluteValue( const double &value )
{ return std::abs( value ); }



void createAlphaMachBoundingInterpolators(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::vector< double > &alphaMachEnvelopeUB,
        const std::vector< double > &alphaMachEnvelopeLB,
        const std::string &outputPath,
        const std::string &outputSubFolder )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //!
    std::map< double, double > map_alphaMachEnvelopeUB, map_alphaMachEnvelopeLB;

    //! Declare and initialize alpha-Mach Interpolator Settings.
    std::shared_ptr< tudat::interpolators::InterpolatorSettings > alphaMachEnvelopeInterpolatorSettings = std::make_shared< tudat::interpolators::InterpolatorSettings >( tudat::interpolators::InterpolatorTypes::linear_interpolator );

    //! Declare and initialize size vector size.
    unsigned long alphaMachEnvelopeUBSize = alphaMachEnvelopeUB.size() / 2;
    unsigned long alphaMachEnvelopeLBSize = alphaMachEnvelopeLB.size() / 2;

    //! Declare alpha-Mach bound vectors.
    Eigen::VectorXd alphaMachEnvelopeUB_alpha( alphaMachEnvelopeUBSize );
    Eigen::VectorXd alphaMachEnvelopeUB_Mach(alphaMachEnvelopeUBSize );

    Eigen::VectorXd alphaMachEnvelopeLB_alpha( alphaMachEnvelopeLBSize );
    Eigen::VectorXd alphaMachEnvelopeLB_Mach(alphaMachEnvelopeLBSize );

    //! Populate vectors that will be used with data map.
    for( unsigned long i = 0; i < alphaMachEnvelopeUBSize; i++ )
    {
        alphaMachEnvelopeUB_Mach( i )  =  alphaMachEnvelopeUB[ i ];
        alphaMachEnvelopeUB_alpha( i ) =  alphaMachEnvelopeUB[ i + alphaMachEnvelopeUBSize ];
    }

    //! Populate vectors that will be used with data map.
    for( unsigned long i = 0; i < alphaMachEnvelopeLBSize; i++ )
    {
        alphaMachEnvelopeLB_Mach( i )  =  alphaMachEnvelopeLB[ i ];
        alphaMachEnvelopeLB_alpha( i ) =  alphaMachEnvelopeLB[ i + alphaMachEnvelopeLBSize ];
    }

    //! Populate data map for interpolator.
    for ( unsigned long i = 0; i < alphaMachEnvelopeUBSize; ++i )
    { map_alphaMachEnvelopeUB[ alphaMachEnvelopeUB_Mach( i ) ] = alphaMachEnvelopeUB_alpha( i ); }

    //! Populate data map for interpolator.
    for ( unsigned long i = 0; i < alphaMachEnvelopeLBSize; ++i )
    { map_alphaMachEnvelopeLB[ alphaMachEnvelopeLB_Mach( i ) ] = alphaMachEnvelopeLB_alpha( i ); }

    std::pair< double, double > alphaMachEnvelopeUBInterpolatorDomainInterval = std::make_pair( alphaMachEnvelopeUB_Mach.minCoeff(), alphaMachEnvelopeUB_Mach.maxCoeff() );
    std::pair< double, double > alphaMachEnvelopeUBInterpolatorRangeInterval  = std::make_pair( alphaMachEnvelopeUB_alpha.minCoeff(), alphaMachEnvelopeUB_alpha.maxCoeff() );
    std::pair< double, double > alphaMachEnvelopeUBInterpolatorBoundaryValues = std::make_pair( alphaMachEnvelopeUB_alpha( 0 ), alphaMachEnvelopeUB_alpha( alphaMachEnvelopeUB_alpha.size() - 1 ) );

    std::pair< double, double > alphaMachEnvelopeLBInterpolatorDomainInterval = std::make_pair( alphaMachEnvelopeLB_Mach.minCoeff(), alphaMachEnvelopeLB_Mach.maxCoeff() );
    std::pair< double, double > alphaMachEnvelopeLBInterpolatorRangeInterval  = std::make_pair( alphaMachEnvelopeLB_alpha.minCoeff(), alphaMachEnvelopeLB_alpha.maxCoeff() );
    std::pair< double, double > alphaMachEnvelopeLBInterpolatorBoundaryValues = std::make_pair( alphaMachEnvelopeLB_alpha( 0 ), alphaMachEnvelopeLB_alpha( alphaMachEnvelopeLB_alpha.size() - 1 ) );

    //! Create alpha-Mach upper bound Interpolator.
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > >  alphaMachEnvelopeUBInterpolator =
            tudat::interpolators::createOneDimensionalInterpolator< double, double >(
                map_alphaMachEnvelopeUB,
                alphaMachEnvelopeInterpolatorSettings,
                alphaMachEnvelopeUBInterpolatorBoundaryValues );

    //! Create alpha-Mach lower bound Interpolator.
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > >  alphaMachEnvelopeLBInterpolator =
            tudat::interpolators::createOneDimensionalInterpolator< double, double >(
                map_alphaMachEnvelopeLB,
                alphaMachEnvelopeInterpolatorSettings,
                alphaMachEnvelopeLBInterpolatorBoundaryValues );

    bislip::Variables::printAlphaMachBounds(
                alphaMachEnvelopeUBInterpolator,
                alphaMachEnvelopeUBInterpolatorDomainInterval,
                alphaMachEnvelopeLBInterpolator,
                alphaMachEnvelopeLBInterpolatorDomainInterval,
                outputPath,
                outputSubFolder );

    //! Set alpha-Mach upper bound Interpolator
    bislipSystems->setAlphaMachEnvelopeUBInterpolator( alphaMachEnvelopeUBInterpolator );

    //! Set alpha-Mach lower bound Interpolator
    bislipSystems->setAlphaMachEnvelopeLBInterpolator( alphaMachEnvelopeLBInterpolator );

}
void printAlphaMachBounds(
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &alphaMachEnvelopeUBInterpolator,
        const std::pair< double, double > &alphaMachEnvelopeUBInterpolatorDomainInterval,
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &alphaMachEnvelopeLBInterpolator,
        const std::pair< double, double > &alphaMachEnvelopeLBInterpolatorDomainInterval,
        const std::string &outputPath,
        const std::string &outputSubFolder )
{
    //! Declare data map to contain vectors of interpolated values.
    std::map< double, Eigen::Vector2d > map_alphaMachBounds;
    Eigen::Vector2d alphaMachBounds;

    //! Loop to populate vectors of interpolated values and then pass to data map.
    //!     Number of evaluations has been arbitrarily selected.
    double pp = 0;
    double domainInterval = alphaMachEnvelopeLBInterpolatorDomainInterval.second - alphaMachEnvelopeLBInterpolatorDomainInterval.first;
    for ( unsigned int i = 0; i < 1001; ++i )
    {
        alphaMachBounds( 0 ) = alphaMachEnvelopeLBInterpolator->interpolate( pp * domainInterval / 1000 );
        alphaMachBounds( 1 ) = alphaMachEnvelopeUBInterpolator->interpolate( pp * domainInterval / 1000 );
        map_alphaMachBounds[ pp * domainInterval / 1000 ] = alphaMachBounds;
        pp += 1;
    }

    //std::cout << "Print out Interpolators and DVs" << std::endl;
    //! Print data maps containing vectors of evaluation of interpolators.
    tudat::input_output::writeDataMapToTextFile( map_alphaMachBounds,
                                                 "alphaMachBounds",
                                                 outputPath + outputSubFolder,
                                                 "",
                                                 std::numeric_limits< double >::digits10,
                                                 std::numeric_limits< double >::digits10,
                                                 "," );
}



void createHeadingErrorDeadBandInterpolator(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::vector< double > &headingErrorDeadBand,
        const std::string &outputPath,
        const std::string &outputSubFolder )
{
    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Coarse interpolator section
    std::map< double, double > map_headingErrorDeadBandInterpolator;

    //! Declare and initialize Interpolator Settings.
    std::shared_ptr< tudat::interpolators::InterpolatorSettings > headingErrorDeadbandInterpolatorSettings = std::make_shared< tudat::interpolators::InterpolatorSettings >( tudat::interpolators::InterpolatorTypes::linear_interpolator );

    //! Declare and initialize size vector size.
    unsigned long deadBandSize = headingErrorDeadBand.size() / 2;

    //! Declare Heading Error Deadband vectors.
    Eigen::VectorXd headingErrorDeadBand_distance( deadBandSize );
    Eigen::VectorXd headingErrorDeadBand_error(deadBandSize );

    //! Populate vectors that will be used with data map.
    for( unsigned long i = 0; i < deadBandSize; i++ )
    {
        headingErrorDeadBand_distance( i ) =  headingErrorDeadBand[ i ];
        headingErrorDeadBand_error( i )    =  headingErrorDeadBand[ i + deadBandSize ];
    }

    //! Populate data map for interpolator.
    for ( unsigned long i = 0; i < deadBandSize; ++i )
    { map_headingErrorDeadBandInterpolator[ headingErrorDeadBand_distance( i ) ] = headingErrorDeadBand_error( i ); }

    // std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > headingErrorDeadBandInterpolator =
    //       tudat::interpolators::createOneDimensionalInterpolator( map_headingErrorDeadBandInterpolator, interpolatorSettings );
    std::pair< double, double > headingErrorDeadBandInterpolatorDomainInterval = std::make_pair( headingErrorDeadBand_distance.minCoeff(), headingErrorDeadBand_distance.maxCoeff() );
    std::pair< double, double > headingErrorDeadBandInterpolatorRangeInterval  = std::make_pair( headingErrorDeadBand_error.minCoeff(), headingErrorDeadBand_error.maxCoeff() );
    std::pair< double, double > headingErrorDeadBandInterpolatorBoundaryValues = std::make_pair( headingErrorDeadBand_error( 0 ), headingErrorDeadBand_error( headingErrorDeadBand_error.size() - 1 ) );

    //! Create Heading Error Deadband Interpolator.
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > >  headingErrorDeadBandInterpolator =
            tudat::interpolators::createOneDimensionalInterpolator< double, double >(
                map_headingErrorDeadBandInterpolator,
                headingErrorDeadbandInterpolatorSettings,
                headingErrorDeadBandInterpolatorBoundaryValues );

    //! Set Heading Error Deadband Interpolator
    bislipSystems->setHeadingErrorDeadBandInterpolator( headingErrorDeadBandInterpolator );

    bislip::Variables::printHeadingErrorDeadBandBounds(
                headingErrorDeadBandInterpolator,
                headingErrorDeadBandInterpolatorDomainInterval,
                outputPath,
                outputSubFolder );

}

void printHeadingErrorDeadBandBounds(
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &headingErrorDeadBandInterpolator,
        const std::pair< double, double > &headingErrorDeadBandInterpolatorDomainInterval,
        const std::string &outputPath,
        const std::string &outputSubFolder )
{
    //! Declare data map to contain vectors of interpolated values.
    std::map< double, Eigen::Vector2d > map_HeadingErrorDeadBandBounds;
    Eigen::Vector2d headingErrorDeadBandBounds;

    //! Loop to populate vectors of interpolated values and then pass to data map.
    //!     Number of evaluations has been arbitrarily selected.
    double pp = 0;
    double domainInterval = headingErrorDeadBandInterpolatorDomainInterval.second - headingErrorDeadBandInterpolatorDomainInterval.first;
    for ( unsigned int i = 0; i < 1001; ++i )
    {
        headingErrorDeadBandBounds( 0 ) = headingErrorDeadBandInterpolator->interpolate( pp * domainInterval / 1000 );
        headingErrorDeadBandBounds( 1 ) = -headingErrorDeadBandBounds( 0 ) ;
        map_HeadingErrorDeadBandBounds[ pp * domainInterval / 1000 ] = headingErrorDeadBandBounds;
        pp += 1;
    }

    //std::cout << "Print out Interpolators and DVs" << std::endl;
    //! Print data maps containing vectors of evaluation of interpolators.
    tudat::input_output::writeDataMapToTextFile( map_HeadingErrorDeadBandBounds,
                                                 "headingErrorDeadBandBounds",
                                                 outputPath + outputSubFolder,
                                                 "",
                                                 std::numeric_limits< double >::digits10,
                                                 std::numeric_limits< double >::digits10,
                                                 "," );
}


void createKourouGuidanceInterpolators(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::vector< double > &kourouAngleofAttackHistory,
        const std::vector< double > &kourouBankAngleHistory,
        const std::string &outputPath,
        const std::string &outputSubFolder )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //!
    std::map< double, double > map_kourouAngleofAttackHistory, map_kourouBankAngleHistory;

    //! Declare and initialize alpha-Mach Interpolator Settings.
    std::shared_ptr< tudat::interpolators::InterpolatorSettings > kourouInterpolatorSettings = std::make_shared< tudat::interpolators::InterpolatorSettings >( tudat::interpolators::InterpolatorTypes::linear_interpolator );

    //! Declare and initialize size vector size.
    unsigned long kourouAngleofAttackHistorySize = kourouAngleofAttackHistory.size() / 2;
    unsigned long kourouBankAngleHistorySize = kourouBankAngleHistory.size() / 2;

    //! Declare alpha-Mach bound vectors.
    Eigen::VectorXd kourouAngleofAttack_alpha( kourouAngleofAttackHistorySize );
    Eigen::VectorXd kourouAngleofAttack_time( kourouAngleofAttackHistorySize );
    Eigen::VectorXd kourouBankAngle_sigma( kourouBankAngleHistorySize );
    Eigen::VectorXd kourouBankAngle_time( kourouBankAngleHistorySize );


    //! Populate vectors that will be used with data map.
    for( unsigned long i = 0; i < kourouAngleofAttackHistorySize; i++ )
    {
        kourouAngleofAttack_time( i )  =  kourouAngleofAttackHistory[ i ];
        kourouAngleofAttack_alpha( i ) =  kourouAngleofAttackHistory[ i + kourouAngleofAttackHistorySize ];
    }

    //! Populate data map for interpolator.
    for ( unsigned long i = 0; i < kourouAngleofAttackHistorySize; ++i )
    { map_kourouAngleofAttackHistory[ kourouAngleofAttack_time( i ) ] = kourouAngleofAttack_alpha( i ); }

    //! Populate vectors that will be used with data map.
    for( unsigned long i = 0; i < kourouBankAngleHistorySize; i++ )
    {
        kourouBankAngle_time( i )  =  kourouBankAngleHistory[ i ];
        kourouBankAngle_sigma( i ) =  kourouBankAngleHistory[ i + kourouBankAngleHistorySize ];
    }

    //! Populate data map for interpolator.
    for ( unsigned long i = 0; i < kourouBankAngleHistorySize; ++i )
    { map_kourouBankAngleHistory[ kourouBankAngle_time( i ) ] = kourouBankAngle_sigma( i ); }


    std::pair< double, double > kourouAngleOfAttackInterpolatorDomainInterval = std::make_pair( kourouAngleofAttack_time.minCoeff(), kourouAngleofAttack_time.maxCoeff() );
    std::pair< double, double > kourouAngleofAttackInterpolatorRangeInterval  = std::make_pair( kourouAngleofAttack_alpha.minCoeff(), kourouAngleofAttack_alpha.maxCoeff() );
    std::pair< double, double > kourouAngleofAttackInterpolatorBoundaryValues = std::make_pair( kourouAngleofAttack_alpha( 0 ), kourouAngleofAttack_alpha( kourouAngleofAttack_alpha.size() - 1 ) );

    std::pair< double, double > kourouBankAngleInterpolatorDomainInterval = std::make_pair( kourouBankAngle_time.minCoeff(), kourouBankAngle_time.maxCoeff() );
    std::pair< double, double > kourouBankAngleInterpolatorRangeInterval  = std::make_pair( kourouBankAngle_sigma.minCoeff(), kourouBankAngle_sigma.maxCoeff() );
    std::pair< double, double > kourouBankAngleInterpolatorBoundaryValues = std::make_pair( kourouBankAngle_sigma( 0 ), kourouBankAngle_sigma( kourouBankAngle_sigma.size() - 1 ) );


    //! Create Heading Error Deadband Interpolator.
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > >  kourouAngleOfAttackInterpolator =
            tudat::interpolators::createOneDimensionalInterpolator< double, double >(
                map_kourouAngleofAttackHistory,
                kourouInterpolatorSettings,
                kourouAngleofAttackInterpolatorBoundaryValues );

    //! Create Heading Error Deadband Interpolator.
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > >  kourouBankAngleInterpolator =
            tudat::interpolators::createOneDimensionalInterpolator< double, double >(
                map_kourouBankAngleHistory,
                kourouInterpolatorSettings,
                kourouBankAngleInterpolatorBoundaryValues );

    bislip::Variables::printKourouGuidanceInterpolators(
                kourouAngleOfAttackInterpolator,
                kourouAngleOfAttackInterpolatorDomainInterval,
                kourouBankAngleInterpolator,
                kourouBankAngleInterpolatorDomainInterval,
                outputPath,
                outputSubFolder );

    //! Set Heading Error Deadband Interpolator
    bislipSystems->setKourouAngleOfAttackInterpolator( kourouAngleOfAttackInterpolator );
    bislipSystems->setKourouBankAngleInterpolator( kourouBankAngleInterpolator );

    //! Turn on flag for validation
    bislipSystems->setValidationFlag( true );

}



void printKourouGuidanceInterpolators(
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &kourouAngleOfAttackInterpolator,
        const std::pair< double, double > &kourouAngleOfAttackInterpolatorDomainInterval,
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &kourouBankAngleInterpolator,
        const std::pair< double, double > &kourouBankAngleInterpolatorDomainInterval,
        const std::string &outputPath,
        const std::string &outputSubFolder )
{
    //! Declare data map to contain vectors of interpolated values.
    std::map< double, Eigen::Vector2d > map_kourouGuidanceEvaluation;
    Eigen::Vector2d KourouGuidance;

    //! Loop to populate vectors of interpolated values and then pass to data map.
    //!     Number of evaluations has been arbitrarily selected.
    double pp = 0;
    double domainInterval = kourouAngleOfAttackInterpolatorDomainInterval.second - kourouAngleOfAttackInterpolatorDomainInterval.first;
    for ( unsigned int i = 0; i < 1001; ++i )
    {
        KourouGuidance( 0 ) = kourouAngleOfAttackInterpolator->interpolate( pp * domainInterval / 1000 );
        KourouGuidance( 1 ) = kourouBankAngleInterpolator->interpolate( pp * domainInterval / 1000 );
        map_kourouGuidanceEvaluation[ pp * domainInterval / 1000 ] = KourouGuidance;
        pp += 1;
    }

    //std::cout << "Print out Interpolators and DVs" << std::endl;
    //! Print data maps containing vectors of evaluation of interpolators.
    tudat::input_output::writeDataMapToTextFile( map_kourouGuidanceEvaluation,
                                                 "kourouGuidanceEvaluation",
                                                 outputPath + outputSubFolder,
                                                 "",
                                                 std::numeric_limits< double >::digits10,
                                                 std::numeric_limits< double >::digits10,
                                                 "," );
}


double computeSumOfEigenVectorXd( const Eigen::VectorXd &vector )
{
    double vectorSum = 0;
    for ( int i = 0; i < vector.size() - 1; i++ )
    { vectorSum += vector( i ); }

    return vectorSum;
}

Eigen::VectorXd computeElementwiseDivisionOfEigenVectorXd(
        const Eigen::VectorXd numerator,
        const Eigen::VectorXd denominator )
{
    Eigen::VectorXd quotient( numerator.size() );
    for ( int i = 0; i <  numerator.size() - 1; i++ )
    {
        quotient( i ) = numerator( i ) / denominator( i );
        if ( std::isnan( quotient( i ) ) == true ) { quotient( i ) = 0.0;}
    }

    return quotient;
}

double computeHeatingRateChapman(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    //std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
    //          bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    //! Extract Vehicle Systems pointer.
    //std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( ) ;

    /*
    const double currentAirDensity = flightConditions->getCurrentDensity( );
    const double currentAirspeed = flightConditions->getCurrentAirspeed( );
    const double airTemperature = flightConditions->getCurrentFreestreamTemperature( );
    const double machNumber = flightConditions->getCurrentMachNumber( );
    const double wallEmissivity = vehicleSystems->getWallEmissivity( );
    const double noseRadius = vehicleSystems->getNoseRadius();

    // Compute adiabatic wall temperature.
    double adiabaticWallTemperature
            = tudat::aerodynamics::computeAdiabaticWallTemperature( airTemperature , machNumber );
*/
    // bislipSystems->setWorkingRadius( bislipSystems->getNoseRadius() );
    if( debugInfo == 1 ){ std::cout << "          Chapman Heat Flux" << std::endl; }
    const double chapmanHeatFlux = bislip::Variables::computeStagnationHeatFlux( bodyMap, vehicleName);
    // double stagnationConvectiveHeatFlux =  C1 * ( currentAirDensity / std::sqrt( noseRadius ) ) * std::pow( currentAirspeed, C2 ) * ( 1.0 - wallTemperature / adiabaticWallTemperature );


    return chapmanHeatFlux;
}

double computeHeatingRate (
        const double &currentAirdensity,
        const double &currentAirspeed,
        const double &C,
        const double &N,
        const double &M)
{
    /*  std::cout << "computeHeatingRate----------------------------------------------------------" << std::endl;

    std::cout << "currentAirdensity: " << currentAirdensity << std::endl;
    std::cout << "currentAirspeed: " << currentAirspeed << std::endl;
    std::cout << "C: " << C << std::endl;
    std::cout << "N: " << N << std::endl;
    std::cout << "M: " << M << std::endl;
    std::cout << "C * std::pow( currentAirdensity, N ) * std::pow( currentAirspeed, M ): " << C * std::pow( currentAirdensity, N ) * std::pow( currentAirspeed, M ) << std::endl;
*/
    return C * std::pow( currentAirdensity, N ) * std::pow( currentAirspeed, M );
}

double computeStagnationConvectiveHeatFlux (
        const double &currentAirdensity,
        const double &currentAirspeed,
        const double &C_s,
        const double &N,
        const double &M,
        const double &adiabaticWallTemperature,
        const double &currentWallTemperature)
{
    double C = C_s * ( 1 - ( currentWallTemperature / adiabaticWallTemperature ) ) ;
    //double stagnationConvectiveHeatFlux = bislip::Variables::computeHeatingRate ( currentAirdensity, currentAirspeed, C, N, M) ;

    /* std::cout << "computeStagnationConvectiveHeatFlux--------------------------------------------" << std::endl;
    std::cout << "currentAirdensity: " << currentAirdensity << std::endl;
    std::cout << "currentAirspeed: " << currentAirspeed << std::endl;
    std::cout << "currentWallTemperature: " << currentWallTemperature << std::endl;
    std::cout << "adiabaticWallTemperature: " << adiabaticWallTemperature << std::endl;
    std::cout << "C: " << C << std::endl;
    std::cout << "C_s: " << C_s << std::endl;
    std::cout << "N: " << N << std::endl;
    std::cout << "M: " << M << std::endl;
    std::cout << "stagnationConvectiveHeatFlux: " << stagnationConvectiveHeatFlux << std::endl;
*/
    return bislip::Variables::computeHeatingRate( currentAirdensity, currentAirspeed, C, N, M );
}

double computeStagnationHeatFlux (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    /*
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Extract Vehicle Systems pointer.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "Extract current conditions" << std::endl; }
    const double currentAirdensity     = flightConditions->getCurrentDensity( );
    const double currentAirspeed       = flightConditions->getCurrentAirspeed( );
    const double currentAirTemperature = flightConditions->getCurrentFreestreamTemperature( );
    const double currentMachNumber     = flightConditions->getCurrentMachNumber( );
    if( debugInfo == 1 ){ std::cout << "     currentAirdensity     = " << currentAirdensity << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentAirspeed       = " << currentAirspeed << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentAirTemperature = " << currentAirTemperature << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentMachNumber     = " << currentMachNumber << std::endl; }

    if( debugInfo == 1 ){ std::cout << "Set Curvature Radius" << std::endl; }
    double curvatureRadius = bislipSystems->getWorkingRadius( );

    if( debugInfo == 1 ){ std::cout << "Extract wall emissivity" << std::endl; }
    const double wallEmissivity = vehicleSystems->getWallEmissivity( );


    const double M = 3.0;
    const double N = 0.5;

    const double C_s = std::pow( curvatureRadius, -0.5 ) * std::pow( 1.225, -N ) * std::pow( 7905, -M ) * ( 1.06584E9 );
    //const double C_s = std::pow( curvatureRadius, -0.5 ) * ( 1.83 ) * ( 1E-8 ) * ( 100 * 100 );


    if( debugInfo == 1 ){ std::cout << "     C_s = " << C_s << std::endl; }



    //! Compute adiabatic wall temperature.
    if( debugInfo == 1 ){ std::cout << "Compute adiabatic wall temperature." << std::endl; }
    double adiabaticWallTemperature =
            tudat::aerodynamics::computeAdiabaticWallTemperature( currentAirTemperature, currentMachNumber );
    if( debugInfo == 1 ){ std::cout << "     adiabaticWallTemperature = " << adiabaticWallTemperature << std::endl; }

    if( debugInfo == 1 ){ std::cout << "Set Heat Transfer Function." << std::endl; }
    std::function< double( const double ) > heatTransferFunction =
            std::bind( &bislip::Variables::computeStagnationConvectiveHeatFlux, currentAirdensity, currentAirspeed, C_s, N, M, adiabaticWallTemperature, std::placeholders::_1 );
  */
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Extract Vehicle Systems pointer.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "Extract wall emissivity" << std::endl; }
    const double wallEmissivity = vehicleSystems->getWallEmissivity( );

    if( debugInfo == 1 ){ std::cout << "Extract current conditions" << std::endl; }
    const double currentAirTemperature = flightConditions->getCurrentFreestreamTemperature( );
    const double currentMachNumber     = flightConditions->getCurrentMachNumber( );
    if( debugInfo == 1 ){ std::cout << "     currentAirTemperature = " << currentAirTemperature << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentMachNumber     = " << currentMachNumber << std::endl; }

    //! Compute adiabatic wall temperature.
    if( debugInfo == 1 ){ std::cout << "Compute adiabatic wall temperature." << std::endl; }
    double adiabaticWallTemperature =
            bislip::Variables::computeAdiabaticWallTemperature( currentAirTemperature, currentMachNumber );
    if( debugInfo == 1 ){ std::cout << "     adiabaticWallTemperature = " << adiabaticWallTemperature << std::endl; }



    if( debugInfo == 1 ){ std::cout << "Get heat transfer function" << std::endl; }
    std::function< double( const double ) > heatTransferFunction = bislip::Variables::getStagnationHeatTransferFunction( bodyMap, vehicleName );

    if( debugInfo == 1 ){ std::cout << "Find equilibrium wall temperature" << std::endl; }
    double equilibriumWallTemperature = bislip::Variables::findEquilibriumWallTemperature( heatTransferFunction, wallEmissivity, adiabaticWallTemperature);


    // std::function< double( const double ) > equilibriumWallTemperatureRootFindingFunction =
    //       std::bind( &bislip::Variables::computeEquilibiumWallTemperatureRootFinder, heatTransferFunction, wallEmissivity, std::placeholders::_1  );

    //   double root = rootFinderBisection( equilibriumWallTemperatureRootFindingFunction, 0.0, adiabaticWallTemperature );
    const double currentAirdensity     = flightConditions->getCurrentDensity( );
    const double currentAirspeed       = flightConditions->getCurrentAirspeed( );


    /*
    double a = 0.0;
    double b = adiabaticWallTemperature;
    double f_a, f_root;
    if ( equilibriumWallTemperatureRootFindingFunction( a ) * equilibriumWallTemperatureRootFindingFunction( b ) < 0.0 )
    {
        root = ( a + b ) / 2.0;

        if( debugInfo == 1 ){ std::cout << "     Starting Stagnation Eq. wall temp. search:" << std::endl; }
        while ( std::abs( equilibriumWallTemperatureRootFindingFunction( root ) ) > 1e-7 )
        {

            f_a       = equilibriumWallTemperatureRootFindingFunction( a );
            f_root    = equilibriumWallTemperatureRootFindingFunction( root );

            if ( f_a * f_root < 0.0 ) { b = root; }
            else { a = root; }
            root = ( a + b ) / 2.0;
            if( debugInfo == 1 ){ std::cout << "     root = " << root << "     |     " << "equilibriumWallTemperatureRootFindingFunction( root )  = " << std::abs( equilibriumWallTemperatureRootFindingFunction( root ) ) << std::endl; }
        }
    }


    */
    bislipSystems->setWallTemperature( equilibriumWallTemperature );

    double heatflux = heatTransferFunction( equilibriumWallTemperature );



    if ( bislipSystems->getWorkingRadius() == 0.8) {
        if( debugInfo == 40 ){ std::cout << heatflux << "," << equilibriumWallTemperature << "," << adiabaticWallTemperature << "," << currentAirdensity << "," << currentAirspeed << std::endl; }
    }

    return heatflux;
}


double findEquilibriumWallTemperature(
        const std::function< double( const double ) > &heatTransferFunction,
        const double &wallEmissivity,
        const double &adiabaticWallTemperature )
{
    std::function< double( const double ) > equilibriumWallTemperatureRootFindingFunction =
            std::bind( &bislip::Variables::computeEquilibiumWallTemperatureRootFinder, heatTransferFunction, wallEmissivity, std::placeholders::_1  );

    return bislip::Variables::rootFinderBisection( equilibriumWallTemperatureRootFindingFunction, 0.0, adiabaticWallTemperature );
}





std::function< double( const double ) > getStagnationHeatTransferFunction(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{

    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "Extract current conditions" << std::endl; }
    const double currentAirdensity     = flightConditions->getCurrentDensity( );
    const double currentAirspeed       = flightConditions->getCurrentAirspeed( );
    const double currentAirTemperature = flightConditions->getCurrentFreestreamTemperature( );
    const double currentMachNumber     = flightConditions->getCurrentMachNumber( );
    if( debugInfo == 1 ){ std::cout << "     currentAirdensity     = " << currentAirdensity << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentAirspeed       = " << currentAirspeed << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentAirTemperature = " << currentAirTemperature << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentMachNumber     = " << currentMachNumber << std::endl; }

    if( debugInfo == 1 ){ std::cout << "Set Curvature Radius" << std::endl; }
    double curvatureRadius = bislipSystems->getWorkingRadius( );

    const double M = 3.0;
    const double N = 0.5;

    const double C_s = std::pow( curvatureRadius, -0.5 ) * std::pow( 1.225, -N ) * std::pow( 7905, -M ) * ( 1.06584E8 );

    if( debugInfo == 1 ){ std::cout << "     C_s = " << C_s << std::endl; }

    //! Compute adiabatic wall temperature.
    if( debugInfo == 1 ){ std::cout << "Compute adiabatic wall temperature." << std::endl; }
    double adiabaticWallTemperature =
            bislip::Variables::computeAdiabaticWallTemperature( currentAirTemperature, currentMachNumber );
    if( debugInfo == 1 ){ std::cout << "     adiabaticWallTemperature = " << adiabaticWallTemperature << std::endl; }


    if( debugInfo == 1 ){ std::cout << "Set Heat Transfer Function." << std::endl; }
    std::function< double( const double ) > heatTransferFunction =
            std::bind( &bislip::Variables::computeStagnationConvectiveHeatFlux, currentAirdensity, currentAirspeed, C_s, N, M, adiabaticWallTemperature, std::placeholders::_1 );

    return heatTransferFunction;
}


//! Compute the adiabatic wall temperature experienced by a vehicle.
double computeAdiabaticWallTemperature(
        const double airTemperature, const double machNumber )
{

    double ratioSpecificHeats = 1.4;
    double recoveryFactor = 0.85;

    double totalTemperature
            = airTemperature * ( 1 + 0.5 * ( ratioSpecificHeats - 1 ) * machNumber * machNumber );

    return airTemperature + recoveryFactor * ( totalTemperature - airTemperature );
}







//! Funtion to compute the equilibrium heat flux experienced by a vehicle
/*double computeEquilibriumHeatflux(
        const std::function< double( const double ) > heatTransferFunction,
        const double wallEmmisivity,
        const double adiabaticWallTemperature,
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    return heatTransferFunction( bislip::Variables::computeEquilibiumWallTemperature( heatTransferFunction, wallEmmisivity, adiabaticWallTemperature, bodyMap, vehicleName ) );
}
*/

double computeFlatPlateHeat (
        const double &currentAirdensity,
        const double &currentAirspeed,
        const double &C_turb_FP_1,
        const double &C_turb_FP_2,
        const double &adiabaticWallTemperature,
        const double &currentWallTemperature)
{
    double C = 0.0;
    double M = 0.0;
    double N = 0.0;


    //! Turbulent Boundary Layer
    if ( currentAirspeed <= 3962.0 )
    {
        M = 3.37;
        N = 0.8;
        C =  C_turb_FP_1 * std::pow( 556 / currentWallTemperature, 1.0 / 4.0 ) * ( 1.0 - 1.11 * ( currentWallTemperature / adiabaticWallTemperature ) );
    }
    else if ( currentAirspeed > 3962.0 )
    {
        M = 3.7;
        N = 0.8;
        C = C_turb_FP_2 * ( 1 - 1.11 * ( currentWallTemperature / adiabaticWallTemperature ) );
    }

    return bislip::Variables::computeHeatingRate ( currentAirdensity, currentAirspeed, C, N, M);;
}

double computeFlatPlateHeatFlux (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Extract Vehicle Systems pointer.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( );

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "Calculating Flat Plate Heat Flux" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Extract current flight conditions" << std::endl; }
    const double currentAirdensity     = flightConditions->getCurrentDensity( );
    const double currentAirspeed       = flightConditions->getCurrentAirspeed( );
    const double currentAirTemperature = flightConditions->getCurrentFreestreamTemperature( );
    const double currentMachNumber     = flightConditions->getCurrentMachNumber( );
    if( debugInfo == 1 ){ std::cout << "     currentAirdensity     = " << currentAirdensity << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentAirspeed       = " << currentAirspeed << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentAirTemperature = " << currentAirTemperature << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentMachNumber     = " << currentMachNumber << std::endl; }


    if( debugInfo == 1 ){ std::cout << "     Extract vehicle parameters" << std::endl; }
    const double wallEmissivity = vehicleSystems->getWallEmissivity( );
    const double phi = bislipSystems->getLocalBodyAngle( );
    const double x_T = bislipSystems->getTransitionDistance( );

    const double sin_phi_term = std::pow( std::sin( phi ), 1.6 ) ;
    const double cos_phi = std::cos( phi );
    const double x_T_term = std::pow( x_T, -1.0 / 5.0 ) ;
    const double C_turb_FP_1 = ( 3.35 ) * sin_phi_term * std::pow( cos_phi, 1.78 ) * x_T_term * 1E-8;
    const double C_turb_FP_2 = ( 2.2 ) * sin_phi_term * std::pow( cos_phi, 2.08 ) * x_T_term * 1E-9;
    if( debugInfo == 1 ){ std::cout << "     C_turb_FP_1 = " << C_turb_FP_1 << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     C_turb_FP_2 = " << C_turb_FP_2 << std::endl; }




    //! Compute adiabatic wall temperature.
    if( debugInfo == 1 ){ std::cout << "     Compute adiabatic wall temperature" << std::endl; }
    double adiabaticWallTemperature
            = tudat::aerodynamics::computeAdiabaticWallTemperature( currentAirTemperature , currentMachNumber );
    if( debugInfo == 1 ){ std::cout << "     adiabaticWallTemperature = " << adiabaticWallTemperature << std::endl; }

    if( debugInfo == 1 ){ std::cout << "     Set heat transfer function" << std::endl; }
    std::function< double( const double ) > heatTransferFunction =
            std::bind( &bislip::Variables::computeFlatPlateHeat, currentAirdensity, currentAirspeed, C_turb_FP_1, C_turb_FP_2, adiabaticWallTemperature, std::placeholders::_1 );

    std::function< double( const double ) > equilibriumWallTemperatureRootFindingFunction =
            std::bind( &bislip::Variables::computeEquilibiumWallTemperatureRootFinder, heatTransferFunction, wallEmissivity, std::placeholders::_1  );

    double a = 0.0; equilibriumWallTemperatureRootFindingFunction( 0.0 );
    double b = adiabaticWallTemperature;
    double root = 0.0;
    double f_a, f_root;
    if ( equilibriumWallTemperatureRootFindingFunction( a ) * equilibriumWallTemperatureRootFindingFunction( b ) < 0.0 )
    {
        root = ( a + b ) / 2.0;

        if( debugInfo == 1 ){ std::cout << "     Starting Flat Plate Eq. wall temp. search:" << std::endl; }
        while ( std::abs( equilibriumWallTemperatureRootFindingFunction( root ) ) > 1e-7 )
        {

            f_a       = equilibriumWallTemperatureRootFindingFunction( a );
            f_root    = equilibriumWallTemperatureRootFindingFunction( root );

            if ( f_a * f_root < 0.0 ) { b = root; }
            else { a = root; }
            root = ( a + b ) / 2.0;
            if( debugInfo == 1 ){ std::cout << "     root = " << root << "     |     " << "equilibriumWallTemperatureRootFindingFunction( root )  = " << std::abs( equilibriumWallTemperatureRootFindingFunction( root ) ) << std::endl; }
        }
    }

    bislipSystems->setWallTemperature( root );

    return heatTransferFunction( root );


    // if( debugInfo == 1 ){ std::cout << "     Determine equilibrium flat plate heat flux" << std::endl; }
    // double flatPlateHeatFlux = bislip::Variables::computeEquilibriumHeatflux( heatTransferFunction, wallEmissivity, adiabaticWallTemperature, bodyMap, vehicleName );

    //return flatPlateHeatFlux;
}

double computeHeatingRateTauber(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "Calculating Tauber Heat Flux for Leading Edges" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "     Calculate Stagnation Heat Flux Component of Tauber Heat Flux" << std::endl; }

    bislipSystems->setWallTemperature( bislipSystems->getTauberWallTempStagnation() );
    //if( debugInfo == 1 ){ std::cout << "     Initial Wall Temperature = " <<  bislipSystems->getTauberWallTempStagnation( ) << std::endl; }

    const double q_dot_s = bislip::Variables::computeStagnationHeatFlux( bodyMap, vehicleName);

    if( debugInfo == 1 ){ std::cout << "     Set Stagnation Heat Flux Component of Tauber Heat Flux" << std::endl; }
    bislipSystems->setTauberHeatFluxStagnation( q_dot_s );
    bislipSystems->setTauberWallTempStagnation( bislipSystems->getWallTemperature() );
    if( debugInfo == 1 ){ std::cout << "          Current Tauber Stagnation Heat Flux = " <<  bislipSystems->getTauberHeatFluxStagnation( ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "          Current Tauber Stagnation Eq. Temp. = " <<  bislipSystems->getTauberWallTempStagnation( ) << std::endl; }

    if( debugInfo == 1 ){ std::cout << "     Calculate Flat Plate Heat Flux Component of Tauber Heat Flux" << std::endl; }
    bislipSystems->setWallTemperature( bislipSystems->getTauberWallTempFlatPlate() );
    //if( debugInfo == 1 ){ std::cout << "     Initial Wall Temperature = " <<  bislipSystems->getTauberWallTempFlatPlate( ) << std::endl; }
    const double q_dot_FP = ( 100 * 100 ) * bislip::Variables::computeFlatPlateHeatFlux( bodyMap, vehicleName);

    if( debugInfo == 1 ){ std::cout << "     Set Flat Plate Heat Flux Component of Tauber Heat Flux" << std::endl; }
    bislipSystems->setTauberHeatFluxFlatPlate( q_dot_FP );
    bislipSystems->setTauberWallTempFlatPlate( bislipSystems->getWallTemperature() );
    if( debugInfo == 1 ){ std::cout << "          Current Tauber Flat Plate Heat Flux = " <<  bislipSystems->getTauberHeatFluxFlatPlate( ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "          Current Tauber Flat Plate Eq. Temp. = " <<  bislipSystems->getTauberWallTempFlatPlate( ) << std::endl; }



    //const double lambda = bislipSystems->getWingSweepAngle();
    const double sin_lambda = std::sin( bislipSystems->getWingSweepAngle() );
    const double cos_lambda = std::cos( bislipSystems->getWingSweepAngle() );

    if( debugInfo == 1 ){ std::cout << "     Calculate Tauber Heat Flux Component of Tauber Heat Flux" << std::endl; }
    double leadingEdgeHeatFlux = std::sqrt( 0.5 * q_dot_s * q_dot_s * cos_lambda * cos_lambda + q_dot_FP * q_dot_FP * sin_lambda * sin_lambda );

    if( debugInfo == 1 ){ std::cout << "     Set Tauber Heat Flux" << std::endl; }
    bislipSystems->setCurrentHeatFluxTauber( leadingEdgeHeatFlux );

    return leadingEdgeHeatFlux;
}

double computeRadiativeHeatFlux(
        const double &wallEmissivity,
        const double &bodyTemperature )
{
    return wallEmissivity * tudat::physical_constants::STEFAN_BOLTZMANN_CONSTANT *
            bodyTemperature * bodyTemperature * bodyTemperature * bodyTemperature;
}

//! Function to compute the equilibrium wall temperature from the heat input and emmisivity
double computeEquilibiumWallTemperatureRootFinder(
        const std::function< double( const double ) > &heatTransferFunction,
        const double &wallEmmisivity,
        const double &wallTemperature )

{

    const double evaluation = heatTransferFunction( wallTemperature ) - bislip::Variables::computeRadiativeHeatFlux( wallEmmisivity, wallTemperature );

    return evaluation;
}

/*
double computeEquilibiumWallTemperature(
        const std::function< double( const double ) > heatTransferFunction,
        const double wallEmmisivity,
        const double adiabaticWallTemperature,
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "Create the object that contains the function who's root needs to be found" << std::endl; }

    //! Create the object that contains the function who's root needs to be found.
    std::shared_ptr< bislip::Variables::EquilibriumTemperatureFunction > equilibriumTemperatureFunction
            = std::make_shared< bislip::Variables::EquilibriumTemperatureFunction >(
                heatTransferFunction, wallEmmisivity, adiabaticWallTemperature, bislipSystems );

    // Compute wall temperature, first try secant method, use bisection as backup if secany unsuccesfull.
    double wallTemperature = TUDAT_NAN;
   *//*
    try
    {
        if( debugInfo == 1 ){ std::cout << "Try Secant Root Finder" << std::endl; }

        tudat::root_finders::SecantRootFinder::TerminationFunction terminationConditionFunction =
                std::bind( &tudat::root_finders::termination_conditions::RootRelativeToleranceTerminationCondition< double >::checkTerminationCondition,
                           std::make_shared< tudat::root_finders::termination_conditions::RootRelativeToleranceTerminationCondition< double > >(),
                           std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );
        tudat::root_finders::SecantRootFinder secant( terminationConditionFunction );
        wallTemperature = secant.execute(
                    equilibriumTemperatureFunction, equilibriumTemperatureFunction->getInitialGuess( ) );
    }
    catch ( std::runtime_error )
    {
    *//*
        try
        {
            if( debugInfo == 1 ){ std::cout << "Try Bisection Root Finder" << std::endl; }

            tudat::root_finders::Bisection::TerminationFunction terminationConditionFunction =
                    std::bind( &tudat::root_finders::termination_conditions::RootRelativeToleranceTerminationCondition< double >::checkTerminationCondition,
                               std::make_shared< tudat::root_finders::termination_conditions::RootRelativeToleranceTerminationCondition< double > >( 1.0e-12, 2000 ),
                               std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5 );
            tudat::root_finders::Bisection bisection( terminationConditionFunction );
            wallTemperature = bisection.execute(
                        equilibriumTemperatureFunction, equilibriumTemperatureFunction->getInitialGuess( ) );
        }
        catch ( std::runtime_error )
        {
            //wallTemperature = bislipSystems->getWallTemperature();
            //std::cout << "could not find equilibrium wall temperature, using previous value." << std::endl;


            throw std::runtime_error( "Error, could not find equilibrium wall temperature....mine" );
        }

   // }
    bislipSystems->setWallTemperature( wallTemperature );
    if( debugInfo == 1 ){ std::cout << "Root Finder Eq. Temperature = " << bislipSystems->getWallTemperature() << std::endl; }

    return wallTemperature;
}
                   */

//! Function to end the ascent phase once the both the flight-path angle and its rate are both negative.
bool customTermination_FlightPathAngleCombination_Ascent(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName)

{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Extract current conditions.
    const double currentLatitude  = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double currentLongitude = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );
    const double currentFlightPathAngle     = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
    const double currentFlightPathAngleRate = bislipSystems->getCurrentFlightPathAngleRate();
    const double initialDistanceToTarget = bislipSystems->getInitialDistanceToTarget();


    //! Calculate angular distance to go.
    const double angularDistanceTravelled = bislip::Variables::computeAngularDistance (
                bodyMap.at( vehicleName )->getBislipSystems()->getInitialLat(),
                bodyMap.at( vehicleName )->getBislipSystems()->getInitialLon(),
                currentLatitude,
                currentLongitude );

    bool stop = false;

    if ( currentFlightPathAngle < 0.0 )
    {
        if ( currentFlightPathAngleRate * currentFlightPathAngle > 0.0 )
        {
            if (  angularDistanceTravelled / bislipSystems->getInitialDistanceToTarget() > bislipSystems->getAscentTerminationDistanceRatio() ) { stop = true; }
        }
    }
    return stop;
}


bool customTermination_FlightPathAngleCombination_Descent(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName)

{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    const double currentFlightPathAngle     = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
    const double currentFlightPathAngleRate = bislipSystems->getCurrentFlightPathAngleRate();

    bool stop = false;

    if ( currentFlightPathAngle > 0.0 )
    {
        if ( currentFlightPathAngleRate * currentFlightPathAngle > 0.0 ) { stop = true; }
    }

    return stop;
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
        if ( startIterator == 0 )
        {
            //    std::cout << "dependentVariable_Violation.size(): " << dependentVariable_Violation.size() << std::endl;
            //    std::cout << "dependentVariable_Violation( " << 0 << " ): " << dependentVariable_Violation( 0 ) << std::endl;

            for ( long i = startIterator + 1; i < endIterator; i++ )
            {
                if ( dependentVariable_TimeHistory( i ) > dependentVariable_TimeHistory( i - 1 ) ) { dependentVariable_Violation( i ) = dependentVariable_TimeHistory( i ) - dependentVariable_TimeHistory( i - 1 ); }
                //        std::cout << "dependentVariable_Violation( " << i << " ): " << dependentVariable_Violation( i ) << std::endl;
            }
            // std::cout << "dependentVariable_Violation: " << dependentVariable_Violation << std::endl;
            // std::cout << "dependentVariable_TimeHistory: " << dependentVariable_TimeHistory << std::endl;
            // std::cout << "dependentVariable_TimeHistory.size(): " << dependentVariable_TimeHistory.size() << std::endl;

        }
        else if ( startIterator == 1 )
        {
            for ( long i = startIterator; i < endIterator + 1; i++ )
            {
                if ( dependentVariable_TimeHistory( i ) < dependentVariable_TimeHistory( i - 1 ) ) { dependentVariable_Violation( i ) = dependentVariable_TimeHistory( i - 1 ) - dependentVariable_TimeHistory( i ); }
                //  std::cout << "dependentVariable_Violation( " << i << " ): " << dependentVariable_Violation( i ) << std::endl;
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

    double dependentVariable_ViolationSum = bislip::Variables::computeSumOfEigenVectorXd( dependentVariable_Violation );

    if ( direct == true ) { penalty = dependentVariable_ViolationSum; }
    else { penalty =  ( maximum_Violation / constraint ) + ( fixedStepSize * dependentVariable_ViolationSum ) / ( tof * constraint ); }
    //std::cout << "penalty = " << penalty << std::endl;
    return penalty;
}


double computeMonotonicPenalty (
        const Eigen::VectorXd &depVar_TimeHistory,
        const std::string &category )
{
    Eigen::VectorXd dependentVariable_Violation( depVar_TimeHistory.size() );
    dependentVariable_Violation = Eigen::VectorXd::Zero( depVar_TimeHistory.size() );
    // std::cout << "                                      category: " << category << std::endl;
    // std::cout << "                                      depVar_TimeHistory.size() = " << depVar_TimeHistory.size() << std::endl;

    double monotonicPenalty = 0.0;
    if ( category == "Increasing" )
    { for ( long i = 1; i < depVar_TimeHistory.size(); i++ ) { if ( depVar_TimeHistory( i ) < depVar_TimeHistory( i - 1 ) ) { monotonicPenalty += depVar_TimeHistory( i - 1 ) - depVar_TimeHistory( i ); } } }
    else if ( category == "Decreasing" )
    {

        for ( long i = 1; i < depVar_TimeHistory.size(); i++ ) {


            if ( depVar_TimeHistory( i ) > depVar_TimeHistory( i - 1 ) ) {

                monotonicPenalty += depVar_TimeHistory( i ) - depVar_TimeHistory( i - 1 );

                // std::cout << "                                      depVar_TimeHistory( "<< i << " ) - depVar_TimeHistory( " << i - 1 <<" ) = " << depVar_TimeHistory( i ) << " - " << depVar_TimeHistory( i - 1 ) << " = " << depVar_TimeHistory( i ) - depVar_TimeHistory( i - 1 ) <<std::endl;
                // std::cout << "                                                                                             monotonicPenalty = " << monotonicPenalty << std::endl;


            } } }

    return monotonicPenalty;
}

double computeCompoundViolationPenalty(
        const Eigen::VectorXd &depVar_TimeHistory,
        const double &constraint,
        const double &propagationStepSize,
        const double &normalizer )
{
/*
    std::cout << "size of vector = " << depVar_TimeHistory.size() <<std::endl;

    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "depVar_TimeHistory = " << depVar_TimeHistory <<std::endl;


    if (depVar_TimeHistory.hasNaN() == true)
    {
        std::cout << "FUCCK!!!!" <<std::endl;

    }
*/

    double maximumValue = bislip::Variables::determineMaximumValueFromEigenVectorXd( depVar_TimeHistory );



  //  std::cout << "maximumValue = " << maximumValue <<std::endl;
  //  std::cout << "constraint   = " << constraint <<std::endl;

    double compoundViolationPenalty = 0.0;

    if( maximumValue > constraint )
    {

        //std::cout << "----------------------------------------------------" << std::endl;

        Eigen::VectorXd constraintViolations = bislip::Variables::computeConstraintViolations( depVar_TimeHistory, constraint );

        double constraintViolationPenalty = bislip::Variables::computeSumOfEigenVectorXd( constraintViolations ) * ( propagationStepSize / normalizer );

    //    std::cout << "constraintViolationPenalty = " << constraintViolationPenalty <<std::endl;


        double maximumConstraintViolationPenalty = ( maximumValue - constraint ) / constraint;


      //  std::cout << "maximumConstraintViolationPenalty = " << maximumConstraintViolationPenalty <<std::endl;


        //  std::cout << "|||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;

        // std::cout << "    startIterator:  " <<  startIterator << std::endl;
        //std::cout << "    endIterator:  " <<  endIterator << std::endl;
        /* std::cout << "    constraint:  " <<  constraint << std::endl;
        std::cout << "    propagationStepSize:  " <<  propagationStepSize << std::endl;
        std::cout << "    normalizer:  " <<  normalizer << std::endl;
        std::cout << "    maximumValue:  " <<  maximumValue << std::endl;
        std::cout << "    constraint:  " <<  constraint << std::endl;
        std::cout << "    constraintViolationPenalty:  " <<  constraintViolationPenalty << std::endl;
        std::cout << "    maximumConstraintViolationPenalty:  " <<  maximumConstraintViolationPenalty << std::endl;
*/


        compoundViolationPenalty = maximumConstraintViolationPenalty + constraintViolationPenalty;
    //    std::cout << "compoundViolationPenalty = " << compoundViolationPenalty <<std::endl;


    }

   // std::cout << "----------------------------------------------------" << std::endl;
    return compoundViolationPenalty;
}


double computeConstraintViolationPenalty(
        const Eigen::VectorXd &depVar_TimeHistory,
        const double &constraint,
        const double &propagationStepSize,
        const double &normalizer )
{
    return bislip::Variables::computeSumOfEigenVectorXd( bislip::Variables::computeConstraintViolations( depVar_TimeHistory, constraint ) ) * ( propagationStepSize / normalizer );
}


//! Function to calculate the violation magnitudes for a particular trajectory hase.
Eigen::VectorXd computeConstraintViolations(
        const Eigen::VectorXd &depVar_TimeHistory,
        const double &constraint )
{

    Eigen::VectorXd depVar_Violation( depVar_TimeHistory.size() );
    depVar_Violation = Eigen::VectorXd::Zero( depVar_TimeHistory.size() );
    //std::cout << "    constraint:  " << constraint << std::endl;

    //! Calculate the violation magnitude
    for ( long i = 0; i < depVar_TimeHistory.size(); i++ )
    {
        if ( depVar_TimeHistory ( i ) > constraint )
        { depVar_Violation( i ) = depVar_TimeHistory( i ) - constraint;


        }
        // std::cout << "    depVar_Violation( " << i  << " ):  " <<  depVar_Violation( i ) << std::endl;

    }

    //! Zero out all other components of vector (different trajectory phase).
    // if ( startIterator == 0 )
    // { for ( long i = endIterator; i < depVar_TimeHistory.size(); i++ ) { depVar_Violation( i ) = 0; } }
    // else
    // { for ( long i = 0; i < endIterator; i++ ) { depVar_Violation( i ) = 0; } }

    return depVar_Violation;
}

double computeDeadbandViolationPenalty(
        const Eigen::VectorXd &depVar_TimeHistory,
        const Eigen::VectorXd &depVar_Constraint,
        const double &propagationStepSize,
        const double &normalizer )
{
    //! Extract Bislip Systems pointer.
    //  std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //int debugInfo = bislipSystems->getDebugInfo();

    // Eigen::VectorXd interpolatorEvaluation( depVar_TimeHistory.size() );
    // interpolatorEvaluation = Eigen::VectorXd::Zero(depVar_TimeHistory.size() );
    Eigen::VectorXd deadbandViolation = depVar_TimeHistory.cwiseAbs() - depVar_Constraint;
    //deadbandViolation = Eigen::VectorXd::Zero(depVar_TimeHistory.size() );

    //if( debugInfo == 1 ){ std::cout << "Calculating Heading Error DeadBand Violation Penalty" << std::endl; }

    for ( long i = 0; i < depVar_TimeHistory.size(); i++ )
    {
        //        interpolatorEvaluation( i ) =
        //              bodyMap.at( vehicleName )->getBislipSystems()->getHeadingErrorDeadBandInterpolator()->interpolate( depVar_InterpolatorArgument( i ) );

        //    deadbandViolation( i ) = std::abs( depVar_TimeHistory( i ) ) - interpolatorEvaluation( i );

        if ( deadbandViolation( i ) < 0 ) { deadbandViolation( i ) = 0; }
    }



    return bislip::Variables::computeSumOfEigenVectorXd( deadbandViolation ) * ( propagationStepSize / normalizer );

}

//! Function to compute the maximum constrain violation from a vector of violations
double determineMaximumValueFromEigenVectorXd(
        const Eigen::VectorXd &vectorXd )
{
    /*
     * std::cout << "           Finding max value of vector" << std::endl;

    if (vectorXd.hasNaN() == true)
    {
        std::cout << "           FUUUUUUUUUUUUUUCK" << std::endl;


    }

    */
    //std::ptrdiff_t index;
    int i;
    double maximumValue = vectorXd.maxCoeff( &i );
    //std::cout << "           Max value of vector has been found = " << maximumValue << std::endl;

    return maximumValue;
}


}; //namespace Variables

                 } // namespace bislip










