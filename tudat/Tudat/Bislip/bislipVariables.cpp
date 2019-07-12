#include <Tudat/Bislip/bislipVariables.h>
#include <Tudat/Astrodynamics/Propagators/bodyMassStateDerivative.h>

namespace bislip { namespace Variables {

/*
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
*/

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



//! https://thispointer.com/how-to-remove-substrings-from-a-string-in-c/
void eraseAllSubStr(
        std::string & mainStr,
        const std::string & toErase)
{
    size_t pos = std::string::npos;

    // Search for the substring in string in a loop untill nothing is found
    while ((pos  = mainStr.find(toErase) )!= std::string::npos)
    {
        // If found then erase it from string
        mainStr.erase(pos, toErase.length());
    }
}

//! https://thispointer.com/how-to-remove-substrings-from-a-string-in-c/
void eraseSubStrings(
        std::string & mainStr,
        const std::vector<std::string> & strList)
{
    // Iterate over the given list of substrings. For each substring call eraseAllSubStr() to
    // remove its all occurrences from main string.
    std::for_each(strList.begin(), strList.end(), std::bind(eraseAllSubStr, std::ref(mainStr), std::placeholders::_1));
}


void printEigenMatrixXdToFile( const Eigen::MatrixXd & matrixToPrint,
                               const std::string &filename,
                               const std::string &outputLocation )
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
                                            outputLocation );
    //}
    //else
    /* {
        tudat::input_output::writeMatrixToFile( matrixToPrint,
                                                "fitness_" + fileSuffix + ".dat",
                                                16,
                                                bislip::Variables::getOutputPath( ) + outputSubFolder);
    } */
}



unsigned int getMillisSincePlayTime(
        const std::chrono::time_point< std::chrono::system_clock > &playTime,
        const std::chrono::time_point< std::chrono::system_clock > &printTime)
{
    //! get start time
    // std::time_t timerStart = std::chrono::system_clock::from_time_t( playTime );
    // std::tm *startTimeTM= std::localtime( &timerStart );

    //! current time
    //    std::time_t timerStop = std::chrono::system_clock::to_time_t( printTime );
    //std::tm *stopTimeTM= std::localtime( &timerStop );
    //auto startTime = std::chrono::system_clock::from_time_t( std::mktime( startTimeTM ) );

    // number of milliseconds between start time and now, ie current time in millis
    return std::chrono::duration_cast< std::chrono::milliseconds >( printTime - playTime ).count();
}
unsigned int millis_since_midnight( )
{
    // current time
    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();

    // get midnight
    time_t tnow = std::chrono::system_clock::to_time_t( now );
    tm *date = std::localtime( &tnow );
    date->tm_hour = 0;
    date->tm_min = 0;
    date->tm_sec = 0;
    auto midnight = std::chrono::system_clock::from_time_t( std::mktime( date ) );

    // number of milliseconds between midnight and now, ie current time in millis
    // The same technique can be used for time since epoch
    return std::chrono::duration_cast<std::chrono::milliseconds>( now - midnight ).count();
}

std::chrono::time_point< std::chrono::system_clock > getDateTime( )
{
    // current date/time based on current system
    //std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();

    //std::time_t playTime = std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() );
    // current time
    std::chrono::time_point<std::chrono::system_clock> playTime = std::chrono::system_clock::now();

    return playTime;

}


std::string convertDateTimeToString( const bool &useLocalTime, const std::chrono::time_point<std::chrono::system_clock> &playTime )
{

    std::tm *ptmNow;
    time_t time = std::chrono::system_clock::to_time_t( playTime );

    if( useLocalTime ){ ptmNow = localtime( &time ); }
    else { ptmNow = gmtime( &time ); }

    std::stringstream currentDateTime;
    currentDateTime << 1900 + ptmNow->tm_year;

    //month
    if ( ptmNow->tm_mon < 9 )
        //Fill in the leading 0 if less than 10
        currentDateTime << "0" << 1 + ptmNow->tm_mon;
    else
        currentDateTime << ( 1 + ptmNow->tm_mon );

    //day
    if ( ptmNow->tm_mday < 10 )
        currentDateTime << "0" << ptmNow->tm_mday << "_";
    else
        currentDateTime <<  ptmNow->tm_mday << "_";

    //hour
    if ( ptmNow->tm_hour < 10 )
        currentDateTime << "0" << ptmNow->tm_hour;
    else
        currentDateTime << ptmNow->tm_hour;

    //min
    if ( ptmNow->tm_min < 10 )
        currentDateTime << "0" << ptmNow->tm_min;
    else
        currentDateTime << ptmNow->tm_min;

    //sec
    if ( ptmNow->tm_sec < 10 )
        currentDateTime << "0" << ptmNow->tm_sec;
    else
        currentDateTime << ptmNow->tm_sec;

    return currentDateTime.str();

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
        const double &maximum,
        const double &initialGuess )
{

    double a = minimum;
    double b = maximum;
    double root = ( a + b ) / 2.0;
    if( ( a < initialGuess ) && ( initialGuess < b ) )
    { root = initialGuess; }

    double previousRoot = 0.0;
    double f_a, f_root;
    if( function( a ) * function( b ) < 0.0 )
    {
        //root = ( a + b ) / 2.0;
        unsigned int p = 0;
        /*
        std::cout << "     Starting Bisection search:" << std::endl;
        std::cout << "          function( " << a << " )  = " << function( a ) << std::endl;
        std::cout << "          function( " << b << " )  = " << function( b ) << std::endl;
        std::cout << "          function( " << initialGuess << " )  = " << function( initialGuess ) << std::endl;
        std::cout << "          function( " << root << " )  = " << function( root ) << std::endl;
         */
        do
        {
            f_a       = function( a );
            f_root    = function( root );

            if ( f_a * f_root < 0.0 ) { b = root; }
            else { a = root; }


            previousRoot = root;
            root = ( a + b ) / 2.0;


            // std::cout << "function( " << root << " )  = " << std::abs( function( root ) ) << "  |  Truncated dif = " <<  std::abs( std::trunc( ( previousRoot - root ) * 1000000000 ) ) / 1000000000 << std::endl;

            if( p > 1000 )
            {

                std::cout << "function( " << root << " )  = " << std::abs( function( root ) ) << "  |  Truncated dif = " <<  std::abs( std::trunc( ( previousRoot - root ) * 1000000000 ) ) / 1000000000 << std::endl;


                //throw std::runtime_error( "Equilibrium wall temperature could not converge" );
            }
            p += 1;

        } while( ( std::abs( function( root ) ) > 1e-9 ) && p < 1000 );//&& ( std::abs( std::trunc( ( previousRoot - root ) * 1000000000 ) ) / 1000000000 > 0.0 ) )

    }
    return root;
}

double goldenSectionSearch(
        const std::function< double( const double ) > &function,
        const double &minimum,
        const double &maximum )
{
    double a = minimum;
    double b = maximum;
    double gr = ( std::sqrt( 5 ) + 1 ) / 2;
    double c = b - ( b - a ) / gr;
    double d = a + ( b - a ) / gr;

    //  std::cout << "          Goldern Section Search" << std::endl;



    //std::cout << "              function( " << a << " ) = " << std::abs( function( a ) ) << std::endl;
    //std::cout << "              function( " << b << " ) = " << std::abs( function( b ) ) << std::endl;

    if( ( function( a ) > 0.0 ) && ( function( b ) > 0.0 ) )
    {
        while ( std::abs( c - d ) > 1e-9 )
        {
            if( function( c ) * function( d ) < 0.0 )
            { b = d; }
            else
            { a = c; }

            c = b - ( b - a ) / gr;
            d = a + ( b - a ) / gr;

            //      std::cout << "              function( " << c << " ) = " << std::abs( function( c ) ) << std::endl;
            //    std::cout << "              function( " << d << " ) = " << std::abs( function( d ) ) << std::endl;
        }
    }
    else if( ( function( a ) < 0.0 ) && ( function( b ) < 0.0 ) )
    {
        while ( std::abs( c - d ) > 1e-9 )
        {
            if( function( c ) * function( d ) > 0.0 )
            { b = d; }
            else
            { a = c; }

            c = b - ( b - a ) / gr;
            d = a + ( b - a ) / gr;

            //  std::cout << "              function( " << c << " ) = " << std::abs( function( c ) ) << std::endl;
            //  std::cout << "              function( " << d << " ) = " << std::abs( function( d ) ) << std::endl;
        }
    }
    return ( b + a ) / 2;
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

double getCurrentAltitude(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    const double currentAltitude = flightConditions->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );

    return currentAltitude;
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

bislip::Parameters::Interpolators passOptimizationParameter (
        const bislip::Parameters::Interpolators &parameter )
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
        const bislip::Parameters::Interpolators &parameter,
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    double currentHeight = 0.0;
    double currentAirspeed = 0.0;

    if( debugInfo == 1 ){ std::cout << "Selecting source of airspeed/altitude" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Current Time  = " << flightConditions->getCurrentTime() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Starting Time = " << bislipSystems->getStartingEpoch() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Guidance Step = " << bislipSystems->getGuidanceStepSize() << std::endl; }


    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "     Selecting initial airspeed/altitude" << std::endl; }

        currentHeight   = bislipSystems->getInitialHeight();
        currentAirspeed = bislipSystems->getInitialAirspeed();
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "     Selecting propagated airspeed/altitude" << std::endl; }

        currentHeight   = flightConditions->getCurrentAltitude();
        currentAirspeed = flightConditions->getCurrentAirspeed();
    }

    double normalizedSpecificEnergy = bislip::Variables::computeNormalizedSpecificEnergy( currentHeight, currentAirspeed, bislipSystems->getMaximumSpecificEnergy() );

    if( normalizedSpecificEnergy > 1.0 ) { normalizedSpecificEnergy = 1.0; }

    if( debugInfo == 1 ){ std::cout << "     Current Height             = " << currentHeight << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Current Airspeed           = " << currentAirspeed << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Maximum Specific Energy    = " << bislipSystems->getMaximumSpecificEnergy() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Normalized Specific Energy = " << normalizedSpecificEnergy << std::endl; }

    if( debugInfo == 1 ){ std::cout << "     Evaluating Interpolator for: " << parameter << std::endl; }
    //! Evaluate interpolator.
    double evaluation = ( bislipSystems->getParameterInterpolator( parameter ) )->interpolate( normalizedSpecificEnergy );
    if( debugInfo == 1 ){ std::cout << "        Evaluation = " << evaluation << std::endl; }


    if( debugInfo == 1 ){ std::cout << "        Extract parameter bounds" << std::endl; }
    //! Declare and assign parameter bounds.
    double lowerBound = ( bislipSystems->getParameterLowerBounds( parameter ) )->interpolate( normalizedSpecificEnergy );
    double upperBound = ( bislipSystems->getParameterUpperBounds( parameter ) )->interpolate( normalizedSpecificEnergy );

    if( debugInfo == 1 ){ std::cout << "            Lower Bound = " << lowerBound << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Upper Bound = " << upperBound << std::endl; }




    if( debugInfo == 1 ){ std::cout << "        Impose parameter bounds" << std::endl; }
    if ( evaluation < lowerBound ){ evaluation = lowerBound; }
    if ( evaluation > upperBound ){ evaluation = upperBound; }

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

    if( debugInfo == 1 ){ std::cout << "            Compute Body-Fixed Thrust Direction" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            --------------------------------" << std::endl; }

    //! Evaluate the Thrust Elevation Angle Interpolator and convert to radians.
    double thrustElevationAngle = bislipSystems->getCurrentThrustElevationAngle();

    if( debugInfo == 1 ){ std::cout << "                Thrust Elevation Angle for Body-Fixed Thrust Direction = " << thrustElevationAngle << std::endl; }

    //! Evaluate the Thrust Azimuth Angle Interpolator and convert to radians.
    double thrustAzimuthAngle = bislipSystems->getCurrentThrustAzimuthAngle();

    if( debugInfo == 1 ){ std::cout << "                Thrust Azimuth Angle for Body-Fixed Thrust Direction   = " << thrustAzimuthAngle << std::endl; }

    //! Declare and calculate the body-fixed thrust direction.
    Eigen::Vector3d bodyFixedThrustDirection;
    bodyFixedThrustDirection( 0 ) = std::cos( thrustElevationAngle ) * std::cos( thrustAzimuthAngle );
    bodyFixedThrustDirection( 1 ) = std::cos( thrustElevationAngle ) * std::sin( thrustAzimuthAngle );
    bodyFixedThrustDirection( 2 ) = -std::sin( thrustElevationAngle );
    if( debugInfo == 1 ){ std::cout << "               Body-Fixed Thrust Direction       = [ " << bodyFixedThrustDirection( 0 ) << " , " << bodyFixedThrustDirection( 1 ) << " , " << bodyFixedThrustDirection( 2 ) << " ]" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Ending Computation of Body-Fixed Total Load" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "            --------------------------------" << std::endl; }

    return bodyFixedThrustDirection;
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
    if( debugInfo == 1 ){ std::cout << "            Starting Compution Body-Fixed Thrust Vector" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            -------------------------------------------" << std::endl; }

    //! Determine the Body-Fixed Thrust Direction Vector.
    //const Eigen::Vector3d BodyFixedThrustDirection = bislipSystems->getCurrentBodyFixedThrustDirection();
    const Eigen::Vector3d BodyFixedThrustDirection = bislip::Variables::computeBodyFixedThrustDirection( bodyMap, vehicleName );
    if( debugInfo == 1 ){ std::cout << "                Body-fixed Thrust Direction Vector = [ " << BodyFixedThrustDirection(0) << " , " << BodyFixedThrustDirection(1) << " , " << BodyFixedThrustDirection(2) << " ]" << std::endl; }

    //! Determine the Thrust Magnitude.
    const double thrustMagnitude = bislip::Variables::computeThrustMagnitude( bodyMap, vehicleName );
    //const double thrustMagnitude = bislipSystems->getCurrentThrustMagnitude();
    if( debugInfo == 1 ){ std::cout << "            Thrust Magnitude for Body-Fixed Thrust Vector = " << thrustMagnitude << std::endl; }


    //! Declare and calculate the body-fixed thrust vector.
    Eigen::Vector3d BodyFixedThrustVector;
    BodyFixedThrustVector( 0 ) = thrustMagnitude * BodyFixedThrustDirection( 0 );
    BodyFixedThrustVector( 1 ) = thrustMagnitude * BodyFixedThrustDirection( 1 );
    BodyFixedThrustVector( 2 ) = thrustMagnitude * BodyFixedThrustDirection( 2 );

    if( debugInfo == 1 ){ std::cout << "                Body-fixed Thrust Load Vector = [ " << BodyFixedThrustVector(0) << " , " << BodyFixedThrustVector(1) << " , " << BodyFixedThrustVector(2) << " ]" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Starting Compution Body-Fixed Thrust Vector" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            -------------------------------------------" << std::endl; }

    return BodyFixedThrustVector;
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

    if( debugInfo == 1 ){ std::cout << "    Starting Computation of Thrust Magnitude" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    --------------------------" << std::endl; }

    const double throttle = bislipSystems->getCurrentThrottleSetting();
    if( debugInfo == 1 ){ std::cout << "        Throttle Setting         = " << throttle << std::endl; }

    const double maxThrustMagnitude = bislipSystems->getMaxThrust();
    if( debugInfo == 1 ){ std::cout << "        Max Thrust Magnitude     = " << maxThrustMagnitude << std::endl; }

    const double engineStatus = double( bislipSystems->getCurrentEngineStatus() );
    if( debugInfo == 1 ){ std::cout << "        Current Engine Status    = " << engineStatus << std::endl; }

    const double currentThrustMagnitude = engineStatus * maxThrustMagnitude * throttle;
    if( debugInfo == 1 ){ std::cout << "        Current Thrust Magnitude = " << currentThrustMagnitude << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Ending Computation of Thrust Magnitude" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    --------------------------" << std::endl; }

    return currentThrustMagnitude;
}

bool determineEngineStatus (
        const double &currentMass,
        const double &landingMass )
//const double &normalizedSpecificEnergy )
{



    bool currentEngineStatus = true;
    if ( currentMass <= landingMass ){ currentEngineStatus = false; }
    //if ( normalizedSpecificEnergy >= 0.95 ){ currentEngineStatus = false; }

    return currentEngineStatus;
}

Eigen::Vector3d computeLocalGravity (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName )
{   
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "         Selecting Source of Values" << std::endl; }
    double currentAltitude;
    double currentLatitude;
    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentAltitude = bislipSystems->getInitialAltitude();
        currentLatitude = bislipSystems->getInitialLat();
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentAltitude = flightConditions->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );
        currentLatitude = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    }

    const double gravitationalParameter = bodyMap.at( centralBodyName )->getGravityFieldModel( )->getGravitationalParameter();
    const double earthRadius = bodyMap.at( centralBodyName )->getShapeModel( )->getAverageRadius();

    if( debugInfo == 1 ){ std::cout << "                Current Time            = " << flightConditions->getCurrentTime() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Altitude        = " << currentAltitude << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Latitude        = " << tudat::unit_conversions::convertDegreesToRadians( currentLatitude ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Gravitational Parameter = " << gravitationalParameter << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Earth Radius            = " << earthRadius << std::endl; }

    const double J2 = 1082.626523e-6;
    const double J3 = 2.532153e-7;
    const double J4 = 1.6109876e-7;

    const double s_latitude    = std::sin( currentLatitude );
    const double s_latitude_sq = s_latitude * s_latitude;
    const double c_latitude    = std::cos( currentLatitude );

    const double g_n_pre = -3.0 * ( gravitationalParameter / ( currentAltitude * currentAltitude ) ) * pow( earthRadius / currentAltitude , 2.0 ) * s_latitude * c_latitude;
    const double g_n_sub_1 = J2;
    const double g_n_sub_2 = ( J3 / 2.0 ) * ( earthRadius / currentAltitude ) * ( 1.0 / s_latitude ) * ( 5.0 * s_latitude_sq - 1.0 );
    const double g_n_sub_3 = ( 5.0 * J4 / 6.0 ) * pow( earthRadius / currentAltitude , 2.0 ) * ( 7.0 * s_latitude_sq - 3.0 );
    const double g_n = g_n_pre * ( g_n_sub_1 + g_n_sub_2 + g_n_sub_3 );

    const double g_d_pre = ( gravitationalParameter / ( currentAltitude * currentAltitude ) );
    const double g_d_sub_1 = ( 3.0 * J2 / 2.0 ) * pow( earthRadius / currentAltitude , 2.0 ) * ( 3.0 * s_latitude_sq - 1.0 );
    const double g_d_sub_2 = 2.0 * J3  * pow( earthRadius / currentAltitude , 3.0 ) * ( s_latitude ) * ( 5.0 * s_latitude_sq - 3.0 );
    const double g_d_sub_3 = ( 5.0 * J4 / 8.0 ) * pow( earthRadius / currentAltitude , 4.0 ) * ( 35.0 * s_latitude_sq * s_latitude_sq - 30.0 * s_latitude_sq + 3.0 );
    const double g_d = g_d_pre * ( 1.0 - g_d_sub_1 - g_d_sub_2 - g_d_sub_3 );

    Eigen::Vector3d localGravity ( 3 );
    localGravity << g_n, 0, g_d;


    //bislipSystems->setCurrentLocalGravityVector( localGravity );

    return localGravity;
}

double computeCurrentDragForce(
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

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){std::cout << "         Extract Full Current Coefficients " << std::endl; }

    const Eigen::Vector6d currentCoefficients = bislipSystems->getFullCurrentCoefficients();

    if( debugInfo == 1 ){ std::cout << "         Selecting Source of Dynamic Pressure" << std::endl; }

    double currentDynamicPressure;
    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentDynamicPressure = bislipSystems->getInitialDynamicPressure();
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    }

    if( debugInfo == 1 ){std::cout << "         Compute Drag Force " << std::endl; }

    const double currentDragForce = ( currentDynamicPressure ) * ( bislipSystems->getReferenceArea( ) ) * currentCoefficients[ 0 ];

    if( debugInfo == 1 ){std::cout << "         Set Drag Force " << std::endl; }
    bislipSystems->setCurrentDragForce( currentDragForce );

    return currentDragForce;


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

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){std::cout << "         Determine Full Current Coefficients " << std::endl; }

    const Eigen::Vector6d currentCoefficients = bislipSystems->getFullCurrentCoefficients();

    if( debugInfo == 1 ){ std::cout << "         Selecting Source of Dynamic Pressure" << std::endl; }

    double currentDynamicPressure;
    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentDynamicPressure = bislipSystems->getInitialDynamicPressure();
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    }


    //std::cout << "Mach Number:                           " << flightConditions->getCurrentMachNumber() << std::endl;
    //std::cout << "currentAngleOfAttack:                  " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentAngleOfAttack() ) << std::endl;
    //std::cout << "Lift Coefficient:                      " << currentCoefficients[ 2 ] << std::endl;
    //std::cout << "Reference Area:                        " << bislipSystems->getReferenceArea( ) << std::endl;
    if( debugInfo == 1 ){std::cout << "         Compute Lift Force " << std::endl; }

    const double currentLiftForce = ( currentDynamicPressure ) * ( bislipSystems->getReferenceArea( ) ) * currentCoefficients[ 2 ];

    if( debugInfo == 1 ){std::cout << "         Set Lift Force " << std::endl; }

    bislipSystems->setCurrentLiftForce( currentLiftForce );

    return currentLiftForce;
}


double optimizeAngleOfAttack(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;


    int debugInfo = bislipSystems->getDebugInfo();



    if( debugInfo == 1 ){ std::cout << "         Selecting Source of Values" << std::endl; }
    double currentHeight;
    double currentAirspeed;
    double currentMachNumber;
    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentHeight     = bislipSystems->getInitialHeight();
        currentAirspeed   = bislipSystems->getInitialAirspeed();
        currentMachNumber = bislipSystems->getInitialMachNumber();
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentHeight     = flightConditions->getCurrentAltitude();
        currentAirspeed   = flightConditions->getCurrentAirspeed();
        currentMachNumber = flightConditions->getCurrentMachNumber();
    }


    if( debugInfo == 1 ){ std::cout << "        Starting Optimization of Angle of Attack" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "        --------------------------" << std::endl; }

    const double normalizedSpecificEnergy    = bislip::Variables::computeNormalizedSpecificEnergy( currentHeight, currentAirspeed, bislipSystems->getMaximumSpecificEnergy() );

    // std::cout << "normalizedSpecificEnergy:      " << normalizedSpecificEnergy << std::endl;

    double currentAngleOfAttack = bislipSystems->getCurrentAngleOfAttack();

    const double maximumMechanicalLoad       = bislipSystems->getMechanicalLoadConstraint();

    if( debugInfo == 1 ){ std::cout << "Extract Angle of Attack Bounds." << std::endl; }
    double angleOfAttackLowerBound_1 = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getParameterLowerBounds( bislip::Parameters::Interpolators::AngleOfAttack ) )->interpolate( normalizedSpecificEnergy ) );
    double angleOfAttackUpperBound_1 = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getParameterUpperBounds( bislip::Parameters::Interpolators::AngleOfAttack ) )->interpolate( normalizedSpecificEnergy ) );
    double angleOfAttackLowerBound_2 = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getAlphaMachEnvelopeLBInterpolator() )->interpolate( currentMachNumber ) );
    double angleOfAttackUpperBound_2 = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getAlphaMachEnvelopeUBInterpolator() )->interpolate( currentMachNumber ) );

    double angleOfAttackLowerBound = std::max( angleOfAttackLowerBound_1, angleOfAttackLowerBound_2 );
    double angleOfAttackUpperBound = std::min( angleOfAttackUpperBound_1, angleOfAttackUpperBound_2 );

    if( debugInfo == 1 ){ std::cout << "    Angle of attack Lower Bound: " << angleOfAttackLowerBound << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    Angle of attack Upper Bound: " << angleOfAttackUpperBound << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Current Angle of Attack     = " << currentAngleOfAttack << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Angle of Attack Lower Bound = " << angleOfAttackLowerBound << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Angle of Attack Upper Bound = " << angleOfAttackUpperBound << std::endl; }


    if( debugInfo == 1 ){ std::cout << "           Create Local Function to Calculate Total Body-Fixed g-load Given an Angle of Attack" << std::endl; }
    std::function< double( const double ) > angleofAttackEvaluation =
            std::bind( &bislip::Variables::angleofAttackEvaluationFunction, bodyMap, vehicleName, std::placeholders::_1 );

    if( debugInfo == 1 ){ std::cout << "           Evaluating the Current Angle of Attack" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "           -----------------------------------------------" << std::endl; }

    const double currentBodyFixedTotal_g_Load_Magnitude = angleofAttackEvaluation( currentAngleOfAttack );
    if( debugInfo == 1 ){ std::cout << "           Evaluation of Current Angle of Attack Complete" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "           -----------------------------------------------" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Mechanical Load Limit        = " << maximumMechanicalLoad << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Mechanical Load      = " << currentBodyFixedTotal_g_Load_Magnitude << std::endl; }


    //if( debugInfo == 1 ){ std::cout << "            Evaluated Throttle Setting:    " << evaluatedThrottleSetting << std::endl;
    if( debugInfo == 1 ){ std::cout << "  " << std::endl; }

    double optimizedAngleofAttack = currentAngleOfAttack;
    double optimizedMechanicalLoad = 0.0;

    // if( bislipSystems->getCurrentEngineStatus() == true )
    // {
    //   if( debugInfo == 1 ){ std::cout << "            Engine is Currently ON" << std::endl; }

    // double noThrustBodyFixedTotal_g_Load_Magnitude    = throttleSettingEvaluation( 0.0 );
    double lowerThrustBodyFixedTotal_g_Load_Magnitude = angleofAttackEvaluation( angleOfAttackLowerBound );
    double upperThrustBodyFixedTotal_g_Load_Magnitude = angleofAttackEvaluation( angleOfAttackUpperBound );

    if( debugInfo == 1 ){ std::cout << "            Lowest Angle of Attack Mechanical Load:    " << lowerThrustBodyFixedTotal_g_Load_Magnitude << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Highest Angle of Attack Mechanical Load:   " << upperThrustBodyFixedTotal_g_Load_Magnitude << std::endl; }

    if( lowerThrustBodyFixedTotal_g_Load_Magnitude > maximumMechanicalLoad )
    {
        optimizedAngleofAttack  = angleOfAttackLowerBound;
        optimizedMechanicalLoad = angleofAttackEvaluation( angleOfAttackLowerBound );
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "            Mechanical Load Limit has NOT been reached" << std::endl; }
        if( debugInfo == 1 ){ std::cout << "                Create Local Function to Search Angle of Attack to Optimize the Mechanical Load" << std::endl; }

        std::function< double( const double ) > angleOfAttackOptimizationFunction =
                std::bind( &bislip::Variables::angleOfAttackOptimizationFunction, bodyMap, vehicleName, std::placeholders::_1 );

        double f_a = angleOfAttackOptimizationFunction( angleOfAttackLowerBound );
        double f_b = angleOfAttackOptimizationFunction( angleOfAttackUpperBound );
        if( debugInfo == 1 ){ std::cout << "                f( " << angleOfAttackLowerBound << " ) = " << f_a << std::endl; }
        if( debugInfo == 1 ){ std::cout << "                f( " << angleOfAttackUpperBound << " ) = " << f_b << std::endl; }

        if( f_a * f_b < 0 )
        {
            if( debugInfo == 1 ){ std::cout << "                Bisection Search" << std::endl; }
            optimizedAngleofAttack = bislip::Variables::rootFinderBisection( angleOfAttackOptimizationFunction, angleOfAttackLowerBound, angleOfAttackUpperBound, currentAngleOfAttack );
            if( debugInfo == 1 ){ std::cout << "                Bisection Search Complete" << std::endl; }
            optimizedMechanicalLoad = angleofAttackEvaluation( optimizedAngleofAttack );
        }
        else
        {
            if( debugInfo == 1 ){ std::cout << "                Golden Ratio Search" << std::endl; }
            optimizedAngleofAttack = bislip::Variables::goldenSectionSearch( angleOfAttackOptimizationFunction, angleOfAttackLowerBound, angleOfAttackUpperBound );
            if( debugInfo == 1 ){ std::cout << "                Golden Ratio Search Complete" << std::endl; }
            optimizedMechanicalLoad = angleofAttackEvaluation( optimizedAngleofAttack );
        }
    }
    //}
    if( debugInfo == 1 ){ std::cout << "            Optimized Angle of Attack = " << optimizedAngleofAttack << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Optimized Mechanical Load = " << optimizedMechanicalLoad << std::endl; }
    if( debugInfo == 1 ){ std::cout << "        Ending Optimization of Angle of Attack" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "        --------------------------" << std::endl; }

    return optimizedAngleofAttack;
}


double angleofAttackEvaluationFunction(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const double &angleOfAttack )

{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Extract Vehicle Systems pointer.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( );

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "                Setting new Angle of Attack to Evaluate its Effect" << std::endl; }
    bislipSystems->setCurrentAngleOfAttack( angleOfAttack );

    const double n = bislip::Variables::computeBodyFixedTotal_g_Load_Magnitude( bodyMap, vehicleName );


    return n;
}

double angleOfAttackOptimizationFunction(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const double &angleOfAttack )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;


    double evaluation = bislipSystems->getMechanicalLoadConstraint() - bislip::Variables::angleofAttackEvaluationFunction( bodyMap, vehicleName, angleOfAttack );

    //std::cout << "      evaluation: " << evaluation << std::endl;

    return evaluation;
}


double determineThrottleSetting(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();
    if( debugInfo == 1 ){ std::cout << "        Starting Determination of Throttle Setting" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "        --------------------------" << std::endl; }

    const double maximumMechanicalLoad       = bislipSystems->getMechanicalLoadConstraint();
    //const double normalizedSpecificEnergy    = bislip::Variables::computeNormalizedSpecificEnergy( flightConditions->getCurrentAltitude(), flightConditions->getCurrentAirspeed(), bislipSystems->getMaximumSpecificEnergy() );

    // std::cout << "normalizedSpecificEnergy:      " << normalizedSpecificEnergy << std::endl;

    double currentThrottleSetting            = 0.0;
    double currentFlightPathAngle            = 0.0;
    double currentFlightPathAngleRate        = 0.0;



    if( std::isnan( bislipSystems->getCurrentThrottleSetting() ) == true )
    {
        if( debugInfo == 1 ){ std::cout << "           Throttle Setting has not been previously set. Assigning a value of '0.5'." << std::endl; }

        currentThrottleSetting     = 0.5;
        bislipSystems->setCurrentThrottleSetting( currentThrottleSetting );
        bislipSystems->setCurrentThrustMagnitude(  bislip::Variables::computeThrustMagnitude( bodyMap, vehicleName ) );
        currentFlightPathAngle     = bislipSystems->getInitialFlightPathAngle();
        currentFlightPathAngleRate = bislip::Variables::computeFlightPathAngleRate( bodyMap, vehicleName, centralBodyName );//  currentThrottleSetting = bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::ThrottleSetting, bodyMap, vehicleName, centralBodyName );
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "           Extracting previous Throttle Setting" << std::endl; }

        currentThrottleSetting     = bislipSystems->getCurrentThrottleSetting();
        currentFlightPathAngle     = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        currentFlightPathAngleRate = bislipSystems->getCurrentFlightPathAngleRate();
    }


    //    const double evaluatedThrottleSetting    = bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::ThrottleSetting, bodyMap, vehicleName );
    //const double throttleSettingLowerBound   = ( bislipSystems->getParameterLowerBounds( bislip::Parameters::Interpolators::ThrottleSetting ) )->interpolate( normalizedSpecificEnergy );
    //const double throttleSettingUpperBound   = ( bislipSystems->getParameterUpperBounds( bislip::Parameters::Interpolators::ThrottleSetting ) )->interpolate( normalizedSpecificEnergy );
    const double throttleSettingLowerBound = ( bislipSystems->getThrottleSettingLimits() ).first;
    const double throttleSettingUpperBound = ( bislipSystems->getThrottleSettingLimits() ).second;
    if( debugInfo == 1 ){ std::cout << "                Current Throttle Setting     = " << currentThrottleSetting << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Throttle Setting Lower Bound = " << throttleSettingLowerBound << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Throttle Setting Upper Bound = " << throttleSettingUpperBound << std::endl; }



    if( debugInfo == 1 ){ std::cout << "           Create Local Function to Calculate Total Body-Fixed g-load Given a Throttle Setting" << std::endl; }
    std::function< double( const double ) > throttleSettingEvaluation =
            std::bind( &bislip::Variables::throttleSettingEvaluationFunction, bodyMap, vehicleName, centralBodyName, std::placeholders::_1 );

    if( debugInfo == 1 ){ std::cout << "           Evaluating the Current Throttle Setting" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "           -----------------------------------------------" << std::endl; }

    const double currentBodyFixedTotal_g_Load_Magnitude = throttleSettingEvaluation( currentThrottleSetting );
    if( debugInfo == 1 ){ std::cout << "           Evaluation of Current Throttle Setting Complete" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "           -----------------------------------------------" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Mechanical Load Limit        = " << maximumMechanicalLoad << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Mechanical Load      = " << currentBodyFixedTotal_g_Load_Magnitude << std::endl; }


    //if( debugInfo == 1 ){ std::cout << "            Evaluated Throttle Setting:    " << evaluatedThrottleSetting << std::endl;
    if( debugInfo == 1 ){ std::cout << "  " << std::endl; }

    double limitedThrottleSetting = currentThrottleSetting;
    double limitedMechanicalLoad = 0.0;

    //    if( ( bislipSystems->getCurrentEngineStatus() == true && currentFlightPathAngle >= 0 ) || ( bislipSystems->getCurrentEngineStatus() == true && currentFlightPathAngle < 0 && currentFlightPathAngleRate > 0 ) )
    //    if( ( bislipSystems->getCurrentEngineStatus() == true && currentFlightPathAngle >= 0 ) || ( bislipSystems->getCurrentEngineStatus() == true && currentFlightPathAngle < 0 && currentFlightPathAngleRate > 0 ) )
    if( bislipSystems->getCurrentEngineStatus() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Engine is Currently ON" << std::endl; }

        double noThrustBodyFixedTotal_g_Load_Magnitude    = throttleSettingEvaluation( 0.0 );
        double lowerThrustBodyFixedTotal_g_Load_Magnitude = throttleSettingEvaluation( throttleSettingLowerBound );
        double upperThrustBodyFixedTotal_g_Load_Magnitude = throttleSettingEvaluation( throttleSettingUpperBound );

        if( debugInfo == 1 ){ std::cout << "            No Thrust Mechanical Load      = " << noThrustBodyFixedTotal_g_Load_Magnitude << std::endl; }
        if( debugInfo == 1 ){ std::cout << "            Lowest Thrust Mechanical Load  = " << lowerThrustBodyFixedTotal_g_Load_Magnitude << std::endl; }
        if( debugInfo == 1 ){ std::cout << "            Highest Thrust Mechanical Load = " << upperThrustBodyFixedTotal_g_Load_Magnitude << std::endl; }

        if( lowerThrustBodyFixedTotal_g_Load_Magnitude >= maximumMechanicalLoad )
        {
            if( debugInfo == 1 ){ std::cout << "            Mechanical Load Limit has been exceeded. Assigning the lowest possible Throttle Setting." << std::endl; }

            limitedThrottleSetting = throttleSettingLowerBound;
            limitedMechanicalLoad  = lowerThrustBodyFixedTotal_g_Load_Magnitude;
        }
        else if( upperThrustBodyFixedTotal_g_Load_Magnitude <= maximumMechanicalLoad )
        {
            if( debugInfo == 1 ){ std::cout << "            Highest possible Throttle Setting does not exceed Mechanical Load Limit. Assigning the highest possible Throttle Setting." << std::endl; }

            limitedThrottleSetting = throttleSettingUpperBound;
            limitedMechanicalLoad  = upperThrustBodyFixedTotal_g_Load_Magnitude;
        }
        else
        {
            if( debugInfo == 1 ){ std::cout << "            Mechanical Load Limit has NOT been reached" << std::endl; }
            if( debugInfo == 1 ){ std::cout << "                Create Local Function to Search Throttle Setting to Maximize the Mechanical Load" << std::endl; }

            std::function< double( const double ) > throttleSettingLimitSearchFunction =
                    std::bind( &bislip::Variables::throttleSettingLimitSearchFunction, bodyMap, vehicleName, centralBodyName, std::placeholders::_1 );

            double f_a = throttleSettingLimitSearchFunction( throttleSettingLowerBound );
            double f_b = throttleSettingLimitSearchFunction( throttleSettingUpperBound );
            if( debugInfo == 1 ){ std::cout << "                f( " << throttleSettingLowerBound << " ) = " << f_a << std::endl; }
            if( debugInfo == 1 ){ std::cout << "                f( " << throttleSettingUpperBound << " ) = " << f_b << std::endl; }

            if( f_a * f_b < 0 )
            {
                if( debugInfo == 1 ){ std::cout << "                Bisection Search" << std::endl; }
                limitedThrottleSetting = bislip::Variables::rootFinderBisection( throttleSettingLimitSearchFunction, throttleSettingLowerBound, throttleSettingUpperBound, currentThrottleSetting );
                if( debugInfo == 1 ){ std::cout << "                Bisection Search Complete" << std::endl; }
                limitedMechanicalLoad = throttleSettingEvaluation( limitedThrottleSetting );
            }
            else
            {
                if( debugInfo == 1 ){ std::cout << "                Golden Ratio Search" << std::endl; }
                limitedThrottleSetting = bislip::Variables::goldenSectionSearch( throttleSettingLimitSearchFunction, throttleSettingLowerBound, throttleSettingUpperBound );
                if( debugInfo == 1 ){ std::cout << "                Golden Ratio Search Complete" << std::endl; }
                limitedMechanicalLoad = throttleSettingEvaluation( limitedThrottleSetting );
            }
        }
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "                Negative Flight-Path Angle with Negative Flight-path Angle Rate" << std::endl; }
        limitedThrottleSetting = 0.0;
        limitedMechanicalLoad = throttleSettingEvaluation( limitedThrottleSetting );
    }
    if( debugInfo == 1 ){ std::cout << "            Limited Throttle Setting:   " << limitedThrottleSetting << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Limited Mechanical Load:    " << limitedMechanicalLoad << std::endl; }
    if( debugInfo == 1 ){ std::cout << "        Ending Determination of Throttle Setting" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "        --------------------------" << std::endl; }

    return limitedThrottleSetting;
}


double throttleSettingEvaluationFunction(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName,
        const double &throttleSetting )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Extract Vehicle Systems pointer.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( );

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "                Setting new Throttle Setting to Evaluate its Effect" << std::endl; }
    bislipSystems->setCurrentThrottleSetting( throttleSetting );

    //if( debugInfo == 1 ){ std::cout << "                Setting new Thrust Magnitude to Evaluate its Effect" << std::endl; }
    //bislipSystems->setCurrentThrustMagnitude( bislip::Variables::computeThrustMagnitude( bodyMap, vehicleName ) );

    //if( debugInfo == 1 ){ std::cout << "                Setting new Thrust Magnitude to Evaluate its Effect" << std::endl; }

    //bislipSystems->setCurrentBodyFixedThrustDirection( bislip::Variables::computeBodyFixedThrustDirection( bodyMap, vehicleName ) );
    /*
    if( debugInfo == 1 ){ std::cout << "                Setting new Control Surfaces to Evaluate their Effect" << std::endl; }
   Eigen::Vector2d controlSurfaceDeflections = bislip::Variables::computeControlSurfaceDeflection( bodyMap, vehicleName );
    bislipSystems->setCurrentBodyFlapAngle( controlSurfaceDeflections( 0 ) );
    bislipSystems->setCurrentElevonAngle( controlSurfaceDeflections( 1 ) );

    vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", bislipSystems->getCurrentBodyFlapAngle() );
    vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonLeft", bislipSystems->getCurrentElevonAngle() );
    vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonRight", bislipSystems->getCurrentElevonAngle() );

    bislipSystems->setCurrentLiftForce( bislip::Variables::computeCurrentLiftForce( bodyMap, vehicleName ) );
    bislipSystems->setCurrentDragForce( bislip::Variables::computeCurrentDragForce( bodyMap, vehicleName ) );


    bislipSystems->setCurrentThrustElevationAngle( bislip::Variables::determineThrustElevationAngle( bodyMap, vehicleName, centralBodyName ) );

*/
    const double n = bislip::Variables::computeBodyFixedTotal_g_Load_Magnitude( bodyMap, vehicleName );
    /*

    double currentFlightPathAngle;
    double currentBankAngle;
    double currentLatitude;
    double currentHeading;
    double currentAltitude;
    double currentAirspeed;
    double currentDynamicPressure;

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentFlightPathAngle = bislipSystems->getInitialFlightPathAngle();
        currentBankAngle       = bislipSystems->getCurrentBankAngle();
        currentLatitude        = bislipSystems->getInitialLat();
        currentHeading         = bislipSystems->getInitialHeadingAngle();
        currentAltitude        = bislipSystems->getInitialAltitude();
        currentAirspeed        = bislipSystems->getInitialAirspeed();
        currentDynamicPressure = bislipSystems->getInitialDynamicPressure();
    }

    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentFlightPathAngle = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        currentBankAngle       = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::bank_angle );
        currentLatitude        = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        currentHeading         = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
        currentAltitude        = flightConditions->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );
        currentAirspeed        = flightConditions->getCurrentAirspeed( );
        currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    }


    //const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    const double currentMass          = bodyMap.at( vehicleName )->getBodyMass();
    const double currentAngleOfAttack = bislipSystems->getCurrentAngleOfAttack();

    const double rotationRateEarth = 7.292115*1E-5;

    //const Eigen::Vector6d currentCoefficients = bislip::Variables::computePartialCurrentCoefficients( bodyMap, vehicleName );
    const double currentLiftForce = bislipSystems->getCurrentLiftForce();
    const double currentDragForce = bislipSystems->getCurrentDragForce();
    const double currentThrustMagnitude = bislip::Variables::computeThrustMagnitude( bodyMap, vehicleName );

    const Eigen::Vector3d aerodynamicFrameTotal_g_Load_Vector =
            ( tudat::reference_frames::getBodyToAirspeedBasedAerodynamicFrameTransformationMatrix( currentAngleOfAttack, 0.0 ) ) * ( bislip::Variables::computeBodyFixedTotal_g_Load_Vector( bodyMap, vehicleName ) );

    const Eigen::Vector3d aerodynamicFrameTotal_Load_Vector = bislip::Variables::computeAerodynamicFrameTotalLoad( bodyMap, vehicleName );
    const double aerodynamicFrameTotal_g_Load_x = aerodynamicFrameTotal_g_Load_Vector( 0 );



  const double n_aero_xg = aerodynamicFrameTotal_g_Load_x * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION;

    const Eigen::Vector2d gravs = bislipSystems->getCurrentLocalGravityVector();
    // const double currentThrustMagnitude      = bislipSystems->getCurrentThrustMagnitude();
    const double currentThrustElevationAngle = bislipSystems->getCurrentThrustElevationAngle();
    const double currentThrustAzimuthAngle   = bislipSystems->getCurrentThrustAzimuthAngle();


    const double s_alpha = std::sin( currentAngleOfAttack );
    const double c_alpha = std::cos( currentAngleOfAttack );
    const double s_delta = std::sin( currentLatitude );
    const double c_delta = std::cos( currentLatitude );
    // const double s_gamma = std::sin( currentFlightPathAngle );
    // const double c_gamma = std::cos( currentFlightPathAngle );
     const double s_sigma = std::sin( currentBankAngle );
     const double c_sigma = std::cos( currentBankAngle );
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
    //!     E = sqrt( B * B + C * C )
    //!     phi = arctan( B, C )

     const double B1 = currentThrustMagnitude * s_thrustAzimuthAngle * c_thrustElevationAngle / currentMass;
     const double C1 = ( currentLiftForce + currentThrustMagnitude * ( s_alpha * c_thrustAzimuthAngle * c_thrustElevationAngle + c_alpha * s_thrustElevationAngle ) ) / currentMass;
     const double E1 = std::sqrt( ( B1 * B1 ) + ( C1 * C1 ) );

     double phi1 = 0;
     if( B1 != 0.0 ) { phi1 = std::atan2( B1 , C1 ); }





     const double B2 = gravs( 0 ) * c_chi + rotationRateEarth * rotationRateEarth * currentAltitude * c_delta * s_delta * c_chi;
     const double C2 = ( currentAirspeed * currentAirspeed / currentAltitude ) - gravs( 1 ) + rotationRateEarth * rotationRateEarth * currentAltitude * c_delta * c_delta;
     const double E2 = std::sqrt( ( B2 * B2 ) + ( C2 * C2 ) );
     double phi2 = 0;
     if( B2 != 0.0 ) { phi2 = std::atan2( B2 , C2 ); }


     const double numerator_sigma_term = std::copysign( E1 , B1 ) * std::cos( currentBankAngle + phi1 );
     const double numerator_other_term = -2 * rotationRateEarth * currentAirspeed * c_delta * s_chi;
     const double numerator = numerator_other_term - numerator_sigma_term;
     const double denominator = std::copysign( E2 , B2 );
     const double arg = numerator / denominator;
     const double acosArg = std::acos( arg );
     const double gamma = acosArg - phi2;
     const double TminusD = aerodynamicFrameTotal_Load_Vector(0);
     const double currentWeight =  (currentMass * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION);

             const double sin_gamma = TminusD / currentWeight;
             const double gamma1 = std::asin(sin_gamma);

     //const double B = rotationRateEarth * rotationRateEarth * currentAltitude * c_delta * c_delta - gravs( 1 );
     //const double C = -( rotationRateEarth * rotationRateEarth * currentAltitude * c_delta * s_delta + gravs( 0 ) ) * c_chi;
     //const double E = std::sqrt( ( B * B ) + ( C * C ) );
     //const double b = rotationRateEarth * rotationRateEarth * currentAltitude * c_delta * c_delta - gravs( 1 );
     //const double a = -( rotationRateEarth * rotationRateEarth * currentAltitude * c_delta * s_delta + gravs( 0 ) ) * c_chi;
     //const double E = std::sqrt( ( a * a ) + ( b * b ) );
     //const double LHS = ( n_aero_xg * currentMass + currentDragForce - currentThrustMagnitude * c_alpha ) / ( E * currentMass );

     //const double arccosArgument = LHS / ( bislip::Variables::determineSignOfValue( a ) ) ;

     //const double cosArgument = std::acos( arccosArgument );
    // const double phi = std::atan2( B, C );
    //const double phi = std::atan( -b / a );

    //double flightPathAngle = cosArgument - phi;


    std::cout << "        Checking value of Flight-Path Angle" << std::endl;
    std::cout << "        ------------------------------------------------" << std::endl;
    std::cout << "              B1                      = " << B1 << std::endl;
    std::cout << "              C1                      = " << C1 << std::endl;
    std::cout << "              E1                      = " << E1 << std::endl;
    std::cout << "              phi1                    = " << phi1 << std::endl;
    std::cout << "              B2                      = " << B2 << std::endl;
    std::cout << "              C2                      = " << C2 << std::endl;
    std::cout << "              E2                      = " << E2 << std::endl;
    std::cout << "              phi2                    = " << phi2 << std::endl;
    std::cout << "              numerator_other_term    = " << numerator_other_term << std::endl;
    std::cout << "              numerator_sigma_term    = " << numerator_sigma_term << std::endl;


    std::cout << "              numerator               = " << numerator << std::endl;
    std::cout << "              denominator             = " << denominator << std::endl;
    std::cout << "              arg                     = " << arg << std::endl;
    std::cout << "              acosArg                 = " << acosArg << std::endl;
    std::cout << "              gamma                   = " << gamma << std::endl;
    std::cout << "              TminusD                 = " << TminusD << std::endl;
    std::cout << "              currentWeight           = " << currentWeight << std::endl;
    std::cout << "              sin_gamma               = " << sin_gamma << std::endl;
    std::cout << "              gamma1                  = " << gamma1 << std::endl;

   std::cout << "              a                      = " << a << std::endl;
    std::cout << "              b                      = " << b << std::endl;
    std::cout << "              E                      = " << E << std::endl;
    std::cout << "              LHS                    = " << LHS << std::endl;
    std::cout << "              n                      = " << n << std::endl;
    std::cout << "              Aero Frame g-Load - x  = " << aerodynamicFrameTotal_g_Load_x << std::endl;
    std::cout << "              n_aero_x * g           = " << n_aero_xg << std::endl;
    std::cout << "              n_aero_x * g * mass    = " << n_aero_xg * currentMass << std::endl;
    std::cout << "              T * c_alpha            = " << c_alpha * currentThrustMagnitude << std::endl;
    std::cout << "              dif                    = " << n_aero_xg * currentMass - c_alpha * currentThrustMagnitude << std::endl;
    std::cout << "              dif2                   = " << c_alpha * currentThrustMagnitude - currentDragForce<< std::endl;
    std::cout << "              currentWeight          = " << currentMass * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION<< std::endl;
    std::cout << "              sin(gamma)             = " << (c_alpha * currentThrustMagnitude - currentDragForce)/(currentMass * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION) << std::endl;
    std::cout << "              currentMass            = " << currentMass << std::endl;
    std::cout << "              currentDragForce       = " << currentDragForce << std::endl;
    std::cout << "              currentThrustMagnitude = " << currentThrustMagnitude << std::endl;
    std::cout << "              arccosArgument         = " << arccosArgument << std::endl;
    std::cout << "              cosArgument            = " << cosArgument << std::endl;
    std::cout << "              phi                    = " << phi << std::endl;
    std::cout << "              flightPathAngle        = " << flightPathAngle << std::endl;
    std::cout << "              currentFlightPathAngle = " << currentFlightPathAngle << std::endl;
*/




    return n;
}

double throttleSettingLimitSearchFunction(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName,
        const double &throttleSetting )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;


    double evaluation = bislipSystems->getMechanicalLoadConstraint() - bislip::Variables::throttleSettingEvaluationFunction( bodyMap, vehicleName, centralBodyName, throttleSetting );

    //std::cout << "      evaluation: " << evaluation << std::endl;

    return evaluation;
}

double determineThrustElevationAngle(
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
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( );

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 2 ){ std::cout << "        Starting Determination of Thrust Elevation Angle" << std::endl; }
    if( debugInfo == 2 ){ std::cout << "        ------------------------------------------------" << std::endl; }

    if( debugInfo == 2 ){ std::cout << "            Extract Relevant Conditions" << std::endl; }
    const std::string currentTrajectoryPhase = bislipSystems->getCurrentTrajectoryPhase();

    double currentFlightPathAngle;
    double currentLatitude;
    double currentHeading;
    double currentAltitude;
    double currentAirspeed;
    double currentDynamicPressure;
    double currentThrustElevationAngle;


    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 2 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentFlightPathAngle      = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialFlightPathAngle() );
        currentLatitude             = bislipSystems->getInitialLat();
        currentHeading              = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialHeadingAngle() );
        currentAltitude             = bislipSystems->getInitialAltitude();
        currentAirspeed             = bislipSystems->getInitialAirspeed();
        currentDynamicPressure      = bislipSystems->getInitialDynamicPressure();
        currentThrustElevationAngle = 0.0;
    }

    else
    {
        if( debugInfo == 2 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentFlightPathAngle      = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        currentLatitude             = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        currentHeading              = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
        currentAltitude             = flightConditions->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );
        currentAirspeed             = flightConditions->getCurrentAirspeed( );
        currentDynamicPressure      = flightConditions->getCurrentDynamicPressure();
        currentThrustElevationAngle = bislipSystems->getCurrentThrustElevationAngle();
    }


    double newThrustElevationAngle;

    const double currentMass                = bodyMap.at( vehicleName )->getBodyMass();
    const double currentAngleOfAttack       = bislipSystems->getCurrentAngleOfAttack();
    //const double currentFlightPathAngleRate = bislipSystems->getCurrentFlightPathAngleRate();
    const double currentBankAngle           = bislipSystems->getCurrentBankAngle();
    const double currentThrustMagnitude     = bislipSystems->getCurrentThrustMagnitude();
    const double currentThrustAzimuthAngle  = bislipSystems->getCurrentThrustAzimuthAngle();
    const double currentLiftForce           = bislipSystems->getCurrentLiftForce();
    const double currentDragForce           = bislipSystems->getCurrentDragForce();
    const Eigen::Vector3d gravs             = bislipSystems->getCurrentLocalGravityVector();
    const Eigen::Vector3d COM               = bislipSystems->getMassReferenceCenter();
    const Eigen::Vector3d COT               = bislipSystems->getThrustReferenceCenter();
    const Eigen::Vector3d MRC               = bislipSystems->getMomentReferenceCenter();
    const double rotationRateEarth          = 7.292115*1E-5;

    //const Eigen::Vector6d currentCoefficients = bislip::Variables::computePartialCurrentCoefficients( bodyMap, vehicleName );
    // if( debugInfo == 2 ){ std::cout << "            Compute Lift Force" << std::endl; }
    // const double currentLift = bislip::Variables::computeCurrentLiftForce( bodyMap, vehicleName );
    // if( debugInfo == 2 ){ std::cout << "            Compute Local Gravity Vector" << std::endl; }
    // const Eigen::Vector2d gravs = bislip::Variables::computeLocalGravity( bodyMap, vehicleName, centralBodyName );


    //const double currentThrottleSetting      = bislipSystems->getCurrentThrottleSetting();
    //tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::ThrustElevationAngle, bodyMap_, vehicleName_ ) ) );




    if( debugInfo == 2 ){ std::cout << "                Current Time                       = " << flightConditions->getCurrentTime() - bislipSystems->getStartingEpoch() << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Mass                       = " << currentMass << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Airspeed                   = " << currentAirspeed << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Dynamic Pressure           = " << currentDynamicPressure << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Angle of Attack            = " << tudat::unit_conversions::convertRadiansToDegrees( currentAngleOfAttack ) << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Lift Force                 = " << currentLiftForce << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Drag Force                 = " << currentDragForce << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Throttle Setting           = " << bislipSystems->getCurrentThrottleSetting() << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Thrust Magnitude           = " << currentThrustMagnitude << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Body-Fixed Total g-Load    = " << bislip::Variables::computeBodyFixedTotal_g_Load_Magnitude( bodyMap, vehicleName ) << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Bank Angle                 = " << currentBankAngle << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Local Gravity Vector               = [ " << gravs( 0 ) << ", " << gravs( 1 ) << ", " << gravs( 2 ) << " ]" << std::endl; }

    const double thrustElevationLowerBound = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getThrustElevationLimits() ).first );
    const double thrustElevationUpperBound = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getThrustElevationLimits() ).second );
    if( debugInfo == 2 ){ std::cout << "                Current Thrust Elevation Angle     = " << tudat::unit_conversions::convertRadiansToDegrees( currentThrustElevationAngle ) << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Thrust Elevation Angle Lower Bound = " << tudat::unit_conversions::convertRadiansToDegrees( thrustElevationLowerBound ) << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Thrust Elevation Angle Upper Bound = " << tudat::unit_conversions::convertRadiansToDegrees( thrustElevationUpperBound ) << std::endl; }


    if( debugInfo == 2 ){ std::cout << "           Create Local Function to Calculate the Flight-Path Angle Rate Given a Thrust Elevation Angle" << std::endl; }
    std::function< double( const double ) > thrustElevationAngleEvaluation =
            std::bind( &bislip::Variables::thrustElevationAngleEvaluationFunction, bodyMap, vehicleName, centralBodyName, std::placeholders::_1 );

    if( debugInfo == 2 ){ std::cout << "           Evaluation of Thrust Elevation Angles" << std::endl; }
    if( debugInfo == 2 ){ std::cout << "           -----------------------------------------------" << std::endl; }
    const double currentFlightPathAngleRate = thrustElevationAngleEvaluation( currentThrustElevationAngle );
    if( debugInfo == 2 ){ std::cout << "                Current Flight-Path Angle          = " << currentFlightPathAngle << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Flight-Path Angle Rate     = " << currentFlightPathAngleRate << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                Current Trajectory Phase           = " << currentTrajectoryPhase << std::endl; }


    if( currentTrajectoryPhase == "Ascent" && currentFlightPathAngle <= 0 )
    {


        /*const double s_alpha = std::sin( currentAngleOfAttack );
    const double c_alpha = std::cos( currentAngleOfAttack );
    const double c_sigma = std::cos( currentBankAngle );
    const double s_delta = std::sin( currentLatitude );
    const double c_delta = std::cos( currentLatitude );
    const double s_gamma = std::sin( currentFlightPathAngle );
    const double c_gamma = std::cos( currentFlightPathAngle );
    const double s_chi = std::sin( currentHeading );
    const double c_chi = std::cos( currentHeading );
    const double s_thrustAzimuthAngle = std::sin( currentThrustAzimuthAngle );
    const double c_thrustAzimuthAngle = std::cos( currentThrustAzimuthAngle );
    //  const double s_thrustElevationAngle = std::sin( currentThrustElevationAngle );
    //const double c_thrustElevationAngle = std::cos( currentThrustElevationAngle );


*/




        const double lowerThrustElevationAngleFlightPathAngleRate = thrustElevationAngleEvaluation( thrustElevationLowerBound );
        const double upperThrustElevationAngleFlightPathAngleRate = thrustElevationAngleEvaluation( thrustElevationUpperBound );


        if( debugInfo == 2 ){ std::cout << "                Lowest Thrust Elevation Angle Flight-Path Angle Rate  = " << lowerThrustElevationAngleFlightPathAngleRate << std::endl; }
        if( debugInfo == 2 ){ std::cout << "                Highest Thrust Elevation Angle Flight-Path Angle Rate = " << upperThrustElevationAngleFlightPathAngleRate << std::endl; }

        std::function< double( const double ) > thrustElevationAngleLimitSearchFunction =
                std::bind( &bislip::Variables::thrustElevationAngleEvaluationFunction, bodyMap, vehicleName, centralBodyName, std::placeholders::_1 );

        double f_a = thrustElevationAngleLimitSearchFunction( thrustElevationLowerBound );
        double f_b = thrustElevationAngleLimitSearchFunction( thrustElevationUpperBound );
        if( debugInfo == 2 ){ std::cout << "                f( " << thrustElevationLowerBound << " ) = " << f_a << std::endl; }
        if( debugInfo == 2 ){ std::cout << "                f( " << thrustElevationUpperBound << " ) = " << f_b << std::endl; }

        double thrustElevationAngleForZeroFlightPathAngle;
        double newFlightPathAngle;



        if( f_a * f_b < 0 )
        {
            if( debugInfo == 2 ){ std::cout << "                Bisection Search" << std::endl; }
            thrustElevationAngleForZeroFlightPathAngle = bislip::Variables::rootFinderBisection( thrustElevationAngleLimitSearchFunction, thrustElevationLowerBound, thrustElevationUpperBound, currentThrustElevationAngle );
            if( debugInfo == 2 ){ std::cout << "                Bisection Search Complete" << std::endl; }
            newFlightPathAngle = thrustElevationAngleEvaluation( thrustElevationAngleForZeroFlightPathAngle );
        }
        else
        {

            if( std::abs( f_a ) < std::abs( f_b ) && std::abs( f_a ) <= std::abs( currentFlightPathAngleRate ) )
            {
                if( debugInfo == 2 ){ std::cout << "                Lower Bound is closer to zero" << std::endl; }

                thrustElevationAngleForZeroFlightPathAngle = thrustElevationLowerBound;
                newFlightPathAngle = thrustElevationAngleEvaluation( thrustElevationAngleForZeroFlightPathAngle );
            }

            else if( std::abs( f_b ) < std::abs( f_a ) && std::abs( f_b ) <= std::abs( currentFlightPathAngleRate ) )
            {
                if( debugInfo == 2 ){ std::cout << "                Upper Bound is closer to zero" << std::endl; }

                thrustElevationAngleForZeroFlightPathAngle = thrustElevationUpperBound;
                newFlightPathAngle = thrustElevationAngleEvaluation( thrustElevationAngleForZeroFlightPathAngle );
            }
            else
            {
                if( debugInfo == 2 ){ std::cout << "                Golden Ratio Search" << std::endl; }
                thrustElevationAngleForZeroFlightPathAngle = bislip::Variables::goldenSectionSearch( thrustElevationAngleLimitSearchFunction, thrustElevationLowerBound, thrustElevationUpperBound );
                if( debugInfo == 2 ){ std::cout << "                Golden Ratio Search Complete" << std::endl; }
                newFlightPathAngle = thrustElevationAngleEvaluation( thrustElevationAngleForZeroFlightPathAngle );
            }
        }

        newThrustElevationAngle = thrustElevationAngleForZeroFlightPathAngle;


        /*
                       const double A1 = 2.0 * rotationRateEarth * currentAirspeed * c_delta * s_chi;
                       const double A2 = ( currentAirspeed * currentAirspeed / currentAltitude ) * c_gamma;
                       const double A3 = rotationRateEarth * rotationRateEarth * currentAltitude * c_delta * ( c_delta * c_gamma + s_gamma * s_delta * c_chi );

                       const double A = A1 + A2 + A3;
                       const double B = currentThrustMagnitude * c_alpha;
    const double C = currentThrustMagnitude * s_alpha * c_thrustAzimuthAngle;
    const double D = gravs( 1 ) * c_gamma - gravs( 0 ) * s_gamma * c_chi;
    const double E = std::sqrt( ( B * B ) + ( C * C ) );

    const double arccosArgument = -( 1 / E ) * ( currentLiftForce + currentMass * ( ( A - D ) / c_sigma ) );
    double arccos = 0;

    if( std::abs( arccosArgument ) >= 1.0 ) { arccos = 0.0; }
    else if( std::isnan( arccosArgument ) == false ){ arccos = std::acos( arccosArgument ); }

    double currentThrustElevationAngle = arccos + std::atan2( B, C );
*/
        /*
    const double a = COT( 0 );
    const double b = COT( 2 );
    const double c = std::sqrt( ( a * a ) + ( b * b ) );
    const double LHS = COM( 0 ) * ( currentDragForce * s_alpha - currentLiftForce * c_alpha ) / currentThrustMagnitude;
    const double arccosArgument = LHS / ( bislip::Variables::determineSignOfValue( a ) * c );
    double arccos = 0.0;
    double phi = 0.0;


    if( ( a > 0 ) && ( b == 0.0 ) ) { phi = tudat::mathematical_constants::PI / 2; }
    else if( ( a < 0 ) && ( b == 0.0 ) ) { phi = -tudat::mathematical_constants::PI / 2; }
    else { phi = std::atan( a / b ); }

    if( arccosArgument >= 1.0 ) { arccos = 0.0; }
    else if( arccosArgument <= -1.0 ) { arccos = tudat::mathematical_constants::PI; }
    else if( std::isnan( arccosArgument ) == false ){ arccos = std::acos( arccosArgument ); }

    double currentThrustElevationAngle = std::acos( arccosArgument ) - phi;


    double x1 = 2 * std::atan2( ( b + std::sqrt( ( a * a ) + ( b * b ) - ( LHS * LHS ) ) ) , ( a + LHS ) );
    double x2 = 2 * std::atan2( ( b - std::sqrt( ( a * a ) + ( b * b ) - ( LHS * LHS ) ) ) , ( a + LHS ) );


*/

        /* if( debugInfo == 2 ){ std::cout << "                a                                 = " << a << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                b                                 = " << b << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                c                                 = " << c << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                LHS                               = " << LHS << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                arccosArgument                    = " << arccosArgument << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                arccos                            = " << arccos << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                phi                               = " << phi << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                std::atan2( b, a )                = " << std::atan2( b, a ) << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                x1                                = " << x1 << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                x2                                = " << x2 << std::endl; }
*/
        // if( debugInfo == 2 ){ std::cout << "                currentThrustElevationAngle       = " << currentThrustElevationAngle << " rad" << std::endl; }


        if( debugInfo == 2 ){ std::cout << "                Calculated Thrust Elevation Angle = " << tudat::unit_conversions::convertRadiansToDegrees( thrustElevationAngleForZeroFlightPathAngle ) << " deg" << std::endl; }
        if( debugInfo == 2 ){ std::cout << "                Calculated Flight-Path Angle Rate = " << newFlightPathAngle << " rad/s" << std::endl; }
        // if( debugInfo == 2 ){ std::cout << "                    Extract Thrust Elevation Angle Bounds" << std::endl; }

        /* const double lowerbound = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getThrustElevationLimits() ).first );
    const double upperbound = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getThrustElevationLimits() ).second );
    if( debugInfo == 2 ){ std::cout << "                        lowerbound                = " << lowerbound << std::endl; }
    if( debugInfo == 2 ){ std::cout << "                        upperbound                = " << upperbound << std::endl; }

    if( debugInfo == 2 ){ std::cout << "                    Impose Thrust Elevation Angle Bounds" << std::endl; }
    if( currentThrustElevationAngle < lowerbound ) { currentThrustElevationAngle = lowerbound; }
    if( currentThrustElevationAngle > upperbound ) { currentThrustElevationAngle = upperbound; }

    if( debugInfo == 2 ){ std::cout << "                Bound Elevation Angle             = " <<  tudat::unit_conversions::convertRadiansToDegrees( currentThrustElevationAngle ) << std::endl; }

*/
        /*
    std::cout << "A =                           " << A << std::endl;
    std::cout << "B =                           " << B << std::endl;
    std::cout << "C =                           " << C << std::endl;
    */

    }
    else
    {
        newThrustElevationAngle = currentThrustElevationAngle;

    }

    if( debugInfo == 2 ){ std::cout << "        Ending Determination of Thrust Elevation Angle" << std::endl; }
    if( debugInfo == 2 ){ std::cout << "        ------------------------------------------------" << std::endl; }

    return newThrustElevationAngle;
}


double thrustElevationAngleEvaluationFunction(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName,
        const double &thrustElevationAngle )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( );


    bislipSystems->setCurrentThrustElevationAngle( thrustElevationAngle );

    const double flightPathAngleRate = bislip::Variables::computeFlightPathAngleRate( bodyMap, vehicleName, centralBodyName );

    return flightPathAngleRate;
}
/*
double thrustElevationAngleSearchFunction(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::string &centralBodyName,
        const double &thrustElevationAngle )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( );


    bislipSystems->setCurrentThrustElevationAngle( thrustElevationAngle );

    const double calculatedFlightPathAngleRate = bislip::Variables::computeFlightPathAngleRate( bodyMap, vehicleName, centralBodyName );



    return flightPathAngleRate;
}
*/





double computeSkipSuppressionLimit(
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
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( );

    int debugInfo = bislipSystems->getDebugInfo();
    if( debugInfo == 1 ){ std::cout << "        Starting Determination of Skip Suppression Limit" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "        ------------------------------------------------" << std::endl; }

    double currentFlightPathAngle;
    double currentLatitude;
    double currentHeading;
    double currentAltitude;
    double currentAirspeed;
    double currentDynamicPressure;

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentFlightPathAngle = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialFlightPathAngle() );
        currentLatitude        = bislipSystems->getInitialLat();
        currentHeading         = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialHeadingAngle() );
        currentAltitude        = bislipSystems->getInitialAltitude();
        currentAirspeed        = bislipSystems->getInitialAirspeed();
        currentDynamicPressure = bislipSystems->getInitialDynamicPressure();
    }

    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentFlightPathAngle = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        currentLatitude        = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        currentHeading         = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
        currentAltitude        = flightConditions->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );
        currentAirspeed        = flightConditions->getCurrentAirspeed( );
        currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    }


    //const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    const double currentMass            = bodyMap.at( vehicleName )->getBodyMass();

    const double rotationRateEarth = 7.292115*1E-5;

    //const Eigen::Vector6d currentCoefficients = bislip::Variables::computePartialCurrentCoefficients( bodyMap, vehicleName );
    //const double currentLift = bislip::Variables::computeCurrentLiftForce( bodyMap, vehicleName );
    const double currentLift = bislipSystems->getCurrentLiftForce();

    if( debugInfo == 1 ){ std::cout << "            Extract Relevant Conditions" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Lift                 = " << currentLift << std::endl; }

    //const Eigen::Vector2d gravs = bislip::Variables::computeLocalGravity( bodyMap, vehicleName, centralBodyName );
    const Eigen::Vector3d gravs = bislipSystems->getCurrentLocalGravityVector();
    if( debugInfo == 1 ){ std::cout << "                Current Local Gravity Vector = [ " << gravs( 0 ) << ", " << gravs( 1 ) << ", " << gravs( 2 ) << " ]" << std::endl; }



    const double currentAngleOfAttack        = bislipSystems->getCurrentAngleOfAttack();
    const double currentThrustMagnitude      = bislipSystems->getCurrentThrustMagnitude();
    const double currentThrustElevationAngle = bislipSystems->getCurrentThrustElevationAngle();
    const double currentThrustAzimuthAngle   = bislipSystems->getCurrentThrustAzimuthAngle();


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
    const double D = ( gravs( 0 ) * s_gamma * c_chi - gravs( 2 ) * c_gamma );
    const double E = std::sqrt( ( B * B ) + ( C * C ) );
    const double arccosArgument = -( D + A ) * ( currentMass / E );
    double cosArgument = 0;

    if( arccosArgument >= 1.0 ) { cosArgument = 0; }
    else if( arccosArgument <= -1.0 ) { cosArgument = tudat::mathematical_constants::PI; }
    else { cosArgument = std::acos( arccosArgument ); }

    //const double argument = -( currentMass / c ) * ( term1 + term2 + term3 - term4 + term5 );
    double limit = cosArgument + std::atan2( B, C );

    //if( limit < 0.0 ) { limit = std::abs( limit ); }
    if( limit > tudat::mathematical_constants::PI / 2 ) { limit = tudat::mathematical_constants::PI / 2; }
    if( limit < -tudat::mathematical_constants::PI / 2 ) { limit = -tudat::mathematical_constants::PI / 2; }


    if( debugInfo == 1 ){ std::cout << "                Skip Suppression Limit = " << limit << std::endl; }


    /*
    if( std::abs( arccosArgument ) > 1.0 )
    {
        std::cout << "      M   = " << currentMass << std::endl;
        std::cout << "      L   = " << currentLift << std::endl;
        std::cout << "      T   = " << currentThrustMagnitude << std::endl;
        std::cout << "      A   = " << A << std::endl;
        std::cout << "      B   = " << B << std::endl;
        std::cout << "      C   = " << C << std::endl;
        std::cout << "      D   = " << D << std::endl;
        std::cout << "      E   = " << E << std::endl;
        std::cout << "      ARG = " << arccosArgument << std::endl;
        std::cout << "      limit = " << limit << std::endl;

        //   limit = 0.0;

    }
*/
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
    return std::abs( limit );
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

    int debugInfo = bislipSystems->getDebugInfo();
    if( debugInfo == 1 ){ std::cout << "        Starting Determination of Flight-Path Angle Rate" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "        ------------------------------------------------" << std::endl; }

    double currentFlightPathAngle;
    double currentLatitude;
    double currentHeading;
    double currentAltitude;
    double currentAirspeed;
    double currentDynamicPressure;

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentFlightPathAngle = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialFlightPathAngle() );
        currentLatitude        = bislipSystems->getInitialLat();
        currentHeading         = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialHeadingAngle() );
        currentAltitude        = bislipSystems->getInitialAltitude();
        currentAirspeed        = bislipSystems->getInitialAirspeed();
        currentDynamicPressure = bislipSystems->getInitialDynamicPressure();
        //currentLift            = bislipSystems->getCurrentLiftForce();
    }

    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentFlightPathAngle = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        currentLatitude        = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        currentHeading         = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
        currentAltitude        = flightConditions->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );
        currentAirspeed        = flightConditions->getCurrentAirspeed( );
        currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
        //currentLift            = bislip::Variables::computeCurrentLiftForce( bodyMap, vehicleName );

    }

    //const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    const double currentMass       = bodyMap.at( vehicleName )->getBodyMass();

    const double rotationRateEarth = 7.292115*1E-5;

    //const Eigen::Vector6d currentCoefficients = bislip::Variables::computePartialCurrentCoefficients( bodyMap, vehicleName );
    //const double currentLift = bislip::Variables::computeCurrentLiftForce( bodyMap, vehicleName );
    const double currentLift = bislipSystems->getCurrentLiftForce();

    //const Eigen::Vector2d gravs = bislip::Variables::computeLocalGravity( bodyMap, vehicleName, centralBodyName );
    const Eigen::Vector3d gravs = bislipSystems->getCurrentLocalGravityVector();

    const double currentAngleOfAttack        = bislipSystems->getCurrentAngleOfAttack();
    const double currentBankAngle            = bislipSystems->getCurrentBankAngle();
    const double currentThrustMagnitude      = bislipSystems->getCurrentThrustMagnitude();
    const double currentThrustElevationAngle = bislipSystems->getCurrentThrustElevationAngle();
    const double currentThrustAzimuthAngle   = bislipSystems->getCurrentThrustAzimuthAngle();

    if( debugInfo == 1 ){ std::cout << "                Current Flight-Path Angle      = " << currentFlightPathAngle << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Latitude               = " << tudat::unit_conversions::convertRadiansToDegrees( currentLatitude ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Heading                = " << tudat::unit_conversions::convertRadiansToDegrees( currentHeading ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Altitude               = " << currentAltitude << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Airspeed               = " << currentAirspeed << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Dynamic Pressure       = " << currentDynamicPressure << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Angle of Attack        = " << tudat::unit_conversions::convertRadiansToDegrees( currentAngleOfAttack ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Bank Angle             = " << tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Thrust Magnitude       = " << currentThrustMagnitude << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Thrust Elevation Angle = " << tudat::unit_conversions::convertRadiansToDegrees( currentThrustElevationAngle ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Thrust Azimuth Angle   = " << tudat::unit_conversions::convertRadiansToDegrees( currentThrustAzimuthAngle ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Lift Force             = " << currentLift << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Mass                   = " << currentMass << std::endl; }



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

    const double F_gamma_1 = currentThrustMagnitude * ( s_thrustAzimuthAngle * c_thrustElevationAngle ) * s_sigma;
    const double F_gamma_2 = ( currentLift + currentThrustMagnitude * ( s_alpha * c_thrustAzimuthAngle * c_thrustElevationAngle  + c_alpha * s_thrustElevationAngle ) ) * c_sigma;
    const double F_gamma_3 = currentMass * ( gravs( 0 ) * s_gamma * c_chi - gravs( 2 ) * c_gamma );

    const double F_gamma = F_gamma_1 + F_gamma_2 + F_gamma_3;

    const double A1 = 2.0 * rotationRateEarth * currentAirspeed * c_delta * s_chi;
    const double A2 = ( currentAirspeed * currentAirspeed / currentAltitude ) * c_gamma;
    const double A3 = rotationRateEarth * rotationRateEarth * currentAltitude * c_delta * ( c_delta * c_gamma + s_gamma * s_delta * c_chi );

    const double A = A1 + A2 + A3;

    const double flightPathAngleRate = ( ( F_gamma / currentMass ) + A ) / currentAirspeed;
    if( debugInfo == 1 ){ std::cout << "                Flight-Path Angle Rate         = " << flightPathAngleRate << std::endl; }

    return flightPathAngleRate;
}



double computeCumulativeCartesianDistanceTravelled(
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



double computeCumulativeAngularDistanceTravelled(
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


double computeCumulativeAngularDistanceTravelledDifference(
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

double computeAngularDistanceCoveredRatio(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    const double currentLatitude  = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double currentLongitude = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );

    double netAngularDisplacement = bislip::Variables::computeAngularDistance( bislipSystems->getInitialLat(), bislipSystems->getInitialLon(),
                                                                               currentLatitude, currentLongitude );

    double angularDistanceCoveredRatio = netAngularDisplacement / bislipSystems->getInitialDistanceToTarget() ;

    if( debugInfo == 1 ){ std::cout << "Angular Distance Covered       = " << netAngularDisplacement << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Angular Distance Covered Ratio = " << angularDistanceCoveredRatio << std::endl; }

    return angularDistanceCoveredRatio;
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
    double currentMachNumber;

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentMachNumber = bislipSystems->getInitialMachNumber();
    }

    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentMachNumber = flightConditions->getCurrentMachNumber();
    }

    Eigen::Vector2d controlSurfaceDeflections( 2 );
    Eigen::VectorXd CmBound ( 2 );

    if( debugInfo == 1 ){ std::cout << "Create primary input for aerodynamic coefficient interface" << std::endl; }
    //! Create primary input for aerodynamic coefficient interface.
    //!     Take care of order of input (this depends on how the coefficients are created)!
    const std::vector < double > aerodynamicCoefficientInput = bislip::Variables::getAerodynamicCoefficientInput( currentAngleOfAttack, currentMachNumber );

    if( debugInfo == 1 ){ std::cout << "Create current body-flap pitch moment function" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    Current Angle of Attack = " << currentAngleOfAttack << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    Current Mach Number     = " << currentMachNumber << std::endl; }
    std::function< double( const double ) > bodyFlapPitchMomentFunction =
            std::bind( &bislip::Variables::computeFullPitchMomentCoefficient, coefficientInterface, aerodynamicCoefficientInput, std::placeholders::_1, 0.0 );

    double delta_bf_min = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getBodyFlapDeflectionLimits() ).first );
    double current_delta_bf = vehicleSystems->getCurrentControlSurfaceDeflection( "BodyFlap" );
    double delta_bf_max = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getBodyFlapDeflectionLimits() ).second );

    CmBound( 0 ) = bodyFlapPitchMomentFunction( delta_bf_min );
    double currentCm = bodyFlapPitchMomentFunction( current_delta_bf );
    CmBound( 1 ) = bodyFlapPitchMomentFunction( delta_bf_max );

    if( debugInfo == 1 ){ std::cout << "     C_m( "<< ( bislipSystems->getBodyFlapDeflectionLimits() ).first <<", 0.0 ) = " << CmBound( 0 ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     C_m( "<<  tudat::unit_conversions::convertRadiansToDegrees( vehicleSystems->getCurrentControlSurfaceDeflection( "BodyFlap" ) ) <<" , 0.0 ) = " << currentCm << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     C_m( "<< ( bislipSystems->getBodyFlapDeflectionLimits() ).second <<" , 0.0 ) = " << CmBound( 1 ) << std::endl; }

    if( debugInfo == 1 ){ std::cout << "     Starting Control Surface Deflection Angle Search" << std::endl; }

    double delta_bf = 0.0;
    if( CmBound( 0 ) * CmBound( 1 ) < 0 )
    {
        if( debugInfo == 1 ){ std::cout << "        Starting BodyFlap Angle search with Bisection Method" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "            Current BodyFlap Angle = " << tudat::unit_conversions::convertRadiansToDegrees( current_delta_bf ) << "     |     " << "C_m( "<<  tudat::unit_conversions::convertRadiansToDegrees( current_delta_bf ) << ", " << 0.0 << " ) = " << std::abs( bodyFlapPitchMomentFunction( current_delta_bf ) ) << std::endl; }

        delta_bf = bislip::Variables::rootFinderBisection( bodyFlapPitchMomentFunction, delta_bf_min, delta_bf_max, current_delta_bf );

        if( debugInfo == 1 ){ std::cout << "            Root Found             = " << tudat::unit_conversions::convertRadiansToDegrees( delta_bf ) << "     |     " << "C_m( "<<  tudat::unit_conversions::convertRadiansToDegrees( delta_bf ) << ", " <<  0.0 << " ) = " << std::abs( bodyFlapPitchMomentFunction( delta_bf ) ) << std::endl; }

        controlSurfaceDeflections << delta_bf, 0.0 ;
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "        Current AoA-Mach combination requires both body-flap and elevons for pitch trim." << std::endl; }

        if( debugInfo == 1 ){ std::cout << "            Starting BodyFlap angle search that minimizes pitch moment with Golden Search" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "                Current BodyFlap Angle = " << tudat::unit_conversions::convertRadiansToDegrees( current_delta_bf ) << "     |     " << "C_m( "<<  tudat::unit_conversions::convertRadiansToDegrees( current_delta_bf ) << ", " << 0.0 << " ) = " << std::abs( bodyFlapPitchMomentFunction( current_delta_bf ) ) << std::endl; }

        delta_bf = bislip::Variables::goldenSectionSearch( bodyFlapPitchMomentFunction, delta_bf_min, delta_bf_max );

        if( debugInfo == 1 ){ std::cout << "                BodyFlap Angle Found   = " << tudat::unit_conversions::convertRadiansToDegrees( delta_bf ) << "     |     " << "C_m( "<<  tudat::unit_conversions::convertRadiansToDegrees( delta_bf ) << ", " <<  0.0 << " ) = " << std::abs( bodyFlapPitchMomentFunction( delta_bf ) ) << std::endl; }

        if( debugInfo == 1 ){ std::cout << "                Create current elevon pitch moment function" << std::endl; }
        std::function< double( const double ) > elevonPitchMomentFunction =
                std::bind( &bislip::Variables::computeFullPitchMomentCoefficient, coefficientInterface, aerodynamicCoefficientInput, delta_bf, std::placeholders::_1 );

        controlSurfaceDeflections( 0 ) = delta_bf;

        const double delta_el_min = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getElevonDeflectionLimits() ).first );
        double current_delta_el = vehicleSystems->getCurrentControlSurfaceDeflection( "ElevonRight" );
        const double delta_el_max = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getElevonDeflectionLimits() ).second );

        CmBound( 0 ) = elevonPitchMomentFunction( delta_el_min );
        currentCm = elevonPitchMomentFunction( current_delta_bf );
        CmBound( 1 ) = elevonPitchMomentFunction( delta_el_max );

        if( debugInfo == 1 ){ std::cout << "                    C_m( "<<  tudat::unit_conversions::convertRadiansToDegrees( delta_bf ) <<", "<<  tudat::unit_conversions::convertRadiansToDegrees( delta_el_min ) <<" ) = " << CmBound( 0 ) << std::endl; }
        if( debugInfo == 1 ){ std::cout << "                    C_m( "<<  tudat::unit_conversions::convertRadiansToDegrees( delta_bf ) <<", "<<  tudat::unit_conversions::convertRadiansToDegrees( current_delta_el ) <<" ) = " << currentCm << std::endl; }
        if( debugInfo == 1 ){ std::cout << "                    C_m( "<<  tudat::unit_conversions::convertRadiansToDegrees( delta_bf ) <<", "<<  tudat::unit_conversions::convertRadiansToDegrees( delta_el_max ) <<" ) = " << CmBound( 1 ) << std::endl; }

        if( CmBound( 0 ) * CmBound( 1 ) < 0 )
        {
            if( debugInfo == 1 ){ std::cout << "                    Starting elevon search:" << std::endl; }

            if( debugInfo == 1 ){ std::cout << "                    Current Elevon Angle = " << tudat::unit_conversions::convertRadiansToDegrees( current_delta_el ) << "     |     " << "C_m( "<<  tudat::unit_conversions::convertRadiansToDegrees( delta_bf ) << ", " <<  tudat::unit_conversions::convertRadiansToDegrees( current_delta_el ) << " ) = " << std::abs( elevonPitchMomentFunction( current_delta_el ) ) << std::endl; }

            double delta_el = bislip::Variables::rootFinderBisection( elevonPitchMomentFunction, delta_el_min, delta_el_max, current_delta_el );

            if( debugInfo == 1 ){ std::cout << "                    Root Found           = " << tudat::unit_conversions::convertRadiansToDegrees( delta_el ) << "     |     " << "C_m( "<<  tudat::unit_conversions::convertRadiansToDegrees( delta_bf ) << ", " <<  tudat::unit_conversions::convertRadiansToDegrees( delta_el ) << " ) = " << std::abs( elevonPitchMomentFunction( delta_el ) ) << std::endl; }

            controlSurfaceDeflections( 1 ) = delta_el;
        }
    }

    if( debugInfo == 1 ){ std::cout << "PA FUERA!" << std::endl; }

    return controlSurfaceDeflections;
}


double computeFullPitchMomentCoefficient(
        const std::shared_ptr< tudat::aerodynamics::AerodynamicCoefficientInterface > &coefficientInterface,
        const std::vector< double > &aerodynamicCoefficientInput,
        const double &bodyFlapDeflection,
        const double &elevonDeflection )
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

    if( debugInfo == 1 ){ std::cout << "    Starting Computation of Full Current Coefficients" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    -------------------------------------------------" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "        Extracting current data" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Selecting source of Mach Number" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Time             = " << flightConditions->getCurrentTime() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Starting Time            = " << bislipSystems->getStartingEpoch() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Guidance Step            = " << bislipSystems->getGuidanceStepSize() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Initial Mass             = " << bislipSystems->getInitialMass() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Mass             = " << bodyMap.at( vehicleName )->getBodyMass() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Initial Airspeed         = " << bislipSystems->getInitialAirspeed() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Initial Speed of Sournd  = " << bislipSystems->getInitialSpeedOfSound() << std::endl; }

    double currentMachNumber = 0.0;

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Mach" << std::endl; }

        currentMachNumber = bislipSystems->getInitialMachNumber();
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Mach" << std::endl; }

        currentMachNumber = flightConditions->getCurrentMachNumber();
    }

    const double currentAngleOfAttack           = bislipSystems->getCurrentAngleOfAttack();
    const double currentBodyFlapAngleDeflection = vehicleSystems->getCurrentControlSurfaceDeflection( "BodyFlap" );
    const double currentElevonAngleDeflection   = vehicleSystems->getCurrentControlSurfaceDeflection( "ElevonLeft" );
    if( debugInfo == 1 ){ std::cout << "            Current Angle Of Attack           = " << tudat::unit_conversions::convertRadiansToDegrees( currentAngleOfAttack ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Current Mach Number               = " << currentMachNumber << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Current BodyFlap Angle Deflection = " << tudat::unit_conversions::convertRadiansToDegrees( currentBodyFlapAngleDeflection ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Current Elevon Angle Deflection   = " << tudat::unit_conversions::convertRadiansToDegrees( currentElevonAngleDeflection ) << std::endl; }

    if( debugInfo == 1 ){ std::cout << "        Create primary input for aerodynamic coefficient interface." << std::endl; }
    //! Create primary input for aerodynamic coefficient interface.
    //!     Take care of order of input (this depends on how the coefficients are created)!
    std::vector< double > aerodynamicCoefficientInput = bislip::Variables::getAerodynamicCoefficientInput( currentAngleOfAttack, currentMachNumber );

    if( debugInfo == 1 ){ std::cout << "        Create control surface input for aerodynamic coefficient interface." << std::endl; }
    //! Create control surface input for aerodynamic coefficient interface.
    //!     Take care of order of input (this depends on how the coefficients are created)!
    std::map< std::string, std::vector< double > > controlSurfaceCoefficientInput = bislip::Variables::getControlSurfaceCoefficientInput( currentAngleOfAttack, currentMachNumber, currentBodyFlapAngleDeflection, currentElevonAngleDeflection );

    if( debugInfo == 1 ){ std::cout << "        Update and retrieve current aerodynamic coefficients." << std::endl; }
    coefficientInterface->updateFullCurrentCoefficients( aerodynamicCoefficientInput, controlSurfaceCoefficientInput );
    if( debugInfo == 1 ){ std::cout << "    Ending Computation of Full Current Coefficients" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    -------------------------------------------------" << std::endl; }
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
                    bislip::Parameters::Interpolators::ThrustElevationAngle,
                    bodyMap,
                    vehicleName ) );*/
    if( debugInfo == 1 ){ std::cout << "    computeBodyFlapCmIncrement -----> Extracting current thrust azimuth angle" << std::endl; }

    double thrustAzimuthAngle = bislipSystems->getCurrentThrustAzimuthAngle();


    /*tudat::unit_conversions::convertDegreesToRadians(
                bislip::Variables::evaluateGuidanceInterpolator (
                    bislip::Parameters::Interpolators::ThrustAzimuthAngle,
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

    if( debugInfo == 1 ){ std::cout << "    Starting Computation of Body-Fixed Total Load" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    ---------------------------------------------------" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "        Calculate the Body-Fixed Aerodynamic Load" << std::endl; }
    Eigen::Vector3d bodyfixedAerodynamicLoad = bislip::Variables::computeBodyFixedAerodynamicLoad( bodyMap, vehicleName );
    if( debugInfo == 1 ){ std::cout << "            Body-fixed Aerodynamic Load = [ " << bodyfixedAerodynamicLoad(0) << " , " << bodyfixedAerodynamicLoad(1) << " , " << bodyfixedAerodynamicLoad(2) << " ]" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "        Calculate the Body-Fixed Thrust Load" << std::endl; }
    Eigen::Vector3d bodyfixedThrustLoad = bislip::Variables::computeBodyFixedThrustVector( bodyMap, vehicleName );
    if( debugInfo == 1 ){ std::cout << "            Body-fixed Thrust Load      = [ " << bodyfixedThrustLoad(0) << " , " << bodyfixedThrustLoad(1) << " , " << bodyfixedThrustLoad(2) << " ]" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "        Calculate the Total Body-Fixed Load" << std::endl; }
    Eigen::Vector3d totalLoadBodyFrame = bodyfixedAerodynamicLoad + bodyfixedThrustLoad;
    if( debugInfo == 1 ){ std::cout << "            Total Body-fixed Load       = [ " << totalLoadBodyFrame(0) << " , " << totalLoadBodyFrame(1) << " , " << totalLoadBodyFrame(2) << " ]" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    Ending Computation of Body-Fixed Total Load" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    ---------------------------------------------------" << std::endl; }

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

    if( debugInfo == 1 ){ std::cout << "            Starting Computation of Body-Fixed Aerodynamic Load" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            ---------------------------------------------------" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Extracting current data" << std::endl; }
    const double currentAngleOfAttack   = bislipSystems->getCurrentAngleOfAttack();

    if( debugInfo == 1 ){ std::cout << "                    Selecting Source of Dynamic Pressure" << std::endl; }
    double currentDynamicPressure = 0.0;

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "                    Selecting Initial Values" << std::endl; }

        currentDynamicPressure = bislipSystems->getInitialDynamicPressure();
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "                    Selecting Propagated Values" << std::endl; }

        currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    }

    //const double airDensity    = flightConditions->getCurrentDensity();
    //const double airSpeed      = flightConditions->getCurrentAirspeed();
    //const double machNumber    = flightConditions->getCurrentMachNumber();
    if( debugInfo == 1 ){ std::cout << "                currentAngleOfAttack   = " << currentAngleOfAttack << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                currentDynamicPressure = " << currentDynamicPressure << std::endl; }

    //! Extract information from Bislip/Vehicle Systems.
    const double referenceArea = bislipSystems->getReferenceValues( )( 0 );
    if( debugInfo == 1 ){ std::cout << "                reference area         = " << referenceArea << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Calculate Aerodynamic Load in the Aerodynamic Frame" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                ---------------------------------------------------" << std::endl; }
    Eigen::Vector3d aerodynamicLoad = bislip::Variables::computeAerodynamicFrameAerodynamicLoad( bodyMap, vehicleName ) ;
    if( debugInfo == 1 ){ std::cout << "                    Aerodynamic Load in the Aerodynamic Frame = [ " << aerodynamicLoad(0) << " , " << aerodynamicLoad(1) << " , " << aerodynamicLoad(2) << " ]" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                ---------------------------------------------------" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Transform aerodynamic load from the Aerodynamic Frame to the Body Frame" << std::endl; }
    Eigen::Vector3d aerodynamicLoadBodyFrame = tudat::reference_frames::getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix( currentAngleOfAttack, 0.0 ) * aerodynamicLoad ;
    if( debugInfo == 1 ){ std::cout << "                    Body-fixed Aerodynamic Load = [ " << aerodynamicLoadBodyFrame(0) << " , " << aerodynamicLoadBodyFrame(1) << " , " << aerodynamicLoadBodyFrame(2) << " ]" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            Ending Computation of Body-Fixed Aerodynamic Load" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            ---------------------------------------------------" << std::endl; }

    return aerodynamicLoadBodyFrame ;
}

Eigen::Vector3d computeAerodynamicFrameAerodynamicLoad(
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

    if( debugInfo == 1 ){ std::cout << "            Starting Computation of Aerodynamic Frame Aerodynamic Load" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            ---------------------------------------------------" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Extracting current data" << std::endl; }
    const double currentAngleOfAttack   = bislipSystems->getCurrentAngleOfAttack();

    if( debugInfo == 1 ){ std::cout << "                    Selecting Source of Dynamic Pressure" << std::endl; }
    double currentDynamicPressure = 0.0;

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "                    Selecting Initial Values" << std::endl; }

        currentDynamicPressure = bislipSystems->getInitialDynamicPressure();
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "                    Selecting Propagated Values" << std::endl; }

        currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    }

    //const double airDensity    = flightConditions->getCurrentDensity();
    //const double airSpeed      = flightConditions->getCurrentAirspeed();
    //const double machNumber    = flightConditions->getCurrentMachNumber();
    if( debugInfo == 1 ){ std::cout << "                Current Angle Of Attack  = " << tudat::unit_conversions::convertRadiansToDegrees( currentAngleOfAttack ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Dynamic Pressure = " << currentDynamicPressure << std::endl; }

    //! Extract information from Bislip/Vehicle Systems.
    const double referenceArea = bislipSystems->getReferenceValues( )( 0 );
    if( debugInfo == 1 ){ std::cout << "                Reference Area           = " << referenceArea << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Calculate Aerodynamic Load in the Aerodynamic Frame" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                ---------------------------------------------------" << std::endl; }
    Eigen::Vector3d aerodynamicLoad = -currentDynamicPressure * referenceArea * ( bislip::Variables::computeFullCurrentCoefficients( bodyMap, vehicleName ) ).segment( 0, 3 ) ;
    if( debugInfo == 1 ){ std::cout << "                    Aerodynamic Load in the Aerodynamic Frame = [ " << aerodynamicLoad( 0 ) << " , " << aerodynamicLoad( 1 ) << " , " << aerodynamicLoad( 2 ) << " ]" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                ---------------------------------------------------" << std::endl; }

    return aerodynamicLoad ;

}


Eigen::Vector3d computeAerodynamicFrameTotalLoad(
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

    if( debugInfo == 1 ){ std::cout << "            Starting Computation of Aerodynamic Frame Total Load" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            ---------------------------------------------------" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Extracting current data" << std::endl; }
    const double currentAngleOfAttack   = bislipSystems->getCurrentAngleOfAttack();

    Eigen::Vector3d aerodynamicFrameTotalLoad =  tudat::reference_frames::getBodyToAirspeedBasedAerodynamicFrameTransformationMatrix( currentAngleOfAttack, 0.0 ) * ( bislip::Variables::computeBodyFixedTotalLoad( bodyMap, vehicleName ) );


    return aerodynamicFrameTotalLoad;

}


Eigen::Vector3d computeAerodynamicFrameTotalAcceleration(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    //! Extract Vehicle Systems pointer.
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems( ) ;

    if( debugInfo == 1 ){ std::cout << "            Starting Computation of Aerodynamic Frame Total Acceleration" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "            ---------------------------------------------------" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "                Extracting relevant data" << std::endl; }
    const double currentMass = bodyMap.at( vehicleName )->getBodyMass();

    return bislip::Variables::computeAerodynamicFrameTotalLoad( bodyMap, vehicleName ) / currentMass;
}

Eigen::Vector3d computePassengerFrameJerk(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName,
        const std::function< double() > massRateFunction)//const double massRate )
//        const std::map< std::string, std::shared_ptr< tudat::basic_astrodynamics::MassRateModel > > &nBodyModel )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;


    int debugInfo = bislipSystems->getDebugInfo();


    //double currentMach                  = flightConditions->getCurrentMachNumber();
    double currentMass                  = bodyMap.at( vehicleName )->getBodyMass();
    double sq_currentMass               = currentMass * currentMass;
    //std::vector< std::string > bodiesToIntegrate;
    //bodiesToIntegrate.push_back( vehicleName );
    //tudat::propagators::BodyMassStateDerivative< double, double > massModel( nBodyModel, bodiesToIntegrate );

    //double massRate                     = massModel.getTotalMassRateForBody( vehicleName );
    double massRate = massRateFunction();
    double currentAngleOfAttack         = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );
    double currentThrustElevationAngle  = bislipSystems->getCurrentThrustElevationAngle();
    double currentThrustAzimuthAngle    = bislipSystems->getCurrentThrustAzimuthAngle();
    double currentThrustMagnitude       = bislipSystems->getCurrentThrustMagnitude();
    double currentDragCoefficient       = ( bislipSystems->getFullCurrentCoefficients() )( 0 );
    double currentLiftCoefficient       = ( bislipSystems->getFullCurrentCoefficients() )( 2 );
    Eigen::Vector3d currentAerodynamicFrameAerodynamicLoad = bislip::Variables::computeAerodynamicFrameAerodynamicLoad( bodyMap, vehicleName );
    double currentDragForce             = std::abs( currentAerodynamicFrameAerodynamicLoad( 0 ) );
    double currentLiftForce             = std::abs( currentAerodynamicFrameAerodynamicLoad( 2 ) );

    const double s_alpha                = std::sin( currentAngleOfAttack );
    const double c_alpha                = std::cos( currentAngleOfAttack );

    const double s_thrustElevationAngle = std::sin( currentThrustElevationAngle );
    const double c_thrustElevationAngle = std::cos( currentThrustElevationAngle );
    const double s_thrustAzimuthAngle   = std::sin( currentThrustAzimuthAngle );
    const double c_thrustAzimuthAngle   = std::cos( currentThrustAzimuthAngle );

    double currentHeight                = flightConditions->getCurrentAltitude();
    double currentAirspeed              = flightConditions->getCurrentAirspeed();
    double currentDensity               = flightConditions->getCurrentDensity();
    double currentFlightPathAngle       = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );

    double propagationStepSize          = bislipSystems->getPropagationStepSize();
    double referenceArea                = bislipSystems->getReferenceArea();


    //double heightRate                   = currentAirspeed * std::sin( currentFlightPathAngle );
    //double machRate                     = ( currentMach - previousMach ) / propagationStepSize;
    //double angleOfAttackRate            = ( currentAngleOfAttack - previousAngleOfAttack ) / propagationStepSize;
    //double densityRate                  = bislip::Variables::computeDensityRate( currentHeight, heightRate, bislipSystems->getDensityParameterMapForJerk() );
    //double airspeedRate                 = bislip::Variables::computeAirspeedRate( bodyMap, vehicleName );
    //double thrustElevationAngleRate     = ( currentThrustElevationAngle - previousThrustElevationAngle ) / propagationStepSize;
    //double thrustAzimuthAngleRate       = ( currentThrustAzimuthAngle - previousThrustAzimuthAngle ) / propagationStepSize;
    //double thrustRate                   = ( currentThrust - previousThrust ) / propagationStepSize;
    //double dragCoefficientRate          = ( currentDragCoefficient - previousDragCoefficient ) / propagationStepSize;
    //double liftCoefficientRate          = ( currentLiftCoefficient - previousLiftCoefficient ) / propagationStepSize;

    //double previousMach                 = previousConditions[ bislip::BislipVehicleSystems::PreviousConditionRates::mach_number ];
    std::map< bislip::BislipVehicleSystems::PreviousConditions, double > previousConditions = bislipSystems->getPreviousConditions( );


    double previousHeight               = previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::height ];
    double previousDensity              = previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::density ];
    double previousAirspeed             = previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::airspeed ];
    double previousAngleOfAttack        = previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::angle_of_attack ];
    double previousThrustElevationAngle = previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::thrust_elevation_angle ];
    double previousThrustAzimuthAngle   = previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::thrust_azimuth_angle ];
    double previousThrustMagnitude      = previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::thrust_magnitude ];
    double previousDragCoefficient      = previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::drag_coefficient ];
    double previousLiftCoefficient      = previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::lift_coefficient ];


    double heightRate               = ( currentHeight - previousHeight ) / propagationStepSize;
    double densityRate              = ( currentDensity - previousDensity ) / propagationStepSize;
    double airspeedRate             = ( currentAirspeed - previousAirspeed ) / propagationStepSize;
    double angleOfAttackRate        = ( currentAngleOfAttack - previousAngleOfAttack ) / propagationStepSize;
    double thrustElevationAngleRate = ( currentThrustElevationAngle - previousThrustElevationAngle ) / propagationStepSize;
    double thrustAzimuthAngleRate   = ( currentThrustAzimuthAngle - previousThrustAzimuthAngle ) / propagationStepSize;
    double thrustRate               = ( currentThrustMagnitude - previousThrustMagnitude ) / propagationStepSize;
    double dragCoefficientRate      = ( currentDragCoefficient - previousDragCoefficient ) / propagationStepSize;
    double liftCoefficientRate      = ( currentLiftCoefficient - previousLiftCoefficient ) / propagationStepSize;


    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::height ]                 = currentHeight;
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::density ]                = currentDensity;
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::airspeed ]               = currentAirspeed;
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::angle_of_attack ]        = currentAngleOfAttack;
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::thrust_elevation_angle ] = currentThrustElevationAngle;
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::thrust_azimuth_angle ]   = currentThrustAzimuthAngle;
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::thrust_magnitude ]       = currentThrustMagnitude;
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::drag_coefficient ]       = currentDragCoefficient;
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::lift_coefficient ]       = currentLiftCoefficient;




    bislipSystems->setPreviousConditions( previousConditions );




//    double dragRate = ( 1.0 / 2.0 ) * referenceArea * currentAirspeed * ( currentDensity * ( currentAirspeed * dragCoefficientRate + 2 * currentDragCoefficient * airspeedRate ) + currentDragCoefficient * currentAirspeed * densityRate );
    double dragRate = ( 1.0 / 2.0 ) * referenceArea * currentAirspeed * ( currentAirspeed * ( currentDensity * dragCoefficientRate + currentDragCoefficient * densityRate ) + 2 * currentDragCoefficient * currentDensity * airspeedRate );

   // double liftRate = ( 1.0 / 2.0 ) * referenceArea * currentAirspeed * ( currentDensity * ( currentAirspeed * liftCoefficientRate + 2 * currentLiftCoefficient * airspeedRate ) + currentLiftCoefficient * currentAirspeed * densityRate );
    double liftRate = ( 1.0 / 2.0 ) * referenceArea * currentAirspeed * ( currentAirspeed * ( currentDensity * liftCoefficientRate + currentLiftCoefficient * densityRate ) + 2 * currentLiftCoefficient * currentDensity * airspeedRate );
/*
    std::cout << "Current Mass                = " << currentMass << std::endl;
    std::cout << "Mass Rate                   = " << massRate << std::endl;
    std::cout << "Current Height              = " << currentHeight << std::endl;
    std::cout << "Height Rate                 = " << heightRate << std::endl;
    std::cout << "Current Density             = " << currentDensity << std::endl;
    std::cout << "Density Rate                = " << densityRate << std::endl;
    std::cout << "Current Airspeed            = " << currentAirspeed << std::endl;
    std::cout << "Airspeed Rate               = " << airspeedRate << std::endl;
   // std::cout << "Propagation Step Size       = " << propagationStepSize << std::endl;

    std::cout << "Current Angle Of Attack     = " << currentAngleOfAttack << std::endl;
   // std::cout << "Previous Angle Of Attack    = " << previousAngleOfAttack << std::endl;
    std::cout << "Angle Of Attack Rate        = " << angleOfAttackRate << std::endl;


    std::cout << "Thrust Elevation Angle Rate = " << thrustElevationAngleRate << std::endl;
    std::cout << "Thrust Azimuth Angle Rate   = " << thrustAzimuthAngleRate << std::endl;
    std::cout << "Current Thrust              = " << currentThrustMagnitude << std::endl;
    //std::cout << "Previous Thrust             = " << previousThrust << std::endl;
    std::cout << "Thrust Rate                 = " << thrustRate << std::endl;
    std::cout << "Current Drag Force          = " << currentDragForce << std::endl;
    std::cout << "Current Lift Force          = " << currentLiftForce << std::endl;
    std::cout << "Drag Rate                   = " << dragRate << std::endl;
    std::cout << "Lift Rate                   = " << liftRate << std::endl;
    std::cout << "Current Drag Coefficient    = " << currentDragCoefficient << std::endl;
    std::cout << "Current Lift Coefficient    = " << currentLiftCoefficient << std::endl;
    //std::cout << "Previous Drag Coefficient   = " << previousDragCoefficient << std::endl;
    //std::cout << "Previous Lift Coefficient   = " << previousLiftCoefficient << std::endl;
    std::cout << "Drag Coefficient Rate       = " << dragCoefficientRate << std::endl;
    std::cout << "Lift Coefficient Rate       = " << liftCoefficientRate << std::endl;
*/

    Eigen::Vector3d passengerFrameJerk;
    double passengerFrameJerk_x11 = referenceArea * currentAirspeed * currentAirspeed * currentDragCoefficient * densityRate * c_alpha / ( 2 * currentMass );
    double passengerFrameJerk_x12 = referenceArea * currentDensity * currentAirspeed * currentDragCoefficient * airspeedRate * c_alpha / currentMass;
    double passengerFrameJerk_x13 = referenceArea * currentDensity * currentAirspeed * currentAirspeed * dragCoefficientRate * c_alpha /  ( 2 * currentMass );
    double passengerFrameJerk_x14 = -referenceArea * currentDensity * currentAirspeed * currentAirspeed * currentDragCoefficient * s_alpha * angleOfAttackRate /  ( 2 * currentMass );
    double passengerFrameJerk_x15 = -referenceArea * currentDensity * currentAirspeed * currentAirspeed * currentDragCoefficient * c_alpha * massRate /  ( 2 * currentMass * currentMass );

    double passengerFrameJerk_x21 = referenceArea * currentAirspeed * currentAirspeed * currentLiftCoefficient * densityRate * s_alpha / ( 2 * currentMass );
    double passengerFrameJerk_x22 = referenceArea * currentDensity * currentAirspeed * currentLiftCoefficient * airspeedRate * s_alpha / currentMass;
    double passengerFrameJerk_x23 = referenceArea * currentDensity * currentAirspeed * currentAirspeed * liftCoefficientRate * s_alpha /  ( 2 * currentMass );
    double passengerFrameJerk_x24 = referenceArea * currentDensity * currentAirspeed * currentAirspeed * currentLiftCoefficient * c_alpha * angleOfAttackRate /  ( 2 * currentMass );
    double passengerFrameJerk_x25 = -referenceArea * currentDensity * currentAirspeed * currentAirspeed * currentLiftCoefficient * s_alpha * massRate /  ( 2 * currentMass * currentMass );

    double passengerFrameJerk_x1 = passengerFrameJerk_x11 + passengerFrameJerk_x12 + passengerFrameJerk_x13 + passengerFrameJerk_x14 + passengerFrameJerk_x15;
    double passengerFrameJerk_x2 = passengerFrameJerk_x21 + passengerFrameJerk_x22 + passengerFrameJerk_x23 + passengerFrameJerk_x24 + passengerFrameJerk_x25;
    double passengerFrameJerk_x3 = ( currentThrustMagnitude * massRate * c_thrustElevationAngle * c_thrustAzimuthAngle / sq_currentMass
                                     - thrustRate * c_thrustElevationAngle * c_thrustAzimuthAngle / currentMass
                                     + currentThrustMagnitude * thrustElevationAngleRate * s_thrustElevationAngle * c_thrustAzimuthAngle / currentMass
                                     + currentThrustMagnitude * thrustAzimuthAngleRate * c_thrustElevationAngle * s_thrustAzimuthAngle / currentMass );


    double passengerFrameJerk_z11 = referenceArea * currentAirspeed * currentAirspeed * currentDragCoefficient * densityRate * s_alpha / ( 2 * currentMass );
    double passengerFrameJerk_z12 = referenceArea * currentDensity * currentAirspeed * currentDragCoefficient * airspeedRate * s_alpha / currentMass;
    double passengerFrameJerk_z13 = referenceArea * currentDensity * currentAirspeed * currentAirspeed * dragCoefficientRate * s_alpha /  ( 2 * currentMass );
    double passengerFrameJerk_z14 = referenceArea * currentDensity * currentAirspeed * currentAirspeed * currentDragCoefficient * c_alpha * angleOfAttackRate /  ( 2 * currentMass );
    double passengerFrameJerk_z15 = -referenceArea * currentDensity * currentAirspeed * currentAirspeed * currentDragCoefficient * s_alpha * massRate /  ( 2 * currentMass * currentMass );

    double passengerFrameJerk_z21 = referenceArea * currentAirspeed * currentAirspeed * currentLiftCoefficient * densityRate * c_alpha / ( 2 * currentMass );
    double passengerFrameJerk_z22 = referenceArea * currentDensity * currentAirspeed * currentLiftCoefficient * airspeedRate * c_alpha / currentMass;
    double passengerFrameJerk_z23 = referenceArea * currentDensity * currentAirspeed * currentAirspeed * liftCoefficientRate * c_alpha /  ( 2 * currentMass );
    double passengerFrameJerk_z24 = -referenceArea * currentDensity * currentAirspeed * currentAirspeed * currentLiftCoefficient * s_alpha * angleOfAttackRate /  ( 2 * currentMass );
    double passengerFrameJerk_z25 = -referenceArea * currentDensity * currentAirspeed * currentAirspeed * currentLiftCoefficient * c_alpha * massRate /  ( 2 * currentMass * currentMass );

    double passengerFrameJerk_z1 = -( passengerFrameJerk_z11 + passengerFrameJerk_z12 + passengerFrameJerk_z13 + passengerFrameJerk_z14 + passengerFrameJerk_z15 );
    double passengerFrameJerk_z2 = passengerFrameJerk_z21 + passengerFrameJerk_z22 + passengerFrameJerk_z23 + passengerFrameJerk_z24 + passengerFrameJerk_z25;
    double passengerFrameJerk_z3 = ( currentMass * currentThrustMagnitude * thrustElevationAngleRate * c_thrustElevationAngle - currentThrustMagnitude * s_thrustElevationAngle * massRate + currentMass * s_thrustElevationAngle * thrustRate ) / ( currentMass * currentMass );

    //double passengerFrameJerk_x1 = currentMass * c_alpha * dragRate - currentDragForce * ( currentMass * angleOfAttackRate * s_alpha + massRate * c_alpha ) / sq_currentMass;
    //double passengerFrameJerk_x2 = currentMass * s_alpha * liftRate + currentLiftForce * ( currentMass * angleOfAttackRate * c_alpha - massRate * s_alpha ) / sq_currentMass;
    //double passengerFrameJerk_x3 = ( currentThrustMagnitude * massRate * c_thrustElevationAngle * c_thrustAzimuthAngle / sq_currentMass
    //                                - thrustRate * c_thrustElevationAngle * c_thrustAzimuthAngle / currentMass
    //                              + currentThrustMagnitude * thrustElevationAngleRate * s_thrustElevationAngle * c_thrustAzimuthAngle / currentMass
         //                            + currentThrustMagnitude * thrustAzimuthAngleRate * c_thrustElevationAngle * s_thrustAzimuthAngle / currentMass );

//double passengerFrameJerk_x31 = ( currentMass * thrustRate - currentThrustMagnitude * massRate ) / sq_currentMass;
    //std::cout << "passengerFrameJerk_x1       = " << passengerFrameJerk_x1 << std::endl;
    //std::cout << "passengerFrameJerk_x2       = " << passengerFrameJerk_x2 << std::endl;
    //std::cout << "passengerFrameJerk_x3       = " << passengerFrameJerk_x3 << std::endl;
    //std::cout << "passengerFrameJerk_x31      = " << passengerFrameJerk_x31 << std::endl;




    passengerFrameJerk( 0 ) = passengerFrameJerk_x1 + passengerFrameJerk_x2 + passengerFrameJerk_x3;


    passengerFrameJerk( 1 ) =
            - currentThrustMagnitude * massRate * c_thrustElevationAngle * s_thrustAzimuthAngle / sq_currentMass
            + thrustRate * c_thrustElevationAngle * s_thrustAzimuthAngle / currentMass
            - currentThrustMagnitude * thrustElevationAngleRate * s_thrustElevationAngle * c_thrustAzimuthAngle / currentMass
            + currentThrustMagnitude * thrustAzimuthAngleRate * c_thrustElevationAngle * c_thrustAzimuthAngle / currentMass;


    passengerFrameJerk( 2 ) = passengerFrameJerk_z1 + passengerFrameJerk_z2 + passengerFrameJerk_z3;


    bislipSystems->setPreviousConditions( previousConditions );

    return passengerFrameJerk;
}


double computeDensityRate(
        const double &height,
        const double &heightRate,
        const std::map < int, Eigen::VectorXd > &densityParameterMap )
{
    //! https://en.wikipedia.org/wiki/Barometric_formula#Density_equations
    double densityRate = 0.0;

    //for( std::map< int, Eigen::VectorXd >::iterator it = densityParameterMap.begin(); it != densityParameterMap.end(); ++it )
    for( int i = 0; i < 7; i++ )
    {

        if( i != 6)
        {
            if( height >= ( densityParameterMap.at( i ) )( 0 ) && height < ( densityParameterMap.at( i + 1 ) )( 0 ) )
            {
                const double h   = ( densityParameterMap.at( i ) )( 0 );
                const double rho = ( densityParameterMap.at( i ) )( 1 );
                const double T   = ( densityParameterMap.at( i ) )( 2 );
                const double L   = ( densityParameterMap.at( i ) )( 3 );
                const double R   = ( densityParameterMap.at( i ) )( 4 );
                const double g_0 = ( densityParameterMap.at( i ) )( 5 );
                const double M   = ( densityParameterMap.at( i ) )( 6 );
                const double d1  = ( g_0 * M ) / ( R * L );
                const double d2  = ( g_0 * M ) / ( R * T );

                if( L != 0.0 )
                {
                    const double den = T + L * ( height - h );

                    densityRate = -T * rho * L * ( d1 + 1 ) * heightRate * std::pow( T / den , d1 ) / ( den * den );
                }
                else
                {
                    densityRate = -d2 * rho * heightRate * std::exp( d2 * ( h - height ) );
                }
            }
        }
        else if( i == 6 )
        {
            if( height >= ( densityParameterMap.at( i ) )( 0 ) )
            {
                const double h   = ( densityParameterMap.at( i ) )( 0 );
                const double rho = ( densityParameterMap.at( i ) )( 1 );
                const double T   = ( densityParameterMap.at( i ) )( 2 );
                const double L   = ( densityParameterMap.at( i ) )( 3 );
                const double R   = ( densityParameterMap.at( i ) )( 4 );
                const double g_0 = ( densityParameterMap.at( i ) )( 5 );
                const double M   = ( densityParameterMap.at( i ) )( 6 );
                const double d1  = ( g_0 * M ) / ( R * L );
                const double den = T + L * ( height - h );

                densityRate = -T * rho * L * ( d1 + 1 ) * heightRate * std::pow( T / den , d1 ) / ( den * den );

            }
        }
    }



    return densityRate;
}


double computeAirspeedRate(
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
    if( debugInfo == 1 ){ std::cout << "        Starting Determination of Airspeed Rate" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "        ------------------------------------------------" << std::endl; }

    double currentFlightPathAngle;
    double currentLatitude;
    double currentHeading;
    double currentAltitude;
    double currentAirspeed;
    double currentDynamicPressure;

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentFlightPathAngle = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialFlightPathAngle() );
        currentLatitude        = bislipSystems->getInitialLat();
        currentHeading         = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialHeadingAngle() );
        currentAltitude        = bislipSystems->getInitialAltitude();
        currentAirspeed        = bislipSystems->getInitialAirspeed();
        currentDynamicPressure = bislipSystems->getInitialDynamicPressure();
        //currentLift            = bislipSystems->getCurrentLiftForce();
    }

    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentFlightPathAngle = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        currentLatitude        = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        currentHeading         = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
        currentAltitude        = flightConditions->getCurrentBodyCenteredBodyFixedState( ).segment( 0, 3 ).norm( );
        currentAirspeed        = flightConditions->getCurrentAirspeed( );
        currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
        //currentLift            = bislip::Variables::computeCurrentLiftForce( bodyMap, vehicleName );

    }

    //const double currentDynamicPressure = flightConditions->getCurrentDynamicPressure();
    const double currentMass       = bodyMap.at( vehicleName )->getBodyMass();

    const double rotationRateEarth = 7.292115*1E-5;

    //const Eigen::Vector6d currentCoefficients = bislip::Variables::computePartialCurrentCoefficients( bodyMap, vehicleName );
    //const double currentLift = bislip::Variables::computeCurrentLiftForce( bodyMap, vehicleName );
    //const double currentLift = bislipSystems->getCurrentLiftForce();
    const double currentDrag = std::abs( ( bislip::Variables::computeAerodynamicFrameAerodynamicLoad( bodyMap, vehicleName ) )( 0 ) );

    //const Eigen::Vector2d gravs = bislip::Variables::computeLocalGravity( bodyMap, vehicleName, centralBodyName );
    const Eigen::Vector3d gravs = bislipSystems->getCurrentLocalGravityVector();

    const double currentAngleOfAttack        = bislipSystems->getCurrentAngleOfAttack();
    const double currentBankAngle            = bislipSystems->getCurrentBankAngle();
    const double currentThrustMagnitude      = bislipSystems->getCurrentThrustMagnitude();
    const double currentThrustElevationAngle = bislipSystems->getCurrentThrustElevationAngle();
    const double currentThrustAzimuthAngle   = bislipSystems->getCurrentThrustAzimuthAngle();

    if( debugInfo == 1 ){ std::cout << "                Current Flight-Path Angle      = " << currentFlightPathAngle << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Latitude               = " << tudat::unit_conversions::convertRadiansToDegrees( currentLatitude ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Heading                = " << tudat::unit_conversions::convertRadiansToDegrees( currentHeading ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Altitude               = " << currentAltitude << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Airspeed               = " << currentAirspeed << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Dynamic Pressure       = " << currentDynamicPressure << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Angle of Attack        = " << tudat::unit_conversions::convertRadiansToDegrees( currentAngleOfAttack ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Bank Angle             = " << tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Thrust Magnitude       = " << currentThrustMagnitude << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Thrust Elevation Angle = " << tudat::unit_conversions::convertRadiansToDegrees( currentThrustElevationAngle ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Thrust Azimuth Angle   = " << tudat::unit_conversions::convertRadiansToDegrees( currentThrustAzimuthAngle ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Drag Force             = " << currentDrag << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                Current Mass                   = " << currentMass << std::endl; }



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

    const double F_V = -currentDrag + currentThrustMagnitude * c_alpha * c_thrustAzimuthAngle * c_thrustElevationAngle
            - currentThrustMagnitude * s_alpha * s_thrustElevationAngle
            - currentMass * ( gravs( 0 ) * c_gamma * c_chi - gravs( 2 ) * s_gamma );

    const double A = rotationRateEarth * rotationRateEarth * currentAltitude * c_delta *  ( c_delta * s_gamma - c_gamma * s_delta * c_chi );

    const double airspeedRate = ( ( F_V / currentMass ) + A );
    if( debugInfo == 1 ){ std::cout << "                Airspeed Rate         = " << airspeedRate << std::endl; }


    return airspeedRate;
}



Eigen::Matrix3d getLocalVerticalToBodyFrameTransformationMatrix(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    //! Extract current conditions.
    if( debugInfo == 1 ){ std::cout << "        Selecting source of values" << std::endl; }

    double currentAngleOfAttack;
    double currentAngleOfSideSlip;
    double currentBankAngle;
    double currentFlightPathAngle;
    double currentHeadingAngle;

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting initial values" << std::endl; }

        currentAngleOfAttack   = bislipSystems->getCurrentAngleOfAttack();
        currentAngleOfSideSlip = 0.0;
        currentBankAngle       = bislipSystems->getCurrentBankAngle();
        currentFlightPathAngle = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialFlightPathAngle() );
        currentHeadingAngle    = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialHeadingAngle() );
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting propagated values" << std::endl; }

        currentAngleOfAttack   = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );
        currentAngleOfSideSlip = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_sideslip );
        currentBankAngle       = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::bank_angle );
        currentFlightPathAngle = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        currentHeadingAngle    = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
    }


    Eigen::Matrix3d localVerticalToBodyFrameTransformationMatrix =
            tudat::reference_frames::getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix( currentAngleOfAttack, currentAngleOfSideSlip ) *
            tudat::reference_frames::getTrajectoryToAerodynamicFrameTransformationMatrix( currentBankAngle ) *
            tudat::reference_frames::getLocalVerticalFrameToTrajectoryTransformationMatrix( currentFlightPathAngle, currentHeadingAngle );


    return localVerticalToBodyFrameTransformationMatrix;
}

Eigen::VectorXd computeNumericalDerivativeOfVector(
        const Eigen::VectorXd &f,
        const Eigen::VectorXd &x )
{
    Eigen::VectorXd numericalDerivative( f.size() );

    //! Left boundary value
    double x2minusx1 = ( x( 1 ) - x( 0 ) );
    double x3minusx1 = ( x( 2 ) - x( 0 ) );
    double x3minusx2 = ( x( 2 ) - x( 1 ) );
    double leftBoundaryNUM = -x2minusx1 * x2minusx1 * f( 2 ) + x3minusx1 * x3minusx1 * f( 1 ) - ( x3minusx1 * x3minusx1 - x2minusx1 * x2minusx1 ) * f( 0 );
    double leftBoundaryDEN =  x2minusx1 * x3minusx1 * x3minusx2;
    numericalDerivative( 0 ) = leftBoundaryNUM / leftBoundaryDEN;

    //! First point AFTER boundary
    numericalDerivative( 1 ) = ( f( 2 ) - f( 0 ) ) / ( 2 * ( x( 2 ) - x( 0 ) ) );

    //! Last point BEFORE boundary
    numericalDerivative( numericalDerivative.size( ) - 2 ) = ( f( f.size() - 1 ) - f( f.size() - 3 ) ) / ( 2 * ( x( x.size() - 1 ) - x( x.size() - 3 ) ) );

    //! Right boundary value
    double xn1minusxn = ( x( 1 ) - x( 0 ) );
    double x3n2minusxn = ( x( 2 ) - x( 0 ) );
    double x3n2minusxn1 = ( x( 2 ) - x( 1 ) );
    double rightBoundaryNUM = x2minusx1 * x2minusx1 * f( f.size() - 3 ) - x3minusx1 * x3minusx1 * f( f.size() - 2 ) + ( x3minusx1 * x3minusx1 - x2minusx1 * x2minusx1 ) * f( f.size() - 1 );
    double rightBoundaryDEN =  x3n2minusxn * x3n2minusxn * x3n2minusxn1;
    numericalDerivative( numericalDerivative.size() - 1 ) = rightBoundaryNUM / rightBoundaryDEN;
    //numericalDerivative( numericalDerivative.size() - 1 ) = ( f( f.size() - 1 ) - f( f.size() - 2 ) ) / ( x( x.size() - 1 ) - x( x.size() - 2 ) );


    //! Central difference for the interior
    for ( unsigned int i = 2; i < numericalDerivative.size( ) - 2; i++ )
    {
        //numericalDerivative( i ) = ( f( i + 1 ) - f( i - 1 ) ) / ( 2 * ( x( i + 1 ) - x( i - 1 ) ) );
        numericalDerivative( i ) = ( -f( i + 2 ) + 8 * f( i + 1 ) - 8 * f( i - 1 ) + f( i - 2 ) ) / ( 12 * ( x( i + 1 ) - x( i - 1 ) ) );

    }

    return numericalDerivative;
}


Eigen::Vector3d computeBodyFixedTotal_g_Load_Vector(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "                Determine Total Body-Fixed Load Vector" << std::endl; }

    Eigen::Vector3d bodyFixedTotal_g_Load_Vector = bislip::Variables::computeBodyFixedTotalLoad( bodyMap, vehicleName );

    if( debugInfo == 1 ){ std::cout << "                Normalizing the Total Body-Fixed Load Vector" << std::endl; }

    //! Determine the Body-Fixed Total Load g-load vector.
    Eigen::Vector3d totalLoadBodyFrame_g_Load_Vector = bodyFixedTotal_g_Load_Vector / ( tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * bodyMap.at( vehicleName )->getBodyMass() );

    return totalLoadBodyFrame_g_Load_Vector;
}

double computeBodyFixedTotal_g_Load_Magnitude (
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "                Determine Total Body-Fixed g-load Vector" << std::endl; }
    Eigen::Vector3d bodyFixedTotal_g_Load_Vector = bislip::Variables::computeBodyFixedTotal_g_Load_Vector( bodyMap, vehicleName );

    if( debugInfo == 1 ){ std::cout << "                Calculating Magnitude of the Total Body-Fixed g-load Vector" << std::endl; }
    double totalLoadBodyFrame_g_Load_Magnitude = ( bodyFixedTotal_g_Load_Vector ).norm( );

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

Eigen::Vector3d computePassengerFrameTotalLoad(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Determine the Passenger-Fixed Total Load g-load vector.
    Eigen::Vector3d passengerFrameTotaLoad = ( bislipSystems->getBodyFrameToPassengerFrameTransformationMatrix() ) * ( bislip::Variables::computeBodyFixedTotalLoad( bodyMap, vehicleName ) );

    return passengerFrameTotaLoad;
}

Eigen::Vector3d computePassengerFrameTotal_g_Load_Vector(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Determine the Passenger-Fixed Total Load g-load vector.
    Eigen::Vector3d totalLoadPassengerFrame_g_Load_Vector = ( bislip::Variables::computePassengerFrameTotalLoad( bodyMap, vehicleName ) ) / ( tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * bodyMap.at( vehicleName )->getBodyMass() );

    return totalLoadPassengerFrame_g_Load_Vector;
}

Eigen::Vector3d computePassengerFrameTotalAcceleration(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Determine the Passenger-Fixed Total Load g-load vector.
    Eigen::Vector3d passengerFrameTotalAcceleration = bislip::Variables::computePassengerFrameTotalLoad( bodyMap, vehicleName ) / bodyMap.at( vehicleName )->getBodyMass();

    return passengerFrameTotalAcceleration;
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
    if( debugInfo == 1 ){ std::cout << "DERP" << std::endl; }
    bislipSystems->setTempBankAngle( newBankAngle );


    /*
    //! Determine if bank angle reversal is required.
    const bool bankReversal = bislip::Variables::determineBankAngleReversal( bodyMap, vehicleName, newBankAngle );

    double returnedBankAngle = newBankAngle;
    //! Impose bank angle reversal on new bank angle.
    if ( bankReversal == true ) { returnedBankAngle = -newBankAngle; }

    if( debugInfo == 1 ){ std::cout << "newBankAngle = " << newBankAngle << "   |   returnedBankAngle = " << returnedBankAngle << "   |   bislipSystems->getBankAngleReversalTrigger( ) = " << bislipSystems->getBankAngleReversalTrigger(  ) << std::endl; }
*/

    double returnedBankAngle = bislip::Variables::returnReversedBankAngle( bodyMap, vehicleName, newBankAngle );

    //if( debugInfo == 1 ){ std::cout << "newBankAngle = " << newBankAngle << "   |   returnedBankAngle = " << returnedBankAngle << "   |   bislipSystems->getBankAngleReversalTrigger( ) = " << bislipSystems->getBankAngleReversalTrigger(  ) << std::endl; }

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

    if( debugInfo == 1 ){ std::cout << "newBankAngle = " << newBankAngle << "   |   returnedBankAngle = " << returnedReversedBankAngle << "   |   bislipSystems->getBankAngleReversalTrigger( ) = " << bislipSystems->getBankAngleReversalTrigger(  ) << std::endl; }

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
    newBankAngle = bislip::Variables::determineSignedBankAngle( bislip::Variables::determineSignOfValue( bislipSystems->getCurrentBankAngle() ), std::abs( bislipSystems->getEvaluatedBankAngle( ) ) );
    //}

    return newBankAngle;
}

int determineSignOfValue ( const double &value )
{
    //! This ensures that the new bank angle has the same sign as in current bank angle.
    int sign = 1;
    if ( value < 0.0 ){ sign = -1; }
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

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "    Starting Computation of Heading Error Deadband" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    ---------------------------------------------------" << std::endl; }

    //! Extract current conditions.
    if( debugInfo == 1 ){ std::cout << "        Selecting source of airspeed/altitude" << std::endl; }

    double currentLatitude;
    double currentLongitude;

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting initial values" << std::endl; }

        currentLatitude  = bislipSystems->getInitialLat();
        currentLongitude = bislipSystems->getInitialLon();
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting propagated values" << std::endl; }

        currentLatitude  = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        currentLongitude = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );
    }


    if( debugInfo == 1 ){ std::cout << "        Calculate Angular Distance To Go" << std::endl; }
    const double angularDistanceToGo_deg = tudat::unit_conversions::convertRadiansToDegrees(
                bislip::Variables::computeAngularDistance (
                    currentLatitude,
                    currentLongitude,
                    bodyMap.at( vehicleName )->getBislipSystems()->getTargetLat(),
                    bodyMap.at( vehicleName )->getBislipSystems()->getTargetLon() ) );

    double headingErrorDeadBand = 0.0;
    if( angularDistanceToGo_deg > bislipSystems->getHeadingErrorDeadBandLowDistanceTrigger() )
    { headingErrorDeadBand = bislipSystems->getHeadingErrorDeadBandCoarseInterpolator()->interpolate( angularDistanceToGo_deg ); }
    else
    { headingErrorDeadBand = bislipSystems->getHeadingErrorDeadBandLowDistanceInterpolator()->interpolate( angularDistanceToGo_deg ); }

    return headingErrorDeadBand;
}

double convertRadiansToDegrees( const double &angleInRadians )
{ return angleInRadians * 180.0 / tudat::mathematical_constants::PI; }

double convertDegreesToRadians( const double &angleInDegrees )
{ return angleInDegrees * tudat::mathematical_constants::PI / 180.0; }

double convertNegativeAnglesInRadiansToPositive( const double &angleInRadians )
{
    double positiveAngle = angleInRadians;
    if( angleInRadians < 0.0 ) { positiveAngle = angleInRadians + 2 * tudat::mathematical_constants::PI; }
    return positiveAngle;
}

double convertNegativeAnglesInDegreesToPositive( const double &angleInDegrees )
{
    double positiveAngle = angleInDegrees;
    if( angleInDegrees < 0.0 ) { positiveAngle = angleInDegrees + 360.0; }
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
        const std::vector< double > &headingErrorDeadBandCoarse,
        const std::vector< double > &headingErrorDeadBandLowDistance,
        const std::string &outputPath,
        const std::string &outputSubFolder )
{
    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Coarse interpolator section
    std::map< double, double > map_headingErrorDeadBandCoarseInterpolator;

    //! Declare and initialize Interpolator Settings.
    std::shared_ptr< tudat::interpolators::InterpolatorSettings > headingErrorDeadBandCoarseInterpolatorSettings = std::make_shared< tudat::interpolators::InterpolatorSettings >( tudat::interpolators::InterpolatorTypes::linear_interpolator );

    //! Declare and initialize size vector size.
    unsigned long deadBandSizeCoarse = headingErrorDeadBandCoarse.size() / 2;

    //! Declare Heading Error Deadband vectors.
    Eigen::VectorXd headingErrorDeadBandCoarse_distance( deadBandSizeCoarse );
    Eigen::VectorXd headingErrorDeadBandCoarse_error(deadBandSizeCoarse );

    //! Populate vectors that will be used with data map.
    for( unsigned long i = 0; i < deadBandSizeCoarse; i++ )
    {
        headingErrorDeadBandCoarse_distance( i ) =  headingErrorDeadBandCoarse[ i ];
        headingErrorDeadBandCoarse_error( i )    =  headingErrorDeadBandCoarse[ i + deadBandSizeCoarse ];
    }

    //! Populate data map for interpolator.
    for ( unsigned long i = 0; i < deadBandSizeCoarse; ++i )
    { map_headingErrorDeadBandCoarseInterpolator[ headingErrorDeadBandCoarse_distance( i ) ] = headingErrorDeadBandCoarse_error( i ); }

    // std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > headingErrorDeadBandCoarseInterpolator =
    //       tudat::interpolators::createOneDimensionalInterpolator( map_headingErrorDeadBandCoarseInterpolator, interpolatorSettings );
    std::pair< double, double > headingErrorDeadBandCoarseInterpolatorDomainInterval = std::make_pair( headingErrorDeadBandCoarse_distance.minCoeff(), headingErrorDeadBandCoarse_distance.maxCoeff() );
    std::pair< double, double > headingErrorDeadBandCoarseInterpolatorRangeInterval  = std::make_pair( headingErrorDeadBandCoarse_error.minCoeff(), headingErrorDeadBandCoarse_error.maxCoeff() );
    std::pair< double, double > headingErrorDeadBandCoarseInterpolatorBoundaryValues = std::make_pair( headingErrorDeadBandCoarse_error( 0 ), headingErrorDeadBandCoarse_error( headingErrorDeadBandCoarse_error.size() - 1 ) );

    //! Create Heading Error Deadband Interpolator.
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > >  headingErrorDeadBandCoarseInterpolator =
            tudat::interpolators::createOneDimensionalInterpolator< double, double >(
                map_headingErrorDeadBandCoarseInterpolator,
                headingErrorDeadBandCoarseInterpolatorSettings,
                headingErrorDeadBandCoarseInterpolatorBoundaryValues );

    //! LowDistance interpolator section
    std::map< double, double > map_headingErrorDeadBandLowDistanceInterpolator;

    //! Declare and initialize Interpolator Settings.
    std::shared_ptr< tudat::interpolators::InterpolatorSettings > headingErrorDeadBandLowDistanceInterpolatorSettings = std::make_shared< tudat::interpolators::InterpolatorSettings >( tudat::interpolators::InterpolatorTypes::linear_interpolator );

    //! Declare and initialize size vector size.
    unsigned long deadBandSizeLowDistance = headingErrorDeadBandLowDistance.size() / 2;

    //! Declare Heading Error Deadband vectors.
    Eigen::VectorXd headingErrorDeadBandLowDistance_distance( deadBandSizeLowDistance );
    Eigen::VectorXd headingErrorDeadBandLowDistance_error(deadBandSizeLowDistance );

    //! Populate vectors that will be used with data map.
    for( unsigned long i = 0; i < deadBandSizeLowDistance; i++ )
    {
        headingErrorDeadBandLowDistance_distance( i ) =  headingErrorDeadBandLowDistance[ i ];
        headingErrorDeadBandLowDistance_error( i )    =  headingErrorDeadBandLowDistance[ i + deadBandSizeLowDistance ];
    }

    //! Populate data map for interpolator.
    for ( unsigned long i = 0; i < deadBandSizeLowDistance; ++i )
    { map_headingErrorDeadBandLowDistanceInterpolator[ headingErrorDeadBandLowDistance_distance( i ) ] = headingErrorDeadBandLowDistance_error( i ); }

    // std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > headingErrorDeadBandLowDistanceInterpolator =
    //       tudat::interpolators::createOneDimensionalInterpolator( map_headingErrorDeadBandLowDistanceInterpolator, interpolatorSettings );
    std::pair< double, double > headingErrorDeadBandLowDistanceInterpolatorDomainInterval = std::make_pair( headingErrorDeadBandLowDistance_distance.minCoeff(), headingErrorDeadBandLowDistance_distance.maxCoeff() );
    std::pair< double, double > headingErrorDeadBandLowDistanceInterpolatorRangeInterval  = std::make_pair( headingErrorDeadBandLowDistance_error.minCoeff(), headingErrorDeadBandLowDistance_error.maxCoeff() );
    std::pair< double, double > headingErrorDeadBandLowDistanceInterpolatorBoundaryValues = std::make_pair( headingErrorDeadBandLowDistance_error( 0 ), headingErrorDeadBandLowDistance_error( headingErrorDeadBandLowDistance_error.size() - 1 ) );

    //! Set trigger to use the Low Distance Interpolator
    bislipSystems->setHeadingErrorDeadBandLowDistanceTrigger( headingErrorDeadBandLowDistanceInterpolatorDomainInterval.second );

    //! Create Heading Error Deadband Interpolator.
    std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > >  headingErrorDeadBandLowDistanceInterpolator =
            tudat::interpolators::createOneDimensionalInterpolator< double, double >(
                map_headingErrorDeadBandLowDistanceInterpolator,
                headingErrorDeadBandLowDistanceInterpolatorSettings,
                headingErrorDeadBandLowDistanceInterpolatorBoundaryValues );

    //! Set Coarse Heading Error Deadband Interpolator
    bislipSystems->setHeadingErrorDeadBandCoarseInterpolator( headingErrorDeadBandCoarseInterpolator );

    //! Set Low Distance Heading Error Deadband Interpolator
    bislipSystems->setHeadingErrorDeadBandLowDistanceInterpolator( headingErrorDeadBandLowDistanceInterpolator );

    bislip::Variables::printHeadingErrorDeadBandBounds(
                headingErrorDeadBandCoarseInterpolator,
                headingErrorDeadBandCoarseInterpolatorDomainInterval,
                headingErrorDeadBandLowDistanceInterpolator,
                headingErrorDeadBandLowDistanceInterpolatorDomainInterval,
                outputPath,
                outputSubFolder );

}

void printHeadingErrorDeadBandBounds(
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &headingErrorDeadBandCoarseInterpolator,
        const std::pair< double, double > &headingErrorDeadBandCoarseInterpolatorDomainInterval,
        const std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > &headingErrorDeadBandLowDistanceInterpolator,
        const std::pair< double, double > &headingErrorDeadBandLowDistanceInterpolatorDomainInterval,
        const std::string &outputPath,
        const std::string &outputSubFolder )
{
    //! Declare data map to contain vectors of interpolated values.
    std::map< double, Eigen::Vector2d > map_HeadingErrorDeadBandBounds;
    Eigen::Vector2d headingErrorDeadBandBounds;

    //! Loop to populate vectors of interpolated values and then pass to data map.
    //!     Number of evaluations has been arbitrarily selected.
    double pp = 0;
    double domainInterval = headingErrorDeadBandLowDistanceInterpolatorDomainInterval.second - headingErrorDeadBandLowDistanceInterpolatorDomainInterval.first;
    for ( unsigned int i = 0; i < 1001; ++i )
    {
        headingErrorDeadBandBounds( 0 ) = -headingErrorDeadBandLowDistanceInterpolator->interpolate( pp * domainInterval / 1000 );
        headingErrorDeadBandBounds( 1 ) = std::abs( headingErrorDeadBandBounds( 0 ) );
        map_HeadingErrorDeadBandBounds[ pp * domainInterval / 1000 ] = headingErrorDeadBandBounds;
        pp += 1;
    }

    pp = 0;
    domainInterval = headingErrorDeadBandCoarseInterpolatorDomainInterval.second - headingErrorDeadBandLowDistanceInterpolatorDomainInterval.second;
    for ( unsigned int i = 0; i < 1001; ++i )
    {
        headingErrorDeadBandBounds( 0 ) = -headingErrorDeadBandCoarseInterpolator->interpolate( headingErrorDeadBandLowDistanceInterpolatorDomainInterval.second + pp * domainInterval / 1000 );
        headingErrorDeadBandBounds( 1 ) = std::abs( headingErrorDeadBandBounds( 0 ) );
        map_HeadingErrorDeadBandBounds[ headingErrorDeadBandLowDistanceInterpolatorDomainInterval.second + pp * domainInterval / 1000 ] = headingErrorDeadBandBounds;
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
    for ( int i = 0; i < vector.size(); i++ )
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
    if( debugInfo == 1 ){ std::cout << "     currentAirSpeed       = " << flightConditions->getCurrentAirspeed() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentSpeedOfSound   = " << flightConditions->getCurrentSpeedOfSound() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentMachNumber     = " << currentMachNumber << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     currentWallTemp       = " << bislipSystems->getWallTemperature() << std::endl; }

    //! Compute adiabatic wall temperature.
    if( debugInfo == 1 ){ std::cout << "Compute adiabatic wall temperature." << std::endl; }
    double adiabaticWallTemperature =
            bislip::Variables::computeAdiabaticWallTemperature( currentAirTemperature, currentMachNumber );
    if( debugInfo == 1 ){ std::cout << "     adiabaticWallTemperature = " << adiabaticWallTemperature << std::endl; }


    if( debugInfo == 1 ){ std::cout << "Get heat transfer function" << std::endl; }
    std::function< double( const double ) > heatTransferFunction = bislip::Variables::getStagnationHeatTransferFunction( bodyMap, vehicleName );

    if( debugInfo == 1 ){ std::cout << "Find equilibrium wall temperature" << std::endl; }
    double equilibriumWallTemperature = bislip::Variables::findEquilibriumWallTemperature( heatTransferFunction, wallEmissivity, adiabaticWallTemperature, bislipSystems->getWallTemperature() );


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
    if( debugInfo == 1 ){ std::cout << "Set equilibrium wall temperature" << std::endl; }

    bislipSystems->setWallTemperature( equilibriumWallTemperature );

    double heatflux = heatTransferFunction( equilibriumWallTemperature );



    if ( bislipSystems->getWorkingRadius() == 0.8) {
        if( debugInfo == 1 ){ std::cout << heatflux << "," << equilibriumWallTemperature << "," << adiabaticWallTemperature << "," << currentAirdensity << "," << currentAirspeed << std::endl; }
    }

    return heatflux;
}


double findEquilibriumWallTemperature(
        const std::function< double( const double ) > &heatTransferFunction,
        const double &wallEmissivity,
        const double &adiabaticWallTemperature,
        const double &currentWallTemperature )
{
    std::function< double( const double ) > equilibriumWallTemperatureRootFindingFunction =
            std::bind( &bislip::Variables::computeEquilibiumWallTemperatureRootFinder, heatTransferFunction, wallEmissivity, std::placeholders::_1  );

    return bislip::Variables::rootFinderBisection( equilibriumWallTemperatureRootFindingFunction, 0.0, adiabaticWallTemperature, currentWallTemperature );
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

    double M = 3.15;
    //if( bislipSystems->getValidationFlag( ) == true ) { M = 3.15; }

    double N = 0.5;

    const double C_s = std::pow( curvatureRadius, -0.5 ) * std::pow( 1.225, -N ) * std::pow( 7905, -M ) * ( 119988627.724096 );
    //const double C_s = std::pow( curvatureRadius, -0.5 ) * ( 1.83E-4 );

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

    //double a = 0.0; equilibriumWallTemperatureRootFindingFunction( 0.0 );
    //double b = adiabaticWallTemperature;
    //double root = 0.0;
    //double f_a, f_root;
    /*
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
            //if( debugInfo == 1 ){ std::cout << "     root = " << root << "     |     " << "equilibriumWallTemperatureRootFindingFunction( root )  = " << std::abs( equilibriumWallTemperatureRootFindingFunction( root ) ) << std::endl; }
        }
    }

    */

    double equilibriumWallTemperature = 0.0;
    if ( equilibriumWallTemperatureRootFindingFunction( 0.0 ) * equilibriumWallTemperatureRootFindingFunction( adiabaticWallTemperature ) < 0.0 )
    {
        if( debugInfo == 1 ){ std::cout << "     Starting Flat Plate Eq. wall temp. search:" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "                Bisection Search" << std::endl; }
        equilibriumWallTemperature = bislip::Variables::rootFinderBisection( equilibriumWallTemperatureRootFindingFunction, 0.0, adiabaticWallTemperature, 10.0 );
        if( debugInfo == 1 ){ std::cout << "                Bisection Search Complete" << std::endl; }
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "     Starting Flat Plate Eq. wall temp. search:" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "                Golden Ratio Search" << std::endl; }
        equilibriumWallTemperature = bislip::Variables::goldenSectionSearch( equilibriumWallTemperatureRootFindingFunction, 0.0, adiabaticWallTemperature );
        if( debugInfo == 1 ){ std::cout << "                Golden Ratio Search Complete" << std::endl; }
    }


    bislipSystems->setWallTemperature( equilibriumWallTemperature );

    return heatTransferFunction( equilibriumWallTemperature );


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

    if( debugInfo == 1 ){ std::cout << "        Calculating Tauber Heat Flux for Leading Edges" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "             Calculate Stagnation Heat Flux Component of Tauber Heat Flux" << std::endl; }

    bislipSystems->setWallTemperature( bislipSystems->getTauberWallTempStagnation() );
    //if( debugInfo == 1 ){ std::cout << "     Initial Wall Temperature = " <<  bislipSystems->getTauberWallTempStagnation( ) << std::endl; }

    const double q_dot_s = bislip::Variables::computeStagnationHeatFlux( bodyMap, vehicleName);

    if( debugInfo == 1 ){ std::cout << "             Set Stagnation Heat Flux Component of Tauber Heat Flux" << std::endl; }
    bislipSystems->setTauberHeatFluxStagnation( q_dot_s );
    bislipSystems->setTauberWallTempStagnation( bislipSystems->getWallTemperature() );
    if( debugInfo == 1 ){ std::cout << "                  Current Tauber Stagnation Heat Flux = " <<  bislipSystems->getTauberHeatFluxStagnation( ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                  Current Tauber Stagnation Eq. Temp. = " <<  bislipSystems->getTauberWallTempStagnation( ) << std::endl; }

    if( debugInfo == 1 ){ std::cout << "             Calculate Flat Plate Heat Flux Component of Tauber Heat Flux" << std::endl; }
    bislipSystems->setWallTemperature( bislipSystems->getTauberWallTempFlatPlate() );
    //if( debugInfo == 1 ){ std::cout << "     Initial Wall Temperature = " <<  bislipSystems->getTauberWallTempFlatPlate( ) << std::endl; }
    const double q_dot_FP = ( 100 * 100 ) * bislip::Variables::computeFlatPlateHeatFlux( bodyMap, vehicleName);

    if( debugInfo == 1 ){ std::cout << "             Set Flat Plate Heat Flux Component of Tauber Heat Flux" << std::endl; }
    bislipSystems->setTauberHeatFluxFlatPlate( q_dot_FP );
    bislipSystems->setTauberWallTempFlatPlate( bislipSystems->getWallTemperature() );
    if( debugInfo == 1 ){ std::cout << "                  Current Tauber Flat Plate Heat Flux = " <<  bislipSystems->getTauberHeatFluxFlatPlate( ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "                  Current Tauber Flat Plate Eq. Temp. = " <<  bislipSystems->getTauberWallTempFlatPlate( ) << std::endl; }



    //const double lambda = bislipSystems->getWingSweepAngle();
    const double sin_lambda = std::sin( bislipSystems->getWingSweepAngle() );
    const double cos_lambda = std::cos( bislipSystems->getWingSweepAngle() );

    if( debugInfo == 1 ){ std::cout << "             Calculate Tauber Heat Flux Component of Tauber Heat Flux" << std::endl; }
    double leadingEdgeHeatFlux = std::sqrt( 0.5 * q_dot_s * q_dot_s * cos_lambda * cos_lambda + q_dot_FP * q_dot_FP * sin_lambda * sin_lambda );

    if( debugInfo == 1 ){ std::cout << "                     Set Tauber Heat Flux" << std::endl; }
    bislipSystems->setCurrentHeatFluxTauber( leadingEdgeHeatFlux );

    return leadingEdgeHeatFlux;
}

double computeRadiativeHeatFlux(
        const double &wallEmissivity,
        const double &bodyTemperature )
{
    return wallEmissivity * tudat::physical_constants::STEFAN_BOLTZMANN_CONSTANT * bodyTemperature * bodyTemperature * bodyTemperature * bodyTemperature;
}

//! Function to compute the equilibrium wall temperature from the heat input and emmisivity
double computeEquilibiumWallTemperatureRootFinder(
        const std::function< double( const double ) > &heatTransferFunction,
        const double &wallEmmisivity,
        const double &wallTemperature )

{

    // const double evaluation = heatTransferFunction( wallTemperature ) - bislip::Variables::computeRadiativeHeatFlux( wallEmmisivity, wallTemperature );

    return heatTransferFunction( wallTemperature ) - bislip::Variables::computeRadiativeHeatFlux( wallEmmisivity, wallTemperature );;
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
        const std::string &vehicleName,
        const std::string &centralBodyName )

{
    //! Extract flight conditions pointer.
    std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > flightConditions = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                bodyMap.at( vehicleName )->getFlightConditions( ) );

    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;

    //! Extract current conditions.
    const double currentLatitude            = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
    const double currentLongitude           = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );
    const double currentFlightPathAngle     = flightConditions->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
    const double currentFlightPathAngleRate = bislip::Variables::computeFlightPathAngleRate( bodyMap, vehicleName, centralBodyName );
    const double initialDistanceToTarget    = bislipSystems->getInitialDistanceToTarget();


    //! Calculate angular distance to go.
    const double angularDistanceTravelled = bislip::Variables::computeAngularDistance (
                bodyMap.at( vehicleName )->getBislipSystems()->getInitialLat(),
                bodyMap.at( vehicleName )->getBislipSystems()->getInitialLon(),
                currentLatitude,
                currentLongitude );

    bool stop = false;

    if ( currentFlightPathAngle < 0.0 )
    {
        if ( currentFlightPathAngleRate < 0.0 )
        {
            if ( ( angularDistanceTravelled / initialDistanceToTarget ) > bislipSystems->getAscentTerminationDistanceRatio() ) { stop = true; }
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
                       const double E_hat = bislip::Variables::computeNormalizedSpecificEnergy( current_height, current_V, bodyMap.at( vehicleName )->getBislipSystems()->getMaximumSpecificEnergy() );

                       //! Evaluate current throttle setting and thrust elevation angle.
                       double throttle = bislip::Variables::evaluateGuidanceInterpolator(
                                       //current_gamma,
                                       bislip::Parameters::Interpolators::ThrottleSetting,
                                       bodyMap.at( vehicleName )->getBislipSystems(),
                                       current_height,
                                       current_V,
                                       bodyMap.at( vehicleName )->getBislipSystems()->getMaximumSpecificEnergy() );

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

        double constraintViolationPenalty = ( constraintViolations.sum() ) * ( propagationStepSize / normalizer );

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


double estimatedFlightPathAngle (
        const Eigen::Vector3d &aerodynamicFrameTotalLoad,
        const double &currentMass )
{





    return std::asin( aerodynamicFrameTotalLoad( 0 ) / ( currentMass * tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );
}


int determineNumberOfHardViolations (
        const Eigen::VectorXd &depVar,
        const double &hardConstraint )
{
    int numberOfViolations = 0;

    for( int i = 0; i < depVar.size(); i++ ) { if( depVar( i ) > hardConstraint ) { numberOfViolations += 1; } }

    return numberOfViolations;
}




bool convertTrajectoryPhaseToBoolean(
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const std::string &vehicleName )
{
    //! Extract Bislip Systems pointer.
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( vehicleName )->getBislipSystems( ) ;


    bool phaseInBoolean = false;

    if( bislipSystems->getCurrentTrajectoryPhase() != "Ascent" )
    {
        phaseInBoolean = true;
    }


    return phaseInBoolean;
}









}; //namespace Variables
                 } // namespace bislip










