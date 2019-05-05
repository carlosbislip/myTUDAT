#include "Space4ErrBody_Headers/Space4ErrBodyProblem.h"

#include "Space4ErrBody_Headers/getPagmoAlgorithm.h"
//#include "Space4ErrBody_Headers/saveOptimizationResults.h"

using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::aerodynamics;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::unit_conversions;
using namespace tudat;
//using namespace tudat_applications;
using namespace pagmo;
using namespace bislip;
using namespace bislip::Variables;
using namespace bislip::Parameters;
using namespace bislip::Utilities;
using namespace boost;
//using namespace tudat_pagmo_applications;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////            Main Program                  //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    //! Create shared pointer for problem input.
    std::shared_ptr< bislip::ProblemInput > problemInput = std::make_shared< bislip::ProblemInput >( );

    //! Get this run's time stamp. Will be used to create a unique output
    //! subfolder name, where all files created by the optimizer will be stored.
    std::string playTime = bislip::Variables::getCurrentDateTime( false );

    problemInput->setPlayTime( playTime );



    std::string primaryFileName;
    std::cout << "Primary File Name: ";
    std::getline ( std::cin, primaryFileName );


    //! Get all inputs.
    const std::vector< std::string > primary                   = bislip::Variables::getDataString ( primaryFileName + ".txt" );
    const std::vector< std::string > ascentParameterList       = bislip::Variables::getDataString ( primary[ 0 ] );
    const std::vector< double > ascentParameterBounds          = bislip::Variables::getDataNumeri ( primary[ 1 ] );
    const std::vector< std::string > descentParameterList      = bislip::Variables::getDataString ( primary[ 2 ] );
    const std::vector< double > descentParameterBounds         = bislip::Variables::getDataNumeri ( primary[ 3 ] );
    const std::vector< std::string > vehicleParameterList      = bislip::Variables::getDataString ( primary[ 4 ] );
    const std::vector< double > vehicleParameterValues         = bislip::Variables::getDataNumeri ( primary[ 5 ] );
    const std::vector< std::string > aeroCoeffFileList         = bislip::Variables::getDataString ( primary[ 6 ] );
    const std::vector< std::string > optimizationSettingsList  = bislip::Variables::getDataString ( primary[ 7 ] );
    const std::vector< double > optimizationSettingsValues     = bislip::Variables::getDataNumeri ( primary[ 8 ] );
    const std::vector< std::string > simulationSettingsList    = bislip::Variables::getDataString ( primary[ 9 ] );
    const std::vector< double > simulationSettingsValues       = bislip::Variables::getDataNumeri ( primary[ 10 ] );
    const std::vector< std::string > outputSettingsList        = bislip::Variables::getDataString ( primary[ 11 ] );
    const std::vector< double > outputSettingsValues           = bislip::Variables::getDataNumeri ( primary[ 12 ] );
    const std::vector< std::string > initialConditionsList     = bislip::Variables::getDataString ( primary[ 13 ] );
    const std::vector< double > initialConditionsValues        = bislip::Variables::getDataNumeri ( primary[ 14 ] );
    const std::vector< std::string > constraintsList           = bislip::Variables::getDataString ( primary[ 15 ] );
    const std::vector< double > constraintsValues              = bislip::Variables::getDataNumeri ( primary[ 16 ] );
    const std::vector< double > headingErrorDeadBand           = bislip::Variables::getDataNumeri ( primary[ 17 ] );
    const std::vector< double > alphaMachEnvelopeLB            = bislip::Variables::getDataNumeri ( primary[ 18 ] );
    const std::vector< double > alphaMachEnvelopeUB            = bislip::Variables::getDataNumeri ( primary[ 19 ] );
    const std::string problemName                              = primary.rbegin()[ 1 ];
    const std::string vehicleName                              = primary.rbegin()[ 0 ];

    problemInput->setProblemName( problemName );
    problemInput->setOutputSettings( outputSettingsValues );
    problemInput->setSimulationSettings( simulationSettingsValues );
    problemInput->setVehicleName( vehicleName );
    problemInput->setVehicleParameterList( vehicleParameterList );
    problemInput->setVehicleParameters( vehicleParameterValues );
    problemInput->setAerodynamicDatabaseFileList( aeroCoeffFileList );
    problemInput->setAscentParameterList( ascentParameterList );
    problemInput->setAscentParameterBounds( ascentParameterBounds );
    problemInput->setDescentParameterList( descentParameterList );
    problemInput->setDescentParameterBounds( descentParameterBounds );
    problemInput->setInitialConditions( initialConditionsValues );
    problemInput->setConstraints( constraintsValues );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            UNPACK INPUT DATA             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //std::cout << "Unpacking data" << std::endl;

    //! Declare and initialize simulation start epoch.
    const double simulationStartEpoch = simulationSettingsValues[ 0 ]; // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000

    //! Declare and initialize maximum simulation time of flight.
    const double max_tof = simulationSettingsValues[ 1 ];

    //! Declare and initialize simulation end epoch.
    const double simulationEndEpoch = simulationStartEpoch + max_tof;

    //! Declare and initialize numerical integration fixed step size.
    const double propagationStepSize = simulationSettingsValues[ 2 ];

    //! Declare and initialize propagation to guidance sampling ratio.
    const double guidanceStepSize = simulationSettingsValues[ 3 ];

    const int debugInfo = int( simulationSettingsValues[ 4 ] );

    //! Declare and initialize number of control nodes.
    //const unsigned long nodesAscent = simulationSettingsValues[ 3 ];
    //   const unsigned long nodesDescent = simulationSettingsValues[ 4 ];

    //! Declare and initialize Reference area
    const double S_ref = vehicleParameterValues[ 0 ]; // m^2

    //! Declare and initialize Reference length: wingspan
    const double b_ref = vehicleParameterValues[ 1 ]; // m

    //! Declare and initialize Reference length: chord
    const double c_ref = vehicleParameterValues[ 2 ]; // m

    const Eigen::Vector3d referenceValues( S_ref, b_ref, c_ref ); // m^2, m, m

    //! Declare and initialize position vector of moment reference center
    const Eigen::Vector3d R_mrc( vehicleParameterValues[ 3 ], vehicleParameterValues[ 4 ], vehicleParameterValues[ 5 ] ); // m

    //! Declare and initialize position vector of center of mass
    const Eigen::Vector3d R_com( vehicleParameterValues[ 6 ], vehicleParameterValues[ 7 ], vehicleParameterValues[ 8 ] ); // m

    //! Declare and initialize position vector of center of thrust
    const Eigen::Vector3d R_cot( vehicleParameterValues[ 9 ], vehicleParameterValues[ 10 ], vehicleParameterValues[ 11 ] ); // m

    //! Declare and initialize initial mass
    const double initialMass_Ascent = vehicleParameterValues[ 12 ]; // kg

    //! Declare and initialize landing mass
    const double landingMass = vehicleParameterValues[ 13 ]; // kg

    //! Declare and initialize specific impulse
    const double specificImpulse = vehicleParameterValues[ 14 ]; // kg

    //! Declare and initialize maximum engine thrust
    const double maxThrust = vehicleParameterValues[ 15 ]; // N

    //! Declare and initialize relative chest forward angle
    const double relativeChestForwardAngle_deg = vehicleParameterValues[ 16 ]; // deg

    //! Declare and initialize vertebral column inclination angle
    const double vertebralColumnInclinationAngle_deg = vehicleParameterValues[ 17 ]; // deg

    //! Declare and initialize the minimum dynamic pressure required for bank angle assignment.
    const double minimumDynamicPressureforControlSurface = vehicleParameterValues[ 18 ];

    //! Declare and initialize nose radius
    const double noseRadius = vehicleParameterValues[ 19 ]; // m

    //! Declare and initialize nose radius
    const double leadingEdgeRadius = vehicleParameterValues[ 20 ]; // m

    //! Declare and initialize starting height
    const double initialHeight = initialConditionsValues[ 2 ]; // m

    //! Declare and initialize initial Mach
    const double Mach_i = initialConditionsValues[ 3 ]; // -

    //! Declare and initialize initial flight-path angle
    const double gamma_i_deg = initialConditionsValues[ 4 ]; // deg

    //! Declare and initialize initial bodyflap angle
    const double bodyflap_i_deg = initialConditionsValues[ 5 ]; // deg
    const double elevon_i_deg   = initialConditionsValues[ 6 ]; // deg




    //! Declare and initialize starting position coordinates.
    const double initialLat_deg = initialConditionsValues[ 0 ];
    const double initialLon_deg = initialConditionsValues[ 1 ];

    //! Declare and initialize final position coordinates and additional termination conditions
    const double targetLat_deg = constraintsValues[ 0 ];
    const double targetLon_deg = constraintsValues[ 1 ];

    //! Declare and initialize various termination conditions
    const double finalDistanceToTarget_deg = constraintsValues[ 2 ];
    const double maximumHeightAllowable    = constraintsValues[ 3 ];
    const double minimumHeightAllowable    = constraintsValues[ 4 ];

    //const double V_UP = constraintsValues[ 5 ];
    const double constraint_MechanicalLoad         = constraintsValues[ 6 ];
    const double constraint_HeatingRate            = constraintsValues[ 7 ];
    const double constraint_DynamicPressure        = constraintsValues[ 8 ];
    const double constraint_PitchMomentCoefficient = constraintsValues[ 9 ];

    //! Convert angles from degrees to radians
    const double initialLat_rad                      = unit_conversions::convertDegreesToRadians( initialLat_deg );
    const double initialLon_rad                      = unit_conversions::convertDegreesToRadians( initialLon_deg );
    const double targetLat_rad                       = unit_conversions::convertDegreesToRadians( targetLat_deg );
    const double targetLon_rad                       = unit_conversions::convertDegreesToRadians( targetLon_deg );
    const double finalDistanceToTarget_rad           = unit_conversions::convertDegreesToRadians( finalDistanceToTarget_deg );
    const double initialFlightPathAngle_rad          = unit_conversions::convertDegreesToRadians( gamma_i_deg );
    const double relativeChestForwardAngle_rad       = unit_conversions::convertDegreesToRadians( relativeChestForwardAngle_deg );
    const double vertebralColumnInclinationAngle_rad = unit_conversions::convertDegreesToRadians( vertebralColumnInclinationAngle_deg );

    //! Pre-define various variables used to determine fitness.
    double targetLat_deg_calc          = initialLat_deg;
    double targetLon_deg_calc          = initialLon_deg;
    double initialDistanceToTarget_rad = bislip::Variables::computeAngularDistance( initialLat_rad, initialLon_rad, targetLat_rad, targetLon_rad );
    double initialDistanceToTarget_deg = unit_conversions::convertRadiansToDegrees( initialDistanceToTarget_rad );

    //! Still working on these
    //const double noseRadius = 3.0;
    //const double leadingEdgeRadius = 1.5;
    const double lambda = unit_conversions::convertDegreesToRadians( 30.0 );
    const double epsilon = 0.7;
    const double x_T = 1E-9;
    const double phi = unit_conversions::convertDegreesToRadians( 3.0 );
    Eigen::Matrix3d bodyFrameToPassengerFrameTransformationMatrix =
            ( bislip::Variables::computeRotationMatrixONE( -( relativeChestForwardAngle_rad + tudat::mathematical_constants::PI ) ) ) * ( bislip::Variables::computeRotationMatrixTHREE( -( vertebralColumnInclinationAngle_rad + tudat::mathematical_constants::PI ) ) );

    if( debugInfo == 1 ){ std::cout << "Creating partial output subfolder filename suffix" << std::endl; }
    //! Create partial output subfolder filename suffix based on optimization
    //! settings and fixed time step.
    const std::string this_run_settings = std::to_string( int( optimizationSettingsValues[ 0 ] ) ) + "_" +
            std::to_string( int( optimizationSettingsValues[ 1 ] ) ) + "_" +
            std::to_string( int( optimizationSettingsValues[ 2 ] ) ) + "_" +
            std::to_string( int( optimizationSettingsValues[ 3 ] ) ) + "_" +
            std::to_string( int( optimizationSettingsValues[ 4 ] ) ) + "_" +
            std::to_string( simulationSettingsValues[ 2 ] );

    //! Create output subfolder filename. Based on arbitrary prefix, the
    //! previously created partial suffix, and the run's time stamp.
    const std::string outputSubFolder = "OUTPUT_" + this_run_settings + "_" + playTime + "/";

    const std::string outputPath = bislip::Variables::getOutputPath( );

    problemInput->setOutputPath( outputPath );
    problemInput->setOutputSubFolder( outputSubFolder );

    //! Determine number of Ascent parameters to vary based of size of ascentParameterBounds vector.
    const unsigned long nodesAscent = simulationSettingsValues.rbegin()[ 1 ];
    const int N = ( ascentParameterList.size() - 6 ) * nodesAscent - 1;
    const int M = ascentParameterBounds.size();

    //! Determine number of Descent parameters to vary based of size of descentParameterBounds vector.
    const unsigned long nodesDescent = simulationSettingsValues.rbegin()[ 0 ];
    const int NN = ( descentParameterList.size() - 2 ) * nodesDescent - 1;
    const int MM = descentParameterBounds.size();

    //! Create vector containing decision vector bounds.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( N + 6 + NN + 2, 0.0 ) );
    //!---------------------------------------   ^ for lb/up  (rows)       ^ for # of parameters


    if( debugInfo == 1 ){ std::cout << "Constructing Bounds Matrix" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Decision Vector Bounds" << std::endl; }
    //! Loop to build the bounds matrix.
    unsigned long p = 0;

    for( int i = 0; i < N; i++ )
    {
        if( i < ( nodesAscent - 1 ) )
        {
            p = 0;
        }
        else if( ( i - ( nodesAscent - 1 ) ) % nodesAscent == 0 )
        {
            p += 2;
        }

        bounds[ 0 ][ i ] = ascentParameterBounds[ p ];
        bounds[ 1 ][ i ] = ascentParameterBounds[ p + 1 ];

        if( debugInfo == 1 )
        {
            std::cout << "bounds[ 0 ][ " << i << " ] = " << ascentParameterBounds[ p ] << std::endl;
            std::cout << "bounds[ 1 ][ " << i << " ] = " << ascentParameterBounds[ p + 1 ] << std::endl;
        }
    }

    bounds[ 0 ][ N ] = ascentParameterBounds[ M - 12 ];
    bounds[ 1 ][ N ] = ascentParameterBounds[ M - 11 ];
    bounds[ 0 ][ N + 1 ] = ascentParameterBounds[ M - 10 ];
    bounds[ 1 ][ N + 1 ] = ascentParameterBounds[ M - 9 ];
    bounds[ 0 ][ N + 2 ] = ascentParameterBounds[ M - 8 ];
    bounds[ 1 ][ N + 2 ] = ascentParameterBounds[ M - 7 ];
    bounds[ 0 ][ N + 3 ] = ascentParameterBounds[ M - 6 ];
    bounds[ 1 ][ N + 3 ] = ascentParameterBounds[ M - 5 ];
    bounds[ 0 ][ N + 4 ] = ascentParameterBounds[ M - 4 ];
    bounds[ 1 ][ N + 4 ] = ascentParameterBounds[ M - 3 ];
    bounds[ 0 ][ N + 5 ] = ascentParameterBounds[ M - 2 ];
    bounds[ 1 ][ N + 5 ] = ascentParameterBounds[ M - 1 ];

    if( debugInfo == 1 )
    {
        std::cout << "bounds[ 0 ][ " << N << " ] = " << ascentParameterBounds[ M - 12 ] << std::endl;
        std::cout << "bounds[ 1 ][ " << N << " ] = " << ascentParameterBounds[ M - 11 ] << std::endl;
        std::cout << "bounds[ 0 ][ " << N + 1 << " ] = " << ascentParameterBounds[ M - 10 ] << std::endl;
        std::cout << "bounds[ 1 ][ " << N + 1 << " ] = " << ascentParameterBounds[ M - 9 ] << std::endl;
        std::cout << "bounds[ 0 ][ " << N + 2 << " ] = " << ascentParameterBounds[ M - 8 ] << std::endl;
        std::cout << "bounds[ 1 ][ " << N + 2 << " ] = " << ascentParameterBounds[ M - 7 ] << std::endl;
        std::cout << "bounds[ 0 ][ " << N + 3 << " ] = " << ascentParameterBounds[ M - 6 ] << std::endl;
        std::cout << "bounds[ 1 ][ " << N + 3 << " ] = " << ascentParameterBounds[ M - 5 ] << std::endl;
        std::cout << "bounds[ 0 ][ " << N + 4 << " ] = " << ascentParameterBounds[ M - 4 ] << std::endl;
        std::cout << "bounds[ 1 ][ " << N + 4 << " ] = " << ascentParameterBounds[ M - 3 ] << std::endl;
        std::cout << "bounds[ 0 ][ " << N + 5 << " ] = " << ascentParameterBounds[ M - 2 ] << std::endl;
        std::cout << "bounds[ 1 ][ " << N + 5 << " ] = " << ascentParameterBounds[ M - 1 ] << std::endl;
    }

    for( int i = 0; i < NN ; i++ )
    {
        if( i < ( nodesDescent -1 ) )
        {
            p = 0;
        }
        else if( ( i - ( nodesDescent - 1 ) ) % nodesDescent == 0 )
        {
            p += 2;
        }

        bounds[ 0 ][ ( N + 6 ) + i ] = descentParameterBounds[ p ];
        bounds[ 1 ][ ( N + 6 ) + i ] = descentParameterBounds[ p + 1 ];

        if( debugInfo == 1 )
        {

            std::cout << "bounds[ 0 ][ " << ( N + 6 ) + i << " ] = " << descentParameterBounds[ p ] << std::endl;
            std::cout << "bounds[ 1 ][ " << ( N + 6 ) + i << " ] = " << descentParameterBounds[ p + 1 ] << std::endl;
        }
    }


    bounds[ 0 ][ ( N + 6 ) + NN ]     = descentParameterBounds[ MM - 4 ];
    bounds[ 1 ][ ( N + 6 ) + NN ]     = descentParameterBounds[ MM - 3 ];
    bounds[ 0 ][ ( N + 6 ) + NN + 1 ] = descentParameterBounds[ MM - 2 ];
    bounds[ 1 ][ ( N + 6 ) + NN + 1 ] = descentParameterBounds[ MM - 1 ];

    if( debugInfo == 1 )
    {
        std::cout << "bounds[ 0 ][ " << ( N + 6 ) + NN << " ] = " << descentParameterBounds[ MM - 4 ] << std::endl;
        std::cout << "bounds[ 1 ][ " << ( N + 6 ) + NN << " ] = " << descentParameterBounds[ MM - 3 ] << std::endl;
        std::cout << "bounds[ 0 ][ " << ( N + 6 ) + NN + 1 << " ] = " << descentParameterBounds[ MM - 2 ] << std::endl;
        std::cout << "bounds[ 1 ][ " << ( N + 6 ) + NN + 1 << " ] = " << descentParameterBounds[ MM - 1 ] << std::endl;

        std::cout << "N = " << N << std::endl;
        std::cout << "M = " << M << std::endl;
        std::cout << "NN = " << NN << std::endl;
        std::cout << "MM = " << MM << std::endl;
        //std::cout << "descentParameterBounds[ MM - 2 ]: " << descentParameterBounds[ MM - 2 ] << std::endl;
        //std::cout << "descentParameterBounds[ MM - 1 ]: " << descentParameterBounds[ MM - 1 ] << std::endl;
        //std::cout << "N + NN + 3 + 1: " << N + NN + 3 + 1 << std::endl;
    }


    problemInput->setDecisionVectorBounds( bounds );

    //! Load spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE PARAMETER BOUNDS          ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Creating Parameter Bounds" << std::endl; }

    //! Declare data maps of strings and pairs to store parameter bounds.
    //!     Strings are used for the parameter names.
    //!     Pairs are used for lower and upper bounds.
    std::map< bislip::Parameters::Optimization, std::pair < double, double > > ascentParameterBoundsMap, descentParameterBoundsMap;

    //! Declare and initialize counter.
    //!     Information is extracted directly from vectors created from input data, which requires the counter for the upper bound value.
    p = 0;

    //! Loop to populate the data map with pairs associated to string values related to Ascent phase.
    for( unsigned long i = 0; i < ascentParameterList.size(); i++)
    {
        if ( ascentParameterList[ i ] == "Angle of Attack" ) { ascentParameterBoundsMap[ bislip::Parameters::Optimization::AngleOfAttack ]               = std::make_pair( ascentParameterBounds[ p ], ascentParameterBounds[ p + 1 ] ); }
        if ( ascentParameterList[ i ] == "Bank Angle" ) { ascentParameterBoundsMap[ bislip::Parameters::Optimization::BankAngle ]                        = std::make_pair( ascentParameterBounds[ p ], ascentParameterBounds[ p + 1 ] ); }
        if ( ascentParameterList[ i ] == "Thrust Elevation Angle" ) { ascentParameterBoundsMap[ bislip::Parameters::Optimization::ThrustElevationAngle ] = std::make_pair( ascentParameterBounds[ p ], ascentParameterBounds[ p + 1 ] ); }
        if ( ascentParameterList[ i ] == "Thrust Azimuth Angle" ) { ascentParameterBoundsMap[ bislip::Parameters::Optimization::ThrustAzimuthAngle ]     = std::make_pair( ascentParameterBounds[ p ], ascentParameterBounds[ p + 1 ] ); }
        if ( ascentParameterList[ i ] == "Throttle Setting" ) { ascentParameterBoundsMap[ bislip::Parameters::Optimization::ThrottleSetting ]            = std::make_pair( ascentParameterBounds[ p ], ascentParameterBounds[ p + 1 ] ); }
        if ( ascentParameterList[ i ] == "Node Interval" ) { ascentParameterBoundsMap[ bislip::Parameters::Optimization::NodeInterval ]                  = std::make_pair( ascentParameterBounds[ p ], ascentParameterBounds[ p + 1 ] ); }
        if ( ascentParameterList[ i ] == "Initial Velocity" ) { ascentParameterBoundsMap[ bislip::Parameters::Optimization::InitialVelocity ]            = std::make_pair( ascentParameterBounds[ p ], ascentParameterBounds[ p + 1 ] ); }
        if ( ascentParameterList[ i ] == "Maximum Velocity" ) { ascentParameterBoundsMap[ bislip::Parameters::Optimization::MaximumVelocity ]            = std::make_pair( ascentParameterBounds[ p ], ascentParameterBounds[ p + 1 ] ); }
        if ( ascentParameterList[ i ] == "Maximum Height" ) { ascentParameterBoundsMap[ bislip::Parameters::Optimization::MaximumHeight ]                = std::make_pair( ascentParameterBounds[ p ], ascentParameterBounds[ p + 1 ] ); }
        if ( ascentParameterList[ i ] == "Additional Mass" ) { ascentParameterBoundsMap[ bislip::Parameters::Optimization::AdditionalMass ]              = std::make_pair( ascentParameterBounds[ p ], ascentParameterBounds[ p + 1 ] ); }
        p += 2;
    }

    problemInput->setAscentParameterBoundsMap( ascentParameterBoundsMap );

    //! Re-initialize counter.
    p = 0;

    //! Loop to populate the data map with pairs associated to string values related to Descent phase.
    for( unsigned long i = 0; i < descentParameterList.size(); i++)
    {
        if ( descentParameterList[ i ] == "Angle of Attack" ) { descentParameterBoundsMap[ bislip::Parameters::Optimization::AngleOfAttack ]               = std::make_pair( descentParameterBounds[ p ], descentParameterBounds[ p + 1 ] ); }
        if ( descentParameterList[ i ] == "Bank Angle" ) { descentParameterBoundsMap[ bislip::Parameters::Optimization::BankAngle ]                        = std::make_pair( descentParameterBounds[ p ], descentParameterBounds[ p + 1 ] ); }
        if ( descentParameterList[ i ] == "Thrust Elevation Angle" ) { descentParameterBoundsMap[ bislip::Parameters::Optimization::ThrustElevationAngle ] = std::make_pair( descentParameterBounds[ p ], descentParameterBounds[ p + 1 ] ); }
        if ( descentParameterList[ i ] == "Thrust Azimuth Angle" ) { descentParameterBoundsMap[ bislip::Parameters::Optimization::ThrustAzimuthAngle ]     = std::make_pair( descentParameterBounds[ p ], descentParameterBounds[ p + 1 ] ); }
        if ( descentParameterList[ i ] == "Throttle Setting" ) { descentParameterBoundsMap[ bislip::Parameters::Optimization::ThrottleSetting ]            = std::make_pair( descentParameterBounds[ p ], descentParameterBounds[ p + 1 ] ); }
        if ( descentParameterList[ i ] == "Node Interval" ) { descentParameterBoundsMap[ bislip::Parameters::Optimization::NodeInterval ]                  = std::make_pair( descentParameterBounds[ p ], descentParameterBounds[ p + 1 ] ); }
        if ( descentParameterList[ i ] == "Initial Velocity" ) { descentParameterBoundsMap[ bislip::Parameters::Optimization::InitialVelocity ]            = std::make_pair( descentParameterBounds[ p ], descentParameterBounds[ p + 1 ] ); }
        if ( descentParameterList[ i ] == "Maximum Velocity" ) { descentParameterBoundsMap[ bislip::Parameters::Optimization::MaximumVelocity ]            = std::make_pair( descentParameterBounds[ p ], descentParameterBounds[ p + 1 ] ); }
        if ( descentParameterList[ i ] == "Maximum Height" ) { descentParameterBoundsMap[ bislip::Parameters::Optimization::MaximumHeight ]                = std::make_pair( descentParameterBounds[ p ], descentParameterBounds[ p + 1 ] ); }
        if ( descentParameterList[ i ] == "Additional Mass" ) { descentParameterBoundsMap[ bislip::Parameters::Optimization::AdditionalMass ]              = std::make_pair( descentParameterBounds[ p ], descentParameterBounds[ p + 1 ] ); }
        p += 2;
    }

    problemInput->setDescentParameterBoundsMap( descentParameterBoundsMap );

    if( debugInfo == 1 ){ std::cout << "Parameter Bounds Created" << std::endl; }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Creating environment" << std::endl; }

    //! Declare and initialize central body name.
    const std::string centralBodyName = "Earth";

    //! Declare and initialize simulation body settings data map.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { centralBodyName },
                                    simulationStartEpoch - 10000.0 * propagationStepSize,
                                    simulationEndEpoch + 10000.0 * propagationStepSize );

    //    std::cout << "Define atmospheric model." << std::endl;
    //! Define atmospheric model.
    bodySettings[ centralBodyName ]->atmosphereSettings = std::make_shared< AtmosphereSettings >(
                nrlmsise00 );

    // std::cout << "Define ephemeris model settings." << std::endl;
    //! Define ephemeris model settings.
    //! This is an acceptable 'cheat' were Earth is placed at the barycenter.
    //! Use only when there arent any third body perturbations (Moon, Sun, etc.)
    bodySettings[ centralBodyName ]->ephemerisSettings =
            std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );

    // std::cout << "Reset ephemeris to J2000." << std::endl;
    //! Reset ephemeris to J2000.
    bodySettings[ centralBodyName ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    //! Define Earth's radius. Using spice here. Is there a way to get a 'radius
    //! field'? Im interested in includig Earth's flattening, yet am unsure how
    //! to properly do it.
    const double radiusEarth = spice_interface::getAverageRadius( centralBodyName );
    double radiusEarth_i = radiusEarth;
    double radiusEartmaximumHeightAllowable = radiusEarth;

    //!---------------------------- Still playing around with this section

    //! https://en.wikipedia.org/wiki/Flattening#Numerical_values_for_planets
    //! https://www.vcalc.com/wiki/vCalc/WGS-84+Earth+flattening+factor
    double flattening_Earth = 1 / 1;//nput_data_.back();
    if ( flattening_Earth != 1 )
    {
        bodySettings[ centralBodyName ]->shapeModelSettings = std::make_shared< OblateSphericalBodyShapeSettings >( radiusEarth, flattening_Earth );

        // Declare variable in which raw result is to be put by Spice function.
        double radii[ 3 ];

        // Call Spice function to retrieve gravitational parameter.
        SpiceInt numberOfReturnedParameters;
        bodvrd_c( "Earth", "RADII", 3, &numberOfReturnedParameters, radii );
        radii[0] = radii [0] * 1e3 ;
        radii[1] = radii [1] * 1e3 ;
        radii[1] = radii [2] * 1e3 ;
        // Compute average and convert from km to m.
        //return unit_conversions::convertKilometersToMeters< double >(
        //           radii[ 0 ] + radii[ 1 ] + radii[ 2 ] ) / 3.0;

        const double radius_Equator =  6378137.0;
        const double radius_pole = 6356752.3142;
        /* const double polar_angle_i_deg = 90 - lat_i_deg;
        const double polar_angle_f_deg = 90 - lat_f_deg;
        const double polar_angle_i_rad = polar_angle_i_deg * mathematical_constants::PI / 180;
        const double polar_angle_f_rad = polar_angle_f_deg * mathematical_constants::PI / 180;

        const double radiusEarth_i_x = radii[ 0 ] * std::cos(polar_angle_i_rad) * std::cos(initialLon_rad);
        const double radiusEarth_i_y = radii[ 1 ] * std::cos(polar_angle_i_rad) * std::sin(initialLon_rad);
        const double radiusEarth_i_z = radii[ 2 ] * std::sin(polar_angle_i_rad);
        radiusEarth_i = std::sqrt(pow(radiusEarth_i_x,2) + pow(radiusEarth_i_y,2) + pow(radiusEarth_i_z,2));
        const double radiusEartmaximumHeightAllowable_x = radii[ 0 ] * std::cos(polar_angle_f_rad) * std::cos(targetLon_rad);
        const double radiusEartmaximumHeightAllowable_y = radii[ 1 ] * std::cos(polar_angle_f_rad) * std::sin(targetLon_rad);
        const double radiusEartmaximumHeightAllowable_z = radii[ 2 ] * std::sin(polar_angle_f_rad);
        radiusEartmaximumHeightAllowable = std::sqrt(pow(radiusEartmaximumHeightAllowable_x,2) + pow(radiusEartmaximumHeightAllowable_y,2) + pow(radiusEartmaximumHeightAllowable_z,2));
*/
    }
    //!--------------------------------------------------------------

    if( debugInfo == 1 ){ std::cout << "Creating Earth object" << std::endl; }

    //! Create Earth object
    tudat::simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    //std::cout << "Create Earth's rotational ephemeris" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Creating vehicle" << std::endl; }

    ///////// Start: Vehicle Aerodynamics Section
    if( debugInfo == 1 ){ std::cout << "           Start: Vehicle Aerodynamics Section" << std::endl; }

    //! Define physical meaning of independent variables
    //!     Angle of attack
    //!     Mach number
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames;
    independentVariableNames.push_back( aerodynamics::angle_of_attack_dependent );
    independentVariableNames.push_back( aerodynamics::mach_number_dependent );

    // Define physical meaning of independent variables for control surface increments, in this case Mach number, angle of attack and control surface deflection
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > controlSurfaceIndependentVariableNames;
    controlSurfaceIndependentVariableNames.push_back( aerodynamics::angle_of_attack_dependent );
    controlSurfaceIndependentVariableNames.push_back( aerodynamics::mach_number_dependent );
    controlSurfaceIndependentVariableNames.push_back( aerodynamics::control_surface_deflection_dependent );


    //! Define list of files for force coefficients. Entry denotes direction.
    //!     0 : x-direction (C ~D~/C ~X~)
    //!     1 : y-direction (C ~S~/C ~Y~)
    //!     2 : z-direction (C ~L~/C ~Z~)
    std::map< int, std::string > forceCoefficientFiles;
    std::map< int, std::string > forceCoefficientFiles_CS_B;
    std::map< int, std::string > forceCoefficientFiles_CS_EL;
    std::map< int, std::string > forceCoefficientFiles_CS_ER;
    std::map< int, std::string > momentCoefficientFiles;
    std::map< int, std::string > momentCoefficientFiles_CS_B;
    std::map< int, std::string > momentCoefficientFiles_CS_EL;
    std::map< int, std::string > momentCoefficientFiles_CS_ER;

    std::string BODYFLAP = "BodyFlap";
    std::string ELEVON_L = "ElevonLeft";
    std::string ELEVON_R = "ElevonRight";

    if( debugInfo == 1 ){ std::cout << "           Assigning files names for Aerodynamic Coefficient Database" << std::endl; }

    // Define list of files for force coefficients.
    forceCoefficientFiles[ 0 ] = aeroCoeffFileList[ 0 ]; // Set drag coefficient file
    forceCoefficientFiles[ 2 ] = aeroCoeffFileList[ 1 ]; // Set lift coefficient file

    // Define list of files for moment coefficients.
    momentCoefficientFiles[ 1 ] = aeroCoeffFileList[ 2 ]; // Set pitch moment coefficient file

    // Define list of files for force coefficients for control surfaces ( Bodyflap )
    forceCoefficientFiles_CS_B[ 0 ] = aeroCoeffFileList[ 3 ];
    forceCoefficientFiles_CS_B[ 2 ] = aeroCoeffFileList[ 4 ];

    // Define list of files for moment coefficients for control surfaces ( Bodyflap )
    momentCoefficientFiles_CS_B[ 1 ] = aeroCoeffFileList[ 5 ];

    // Define list of files for force coefficients for control surfaces ( Elevon - Left )
    forceCoefficientFiles_CS_EL[ 0 ] = aeroCoeffFileList[ 6 ];
    forceCoefficientFiles_CS_EL[ 2 ] = aeroCoeffFileList[ 7 ];

    // Define list of files for moment coefficients for control surfaces ( Elevon - Left )
    momentCoefficientFiles_CS_EL[ 1 ] = aeroCoeffFileList[ 8 ];

    // Define list of files for force coefficients for control surfaces ( Elevon - Right )
    forceCoefficientFiles_CS_ER[ 0 ] = aeroCoeffFileList[ 9 ];
    forceCoefficientFiles_CS_ER[ 2 ] = aeroCoeffFileList[ 10 ];

    // Define list of files for moment coefficients for control surfaces ( Elevon - Right )
    momentCoefficientFiles_CS_ER[ 1 ] = aeroCoeffFileList[ 11 ];

    //! Define reference frame in which the loaded coefficients are defined.
    //! Have to get some more background info here to properly understand it.
    bool areCoefficientsInAerodynamicFrame = true;
    bool areCoefficientsInNegativeAxisDirection = true;

    if( debugInfo == 1 ){ std::cout << "           Load and parse coefficient files; create coefficient settings" << std::endl; }

    //! Load and parse coefficient files; create coefficient settings.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            tudat::simulation_setup::readTabulatedAerodynamicCoefficientsFromFiles(
                forceCoefficientFiles,
                momentCoefficientFiles,
                b_ref,
                S_ref,
                c_ref,
                R_mrc,
                independentVariableNames,
                areCoefficientsInAerodynamicFrame,
                areCoefficientsInNegativeAxisDirection );

    if( debugInfo == 1 ){ std::cout << "           Add settings for control surface increments to main aerodynamic coefficients" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "           --------> Bodyflap" << std::endl; }

    //! Add settings for control surface increments to main aerodynamic coefficients
    aerodynamicCoefficientSettings->setControlSurfaceSettings(
                tudat::simulation_setup::readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                    forceCoefficientFiles_CS_B,
                    momentCoefficientFiles_CS_B,
                    controlSurfaceIndependentVariableNames), BODYFLAP );

    if( debugInfo == 1 ){ std::cout << "           --------> Left Elevon" << std::endl; }

    aerodynamicCoefficientSettings->setControlSurfaceSettings(
                tudat::simulation_setup::readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                    forceCoefficientFiles_CS_EL,
                    momentCoefficientFiles_CS_EL,
                    controlSurfaceIndependentVariableNames), ELEVON_L );

    if( debugInfo == 1 ){ std::cout << "           --------> Right Elevon" << std::endl; }

    aerodynamicCoefficientSettings->setControlSurfaceSettings(
                tudat::simulation_setup::readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                    forceCoefficientFiles_CS_ER,
                    momentCoefficientFiles_CS_ER,
                    controlSurfaceIndependentVariableNames), ELEVON_R );

    ///////// End: Vehicle Aerodynamics Section
    if( debugInfo == 1 ){ std::cout << "           End: Vehicle Aerodynamics Section" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "Create vehicle object" << std::endl; }
    //! Create vehicle objects.
    bodyMap[ vehicleName ] = std::make_shared< simulation_setup::Body >( );

    //! Create vehicle systems and initialize by setting the landing (dry) Mass.
    std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared< system_models::VehicleSystems >( landingMass );
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = std::make_shared< bislip::BislipVehicleSystems >( landingMass );

    //! Pass shared pointer of vehicle and bislip systems.
    bodyMap.at( vehicleName )->setVehicleSystems( vehicleSystems );
    bodyMap.at( vehicleName )->setBislipSystems( bislipSystems );

    //vehicleSystems->setInitialMass( initialMass_Ascent );

    //! Set nose radius.
    vehicleSystems->setNoseRadius( noseRadius );

    //! Set wall emmisivity.
    vehicleSystems->setWallEmissivity( epsilon );

    bislipSystems->setDebugInfo( debugInfo );

    //! Set wing sweep angle.
    bislipSystems->setWingSweepAngle( lambda );

    //! Set nose radius.
    bislipSystems->setNoseRadius( noseRadius );

    //! Set leading edge radius.
    bislipSystems->setLeadingEdgeRadius( leadingEdgeRadius );

    //! Set local body angle.
    bislipSystems->setLocalBodyAngle( phi );

    //! Set transition length.
    bislipSystems->setTransitionDistance( x_T );

    bislipSystems->setInitialLat( initialLat_rad );
    bislipSystems->setInitialLon( initialLon_rad );
    bislipSystems->setInitialCoordinates( std::make_pair( bislipSystems->getInitialLat( ), bislipSystems->getInitialLon( ) ) );
    bislipSystems->setTargetLat( targetLat_rad );
    bislipSystems->setTargetLon( targetLon_rad );
    bislipSystems->setTargetCoordinates( std::make_pair( bislipSystems->getTargetLat( ), bislipSystems->getTargetLon( ) ) );
    bislipSystems->setInitialDistanceToTarget( initialDistanceToTarget_rad );
    bislipSystems->setFinalDistanceToTarget( finalDistanceToTarget_rad );
    bislipSystems->setStartingEpoch( simulationStartEpoch );
    //bislipSystems->setSamplingRatio( samplingRatio );
    bislipSystems->setReferenceArea( S_ref );
    bislipSystems->setReferenceValues( referenceValues );
    bislipSystems->setMassReferenceCenter( R_com );
    bislipSystems->setMomentReferenceCenter( R_mrc );
    bislipSystems->setThrustReferenceCenter( R_cot );
    bislipSystems->setRelativeChestForwardAngle( relativeChestForwardAngle_rad );
    bislipSystems->setVertebralColumnInclinationAngle( vertebralColumnInclinationAngle_rad );
    bislipSystems->setBodyFrameToPassengerFrameTransformationMatrix( bodyFrameToPassengerFrameTransformationMatrix );

    bislipSystems->setMinimumDynamicPressureforControlSurface( minimumDynamicPressureforControlSurface );
    bislipSystems->setAverageEarthRadius( radiusEarth );
    bislipSystems->setBodyFlapDeflectionLimits( std::make_pair( -20, 30 ) );
    bislipSystems->setElevonDeflectionLimits( std::make_pair( -40, 40 ) );

    bislipSystems->setBankAngleReversalTimepoint( simulationStartEpoch );
    bislipSystems->setValidationFlag( false );





    //! Set vehicle aerodynamic coefficients.
    bodyMap.at( vehicleName )->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, vehicleName ) );

    //std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface >  > controlSurfaceList;
    //controlSurfaceList[ ELEVON_L ] = controlSurfaceInterface;

    //! Set vehicle system models for aerodynamic coefficient increments.
    //bodyMap[ vehicleName ]->setVehicleSystems( systemsModels );
    //bodyMap[ vehicleName ]->getAerodynamicCoefficientInterface( )->setControlSurfaceIncrements( controlSurfaceList );

    //! Finalize body creation. Not entirely sure what this does. I believe it
    //! may place the body (Earth) at the barycenter. Probably another hack to
    //! facilitate bringing in everything into the inertial frame.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    bodyMap.at( vehicleName )->setEphemeris(
                std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d  > >( ), centralBodyName ) );
    /*
     *Thr problem with the mass propagation can be solved (-ish), by turning the throw.... of this error into a std::cerr (printing the error, but not terminating the program). Alternatively, it should also work if you use:
     *
     *   bodyMap[ "Spaceplane" ]->setEphemeris(
     *               std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
     *                   std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d  > >( ),   "Earth" ) );
     *
     *What I think is going on is the following:
     * - When setting up the propagation, the code checks if all required environment models are available.
     * - The mass propagation requires the current state of the spacecraft
     * - The current state of the spacecraft requires that either the spacecraft has an ephemeris, or this state is propagated.
     * - When propagating the mass and state, the state is propagated. However, when making this check for the mass propagation, the code does not know this, resulting in this error.
     *
     *Either one of the fixes I give above should work. Let me know if it works out,
     */


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE KNOWN STATES              ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Creating known states" << std::endl;
    if( debugInfo == 1 ){ std::cout << "Creating (partial) known states" << std::endl; }

    //! Impose constraints on first and last energy nodesAscent
    //const double mu = spice_interface::getBodyGravitationalParameter( centralBodyName );
    //a = 301.7;//NRLMSISE00Atmosphere::getSpeedOfSound( R_E + height( 0 ), 0, 0, simulationStartEpoch );
    //double V_i = a * Mach_i;
    //double V_DN = 0.99 * sqrt( mu / ( R_E + maximumHeightAllowable ) );

    //! Set spherical elements for vehicle's initial state. Arbitrarily chosen to be in Earth-Fixed frame.
    //!     Transformation to Inertial Frame is done within fitness function, as the initial velocity is a
    //!     decision variable.
    Eigen::Vector6d initialState_spherical;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::radiusIndex )       = radiusEarth + initialHeight;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::latitudeIndex )     = initialLat_rad;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::longitudeIndex )    = initialLon_rad;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::speedIndex )        = 0.0;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::flightPathIndex )   = initialFlightPathAngle_rad;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::headingAngleIndex ) = 0.0;

    problemInput->setInitialState_Spherical( initialState_spherical );

    //! Create Earth's rotational ephemeris.
    std::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodyMap.at( centralBodyName )->getRotationalEphemeris( );

    problemInput->setEarthRotationalEphemeris( earthRotationalEphemeris );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Setting Accelerations" << std::endl; }

    //! Declare acceleration data map.
    SelectedAccelerationMap accelerationSettingsMap;

    //! Define gravitational model.
    //!     The central body acts this force on the vehicle.
    //!     Arbitrary maximum degree/order. Equivalent functionality to Cartesian with corresponding maximum degree/order.
    accelerationSettingsMap[ vehicleName ][ centralBodyName ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    //! Define aerodynamic accelerations.
    //!     The atmosphere of the central body acts this force on the vehicle.
    accelerationSettingsMap[ vehicleName ][ centralBodyName ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );

    //! Declare and initialize thrust guidance pointer.
    std::shared_ptr< bislip::MyGuidance > ThrustGuidance = std::make_shared< bislip::MyGuidance >(
                bodyMap,
                vehicleName,
                centralBodyName,
                bislipSystems->getStartingEpoch() );

    //! Declare and initialize thrust guidance direction guidance settings.
    std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            std::make_shared< ThrustDirectionGuidanceSettings >( thrust_direction_from_existing_body_orientation, centralBodyName );

    //! Declare and initialize thrust magnitude function.
    std::function< double( const double ) > thrustMagnitudeFunction =
            std::bind( &bislip::MyGuidance::getCurrentThrustMagnitude, ThrustGuidance );

    //! Declare and initialize thrust specific impulse function.
    std::function< double( const double ) > specificImpulseFunction =
            std::bind( &bislip::MyGuidance::getCurrentSpecificImpulse, ThrustGuidance );

    //! Declare and initialize engine status function.
    std::function< bool( const double ) > isEngineOnFunction =
            std::bind( &bislip::MyGuidance::getCurrentEngineStatus, ThrustGuidance );

    //! Declare and initialize thrust direction vector function.
    std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection =
            std::bind( &bislip::MyGuidance::getCurrentBodyFixedThrustDirection, ThrustGuidance );

    //! Declare and initialize thrust guidance update function.
    std::function< void( const double ) > customThrustResetFunction =
            std::bind( &bislip::MyGuidance::updateGuidance, ThrustGuidance, std::placeholders::_1 );

    //! Declare and initialize thrust magnitude settings.
    std::shared_ptr< ThrustMagnitudeSettings > thrustMagnitudeSettings =
            std::make_shared< FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction,
                specificImpulseFunction,
                isEngineOnFunction,
                bodyFixedThrustDirection,
                customThrustResetFunction );

    //! Pass engine capabilities to vehicle systems.
    bislipSystems->setMaxThrust( maxThrust );
    bislipSystems->setSpecificImpulse( specificImpulse );

    //! Define thrust acceleration settings.
    //!     The vehicle acts this force on itself.
    accelerationSettingsMap[ vehicleName ][ vehicleName ].push_back( std::make_shared< ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );

    //std::cout << "Accelerations Set" << std::endl;
    if( debugInfo == 1 ){ std::cout << "Accelerations Set" << std::endl; }

    //! Declare vectors that will contain the names of relevant bodies.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;

    //! Define bodies that will be propagated. Only 1.
    bodiesToIntegrate.push_back( vehicleName );
    problemInput->setBodiesToIntegrate( bodiesToIntegrate );

    //! Define central bodies. Only 1, Earth.
    centralBodies.push_back( centralBodyName );
    problemInput->setCentralBodies( centralBodies );

    //! Set acceleration models.
    basic_astrodynamics::AccelerationMap accelerationMap =
            createAccelerationModelsMap(
                bodyMap,
                accelerationSettingsMap,
                bodiesToIntegrate,
                centralBodies );

    problemInput->setAccelerationMap( accelerationMap );

    //! **************************************************************************************
    //! ********** Aerodynamic guidance is set AFTER the accelerations and BEFORE propagating.
    //! **************************************************************************************

    if( debugInfo == 1 ){ std::cout << "Declare and assign aerodynamic guidance functions" << std::endl; }

    //! Declare and assign aerodynamic guidance functions.
    std::shared_ptr< bislip::MyGuidance > AeroGuidance = std::make_shared< bislip::MyGuidance >(
                bodyMap,
                vehicleName,
                centralBodyName,
                bislipSystems->getStartingEpoch() );

    //! Set Guidance angle functions.
    setGuidanceAnglesFunctions( AeroGuidance, bodyMap.at( vehicleName ) );

    //! **************************************************************************************
    //! ********** Create alpha-Mach interpolators delineating aerodyanmic data range.
    //! **************************************************************************************

    if( debugInfo == 1 ){ std::cout << "     Create alpha-Mach interpolators delineating aerodynamic data range" << std::endl; }

    bislip::Variables::createAlphaMachBoundingInterpolators( bodyMap, vehicleName, alphaMachEnvelopeUB, alphaMachEnvelopeLB, outputPath, outputSubFolder );

    //! **************************************************************************************
    //! ********** Create heading error deadband interpolators used for bank angle reversals.
    //! **************************************************************************************

    if( debugInfo == 1 ){ std::cout << "     Create heading error deadband interpolators used for bank angle reversals" << std::endl; }

    bislip::Variables::createHeadingErrorDeadBandInterpolator( bodyMap, vehicleName, headingErrorDeadBand, outputPath, outputSubFolder);

    //! **************************************************************************************
    //! ********** Create aerodynamic guidance interpolators used for HORUS Validation.
    //! **************************************************************************************

    if ( ( problemName.find( "HORUS" ) != std::string::npos ) && ( problemName.find( "Validation" ) != std::string::npos ) && ( problemName.find( "Kourou" ) != std::string::npos ) )
    {
        bislipSystems->setValidationFlag( true );
        bislipSystems->setMaxThrust( 0.0 );

        bislipSystems->setMaxThrust( 0.0 );




        if( debugInfo == 1 ){ std::cout << "     Create HORUS Re-entry to Kourou guidance interpolators" << std::endl; }

        const std::vector< double > kourouAngleofAttackHistory = bislip::Variables::getDataNumeri ( primary[ 20 ] );
        const std::vector< double > kourouBankAngleHistory     = bislip::Variables::getDataNumeri ( primary[ 21 ] );

        bislip::Variables::createKourouGuidanceInterpolators( bodyMap, vehicleName, kourouAngleofAttackHistory, kourouBankAngleHistory, outputPath, outputSubFolder );
    }

    if( debugInfo == 1 ){ std::cout << "Vehicle guidance is set" << std::endl; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          CREATE LIST OF DEPENDENT VARIABLES        ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Creating list of dependent variables" << std::endl; }

    //! Create vector that will contian the list of dependent variables to save/output.
    //!     The file that prints out the text saying what is saved has been modified
    //!     to mute this message. During optimization this would clutter the
    //!     terminal screen and is incredibly annoying.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dep_varList;

    //! Add dependent variables to list.
    //! 0-2
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_cartesian_position,
                    vehicleName,
                    centralBodyName ) );
    //! 3-5
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_spherical_position,
                    vehicleName,
                    centralBodyName ) );
    //! 6
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::heading_angle,
                    centralBodyName) );
    //! 7
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::flight_path_angle,
                    centralBodyName) );
    //! 8
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::angle_of_attack,
                    centralBodyName) );
    //! 9
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::angle_of_sideslip,
                    centralBodyName) );
    //! 10
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::bank_angle,
                    centralBodyName) );
    //! 11
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicleName,
                    centralBodyName ) );
    //! 12
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    mach_number_dependent_variable,
                    vehicleName ) );
    //! 13
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    airspeed_dependent_variable,
                    vehicleName ) );
    //! 14
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_aerodynamic_g_load_variable,
                    vehicleName,
                    centralBodyName) );
    //! 15-17
    dep_varList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic,
                    vehicleName,
                    centralBodyName,
                    false,
                    -1 )); //! false prints vector components. -1 prints all components.
    //! 18-20
    dep_varList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    spherical_harmonic_gravity,
                    vehicleName,
                    centralBodyName,
                    false,
                    -1 )); //! false prints vector components. -1 prints all components
    //! 21-23
    dep_varList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration,
                    vehicleName,
                    vehicleName,
                    false,
                    -1 )); //! false prints vector components. -1 prints all components.
    //! 24
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_dynamic_pressure_dependent_variable,
                    vehicleName,
                    centralBodyName) );
    //! 25
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_aerodynamic_heat_rate_dependent_variable,
                    vehicleName,
                    centralBodyName) );
    //! 26
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_mass_rate_dependent_variables,
                    vehicleName ) );
    //! 27
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    current_mass,
                    vehicleName ) );
    //! 28
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    specific_energy,
                    vehicleName,
                    centralBodyName ) );
    //! 29
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    normalized_specific_energy,
                    vehicleName,
                    centralBodyName ) );
    //! 30
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    evaluated_throttle_setting,
                    vehicleName,
                    centralBodyName ) );
    //! 31
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    evaluated_thrust_elevation_angle,
                    vehicleName,
                    centralBodyName ) );
    //! 32
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    evaluated_thrust_azimuth_angle,
                    vehicleName,
                    centralBodyName ) );
    //! 33
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    evaluated_angle_of_attack,
                    vehicleName,
                    centralBodyName ) );
    //! 34
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    evaluated_bank_angle,
                    vehicleName,
                    centralBodyName ) );
    //! 35
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    engine_status,
                    vehicleName,
                    centralBodyName ) );
    //! 36
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    angular_distance_traveled,
                    vehicleName,
                    centralBodyName ) );
    //! 37
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    angular_distance_to_go,
                    vehicleName,
                    centralBodyName ) );
    //! 38
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heading_to_target,
                    vehicleName,
                    centralBodyName ) );
    //! 39
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heading_error,
                    vehicleName,
                    centralBodyName ) );
    //! 40
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heat_flux_tauber,
                    vehicleName,
                    centralBodyName ) );
    //! 41-43
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_thrust_vector,
                    vehicleName,
                    centralBodyName ) );
    //! 44
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    bending_moment,
                    vehicleName,
                    centralBodyName ) );
    //! 45-46
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_gravity,
                    vehicleName,
                    centralBodyName ) );
    //! 47
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    skip_suppression_limit,
                    vehicleName,
                    centralBodyName ) );
    //! 48
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    control_surface_deflection_dependent_variable,
                    vehicleName,
                    BODYFLAP ) );
    //! 49
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    bodyflap_deflection_moment_coefficient_increment,
                    vehicleName,
                    centralBodyName ) );
    //! 50
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    bodyflap_deflection_moment_coefficient_increment_dif,
                    vehicleName,
                    centralBodyName ) );
    //! 51-53
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_total_load_vector,
                    vehicleName,
                    centralBodyName ) );
    //! 54-56
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_total_g_load_vector,
                    vehicleName,
                    centralBodyName ) );
    //! 57
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_total_g_load_magnitude,
                    vehicleName,
                    centralBodyName ) );
    //! 58-60
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_aero_load_vector,
                    vehicleName,
                    centralBodyName ) );
    //! 61-63
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_force_coefficients_dependent_variable,
                    vehicleName,
                    centralBodyName ) );
    //! 64-66
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_moment_coefficients_dependent_variable,
                    vehicleName,
                    centralBodyName ) );
    //! 67
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heat_flux_chapman,
                    vehicleName,
                    centralBodyName ) );
    //! 68
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_density_dependent_variable,
                    vehicleName,
                    centralBodyName ) );
    //! 69-71
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    passenger_fixed_total_g_load_vector,
                    vehicleName,
                    centralBodyName ) );
    //! 72
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    commanded_throttle_setting,
                    vehicleName,
                    centralBodyName ) );
    //! 73
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    commanded_thrust_elevation_angle,
                    vehicleName,
                    centralBodyName ) );
    //! 74
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    commanded_thrust_azimuth_angle,
                    vehicleName,
                    centralBodyName ) );
    //! 75
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    commanded_angle_of_attack,
                    vehicleName,
                    centralBodyName ) );
    //! 76
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    commanded_bank_angle,
                    vehicleName,
                    centralBodyName ) );
    //! 77
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    current_lift_magnitude,
                    vehicleName,
                    centralBodyName ) );
    //! 78
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    current_heading_error_deadband,
                    vehicleName,
                    centralBodyName ) );
    //! 79
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    temp_bank_angle,
                    vehicleName,
                    centralBodyName ) );
    //! 80
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    reversal_conditional,
                    vehicleName,
                    centralBodyName ) );
    //! 81
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    bank_reversal_trigger,
                    vehicleName,
                    centralBodyName ) );
    //! 82
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    wall_temperature_chapman,
                    vehicleName,
                    centralBodyName ) );
    //! 83
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    wall_temperature_tauber_stagnation,
                    vehicleName,
                    centralBodyName ) );
    //! 84
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    wall_temperature_tauber_flatplate,
                    vehicleName,
                    centralBodyName ) );
    //! 85
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heat_flux_tauber_stagnation,
                    vehicleName,
                    centralBodyName ) );
    //! 86
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heat_flux_tauber_flatplate,
                    vehicleName,
                    centralBodyName ) );
    //! 87
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    cumulative_angular_distance_travelled,
                    vehicleName,
                    centralBodyName ) );
    //! 88
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    groundtrack_difference,
                    vehicleName,
                    centralBodyName ) );
    //! 89
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    time_of_flight,
                    vehicleName,
                    centralBodyName ) );
    //! 90
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    flight_path_angle_rate,
                    vehicleName,
                    centralBodyName ) );
    //! 91
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    cumulative_cartesian_distance_travelled,
                    vehicleName,
                    centralBodyName ) );
    //! 92
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    control_surface_deflection_dependent_variable,
                    vehicleName,
                    ELEVON_L ) );
    //! 93
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    thrust_force_magnitude,
                    vehicleName,
                    centralBodyName ) );
    //! 94
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    speed_of_sound,
                    vehicleName,
                    centralBodyName ) );
    //! 95
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    adiabatic_wall_temperature,
                    vehicleName,
                    centralBodyName ) );



    /*dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    stagnation_point_heat_flux_dependent_variable,
                    vehicleName,
                    centralBodyName) );*/

    //! Create and initialize object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dep_varList );

    problemInput->setDependentVariablesToSave( dependentVariablesToSave );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE TERMINATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Creating Termination Settings" << std::endl; }

    //! Define CUSTOM termination settings.
    std::shared_ptr< PropagationTerminationSettings > customTermination_FlightPathAngleCombination_Ascent =
            std::make_shared< PropagationCustomTerminationSettings >(
                boost::bind( &bislip::Variables::customTermination_FlightPathAngleCombination_Ascent, bodyMap, vehicleName ) );

    std::shared_ptr< PropagationTerminationSettings > customTermination_FlightPathAngleCombination_Descent =
            std::make_shared< PropagationCustomTerminationSettings >(
                boost::bind( &bislip::Variables::customTermination_FlightPathAngleCombination_Descent, bodyMap, vehicleName ) );


    //! Define dependent variable termination settings.
    std::shared_ptr< PropagationTerminationSettings > thrustTerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration,
                    vehicleName,
                    vehicleName,
                    true ), bislipSystems->getMaxThrust() , false );

    std::shared_ptr< PropagationTerminationSettings > angularDistanceToGo_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    angular_distance_to_go,
                    vehicleName), tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getFinalDistanceToTarget() ), true );

    std::shared_ptr< PropagationTerminationSettings > netAngularDistanceTravelled_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    angular_distance_traveled,
                    vehicleName), tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getInitialDistanceToTarget() ), false );

    std::shared_ptr< PropagationTerminationSettings > mass_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    current_mass,
                    vehicleName), bislipSystems->getLandingMass(), true );

    std::shared_ptr< PropagationTerminationSettings > E_hat_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    normalized_specific_energy,
                    vehicleName), 1.0 , false );

    std::shared_ptr< PropagationTerminationSettings > q_dyn_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_dynamic_pressure_dependent_variable,
                    vehicleName,
                    centralBodyName), constraint_DynamicPressure * 1.5 , false );

    std::shared_ptr< PropagationTerminationSettings > q_dot_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_aerodynamic_heat_rate_dependent_variable,
                    vehicleName,
                    centralBodyName), constraint_HeatingRate , false );

    std::shared_ptr< PropagationTerminationSettings > maximumHeightAllowable_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicleName,
                    centralBodyName), maximumHeightAllowable , false );

    std::shared_ptr< PropagationTerminationSettings > maximumAirspeed_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    airspeed_dependent_variable,
                    vehicleName,
                    centralBodyName), maximumHeightAllowable , false );

    std::shared_ptr< PropagationTerminationSettings > initialHeight_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicleName,
                    centralBodyName), initialHeight * 9.0 / 10.0, true );

    std::shared_ptr< PropagationTerminationSettings > minimumHeightAllowable_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicleName,
                    //              centralBodyName), 0 , true );
                    centralBodyName), minimumHeightAllowable , true );



    std::shared_ptr< PropagationTerminationSettings > total_aero_g_load_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_aerodynamic_g_load_variable,
                    vehicleName,
                    centralBodyName), constraint_MechanicalLoad * 1.5 , false );

    std::shared_ptr< PropagationTerminationSettings > maximumTimeOfFlight_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    time_of_flight,
                    vehicleName,
                    centralBodyName), max_tof , false );

    std::shared_ptr< PropagationTerminationSettings > negativeFlightPathAngleRate_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    flight_path_angle_rate,
                    vehicleName,
                    centralBodyName), -0.000001 , true );


    //! Create list of terminations setting for ASCENT.
    std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettingsList_Ascent;
    terminationSettingsList_Ascent.push_back( customTermination_FlightPathAngleCombination_Ascent );
    //terminationSettingsList_Ascent.push_back( E_hat_TerminationSettings );
    //terminationSettingsList_Ascent.push_back( angularDistanceToGo_TerminationSettings );
    terminationSettingsList_Ascent.push_back( netAngularDistanceTravelled_TerminationSettings );
    terminationSettingsList_Ascent.push_back( initialHeight_TerminationSettings );
    //terminationSettingsList_Ascent.push_back( total_aero_g_load_TerminationSettings );
    // terminationSettingsList_Ascent.push_back( q_dyn_TerminationSettings );
    //terminationSettingsList.push_back( q_dot_TerminationSettings );
    //terminationSettingsList_Ascent.push_back( maximumHeightAllowable_TerminationSettings );
    //terminationSettingsList_Ascent.push_back( minimumHeightAllowable_TerminationSettings );
    //terminationSettingsList_Ascent.push_back( negativeFlightPathAngle_TerminationSettings );
    //terminationSettingsList_Ascent.push_back( negativeFlightPathAngleRate_TerminationSettings );
    terminationSettingsList_Ascent.push_back( maximumTimeOfFlight_TerminationSettings );


    //! Finalize terminations setting for ASCENT.
    std::shared_ptr< PropagationTerminationSettings > terminationSettings_Ascent = std::make_shared<
            propagators::PropagationHybridTerminationSettings >( terminationSettingsList_Ascent, true );

    problemInput->setAscentTerminationSettings( terminationSettings_Ascent );

    //! Create list of terminations setting for DESCENT.
    std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettingsList_Descent;
    // terminationSettingsList_Descent.push_back( customTermination );
    //terminationSettingsList_Descent.push_back( thrustTerminationSettings );
    terminationSettingsList_Descent.push_back( angularDistanceToGo_TerminationSettings );
    terminationSettingsList_Descent.push_back( netAngularDistanceTravelled_TerminationSettings );
    //terminationSettingsList_Descent.push_back( total_aero_g_load_TerminationSettings );
    //terminationSettingsList_Descent.push_back( mass_TerminationSettings );
    //terminationSettingsList_Descent.push_back( E_hat_TerminationSettings );
    //terminationSettingsList_Descent.push_back( q_dyn_TerminationSettings );
    //terminationSettingsList_Descent.push_back( q_dot_TerminationSettings );
    //    terminationSettingsList_Descent.push_back( height_UP_TerminationSettings );
    terminationSettingsList_Descent.push_back( minimumHeightAllowable_TerminationSettings );
    //terminationSettingsList_Descent.push_back( positiveFlightPathAngle_TerminationSettings );
    terminationSettingsList_Descent.push_back( maximumTimeOfFlight_TerminationSettings );

    //! Finalize terminations setting for DESCENT.
    std::shared_ptr< PropagationTerminationSettings > terminationSettings_Descent = std::make_shared<
            propagators::PropagationHybridTerminationSettings >( terminationSettingsList_Descent, true );

    problemInput->setDescentTerminationSettings( terminationSettings_Descent );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE MASS RATE SETTINGS: ASCENT            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Calculate initial heading angle to target: https://www.movable-type.co.uk/scripts/latlong.html
    double chi_i_deg_calc = bislip::Variables::convertNegativeAnglesToPositive( unit_conversions::convertRadiansToDegrees( bislip::Variables::computeHeadingToTarget( initialLat_rad , initialLon_rad , targetLat_rad , targetLon_rad ) ) );

    //! Calculate ground distance to cover using Haversine formula and Sperical Law of Cosines: https://www.movable-type.co.uk/scripts/latlong.html
    const double a = std::sin( (targetLat_rad - initialLat_rad) / 2) * std::sin( (targetLat_rad - initialLat_rad) / 2) + std::cos( initialLat_rad ) * std::cos( initialLon_rad ) * std::sin( (targetLon_rad - initialLon_rad) / 2) * std::sin( (targetLon_rad - initialLon_rad) / 2);
    const double d_angular = 2 * std::atan2( std::sqrt(a) , std::sqrt(1 - a) );
    const double d_spherical_law_cosines = unit_conversions::convertRadiansToDegrees( bislip::Variables::computeAngularDistance( initialLat_rad , initialLon_rad , targetLat_rad , targetLon_rad ) );
    //std::acos( std::sin(initialLat_rad) * std::sin(targetLat_rad) + std::cos(initialLat_rad) * std::cos(targetLat_rad) * std::cos(targetLon_rad-initialLon_rad) ) * spice_interface::getAverageRadius( "Earth" );

    //! Print to Terminal screen.
    std::cout << " " << std::endl;
    std::cout << problemName << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Parameters to vary during Ascent" << std::endl;
    p=0;
    for( int i = 0; i < int( ascentParameterList.size() ) ; i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.1f %-2s %-8.1f %-2s\n") % ascentParameterList[i] % "[" % ascentParameterBounds[p] % "," % ascentParameterBounds[p+1] % "]";
        p+=2;
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Parameters to vary during Descent" << std::endl;
    p=0;
    for( int i = 0; i < int( descentParameterList.size() ) ; i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.1f %-2s %-8.1f %-2s\n") % descentParameterList[i] % "[" % descentParameterBounds[p] % "," % descentParameterBounds[p+1] % "]";
        p+=2;
    }

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Vehicle Parameters" << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[0] << "      " <<  vehicleParameterValues[0] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[1] << "      " <<  vehicleParameterValues[1] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[2] << "      " <<  vehicleParameterValues[2] << std::endl;
    std::cout << std::setw(30) << "R_mrc" << "     [ " << vehicleParameterValues[3] << " , " << vehicleParameterValues[4] << " , " << vehicleParameterValues[5] << " ]" <<  std::endl;
    std::cout << std::setw(30) << "R_com" << "     [ " << vehicleParameterValues[6] << " , " << vehicleParameterValues[7] << " , " << vehicleParameterValues[8] << " ]" <<  std::endl;
    std::cout << std::setw(30) << "R_cot" << "     [ " << vehicleParameterValues[9] << " , " << vehicleParameterValues[10] << " , " << vehicleParameterValues[11] << " ]" <<  std::endl;
    std::cout << std::setw(30) << vehicleParameterList[12] << "      " <<  vehicleParameterValues[12] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[13] << "      " <<  vehicleParameterValues[13] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[14] << "      " <<  vehicleParameterValues[14] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[15] << "      " <<  vehicleParameterValues[15] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[16] << "      " <<  vehicleParameterValues[16] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[17] << "      " <<  vehicleParameterValues[17] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[18] << "      " <<  vehicleParameterValues[18] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[19] << "      " <<  vehicleParameterValues[19] << std::endl;
    std::cout << std::setw(30) << vehicleParameterList[20] << "      " <<  vehicleParameterValues[20] << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Aerodynamic Coefficients File Names" << std::endl;
    for( int i = 0; i < int( aeroCoeffFileList.size() ) ; i++)
    {
        std::cout << "File " << i + 1 << ":      " <<  aeroCoeffFileList[ i ] << std::endl;
    }
    //  std::cout << "Coefficient File: " << "      " <<  aeroCoeffFileList[1] << std::endl;
    //  std::cout << "Coefficient File: " << "      " <<  aeroCoeffFileList[2] << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Optimizations Settings" << std::endl;
    for( int i = 0; i < int( optimizationSettingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.0f\n") % optimizationSettingsList[i] % "" % optimizationSettingsValues[i];
    }

    const int index = int( optimizationSettingsValues[ 0 ] );
    std::string algo_method;
    switch( index )
    {
    case 0:
    { algo_method = "NSGA2"; break; }
    case 1:
    { algo_method = "MOEAD"; break; }
    case 2:
    { algo_method = "IHS"; break; }
    }


    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Simulation Settings" << std::endl;
    for( int i = 0; i < int( simulationSettingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.0f\n") % simulationSettingsList[i] % "" % simulationSettingsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Output Settings" << std::endl;
    for( int i = 0; i < int( outputSettingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.2f\n") % outputSettingsList[i] % "" % outputSettingsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Initial Conditions" << std::endl;
    for( int i = 0; i < int( initialConditionsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.2f\n") % initialConditionsList[i] % "" % initialConditionsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Termination Conditions & Constraints" << std::endl;
    for( int i = 0; i < int( constraintsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.2f\n") % constraintsList[i] % "" % constraintsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Ground distance to cover " << std::endl;
    std::cout << "  Haversine Formula:               " << unit_conversions::convertRadiansToDegrees( d_angular ) << " degrees." << std::endl;
    std::cout << "  Spherical Law of Cosines:        " << d_spherical_law_cosines << " degrees." << std::endl;
    std::cout << "                                   " << d_spherical_law_cosines * radiusEarth / 1E3 << " km." << std::endl;
    std::cout << "  Initial heading angle to Target: " << chi_i_deg_calc << " degrees. Calculated."  << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Output subfolder: " << std::endl;
    std::cout << "     '" << outputSubFolder <<"'" << std::endl;
    std::cout << " DV   " << "  Lower Boundary " << "  Upper Boundary " << std::endl;
    for( int i = 0; i < ( N + 6 ) + NN + 2; i++ )
    {
        std::cout << std::fixed << std::setprecision(10) <<
                     std::setw(6) << i <<
                     std::setw(17) << bounds[ 0 ][ i ] <<
                     std::setw(17) << bounds[ 1 ][ i ] << std::endl;
    }


    //*********************************************************************************************************************

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PAGMO PROBLEM TO COMPUTE FITNESS            ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bislipSystems->setPropagationStepSize( propagationStepSize );
    bislipSystems->setGuidanceStepSize( guidanceStepSize );


    if( debugInfo == 1 ){ std::cout << "Creating Pagmo Problem" << std::endl; }

    //! Create object to compute the problem fitness.
    pagmo::problem prob{ Space4ErrBodyProblem( problemInput, bodyMap ) };

    //! Retrieve algorithm. Three options available in the following function:
    //!        getMultiObjectiveAlgorithm
    //!               case 0 --> nsga2
    //!               case 1 --> moead
    //!               case 2 --> ihs
    //! Selection is currently arbitrary. moead.hpp should be modified such that
    //! the points are generated with a Halton sequence.

    pagmo::algorithm algo{ getPagmoAlgorithm( int( optimizationSettingsValues[ 0 ] ) ) };

    //! Assign population size.
    pagmo::population::size_type populationSize = int( optimizationSettingsValues[ 1 ] );

    //! Assign archipelago size.
    //pagmo::archipelago::size_type archipelagoSize = int( optimizationSettingsValues[ 2 ] );

    //! Assign population per archipelago
    //const int pop_per_archi = populationSize/archipelagoSize;

    //! Assign number of evolutions
    const int evolutions =  int( optimizationSettingsValues[ 3 ] );

    if( debugInfo == 1 ){ std::cout << "Set Seed for reproducible results" << std::endl; }

    //! Set seed for reproducible results.
    pagmo::random_device::set_seed( int( optimizationSettingsValues[ 4 ] ) );

    //! Assign number of evolutions
    pagmo::population::size_type monteCarloPopulation =  int( optimizationSettingsValues[ 5 ] );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PAGMO PROBLEM                               ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Create an island with populationSize individuals
    if ( monteCarloPopulation != 0 )
    {
        std::cout << "Running a Monte Carlo population of " << monteCarloPopulation << " individuals" << std::endl;

        bislipSystems->setPropagationStepSize( propagationStepSize * 10 );
        bislipSystems->setGuidanceStepSize( guidanceStepSize * 10 );

        pagmo::island pagmoIsland{ algo, prob, monteCarloPopulation };

        std::cout << monteCarloPopulation << " individuals completed" << std::endl;

        //! Write original (unevolved) population and fitness to file
        if  ( int( outputSettingsValues[ 0 ] ) == 1 )
        {
            //!
            std::string this_run_settings_1 = std::to_string( 0 ) + "_" + this_run_settings;

            Eigen::MatrixXd monteCarloPopulationMatrix = bislip::Variables::convertVectorOfVectorsDoubleToEigenMatrixXd( pagmoIsland.get_population( ).get_x( ) );

            Eigen::MatrixXd monteCarloFitnessMatrix = bislip::Variables::convertVectorOfVectorsDoubleToEigenMatrixXd( pagmoIsland.get_population( ).get_f( ) );

            //! Population
            bislip::Variables::printEigenMatrixXdToFile( monteCarloPopulationMatrix,
                                                         "monteCarloPopulation_" + this_run_settings_1 + ".dat",
                                                         outputSubFolder );
            //! Fitness
            bislip::Variables::printEigenMatrixXdToFile( monteCarloFitnessMatrix,
                                                         "monteCarloFitness_" + this_run_settings_1 + ".dat",
                                                         outputSubFolder );
        }
    }
    else
    {
        pagmo::island pagmoIsland{ algo, prob, populationSize };

        std::cout << " " << std::endl;
        std::cout << "PAGMO island created" << std::endl;
        std::cout << " " << std::endl;


        Eigen::MatrixXd finalPopulationMatrix( pagmoIsland.get_population( ).get_x( ).size( ), pagmoIsland.get_population( ).get_x( ).at( 0 ).size( ) );
        Eigen::MatrixXd finalFitnessMatrix( pagmoIsland.get_population( ).get_f( ).size( ), pagmoIsland.get_population( ).get_f( ).at( 0 ).size( ) );

        for( int k = 1; k < evolutions + 1; k++ )
        {
            pagmoIsland.evolve( );
            while( pagmoIsland.status( ) != pagmo::evolve_status::idle &&
                   pagmoIsland.status( ) != pagmo::evolve_status::idle_error )
            {
                pagmoIsland.wait( );
            }
            pagmoIsland.wait_check( ); // Raises errors

            if  ( int( outputSettingsValues[ 0 ] ) == 1 )
            {
                //if ( k == evolutions )
                //{
                //!
                std::string this_run_settings_1 = std::to_string( k ) + "_" + this_run_settings;

                //! Write current iteration results to file

                finalPopulationMatrix = bislip::Variables::convertVectorOfVectorsDoubleToEigenMatrixXd( pagmoIsland.get_population( ).get_x( ) );

                finalFitnessMatrix = bislip::Variables::convertVectorOfVectorsDoubleToEigenMatrixXd( pagmoIsland.get_population( ).get_f( ) );

                //! Population
                bislip::Variables::printEigenMatrixXdToFile( finalPopulationMatrix,
                                                             "population_" + this_run_settings_1 + ".dat",
                                                             outputSubFolder );
                //! Fitness
                bislip::Variables::printEigenMatrixXdToFile( finalFitnessMatrix,
                                                             "fitness_" + this_run_settings_1 + ".dat",
                                                             outputSubFolder );

                ///}
            }

            std::cout<< "  " << std::endl;
            std::cout<< k    << std::endl;
            std::cout<< "  " << std::endl;

        }


    }




    //! Instantiate an archipelago
    /*  archipelago archi{archipelagoSize, algo, prob, pop_per_archi};

    //! Evolve
    archi.evolve(evolutions);

    //! Wait for the evolutions to be finished
    archi.wait_check();// Raises errors,

    for (const auto &isl : archi) {
        std::cout << pagmoIsland.get_population().champion_f()[0] << '\n';
    }
    for (const auto &isl : archi) {
        std::cout << pagmoIsland.get_population().champion_x()[0] << '\n';
    }
    for (const auto &isl : archi) {
     printPopulationToFile( pagmoIsland.get_population().champion_f()[0],
             "mo_AMSIAD_BALLISTIC_" ,
             false );
     printPopulationToFile( pagmoIsland.get_population().champion_x()[0],
             "mo_AMSIAD_BALLISTIC_" ,
             true );
    }
*/
















}

