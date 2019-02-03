#include "Space4ErrBody_Headers/Space4ErrBodyProblem.h"
#include "Space4ErrBody_Headers/getPagmoAlgorithm.h"
#include "Space4ErrBody_Headers/saveOptimizationResults.h"

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
using namespace tudat_applications;
using namespace pagmo;
using namespace bislip;
using namespace bislip::Variables;
using namespace bislip::Parameters;
using namespace boost;
using namespace tudat_pagmo_applications;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////            Main Program                  //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    //! Get this run's time stamp. Will be used to create a unique output
    //! subfolder name, where all files created by the optimizer will be stored.
    std::string play_time = bislip::Variables::getCurrentDateTime( false );

    //! Get all inputs.
    const std::vector< std::string > primary                   = bislip::Variables::getDataString ( "primary.txt" );
    const std::vector< std::string > parameterList_Ascent      = bislip::Variables::getDataString ( primary[ 0 ] );
    const std::vector< double > parameterBounds_Ascent         = bislip::Variables::getDataNumeri ( primary[ 1 ] );
    const std::vector< std::string > parameterList_Descent     = bislip::Variables::getDataString ( primary[ 2 ] );
    const std::vector< double > parameterBounds_Descent        = bislip::Variables::getDataNumeri ( primary[ 3 ] );
    const std::vector< std::string > vehicleParameterList      = bislip::Variables::getDataString ( primary[ 4 ] );
    const std::vector< double > vehicleParameterValues         = bislip::Variables::getDataNumeri ( primary[ 5 ] );
    const std::vector< std::string > aeroCoeffFileList         = bislip::Variables::getDataString ( primary[ 6 ] );
    const std::vector< std::string > optimization_settingsList = bislip::Variables::getDataString ( primary[ 7 ] );
    const std::vector< double > optimization_settingsValues    = bislip::Variables::getDataNumeri ( primary[ 8 ] );
    const std::vector< std::string > simulation_settingsList   = bislip::Variables::getDataString ( primary[ 9 ] );
    const std::vector< double > simulation_settingsValues      = bislip::Variables::getDataNumeri ( primary[ 10 ] );
    const std::vector< std::string > output_settingsList       = bislip::Variables::getDataString ( primary[ 11 ] );
    const std::vector< double > output_settingsValues          = bislip::Variables::getDataNumeri ( primary[ 12 ] );
    const std::vector< std::string > initialConditionsList     = bislip::Variables::getDataString ( primary[ 13 ] );
    const std::vector< double > initialConditionsValues        = bislip::Variables::getDataNumeri ( primary[ 14 ] );
    const std::vector< std::string > terminationConditionsList = bislip::Variables::getDataString ( primary[ 15 ] );
    const std::vector< double > terminationConditionsValues    = bislip::Variables::getDataNumeri ( primary[ 16 ] );
    const std::string problem_name                             = primary.rbegin()[ 1 ];
    const std::string vehicleName                              = primary.rbegin()[ 0 ];

    //! Set seed for reproducible results.
    pagmo::random_device::set_seed( int( optimization_settingsValues[ 4 ] ) );

    //! Create partial output subfolder filename suffix based on optimization
    //! settings and fixed time step.
    const std::string this_run_settings = std::to_string( int( optimization_settingsValues[ 0 ] ) ) + "_" +
            std::to_string( int( optimization_settingsValues[ 1 ] ) ) + "_" +
            std::to_string( int( optimization_settingsValues[ 2 ] ) ) + "_" +
            std::to_string( int( optimization_settingsValues[ 3 ] ) ) + "_" +
            std::to_string( int( optimization_settingsValues[ 4 ] ) ) + "_" +
            std::to_string( simulation_settingsValues[ 3 ] );

    //! Create output subfolder filename. Based on arbitrary prefix, the
    //! previously created partial suffix, and the run's time stamp.
    const std::string outputSubFolder = "OUTPUT_" + this_run_settings + "_" + play_time + "/";

    //! Determine number of Ascent parameters to vary based of size of parameterBounds_Ascent vector.
    const unsigned long nodesAscent = simulation_settingsValues[ 3 ];
    const int N = ( parameterList_Ascent.size() - 4 ) * nodesAscent - 1;
    const int M = parameterBounds_Ascent.size();

    //! Determine number of Descent parameters to vary based of size of parameterBounds_Descent vector.
    const unsigned long nodesDescent = simulation_settingsValues[ 4 ];
    const int NN = ( parameterList_Descent.size() - 1 ) * nodesDescent - 1;
    const int MM = parameterBounds_Descent.size();

    //! Create vector containing decision vector bounds.
    std::vector< std::vector< double > > bounds( 2, std::vector< double >( N + 4 + NN + 1, 0.0 ) );
    //!---------------------------------------   ^ for lb/up  (rows)       ^ for # of parameters

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

        bounds[ 0 ][ i ] = parameterBounds_Ascent[ p ];
        bounds[ 1 ][ i ] = parameterBounds_Ascent[ p + 1 ];
        std::cout << "bounds[ 0 ][ " << i << " ] = " << parameterBounds_Ascent[ p ] << std::endl;
        std::cout << "bounds[ 1 ][ " << i << " ] = " << parameterBounds_Ascent[ p + 1 ] << std::endl;


    }

    bounds[ 0 ][ N ] = parameterBounds_Ascent[ M - 8 ];
    bounds[ 1 ][ N ] = parameterBounds_Ascent[ M - 7 ];
    bounds[ 0 ][ N + 1 ] = parameterBounds_Ascent[ M - 6 ];
    bounds[ 1 ][ N + 1 ] = parameterBounds_Ascent[ M - 5 ];
    bounds[ 0 ][ N + 2 ] = parameterBounds_Ascent[ M - 4 ];
    bounds[ 1 ][ N + 2 ] = parameterBounds_Ascent[ M - 3 ];
    bounds[ 0 ][ N + 3 ] = parameterBounds_Ascent[ M - 2 ];
    bounds[ 1 ][ N + 3 ] = parameterBounds_Ascent[ M - 1 ];


    std::cout << "bounds[ 0 ][ " << N << " ] = " << parameterBounds_Ascent[ M - 8 ] << std::endl;
    std::cout << "bounds[ 1 ][ " << N << " ] = " << parameterBounds_Ascent[ M - 7 ] << std::endl;
    std::cout << "bounds[ 0 ][ " << N + 1 << " ] = " << parameterBounds_Ascent[ M - 6 ] << std::endl;
    std::cout << "bounds[ 1 ][ " << N + 1 << " ] = " << parameterBounds_Ascent[ M - 5 ] << std::endl;
    std::cout << "bounds[ 0 ][ " << N + 2 << " ] = " << parameterBounds_Ascent[ M - 4 ] << std::endl;
    std::cout << "bounds[ 1 ][ " << N + 2 << " ] = " << parameterBounds_Ascent[ M - 3 ] << std::endl;
    std::cout << "bounds[ 0 ][ " << N + 3 << " ] = " << parameterBounds_Ascent[ M - 2 ] << std::endl;
    std::cout << "bounds[ 1 ][ " << N + 3 << " ] = " << parameterBounds_Ascent[ M - 1 ] << std::endl;


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

        bounds[ 0 ][ ( N + 4 ) + i ] = parameterBounds_Descent[ p ];
        bounds[ 1 ][ ( N + 4 ) + i ] = parameterBounds_Descent[ p + 1 ];
        std::cout << "bounds[ 0 ][ " << ( N + 4 ) + i << " ] = " << parameterBounds_Descent[ p ] << std::endl;
        std::cout << "bounds[ 1 ][ " << ( N + 4 ) + i << " ] = " << parameterBounds_Descent[ p + 1 ] << std::endl;
    }

    bounds[ 0 ][ N + NN + 4 ] = parameterBounds_Descent[ MM - 2 ];
    bounds[ 1 ][ N + NN + 4 ] = parameterBounds_Descent[ MM - 1 ];
    std::cout << "bounds[ 0 ][ " << ( N + 4 ) + NN << " ] = " << parameterBounds_Descent[ MM - 2 ] << std::endl;
    std::cout << "bounds[ 1 ][ " << ( N + 4 ) + NN << " ] = " << parameterBounds_Descent[ MM - 1 ] << std::endl;

    std::cout << "N = " << N << std::endl;
    std::cout << "M = " << M << std::endl;
    std::cout << "NN: " << NN << std::endl;
    std::cout << "MM: " << MM << std::endl;
    //std::cout << "parameterBounds_Descent[ MM - 2 ]: " << parameterBounds_Descent[ MM - 2 ] << std::endl;
    //std::cout << "parameterBounds_Descent[ MM - 1 ]: " << parameterBounds_Descent[ MM - 1 ] << std::endl;
    //std::cout << "N + NN + 3 + 1: " << N + NN + 3 + 1 << std::endl;

    //! Load spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            UNPACK INPUT DATA             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //std::cout << "Unpacking data" << std::endl;

    //! Declare and initialize simulation start epoch.
    const double simulationStartEpoch = simulation_settingsValues[ 0 ]; // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000

    //! Declare and initialize simulation end epoch.
    const double simulationEndEpoch = simulationStartEpoch + simulation_settingsValues[ 1 ];

    //! Declare and initialize numerical integration fixed step size.
    const double fixedStepSize = simulation_settingsValues[ 2 ];

    //! Declare and initialize number of control nodes.
    //const unsigned long nodesAscent = simulation_settingsValues[ 3 ];
    //   const unsigned long nodesDescent = simulation_settingsValues[ 4 ];

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

    //! Declare and initialize starting height
    const double h_i = initialConditionsValues[ 2 ]; // m

    //! Declare and initialize initial Mach
    const double Mach_i = initialConditionsValues[ 3 ]; // -

    //! Declare and initialize initial flight-path angle
    const double gamma_i_deg = initialConditionsValues[ 4 ]; // deg

    //! Declare and initialize starting position coordinates.
    const double initialLat_deg = initialConditionsValues[ 0 ];
    const double initialLon_deg = initialConditionsValues[ 1 ];

    //! Declare and initialize final position coordinates and additional termination conditions
    const double targetLat_deg = terminationConditionsValues[ 0 ];
    const double targetLon_deg = terminationConditionsValues[ 1 ];

    //! Declare and initialize various termination conditions
    const double finalDistanceToTarget_deg = terminationConditionsValues[ 2 ];
    //const double h_UP = terminationConditionsValues[ 3 ];
    const double h_DN = terminationConditionsValues[ 4 ];
    //const double V_UP = terminationConditionsValues[ 5 ];
    const double constraint_MechanicalLoad = terminationConditionsValues[ 6 ];
    const double constraint_HeatingRate = terminationConditionsValues[ 7 ];
    const double constraint_DynamicPressure = terminationConditionsValues[ 8 ];

    //! Still working on these
    const double R_N = 3.0;
    const double lambda = unit_conversions::convertDegreesToRadians( 30.0 );
    const double epsilon = 0.7;
    const double x_T = 1E-9;
    const double phi = unit_conversions::convertDegreesToRadians( 3.0 );

    //! Convert angles from degrees to radians
    const double initialLat_rad             = unit_conversions::convertDegreesToRadians( initialLat_deg );
    const double initialLon_rad             = unit_conversions::convertDegreesToRadians( initialLon_deg );
    const double targetLat_rad              = unit_conversions::convertDegreesToRadians( targetLat_deg );
    const double targetLon_rad              = unit_conversions::convertDegreesToRadians( targetLon_deg );
    const double finalDistanceToTarget_rad  = unit_conversions::convertDegreesToRadians( finalDistanceToTarget_deg );
    const double initialFlightPathAngle_rad = unit_conversions::convertDegreesToRadians( gamma_i_deg );

    //! Pre-define various variables used to determine fitness.
    double targetLat_deg_calc          = initialLat_deg;
    double targetLon_deg_calc          = initialLon_deg;
    double initialDistanceToTarget_rad = bislip::Variables::computeAngularDistance( initialLat_rad, initialLon_rad, targetLat_rad, targetLon_rad );
    double initialDistanceToTarget_deg = unit_conversions::convertRadiansToDegrees( initialDistanceToTarget_rad );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE PARAMETER BOUNDS          ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // std::cout << "Create Parameter Bounds" << std::endl;

    //! Declare data maps of strings and pairs to store parameter bounds.
    //!     Strings are used for the parameter names.
    //!     Pairs are used for lower and upper bounds.
    std::map< bislip::Parameters::Optimization, std::pair < double, double > > Bounds_Ascent, Bounds_Descent;

    //! Declare and initialize counter.
    //!     Information is extracted directly from vectors created from input data, which requires the counter for the upper bound value.
    p = 0;

    //! Loop to populate the data map with pairs associated to string values related to Ascent phase.
    for( unsigned long i = 0; i < parameterList_Ascent.size(); i++)
    {
        if ( parameterList_Ascent[ i ] == "Angle of Attack" ) { Bounds_Ascent[ bislip::Parameters::Optimization::AngleOfAttack ] = std::make_pair( parameterBounds_Ascent[ p ], parameterBounds_Ascent[ p + 1 ] ); }
        if ( parameterList_Ascent[ i ] == "Bank Angle" ) { Bounds_Ascent[ bislip::Parameters::Optimization::BankAngle ] = std::make_pair( parameterBounds_Ascent[ p ], parameterBounds_Ascent[ p + 1 ] ); }
        if ( parameterList_Ascent[ i ] == "Thrust Elevation Angle" ) { Bounds_Ascent[ bislip::Parameters::Optimization::ThrustElevationAngle ] = std::make_pair( parameterBounds_Ascent[ p ], parameterBounds_Ascent[ p + 1 ] ); }
        if ( parameterList_Ascent[ i ] == "Thrust Azimuth Angle" ) { Bounds_Ascent[ bislip::Parameters::Optimization::ThrustAzimuthAngle ] = std::make_pair( parameterBounds_Ascent[ p ], parameterBounds_Ascent[ p + 1 ] ); }
        if ( parameterList_Ascent[ i ] == "Throttle Setting" ) { Bounds_Ascent[ bislip::Parameters::Optimization::ThrottleSetting ] = std::make_pair( parameterBounds_Ascent[ p ], parameterBounds_Ascent[ p + 1 ] ); }
        if ( parameterList_Ascent[ i ] == "Node Interval" ) { Bounds_Ascent[ bislip::Parameters::Optimization::NodeInterval ] = std::make_pair( parameterBounds_Ascent[ p ], parameterBounds_Ascent[ p + 1 ] ); }
        if ( parameterList_Ascent[ i ] == "Initial Velocity" ) { Bounds_Ascent[ bislip::Parameters::Optimization::InitialVelocity ] = std::make_pair( parameterBounds_Ascent[ p ], parameterBounds_Ascent[ p + 1 ] ); }
        if ( parameterList_Ascent[ i ] == "Maximum Velocity" ) { Bounds_Ascent[ bislip::Parameters::Optimization::MaximumVelocity ] = std::make_pair( parameterBounds_Ascent[ p ], parameterBounds_Ascent[ p + 1 ] ); }
        if ( parameterList_Ascent[ i ] == "Maximum Height" ) { Bounds_Ascent[ bislip::Parameters::Optimization::MaximumHeight ] = std::make_pair( parameterBounds_Ascent[ p ], parameterBounds_Ascent[ p + 1 ] ); }
        if ( parameterList_Ascent[ i ] == "Additional Mass" ) { Bounds_Ascent[ bislip::Parameters::Optimization::AdditionalMass ] = std::make_pair( parameterBounds_Ascent[ p ], parameterBounds_Ascent[ p + 1 ] ); }
        p += 2;
    }

    //! Re-initialize counter.
    p = 0;

    //! Loop to populate the data map with pairs associated to string values related to Descent phase.
    for( unsigned long i = 0; i < parameterList_Descent.size(); i++)
    {
        if ( parameterList_Descent[ i ] == "Angle of Attack" ) { Bounds_Descent[ bislip::Parameters::Optimization::AngleOfAttack ] = std::make_pair( parameterBounds_Descent[ p ], parameterBounds_Descent[ p + 1 ] ); }
        if ( parameterList_Descent[ i ] == "Bank Angle" ) { Bounds_Descent[ bislip::Parameters::Optimization::BankAngle ] = std::make_pair( parameterBounds_Descent[ p ], parameterBounds_Descent[ p + 1 ] ); }
        if ( parameterList_Descent[ i ] == "Thrust Elevation Angle" ) { Bounds_Descent[ bislip::Parameters::Optimization::ThrustElevationAngle ] = std::make_pair( parameterBounds_Descent[ p ], parameterBounds_Descent[ p + 1 ] ); }
        if ( parameterList_Descent[ i ] == "Thrust Azimuth Angle" ) { Bounds_Descent[ bislip::Parameters::Optimization::ThrustAzimuthAngle ] = std::make_pair( parameterBounds_Descent[ p ], parameterBounds_Descent[ p + 1 ] ); }
        if ( parameterList_Descent[ i ] == "Throttle Setting" ) { Bounds_Descent[ bislip::Parameters::Optimization::ThrottleSetting ] = std::make_pair( parameterBounds_Descent[ p ], parameterBounds_Descent[ p + 1 ] ); }
        if ( parameterList_Descent[ i ] == "Node Interval" ) { Bounds_Descent[ bislip::Parameters::Optimization::NodeInterval ] = std::make_pair( parameterBounds_Descent[ p ], parameterBounds_Descent[ p + 1 ] ); }
        if ( parameterList_Descent[ i ] == "Initial Velocity" ) { Bounds_Descent[ bislip::Parameters::Optimization::InitialVelocity ] = std::make_pair( parameterBounds_Descent[ p ], parameterBounds_Descent[ p + 1 ] ); }
        if ( parameterList_Descent[ i ] == "Maximum Velocity" ) { Bounds_Descent[ bislip::Parameters::Optimization::MaximumVelocity ] = std::make_pair( parameterBounds_Descent[ p ], parameterBounds_Descent[ p + 1 ] ); }
        if ( parameterList_Descent[ i ] == "Maximum Height" ) { Bounds_Descent[ bislip::Parameters::Optimization::MaximumHeight ] = std::make_pair( parameterBounds_Descent[ p ], parameterBounds_Descent[ p + 1 ] ); }
        if ( parameterList_Descent[ i ] == "Additional Mass" ) { Bounds_Descent[ bislip::Parameters::Optimization::AdditionalMass ] = std::make_pair( parameterBounds_Descent[ p ], parameterBounds_Descent[ p + 1 ] ); }
        p += 2;
    }



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //std::cout << "Creating environment" << std::endl;

    //! Declare and initialize central body name.
    const std::string centralBodyName = "Earth";

    //! Declare and initialize simulation body settings data map.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { centralBodyName },
                                    simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );

    //  std::cout << "Define atmospheric model." << std::endl;
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
    const double radius_Earth = spice_interface::getAverageRadius( centralBodyName );
    double radius_Earth_i = radius_Earth;
    double radius_Earth_UP = radius_Earth;

    //!---------------------------- Still playing around with this section

    //! https://en.wikipedia.org/wiki/Flattening#Numerical_values_for_planets
    //! https://www.vcalc.com/wiki/vCalc/WGS-84+Earth+flattening+factor
    double flattening_Earth = 1 / 1;//nput_data_.back();
    if ( flattening_Earth != 1 )
    {
        bodySettings[ centralBodyName ]->shapeModelSettings = std::make_shared< OblateSphericalBodyShapeSettings >( radius_Earth, flattening_Earth );

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

        const double radius_Earth_i_x = radii[ 0 ] * std::cos(polar_angle_i_rad) * std::cos(lon_i_rad);
        const double radius_Earth_i_y = radii[ 1 ] * std::cos(polar_angle_i_rad) * std::sin(lon_i_rad);
        const double radius_Earth_i_z = radii[ 2 ] * std::sin(polar_angle_i_rad);
        radius_Earth_i = std::sqrt(pow(radius_Earth_i_x,2) + pow(radius_Earth_i_y,2) + pow(radius_Earth_i_z,2));
        const double radius_Earth_UP_x = radii[ 0 ] * std::cos(polar_angle_f_rad) * std::cos(lon_f_rad);
        const double radius_Earth_UP_y = radii[ 1 ] * std::cos(polar_angle_f_rad) * std::sin(lon_f_rad);
        const double radius_Earth_UP_z = radii[ 2 ] * std::sin(polar_angle_f_rad);
        radius_Earth_UP = std::sqrt(pow(radius_Earth_UP_x,2) + pow(radius_Earth_UP_y,2) + pow(radius_Earth_UP_z,2));
*/
    }
    //!--------------------------------------------------------------

    // std::cout << "Create Earth object" << std::endl;
    //! Create Earth object
    tudat::simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    //std::cout << "Create Earth's rotational ephemeris" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Creating vehicle" << std::endl;

    //! Assign coefficient file names.
    //const std::string dragCoefficients = aeroCoeffFileList_[0];
    //const std::string liftCoefficients = aeroCoeffFileList_[1];
    //const std::string dragControlSurfaceForceCoefficients = "HORUS_CD_CS.txt";
    //const std::string liftControlSurfaceForceCoefficients = "HORUS_CL_CS.txt";
    //const std::string momentCoefficients = aeroCoeffFileList_[2];
    //const std::string controlSurfaceMomentCoefficients = "HORUS_CM_CS.txt";

    ///////// Start: Vehicle Aerodynamics Section

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


    // forceCoefficientFiles_CS_B[ 0 ] = dragControlSurfaceForceCoefficients; // Set drag coefficient increment file
    //forceCoefficientFiles_CS_B[ 2 ] = liftControlSurfaceForceCoefficients; // Set lift coefficient increment file
    //momentCoefficientFiles_CS_B[ 1 ] = controlSurfaceMomentCoefficients; // Set pitch moment coefficient increment file

    //       std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > controlSurfaceInterface =
    //              std::make_shared< CustomControlSurfaceIncrementAerodynamicInterface >(
    //                   &dummyControlIncrements,
    //                   boost::assign::list_of( angle_of_attack_dependent )( control_surface_deflection_dependent ) );
    //       std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface >  > controlSurfaceList;
    //       controlSurfaceList[ "bodyflap" ] = controlSurfaceInterface;

    //! Define reference frame in which the loaded coefficients are defined.
    //! Have to get some more background info here to properly understand it.
    bool areCoefficientsInAerodynamicFrame = true;
    bool areCoefficientsInNegativeAxisDirection = true;


    //! Load and parse coefficient files; create coefficient settings.
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            readTabulatedAerodynamicCoefficientsFromFiles(
                forceCoefficientFiles,
                momentCoefficientFiles,
                b_ref,
                S_ref,
                c_ref,
                R_mrc,
                independentVariableNames,
                areCoefficientsInAerodynamicFrame,
                areCoefficientsInNegativeAxisDirection );

    // Add settings for control surface increments to main aerodynamic coefficients
    aerodynamicCoefficientSettings->setControlSurfaceSettings(
                readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                    forceCoefficientFiles_CS_B,
                    momentCoefficientFiles_CS_B,
                    controlSurfaceIndependentVariableNames), BODYFLAP );
    /*
     aerodynamicCoefficientSettings->setControlSurfaceSettings(
                 readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                     forceCoefficientFiles_CS_EL,
                     momentCoefficientFiles_CS_EL,
                     controlSurfaceIndependentVariableNames), ELEVON_L );
     aerodynamicCoefficientSettings->setControlSurfaceSettings(
                 readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
                     forceCoefficientFiles_CS_ER,
                     momentCoefficientFiles_CS_ER,
                     controlSurfaceIndependentVariableNames), ELEVON_R );
*/
    //bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
    //            boost::lambda::constant( constantAngleOfAttack ) );

    ///////// End: Vehicle Aerodynamics Section

    //! Create vehicle objects.
    bodyMap[ vehicleName ] = std::make_shared< simulation_setup::Body >( );

    //! Set body Mass.
    //bodyMap.at( vehicleName )->setConstantBodyMass( initialMass_Ascent );

    //! Create vehicle systems and initialize by setting the landing (dry) Mass.
    std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared< system_models::VehicleSystems >( landingMass );
    std::shared_ptr< bislip::VehicleSystems > bislipSystems = std::make_shared< bislip::VehicleSystems >( landingMass );

    //!
    bodyMap.at( vehicleName )->setVehicleSystems( vehicleSystems );
    bodyMap.at( vehicleName )->setBislipSystems( bislipSystems );

    //vehicleSystems->setInitialMass( initialMass_Ascent );

    //! Set nose radius.
    vehicleSystems->setNoseRadius( R_N );

    //! Set wing sweep angle.
    bislipSystems->setWingSweepAngle( lambda );

    //! Set wall emmisivity.
    vehicleSystems->setWallEmissivity( epsilon );

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
    bislipSystems->setReferenceArea( S_ref );
    bislipSystems->setReferenceValues( referenceValues );
    bislipSystems->setMassReferenceCenter( R_com );
    bislipSystems->setMomentReferenceCenter( R_mrc);

    vehicleSystems->setCurrentControlSurfaceDeflection( BODYFLAP, 0 );


    //! Set vehicle aerodynamic coefficients.
    bodyMap.at( vehicleName )->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface(
                    aerodynamicCoefficientSettings,
                    vehicleName ) );

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

    //! Impose constraints on first and last energy nodesAscent
    //const double mu = spice_interface::getBodyGravitationalParameter( centralBodyName );
    //a = 301.7;//NRLMSISE00Atmosphere::getSpeedOfSound( R_E + height( 0 ), 0, 0, simulationStartEpoch );
    //double V_i = a * Mach_i;
    //double V_DN = 0.99 * sqrt( mu / ( R_E + h_UP ) );

    //! Set spherical elements for vehicle's initial state. Arbitrarily chosen to be in Earth-Fixed frame.
    //!     Transformation to Inertial Frame is done within fitness function, as the initial velocity is a
    //!     decision variable.
    Eigen::Vector6d initialState_spherical;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::radiusIndex )       = radius_Earth + h_i;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::latitudeIndex )     = initialLat_rad;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::longitudeIndex )    = initialLon_rad;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::speedIndex )        = 0.0;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::flightPathIndex )   = initialFlightPathAngle_rad;
    initialState_spherical( tudat::orbital_element_conversions::SphericalOrbitalStateElementIndices::headingAngleIndex ) = bislip::Variables::computeHeadingToTarget( initialLat_rad,initialLon_rad,targetLat_rad,targetLon_rad );

    //! Create Earth's rotational ephemeris.
    std::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodyMap.at( centralBodyName )->getRotationalEphemeris( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Setting Accelerations" << std::endl;

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
                centralBodyName);

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
    //vehicleSystems->setMaxThrust( maxThrust );
    //vehicleSystems->setSpecificImpulse( specificImpulse );
    bislipSystems->setMaxThrust( maxThrust );
    bislipSystems->setSpecificImpulse( specificImpulse );

    //! Define thrust acceleration settings.
    //!     The vehicle acts this force on itself.
    accelerationSettingsMap[ vehicleName ][ vehicleName ].push_back( std::make_shared< ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );

    //std::cout << "Accelerations Set" << std::endl;

    //! Declare vectors that will contain the names of relevant bodies.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;

    //! Define bodies that will be propagated. Only 1.
    bodiesToIntegrate.push_back( vehicleName );

    //! Define central bodies. Only 1, Earth.
    centralBodies.push_back( centralBodyName );

    //! Set acceleration models.
    basic_astrodynamics::AccelerationMap accelerationsMap =
            createAccelerationModelsMap(
                bodyMap,
                accelerationSettingsMap,
                bodiesToIntegrate,
                centralBodies );

    //! **************************************************************************************
    //! ********** Aerodynamic guidance is set AFTER the accelerations and BEFORE propagating.
    //! **************************************************************************************

    //! Declare and assign aerodynamic guidance functions.
    std::shared_ptr< bislip::MyGuidance > AeroGuidance = std::make_shared< bislip::MyGuidance >(
                bodyMap,
                vehicleName,
                centralBodyName );

    //! Set Guidance angle functions.
    setGuidanceAnglesFunctions( AeroGuidance, bodyMap.at( vehicleName ) );

    //! Define constant orientation
    //double constantAngleOfAttack = unit_conversions::convertDegreesToRadians( 10 );
    //double constantBankAngle = unit_conversions::convertDegreesToRadians( 0 );
    //double constantSideSlipeAn    gle = unit_conversions::convertDegreesToRadians( 0 );
    //bodyMap.at( vehicleName )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
    //            [ = ]( ){ return constantAngleOfAttack; }, [ = ]( ){ return constantSideSlipeAngle; }, [ = ]( ){ return constantBankAngle; } );
    //std::cout << "Creating vehicle: Guidance is set" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          CREATE LIST OF DEPENDENT VARIABLES        ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Creating list of dependent variables" << std::endl;

    //! Create vector that will contian the list of dependent variables to save/output.
    //!     The file that prints out the text saying what is saved has been modified
    //!     to mute this message. During optimization this would clutter the
    //!     terminal screen and is incredibly annoying.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dep_varList;

    //! Add dependent variables to list.
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_cartesian_position,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_relative_spherical_position,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::heading_angle,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::flight_path_angle,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::angle_of_attack,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::angle_of_sideslip,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::bank_angle,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    mach_number_dependent_variable,
                    vehicleName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    airspeed_dependent_variable,
                    vehicleName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_aerodynamic_g_load_variable,
                    vehicleName,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic,
                    vehicleName,
                    centralBodyName,
                    false,
                    -1 )); //! false prints vector components. -1 prints all components.
    dep_varList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    spherical_harmonic_gravity,
                    vehicleName,
                    centralBodyName,
                    false,
                    -1 )); //! false prints vector components. -1 prints all components
    dep_varList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration,
                    vehicleName,
                    vehicleName,
                    false,
                    -1 )); //! false prints vector components. -1 prints all components.
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_dynamic_pressure_dependent_variable,
                    vehicleName,
                    centralBodyName) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_aerodynamic_heat_rate_dependent_variable,
                    vehicleName,
                    centralBodyName) );
    /*dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    stagnation_point_heat_flux_dependent_variable,
                    vehicleName,
                    centralBodyName) );*/
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_mass_rate_dependent_variables,
                    vehicleName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    current_mass,
                    vehicleName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    specific_energy,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    normalized_specific_energy,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    evaluated_throttle_setting,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    evaluated_thrust_elevation_angle,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    evaluated_thrust_azimuth_angle,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    evaluated_angle_of_attack,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    evaluated_bank_angle,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    engine_status,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    angular_distance_traveled,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    angular_distance_to_go,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heading_to_target,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heading_error,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    heat_rate_leading_edge,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_thrust_vector,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    bending_moment,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_gravity,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    equilibrium_glide_limit,
                    vehicleName,
                    centralBodyName ) );
    dep_varList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    increment_Cm_bodyflap,
                    vehicleName,
                    centralBodyName ) );



    //! Create and initialize object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dep_varList );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE TERMINATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Create Termination Settings" << std::endl;

    //! Define CUSTOM termination settings.
    // std::shared_ptr< PropagationTerminationSettings > customTermination =
    //         std::make_shared< PropagationCustomTerminationSettings >(
    //             boost::bind( &bislip::StopOrNot, bodyMap, vehicleName, vehicleParameterValues_, terminationConditionsValues_ ) );

    //! Define dependent variable termination settings.
    std::shared_ptr< PropagationTerminationSettings > thrustTerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleAccelerationDependentVariableSaveSettings >(
                    thrust_acceleration,
                    vehicleName,
                    vehicleName,
                    true ),bislipSystems->getMaxThrust(), false );

    std::shared_ptr< PropagationTerminationSettings > d_to_go_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    angular_distance_to_go,
                    vehicleName),bislipSystems->getFinalDistanceToTarget(), true );

    std::shared_ptr< PropagationTerminationSettings > d_traveled_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    angular_distance_traveled,
                    vehicleName),bislipSystems->getInitialDistanceToTarget(), false );

    std::shared_ptr< PropagationTerminationSettings > mass_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< tudat::propagators::SingleDependentVariableSaveSettings >(
                    current_mass,
                    vehicleName),bislipSystems->getLandingMass(), true );

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
                    centralBodyName), constraint_DynamicPressure*1.5 , false );

    std::shared_ptr< PropagationTerminationSettings > q_dot_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    local_aerodynamic_heat_rate_dependent_variable,
                    vehicleName,
                    centralBodyName), constraint_HeatingRate , false );

    // std::shared_ptr< PropagationTerminationSettings > height_UP_TerminationSettings =
    //         std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
    //             std::make_shared< SingleDependentVariableSaveSettings >(
    //                 altitude_dependent_variable,
    //                 vehicleName,
    //                 centralBodyName), h_UP , false );

    std::shared_ptr< PropagationTerminationSettings > height_DN_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable,
                    vehicleName,
                    centralBodyName), h_DN , true );

    std::shared_ptr< PropagationTerminationSettings > pos_flight_path_angle_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::flight_path_angle,
                    centralBodyName) , unit_conversions::convertDegreesToRadians( 0.0 ) , false );

    std::shared_ptr< PropagationTerminationSettings > neg_flight_path_angle_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    vehicleName,
                    reference_frames::AerodynamicsReferenceFrameAngles::flight_path_angle,
                    centralBodyName) , unit_conversions::convertDegreesToRadians( 0.0 ) , true );

    std::shared_ptr< PropagationTerminationSettings > total_aero_g_load_TerminationSettings =
            std::make_shared< tudat::propagators::PropagationDependentVariableTerminationSettings >(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_aerodynamic_g_load_variable,
                    vehicleName,
                    centralBodyName), constraint_MechanicalLoad*1.5 , false );

    //! Create list of terminations setting for ASCENT.
    std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettingsList_Ascent;
    // terminationSettingsList.push_back( customTermination );
    //terminationSettingsList.push_back( thrustTerminationSettings );
    //terminationSettingsList_Ascent.push_back( d_to_go_TerminationSettings );
    //terminationSettingsList_Ascent.push_back( d_traveled_TerminationSettings );
    //terminationSettingsList_Ascent.push_back( total_aero_g_load_TerminationSettings );
    //terminationSettingsList.push_back( E_hat_TerminationSettings );
    // terminationSettingsList_Ascent.push_back( q_dyn_TerminationSettings );
    //terminationSettingsList.push_back( q_dot_TerminationSettings );
    //    terminationSettingsList.push_back( height_UP_TerminationSettings );
    //terminationSettingsList_Ascent.push_back( height_DN_TerminationSettings );
    terminationSettingsList_Ascent.push_back( neg_flight_path_angle_TerminationSettings );

    //! Finalize terminations setting for ASCENT.
    std::shared_ptr< PropagationTerminationSettings > terminationSettings_Ascent = std::make_shared<
            propagators::PropagationHybridTerminationSettings >( terminationSettingsList_Ascent, true );

    //! Create list of terminations setting for DESCENT.
    std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettingsList_Descent;
    // terminationSettingsList_Descent.push_back( customTermination );
    //terminationSettingsList_Descent.push_back( thrustTerminationSettings );
    terminationSettingsList_Descent.push_back( d_to_go_TerminationSettings );
    terminationSettingsList_Descent.push_back( d_traveled_TerminationSettings );
    //terminationSettingsList_Descent.push_back( total_aero_g_load_TerminationSettings );
    //terminationSettingsList_Descent.push_back( mass_TerminationSettings );
    //terminationSettingsList_Descent.push_back( E_hat_TerminationSettings );
    //terminationSettingsList_Descent.push_back( q_dyn_TerminationSettings );
    //terminationSettingsList_Descent.push_back( q_dot_TerminationSettings );
    //    terminationSettingsList_Descent.push_back( height_UP_TerminationSettings );
    terminationSettingsList_Descent.push_back( height_DN_TerminationSettings );
    //terminationSettingsList_Descent.push_back( pos_flight_path_angle_TerminationSettings );

    //! Finalize terminations setting for DESCENT.
    std::shared_ptr< PropagationTerminationSettings > terminationSettings_Descent = std::make_shared<
            propagators::PropagationHybridTerminationSettings >( terminationSettingsList_Descent, true );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE MASS RATE SETTINGS: ASCENT            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    //std::cout << "Create mass rate models" << std::endl;

    //! Declare and initialize mass rate model settings.
    std::shared_ptr< MassRateModelSettings > massRateModelSettings =
            std::make_shared< FromThrustMassModelSettings >( true );

    //! Declare and initialize mass rate model settings.
    std::map< std::string, std::shared_ptr< MassRateModel > > massRateModels;
    massRateModels[ vehicleName ] = createMassRateModel(
                vehicleName, massRateModelSettings, bodyMap, accelerationsMap );

    //! Create settings for propagating the mass of the vehicle.
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back( vehicleName );

    //! Declare and initialize starting mass of vehicle.
    Eigen::VectorXd initialBodyMasses_Ascent( 1 );
    initialBodyMasses_Ascent( 0 ) = initialMass_Ascent;

    //! Declare and initialize mass propagation settings.
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings_Ascent =
            std::make_shared< MassPropagatorSettings< double > >(
                bodiesWithMassToPropagate, massRateModels, initialBodyMasses_Ascent, terminationSettings_Ascent, dependentVariablesToSave );

*/
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PAGMO PROBLEM TO COMPUTE FITNESS            ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Create object to compute the problem fitness.
    pagmo::problem prob{ Space4ErrBodyProblem( bounds,
                                               problem_name,
                                               vehicleName,
                                               parameterList_Ascent,
                                               parameterBounds_Ascent,
                                               parameterList_Descent,
                                               parameterBounds_Descent,
                                               vehicleParameterValues,
                                               aeroCoeffFileList,
                                               simulation_settingsValues,
                                               initialConditionsValues,
                                               terminationConditionsValues,
                                               output_settingsValues,
                                               outputSubFolder,
                                               initialState_spherical,
                                               centralBodies,
                                               bodiesToIntegrate,
                                               Bounds_Ascent,
                                               Bounds_Descent,
                                               bodyMap,
                                               accelerationsMap,
                                               earthRotationalEphemeris,
                                               dependentVariablesToSave,
                                               terminationSettings_Ascent,
                                               terminationSettings_Descent ) };

    //! Retrieve algorithm. Three options available in the following function:
    //!        getMultiObjectiveAlgorithm
    //!               case 0 --> nsga2
    //!               case 1 --> moead
    //!               case 2 --> ihs
    //! Selection is currently arbitrary. moead.hpp should be modified such that
    //! the points are generated with a Halton sequence.

    pagmo::algorithm algo{ getPagmoAlgorithm( int( optimization_settingsValues[ 0 ] ) ) };

    const int index = int( optimization_settingsValues[ 0 ] );
    std::string algo_method;
    switch( index )
    {
    case 0:
    {
        algo_method = "NSGA2";
        break;
    }
    case 1:
    {
        algo_method = "MOEAD";
        break;
    }
    case 2:
    {
        algo_method = "IHS";
        break;
    }
    }

    //! Assign population size.
    pagmo::population::size_type populationSize = int( optimization_settingsValues[ 1 ] );

    //! Assign archipelago size.
    //pagmo::archipelago::size_type archipelagoSize = int( optimization_settingsValues[ 2 ] );

    //! Assign population per archipelago
    //const int pop_per_archi = populationSize/archipelagoSize;

    //! Assign number of evolutions
    const int evolutions =  int( optimization_settingsValues[ 3 ] );

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    /*
    const int index = int( optimization_settingsValues[ 0 ] );
    std::string algo_method;
    switch( index )
    {
    case 0:
    {
        algo_method = "NSGA2";
        break;
    }
    case 1:
    {
        algo_method = "MOEAD";
        break;
    }
    case 2:
    {
        algo_method = "IHS";
        break;
    }
    }
    */
    //! Convert angles from degrees to radians
    const double lat_i_rad = unit_conversions::convertDegreesToRadians( initialConditionsValues[ 0 ] );
    const double lon_i_rad = unit_conversions::convertDegreesToRadians( initialConditionsValues[ 1 ] );
    const double lat_f_rad = unit_conversions::convertDegreesToRadians( terminationConditionsValues[ 0 ] );
    const double lon_f_rad = unit_conversions::convertDegreesToRadians( terminationConditionsValues[ 1 ] );

    //! Calculate initial heading angle: https://www.movable-type.co.uk/scripts/latlong.html
    double chi_i_deg_calc = unit_conversions::convertRadiansToDegrees( bislip::Variables::computeHeadingToTarget( lat_i_rad , lon_i_rad , lat_f_rad , lon_f_rad ) );

    //! If heading angle is negative, this may help visualize it.
    if (chi_i_deg_calc < 0)
    {
        chi_i_deg_calc = 360 + chi_i_deg_calc;
    }

    //! Calculate ground distance to cover using Haversine formula and Sperical Law of Cosines: https://www.movable-type.co.uk/scripts/latlong.html
    const double a = std::sin( (lat_f_rad - lat_i_rad) / 2) * std::sin( (lat_f_rad - lat_i_rad) / 2) + std::cos( lat_i_rad ) * std::cos( lon_i_rad ) * std::sin( (lon_f_rad - lon_i_rad) / 2) * std::sin( (lon_f_rad - lon_i_rad) / 2);
    const double d_angular = 2 * std::atan2( std::sqrt(a) , std::sqrt(1 - a) );
    const double d_spherical_law_cosines = unit_conversions::convertRadiansToDegrees( bislip::Variables::computeAngularDistance( lat_i_rad , lon_i_rad , lat_f_rad , lon_f_rad ) );
    //std::acos( std::sin(lat_i_rad) * std::sin(lat_f_rad) + std::cos(lat_i_rad) * std::cos(lat_f_rad) * std::cos(lon_f_rad-lon_i_rad) ) * spice_interface::getAverageRadius( "Earth" );

    //! Print to Terminal screen.
    std::cout << " " << std::endl;
    std::cout << problem_name << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Parameters to vary during Ascent" << std::endl;
    p=0;
    for( int i = 0; i < int( parameterList_Ascent.size() ) ; i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.1f %-2s %-8.1f %-2s\n") % parameterList_Ascent[i] % "[" % parameterBounds_Ascent[p] % "," % parameterBounds_Ascent[p+1] % "]";
        p+=2;
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Parameters to vary during Descent" << std::endl;
    p=0;
    for( int i = 0; i < int( parameterList_Descent.size() ) ; i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.1f %-2s %-8.1f %-2s\n") % parameterList_Descent[i] % "[" % parameterBounds_Descent[p] % "," % parameterBounds_Descent[p+1] % "]";
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
    for( int i = 0; i < int( optimization_settingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.0f\n") % optimization_settingsList[i] % "" % optimization_settingsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Simulation Settings" << std::endl;
    for( int i = 0; i < int( simulation_settingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.0f\n") % simulation_settingsList[i] % "" % simulation_settingsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Output Settings" << std::endl;
    for( int i = 0; i < int( output_settingsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.0f\n") % output_settingsList[i] % "" % output_settingsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Initial Conditions" << std::endl;
    for( int i = 0; i < int( initialConditionsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.2f\n") % initialConditionsList[i] % "" % initialConditionsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Termination Conditions" << std::endl;
    for( int i = 0; i < int( terminationConditionsList.size() ); i++)
    {
        std::cout << boost::format("%-25s %-2s %-8.2f\n") % terminationConditionsList[i] % "" % terminationConditionsValues[i];
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Ground distance to cover " << std::endl;
    std::cout << "  Haversine Formula:        " << unit_conversions::convertRadiansToDegrees( d_angular ) << " degrees." << std::endl;
    std::cout << "  Spherical Law of Cosines: " << d_spherical_law_cosines << " degrees." << std::endl;
    std::cout << "Heading angle to Target:    " << chi_i_deg_calc << " degrees. Calculated."  << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Output subfolder: " << std::endl;
    std::cout << "     '" << outputSubFolder <<"'" << std::endl;
    std::cout << " DV   " << "  Lower Boundary " << "  Upper Boundary " << std::endl;
    for( int i = 0; i < N + NN + 3 + 1; i++ )
    {
        std::cout << std::fixed << std::setprecision(10) <<
                     std::setw(6) << i <<
                     std::setw(17) << bounds[ 0 ][ i ] <<
                     std::setw(17) << bounds[ 1 ][ i ] << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             EVALUATE PAGMO PROBLEM                             ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //! Create an island with populationSize individuals
    pagmo::island isl{algo, prob, populationSize};

    //! Instantiate an archipelago
    /*  archipelago archi{archipelagoSize, algo, prob, pop_per_archi};

    //! Evolve
    archi.evolve(evolutions);

    //! Wait for the evolutions to be finished
    archi.wait_check();// Raises errors,

    for (const auto &isl : archi) {
        std::cout << isl.get_population().champion_f()[0] << '\n';
    }
    for (const auto &isl : archi) {
        std::cout << isl.get_population().champion_x()[0] << '\n';
    }
    for (const auto &isl : archi) {
     printPopulationToFile( isl.get_population().champion_f()[0],
             "mo_AMSIAD_BALLISTIC_" ,
             false );
     printPopulationToFile( isl.get_population().champion_x()[0],
             "mo_AMSIAD_BALLISTIC_" ,
             true );
    }
*/
    //! Write original (unevolved) population to file
    if  ( int( output_settingsValues[ 0 ] ) == 1 )
    {
        //   printPopulationToFile( isl_pert.get_population( ).get_x( ), this_run_settings_1,
        //                          outputSubFolder, false );
        //   printPopulationToFile( isl_pert.get_population( ).get_f( ), this_run_settings_1,
        //                          outputSubFolder, true );
    }

    for( int i = 0; i < evolutions; i++ )
    {
        isl.evolve( );
        while( isl.status( ) != pagmo::evolve_status::idle &&
               isl.status( ) != pagmo::evolve_status::idle_error )
        {
            isl.wait( );
        }
        isl.wait_check( ); // Raises errors

        if  ( int( output_settingsValues[ 0 ] ) == 1 )
        {

            std::string this_run_settings_1 = std::to_string( i ) + "_" +
                    this_run_settings;

            //! Write current iteration results to file
            printPopulationToFile( isl.get_population( ).get_x( ),
                                   this_run_settings_1,
                                   outputSubFolder,
                                   false );
            printPopulationToFile( isl.get_population( ).get_f( ),
                                   this_run_settings_1,
                                   outputSubFolder,
                                   true ); //! true boolean is for the printPopulationToFile to print with a different prefix: fitness_
        }

        std::cout<<i<<std::endl;

    }

}

