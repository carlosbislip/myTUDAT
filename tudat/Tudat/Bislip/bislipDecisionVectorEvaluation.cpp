#include <Tudat/Bislip/bislipVariables.h>
#include <Tudat/Bislip/bislipProblemInput.h>
#include <Tudat/Bislip/updateGuidance.h>
//#include <../pagmo2/include/pagmo/pagmo.hpp>
//#include <../pagmo2/include/pagmo/utils/multi_objective.hpp>

#include <pagmo/pagmo.hpp>
#include <pagmo/types.hpp>
#include <pagmo/utils/multi_objective.hpp>

//#include <Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h>
//#include <Tudat/External/SpiceInterface/spiceInterface.h>

//#include <Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h>
namespace bislip { namespace Variables {

//! Evaluate parameter vector
void decisionVectorEvaluation( const std::vector< double > &x,
                               const std::shared_ptr< bislip::ProblemInput > &problemInput,
                               const tudat::simulation_setup::NamedBodyMap& bodyMap,
                               Eigen::MatrixXd &depVarTimeHistoryMatrix,
                               int &rowsAscent )
{

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    using namespace tudat::mathematical_constants;
    using namespace tudat::input_output;
    using namespace tudat::unit_conversions;
    using namespace tudat::reference_frames;
    using namespace tudat;
    //using namespace tudat_applications;
    using namespace tudat::aerodynamics;
    using namespace bislip;
    //using namespace std::this_thread; // sleep_for, sleep_until
    //using namespace std::chrono; // nanoseconds, system_clock, seconds

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            UNPACK INPUT DATA             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const std::string vehicleName =  problemInput->getVehicleName();

    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap.at( vehicleName )->getVehicleSystems();

    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at(  vehicleName )->getBislipSystems();

    bislipSystems->setStartingMillis( bislip::Variables::millis_since_midnight() );

    bislipSystems->setDebugInfo( problemInput->getSimulationSettings()[ 9 ] );

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 10 )
    {
        std::cout << "Printing decision vector to screen" << std::endl;
        for( int i = 0; i < int( x.size( ) ); i++ ) { std::cout << "x[ " << i << " ] = " << x[ i ] << std::endl; }
    }

    if( debugInfo == 10 ){std::cout << "Unpacking data" << std::endl; }

    const std::string centralBodyName = problemInput->getCentralBodies()[ 0 ];

    //! Various parameters
    const double radiusEarth = tudat::spice_interface::getAverageRadius( centralBodyName );
    const double g0 = tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION;
    const double pi = tudat::mathematical_constants::PI;
    const double omega_E = 7.292115*1E-5;
    const double J2 = 1082.626523e-6;
    const double J3 = 2.532153e-7;
    const double J4 = 1.6109876e-7;

    //! Declare and initialize simulation start epoch.
    const double simulationStartEpoch = bislipSystems->getStartingEpoch();

    //! Declare and initialize simulation end epoch.
    const double simulationEndEpoch = simulationStartEpoch + problemInput->getSimulationSettings()[ 6 ];

    //! Declare and initialize numerical integration fixed step size.
    const double propagationStepSize = problemInput->getSimulationSettings()[ 7 ];

    const double guidanceStepSize = problemInput->getSimulationSettings()[ 8 ];

    //! Declare and initialize number of control nodes.
    const unsigned long nodesAscent = problemInput->getSimulationSettings().rbegin()[ 1 ];
    const unsigned long nodesDescent = problemInput->getSimulationSettings().rbegin()[ 0 ];

    if( debugInfo == 10 ){std::cout << "nodesAscent = " << nodesAscent <<std::endl; }
    if( debugInfo == 10 ){std::cout << "nodesDescent = " << nodesDescent <<std::endl; }


    //! Declare and initialize position vector of moment reference center
    const Eigen::Vector3d R_mrc( problemInput->getVehicleParameters()[ 3 ], problemInput->getVehicleParameters()[ 4 ], problemInput->getVehicleParameters()[ 5 ] ); // m

    //! Declare and initialize position vector of center of mass
    const Eigen::Vector3d R_com( problemInput->getVehicleParameters()[ 6 ], problemInput->getVehicleParameters()[ 7 ], problemInput->getVehicleParameters()[ 8 ] ); // m

    //! Declare and initialize position vector of center of thrust
    const Eigen::Vector3d R_cot( problemInput->getVehicleParameters()[ 9 ], problemInput->getVehicleParameters()[ 10 ], problemInput->getVehicleParameters()[ 11 ] ); // m

    //! Declare and initialize initial mass
    double initialMass_Ascent = problemInput->getVehicleParameters()[ 12 ]; // kg

    double dryMass = problemInput->getVehicleParameters()[ 13 ]; // kg

    //! Declare and initialize starting height
    const double initialHeight_Ascent = problemInput->getInitialConditions()[ 2 ]; // m

    //! Declare and initialize starting position coordinates.
    const double initialLat_deg = problemInput->getInitialConditions()[ 0 ];
    const double initialLon_deg = problemInput->getInitialConditions()[ 1 ];

    //! Declare and initialize final position coordinates and additional termination conditions
    const double targetLat_deg = problemInput->getConstraints()[ 0 ];
    const double targetLon_deg = problemInput->getConstraints()[ 1 ];

    //! Declare and initialize various constraints
    const double finalDistanceToTarget_deg         = problemInput->getConstraints()[ 2 ];
    const double objectiveHeight_Descent           = problemInput->getConstraints()[ 4 ];
    //const double constraint_MechanicalLoad         = problemInput->getConstraints()[ 6 ];
    const double constraint_ChapmanHeatFlux        = problemInput->getConstraints()[ 7 ];
    const double constraint_DynamicPressure        = problemInput->getConstraints()[ 8 ];
    const double constraint_PitchMomentCoefficient = problemInput->getConstraints()[ 9 ];
    const double constraint_BendingMoment          = problemInput->getConstraints()[ 10 ];
    const double constraint_PassengerPosZLoad      = problemInput->getConstraints()[ 11 ];
    const double constraint_PassengerNegZLoad      = problemInput->getConstraints()[ 12 ];
    const double constraint_PassengerJerk          = problemInput->getConstraints()[ 13 ];

    //! Convert angles from degrees to radians
    const double initialLat_rad             = unit_conversions::convertDegreesToRadians( initialLat_deg );
    const double initialLon_rad             = unit_conversions::convertDegreesToRadians( initialLon_deg );
    const double targetLat_rad              = unit_conversions::convertDegreesToRadians( targetLat_deg );
    const double targetLon_rad              = unit_conversions::convertDegreesToRadians( targetLon_deg );

    //! Pre-define various variables used to determine fitness.
    double initialDistanceToTarget_rad = bislip::Variables::computeAngularDistance( initialLat_rad, initialLon_rad, targetLat_rad, targetLon_rad );
    double initialDistanceToTarget_deg = unit_conversions::convertRadiansToDegrees( initialDistanceToTarget_rad );
    double initialDistanceToTarget_m   = radiusEarth * initialDistanceToTarget_rad;


    if( debugInfo == 10 ){ std::cout << "Reset non-static placeholders." << std::endl; }
    //! Reset value placeholders that store non-static values.
    bislipSystems->setCurrentFlightPathAngleRate( 0.0 );
    bislipSystems->setCurrentBankAngle( 0.0 );
    bislipSystems->setTempBankAngle( 0.0 );
    bislipSystems->setReversalConditional( 0.0 );
    bislipSystems->setBankAngleReversalTrigger( false );
    bislipSystems->setLowDistanceReversalCompleted( false );
    bislipSystems->setWallTemperature( 0.0 );
    bislipSystems->setChapmanWallTemp( 0.0 );
    bislipSystems->setTauberWallTempStagnation( 0.0 );
    bislipSystems->setTauberWallTempFlatPlate( 0.0 );
    bislipSystems->setCurrentFlightPathAngleRate( 0.0 );
    bislipSystems->setCurrentThrottleSetting( NAN );
    bislipSystems->setCurrentBodyFlapAngle( tudat::unit_conversions::convertDegreesToRadians( problemInput->getInitialConditions()[ 5 ] ) );
    bislipSystems->setCurrentElevonAngle( tudat::unit_conversions::convertDegreesToRadians( problemInput->getInitialConditions()[ 6 ] ) );
    bislipSystems->setCumulativeDistanceTravelled( 0.0 );
    bislipSystems->setCumulativeAngularDistanceTravelled( 0.0 );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE NODAL STRUCTURE             /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){std::cout << "Creating nodal structure" << std::endl; }

    //! Declare and allocate vectors of interest for Ascent phase.
    Eigen::VectorXd xn_interval_Ascent( nodesAscent - 1 );
    Eigen::VectorXd xn_Ascent( nodesAscent );
    Eigen::VectorXd alpha_deg_Ascent( nodesAscent );
    Eigen::VectorXd sigma_deg_Ascent( nodesAscent );
    Eigen::VectorXd eps_T_deg_Ascent( nodesAscent );
    Eigen::VectorXd phi_T_deg_Ascent( nodesAscent );
    Eigen::VectorXd throttle_Ascent( nodesAscent );

    //! Declare and allocate vectors of interest for Descent phase.
    Eigen::VectorXd xn_interval_Descent( nodesDescent - 1 );
    Eigen::VectorXd xn_Descent( nodesDescent );
    Eigen::VectorXd alpha_deg_Descent( nodesDescent );
    Eigen::VectorXd sigma_deg_Descent( nodesDescent );
    Eigen::VectorXd eps_T_deg_Descent( nodesDescent );
    Eigen::VectorXd phi_T_deg_Descent( nodesDescent );
    Eigen::VectorXd throttle_Descent( nodesDescent );

    if( debugInfo == 10 ){std::cout << "     Re-allocate Ascent DVs" << std::endl; }

    //! Re-allocate decision vector values into workable vectors for Ascent phase.
    for( unsigned int i = 0; i < nodesAscent; i++ )
    {
        if ( i < ( nodesAscent - 1 ) ){ xn_interval_Ascent( i ) = x[ i ]; }
        alpha_deg_Ascent( i ) = x[ i + 1 * nodesAscent - 1 ];
        sigma_deg_Ascent( i ) = x[ i + 2 * nodesAscent - 1 ];
        eps_T_deg_Ascent( i ) = x[ i + 3 * nodesAscent - 1 ];
        phi_T_deg_Ascent( i ) = x[ i + 4 * nodesAscent - 1 ];
        throttle_Ascent( i )  = x[ i + 5 * nodesAscent - 1 ];
    }

    //! Declare and initialize number of parameters exclusive to Ascent phase.
    //    const unsigned long N = ( problemInput->getAscentParameterList().size() - 6 ) * nodesAscent - 1;
    const unsigned long N = problemInput->getAscentParameterList().size();// - 6 ) * nodesAscent - 1;

    //! Declare and initialize various parameters common to the entire trajectory.
    double initialFlightPathAngle               = x[ N - 7 ];
    const double initialLaunchHeadingAngle      = x[ N - 6 ];
    double initialAirspeed_Ascent               = x[ N - 5 ];
    const double objectiveAirspeed_Ascent       = x[ N - 4 ];
    const double objectiveHeight_Ascent         = x[ N - 3 ];
    double additionalMass                       = x[ N - 2 ];
    const double ascentTerminationDistanceRatio = x[ N - 1 ];




    bislipSystems->setInitialHeight( initialHeight_Ascent );
    bislipSystems->setInitialAltitude( initialHeight_Ascent + tudat::spice_interface::getAverageRadius( centralBodyName ) );

    bislipSystems->setInitialValueFlag( true );
    if( bislipSystems->getValidationFlag( ) == true )
    {

        if( debugInfo == 10 ){ std::cout << "    Setting Trajectory Phase" << std::endl; }
        bislipSystems->setCurrentTrajectoryPhase( "Descent" );
        if( debugInfo == 10 ){std::cout << "         Trajectory Phase         = " << bislipSystems->getCurrentTrajectoryPhase() << std::endl; }

        initialFlightPathAngle = bislipSystems->getInitialFlightPathAngle();

    }
    else
    {
        if( debugInfo == 10 ){ std::cout << "    Setting Trajectory Phase" << std::endl; }
        bislipSystems->setCurrentTrajectoryPhase( "Ascent" );
        if( debugInfo == 10 ){std::cout << "         Trajectory Phase         = " << bislipSystems->getCurrentTrajectoryPhase() << std::endl; }

        if( debugInfo == 10 ){std::cout << "    Identifying additional initial conditions" << std::endl; }

        if( bislipSystems->getInitialMachNumber() > 0.0 && bislipSystems->getInitialHeight() > 0.0 )
        {

            if( debugInfo == 10 ){std::cout << "    Vehicle initialized with Speed and Height" << std::endl; }

            initialAirspeed_Ascent = bislipSystems->getInitialSpeedOfSound() * bislipSystems->getInitialMachNumber();
            bislipSystems->setInitialAirspeed( initialAirspeed_Ascent );
        }
        else if( bislipSystems->getInitialMachNumber() == 0.0 && bislipSystems->getInitialHeight() != 0.0 )
        {
            if( debugInfo == 10 ){std::cout << "    Vehicle initialized at rest and with height---> Initial Airspeed is an optimization parameter" << std::endl; }

            bislipSystems->setInitialMachNumber( initialAirspeed_Ascent / bislipSystems->getInitialSpeedOfSound() );
        }
        else if( bislipSystems->getInitialMachNumber() == 0.0 && bislipSystems->getInitialHeight() == 0.0 )
        {
            if( debugInfo == 10 ){std::cout << "    Vehicle initialized at rest and on the surface" << std::endl; }

            bislipSystems->setInitialMachNumber( initialAirspeed_Ascent / bislipSystems->getInitialSpeedOfSound() );
        }
    }

    bislipSystems->setInitialDynamicPressure( bislipSystems->getInitialDensity() * initialAirspeed_Ascent * initialAirspeed_Ascent / 2 );


    if( debugInfo == 10 ){std::cout << "         Initial Height           = " << bislipSystems->getInitialHeight() << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Initial Altitude         = " << bislipSystems->getInitialAltitude() << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Initial Mach Number      = " << bislipSystems->getInitialMachNumber() << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Initial Airspeed         = " << bislipSystems->getInitialAirspeed() << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Initial Speed of Sound   = " << bislipSystems->getInitialSpeedOfSound() << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Initial Dynamic Pressure = " << bislipSystems->getInitialDynamicPressure() << std::endl; }
    const double goaldelV_Ascent  = objectiveAirspeed_Ascent - initialAirspeed_Ascent;


    //! Set body Mass.
    if( bislipSystems->getValidationFlag() == true )
    {
        initialMass_Ascent = vehicleSystems->getDryMass();
        bislipSystems->setInitialMass( vehicleSystems->getDryMass() );
        additionalMass = 0.0;
    }
    else
    {
        if( additionalMass > 0.0 )
        {
            double increasedVehicleDryMass = dryMass * ( 1.0 + 0.3*( ( additionalMass / initialMass_Ascent ) ) );
            vehicleSystems->setDryMass( increasedVehicleDryMass );
        }
    }

    if( debugInfo == 10 ){std::cout << "         Initial Launch Heading    = " << initialLaunchHeadingAngle << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Initial Flight-Path Angle = " << initialFlightPathAngle << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Initial Launch Airspeed   = " << initialAirspeed_Ascent << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Objective Airspeed        = " << objectiveAirspeed_Ascent << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Objective Height          = " << objectiveHeight_Ascent << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Original Mass - Ascent    = " << initialMass_Ascent << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Additional Mass           = " << additionalMass << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Initial Mass - Ascent     = " << initialMass_Ascent + additionalMass << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Ascent Termination Ratio  = " << ascentTerminationDistanceRatio << std::endl; }

    bodyMap.at( vehicleName )->setConstantBodyMass( initialMass_Ascent + additionalMass );
    bislipSystems->setInitialMass( initialMass_Ascent + additionalMass );
    bislipSystems->setInitialHeight( initialHeight_Ascent );
    bislipSystems->setInitialAltitude( initialHeight_Ascent + tudat::spice_interface::getAverageRadius( centralBodyName ) );
    bislipSystems->setAscentTerminationDistanceRatio( ascentTerminationDistanceRatio );

    if( debugInfo == 10 ){std::cout << "     Re-allocate Descent DVs" << std::endl; }

    //! Re-allocate decision vector values into workable vectors for Descent phase.
    for( unsigned int i = 0; i < nodesDescent; i++ )
    {
        if ( i < ( nodesDescent - 1) ){ xn_interval_Descent( i ) = x[ ( N + 0 ) + i ]; }
        alpha_deg_Descent( i ) = x[ ( N + 0 ) + i + 1 * nodesDescent - 1 ];
        sigma_deg_Descent( i ) = x[ ( N + 0 ) + i + 2 * nodesDescent - 1 ];
        eps_T_deg_Descent( i ) = x[ ( N + 0 ) + i + 3 * nodesDescent - 1 ];
        phi_T_deg_Descent( i ) = x[ ( N + 0 ) + i + 4 * nodesDescent - 1 ];
        throttle_Descent( i )  = x[ ( N + 0 ) + i + 5 * nodesDescent - 1 ];
    }

    //! Declare and initialize number of parameters exclusive to Descent phase.
    //const unsigned long NN = ( problemInput->getDescentParameterList().size() - 2 ) * nodesDescent - 1;
    const unsigned long NN = problemInput->getDescentParameterList().size();// - 2 ) * nodesDescent - 1;

    //! Declare and initialize last parameter common to the entire trajectory.
    //!     Final velocity.
    //!     Skip suppression timing trigger.
    const double objectiveAirspeed_Descent    = x.rbegin()[ 2 ];
    const double skipSuppressionTimingTrigger = x.rbegin()[ 1 ];
    const double constraint_MechanicalLoad    = x.rbegin()[ 0 ];

    bislipSystems->setMechanicalLoadConstraint( constraint_MechanicalLoad );

    /*
    if( debugInfo == 10 ){std::cout << "alpha_deg_Descent            = " << alpha_deg_Descent << std::endl; }
    if( debugInfo == 10 ){std::cout << "sigma_deg_Descent            = " << sigma_deg_Descent << std::endl; }
    if( debugInfo == 10 ){std::cout << "eps_T_deg_Descent            = " << eps_T_deg_Descent << std::endl; }
    if( debugInfo == 10 ){std::cout << "phi_T_deg_Descent            = " << phi_T_deg_Descent << std::endl; }
    if( debugInfo == 10 ){std::cout << "throttle_Descent             = " << throttle_Descent << std::endl; }
    if( debugInfo == 10 ){std::cout << "Objective Airspeed Descent   = " << objectiveAirspeed_Descent << std::endl; }
    if( debugInfo == 10 ){std::cout << "Skip Suppresion Trigger Time = " << skipSuppressionTimingTrigger << std::endl; }
    */
    bislipSystems->setSkipSuppressionTimingTrigger( skipSuppressionTimingTrigger );

    const double goaldelV_Descent = std::abs( objectiveAirspeed_Descent - objectiveAirspeed_Ascent );

    if( debugInfo == 10 ){std::cout << "Create vector of node locations for ascent" << std::endl; }
    //! Create vector of node locations for Ascent phase.
    xn_Ascent( 0 ) = 0;
    for( unsigned int i = 0; i < nodesAscent - 1; i++ )
    {
        xn_Ascent( i + 1 )        = xn_Ascent( i ) + xn_interval_Ascent( i ) / xn_interval_Ascent.sum();
    }

    if( debugInfo == 10 ){std::cout << "Create vector of node locations for descent" << std::endl; }
    //! Create vector of node locations for Descent phase.
    xn_Descent( 0 ) = 0;//xn_Ascent( xn_Ascent.size() - 1 ) ;//xn_interval_Descent.sum();
    for( unsigned int i = 0; i < nodesDescent - 1; i++ )
    {
        xn_Descent( i + 1 )        = xn_Descent( i ) + xn_interval_Descent( i ) / xn_interval_Descent.sum();
    }


    //  std::cout << "xn_interval_Descent =  " << xn_interval_Descent << std::endl;
    // std::cout << "throttle_Descent =  " << throttle_Descent << std::endl;
    // std::cout << "objectiveAirspeed_Descent =  " << objectiveAirspeed_Descent << std::endl;

    // xn_Descent = xn_Descent.reverse().eval();
    /*

    std::cout << "x =  " << std::endl;
    for( int i = 0; i < int( x.size( ) ); i++ )
    {
        std::cout << i << "   " << x[i] << std::endl;
    }

 *
 * std::cout << "-------" << std::endl;
    std::cout << "xn_interval_Ascent =  " << xn_interval_Ascent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "xn_Ascent=  " << xn_Ascent<< std::endl;
    std::cout << "-------" << std::endl;

    std::cout << "alpha_deg_Ascent =  " << alpha_deg_Ascent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "eps_T_deg =  " << eps_T_deg << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "throttle =  " << throttle << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "xn_interval_Descent =  " << xn_interval_Descent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "xn_Descent=  " << xn_Descent<< std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "alpha_deg_Descent =  " << alpha_deg_Descent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "sigma_deg_Descent =  " << sigma_deg_Descent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "initialAirspeed_Ascent =  " << initialAirspeed_Ascent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "objectiveAirspeed_Ascent =  " << objectiveAirspeed_Ascent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "objectiveAirspeed_Descent =  " << objectiveAirspeed_Descent << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "objectiveHeight_Ascent =  " << objectiveHeight_Ascent << std::endl;
    std::cout << "-------" << std::endl;
*/
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            COMPLETE KNOWN ASCENT STATES            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){std::cout << "Complete known ascent states" << std::endl; }

    Eigen::Vector6d initialStateAscent = problemInput->getInitialState_Spherical();
    initialStateAscent( 3 ) = initialAirspeed_Ascent;
    initialStateAscent( 4 ) = tudat::unit_conversions::convertDegreesToRadians( initialFlightPathAngle );
    initialStateAscent( 5 ) = tudat::unit_conversions::convertDegreesToRadians( initialLaunchHeadingAngle );

    bislipSystems->setInitialHeadingAngle( initialLaunchHeadingAngle );
    bislipSystems->setInitialFlightPathAngle( initialFlightPathAngle );

    //! Two things are being done here
    //!     Converting state vector from spherical to Cartesian elements
    //!     Transforming state vector from Earth-Fixed frame to Inertial frame.
    Eigen::Vector6d bodyCenteredBodyFixedState =
            tudat::orbital_element_conversions::convertSphericalOrbitalToCartesianState( initialStateAscent );

    //! Initizalize (via resetting) previous coordinate placeholder
    bislipSystems->setPreviousCoordinates( initialStateAscent.segment( 1, 2 ) );
    bislipSystems->setPreviousCartesianCoordinates( bodyCenteredBodyFixedState.segment( 0, 3 ) );

    const Eigen::Vector6d systemInitialState_Ascent =
            tudat::ephemerides::transformStateToGlobalFrame(
                bodyCenteredBodyFixedState, simulationStartEpoch, problemInput->getEarthRotationalEphemeris() );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE INTERPOLATORS             ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){ std::cout << "Creating interpolators" << std::endl; }

    //! Calculate initial, maximum, and final specific energy levels.
    const double E_i = bislip::Variables::computeSpecificEnergy( initialHeight_Ascent, initialAirspeed_Ascent );
    const double E_max = bislip::Variables::computeSpecificEnergy( objectiveHeight_Ascent, objectiveAirspeed_Ascent );
    const double E_f = bislip::Variables::computeSpecificEnergy( objectiveHeight_Descent, objectiveAirspeed_Descent );

    //! Set maximum Energy level.
    bislipSystems->setMaximumSpecificEnergy( E_max );

    //! Normalize initial, maximum, and final specific energy levels.
    const double E_hat_i   = bislip::Variables::computeNormalizedSpecificEnergy( initialHeight_Ascent, initialAirspeed_Ascent, E_max );
    const double E_hat_max = bislip::Variables::computeNormalizedSpecificEnergy( objectiveHeight_Ascent, objectiveAirspeed_Ascent, E_max );
    const double E_hat_f   = bislip::Variables::computeNormalizedSpecificEnergy( objectiveHeight_Descent, objectiveAirspeed_Descent, E_max );

    if( debugInfo == 10 ){ std::cout << "E_hat_i   = " << E_hat_i << std::endl; }
    if( debugInfo == 10 ){ std::cout << "E_hat_max = " << E_hat_max << std::endl; }
    if( debugInfo == 10 ){ std::cout << "E_hat_f   = " << E_hat_f << std::endl; }

    //! Map normalized specific energy levels to control node locations.
    Eigen::VectorXd E_mapped_Ascent( xn_Ascent.size() );
    Eigen::VectorXd E_mapped_Descent( xn_Descent.size() );

    E_mapped_Ascent  = ( ( E_max - E_i ) * xn_Ascent.array() + E_i ) / E_max;

    if( bislipSystems->getValidationFlag() == true )
    {
        if( debugInfo == 10 ){ std::cout << "Validation Case ---> Overwriting ascent info with descent info "<< std::endl; }

        xn_Ascent = xn_Descent;
        E_mapped_Ascent  = ( ( E_max - E_f ) * xn_Descent.array() + E_f ) / E_max;

        alpha_deg_Ascent = alpha_deg_Descent;
        sigma_deg_Ascent = sigma_deg_Descent;
        eps_T_deg_Ascent = eps_T_deg_Descent;
        phi_T_deg_Ascent = phi_T_deg_Descent;
        throttle_Ascent  = throttle_Descent;
    }
    else { E_mapped_Ascent  = ( ( E_max - E_i ) * xn_Ascent.array() + E_i ) / E_max; }

    E_mapped_Descent = ( ( E_max - E_f ) * xn_Descent.array() + E_f ) / E_max;

    //E_mapped_Ascent = ( ( E_max - E_i ) * xn_Ascent.array() + E_i ) ;
    //E_mapped_Descent = ( ( E_max - E_f ) * xn_Descent.array() + E_f ) ;
    //E_mapped_Descent = E_mapped_Descent.reverse().eval();
    //std::cout << "xn_Ascent=  " << xn_Ascent<< std::endl;
    //std::cout << "xn_Descent=  " << xn_Descent<< std::endl;
    //std::cout << "E_mapped_Ascent=  " << E_mapped_Ascent<< std::endl;
    //std::cout << "E_mapped_Descent=  " << E_mapped_Descent<< std::endl;

    //! Declare data maps used for decision vector values related to Ascent phase.
    std::map< double, double > map_alpha_deg_Ascent, map_eps_T_deg_Ascent, map_phi_T_deg_Ascent, map_throttle_Ascent, map_sigma_deg_Ascent;
    std::map< double, double > map_alpha_deg_Ascent_LB, map_eps_T_deg_Ascent_LB, map_phi_T_deg_Ascent_LB, map_throttle_Ascent_LB, map_sigma_deg_Ascent_LB;
    std::map< double, double > map_alpha_deg_Ascent_UB, map_eps_T_deg_Ascent_UB, map_phi_T_deg_Ascent_UB, map_throttle_Ascent_UB, map_sigma_deg_Ascent_UB;
    std::map< double, Eigen::VectorXd > map_DV_mapped_Ascent;
    Eigen::VectorXd DV_mapped_Ascent ( 6 );


    std::map< bislip::Parameters::Interpolators, Eigen::MatrixXd > ascentParameterBoundsMap = problemInput->getAscentParameterBoundsMap();


    if( debugInfo == 10 ){ std::cout << "Mapping Ascent DVs" << std::endl; }
    //! Associate decision vector values to mapped normalized specific energy levels within data maps.
    for ( unsigned int i = 0; i < E_mapped_Ascent.size( ); ++i )
    {
        map_alpha_deg_Ascent[ E_mapped_Ascent( i ) ] = DV_mapped_Ascent( 0 ) = alpha_deg_Ascent( i );
        map_sigma_deg_Ascent[ E_mapped_Ascent( i ) ] = DV_mapped_Ascent( 1 ) = sigma_deg_Ascent( i );
        map_eps_T_deg_Ascent[ E_mapped_Ascent( i ) ] = DV_mapped_Ascent( 2 ) = eps_T_deg_Ascent( i );
        map_phi_T_deg_Ascent[ E_mapped_Ascent( i ) ] = DV_mapped_Ascent( 3 ) = phi_T_deg_Ascent( i );
        map_throttle_Ascent[ E_mapped_Ascent( i ) ]  = DV_mapped_Ascent( 4 ) = throttle_Ascent( i );
        DV_mapped_Ascent( 5 ) = xn_Ascent( i );
        map_DV_mapped_Ascent[ E_mapped_Ascent( i ) ] = DV_mapped_Ascent;

        map_alpha_deg_Ascent_LB[ E_mapped_Ascent( i ) ] = ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( i , 0 );
        map_sigma_deg_Ascent_LB[ E_mapped_Ascent( i ) ] = ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( i , 0 );
        map_eps_T_deg_Ascent_LB[ E_mapped_Ascent( i ) ] = ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( i , 0 );
        map_phi_T_deg_Ascent_LB[ E_mapped_Ascent( i ) ] = ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( i , 0 );
        map_throttle_Ascent_LB[ E_mapped_Ascent( i ) ]  = ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( i , 0 );

        map_alpha_deg_Ascent_UB[ E_mapped_Ascent( i ) ] = ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( i , 1 );
        map_sigma_deg_Ascent_UB[ E_mapped_Ascent( i ) ] = ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( i , 1 );
        map_eps_T_deg_Ascent_UB[ E_mapped_Ascent( i ) ] = ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( i , 1 );
        map_phi_T_deg_Ascent_UB[ E_mapped_Ascent( i ) ] = ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( i , 1 );
        map_throttle_Ascent_UB[ E_mapped_Ascent( i ) ]  = ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( i , 1 );
    }

    //! Declare data maps used for decision vector values related to Descent phase.
    std::map< double, double > map_alpha_deg_Descent, map_eps_T_deg_Descent, map_phi_T_deg_Descent, map_throttle_Descent, map_sigma_deg_Descent;
    std::map< double, double > map_alpha_deg_Descent_LB, map_eps_T_deg_Descent_LB, map_phi_T_deg_Descent_LB, map_throttle_Descent_LB, map_sigma_deg_Descent_LB;
    std::map< double, double > map_alpha_deg_Descent_UB, map_eps_T_deg_Descent_UB, map_phi_T_deg_Descent_UB, map_throttle_Descent_UB, map_sigma_deg_Descent_UB;
    std::map< double, Eigen::VectorXd > map_DV_mapped_Descent;
    Eigen::VectorXd DV_mapped_Descent ( 6 );

    std::map< bislip::Parameters::Interpolators, Eigen::MatrixXd > descentParameterBoundsMap = problemInput->getDescentParameterBoundsMap();



    if( debugInfo == 10 ){ std::cout << "Mapping Descent DVs" << std::endl; }
    //! Associate decision vector values to mapped normalized specific energy levels within data maps.
    for ( unsigned int i = 0; i < E_mapped_Descent.size( ); ++i )
    {
        map_alpha_deg_Descent[ E_mapped_Descent( i ) ] = DV_mapped_Descent( 0 ) = alpha_deg_Descent( i );
        map_sigma_deg_Descent[ E_mapped_Descent( i ) ] = DV_mapped_Descent( 1 ) = sigma_deg_Descent( i );
        map_eps_T_deg_Descent[ E_mapped_Descent( i ) ] = DV_mapped_Descent( 2 ) = eps_T_deg_Descent( i );
        map_phi_T_deg_Descent[ E_mapped_Descent( i ) ] = DV_mapped_Descent( 3 ) = phi_T_deg_Descent( i );
        map_throttle_Descent[ E_mapped_Descent( i ) ]  = DV_mapped_Descent( 4 ) = throttle_Descent( i );
        DV_mapped_Descent( 5 ) = xn_Descent( i );
        map_DV_mapped_Descent[ E_mapped_Descent( i ) ] = DV_mapped_Descent;


        map_alpha_deg_Descent_LB[ E_mapped_Ascent( i ) ] = descentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( i , 0 );
        map_sigma_deg_Descent_LB[ E_mapped_Ascent( i ) ] = descentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( i , 0 );
        map_eps_T_deg_Descent_LB[ E_mapped_Ascent( i ) ] = descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( i , 0 );
        map_phi_T_deg_Descent_LB[ E_mapped_Ascent( i ) ] = descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( i , 0 );
        map_throttle_Descent_LB[ E_mapped_Ascent( i ) ]  = descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( i , 0 );

        map_alpha_deg_Descent_UB[ E_mapped_Ascent( i ) ] = descentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( i , 1 );
        map_sigma_deg_Descent_UB[ E_mapped_Ascent( i ) ] = descentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( i , 1 );
        map_eps_T_deg_Descent_UB[ E_mapped_Ascent( i ) ] = descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( i , 1 );
        map_phi_T_deg_Descent_UB[ E_mapped_Ascent( i ) ] = descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( i , 1 );
        map_throttle_Descent_UB[ E_mapped_Ascent( i ) ]  = descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( i , 1 );

    }


    if( debugInfo == 10 ){ std::cout << "Creating Optimization Interpolators' Settings" << std::endl; }

    std::shared_ptr< tudat::interpolators::InterpolatorSettings > interpolatorSettings = std::make_shared< tudat::interpolators::InterpolatorSettings >( tudat::interpolators::hermite_spline_interpolator );
    std::shared_ptr< tudat::interpolators::InterpolatorSettings > interpolatorBoundariesSettings = std::make_shared< tudat::interpolators::InterpolatorSettings >( tudat::interpolators::linear_interpolator );
    std::map< bislip::Parameters::Interpolators, std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > > Interpolators_Ascent, Interpolators_Descent;
    std::map< bislip::Parameters::Interpolators, std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > > Interpolators_Ascent_LB, Interpolators_Ascent_UB, Interpolators_Descent_LB, Interpolators_Descent_UB;

    if( debugInfo == 10 ){ std::cout << "Creating Ascent Interpolators" << std::endl; }

    Interpolators_Ascent[ bislip::Parameters::Interpolators::AngleOfAttack ]        = bislip::Variables::createOneDimensionalHermiteInterpolator( alpha_deg_Ascent, E_mapped_Ascent, map_alpha_deg_Ascent, interpolatorSettings );
    Interpolators_Ascent[ bislip::Parameters::Interpolators::BankAngle ]            = bislip::Variables::createOneDimensionalHermiteInterpolator( sigma_deg_Ascent, E_mapped_Ascent, map_sigma_deg_Ascent, interpolatorSettings );
    Interpolators_Ascent[ bislip::Parameters::Interpolators::ThrustElevationAngle ] = bislip::Variables::createOneDimensionalHermiteInterpolator( eps_T_deg_Ascent, E_mapped_Ascent, map_eps_T_deg_Ascent, interpolatorSettings );
    Interpolators_Ascent[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ]   = bislip::Variables::createOneDimensionalHermiteInterpolator( phi_T_deg_Ascent, E_mapped_Ascent, map_phi_T_deg_Ascent, interpolatorSettings );
    Interpolators_Ascent[ bislip::Parameters::Interpolators::ThrottleSetting ]      = bislip::Variables::createOneDimensionalHermiteInterpolator( throttle_Ascent, E_mapped_Ascent, map_throttle_Ascent, interpolatorSettings );

    if( debugInfo == 10 ){ std::cout << "Creating Ascent Interpolators Boundaries" << std::endl; }

    std::pair< double, double > angleOfAttackAscentLBInterpolatorBoundaryValues        = std::make_pair( ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( 0 , 0 ), ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( E_mapped_Ascent.size( ) - 1 , 0 ) );
    std::pair< double, double > bankAngleAscentLBInterpolatorBoundaryValues            = std::make_pair( ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( 0 , 0 ), ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( E_mapped_Ascent.size( ) - 1 , 0 ) );
    std::pair< double, double > thrustElevationAngleAscentLBInterpolatorBoundaryValues = std::make_pair( ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( 0 , 0 ), ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( E_mapped_Ascent.size( ) - 1 , 0 ) );
    std::pair< double, double > thrustAzimuthAngleAscentLBInterpolatorBoundaryValues   = std::make_pair( ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( 0 , 0 ), ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( E_mapped_Ascent.size( ) - 1 , 0 ) );
    std::pair< double, double > throttleSettingAscentLBInterpolatorBoundaryValues      = std::make_pair( ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( 0 , 0 ), ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( E_mapped_Ascent.size( ) - 1 , 0 ) );

    Interpolators_Ascent_LB[ bislip::Parameters::Interpolators::AngleOfAttack ]        = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_alpha_deg_Ascent_LB, interpolatorBoundariesSettings, angleOfAttackAscentLBInterpolatorBoundaryValues );
    Interpolators_Ascent_LB[ bislip::Parameters::Interpolators::BankAngle ]            = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_sigma_deg_Ascent_LB, interpolatorBoundariesSettings, bankAngleAscentLBInterpolatorBoundaryValues );
    Interpolators_Ascent_LB[ bislip::Parameters::Interpolators::ThrustElevationAngle ] = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_eps_T_deg_Ascent_LB, interpolatorBoundariesSettings, thrustElevationAngleAscentLBInterpolatorBoundaryValues );
    Interpolators_Ascent_LB[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ]   = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_phi_T_deg_Ascent_LB, interpolatorBoundariesSettings, thrustAzimuthAngleAscentLBInterpolatorBoundaryValues );
    Interpolators_Ascent_LB[ bislip::Parameters::Interpolators::ThrottleSetting ]      = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_throttle_Ascent_LB, interpolatorBoundariesSettings, throttleSettingAscentLBInterpolatorBoundaryValues );

    std::pair< double, double > angleOfAttackAscentUBInterpolatorBoundaryValues        = std::make_pair( ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( 0 , 1 ), ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( E_mapped_Ascent.size( ) - 1 , 1 ) );
    std::pair< double, double > bankAngleAscentUBInterpolatorBoundaryValues            = std::make_pair( ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( 0 , 1 ), ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( E_mapped_Ascent.size( ) - 1 , 1 ) );
    std::pair< double, double > thrustElevationAngleAscentUBInterpolatorBoundaryValues = std::make_pair( ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( 0 , 1 ), ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( E_mapped_Ascent.size( ) - 1 , 1 ) );
    std::pair< double, double > thrustAzimuthAngleAscentUBInterpolatorBoundaryValues   = std::make_pair( ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( 0 , 1 ), ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( E_mapped_Ascent.size( ) - 1 , 1 ) );
    std::pair< double, double > throttleSettingAscentUBInterpolatorBoundaryValues      = std::make_pair( ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( 0 , 1 ), ascentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( E_mapped_Ascent.size( ) - 1 , 1 ) );

    Interpolators_Ascent_UB[ bislip::Parameters::Interpolators::AngleOfAttack ]        = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_alpha_deg_Ascent_UB, interpolatorBoundariesSettings, angleOfAttackAscentUBInterpolatorBoundaryValues );
    Interpolators_Ascent_UB[ bislip::Parameters::Interpolators::BankAngle ]            = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_sigma_deg_Ascent_UB, interpolatorBoundariesSettings, bankAngleAscentUBInterpolatorBoundaryValues );
    Interpolators_Ascent_UB[ bislip::Parameters::Interpolators::ThrustElevationAngle ] = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_eps_T_deg_Ascent_UB, interpolatorBoundariesSettings, thrustElevationAngleAscentUBInterpolatorBoundaryValues );
    Interpolators_Ascent_UB[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ]   = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_phi_T_deg_Ascent_UB, interpolatorBoundariesSettings, thrustAzimuthAngleAscentUBInterpolatorBoundaryValues );
    Interpolators_Ascent_UB[ bislip::Parameters::Interpolators::ThrottleSetting ]      = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_throttle_Ascent_UB, interpolatorBoundariesSettings, throttleSettingAscentUBInterpolatorBoundaryValues );


    if( debugInfo == 10 ){ std::cout << "Creating Descent Interpolators" << std::endl; }
    //! Declare and initialize interpolators for Descent phase.

    if( debugInfo == 10 ){ std::cout << "   Angle of Attack" << std::endl; }
    Interpolators_Descent[ bislip::Parameters::Interpolators::AngleOfAttack ]        = bislip::Variables::createOneDimensionalHermiteInterpolator( alpha_deg_Descent, E_mapped_Descent, map_alpha_deg_Descent, interpolatorSettings );
    if( debugInfo == 10 ){ std::cout << "   Bank Angle" << std::endl; }
    Interpolators_Descent[ bislip::Parameters::Interpolators::BankAngle ]            = bislip::Variables::createOneDimensionalHermiteInterpolator( sigma_deg_Descent, E_mapped_Descent, map_sigma_deg_Descent, interpolatorSettings );
    if( debugInfo == 10 ){ std::cout << "   Thrust Elevation Angle" << std::endl; }
    Interpolators_Descent[ bislip::Parameters::Interpolators::ThrustElevationAngle ] = bislip::Variables::createOneDimensionalHermiteInterpolator( eps_T_deg_Descent, E_mapped_Descent, map_eps_T_deg_Descent, interpolatorSettings );
    if( debugInfo == 10 ){ std::cout << "   Thrust Azimuth Angle" << std::endl; }
    Interpolators_Descent[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ]   = bislip::Variables::createOneDimensionalHermiteInterpolator( phi_T_deg_Descent, E_mapped_Descent, map_phi_T_deg_Descent, interpolatorSettings );
    if( debugInfo == 10 ){ std::cout << "   Throttle Setting" << std::endl; }
    Interpolators_Descent[ bislip::Parameters::Interpolators::ThrottleSetting ]      = bislip::Variables::createOneDimensionalHermiteInterpolator( throttle_Descent, E_mapped_Descent, map_throttle_Descent, interpolatorSettings );


    if( debugInfo == 10 ){ std::cout << "Creating Descent Interpolators Boundaries" << std::endl; }

    std::pair< double, double > angleOfAttackDescentLBInterpolatorBoundaryValues        = std::make_pair( descentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( 0 , 0 ), descentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( E_mapped_Descent.size( ) - 1 , 0 ) );
    std::pair< double, double > bankAngleDescentLBInterpolatorBoundaryValues            = std::make_pair( descentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( 0 , 0 ), descentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( E_mapped_Descent.size( ) - 1 , 0 ) );
    std::pair< double, double > thrustElevationAngleDescentLBInterpolatorBoundaryValues = std::make_pair( descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( 0 , 0 ), descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( E_mapped_Descent.size( ) - 1 , 0 ) );
    std::pair< double, double > thrustAzimuthAngleDescentLBInterpolatorBoundaryValues   = std::make_pair( descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( 0 , 0 ), descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( E_mapped_Descent.size( ) - 1 , 0 ) );
    std::pair< double, double > throttleSettingDescentLBInterpolatorBoundaryValues      = std::make_pair( descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( 0 , 0 ), descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( E_mapped_Descent.size( ) - 1 , 0 ) );

    Interpolators_Descent_LB[ bislip::Parameters::Interpolators::AngleOfAttack ]        = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_alpha_deg_Descent_LB, interpolatorBoundariesSettings, angleOfAttackDescentLBInterpolatorBoundaryValues );
    Interpolators_Descent_LB[ bislip::Parameters::Interpolators::BankAngle ]            = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_sigma_deg_Descent_LB, interpolatorBoundariesSettings, bankAngleDescentLBInterpolatorBoundaryValues );
    Interpolators_Descent_LB[ bislip::Parameters::Interpolators::ThrustElevationAngle ] = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_eps_T_deg_Descent_LB, interpolatorBoundariesSettings, thrustElevationAngleDescentLBInterpolatorBoundaryValues );
    Interpolators_Descent_LB[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ]   = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_phi_T_deg_Descent_LB, interpolatorBoundariesSettings, thrustAzimuthAngleDescentLBInterpolatorBoundaryValues );
    Interpolators_Descent_LB[ bislip::Parameters::Interpolators::ThrottleSetting ]      = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_throttle_Descent_LB, interpolatorBoundariesSettings, throttleSettingDescentLBInterpolatorBoundaryValues );

    std::pair< double, double > angleOfAttackDescentUBInterpolatorBoundaryValues        = std::make_pair( descentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( 0 , 1 ), descentParameterBoundsMap.at( bislip::Parameters::Interpolators::AngleOfAttack ).coeff( E_mapped_Descent.size( ) - 1 , 1 ) );
    std::pair< double, double > bankAngleDescentUBInterpolatorBoundaryValues            = std::make_pair( descentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( 0 , 1 ), descentParameterBoundsMap.at( bislip::Parameters::Interpolators::BankAngle ).coeff( E_mapped_Descent.size( ) - 1 , 1 ) );
    std::pair< double, double > thrustElevationAngleDescentUBInterpolatorBoundaryValues = std::make_pair( descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( 0 , 1 ), descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustElevationAngle ).coeff( E_mapped_Descent.size( ) - 1 , 1 ) );
    std::pair< double, double > thrustAzimuthAngleDescentUBInterpolatorBoundaryValues   = std::make_pair( descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( 0 , 1 ), descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrustAzimuthAngle ).coeff( E_mapped_Descent.size( ) - 1 , 1 ) );
    std::pair< double, double > throttleSettingDescentUBInterpolatorBoundaryValues      = std::make_pair( descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( 0 , 1 ), descentParameterBoundsMap.at( bislip::Parameters::Interpolators::ThrottleSetting ).coeff( E_mapped_Descent.size( ) - 1 , 1 ) );

    Interpolators_Descent_UB[ bislip::Parameters::Interpolators::AngleOfAttack ]        = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_alpha_deg_Descent_UB, interpolatorBoundariesSettings, angleOfAttackDescentUBInterpolatorBoundaryValues );
    Interpolators_Descent_UB[ bislip::Parameters::Interpolators::BankAngle ]            = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_sigma_deg_Descent_UB, interpolatorBoundariesSettings, bankAngleDescentUBInterpolatorBoundaryValues );
    Interpolators_Descent_UB[ bislip::Parameters::Interpolators::ThrustElevationAngle ] = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_eps_T_deg_Descent_UB, interpolatorBoundariesSettings, thrustElevationAngleDescentUBInterpolatorBoundaryValues );
    Interpolators_Descent_UB[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ]   = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_phi_T_deg_Descent_UB, interpolatorBoundariesSettings, thrustAzimuthAngleDescentUBInterpolatorBoundaryValues );
    Interpolators_Descent_UB[ bislip::Parameters::Interpolators::ThrottleSetting ]      = tudat::interpolators::createOneDimensionalInterpolator< double, double >( map_throttle_Descent_UB, interpolatorBoundariesSettings, throttleSettingDescentUBInterpolatorBoundaryValues );

    //! Declare vectors containing interpolated values.
    Eigen::VectorXd interpolated_values_Ascent( 5 ), interpolated_values_Ascent_LB( 5 ), interpolated_values_Ascent_UB( 5 );
    Eigen::VectorXd interpolated_values_Descent( 5 ), interpolated_values_Descent_LB( 5 ), interpolated_values_Descent_UB( 5 );

    //! Declare data map to contain vectors of interpolated values.
    std::map< double, Eigen::VectorXd > evaluatedInterpolatorsAscent, evaluatedInterpolatorsAscent_LB, evaluatedInterpolatorsAscent_UB;
    std::map< double, Eigen::VectorXd > evaluatedInterpolatorsDescent, evaluatedInterpolatorsDescent_LB, evaluatedInterpolatorsDescent_UB;

    //! Declare evaluation variable.
    double eval;

    if( debugInfo == 10 ){ std::cout << "Evaluate Interpolators" << std::endl; }
    //! Loop to populate vectors of interpolated values and then pass to data map.
    //!     Number of evaluations has been arbitrarily selected.
    double pp = 0;
    for ( unsigned int i = 0; i < 1001; ++i )
    {
        eval = pp * E_mapped_Ascent.maxCoeff() / 1000;
        interpolated_values_Ascent( 0 ) = ( Interpolators_Ascent[ bislip::Parameters::Interpolators::AngleOfAttack ] )->interpolate( eval );
        interpolated_values_Ascent( 1 ) = ( Interpolators_Ascent[ bislip::Parameters::Interpolators::BankAngle ] )->interpolate( eval );
        interpolated_values_Ascent( 2 ) = ( Interpolators_Ascent[ bislip::Parameters::Interpolators::ThrustElevationAngle ] )->interpolate( eval );
        interpolated_values_Ascent( 3 ) = ( Interpolators_Ascent[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ] )->interpolate( eval );
        interpolated_values_Ascent( 4 ) = ( Interpolators_Ascent[ bislip::Parameters::Interpolators::ThrottleSetting ] )->interpolate( eval );

        evaluatedInterpolatorsAscent[ eval ] = interpolated_values_Ascent;

        interpolated_values_Ascent_LB( 0 ) = ( Interpolators_Ascent_LB[ bislip::Parameters::Interpolators::AngleOfAttack ] )->interpolate( eval );
        interpolated_values_Ascent_LB( 1 ) = ( Interpolators_Ascent_LB[ bislip::Parameters::Interpolators::BankAngle ] )->interpolate( eval );
        interpolated_values_Ascent_LB( 2 ) = ( Interpolators_Ascent_LB[ bislip::Parameters::Interpolators::ThrustElevationAngle ] )->interpolate( eval );
        interpolated_values_Ascent_LB( 3 ) = ( Interpolators_Ascent_LB[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ] )->interpolate( eval );
        interpolated_values_Ascent_LB( 4 ) = ( Interpolators_Ascent_LB[ bislip::Parameters::Interpolators::ThrottleSetting ] )->interpolate( eval );

        evaluatedInterpolatorsAscent_LB[ eval ] = interpolated_values_Ascent_LB;

        interpolated_values_Ascent_UB( 0 ) = ( Interpolators_Ascent_UB[ bislip::Parameters::Interpolators::AngleOfAttack ] )->interpolate( eval );
        interpolated_values_Ascent_UB( 1 ) = ( Interpolators_Ascent_UB[ bislip::Parameters::Interpolators::BankAngle ] )->interpolate( eval );
        interpolated_values_Ascent_UB( 2 ) = ( Interpolators_Ascent_UB[ bislip::Parameters::Interpolators::ThrustElevationAngle ] )->interpolate( eval );
        interpolated_values_Ascent_UB( 3 ) = ( Interpolators_Ascent_UB[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ] )->interpolate( eval );
        interpolated_values_Ascent_UB( 4 ) = ( Interpolators_Ascent_UB[ bislip::Parameters::Interpolators::ThrottleSetting ] )->interpolate( eval );

        evaluatedInterpolatorsAscent_UB[ eval ] = interpolated_values_Ascent_UB;

        eval = pp * E_mapped_Descent.maxCoeff() / 1000;
        interpolated_values_Descent( 0 ) = ( Interpolators_Descent[ bislip::Parameters::Interpolators::AngleOfAttack ] )->interpolate( eval );
        interpolated_values_Descent( 1 ) = ( Interpolators_Descent[ bislip::Parameters::Interpolators::BankAngle ] )->interpolate( eval );
        interpolated_values_Descent( 2 ) = ( Interpolators_Descent[ bislip::Parameters::Interpolators::ThrustElevationAngle ] )->interpolate( eval );
        interpolated_values_Descent( 3 ) = ( Interpolators_Descent[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ] )->interpolate( eval );
        interpolated_values_Descent( 4 ) = ( Interpolators_Descent[ bislip::Parameters::Interpolators::ThrottleSetting ] )->interpolate( eval );

        evaluatedInterpolatorsDescent[ eval ] = interpolated_values_Descent;

        interpolated_values_Descent_LB( 0 ) = ( Interpolators_Descent_LB[ bislip::Parameters::Interpolators::AngleOfAttack ] )->interpolate( eval );
        interpolated_values_Descent_LB( 1 ) = ( Interpolators_Descent_LB[ bislip::Parameters::Interpolators::BankAngle ] )->interpolate( eval );
        interpolated_values_Descent_LB( 2 ) = ( Interpolators_Descent_LB[ bislip::Parameters::Interpolators::ThrustElevationAngle ] )->interpolate( eval );
        interpolated_values_Descent_LB( 3 ) = ( Interpolators_Descent_LB[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ] )->interpolate( eval );
        interpolated_values_Descent_LB( 4 ) = ( Interpolators_Descent_LB[ bislip::Parameters::Interpolators::ThrottleSetting ] )->interpolate( eval );

        evaluatedInterpolatorsDescent_LB[ eval ] = interpolated_values_Descent_LB;

        interpolated_values_Descent_UB( 0 ) = ( Interpolators_Descent_UB[ bislip::Parameters::Interpolators::AngleOfAttack ] )->interpolate( eval );
        interpolated_values_Descent_UB( 1 ) = ( Interpolators_Descent_UB[ bislip::Parameters::Interpolators::BankAngle ] )->interpolate( eval );
        interpolated_values_Descent_UB( 2 ) = ( Interpolators_Descent_UB[ bislip::Parameters::Interpolators::ThrustElevationAngle ] )->interpolate( eval );
        interpolated_values_Descent_UB( 3 ) = ( Interpolators_Descent_UB[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ] )->interpolate( eval );
        interpolated_values_Descent_UB( 4 ) = ( Interpolators_Descent_UB[ bislip::Parameters::Interpolators::ThrottleSetting ] )->interpolate( eval );

        evaluatedInterpolatorsDescent_UB[ eval ] = interpolated_values_Descent_UB;

        pp += 1;
    }

    if( debugInfo == 10 ){ std::cout << "Interpolators Evaluated." << std::endl; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           ASSIGN INTERPOLATORS & BOUNDS: ASCENT           /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){ std::cout << "Pass Ascent Interpolators to vehicle systems" << std::endl;}

    //! Pass Ascent phase interpolators to vehicle systems.
    bislipSystems->setParameterInterpolator( Interpolators_Ascent );

    if( debugInfo == 10 ){ std::cout << "Pass Ascent Bounds to vehicle systems" << std::endl;}

    //! Pass parameter bounds to vehicle systems.
    //bislipSystems->setParameterBounds( problemInput->getAscentParameterBoundsMap() );
    bislipSystems->setParameterLowerBounds( Interpolators_Ascent_LB);
    bislipSystems->setParameterUpperBounds( Interpolators_Ascent_UB );


    double simulationEndEpoch_Ascent = simulationStartEpoch;
    Eigen::VectorXd dependentVariableFinalState_Ascent;
    double finalAirspeed_Ascent = initialAirspeed_Ascent;
    double finalMass_Ascent = initialMass_Ascent;
    double initialMass_Descent = initialMass_Ascent;
    double finalSpecificEnergy_Ascent = E_i;

    std::map< double, Eigen::VectorXd > propagationTimeHistoryMap_Ascent;
    std::map< double, Eigen::VectorXd > depVarTimeHistoryMap_Ascent;

    //if( E_i < 2 * E_f )
    // {

    //if( debugInfo == 10 ){ std::cout << "Initial Specific Energy is LOWER than Final Specific Energy" << std::endl; }

    //( bislipSystems->getParameterInterpolator( bislip::Parameters::Optimization::ThrustElevationAngle ) )->interpolate( bislip::Variables::computeNormalizedSpecificEnergy( currentAltitude, currentAirspeed, bislipSystems->getMaximumSpecificEnergy() ) );



    if( debugInfo == 10 ){ std::cout << "Setting various initial values." << std::endl; }

    if( debugInfo == 10 ){std::cout << "    Evaluating Angle of Attack Interpolator" << std::endl; }
    if( bislipSystems->getValidationFlag() == true )
    {
        if( debugInfo == 10 ){std::cout << "         Validation Interpolator" << std::endl; }
        bislipSystems->setCurrentAngleOfAttack( tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getKourouAngleOfAttackInterpolator( ) )->interpolate( 0 ) ) );
    }
    else
    {
        if( debugInfo == 10 ){std::cout << "         Optimization Interpolator" << std::endl; }
        bislipSystems->setCurrentAngleOfAttack( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::AngleOfAttack, bodyMap, vehicleName, centralBodyName ) ) );
    }
    if( debugInfo == 10 ){std::cout << "         Initial Angle of Attack = " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentAngleOfAttack() ) << std::endl; }
    /*
    if( debugInfo == 10 ){std::cout << "    Hardcoding Bank Angle" << std::endl; }
    bislipSystems->setCurrentBankAngle( 0.0 );
    if( debugInfo == 10 ){std::cout << "         Initial Bank Angle = " << bislipSystems->getCurrentBankAngle() << std::endl; }

    if( debugInfo == 10 ){std::cout << "    Repeating Flight-Path Angle" << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Initial Flight-Path Angle = " << bislipSystems->getInitialFlightPathAngle() << std::endl; }

    if( debugInfo == 10 ){std::cout << "    Repeating Heading Angle" << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Initial Heading Angle = " << bislipSystems->getInitialHeading() << std::endl; }

    if( debugInfo == 10 ){std::cout << "    Evaluating Thrust Azimuth Angle Interpolator" << std::endl; }
    bislipSystems->setCurrentThrustAzimuthAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::ThrustAzimuthAngle, bodyMap, vehicleName ) ) );
    if( debugInfo == 10 ){std::cout << "         Initial Thrust Azimuth Angle   = " << bislipSystems->getCurrentThrustAzimuthAngle() << std::endl; }

    if( debugInfo == 10 ){std::cout << "    Repeating Mass Variables" << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Current mass = " << bodyMap.at( vehicleName )->getBodyMass() << std::endl; }
    if( debugInfo == 10 ){std::cout << "         Dry mass     = " << vehicleSystems->getDryMass() << std::endl; }

    if( debugInfo == 10 ){std::cout << "    Evaluating Engine Status Function" << std::endl; }
    bislipSystems->setCurrentEngineStatus( bislip::Variables::determineEngineStatus( bodyMap.at( vehicleName )->getBodyMass(), vehicleSystems->getDryMass() ) );
    if( debugInfo == 10 ){std::cout << "         Initial Engine Status    = " << bislipSystems->getCurrentEngineStatus() << std::endl; }

    if( debugInfo == 10 ){std::cout << "    Evaluating Throttle Setting Interpolator" << std::endl; }
    bislipSystems->setCurrentThrottleSetting( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::ThrottleSetting, bodyMap, vehicleName ) );
    if( debugInfo == 10 ){std::cout << "         Initial Throttle Setting = " << bislipSystems->getCurrentThrottleSetting() << std::endl; }

    if( debugInfo == 10 ){std::cout << "    Evaluating Thrust Magnitude Function" << std::endl; }
    bislipSystems->setCurrentThrustMagnitude( bislip::Variables::computeThrustMagnitude( bodyMap, vehicleName ) );
    if( debugInfo == 10 ){std::cout << "         Initial Thrust Magnitude = " << bislipSystems->getCurrentThrustMagnitude() << std::endl; }

    if( debugInfo == 10 ){std::cout << "    Evaluating Thrust Elevation Angle Interpolator" << std::endl; }
    bislipSystems->setCurrentThrustElevationAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::ThrustElevationAngle, bodyMap, vehicleName ) ) );
    if( debugInfo == 10 ){std::cout << "         Initial Thrust Elevation Angle = " << bislipSystems->getCurrentThrustElevationAngle() << std::endl; }

    if( debugInfo == 10 ){std::cout << "    Evaluating Body-Fixed Thrust Direction Function" << std::endl; }
    bislipSystems->setCurrentBodyFixedThrustDirection( bislip::Variables::computeBodyFixedThrustDirection( bodyMap, vehicleName ) );
    if( debugInfo == 10 ){std::cout << "         Initial Thrust Direction = [ " << ( bislipSystems->getCurrentBodyFixedThrustDirection() )( 0 ) << " , " << ( bislipSystems->getCurrentBodyFixedThrustDirection() )( 1 ) << " , " << ( bislipSystems->getCurrentBodyFixedThrustDirection() )( 2 ) << " ]" << std::endl; }
*/
    bislipSystems->setCurrentThrustElevationAngle( 0.0 );
    bislipSystems->setCurrentThrustAzimuthAngle( 0.0 );
    bislipSystems->setCurrentBankAngle( 0.0 );

    if( debugInfo == 10 ){ std::cout << "    Setting Initial Control Surface Deflections" << std::endl; }
    vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", bislipSystems->getCurrentBodyFlapAngle() );
    vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonLeft", bislipSystems->getCurrentElevonAngle() );
    vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonRight", bislipSystems->getCurrentElevonAngle() );
    if( debugInfo == 10 ){ std::cout << "         Initial BodyFlap Deflection = " << vehicleSystems->getCurrentControlSurfaceDeflection( "BodyFlap" ) << std::endl; }
    if( debugInfo == 10 ){ std::cout << "         Initial Elevon Deflection   = " << vehicleSystems->getCurrentControlSurfaceDeflection( "ElevonLeft" ) << std::endl; }

    if( debugInfo == 10 ){ std::cout << "    Calculate Full Current Coefficients" << std::endl; }
    bislipSystems->setFullCurrentCoefficients( bislip::Variables::computeFullCurrentCoefficients( bodyMap, vehicleName ) );
    if( debugInfo == 10 ){ std::cout << "         Preliminary Coefficient Vector = [ " << bislipSystems->getFullCurrentCoefficients()( 0 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 1 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 2 ) << bislipSystems->getFullCurrentCoefficients()( 3 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 4 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 5 ) << " ]" << std::endl; }

    if( debugInfo == 10 ){ std::cout << "    Calculate Preliminary Aerodynamic Frame Aerodyamic Load" << std::endl; }
    Eigen::Vector3d aerodynamicFrameAerodynamicLoadVector = bislip::Variables::computeAerodynamicFrameAerodynamicLoad( bodyMap, vehicleName );
    if( debugInfo == 10 ){ std::cout << "         Preliminary Aerodyamic Load = [ " << aerodynamicFrameAerodynamicLoadVector( 0 ) << ", " << aerodynamicFrameAerodynamicLoadVector( 1 ) << ", " << aerodynamicFrameAerodynamicLoadVector( 2 ) << " ]" << std::endl; }

    if( debugInfo == 10 ){ std::cout << "    Set Preliminary Drag Force" << std::endl; }
    bislipSystems->setCurrentDragForce( -aerodynamicFrameAerodynamicLoadVector( 0 ) );
    if( debugInfo == 10 ){ std::cout << "         Preliminary Drag Force = " << bislipSystems->getCurrentDragForce() << std::endl; }

    if( debugInfo == 10 ){ std::cout << "    Set Preliminary Lift Force" << std::endl; }
    bislipSystems->setCurrentLiftForce( -aerodynamicFrameAerodynamicLoadVector( 2 ) );
    if( debugInfo == 10 ){ std::cout << "         Preliminary Lift Force = " << bislipSystems->getCurrentLiftForce() << std::endl; }

    if( debugInfo == 10 ){ std::cout << "    Evaluating Local Gravity Vector Function" << std::endl; }
    bislipSystems->setCurrentLocalGravityVector( bislip::Variables::computeLocalGravity( bodyMap, vehicleName, centralBodyName ) );
    if( debugInfo == 10 ){ std::cout << "         Initial Local Gravity Vector = [ " << bislipSystems->getCurrentLocalGravityVector()( 0 ) << ", " << bislipSystems->getCurrentLocalGravityVector()( 1 ) << " ]" << std::endl; }

    if( debugInfo == 10 ){ std::cout << "    Create caller to Evaluate Guidance Functions" << std::endl; }
    bislip::MyGuidance initialEvaluation( bodyMap, vehicleName, centralBodyName, simulationStartEpoch );

    if( debugInfo == 10 ){ std::cout << "    Evaluating All Guidance Functions to re-Determine All Relevant Initial Values" << std::endl; }
    if( bislipSystems->getValidationFlag() == true )
    { initialEvaluation.evaluateValidationGuidanceFunctions( bislipSystems, vehicleSystems, "Descent", 0.0 ); }
    else { initialEvaluation.evaluateGuidanceFunctions( bislipSystems, vehicleSystems, "Ascent", 0.0 ); }

    if( debugInfo == 10 ){ std::cout << "    Calculate Full Current Coefficients" << std::endl; }
    bislipSystems->setFullCurrentCoefficients( bislip::Variables::computeFullCurrentCoefficients( bodyMap, vehicleName ) );
    if( debugInfo == 10 ){ std::cout << "         Initial Coefficient Vector = [ " << bislipSystems->getFullCurrentCoefficients()( 0 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 1 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 2 ) << bislipSystems->getFullCurrentCoefficients()( 3 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 4 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 5 ) << " ]" << std::endl; }

    if( debugInfo == 10 ){ std::cout << "    Calculate Initial Aerodynamic Frame Aerodyamic Load" << std::endl; }
    aerodynamicFrameAerodynamicLoadVector = bislip::Variables::computeAerodynamicFrameAerodynamicLoad( bodyMap, vehicleName );
    if( debugInfo == 10 ){ std::cout << "         Initial Aerodyamic Load = [ " << aerodynamicFrameAerodynamicLoadVector( 0 ) << ", " << aerodynamicFrameAerodynamicLoadVector( 1 ) << ", " << aerodynamicFrameAerodynamicLoadVector( 2 ) << " ]" << std::endl; }

    if( debugInfo == 10 ){ std::cout << "    Set Initial Drag Force" << std::endl; }
    bislipSystems->setCurrentDragForce( -aerodynamicFrameAerodynamicLoadVector( 0 ) );
    if( debugInfo == 10 ){std::cout << "         Initial Drag Force = " << bislipSystems->getCurrentDragForce() << std::endl; }

    if( debugInfo == 10 ){ std::cout << "    Set Initial Lift Force" << std::endl; }
    bislipSystems->setCurrentLiftForce( -aerodynamicFrameAerodynamicLoadVector( 2 ) );
    if( debugInfo == 10 ){std::cout << "         Initial Lift Force = " << bislipSystems->getCurrentLiftForce() << std::endl; }

    if( debugInfo == 10 ){ std::cout << "    Set Initial Flight-Path Angle Rate" << std::endl; }
    bislipSystems->setCurrentFlightPathAngleRate( bislip::Variables::computeFlightPathAngleRate( bodyMap, vehicleName, centralBodyName ) );
    if( debugInfo == 10 ){std::cout << "         Initial Flight-Path Angle Rate = " << bislipSystems->getCurrentFlightPathAngleRate() << std::endl; }

    std::map< bislip::BislipVehicleSystems::PreviousConditions, double > previousConditions;
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::height ]                 = bislipSystems->getInitialHeight( );
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::density ]                = bislipSystems->getInitialDensity( );
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::airspeed ]               = bislipSystems->getInitialAirspeed( );
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::angle_of_attack ]        = bislipSystems->getCurrentAngleOfAttack( );
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::thrust_elevation_angle ] = bislipSystems->getCurrentThrustElevationAngle( );
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::thrust_azimuth_angle ]   = bislipSystems->getCurrentThrustAzimuthAngle( );
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::thrust_magnitude ]       = bislipSystems->getCurrentThrustMagnitude( );
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::drag_coefficient ]       = bislipSystems->getFullCurrentCoefficients()( 0 );
    previousConditions[ bislip::BislipVehicleSystems::PreviousConditions::lift_coefficient ]       = bislipSystems->getFullCurrentCoefficients()( 2 );

    bislipSystems->setInitialValueFlag( false );
    bislipSystems->setPreviousConditions( previousConditions );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE MASS RATE SETTINGS: ASCENT            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){std::cout << "Create mass rate models: Ascent" << std::endl;}

    //! Declare and initialize mass rate model settings.
    std::shared_ptr< MassRateModelSettings > massRateModelSettings =
            std::make_shared< FromThrustMassModelSettings >( true );

    //! Declare and initialize mass rate model settings.
    std::map< std::string, std::shared_ptr< MassRateModel > > massRateModels;
    massRateModels[  vehicleName ] = createMassRateModel(
                vehicleName, massRateModelSettings, bodyMap, problemInput->getAccelerationMap() );

    //! Create settings for propagating the mass of the vehicle.
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back( vehicleName );

    //! Declare and initialize starting mass of vehicle.
    Eigen::VectorXd initialBodyMasses_Ascent( 1 );
    initialBodyMasses_Ascent( 0 ) = initialMass_Ascent + additionalMass;

    //! Declare and initialize mass propagation settings.
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings_Ascent =
            std::make_shared< MassPropagatorSettings< double > >(
                bodiesWithMassToPropagate, massRateModels, initialBodyMasses_Ascent, problemInput->getAscentTerminationSettings(), problemInput->getDependentVariablesToSave() );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE PROPAGATION SETTINGS: ASCENT              ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){std::cout << "Create Propagation Settings: Ascent" << std::endl;}

    //! Create translational propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings_Ascent =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( problemInput->getCentralBodies(),
              problemInput->getAccelerationMap(),
              problemInput->getBodiesToIntegrate(),
              systemInitialState_Ascent,
              problemInput->getAscentTerminationSettings(),
              cowell,
              problemInput->getDependentVariablesToSave() );

    //! Declare and initialize list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector_Ascent;
    propagatorSettingsVector_Ascent.push_back( translationalPropagatorSettings_Ascent );
    propagatorSettingsVector_Ascent.push_back( massPropagatorSettings_Ascent );

    //! Declare and initialize propagation settings for both mass and translational dynamics.
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings_Ascent =
            std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector_Ascent, problemInput->getAscentTerminationSettings(), problemInput->getDependentVariablesToSave() );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE INTEGRATION SETTINGS: ASCENT           ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){ std::cout << "Create Integration Settings: Ascent" << std::endl;}

    //! Declare and initialize integrator settings.
    std::shared_ptr< IntegratorSettings<  > > integratorSettings_Ascent =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationStartEpoch,
              propagationStepSize, 1, false );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          PROPAGATE TRAJECTORY: ASCENT         /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){ std::cout << "Starting Ascent propagation" << std::endl; }

    //! Propagate trajectory.
    SingleArcDynamicsSimulator< double > ascentSimulation(
                bodyMap,
                integratorSettings_Ascent,
                propagatorSettings_Ascent );

    if( debugInfo == 10 ){ std::cout << "Ascent Propagation done" << std::endl; }

    if( debugInfo == 10 ){ std::cout << "Retrieve ASCENT Propagation and Dependent Variable Dependent Variable Time History maps" << std::endl; }
    propagationTimeHistoryMap_Ascent = ascentSimulation.getEquationsOfMotionNumericalSolution( );
    depVarTimeHistoryMap_Ascent      = ascentSimulation.getDependentVariableHistory( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////          RETRIEVE ASCENT-DESCENT PROPAGATION LINKING DATA          //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Retrieve final Ascent epoch
    simulationEndEpoch_Ascent = ( ascentSimulation.getEquationsOfMotionNumericalSolution().rbegin() )->first;

    if( debugInfo == 10 ){ std::cout << "Ascent Phase Initial Epoch = " << simulationStartEpoch << std::endl; }
    if( debugInfo == 10 ){ std::cout << "Ascent Phase Final Epoch   = " << simulationEndEpoch_Ascent << std::endl; }

    //! Retrieve final Ascent state, equal to initial Descent state
    Eigen::Vector6d systemInitialState_Descent = ( ( ascentSimulation.getEquationsOfMotionNumericalSolution().rbegin() )->second ).segment( 0, 6 );

    if( debugInfo == 10 ){ std::cout << "Ascent Phase Final State = [ " << systemInitialState_Descent( 0 ) << " , ";
        std::cout << systemInitialState_Descent( 1 ) << " , " << systemInitialState_Descent( 2 ) << " , " ;
        std::cout << systemInitialState_Descent( 3 ) << " , " << systemInitialState_Descent( 4 ) << " , " ;
        std::cout << systemInitialState_Descent( 5 ) << " ]" <<std::endl; }


    //! Retrieve dependent variables of final Ascent state
    dependentVariableFinalState_Ascent = ( ascentSimulation.getDependentVariableHistory( ).rbegin() )->second;

    //if( debugInfo == 10 ){ std::cout << "Ascent Phase Final DepVar = " << dependentVariableFinalState_Ascent << std::endl; }

    //! Retrieve final Ascent mass, equal to initial Descent mass
    finalAirspeed_Ascent       = dependentVariableFinalState_Ascent[ 13 ];
    finalMass_Ascent           = dependentVariableFinalState_Ascent[ 27 ];
    initialMass_Descent        = finalMass_Ascent;
    finalSpecificEnergy_Ascent = dependentVariableFinalState_Ascent[ 28 ];

    //  }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            COMPLETE KNOWN DESCENT STATES            ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( debugInfo == 10 ){ std::cout << "Complete Known Descent States" << std::endl; }
    if( debugInfo == 10 ){ std::cout << "   Nothing done here, as the initial Descent State is the final Ascent State" << std::endl; }

    /*
    if( E_i > 2 * E_f )
    {

        if( debugInfo == 10 ){ std::cout << "Initial Specific Energy is HIGHER than Final Specific Energy" << std::endl; }

        bislipSystems->setCurrentThrustElevationAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustElevationAngle, bodyMap, vehicleName ) ) );
        if( debugInfo == 10 ){std::cout << "     Initial Thrust Elevation Angel = " << bislipSystems->getCurrentThrustElevationAngle() << std::endl; }

        bislipSystems->setCurrentThrustAzimuthAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustAzimuthAngle, bodyMap, vehicleName ) ) );
        if( debugInfo == 10 ){std::cout << "     Initial Thrust Azimuth Angel   = " << bislipSystems->getCurrentThrustAzimuthAngle() << std::endl; }

        bislipSystems->setCurrentThrottleSetting( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrottleSetting, bodyMap, vehicleName ) );
        if( debugInfo == 10 ){std::cout << "     Initial Throttle Setting       = " << bislipSystems->getCurrentThrottleSetting() << std::endl; }

        if( debugInfo == 10 ){std::cout << "     Current mass = " << bodyMap.at( vehicleName )->getBodyMass() << std::endl; }
        if( debugInfo == 10 ){std::cout << "     Dry mass     = " << vehicleSystems->getDryMass() << std::endl; }

        bislipSystems->setCurrentEngineStatus( bislip::Variables::determineEngineStatus( bodyMap.at( vehicleName )->getBodyMass(), vehicleSystems->getDryMass() ) );
        if( debugInfo == 10 ){std::cout << "     Initial Engine Status    = " << bislipSystems->getCurrentEngineStatus() << std::endl; }

        bislipSystems->setCurrentBodyFixedThrustDirection( bislip::Variables::computeBodyFixedThrustDirection( bodyMap, vehicleName ) );
        if( debugInfo == 10 ){std::cout << "     Initial Thrust Direction = " << bislipSystems->getCurrentBodyFixedThrustDirection() << std::endl; }

        bislipSystems->setCurrentThrustMagnitude( bislip::Variables::computeThrustMagnitude( bodyMap, vehicleName ) );
        if( debugInfo == 10 ){std::cout << "     Initial Thrust Magnitude = " << bislipSystems->getCurrentThrustMagnitude() << std::endl; }



        systemInitialState_Descent = systemInitialState_Ascent;
        initialMass_Descent = initialMass_Ascent + additionalMass;
    }


*/





    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE PROPAGATION SETTINGS: DESCENT              ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){ std::cout << "Create Propagation Settings - Descent" << std::endl; }

    if( debugInfo == 10 ){ std::cout << "   Create Translational Propagation Settings" << std::endl; }
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings_Descent =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( problemInput->getCentralBodies(),
              problemInput->getAccelerationMap(),
              problemInput->getBodiesToIntegrate(),
              systemInitialState_Descent,
              problemInput->getDescentTerminationSettings(),
              cowell,
              problemInput->getDependentVariablesToSave() );

    if( debugInfo == 10 ){ std::cout << "Initial Ascent Mass  = " << initialMass_Ascent + additionalMass << std::endl; }
    if( debugInfo == 10 ){ std::cout << "Initial Descent Mass = " << initialMass_Descent << std::endl; }

    if( debugInfo == 10 ){ std::cout << "       Declare and initialize starting mass of vehicle for Descent phase" << std::endl; }
    Eigen::VectorXd initialBodyMasses_Descent( 1 );
    initialBodyMasses_Descent( 0 ) = initialMass_Descent;

    if( debugInfo == 10 ){ std::cout << "   Create Mass Propagation Settings" << std::endl; }
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings_Descent =
            std::make_shared< MassPropagatorSettings< double > >(
                bodiesWithMassToPropagate, massRateModels, initialBodyMasses_Descent, problemInput->getDescentTerminationSettings(), problemInput->getDependentVariablesToSave() );

    if( debugInfo == 10 ){ std::cout << "   Create List of Propagation Settings" << std::endl; }
    //! Declare and initialize list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector_Descent;
    propagatorSettingsVector_Descent.push_back( translationalPropagatorSettings_Descent );
    propagatorSettingsVector_Descent.push_back( massPropagatorSettings_Descent );

    if( debugInfo == 10 ){ std::cout << "   Create Propagation Settings (for the list created above)" << std::endl; }
    //! Declare and initialize propagation settings for both mass and translational dynamics.
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings_Descent =
            std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector_Descent, problemInput->getDescentTerminationSettings(), problemInput->getDependentVariablesToSave() );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE INTEGRATION SETTINGS: DESCENT           //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){ std::cout << "Create Integration Settings: Descent" << std::endl; }

    std::shared_ptr< IntegratorSettings<  > > integratorSettings_Descent =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationEndEpoch_Ascent,
              propagationStepSize, 1, false );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           ASSIGN INTERPOLATORS & BOUNDS: DESCENT           ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){ std::cout << "Pass Descent Interpolators to Vehicle Systems" << std::endl; }
    bislipSystems->setParameterInterpolator( Interpolators_Descent );

    if( debugInfo == 10 ){ std::cout << "Pass Descent Bounds to Vehicle Systems" << std::endl; }


    if( debugInfo == 10 ){ std::cout << "Pass Parameter Bounds to Vehicle Systems" << std::endl; }
    bislipSystems->setParameterLowerBounds( Interpolators_Descent_LB);
    bislipSystems->setParameterUpperBounds( Interpolators_Descent_UB );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          PROPAGATE TRAJECTORY: DESCENT         ////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){ std::cout << "Starting Descent propagation" << std::endl; }

    bislipSystems->setCurrentTrajectoryPhase( "Descent" );

    //! Propagate trajectory.
    SingleArcDynamicsSimulator< double > descentSimulation(
                bodyMap,
                integratorSettings_Descent,
                propagatorSettings_Descent );

    if( debugInfo == 10 ){ std::cout << "Descent Propagation done" << std::endl; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          RETRIEVE DATA MAPS          //////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 10 ){ std::cout << "Merge ASCENT and DESCENT Time History Maps" << std::endl; }

    if( debugInfo == 10 ){ std::cout << "   Declare Total Propagation Time History Map and Initialize with ASCENT Propagation Time History Map" << std::endl; }
    std::map< double, Eigen::VectorXd > propagationTimeHistoryMap = propagationTimeHistoryMap_Ascent;

    if( debugInfo == 10 ){ std::cout << "   Declare Total Dependent Variable Time History Map and Initialize with ASCENT Dependent Variable Time History Map" << std::endl; }
    std::map< double, Eigen::VectorXd > depVarTimeHistoryMap = depVarTimeHistoryMap_Ascent;

    if( debugInfo == 10 ){ std::cout << "   Retrieve DESCENT Propagation and Dependent Variable Dependent Variable Time History maps" << std::endl; }
    const std::map< double, Eigen::VectorXd > propagationTimeHistoryMap_Descent = descentSimulation.getEquationsOfMotionNumericalSolution( );
    const std::map< double, Eigen::VectorXd > depVarTimeHistoryMap_Descent = descentSimulation.getDependentVariableHistory( );

    if( debugInfo == 10 ){ std::cout << "   Append DESCENT Propagation Time History Map onto TOTAL Propagation Time History Map" << std::endl; }
    propagationTimeHistoryMap.insert( propagationTimeHistoryMap_Descent.begin(), propagationTimeHistoryMap_Descent.end() );

    if( debugInfo == 10 ){ std::cout << "   Append DESCENT Dependent Variable Time History Map onto TOTAL Dependent Variable Time History Map" << std::endl; }
    depVarTimeHistoryMap.insert( depVarTimeHistoryMap_Descent.begin(), depVarTimeHistoryMap_Descent.end() );


    if( debugInfo == 10 ){ std::cout << "Extract Data Map keys (Epochs) and place values into Epoch Vector" << std::endl; }
    //! Extract time vector ( Time History keys )
    std::vector< double > epochVector;
    for ( auto const& element : depVarTimeHistoryMap ) { epochVector.push_back( element.first ); }

    if( debugInfo == 10 ){ std::cout << "   Determine size of Epoch Vector (check)" << std::endl; }

    int epochVectorSize = 0;
    for( std::map< double, Eigen::VectorXd >::iterator it = depVarTimeHistoryMap.begin(); it != depVarTimeHistoryMap.end(); ++it)
    { epochVectorSize += 1; }

    //! Extract final epoch
    const double simulationEndEpoch_Descent =  ( propagationTimeHistoryMap_Descent.rbegin() )->first;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        CLEAN UP AND CALCULATE FITNESS                //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Initialize variables used objectives functions.
    double finalDistanceToTarget_deg_calc = initialDistanceToTarget_deg;
    //std::cout << "Cleaning up and extracting Time Histories" << std::endl;
    double objectiveHeight_Ascent_calc = objectiveHeight_Ascent;
    //double objectiveHeight_Ascent_calc = objectiveHeight_Ascent;
    double objectiveAirspeed_Ascent_calc = objectiveAirspeed_Ascent;
    double objectiveAirspeed_Descent_calc = objectiveAirspeed_Descent;
    //double tof = problemInput->getSimulationSettings()[ 1 ];

    //double maximum_DynamicPressure = constraint_DynamicPressure;
    //double maximum_MechanicalLoad = constraint_MechanicalLoad;
    //double finalMass_Descent = 0;
    const double Isp = bislipSystems->getSpecificImpulse();
    const double theoreticaldelV_Ascent  = g0 * Isp * log( bislipSystems->getInitialMass() / bodyMap.at(  vehicleName )->getVehicleSystems()->getDryMass() );



    if( initialMass_Descent < bodyMap.at(  vehicleName )->getVehicleSystems()->getDryMass() ) { initialMass_Descent = bodyMap.at(  vehicleName )->getVehicleSystems()->getDryMass(); }
    const double theoreticaldelV_Descent = g0 * Isp * log( initialMass_Descent / bodyMap.at(  vehicleName )->getVehicleSystems()->getDryMass() );

    //! Calculate Time of Flight
    const double timeOfFlight_Ascent = simulationEndEpoch_Ascent - simulationStartEpoch;
    const double timeOfFlight_Descent = simulationEndEpoch_Descent - simulationEndEpoch_Ascent;


    //! Declare number of elements in time history.
    const int rowsTotal = depVarTimeHistoryMap.size();
    rowsAscent = depVarTimeHistoryMap_Ascent.size();
    // long rowsDescent = depVarTimeHistoryMap_Descent.size();
    const int columns = ( ( descentSimulation.getDependentVariableHistory( ).begin() )->second ).size();


    if( debugInfo == 10 ){ std::cout << "Convert Dependent Variable Map to Matrix" << std::endl; }
    depVarTimeHistoryMatrix.resize( depVarTimeHistoryMap.size(), columns );
    depVarTimeHistoryMatrix = Eigen::MatrixXd::Zero( depVarTimeHistoryMap.size(), columns );

    if( debugInfo == 10 ){ std::cout << "   depVarTimeHistoryMap.size():  " << depVarTimeHistoryMap.size() << std::endl; }
    if( debugInfo == 10 ){ std::cout << "   depVarTimeHistoryMatrix.size():  " << depVarTimeHistoryMatrix.size() << std::endl; }


    for ( unsigned long i = 0; i < depVarTimeHistoryMap.size(); i++ )
    { depVarTimeHistoryMatrix.row( i ) = depVarTimeHistoryMap.at( epochVector[ i ] ); }

    if( debugInfo == 10 ){ std::cout << "Extract various columns form the Dependent Variable Matrix" << std::endl; }
    Eigen::VectorXd depVar_HeadingAngle                        = depVarTimeHistoryMatrix.col( 6 );
    Eigen::VectorXd depVar_TimeOfFlight                        = depVarTimeHistoryMatrix.col( 90 );
    Eigen::VectorXd depVar_PassengerFrameTotalAcceleration_x   = depVarTimeHistoryMatrix.col( 112 );
    Eigen::VectorXd depVar_PassengerFrameTotalAcceleration_y   = depVarTimeHistoryMatrix.col( 113 );
    Eigen::VectorXd depVar_PassengerFrameTotalAcceleration_z   = depVarTimeHistoryMatrix.col( 114 );

    if( debugInfo == 10 ){ std::cout << "   Convert Heading Angles to Positive Angles" << std::endl; }
    Eigen::VectorXd depVar_HeadingAnglePositive( rowsTotal );
    for ( int i = 0; i < rowsTotal; ++i )
    { depVar_HeadingAnglePositive( i ) = bislip::Variables::convertNegativeAnglesInDegreesToPositive( depVar_HeadingAngle( i ) ); }
    if( debugInfo == 10 ){ std::cout << "   Replace Heading Angle Vector with Positive Heading Angle Vector" << std::endl; }
    depVarTimeHistoryMatrix.col( 6 ) = depVar_HeadingAnglePositive;


    if( debugInfo == 10 ){ std::cout << "   Numerically Calculate Jerk" << std::endl; }
    //Eigen::VectorXd depVar_PassengerFrame_Jerk_x( rowsTotal );
    Eigen::VectorXd jerk_x = bislip::Variables::computeNumericalDerivativeOfVector( depVar_PassengerFrameTotalAcceleration_x, depVar_TimeOfFlight );
    //Eigen::VectorXd depVar_PassengerFrame_Jerk_y( rowsTotal );
    Eigen::VectorXd jerk_y = bislip::Variables::computeNumericalDerivativeOfVector( depVar_PassengerFrameTotalAcceleration_y, depVar_TimeOfFlight );
    //Eigen::VectorXd depVar_PassengerFrame_Jerk_z( rowsTotal );
    Eigen::VectorXd jerk_z = bislip::Variables::computeNumericalDerivativeOfVector( depVar_PassengerFrameTotalAcceleration_z, depVar_TimeOfFlight );

    if( debugInfo == 10 ){ std::cout << "   Resize Time History Matrix to Append Jerk Vector" << std::endl; }
    depVarTimeHistoryMatrix.conservativeResize( depVarTimeHistoryMatrix.rows(), depVarTimeHistoryMatrix.cols() + 3 );

    if( debugInfo == 10 ){ std::cout << "   Append Jerk Vectors to Time History Matrix" << std::endl; }
    depVarTimeHistoryMatrix.col( depVarTimeHistoryMatrix.cols() - 3 ) = jerk_x;
    depVarTimeHistoryMatrix.col( depVarTimeHistoryMatrix.cols() - 2 ) = jerk_y;
    depVarTimeHistoryMatrix.col( depVarTimeHistoryMatrix.cols() - 1 ) = jerk_z;



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PRINT SIMULATION OUTPUT TO FILE               //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( problemInput->getEvolutionEvaluationFlag() != true )
    {
        if( debugInfo == 10 ){ std::cout << "Get this run's file name suffix" << std::endl; }
        std::string simulation_file_name_suffix = problemInput->getRerunFileNameSuffix();

        if( debugInfo == 10 ){ std::cout << "Create file names" << std::endl; }
        std::string complete_file_name_Prop                     = "propagationHistory_" + simulation_file_name_suffix + ".dat";
        std::string complete_file_name_DepVar                   = "dependentVariables_" + simulation_file_name_suffix + ".dat";
        std::string complete_file_name_interpolators_Ascent     = "evaluatedInterpolatorsAscent_" + simulation_file_name_suffix + ".dat";
        std::string complete_file_name_interpolators_Ascent_LB  = "evaluatedInterpolatorsAscent_LB_" + simulation_file_name_suffix + ".dat";
        std::string complete_file_name_interpolators_Ascent_UB  = "evaluatedInterpolatorsAscent_UB_" + simulation_file_name_suffix + ".dat";
        std::string complete_file_name_map_DV_mapped_Ascent     = "map_DV_mapped_Ascent_" + simulation_file_name_suffix + ".dat";
        std::string complete_file_name_interpolators_Descent    = "evaluatedInterpolatorsDescent_" + simulation_file_name_suffix + ".dat";
        std::string complete_file_name_interpolators_Descent_LB = "evaluatedInterpolatorsDescent_LB_" + simulation_file_name_suffix + ".dat";
        std::string complete_file_name_interpolators_Descent_UB = "evaluatedInterpolatorsDescent_UB_" + simulation_file_name_suffix + ".dat";
        std::string complete_file_name_map_DV_mapped_Descent    = "map_DV_mapped_Descent_" + simulation_file_name_suffix + ".dat";

        if( debugInfo == 10 ){ std::cout << "Saving Propagation" << std::endl; }
        if( debugInfo == 10 ){ std::cout << "propagationTimeHistoryMap.size() = " << propagationTimeHistoryMap.size() << std::endl; }

        //! Write propagation history to file.
        tudat::input_output::writeDataMapToTextFile( propagationTimeHistoryMap,
                                                     complete_file_name_Prop,
                                                     problemInput->getOutputPath() + problemInput->getOutputSubFolder() + "/" + "propagationHistory" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        if( debugInfo == 10 ){ std::cout << "Propagation Saved" << std::endl; }

        if( debugInfo == 10 ){ std::cout << "Saving Dependent Variables" << std::endl; }
        if( debugInfo == 10 ){ std::cout << "depVarTimeHistoryMap.size() = " << depVarTimeHistoryMap.size() << std::endl; }

        std::map< double, Eigen::VectorXd > depVarTimeHistoryMapOutput;
        for ( int i = 0; i < rowsTotal; i++ )
        { depVarTimeHistoryMapOutput[ epochVector[ i ] ] = depVarTimeHistoryMatrix.row( i ); }

        //! Write dependant variable to file.
        tudat::input_output::writeDataMapToTextFile( depVarTimeHistoryMapOutput,
                                                     complete_file_name_DepVar,
                                                     problemInput->getOutputPath() + problemInput->getOutputSubFolder() + "/" + "dependentVariables" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );

        if( debugInfo == 10 ){ std::cout << "Dependent Variables Saved" << std::endl; }

        if( debugInfo == 10 ){ std::cout << "Saving Evaluated Interpolators - Ascent" << std::endl; }

        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsAscent,
                                                     complete_file_name_interpolators_Ascent,
                                                     problemInput->getOutputPath() + problemInput->getOutputSubFolder() + "/" + "evaluatedInterpolatorsAscent" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );

        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsAscent_LB,
                                                     complete_file_name_interpolators_Ascent_LB,
                                                     problemInput->getOutputPath() + problemInput->getOutputSubFolder() + "/" + "evaluatedInterpolatorsAscent_LB" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsAscent_UB,
                                                     complete_file_name_interpolators_Ascent_UB,
                                                     problemInput->getOutputPath() + problemInput->getOutputSubFolder() + "/" + "evaluatedInterpolatorsAscent_UB" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );

        if( debugInfo == 10 ){ std::cout << "Evaluated Interpolators - Ascent Saved" << std::endl; }

        if( debugInfo == 10 ){ std::cout << "Saving Evaluated Interpolators - Descent" << std::endl; }
        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsDescent,
                                                     complete_file_name_interpolators_Descent,
                                                     problemInput->getOutputPath() + problemInput->getOutputSubFolder() + "/" + "evaluatedInterpolatorsDescent" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsDescent_LB,
                                                     complete_file_name_interpolators_Descent_LB,
                                                     problemInput->getOutputPath() + problemInput->getOutputSubFolder() + "/" + "evaluatedInterpolatorsDescent_LB" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsDescent_UB,
                                                     complete_file_name_interpolators_Descent_UB,
                                                     problemInput->getOutputPath() + problemInput->getOutputSubFolder() + "/" + "evaluatedInterpolatorsDescent_UB" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );


        if( debugInfo == 10 ){ std::cout << "Evaluated Interpolators - Descent Saved" << std::endl; }

        if( debugInfo == 10 ){ std::cout << "Saving Mapped Nodal Values - Ascent" << std::endl; }
        tudat::input_output::writeDataMapToTextFile( map_DV_mapped_Ascent,
                                                     complete_file_name_map_DV_mapped_Ascent,
                                                     problemInput->getOutputPath() + problemInput->getOutputSubFolder() + "/" + "map_DV_mapped_Ascent" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        if( debugInfo == 10 ){ std::cout << "Decision Vector - Ascent Saved" << std::endl; }

        if( debugInfo == 10 ){ std::cout << "Saving Mapped Nodal Values - Descent" << std::endl; }
        tudat::input_output::writeDataMapToTextFile( map_DV_mapped_Descent,
                                                     complete_file_name_map_DV_mapped_Descent,
                                                     problemInput->getOutputPath() + problemInput->getOutputSubFolder() + "/" + "map_DV_mapped_Descent" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        if( debugInfo == 10 ){ std::cout << "Decision Vector - Descent Saved" << std::endl; }
    }


    //if( debugInfo == 10 ){ std::cout << "Printing to terminal" << std::endl; }
    //int individualNumber = problemInput->getIndividualNumber( ) + 1;
    //problemInput->setIndividualNumber( individualNumber );

    /*
    std::cout  << "Individual " << individualNumber << " | " ;
    for( unsigned int i = 0; i < fitness.size() - 1; i++)
    {
        std::cout  << fitness[ i ] << " | " ;
    }
    //std::cout  << fitness[ fitness.size() - 1 ] << "  ||  " << rowsAscent << "  ||  " << tof << "  ||  " << " Theoretical delV = " << theoreticaldelV_Ascent << "  ||  " << " Goal delV = " << goaldelV_Ascent << "  ||  " << " Actual delV = " << actualdelV_Ascent << std::endl;
    //    std::cout  << fitness[ fitness.size() - 1 ] << "  ||  " << rowsAscent << "  ||  " << tof <<  std::endl;
    //std::cout  << fitness[ fitness.size() - 1 ] << " ||  " << timeOfFlight << " ||  " << saveTrajectoryOutput << " ||  " << rowsAscent << std::endl;
    std::cout  << fitness[ fitness.size() - 1 ] << " ||  " << timeOfFlight << " ||  " << saveTrajectoryOutput << std::endl;

*/
    //  }

    //std::cout << "objectiveHeight_Ascent: " << objectiveHeight_Ascent << std::endl;
    //std::cout << "objectiveHeight_Ascent_calc: " << objectiveHeight_Ascent_calc << std::endl;



    // std::cout <<dynamicsSimulator.getPropagationTerminationReason( )->getPropagationTerminationReason( ) << std::endl;
    //std::this_thread::sleep_for( std::chrono::nanoseconds( 10 ) );
    //using namespace std::this_thread; // sleep_for, sleep_until
    //using namespace std::chrono; // nanoseconds, system_clock, seconds

    //return decisionVectorEvaluationOutput;
    if( debugInfo == 10 ){ std::cout << "Exiting Decision Vector Evaluator" << std::endl; }

}


void reruns( const std::shared_ptr< bislip::ProblemInput > &problemInput,
             const tudat::simulation_setup::NamedBodyMap& bodyMap,
             const int &topIndividuals )
{
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap.at( problemInput->getVehicleName() )->getBislipSystems();

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 10 ){ std::cout << "Set the Evolution Evaluation Flag to FALSE. A FALSE flag triggers the printing." << std::endl; }
    problemInput->setEvolutionEvaluationFlag( false );

    if( debugInfo == 10 ){ std::cout << "Extract population/fitness maps" << std::endl; }
    std::map < std::string, Eigen::VectorXd > fitnessMap = problemInput->getFitness();
    std::map < std::string, Eigen::VectorXd > populationMap = problemInput->getPopulation();

    const int pureFitnessRows = fitnessMap.size();
    const int pureFitnessColumns = ( fitnessMap.begin()->second ).size() - 2;
    std::vector< std::vector< double > > pureFitness( pureFitnessRows, std::vector< double >( pureFitnessColumns, 0.0 ) );

    if( debugInfo == 10 ){ std::cout << "Pull out population/fitness keys (file name tags)" << std::endl; }
    std::vector< std::string > suffixes;
    for( auto const& element : fitnessMap ) { suffixes.push_back( element.first ); }

    if( debugInfo == 10 ){ std::cout << "Convert relevant fitness map content to STL vector of STL vectors" << std::endl; }
    for( int i = 0; i < pureFitnessRows; i++ )
    {
        for( int j = 0; j < pureFitnessColumns; j++ )
        {
            pureFitness[ i ][ j ] = ( fitnessMap.at( suffixes[ i ] ) )( j + 2 );
        }
    }

    if( debugInfo == 10 ){ std::cout << "Sort fitness vector" << std::endl; }
    std::vector< std::vector< double >::size_type > ranking = pagmo::select_best_N_mo( pureFitness, topIndividuals );

    Eigen::VectorXd rerunDecisionVectorEigen;
    std::vector< double > rerunDecisionVectorSTL;
    std::string rerunSuffix;
    const int purePopulationColumns = ( populationMap.begin()->second ).size() - 2;

    std::cout<< "Re-evaluate top " << topIndividuals << " individuals" << std::endl;
    for( int i = 0; i < topIndividuals; i++ )
    {
        rerunSuffix = suffixes[ ranking[ i ] ];

        std::cout<< "   Individual " << int( populationMap.at( rerunSuffix )( 1 ) ) << " : " << rerunSuffix << " | " ;
        for( unsigned int j = 0; j < pureFitnessColumns - 1; j++)
        {
            std::cout  << pureFitness[ ranking[ i ] ][ j ] << " | " ;
        }
        std::cout  << pureFitness[ ranking[ i ] ][ pureFitnessColumns - 1 ] << std::endl;

        problemInput->setRerunFileNameSuffix( rerunSuffix );
        rerunDecisionVectorEigen = ( populationMap.at( rerunSuffix ) ).tail( purePopulationColumns );
        rerunDecisionVectorSTL = tudat::utilities::convertEigenVectorToStlVector( rerunDecisionVectorEigen );

        Eigen::MatrixXd depVarTimeHistoryMatrix;
        int rowsAscent;
        bislip::Variables::decisionVectorEvaluation( rerunDecisionVectorSTL, problemInput, bodyMap, depVarTimeHistoryMatrix, rowsAscent );

        //! Flag the individuals as being printed
        if( debugInfo == 10 ){ std::cout << "Flag individual as being printed" << std::endl; }
        populationMap.at( rerunSuffix )( 0 ) = 1;
        fitnessMap.at( rerunSuffix )( 0 ) = 1;

    }

    if( debugInfo == 10 ){ std::cout << "Reset maps with the printed flags" << std::endl; }
    problemInput->setPopulation( populationMap );
    problemInput->setFitness( fitnessMap );
    std::cout<< "  " << std::endl;

}




}; //namespace Variables
                 } // namespace bislip










