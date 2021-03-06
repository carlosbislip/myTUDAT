/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved.
 * qqwqwqw
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120522    A. Ronse          First creation of code.
 *
 *    References
 *      Williams, Dr. David R., "Moon Fact Sheet", NASA (National Space Science Data Center),
 *         http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html, last accessed: 22 May 2012
 *
 *    Notes
 *
 */

//! Mine
#include "Space4ErrBodyProblem.h"

Space4ErrBodyProblem::Space4ErrBodyProblem(
        const std::shared_ptr< bislip::ProblemInput > &problemInput,
        const tudat::simulation_setup::NamedBodyMap& bodyMap ):
    problemInput_( problemInput ),
    bodyMap_( bodyMap )
{ }

//! Descriptive name of the problem
std::string Space4ErrBodyProblem::get_name() const
{
    return problemInput_->getProblemName();
}

//! Get bounds
std::pair< std::vector< double >, std::vector< double > > Space4ErrBodyProblem::get_bounds() const
{
    return { problemInput_->getDecisionVectorBounds()[ 0 ], problemInput_->getDecisionVectorBounds()[ 1 ] };
}

std::vector< double > Space4ErrBodyProblem::fitness( const std::vector< double > &x )  const
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
    ///
    const std::string vehicleName =  problemInput_->getVehicleName();

    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap_.at( vehicleName )->getVehicleSystems();

    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at(  vehicleName )->getBislipSystems();

    bislipSystems->setStartingMillis( bislip::Variables::millis_since_midnight() );

    bislipSystems->setDebugInfo( problemInput_->getSimulationSettings()[ 9 ] );

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 )
    {
        std::cout << "Printing decision vector to screen" << std::endl;
        std::cout << "x =  " << std::endl;
        for( int i = 0; i < int( x.size( ) ); i++ ) { std::cout << "x[ " << i << " ] = " << x[ i ] << std::endl; }
    }

    if( debugInfo == 1 ){std::cout << "Unpacking data" << std::endl; }

    const std::string centralBodyName = problemInput_->getCentralBodies()[ 0 ];

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
    const double simulationEndEpoch = simulationStartEpoch + problemInput_->getSimulationSettings()[ 6 ];

    //! Declare and initialize numerical integration fixed step size.
    const double propagationStepSize = problemInput_->getSimulationSettings()[ 7 ];

    const double guidanceStepSize = problemInput_->getSimulationSettings()[ 8 ];

    //! Declare and initialize number of control nodes.
    const unsigned long nodesAscent = problemInput_->getSimulationSettings().rbegin()[ 1 ];
    const unsigned long nodesDescent = problemInput_->getSimulationSettings().rbegin()[ 0 ];

    if( debugInfo == 1 ){std::cout << "nodesAscent = " << nodesAscent <<std::endl; }
    if( debugInfo == 1 ){std::cout << "nodesDescent = " << nodesDescent <<std::endl; }


    //! Declare and initialize position vector of moment reference center
    const Eigen::Vector3d R_mrc( problemInput_->getVehicleParameters()[ 3 ], problemInput_->getVehicleParameters()[ 4 ], problemInput_->getVehicleParameters()[ 5 ] ); // m

    //! Declare and initialize position vector of center of mass
    const Eigen::Vector3d R_com( problemInput_->getVehicleParameters()[ 6 ], problemInput_->getVehicleParameters()[ 7 ], problemInput_->getVehicleParameters()[ 8 ] ); // m

    //! Declare and initialize position vector of center of thrust
    const Eigen::Vector3d R_cot( problemInput_->getVehicleParameters()[ 9 ], problemInput_->getVehicleParameters()[ 10 ], problemInput_->getVehicleParameters()[ 11 ] ); // m

    //! Declare and initialize initial mass
    double initialMass_Ascent = problemInput_->getVehicleParameters()[ 12 ]; // kg

    double dryMass = problemInput_->getVehicleParameters()[ 13 ]; // kg

    //! Declare and initialize starting height
    const double initialHeight_Ascent = problemInput_->getInitialConditions()[ 2 ]; // m

    //! Declare and initialize starting position coordinates.
    const double initialLat_deg = problemInput_->getInitialConditions()[ 0 ];
    const double initialLon_deg = problemInput_->getInitialConditions()[ 1 ];

    //! Declare and initialize final position coordinates and additional termination conditions
    const double targetLat_deg = problemInput_->getConstraints()[ 0 ];
    const double targetLon_deg = problemInput_->getConstraints()[ 1 ];

    //! Declare and initialize various constraints
    const double finalDistanceToTarget_deg         = problemInput_->getConstraints()[ 2 ];
    const double objectiveHeight_Descent           = problemInput_->getConstraints()[ 4 ];
    //const double constraint_MechanicalLoad         = problemInput_->getConstraints()[ 6 ];
    const double constraint_ChapmanHeatFlux        = problemInput_->getConstraints()[ 7 ];
    const double constraint_DynamicPressure        = problemInput_->getConstraints()[ 8 ];
    const double constraint_PitchMomentCoefficient = problemInput_->getConstraints()[ 9 ];
    const double constraint_BendingMoment          = problemInput_->getConstraints()[ 10 ];

    //! Convert angles from degrees to radians
    const double initialLat_rad             = unit_conversions::convertDegreesToRadians( initialLat_deg );
    const double initialLon_rad             = unit_conversions::convertDegreesToRadians( initialLon_deg );
    const double targetLat_rad              = unit_conversions::convertDegreesToRadians( targetLat_deg );
    const double targetLon_rad              = unit_conversions::convertDegreesToRadians( targetLon_deg );

    //! Pre-define various variables used to determine fitness.
    double targetLat_deg_calc          = initialLat_deg;
    double targetLon_deg_calc          = initialLon_deg;
    double initialDistanceToTarget_rad = bislip::Variables::computeAngularDistance( initialLat_rad, initialLon_rad, targetLat_rad, targetLon_rad );
    double initialDistanceToTarget_deg = unit_conversions::convertRadiansToDegrees( initialDistanceToTarget_rad );
    double initialDistanceToTarget_m   = radiusEarth * initialDistanceToTarget_rad;


    if( debugInfo == 1 ){ std::cout << "Reset non-static placeholders." << std::endl; }
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
    bislipSystems->setCurrentBodyFlapAngle( tudat::unit_conversions::convertDegreesToRadians( problemInput_->getInitialConditions()[ 5 ] ) );
    bislipSystems->setCurrentElevonAngle( tudat::unit_conversions::convertDegreesToRadians( problemInput_->getInitialConditions()[ 6 ] ) );
    bislipSystems->setCumulativeDistanceTravelled( 0.0 );
    bislipSystems->setCumulativeAngularDistanceTravelled( 0.0 );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE NODAL STRUCTURE             /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){std::cout << "Creating nodal structure" << std::endl; }

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

    if( debugInfo == 1 ){std::cout << "     Re-allocate Ascent DVs" << std::endl; }

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
    //    const unsigned long N = ( problemInput_->getAscentParameterList().size() - 6 ) * nodesAscent - 1;
    const unsigned long N = problemInput_->getAscentParameterList().size();// - 6 ) * nodesAscent - 1;

    //! Declare and initialize various parameters common to the entire trajectory.
    double initialFlightPathAngle               = x[ N - 7 ];
    const double initialLaunchHeading           = x[ N - 6 ];
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

        if( debugInfo == 1 ){ std::cout << "    Setting Trajectory Phase" << std::endl; }
        bislipSystems->setCurrentTrajectoryPhase( "Descent" );
        if( debugInfo == 1 ){std::cout << "         Trajectory Phase         = " << bislipSystems->getCurrentTrajectoryPhase() << std::endl; }

        initialFlightPathAngle = bislipSystems->getInitialFlightPathAngle();

    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "    Setting Trajectory Phase" << std::endl; }
        bislipSystems->setCurrentTrajectoryPhase( "Ascent" );
        if( debugInfo == 1 ){std::cout << "         Trajectory Phase         = " << bislipSystems->getCurrentTrajectoryPhase() << std::endl; }

        if( debugInfo == 1 ){std::cout << "    Identifying additional initial conditions" << std::endl; }

        if( bislipSystems->getInitialMachNumber() > 0.0 && bislipSystems->getInitialHeight() > 0.0 )
        {

            if( debugInfo == 1 ){std::cout << "    Vehicle initialized with Speed and Height" << std::endl; }

            initialAirspeed_Ascent = bislipSystems->getInitialSpeedOfSound() * bislipSystems->getInitialMachNumber();
            bislipSystems->setInitialAirspeed( initialAirspeed_Ascent );
        }
        else if( bislipSystems->getInitialMachNumber() == 0.0 && bislipSystems->getInitialHeight() != 0.0 )
        {
            if( debugInfo == 1 ){std::cout << "    Vehicle initialized at rest and with height---> Initial Airspeed is an optimization parameter" << std::endl; }

            bislipSystems->setInitialMachNumber( initialAirspeed_Ascent / bislipSystems->getInitialSpeedOfSound() );
        }
        else if( bislipSystems->getInitialMachNumber() == 0.0 && bislipSystems->getInitialHeight() == 0.0 )
        {
            if( debugInfo == 1 ){std::cout << "    Vehicle initialized at rest and on the surface" << std::endl; }

            bislipSystems->setInitialMachNumber( initialAirspeed_Ascent / bislipSystems->getInitialSpeedOfSound() );
        }
    }

    bislipSystems->setInitialDynamicPressure( bislipSystems->getInitialDensity() * initialAirspeed_Ascent * initialAirspeed_Ascent / 2 );


    if( debugInfo == 1 ){std::cout << "         Initial Height           = " << bislipSystems->getInitialHeight() << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Initial Altitude         = " << bislipSystems->getInitialAltitude() << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Initial Mach Number      = " << bislipSystems->getInitialMachNumber() << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Initial Airspeed         = " << bislipSystems->getInitialAirspeed() << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Initial Speed of Sound   = " << bislipSystems->getInitialSpeedOfSound() << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Initial Dynamic Pressure = " << bislipSystems->getInitialDynamicPressure() << std::endl; }
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
            double increasedVehicleDryMass = vehicleSystems->getDryMass() * ( 1.0 + 0.3*( ( additionalMass / initialMass_Ascent ) ) );
            vehicleSystems->setDryMass( increasedVehicleDryMass );
        }
    }

    if( debugInfo == 1 ){std::cout << "         Initial Launch Heading    = " << initialLaunchHeading << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Initial Flight-Path Angle = " << initialFlightPathAngle << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Initial Launch Airspeed   = " << initialAirspeed_Ascent << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Objective Airspeed        = " << objectiveAirspeed_Ascent << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Objective Height          = " << objectiveHeight_Ascent << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Initial Mass - Ascent     = " << initialMass_Ascent << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Additional Mass           = " << additionalMass << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Ascent Termination Ratio  = " << ascentTerminationDistanceRatio << std::endl; }

    bodyMap_.at( vehicleName )->setConstantBodyMass( initialMass_Ascent + additionalMass );
    bislipSystems->setInitialMass( initialMass_Ascent + additionalMass );
    bislipSystems->setInitialHeight( initialHeight_Ascent );
    bislipSystems->setInitialAltitude( initialHeight_Ascent + tudat::spice_interface::getAverageRadius( centralBodyName ) );
    bislipSystems->setAscentTerminationDistanceRatio( ascentTerminationDistanceRatio );

    if( debugInfo == 1 ){std::cout << "     Re-allocate Descent DVs" << std::endl; }

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
    //const unsigned long NN = ( problemInput_->getDescentParameterList().size() - 2 ) * nodesDescent - 1;
    const unsigned long NN = problemInput_->getDescentParameterList().size();// - 2 ) * nodesDescent - 1;

    //! Declare and initialize last parameter common to the entire trajectory.
    //!     Final velocity.
    //!     Skip suppression timing trigger.
    const double objectiveAirspeed_Descent    = x.rbegin()[ 2 ];
    const double skipSuppressionTimingTrigger = x.rbegin()[ 1 ];
    const double constraint_MechanicalLoad    = x.rbegin()[ 0 ];

    bislipSystems->setMechanicalLoadConstraint( constraint_MechanicalLoad );

    if( debugInfo == 1 ){std::cout << "alpha_deg_Descent            = " << alpha_deg_Descent << std::endl; }
    if( debugInfo == 1 ){std::cout << "sigma_deg_Descent            = " << sigma_deg_Descent << std::endl; }
    if( debugInfo == 1 ){std::cout << "eps_T_deg_Descent            = " << eps_T_deg_Descent << std::endl; }
    if( debugInfo == 1 ){std::cout << "phi_T_deg_Descent            = " << phi_T_deg_Descent << std::endl; }
    if( debugInfo == 1 ){std::cout << "throttle_Descent             = " << throttle_Descent << std::endl; }
    if( debugInfo == 1 ){std::cout << "Objective Airspeed Descent   = " << objectiveAirspeed_Descent << std::endl; }
    if( debugInfo == 1 ){std::cout << "Skip Suppresion Trigger Time = " << skipSuppressionTimingTrigger << std::endl; }

    bislipSystems->setSkipSuppressionTimingTrigger( skipSuppressionTimingTrigger );

    const double goaldelV_Descent = std::abs( objectiveAirspeed_Descent - objectiveAirspeed_Ascent );

    if( debugInfo == 1 ){std::cout << "Create vector of node locations for ascent" << std::endl; }
    //! Create vector of node locations for Ascent phase.
    xn_Ascent( 0 ) = 0;
    for( unsigned int i = 0; i < nodesAscent - 1; i++ )
    {
        xn_Ascent( i + 1 )        = xn_Ascent( i ) + xn_interval_Ascent( i ) / xn_interval_Ascent.sum();
    }

    if( debugInfo == 1 ){std::cout << "Create vector of node locations for descent" << std::endl; }
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

    if( debugInfo == 1 ){std::cout << "Complete known ascent states" << std::endl; }

    Eigen::Vector6d initialStateAscent = problemInput_->getInitialState_Spherical();
    initialStateAscent( 3 ) = initialAirspeed_Ascent;
    initialStateAscent( 4 ) = tudat::unit_conversions::convertDegreesToRadians( initialFlightPathAngle );
    initialStateAscent( 5 ) = tudat::unit_conversions::convertDegreesToRadians( initialLaunchHeading );

    bislipSystems->setInitialHeading( initialLaunchHeading );
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
                bodyCenteredBodyFixedState, simulationStartEpoch, problemInput_->getEarthRotationalEphemeris() );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE INTERPOLATORS             ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Creating interpolators" << std::endl; }

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

    if( debugInfo == 1 ){ std::cout << "E_hat_i   = " << E_hat_i << std::endl; }
    if( debugInfo == 1 ){ std::cout << "E_hat_max = " << E_hat_max << std::endl; }
    if( debugInfo == 1 ){ std::cout << "E_hat_f   = " << E_hat_f << std::endl; }

    //! Map normalized specific energy levels to control node locations.
    Eigen::VectorXd E_mapped_Ascent( xn_Ascent.size() );
    Eigen::VectorXd E_mapped_Descent( xn_Descent.size() );

    E_mapped_Ascent  = ( ( E_max - E_i ) * xn_Ascent.array() + E_i ) / E_max;

    if( bislipSystems->getValidationFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "Validation Case ---> Overwriting ascent info with descent info "<< std::endl; }

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


    std::map< bislip::Parameters::Interpolators, Eigen::MatrixXd > ascentParameterBoundsMap = problemInput_->getAscentParameterBoundsMap();


    if( debugInfo == 1 ){ std::cout << "Mapping Ascent DVs" << std::endl; }
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

    std::map< bislip::Parameters::Interpolators, Eigen::MatrixXd > descentParameterBoundsMap = problemInput_->getDescentParameterBoundsMap();



    if( debugInfo == 1 ){ std::cout << "Mapping Descent DVs" << std::endl; }
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


    if( debugInfo == 1 ){ std::cout << "Creating Optimization Interpolators' Settings" << std::endl; }

    std::shared_ptr< tudat::interpolators::InterpolatorSettings > interpolatorSettings = std::make_shared< tudat::interpolators::InterpolatorSettings >( tudat::interpolators::hermite_spline_interpolator );
    std::shared_ptr< tudat::interpolators::InterpolatorSettings > interpolatorBoundariesSettings = std::make_shared< tudat::interpolators::InterpolatorSettings >( tudat::interpolators::linear_interpolator );
    std::map< bislip::Parameters::Interpolators, std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > > Interpolators_Ascent, Interpolators_Descent;
    std::map< bislip::Parameters::Interpolators, std::shared_ptr< tudat::interpolators::OneDimensionalInterpolator< double, double > > > Interpolators_Ascent_LB, Interpolators_Ascent_UB, Interpolators_Descent_LB, Interpolators_Descent_UB;

    if( debugInfo == 1 ){ std::cout << "Creating Ascent Interpolators" << std::endl; }

    Interpolators_Ascent[ bislip::Parameters::Interpolators::AngleOfAttack ]        = bislip::Variables::createOneDimensionalHermiteInterpolator( alpha_deg_Ascent, E_mapped_Ascent, map_alpha_deg_Ascent, interpolatorSettings );
    Interpolators_Ascent[ bislip::Parameters::Interpolators::BankAngle ]            = bislip::Variables::createOneDimensionalHermiteInterpolator( sigma_deg_Ascent, E_mapped_Ascent, map_sigma_deg_Ascent, interpolatorSettings );
    Interpolators_Ascent[ bislip::Parameters::Interpolators::ThrustElevationAngle ] = bislip::Variables::createOneDimensionalHermiteInterpolator( eps_T_deg_Ascent, E_mapped_Ascent, map_eps_T_deg_Ascent, interpolatorSettings );
    Interpolators_Ascent[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ]   = bislip::Variables::createOneDimensionalHermiteInterpolator( phi_T_deg_Ascent, E_mapped_Ascent, map_phi_T_deg_Ascent, interpolatorSettings );
    Interpolators_Ascent[ bislip::Parameters::Interpolators::ThrottleSetting ]      = bislip::Variables::createOneDimensionalHermiteInterpolator( throttle_Ascent, E_mapped_Ascent, map_throttle_Ascent, interpolatorSettings );

    if( debugInfo == 1 ){ std::cout << "Creating Ascent Interpolators Boundaries" << std::endl; }

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


    if( debugInfo == 1 ){ std::cout << "Creating Descent Interpolators" << std::endl; }
    //! Declare and initialize interpolators for Descent phase.

    if( debugInfo == 1 ){ std::cout << "   Angle of Attack" << std::endl; }
    Interpolators_Descent[ bislip::Parameters::Interpolators::AngleOfAttack ]        = bislip::Variables::createOneDimensionalHermiteInterpolator( alpha_deg_Descent, E_mapped_Descent, map_alpha_deg_Descent, interpolatorSettings );
    if( debugInfo == 1 ){ std::cout << "   Bank Angle" << std::endl; }
    Interpolators_Descent[ bislip::Parameters::Interpolators::BankAngle ]            = bislip::Variables::createOneDimensionalHermiteInterpolator( sigma_deg_Descent, E_mapped_Descent, map_sigma_deg_Descent, interpolatorSettings );
    if( debugInfo == 1 ){ std::cout << "   Thrust Elevation Angle" << std::endl; }
    Interpolators_Descent[ bislip::Parameters::Interpolators::ThrustElevationAngle ] = bislip::Variables::createOneDimensionalHermiteInterpolator( eps_T_deg_Descent, E_mapped_Descent, map_eps_T_deg_Descent, interpolatorSettings );
    if( debugInfo == 1 ){ std::cout << "   Thrust Azimuth Angle" << std::endl; }
    Interpolators_Descent[ bislip::Parameters::Interpolators::ThrustAzimuthAngle ]   = bislip::Variables::createOneDimensionalHermiteInterpolator( phi_T_deg_Descent, E_mapped_Descent, map_phi_T_deg_Descent, interpolatorSettings );
    if( debugInfo == 1 ){ std::cout << "   Throttle Setting" << std::endl; }
    Interpolators_Descent[ bislip::Parameters::Interpolators::ThrottleSetting ]      = bislip::Variables::createOneDimensionalHermiteInterpolator( throttle_Descent, E_mapped_Descent, map_throttle_Descent, interpolatorSettings );


    if( debugInfo == 1 ){ std::cout << "Creating Descent Interpolators Boundaries" << std::endl; }

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

    if( debugInfo == 1 ){ std::cout << "Evaluate Interpolators" << std::endl; }
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

    if( debugInfo == 1 ){ std::cout << "Interpolators Evaluated." << std::endl; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           ASSIGN INTERPOLATORS & BOUNDS: ASCENT           /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Pass Ascent Interpolators to vehicle systems" << std::endl;}

    //! Pass Ascent phase interpolators to vehicle systems.
    bislipSystems->setParameterInterpolator( Interpolators_Ascent );

    if( debugInfo == 1 ){ std::cout << "Pass Ascent Bounds to vehicle systems" << std::endl;}

    //! Pass parameter bounds to vehicle systems.
    //bislipSystems->setParameterBounds( problemInput_->getAscentParameterBoundsMap() );
    bislipSystems->setParameterLowerBounds( Interpolators_Ascent_LB);
    bislipSystems->setParameterUpperBounds( Interpolators_Ascent_UB );


    double simulationEndEpoch_Ascent = simulationStartEpoch;
    Eigen::Vector6d systemInitialState_Descent;
    Eigen::VectorXd dependentVariableFinalState_Ascent;
    double finalAirspeed_Ascent = initialAirspeed_Ascent;
    double finalMass_Ascent = initialMass_Ascent;
    double initialMass_Descent = initialMass_Ascent;
    double finalSpecificEnergy_Ascent = E_i;

    std::map< double, Eigen::VectorXd > propagationTimeHistoryMap_Ascent;
    std::map< double, Eigen::VectorXd > depVarTimeHistoryMap_Ascent;

    //if( E_i < 2 * E_f )
    // {

    //if( debugInfo == 1 ){ std::cout << "Initial Specific Energy is LOWER than Final Specific Energy" << std::endl; }

    //( bislipSystems->getParameterInterpolator( bislip::Parameters::Optimization::ThrustElevationAngle ) )->interpolate( bislip::Variables::computeNormalizedSpecificEnergy( currentAltitude, currentAirspeed, bislipSystems->getMaximumSpecificEnergy() ) );



    if( debugInfo == 1 ){ std::cout << "Setting various initial values." << std::endl; }

    if( debugInfo == 1 ){std::cout << "    Evaluating Angle of Attack Interpolator" << std::endl; }
    if( bislipSystems->getValidationFlag() == true )
    {
        if( debugInfo == 1 ){std::cout << "         Validation Interpolator" << std::endl; }
        bislipSystems->setCurrentAngleOfAttack( tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getKourouAngleOfAttackInterpolator( ) )->interpolate( 0 ) ) );
    }
    else
    {
        if( debugInfo == 1 ){std::cout << "         Optimization Interpolator" << std::endl; }
        bislipSystems->setCurrentAngleOfAttack( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::AngleOfAttack, bodyMap_, vehicleName, centralBodyName ) ) );
    }
    if( debugInfo == 1 ){std::cout << "         Initial Angle of Attack = " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentAngleOfAttack() ) << std::endl; }
    /*
    if( debugInfo == 1 ){std::cout << "    Hardcoding Bank Angle" << std::endl; }
    bislipSystems->setCurrentBankAngle( 0.0 );
    if( debugInfo == 1 ){std::cout << "         Initial Bank Angle = " << bislipSystems->getCurrentBankAngle() << std::endl; }

    if( debugInfo == 1 ){std::cout << "    Repeating Flight-Path Angle" << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Initial Flight-Path Angle = " << bislipSystems->getInitialFlightPathAngle() << std::endl; }

    if( debugInfo == 1 ){std::cout << "    Repeating Heading Angle" << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Initial Heading Angle = " << bislipSystems->getInitialHeading() << std::endl; }

    if( debugInfo == 1 ){std::cout << "    Evaluating Thrust Azimuth Angle Interpolator" << std::endl; }
    bislipSystems->setCurrentThrustAzimuthAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::ThrustAzimuthAngle, bodyMap_, vehicleName ) ) );
    if( debugInfo == 1 ){std::cout << "         Initial Thrust Azimuth Angle   = " << bislipSystems->getCurrentThrustAzimuthAngle() << std::endl; }

    if( debugInfo == 1 ){std::cout << "    Repeating Mass Variables" << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Current mass = " << bodyMap_.at( vehicleName )->getBodyMass() << std::endl; }
    if( debugInfo == 1 ){std::cout << "         Dry mass     = " << vehicleSystems->getDryMass() << std::endl; }

    if( debugInfo == 1 ){std::cout << "    Evaluating Engine Status Function" << std::endl; }
    bislipSystems->setCurrentEngineStatus( bislip::Variables::determineEngineStatus( bodyMap_.at( vehicleName )->getBodyMass(), vehicleSystems->getDryMass() ) );
    if( debugInfo == 1 ){std::cout << "         Initial Engine Status    = " << bislipSystems->getCurrentEngineStatus() << std::endl; }

    if( debugInfo == 1 ){std::cout << "    Evaluating Throttle Setting Interpolator" << std::endl; }
    bislipSystems->setCurrentThrottleSetting( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::ThrottleSetting, bodyMap_, vehicleName ) );
    if( debugInfo == 1 ){std::cout << "         Initial Throttle Setting = " << bislipSystems->getCurrentThrottleSetting() << std::endl; }

    if( debugInfo == 1 ){std::cout << "    Evaluating Thrust Magnitude Function" << std::endl; }
    bislipSystems->setCurrentThrustMagnitude( bislip::Variables::computeThrustMagnitude( bodyMap_, vehicleName ) );
    if( debugInfo == 1 ){std::cout << "         Initial Thrust Magnitude = " << bislipSystems->getCurrentThrustMagnitude() << std::endl; }

    if( debugInfo == 1 ){std::cout << "    Evaluating Thrust Elevation Angle Interpolator" << std::endl; }
    bislipSystems->setCurrentThrustElevationAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::ThrustElevationAngle, bodyMap_, vehicleName ) ) );
    if( debugInfo == 1 ){std::cout << "         Initial Thrust Elevation Angle = " << bislipSystems->getCurrentThrustElevationAngle() << std::endl; }

    if( debugInfo == 1 ){std::cout << "    Evaluating Body-Fixed Thrust Direction Function" << std::endl; }
    bislipSystems->setCurrentBodyFixedThrustDirection( bislip::Variables::computeBodyFixedThrustDirection( bodyMap_, vehicleName ) );
    if( debugInfo == 1 ){std::cout << "         Initial Thrust Direction = [ " << ( bislipSystems->getCurrentBodyFixedThrustDirection() )( 0 ) << " , " << ( bislipSystems->getCurrentBodyFixedThrustDirection() )( 1 ) << " , " << ( bislipSystems->getCurrentBodyFixedThrustDirection() )( 2 ) << " ]" << std::endl; }
*/
    bislipSystems->setCurrentThrustElevationAngle( 0.0 );
    bislipSystems->setCurrentThrustAzimuthAngle( 0.0 );
    bislipSystems->setCurrentBankAngle( 0.0 );

    if( debugInfo == 1 ){ std::cout << "    Setting Initial Control Surface Deflections" << std::endl; }
    vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", bislipSystems->getCurrentBodyFlapAngle() );
    vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonLeft", bislipSystems->getCurrentElevonAngle() );
    vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonRight", bislipSystems->getCurrentElevonAngle() );
    if( debugInfo == 1 ){ std::cout << "         Initial BodyFlap Deflection = " << vehicleSystems->getCurrentControlSurfaceDeflection( "BodyFlap" ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "         Initial Elevon Deflection   = " << vehicleSystems->getCurrentControlSurfaceDeflection( "ElevonLeft" ) << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Calculate Full Current Coefficients" << std::endl; }
    bislipSystems->setFullCurrentCoefficients( bislip::Variables::computeFullCurrentCoefficients( bodyMap_, vehicleName ) );
    if( debugInfo == 1 ){ std::cout << "         Preliminary Coefficient Vector = [ " << bislipSystems->getFullCurrentCoefficients()( 0 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 1 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 2 ) << bislipSystems->getFullCurrentCoefficients()( 3 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 4 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 5 ) << " ]" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Calculate Preliminary Aerodynamic Frame Aerodyamic Load" << std::endl; }
    Eigen::Vector3d aerodynamicFrameAerodynamicLoadVector = bislip::Variables::computeAerodynamicFrameAerodynamicLoad( bodyMap_, vehicleName );
    if( debugInfo == 1 ){ std::cout << "         Preliminary Aerodyamic Load = [ " << aerodynamicFrameAerodynamicLoadVector( 0 ) << ", " << aerodynamicFrameAerodynamicLoadVector( 1 ) << ", " << aerodynamicFrameAerodynamicLoadVector( 2 ) << " ]" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Set Preliminary Drag Force" << std::endl; }
    bislipSystems->setCurrentDragForce( -aerodynamicFrameAerodynamicLoadVector( 0 ) );
    if( debugInfo == 1 ){ std::cout << "         Preliminary Drag Force = " << bislipSystems->getCurrentDragForce() << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Set Preliminary Lift Force" << std::endl; }
    bislipSystems->setCurrentLiftForce( -aerodynamicFrameAerodynamicLoadVector( 2 ) );
    if( debugInfo == 1 ){ std::cout << "         Preliminary Lift Force = " << bislipSystems->getCurrentLiftForce() << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Evaluating Local Gravity Vector Function" << std::endl; }
    bislipSystems->setCurrentLocalGravityVector( bislip::Variables::computeLocalGravity( bodyMap_, vehicleName, centralBodyName ) );
    if( debugInfo == 1 ){ std::cout << "         Initial Local Gravity Vector = [ " << bislipSystems->getCurrentLocalGravityVector()( 0 ) << ", " << bislipSystems->getCurrentLocalGravityVector()( 1 ) << " ]" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Create caller to Evaluate Guidance Functions" << std::endl; }
    bislip::MyGuidance initialEvaluation( bodyMap_, vehicleName, centralBodyName, simulationStartEpoch );

    if( debugInfo == 1 ){ std::cout << "    Evaluating All Guidance Functions to re-Determine All Relevant Initial Values" << std::endl; }
    if( bislipSystems->getValidationFlag() == true )
    { initialEvaluation.evaluateValidationGuidanceFunctions( bislipSystems, vehicleSystems, "Descent", 0.0 ); }
    else { initialEvaluation.evaluateGuidanceFunctions( bislipSystems, vehicleSystems, "Ascent", 0.0 ); }

    if( debugInfo == 1 ){ std::cout << "    Calculate Full Current Coefficients" << std::endl; }
    bislipSystems->setFullCurrentCoefficients( bislip::Variables::computeFullCurrentCoefficients( bodyMap_, vehicleName ) );
    if( debugInfo == 1 ){ std::cout << "         Initial Coefficient Vector = [ " << bislipSystems->getFullCurrentCoefficients()( 0 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 1 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 2 ) << bislipSystems->getFullCurrentCoefficients()( 3 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 4 ) << ", " << bislipSystems->getFullCurrentCoefficients()( 5 ) << " ]" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Calculate Initial Aerodynamic Frame Aerodyamic Load" << std::endl; }
    aerodynamicFrameAerodynamicLoadVector = bislip::Variables::computeAerodynamicFrameAerodynamicLoad( bodyMap_, vehicleName );
    if( debugInfo == 1 ){ std::cout << "         Initial Aerodyamic Load = [ " << aerodynamicFrameAerodynamicLoadVector( 0 ) << ", " << aerodynamicFrameAerodynamicLoadVector( 1 ) << ", " << aerodynamicFrameAerodynamicLoadVector( 2 ) << " ]" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Set Initial Drag Force" << std::endl; }
    bislipSystems->setCurrentDragForce( -aerodynamicFrameAerodynamicLoadVector( 0 ) );
    if( debugInfo == 1 ){std::cout << "         Initial Drag Force = " << bislipSystems->getCurrentDragForce() << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Set Initial Lift Force" << std::endl; }
    bislipSystems->setCurrentLiftForce( -aerodynamicFrameAerodynamicLoadVector( 2 ) );
    if( debugInfo == 1 ){std::cout << "         Initial Lift Force = " << bislipSystems->getCurrentLiftForce() << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Set Initial Flight-Path Angle Rate" << std::endl; }
    bislipSystems->setCurrentFlightPathAngleRate( bislip::Variables::computeFlightPathAngleRate( bodyMap_, vehicleName, centralBodyName ) );
    if( debugInfo == 1 ){std::cout << "         Initial Flight-Path Angle Rate = " << bislipSystems->getCurrentLiftForce() << std::endl; }



    bislipSystems->setInitialValueFlag( false );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE MASS RATE SETTINGS: ASCENT            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){std::cout << "Create mass rate models: Ascent" << std::endl;}

    //! Declare and initialize mass rate model settings.
    std::shared_ptr< MassRateModelSettings > massRateModelSettings =
            std::make_shared< FromThrustMassModelSettings >( true );

    //! Declare and initialize mass rate model settings.
    std::map< std::string, std::shared_ptr< MassRateModel > > massRateModels;
    massRateModels[  vehicleName ] = createMassRateModel(
                vehicleName, massRateModelSettings, bodyMap_, problemInput_->getAccelerationMap() );

    //! Create settings for propagating the mass of the vehicle.
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back(  vehicleName );

    //! Declare and initialize starting mass of vehicle.
    Eigen::VectorXd initialBodyMasses_Ascent( 1 );
    initialBodyMasses_Ascent( 0 ) = initialMass_Ascent + additionalMass;

    //! Declare and initialize mass propagation settings.
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings_Ascent =
            std::make_shared< MassPropagatorSettings< double > >(
                bodiesWithMassToPropagate, massRateModels, initialBodyMasses_Ascent, problemInput_->getAscentTerminationSettings(), problemInput_->getDependentVariablesToSave() );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE PROPAGATION SETTINGS: ASCENT              ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){std::cout << "Create Propagation Settings: Ascent" << std::endl;}

    //! Create translational propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings_Ascent =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( problemInput_->getCentralBodies(),
              problemInput_->getAccelerationMap(),
              problemInput_->getBodiesToIntegrate(),
              systemInitialState_Ascent,
              problemInput_->getAscentTerminationSettings(),
              cowell,
              problemInput_->getDependentVariablesToSave() );

    //! Declare and initialize list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector_Ascent;
    propagatorSettingsVector_Ascent.push_back( translationalPropagatorSettings_Ascent );
    propagatorSettingsVector_Ascent.push_back( massPropagatorSettings_Ascent );

    //! Declare and initialize propagation settings for both mass and translational dynamics.
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings_Ascent =
            std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector_Ascent, problemInput_->getAscentTerminationSettings(), problemInput_->getDependentVariablesToSave() );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE INTEGRATION SETTINGS: ASCENT           ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Create Integration Settings: Ascent" << std::endl;}

    //! Declare and initialize integrator settings.
    std::shared_ptr< IntegratorSettings<  > > integratorSettings_Ascent =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationStartEpoch,
              propagationStepSize, 1, false );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          PROPAGATE TRAJECTORY: ASCENT         /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    if( debugInfo == 1 ){ std::cout << "Starting Ascent propagation" << std::endl; }

    //! Propagate trajectory.
    SingleArcDynamicsSimulator< double > ascentSimulation(
                bodyMap_,
                integratorSettings_Ascent,
                propagatorSettings_Ascent );

    if( debugInfo == 1 ){ std::cout << "Ascent Propagation done" << std::endl; }


    //! Retrieve propagation and dependent variable maps - Ascent
    propagationTimeHistoryMap_Ascent = ascentSimulation.getEquationsOfMotionNumericalSolution( );
    depVarTimeHistoryMap_Ascent = ascentSimulation.getDependentVariableHistory( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////          RETRIEVE ASCENT-DESCENT PROPAGATION LINKING DATA          //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Retrieve final Ascent epoch
    simulationEndEpoch_Ascent = ( ascentSimulation.getEquationsOfMotionNumericalSolution().rbegin() )->first;

    if( debugInfo == 1 ){ std::cout << "Ascent Phase Initial Epoch = " << simulationStartEpoch << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Ascent Phase Final Epoch   = " << simulationEndEpoch_Ascent << std::endl; }

    //! Retrieve final Ascent state, equal to initial Descent state
    systemInitialState_Descent = ( ascentSimulation.getEquationsOfMotionNumericalSolution().rbegin() )->second;

    if( debugInfo == 1 ){ std::cout << "Ascent Phase Final State = " << systemInitialState_Descent << std::endl; }

    //! Retrieve dependent variables of final Ascent state
    dependentVariableFinalState_Ascent = ( ascentSimulation.getDependentVariableHistory( ).rbegin() )->second;

    if( debugInfo == 1 ){ std::cout << "Ascent Phase Final DepVar = " << dependentVariableFinalState_Ascent << std::endl; }

    //! Retrieve final Ascent mass, equal to initial Descent mass
    finalAirspeed_Ascent       = dependentVariableFinalState_Ascent[ 13 ];
    finalMass_Ascent           = dependentVariableFinalState_Ascent[ 27 ];
    initialMass_Descent        = finalMass_Ascent;
    finalSpecificEnergy_Ascent = dependentVariableFinalState_Ascent[ 28 ];

    //  }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            COMPLETE KNOWN DESCENT STATES            ///////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( debugInfo == 1 ){ std::cout << "Complete Known Descent States" << std::endl; }

    /*
    if( E_i > 2 * E_f )
    {

        if( debugInfo == 1 ){ std::cout << "Initial Specific Energy is HIGHER than Final Specific Energy" << std::endl; }

        bislipSystems->setCurrentThrustElevationAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustElevationAngle, bodyMap_, vehicleName ) ) );
        if( debugInfo == 1 ){std::cout << "     Initial Thrust Elevation Angel = " << bislipSystems->getCurrentThrustElevationAngle() << std::endl; }

        bislipSystems->setCurrentThrustAzimuthAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustAzimuthAngle, bodyMap_, vehicleName ) ) );
        if( debugInfo == 1 ){std::cout << "     Initial Thrust Azimuth Angel   = " << bislipSystems->getCurrentThrustAzimuthAngle() << std::endl; }

        bislipSystems->setCurrentThrottleSetting( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrottleSetting, bodyMap_, vehicleName ) );
        if( debugInfo == 1 ){std::cout << "     Initial Throttle Setting       = " << bislipSystems->getCurrentThrottleSetting() << std::endl; }

        if( debugInfo == 1 ){std::cout << "     Current mass = " << bodyMap_.at( vehicleName )->getBodyMass() << std::endl; }
        if( debugInfo == 1 ){std::cout << "     Dry mass     = " << vehicleSystems->getDryMass() << std::endl; }

        bislipSystems->setCurrentEngineStatus( bislip::Variables::determineEngineStatus( bodyMap_.at( vehicleName )->getBodyMass(), vehicleSystems->getDryMass() ) );
        if( debugInfo == 1 ){std::cout << "     Initial Engine Status    = " << bislipSystems->getCurrentEngineStatus() << std::endl; }

        bislipSystems->setCurrentBodyFixedThrustDirection( bislip::Variables::computeBodyFixedThrustDirection( bodyMap_, vehicleName ) );
        if( debugInfo == 1 ){std::cout << "     Initial Thrust Direction = " << bislipSystems->getCurrentBodyFixedThrustDirection() << std::endl; }

        bislipSystems->setCurrentThrustMagnitude( bislip::Variables::computeThrustMagnitude( bodyMap_, vehicleName ) );
        if( debugInfo == 1 ){std::cout << "     Initial Thrust Magnitude = " << bislipSystems->getCurrentThrustMagnitude() << std::endl; }



        systemInitialState_Descent = systemInitialState_Ascent;
        initialMass_Descent = initialMass_Ascent + additionalMass;
    }


*/





    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE PROPAGATION SETTINGS: DESCENT              ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Create Propagation Settings - Descent" << std::endl; }

    //! Create translational propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings_Descent =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( problemInput_->getCentralBodies(),
              problemInput_->getAccelerationMap(),
              problemInput_->getBodiesToIntegrate(),
              systemInitialState_Descent,
              problemInput_->getDescentTerminationSettings(),
              cowell,
              problemInput_->getDependentVariablesToSave() );

    //! Declare and initialize starting mass of vehicle.
    Eigen::VectorXd initialBodyMasses_Descent( 1 );
    initialBodyMasses_Descent( 0 ) = initialMass_Descent;

    if( debugInfo == 1 ){ std::cout << "Initial Ascent Mass  = " << initialMass_Ascent << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Initial Descent Mass = " << initialMass_Descent << std::endl; }

    /*
    //! Create settings for propagating the mass of the vehicle.
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back(  vehicleName );

    //! Declare and initialize mass rate model settings.
    std::shared_ptr< MassRateModelSettings > massRateModelSettings =
            std::make_shared< FromThrustMassModelSettings >( true );

    //! Declare and initialize mass rate model settings.
    std::map< std::string, std::shared_ptr< MassRateModel > > massRateModels;
    massRateModels[  vehicleName ] = createMassRateModel(
                vehicleName, massRateModelSettings, bodyMap_, problemInput_->getAccelerationMap() );
*/
    //! Declare and initialize mass propagation settings.
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings_Descent =
            std::make_shared< MassPropagatorSettings< double > >(
                bodiesWithMassToPropagate, massRateModels, initialBodyMasses_Descent, problemInput_->getDescentTerminationSettings(), problemInput_->getDependentVariablesToSave() );

    //! Declare and initialize list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector_Descent;
    propagatorSettingsVector_Descent.push_back( translationalPropagatorSettings_Descent );
    propagatorSettingsVector_Descent.push_back( massPropagatorSettings_Descent );

    //! Declare and initialize propagation settings for both mass and translational dynamics.
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings_Descent =
            std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector_Descent, problemInput_->getDescentTerminationSettings(), problemInput_->getDependentVariablesToSave() );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE INTEGRATION SETTINGS: DESCENT           //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Create Integration Settings: Descent" << std::endl; }

    std::shared_ptr< IntegratorSettings<  > > integratorSettings_Descent =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationEndEpoch_Ascent,
              propagationStepSize, 1, false );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           ASSIGN INTERPOLATORS & BOUNDS: DESCENT           ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Pass Descent Interpolators to vehicle systems" << std::endl; }

    //! Pass Descent phase interpolators to vehicle systems.
    bislipSystems->setParameterInterpolator( Interpolators_Descent );

    if( debugInfo == 1 ){ std::cout << "Pass Descent Bounds to vehicle systems" << std::endl; }

    //! Pass parameter bounds to vehicle systems.
    //bislipSystems->setParameterBounds( problemInput_->getDescentParameterBoundsMap()  );
    bislipSystems->setParameterLowerBounds( Interpolators_Descent_LB);
    bislipSystems->setParameterUpperBounds( Interpolators_Descent_UB );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          PROPAGATE TRAJECTORY: DESCENT         ////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Starting Descent propagation" << std::endl; }

    bislipSystems->setCurrentTrajectoryPhase( "Descent" );

    //! Propagate trajectory.
    SingleArcDynamicsSimulator< double > descentSimulation(
                bodyMap_,
                integratorSettings_Descent,
                propagatorSettings_Descent );

    if( debugInfo == 1  ){ std::cout << "Descent Propagation done" << std::endl; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          RETRIEVE DATA MAPS          //////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( debugInfo == 1 ){ std::cout << "Retrieve data maps" << std::endl; }


    //! Retrieve propagation and dependent variable maps - Descent
    const std::map< double, Eigen::VectorXd > propagationTimeHistoryMap_Descent = descentSimulation.getEquationsOfMotionNumericalSolution( );
    const std::map< double, Eigen::VectorXd > depVarTimeHistoryMap_Descent = descentSimulation.getDependentVariableHistory( );

    if( debugInfo == 1 ){ std::cout << "Merge Propagation History Ascent & Descent data maps" << std::endl; }

    //! Merge propagation data maps
    std::map< double, Eigen::VectorXd > propagationTimeHistoryMap = propagationTimeHistoryMap_Ascent;
    //propagationTimeHistoryMap.insert( propagationTimeHistoryMap_Ascent.begin(), propagationTimeHistoryMap_Ascent.end() );
    propagationTimeHistoryMap.insert( propagationTimeHistoryMap_Descent.begin(), propagationTimeHistoryMap_Descent.end() );

    if( debugInfo == 1 ){ std::cout << "Merge Dependent Variable Ascent & Descent data maps" << std::endl; }

    //! Merge dependent variable data maps
    std::map< double, Eigen::VectorXd > depVarTimeHistoryMap = depVarTimeHistoryMap_Ascent;
    //depVarTimeHistoryMap.insert( depVarTimeHistoryMap_Ascent.begin(), depVarTimeHistoryMap_Ascent.end() );
    depVarTimeHistoryMap.insert( depVarTimeHistoryMap_Descent.begin(), depVarTimeHistoryMap_Descent.end() );





    /*
    std::cout << "depVarTimeHistoryMap_Ascent.size(): " << depVarTimeHistoryMap_Ascent.size() << std::endl;
    std::cout << "depVarTimeHistoryMap_Descent.size(): " << depVarTimeHistoryMap_Descent.size() << std::endl;
    std::cout << "depVarTimeHistoryMap.size(): " << depVarTimeHistoryMap.size() << std::endl;

    std::cout << "( depVarTimeHistoryMap_Ascent.begin() )->first: " << ( depVarTimeHistoryMap_Ascent.begin() )->first << std::endl;
    std::cout << "( depVarTimeHistoryMap_Ascent.end() )->first: " << ( depVarTimeHistoryMap_Ascent.end() )->first << std::endl;
    std::cout << "( depVarTimeHistoryMap_Ascent.rbegin() )->first: " << ( depVarTimeHistoryMap_Ascent.rbegin() )->first << std::endl;

    std::cout << "( depVarTimeHistoryMap_Descent.begin() )->first: " << ( depVarTimeHistoryMap_Descent.begin() )->first << std::endl;
    std::cout << "( depVarTimeHistoryMap_Descent.end() )->first: " << ( depVarTimeHistoryMap_Descent.end() )->first << std::endl;
    std::cout << "( depVarTimeHistoryMap_Descent.rbegin() )->first: " << ( depVarTimeHistoryMap_Descent.rbegin() )->first << std::endl;

*/
    if( debugInfo == 1 ){ std::cout << "Determine size of time vector" << std::endl; }

    //! Extract time vector ( Time History keys )
    int timeVectorSize = 0;
    for( std::map< double, Eigen::VectorXd >::iterator it = depVarTimeHistoryMap.begin(); it != depVarTimeHistoryMap.end(); ++it)
    { timeVectorSize += 1; }

    if( debugInfo == 1 ){ std::cout << "Extract Data Map keys and place time values into time vector" << std::endl; }

    //  std::vector<std::string> extract_keys(std::map<std::string, std::string> const& input_map) {
    std::vector< double > timeVector;
    for ( auto const& element : depVarTimeHistoryMap ) { timeVector.push_back( element.first ); }


    // if( debugInfo == 1 ){ for ( unsigned long i = 0; i < timeVector.size(); i++ ) { std::cout << "timeVector[ " << i << " ] = " << timeVector[ i ] << std::endl; } }




    //      return retval;
    // }



    //Eigen::VectorXd timeVector( timeVectorSize );
    //int p = 0;
    //for( std::map< double, Eigen::VectorXd >::iterator it = depVarTimeHistoryMap.begin(); it != depVarTimeHistoryMap.end(); ++it)
    // for ( int i = 1; i < retval.size(); i++ )
    // { timeVector( p ) = it->first; p += 1; std::cout << "timeVector( p ):  " << timeVector( p ) << std::endl; }


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
    //double tof = problemInput_->getSimulationSettings()[ 1 ];

    //double maximum_DynamicPressure = constraint_DynamicPressure;
    //double maximum_MechanicalLoad = constraint_MechanicalLoad;
    //double finalMass_Descent = 0;
    const double Isp = bislipSystems->getSpecificImpulse();
    const double theoreticaldelV_Ascent  = g0 * Isp * log( bislipSystems->getInitialMass() / bodyMap_.at(  vehicleName )->getVehicleSystems()->getDryMass() );



    if( initialMass_Descent < bodyMap_.at(  vehicleName )->getVehicleSystems()->getDryMass() ) { initialMass_Descent = bodyMap_.at(  vehicleName )->getVehicleSystems()->getDryMass(); }
    const double theoreticaldelV_Descent = g0 * Isp * log( initialMass_Descent / bodyMap_.at(  vehicleName )->getVehicleSystems()->getDryMass() );

    //! Retrieve results
    //! Maybe this command would allow the use of certain results to guide
    //! further propagations within the same run. Sort of stitching them together.
    // std::map< double, Eigen::VectorXd > integrationResult =
    //        dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    //! Extract Final Epoch
    //! Trying to find the actual END time of the simulation via key extraction.
    //double simulationEndEpoch_calc = ( --dynamicsSimulator.getEquationsOfMotionNumericalSolution().end() )->first;

    //! Extracting the contents of the last timestep via extracted key
    //const Eigen::Vector6d systemFinalState = integrationResult[simulationEndEpoch_calc];

    //! Told by Dominic that this here gives the final epoch directly
    //const double simulationEndEpoch_calc =  ( dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin() )->first;

    //! Told by Dominic that this here gives the final state directly
    //Eigen::Vector6d systemFinalState =  ( dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin() )->second;

    //! Calculate Time of Flight
    const double timeOfFlight_Ascent = simulationEndEpoch_Ascent - simulationStartEpoch;
    const double timeOfFlight_Descent = simulationEndEpoch_Descent - simulationEndEpoch_Ascent;
    //const double timeOfFlight = simulationEndEpoch_Descent - simulationStartEpoch;

    //! Transform final state from Inertial Frame to Earth-Fixed Frame
    //Eigen::Vector6d systemFinalState_EARTobjectiveHeight_AscentIXED = transformStateToTargetFrame( systemFinalState, simulationEndEpoch_calc, earthRotationalEphemeris );

    //! Transform expected final state to Inertial frame. This is where the
    //! destination is located at the termination of the simulation.
    //! Transformation is required because coordinates are Earth-fixed, similar
    //! to input values.
    //systemFinalStateGOAL = transformStateToGlobalFrame(systemFinalStateGOAL,simulationEndEpoch_calc,earthRotationalEphemeris );

    //! Calculate latitude and longitude of the GOAL state: Inertial Frame.
    //const double altitude_f_GOAL_calc = std::sqrt( pow(systemFinalStateGOAL[0],2) +
    //       pow(systemFinalStateGOAL[1],2) + pow(systemFinalStateGOAL[2],2) ) ;
    //double lon_f_rad_GOAL_calc = std::atan2(systemFinalStateGOAL[1] , systemFinalStateGOAL[0]);
    //double lat_f_rad_GOAL_calc = std::asin(systemFinalStateGOAL[2] / altitude_f_GOAL_calc) ;

    //! Convert lat/lon of GOAL state to degrees: Inertial Frame
    //const double lon_f_deg_GOAL_calc = lon_f_rad_GOAL_calc * 180 / mathematical_constants::PI;
    //const double lat_f_deg_GOAL_calc = lat_f_rad_GOAL_calc * 180 / mathematical_constants::PI;



    //! Extract Dependent variables of final state.
    // const std::map< double, Eigen::VectorXd > depVarTimeHistoryMap_Ascent = dynamicsSimulator.getDependentVariableHistory( );
    // std::cout << "depVar_FINALSTATE: " << depVar_FINALSTATE << std::endl;

    //! Declare number of elements in time history.
    const int rowsTotal = depVarTimeHistoryMap.size();
    const int rowsAscent = depVarTimeHistoryMap_Ascent.size();
    // long rowsDescent = depVarTimeHistoryMap_Descent.size();
    const int columns = ( ( descentSimulation.getDependentVariableHistory( ).begin() )->second ).size();


    if( debugInfo == 1 ){ std::cout << "Convert Dependent Variable Map to Matrix" << std::endl; }
    Eigen::MatrixXd depVarTimeHistoryMatrix( depVarTimeHistoryMap.size(), columns );
    depVarTimeHistoryMatrix = Eigen::MatrixXd::Zero( depVarTimeHistoryMap.size(), columns );

    if( debugInfo == 1 ){ std::cout << "depVarTimeHistoryMap.size():  " << depVarTimeHistoryMap.size() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "depVarTimeHistoryMatrix.size():  " << depVarTimeHistoryMatrix.size() << std::endl; }



    for ( unsigned long i = 0; i < depVarTimeHistoryMap.size(); i++ )
    {
        //  if( debugInfo == 1 ){ std::cout << "timeVector[ " << i << " ] = " << timeVector[ i ] << std::endl;  }
        // if( debugInfo == 1 ){ std::cout << "( depVarTimeHistoryMap.begin() )->first: " << ( depVarTimeHistoryMap )->first << std::endl; }


        depVarTimeHistoryMatrix.row( i ) = depVarTimeHistoryMap.at( timeVector[ i ] );

    }
    // std::cout << "depVarTimeHistoryMatrix.size():  " << depVarTimeHistoryMatrix.size() << std::endl;

    //std::cout << "Aqui!!!!" << std::endl;

    //  rows = depVarTimeHistoryMap_Descent.size();
    //  columns = ( ( descentSimulation.getDependentVariableHistory( ).begin() )->second ).size();

    if( debugInfo == 1 ){ std::cout << "Extract various variable time histories." << std::endl; }

    Eigen::VectorXd depVar_BodyFixedX                  = depVarTimeHistoryMatrix.col( 0 );
    Eigen::VectorXd depVar_BodyFixedY                  = depVarTimeHistoryMatrix.col( 1 );
    Eigen::VectorXd depVar_BodyFixedZ                  = depVarTimeHistoryMatrix.col( 2 );
    Eigen::VectorXd depVar_FlightPathAngle             = depVarTimeHistoryMatrix.col( 7 );
    Eigen::VectorXd depVar_BankAngle                   = depVarTimeHistoryMatrix.col( 10 );
    Eigen::VectorXd depVar_Height                      = depVarTimeHistoryMatrix.col( 11 );
    Eigen::VectorXd depVar_Airspeed                    = depVarTimeHistoryMatrix.col( 13 );
    Eigen::VectorXd depVar_DynamicPressure             = depVarTimeHistoryMatrix.col( 24 );
    Eigen::VectorXd depVar_CurrentMass                 = depVarTimeHistoryMatrix.col( 27 );
    Eigen::VectorXd depVar_SpecificEnergy              = depVarTimeHistoryMatrix.col( 28 );
    Eigen::VectorXd depVar_NormalizedSpecificEnergy    = depVarTimeHistoryMatrix.col( 29 );
    Eigen::VectorXd depVar_AngularDistanceTravelled    = depVarTimeHistoryMatrix.col( 36 );
    Eigen::VectorXd depVar_AngularDistanceToGo         = depVarTimeHistoryMatrix.col( 37 );
    Eigen::VectorXd depVar_HeadingToTarget             = depVarTimeHistoryMatrix.col( 38 );
    Eigen::VectorXd depVar_HeadingError                = depVarTimeHistoryMatrix.col( 39 );
    Eigen::VectorXd depVar_BendingMoment               = depVarTimeHistoryMatrix.col( 44 );
    Eigen::VectorXd depVar_TotalBodyFixed_Z_Component  = depVarTimeHistoryMatrix.col( 56 );
    Eigen::VectorXd depVar_BodyFixedMechanicalLoadMag  = depVarTimeHistoryMatrix.col( 57 );
    Eigen::VectorXd depVar_PitchMomentCoefficient      = depVarTimeHistoryMatrix.col( 66 );
    Eigen::VectorXd depVar_HeatFluxChapman             = depVarTimeHistoryMatrix.col( 67 );
    Eigen::VectorXd depVar_CurrentHeadingErrorDeadband = depVarTimeHistoryMatrix.col( 78 );
    Eigen::VectorXd depVar_BankReversalTrigger         = depVarTimeHistoryMatrix.col( 81 );
    Eigen::VectorXd depVar_GroundtrackCovered          = depVarTimeHistoryMatrix.col( 87 );
    Eigen::VectorXd depVar_GroundtrackDifference       = depVarTimeHistoryMatrix.col( 88 );
    Eigen::VectorXd depVar_TimeOfFlight                = depVarTimeHistoryMatrix.col( 89 );
    Eigen::VectorXd depVar_FlightPathAngleRate         = depVarTimeHistoryMatrix.col( 90 );
    Eigen::VectorXd depVar_CumulativeDistanceTravelled = depVarTimeHistoryMatrix.col( 91 );


    if( debugInfo == 1 ){ std::cout << "  " << std::endl; }
    if( debugInfo == 1 ){ std::cout << "  " << std::endl; }
    if( debugInfo == 1 ){ std::cout << "angularDistanceToGo_TerminationSettings: " << depVar_AngularDistanceToGo( depVar_AngularDistanceToGo.size() - 1 ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "netAngularDistanceTravelled_TerminationSettings: " << depVar_AngularDistanceTravelled( depVar_AngularDistanceTravelled.size() - 1 ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "minimumHeightAllowable_TerminationSettings: " << depVar_Height( depVar_Height.size() - 1 ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "maximumTimeOfFlight_TerminationSettings: " << depVar_TimeOfFlight( depVar_TimeOfFlight.size() - 1 ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "  " << std::endl; }
    if( debugInfo == 1 ){ std::cout << "  " << std::endl; }



    //! Extract Dependent variables of final state.
    const Eigen::VectorXd dependentVariableFinalState_Descent = ( descentSimulation.getDependentVariableHistory( ).rbegin() )->second;

    //const double altitude_f_calc = depVar_FINALSTATE[3];
    const double targetLat_rad_calc  = depVarTimeHistoryMatrix( depVarTimeHistoryMatrix.rows() - 1, 4 );
    const double targetLon_rad_calc  = depVarTimeHistoryMatrix( depVarTimeHistoryMatrix.rows() - 1, 5 );
    //const double E_hat_f_calc        = depVarTimeHistoryMatrix( depVarTimeHistoryMatrix.rows() - 1, 29 );
    std::ptrdiff_t ind_maximum_height;
    objectiveHeight_Ascent_calc      = depVar_Height.maxCoeff( &ind_maximum_height );
    objectiveAirspeed_Ascent_calc    = depVar_Airspeed( ind_maximum_height );
    objectiveHeight_Ascent_calc      = depVarTimeHistoryMatrix( depVarTimeHistoryMatrix.rows() - 1, 11 );
    objectiveAirspeed_Descent_calc   = depVarTimeHistoryMatrix( depVarTimeHistoryMatrix.rows() - 1, 13 );

    std::ptrdiff_t index_MaximumNormalizedSpecificEnergy;
    std::ptrdiff_t index_MaximumAirspeed;
    std::ptrdiff_t index_MaximumMechanicalLoad;
    double maximum_NormalizedSpecificEnergy = E_hat_max;
    double maximumAirspeed = initialAirspeed_Ascent;
    double maximumMechanicalLoad = 0.0;

    maximum_NormalizedSpecificEnergy  = depVar_NormalizedSpecificEnergy.maxCoeff( &index_MaximumNormalizedSpecificEnergy );
    maximumAirspeed = depVar_Airspeed.maxCoeff( &index_MaximumAirspeed );
    maximumMechanicalLoad = depVar_BodyFixedMechanicalLoadMag.maxCoeff( &index_MaximumMechanicalLoad );

    const double actualdelV_Ascent  = maximumAirspeed - initialAirspeed_Ascent;
    const double actualdelV_Descent = std::abs( objectiveAirspeed_Descent - maximumAirspeed );
    const double finalMass_Descent  = dependentVariableFinalState_Descent[ 27 ];

    //std::ptrdiff_t index_MaximumDynamicPressure;
    //maximum_DynamicPressure               = depVar_DynamicPressure.maxCoeff( &index_MaximumDynamicPressure );


    //std::cout << "lat_f_rad_calc: " << lat_f_rad_calc << std::endl;
    //std::cout << "lon_f_rad_calc: " << lon_f_rad_calc << std::endl;
    //std::cout << "objectiveHeight_Ascent_calc     : " << objectiveHeight_Ascent_calc << std::endl;
    //std::cout << "objectiveHeight_Ascent_calc     : " << objectiveHeight_Ascent_calc << std::endl;

    //const double altitude_f_calc = std::sqrt( pow(systemFinalState_EARTobjectiveHeight_AscentIXED[0],2) +
    //       pow(systemFinalState_EARTobjectiveHeight_AscentIXED[1],2) + pow(systemFinalState_EARTobjectiveHeight_AscentIXED[2],2) ) ;
    //const double lon_f_rad_calc = std::atan2(systemFinalState_EARTobjectiveHeight_AscentIXED[1] , systemFinalState_EARTobjectiveHeight_AscentIXED[0]);
    //const double lat_f_rad_calc = std::asin(systemFinalState_EARTobjectiveHeight_AscentIXED[2] / altitude_f_calc) ;

    //! Convert coordinates of final state to degrees: Earth-Fixed Frame
    targetLat_deg_calc = unit_conversions::convertRadiansToDegrees( targetLat_rad_calc );
    targetLon_deg_calc = unit_conversions::convertRadiansToDegrees( targetLon_rad_calc );

    //! Calculate angular distance of final state from target coordinates.
    // const double d_rad = bislip::Variables::computeAngularDistance( lat_f_rad_calc, lon_f_rad_calc, vehicleSystems->getTargetLat(), vehicleSystems->getTargetLon() );


    const double timeOfFlight                          = dependentVariableFinalState_Descent[ 89 ];
    const double finalAirspeed_Descent                 = dependentVariableFinalState_Descent[ 13 ];
    const double finalSpecificEnergy_Descent           = dependentVariableFinalState_Descent[ 28 ];
    const double finalNormalizedSpecificEnergy_Descent = dependentVariableFinalState_Descent[ 29 ];
    const double finalAngularDistanceTravelled_deg     = dependentVariableFinalState_Descent[ 36 ];
    const double finalAngularDistanceToGo_deg          = dependentVariableFinalState_Descent[ 37 ];
    const double finalHeight                           = dependentVariableFinalState_Descent[ 11 ];
    const double finalMechanicalLoad                   = dependentVariableFinalState_Descent[ 57 ];
    //const double d_deg = unit_conversions::convertRadiansToDegrees( d_rad );

    //! Calculate offset of final state from GOAL state: Earth-Fixed Frame
    //const double dif_lat_rad = targetLat_rad - lat_f_rad_calc;
    //const double dif_lon_rad = targetLon_rad - lon_f_rad_calc;

    //! Calculate offsets of final state in degrees: Earth-Fixed Frame
    const double dif_lat_deg = targetLat_deg - targetLat_deg_calc;
    const double dif_lon_deg = targetLon_deg - targetLon_deg_calc;

    //! Calculate "norm" of offsets. This is an arbitrary function I have
    //! implemented to pass on as an 'objective function'. It relates the
    //! differences such that when minimizing the offsets, there is an additional
    //! unsigned value that always goes to zero. Most definitely unsure about how
    //!  'proper' it is, yet is what works for the current BALLISTIC case.
    const double dif_norm = std::sqrt( ( dif_lat_deg * dif_lat_deg ) + ( dif_lon_deg * dif_lon_deg ) );

    //! Calculate offset from maximum elevation.
    const double dif_objectiveHeight_Ascent = objectiveHeight_Ascent - objectiveHeight_Ascent_calc;

    const double dif_objectiveAirspeed_Ascent = abs( objectiveAirspeed_Ascent - objectiveAirspeed_Ascent_calc );

    //! Calculate offset from final elevation.
    //const double dif_objectiveHeight_Ascent = objectiveHeight_Ascent - objectiveHeight_Ascent_calc;
    const double dif_objectiveAirspeed_Descent = abs( objectiveAirspeed_Descent - objectiveAirspeed_Descent_calc );
    //std::cout << " rowsAscent: " <<  rowsAscent<< std::endl;
    //std::cout << " rowsTotal: " <<  rowsTotal<< std::endl;
    //std::cout << " rowsTotal - rowsAscent: " <<  rowsTotal - rowsAscent << std::endl;




    //const double penaltyFinalNodeMagnitudeAscent                  = 1000 * ( xn_Ascent( xn_Ascent.size() - 1 ) - 1 );
    //const double penaltyInitialNodeMagnitudeDescent               = 1000 * ( xn_Descent( 0 ) - 1 );
    //const double penaltyFinalNodeMagnitudeDescent                 = 1000 * abs( xn_Descent( xn_Descent.size() - 1 ) );

    //const double penaltyFinalMappedNormalizedSpecificEnergyAscent    = 1000 * ( E_mapped_Ascent( E_mapped_Ascent.size() - 1 ) - 1 );
    //const double penaltyInitialMappedNormalizedSpecificEnergyDescent = 1000 * ( E_mapped_Descent( 0 ) - 1 );
    //const double penaltyFinalMappedNormalizedSpecificEnergyDescent   = abs( E_mapped_Descent(  E_mapped_Descent.size() - 1 ) ) - E_mapped_Descent(  E_mapped_Descent.size() - 1 );
    //const double penaltyMappedNormalizedSpecificEnergyInterface      = abs( E_mapped_Descent( 0 ) - E_mapped_Ascent( E_mapped_Ascent.size() - 1 ) );


    //std::cout << "rowsAscent = " << rowsAscent << std::endl;
    //std::cout << "rowsTotal = " << rowsTotal << std::endl;
    //std::cout << "rowsTotal - rowsAscent = " << rowsTotal - rowsAscent  << std::endl;
    /*
    if( debugInfo == 1 ){ std::cout << "Penalty: No Flight" << std::endl; }
    const double penaltyNoFlight                                  = ( 1.0 - maximum_NormalizedSpecificEnergy ) * ( 1.0 - maximum_NormalizedSpecificEnergy );

    if( debugInfo == 1 ){ std::cout << "Cost: Cumulative Cartesian Distance Travelled" << std::endl; }
    const double costCumulativeDistanceTravelled = depVar_CumulativeDistanceTravelled( rowsTotal - 1 ) / ( radiusEarth * bislipSystems->getInitialDistanceToTarget( ) );
    // std::cout << "costCumulativeDistanceTravelled: " <<  costCumulativeDistanceTravelled << std::endl;

    if( debugInfo == 1 ){ std::cout << "Penalty: Angular Distance Travelled" << std::endl; }
    //if( debugInfo == 5 ){ std::cout << "depVar_AngularDistanceTravelled" << depVar_AngularDistanceTravelled << std::endl; }
    const double penaltyGroundtrack = bislip::Variables::computeSumOfEigenVectorXd( bislip::Variables::computeElementwiseDivisionOfEigenVectorXd( depVar_GroundtrackDifference, depVar_AngularDistanceTravelled ) ) / 180.0;
    // std::cout << "penaltyGroundtrack: " <<  penaltyGroundtrack << std::endl;

    if( debugInfo == 1 ){ std::cout << "Penalty: Monotonic - Approach" << std::endl; }
    const double penaltyMonotonicApproach = bislip::Variables::computeMonotonicPenalty( depVar_AngularDistanceToGo, "Decreasing" );
*/
    if( debugInfo == 1 ){ std::cout << "Penalty: Angular Distance To Go" << std::endl; }
    const double penaltyDistanceToGo = 1000 * depVar_AngularDistanceToGo( rowsTotal - 1 ) / initialDistanceToTarget_deg;
    //std::cout << "finalAngularDistanceToGo_deg: " <<  depVar_AngularDistanceToGo( rowsTotal - 1 )  << std::endl;
    //std::cout << "penaltyDistanceToGo: " <<  penaltyDistanceToGo << std::endl;
    /*
    if( debugInfo == 1 ){ std::cout << "Penalty: Bank Angle Reversals" << std::endl; }
    double penaltyExcessiveBankReversals = 100000.0 * bislip::Variables::computeSumOfEigenVectorXd( depVar_BankReversalTrigger );
*/


    double costMaximumNormalizedSpecificEnergyAscent = 0.0;
    double penaltyMonotonicEnergyStateAscent = 0.0;
    double penaltydelV_Ascent = 0.0;
    double costFuelMassAscent = 0.0;
    double penaltyFuelMassAscentTEMP = 0.0;
    double penaltyFuelMassAscent = 0.0;

    double costHeatLoadChapmanAscent = 0.0;
    double penaltyHeatFluxChapmanAscent = 0.0;
    double penaltyDynamicPressureAscent = 0.0;
    double penaltyMechanicalLoadAscent = 0.0;
    double penaltyBendingMomentAscent = 0.0;
    double penaltyFlightPathAngleAscent = 0.0;
    double penaltyFlightPathAngleRateAscent = 0.0;
    double penaltyPitchMomentCoefficientAscent = 0.0;
    double penaltyHeadingErrorAscent = 0.0;


    if( debugInfo == 1  ){ std::cout << "Cost: Maximum Normalized Specific Energy - Ascent" << std::endl; }
    costMaximumNormalizedSpecificEnergyAscent     = std::abs( E_mapped_Ascent( E_mapped_Ascent.size() - 1 ) - maximum_NormalizedSpecificEnergy );
    /*
    if( debugInfo == 1  ){ std::cout << "Penalty: Monotonic Energy-State - Ascent" << std::endl; }
    penaltyMonotonicEnergyStateAscent = bislip::Variables::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.head( rowsAscent ), "Increasing" );

    if( debugInfo == 1  ){ std::cout << "Penalty: fitness-V Ascent" << std::endl; }
    penaltydelV_Ascent                               = std::abs( goaldelV_Ascent - actualdelV_Ascent ) / theoreticaldelV_Ascent;

    if( debugInfo == 1  ){ std::cout << "Cost: Fuel Mass - Ascent" << std::endl; }
    costFuelMassAscent                           = std::abs( initialMass_Ascent - finalMass_Ascent ) / initialMass_Ascent;

    if( debugInfo == 1  ){ std::cout << "Penalty: Fuel Mass - Ascent" << std::endl; }
    penaltyFuelMassAscentTEMP = ( finalMass_Ascent / initialMass_Ascent ) * ( 1.0 - std::exp( ( -1.0 ) * ( finalAirspeed_Ascent - std::sqrt( 2 * ( finalSpecificEnergy_Ascent - g0 * objectiveHeight_Ascent ) ) ) / ( g0 * Isp ) ) );
    if( std::isnan( penaltyFuelMassAscentTEMP ) == true ) { penaltyFuelMassAscentTEMP = 1e10; }
    penaltyFuelMassAscent = penaltyFuelMassAscentTEMP;

*/

    if( debugInfo == 1  ){ std::cout << "Cost: Chapman Heat Load - Ascent" << std::endl; }
    costHeatLoadChapmanAscent = bislip::Variables::computeSumOfEigenVectorXd( depVar_HeatFluxChapman.head( rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
    if( debugInfo == 1  ){ std::cout << "     costHeatLoadChapmanAscent: " << costHeatLoadChapmanAscent << std::endl; }
    //std::cout << "costHeatLoadChapmanAscent: " << costHeatLoadChapmanAscent << std::endl;

    if( debugInfo == 1  ){ std::cout << "Penalty: Chapman Heat Flux - Ascent" << std::endl; }
    penaltyHeatFluxChapmanAscent = bislip::Variables::computeCompoundViolationPenalty( depVar_HeatFluxChapman.head( rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Ascent * constraint_ChapmanHeatFlux );
    //std::cout << "penaltyHeatFluxChapmanAscent: " << penaltyHeatFluxChapmanAscent << std::endl;

    if( debugInfo == 1  ){ std::cout << "Penalty: Dynamic Pressure - Ascent" << std::endl; }
    penaltyDynamicPressureAscent = bislip::Variables::computeCompoundViolationPenalty( depVar_DynamicPressure.head( rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Ascent * constraint_DynamicPressure );
    //std::cout << "penaltyDynamicPressureAscent: " << penaltyDynamicPressureAscent << std::endl;

    if( debugInfo == 1  ){ std::cout << "Penalty: Mechanical Load - Ascent" << std::endl; }
    penaltyMechanicalLoadAscent = bislip::Variables::computeCompoundViolationPenalty( depVar_BodyFixedMechanicalLoadMag.head( rowsAscent ), constraint_MechanicalLoad, propagationStepSize, timeOfFlight_Ascent * constraint_MechanicalLoad );
    //std::cout << " penaltyMechanicalLoadAscent:   " <<  penaltyMechanicalLoadAscent << std::endl;

    if( debugInfo == 1  ){ std::cout << "Penalty: Bending Moment - Ascent" << std::endl; }
    penaltyBendingMomentAscent = bislip::Variables::computeCompoundViolationPenalty( depVar_BendingMoment.head( rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Ascent * constraint_BendingMoment );
    //std::cout << "penaltyBendingMomentAscent: " << penaltyBendingMomentAscent << std::endl;


    /*

    if( debugInfo == 1  ){ std::cout << "Penalty: Flight Path - Ascent" << std::endl; }
    penaltyFlightPathAngleAscent = bislip::Variables::computeConstraintViolationPenalty( ( -1.0 ) * depVar_FlightPathAngle.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );

    if( debugInfo == 1  ){ std::cout << "Penalty: Flight Path Rate - Ascent" << std::endl; }
    penaltyFlightPathAngleRateAscent = bislip::Variables::computeConstraintViolationPenalty( ( -1.0 ) * depVar_FlightPathAngleRate.head( rowsAscent ), 0.0, propagationStepSize, 360.0 );

    if( debugInfo == 1  ){ std::cout << "Penalty: Pitch Moment Coefficient - Ascent" << std::endl; }
    penaltyPitchMomentCoefficientAscent = bislip::Variables::computeConstraintViolationPenalty( depVar_PitchMomentCoefficient.head( rowsAscent ).cwiseAbs(), 0.0, propagationStepSize, tof );
    //std::cout << " penaltyPitchMomentCoefficientAscent:    " <<  penaltyPitchMomentCoefficientAscent << std::endl;

    if( debugInfo == 1  ){ std::cout << "Penalty: Heading Error - Ascent" << std::endl; }
    penaltyHeadingErrorAscent = bislip::Variables::computeDeadbandViolationPenalty(
                depVar_HeadingError.head( rowsAscent ), depVar_CurrentHeadingErrorDeadband.head( rowsAscent ), propagationStepSize, tof );
    //std::cout << "penaltyHeadingErrorAscent: " << penaltyHeatFluxChapmanAscent << std::endl;
*/
    /*

    if( debugInfo == 1  ){ std::cout << "Penalty: Maximum Normalized Specific Energy - Descent" << std::endl; }
    const double costMaximumNormalizedSpecificEnergyDescent    = std::abs( finalNormalizedSpecificEnergy_Descent - E_mapped_Descent( E_mapped_Descent.size() - 1 ) );

    if( debugInfo == 1  ){ std::cout << "Penalty: Monotonic Energy-State Descent" << std::endl; }
    const double penaltyMonotonicEnergyStateDescent = bislip::Variables::computeMonotonicPenalty( depVar_NormalizedSpecificEnergy.tail( rowsTotal - rowsAscent ), "Decreasing" );

    if( debugInfo == 1  ){ std::cout << "Penalty: fitness-V Descent" << std::endl; }
    double penaltydelV_Descent                              = std::abs( goaldelV_Descent - actualdelV_Descent ) / theoreticaldelV_Descent;
    if( std::isinf( penaltydelV_Descent ) == true ) { penaltydelV_Descent = 0;}
    // std::cout << " penaltydelV_Descent:   " <<  penaltydelV_Descent << std::endl;

    if( debugInfo == 1  ){ std::cout << "Cost: Consumed Mass - Descent" << std::endl; }
    const double costFuelMassDescent                          = std::abs( finalMass_Ascent - finalMass_Descent ) / initialMass_Ascent;

    if( debugInfo == 1  ){ std::cout << "Penalty: Fuel Mass - Descent" << std::endl; }
    double penaltyFuelMassDescentTEMP = ( finalMass_Descent / initialMass_Descent ) * ( 1.0 - std::exp( ( -1.0 ) * ( finalAirspeed_Descent - std::sqrt( 2 * ( finalSpecificEnergy_Descent - g0 * objectiveHeight_Descent ) ) ) / ( g0 * Isp ) ) );
    if( std::isnan( penaltyFuelMassDescentTEMP ) == true ) { penaltyFuelMassAscentTEMP = 1e10; }
    const double penaltyFuelMassDescent = penaltyFuelMassDescentTEMP;

*/

    if( debugInfo == 1  ){ std::cout << "Cost: Chapman Heat Load - Descent" << std::endl; }
    const double costHeatLoadChapmanDescent = bislip::Variables::computeSumOfEigenVectorXd( depVar_HeatFluxChapman.tail( rowsTotal - rowsAscent ) * propagationStepSize ) / ( timeOfFlight_Descent * constraint_ChapmanHeatFlux );
    //std::cout << " costHeatLoadChapmanDescent:    " <<  costHeatLoadChapmanDescent << std::endl;

    if( debugInfo == 1  ){ std::cout << "Penalty: Chapman Heat Flux - Descent" << std::endl; }
    const double penaltyHeatFluxChapmanDescent = bislip::Variables::computeCompoundViolationPenalty( depVar_HeatFluxChapman.tail( rowsTotal - rowsAscent ), constraint_ChapmanHeatFlux, propagationStepSize, timeOfFlight_Descent * constraint_ChapmanHeatFlux );
    //std::cout << "penaltyHeatFluxChapmanDescent: " << penaltyHeatFluxChapmanDescent << std::endl;

    if( debugInfo == 1  ){ std::cout << "Penalty: Dynamic Pressure - Descent" << std::endl; }
    const double penaltyDynamicPressureDescent = bislip::Variables::computeCompoundViolationPenalty( depVar_DynamicPressure.tail( rowsTotal - rowsAscent ), constraint_DynamicPressure, propagationStepSize, timeOfFlight_Descent * constraint_DynamicPressure );
    //std::cout << "penaltyDynamicPressureDescent: " << penaltyDynamicPressureDescent << std::endl;

    if( debugInfo == 1  ){ std::cout << "Penalty: Mechanical Load - Descent" << std::endl; }
    const double penaltyMechanicalLoadDescent = bislip::Variables::computeCompoundViolationPenalty( depVar_BodyFixedMechanicalLoadMag.tail( rowsTotal - rowsAscent ), constraint_MechanicalLoad, propagationStepSize, timeOfFlight_Descent * constraint_MechanicalLoad );
    //std::cout << " penaltyMechanicalLoadDescent:  " <<  penaltyMechanicalLoadDescent << std::endl;

    if( debugInfo == 1  ){ std::cout << "Penalty: Bending Moment - Descent" << std::endl; }
    const double penaltyBendingMomentDescent =  bislip::Variables::computeCompoundViolationPenalty( depVar_BendingMoment.tail( rowsTotal - rowsAscent ), constraint_BendingMoment, propagationStepSize, timeOfFlight_Descent * constraint_BendingMoment );
    //std::cout << "penaltyBendingMomentDescent: " << penaltyBendingMomentDescent << std::endl;

    /*
    if( debugInfo == 1  ){ std::cout << "Penalty: Flight Path - Descent" << std::endl; }
    const double penaltyFlightPathAngleDescent = bislip::Variables::computeConstraintViolationPenalty( depVar_FlightPathAngle.tail( rowsTotal - rowsAscent ), 0.0, propagationStepSize, 360.0 );

    if( debugInfo == 1  ){ std::cout << "Penalty: Pitch Moment Coefficient - Descent" << std::endl; }
    const double penaltyPitchMomentCoefficientDescent = bislip::Variables::computeConstraintViolationPenalty( depVar_PitchMomentCoefficient.tail( rowsTotal - rowsAscent ).cwiseAbs(), 0.0, propagationStepSize, tof );

    if( debugInfo == 1  ){ std::cout << "Penalty: Heading Error - Descent" << std::endl; }
    double penaltyHeadingErrorDescent =  bislip::Variables::computeDeadbandViolationPenalty(
                depVar_HeadingError.segment( rowsAscent, rowsTotal - rowsAscent ), depVar_CurrentHeadingErrorDeadband.tail( rowsTotal - rowsAscent ), propagationStepSize, tof );
*/


    if( debugInfo == 1  ){ std::cout << "Populating Fitness Vector" << std::endl; }

    //! Assign values to Fitness vector! At the moment these are all 'objective
    //! functions'. To modify this I have change the header file and define how
    //! many elements are objective functions, equality contraints, and
    //! inequality constraints. This vector here must contain them is that exact order: nOF, nEC, nIC.
    std::vector< double > fitness;

    //fitness.push_back( penaltyNoFlight + costCumulativeDistanceTravelled + penaltyGroundtrack +
    //                 penaltyMonotonicApproach + penaltyDistanceToGo + penaltyExcessiveBankReversals );

    fitness.push_back( penaltyDistanceToGo );

    //fitness.push_back( costMaximumNormalizedSpecificEnergyAscent + penaltyMonotonicEnergyStateAscent +
    //                 penaltydelV_Ascent + costFuelMassAscent + penaltyFuelMassAscent );


    fitness.push_back( costHeatLoadChapmanAscent + penaltyHeatFluxChapmanAscent + penaltyDynamicPressureAscent +
                       penaltyMechanicalLoadAscent + penaltyBendingMomentAscent );


    //fitness.push_back( penaltyFlightPathAngleAscent + penaltyFlightPathAngleRateAscent +
    //                 penaltyPitchMomentCoefficientAscent + penaltyHeadingErrorAscent );

    /*
    fitness.push_back( costMaximumNormalizedSpecificEnergyDescent + penaltyMonotonicEnergyStateDescent +
                     penaltydelV_Descent + costFuelMassDescent + penaltyFuelMassDescent );

    if( debugInfo == 1 ){ std::cout << "costMaximumNormalizedSpecificEnergyDescent = " << costMaximumNormalizedSpecificEnergyDescent << std::endl; }
    if( debugInfo == 1 ){ std::cout << "penaltyMonotonicEnergyStateDescent = " << penaltyMonotonicEnergyStateDescent << std::endl; }
    if( debugInfo == 1 ){ std::cout << "penaltydelV_Descent = " << penaltydelV_Descent << std::endl; }
    if( debugInfo == 1 ){ std::cout << "costFuelMassDescent = " << costFuelMassDescent << std::endl; }
    if( debugInfo == 1 ){ std::cout << "penaltyFuelMassDescent = " << penaltyFuelMassDescent << std::endl; }
    if( debugInfo == 1 ){ std::cout << " " << std::endl; }
    if( debugInfo == 1 ){ std::cout << " " << std::endl; }
    if( debugInfo == 1 ){ std::cout << "initialMass_Descent = " << initialMass_Descent << std::endl; }
    if( debugInfo == 1 ){ std::cout << "bodyMap_.at(  vehicleName )->getVehicleSystems()->getDryMass() = " << bodyMap_.at(  vehicleName )->getVehicleSystems()->getDryMass() << std::endl; }
*/


    fitness.push_back( costHeatLoadChapmanDescent + penaltyHeatFluxChapmanDescent + penaltyDynamicPressureDescent +
                       penaltyMechanicalLoadDescent + penaltyBendingMomentDescent );


    //fitness.push_back( penaltyFlightPathAngleDescent + penaltyPitchMomentCoefficientDescent +
    //                 penaltyHeadingErrorDescent );


    /*
    fitness.push_back( penaltyNoFlight  );
    fitness.push_back( costCumulativeDistanceTravelled );
    fitness.push_back( penaltyGroundtrack );
    fitness.push_back( penaltyMonotonicApproach );
    fitness.push_back( penaltyDistanceToGo );
    fitness.push_back( penaltyExcessiveBankReversals );


    fitness.push_back( costMaximumNormalizedSpecificEnergyAscent );
    fitness.push_back( penaltyMonotonicEnergyStateAscent );
    fitness.push_back( penaltydelV_Ascent );
    fitness.push_back( costFuelMassAscent );
    fitness.push_back( penaltyFuelMassAscent );


    fitness.push_back( costHeatLoadChapmanAscent );
    fitness.push_back( penaltyHeatFluxChapmanAscent );
    fitness.push_back( penaltyDynamicPressureAscent );
    fitness.push_back( penaltyMechanicalLoadAscent );
    fitness.push_back( penaltyBendingMomentAscent );


    fitness.push_back( penaltyFlightPathAngleAscent );
    fitness.push_back( penaltyFlightPathAngleRateAscent );
    fitness.push_back( penaltyPitchMomentCoefficientAscent );
    fitness.push_back( penaltyHeadingErrorAscent );


    fitness.push_back( costMaximumNormalizedSpecificEnergyDescent );
    fitness.push_back( penaltyMonotonicEnergyStateDescent );
    fitness.push_back( penaltydelV_Descent );
    fitness.push_back( costFuelMassDescent );
    fitness.push_back( penaltyFuelMassDescent );


    fitness.push_back( costHeatLoadChapmanDescent );
    fitness.push_back( penaltyHeatFluxChapmanDescent );
    fitness.push_back( penaltyDynamicPressureDescent );
    fitness.push_back( penaltyMechanicalLoadDescent );
    fitness.push_back( penaltyBendingMomentDescent );


    fitness.push_back( penaltyFlightPathAngleDescent );
    fitness.push_back( penaltyPitchMomentCoefficientDescent );
    fitness.push_back( penaltyHeadingErrorDescent );

*/

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        EVALUATE CRITERIA TO PRINT DATA               //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd hardConstraintViolation( 10 );
    hardConstraintViolation = Eigen::VectorXd::Zero( hardConstraintViolation.size() );

    if( bislipSystems->getValidationFlag() != true )
    {
        if( depVar_TimeOfFlight( rowsTotal - 1 ) < problemInput_->getHardConstraints()[ 0 ] ) { hardConstraintViolation( 0 ) = 1.0; }
        hardConstraintViolation( 1 ) = bislip::Variables::determineNumberOfHardViolations( depVar_HeadingError, problemInput_->getHardConstraints()[ 1 ] );
        if( depVar_AngularDistanceToGo( rowsTotal - 1 ) > problemInput_->getHardConstraints()[ 2 ] ) { hardConstraintViolation( 2 ) = 1.0; }
        hardConstraintViolation( 3 ) = bislip::Variables::determineNumberOfHardViolations( depVar_DynamicPressure, problemInput_->getHardConstraints()[ 3 ] );
        hardConstraintViolation( 4 ) = bislip::Variables::determineNumberOfHardViolations( depVar_BendingMoment, problemInput_->getHardConstraints()[ 4 ] );
        hardConstraintViolation( 5 ) = bislip::Variables::determineNumberOfHardViolations( depVar_HeatFluxChapman, problemInput_->getHardConstraints()[ 5 ] );
        hardConstraintViolation( 6 ) = bislip::Variables::determineNumberOfHardViolations( depVar_BodyFixedMechanicalLoadMag, problemInput_->getHardConstraints()[ 6 ] );
        if( depVar_Height( rowsTotal - 1 ) > problemInput_->getHardConstraints()[ 7 ] ) { hardConstraintViolation( 7 ) = 1.0; }
        if( depVar_NormalizedSpecificEnergy( rowsTotal - 1 ) > problemInput_->getHardConstraints()[ 8 ] ) { hardConstraintViolation( 8 ) = 1.0; }
        hardConstraintViolation( 9 ) = 0;//bislip::Variables::determineNumberOfHardViolations( depVar_NormalizedSpecificEnergy, 1.0 );
    }

    bool saveTrajectoryOutput = false;
    int sumOfViolations = hardConstraintViolation.sum();
    if( sumOfViolations == 0 ) { saveTrajectoryOutput = true; }


    for( unsigned int i = 0; i < hardConstraintViolation.size() - 1; i++)
    {
        std::cout  << hardConstraintViolation[ i ] << " | " ;
    }
    std::cout  << hardConstraintViolation[ hardConstraintViolation.size() - 1 ] <<  " || " <<  sumOfViolations << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PRINT SIMULATION OUTPUT TO FILE               //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Get time stamp for this specific simulation. This avoids overwriting the
    //! file if another individual with the same properties shows up in other
    //! evolutions.

    std::chrono::time_point< std::chrono::system_clock > printTime = bislip::Variables::getDateTime( );

    std::string printTimeString = bislip::Variables::convertDateTimeToString( false, printTime );

    unsigned int millisSincePlayTime = std::chrono::duration_cast< std::chrono::milliseconds >( printTime - ( problemInput_->getPlayTimePair() ).first ).count();

    //! Create unique filename that cannot be overwritten due to the timestamp.
    std::string simulation_file_name_suffix =
            ( problemInput_->getPlayTimePair() ).second + "_" +
            printTimeString + "_" +
            std::to_string( millisSincePlayTime ) + "_" +
            std::to_string( timeOfFlight ) + "_" +
            std::to_string( finalHeight ) + "_" +
            std::to_string( targetLat_deg_calc ) + "_" +
            std::to_string( targetLon_deg_calc ) + "_" +
            std::to_string( finalAngularDistanceToGo_deg ) + "_" +
            std::to_string( maximumMechanicalLoad ) + "_" +
            std::to_string( finalMechanicalLoad );

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


    if( debugInfo == 1 ){ std::cout << "Printing to files" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Output Settings Values 1: " << problemInput_->getOutputSettings()[ 1 ] << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Output Settings Values 2: " << problemInput_->getOutputSettings()[ 2 ] << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Output Settings Values 3: " << problemInput_->getOutputSettings()[ 3 ] << std::endl; }


    //! Convert decision vector from std::vector< double > to Eigen::VectorXd
    const Eigen::VectorXd x_Eigen = Eigen::Map< const Eigen::VectorXd, Eigen::Unaligned >( x.data(), x.size() );

    //! Create new Eigen::VectorXd to include flag that the individual was printed or not
    Eigen::VectorXd x_EigenFlagged( x_Eigen.size() + 1 );
    x_EigenFlagged << 0, x_Eigen;

    std::map < std::string, Eigen::VectorXd > populationMap = problemInput_->getPopulation();
    populationMap[ simulation_file_name_suffix ] = x_EigenFlagged ;
    problemInput_->setPopulation( populationMap );

    //! Convert fitness vector from std::vector< double > to Eigen::VectorXd
    const Eigen::VectorXd fitness_Eigen = Eigen::Map< const Eigen::VectorXd, Eigen::Unaligned >( fitness.data(), fitness.size() );

    //! Create new Eigen::VectorXd to include flag that the individual was printed or not
    Eigen::VectorXd fitness_EigenFlagged( fitness_Eigen.size() + 1 );
    fitness_EigenFlagged << 0, fitness_Eigen;

    std::map < std::string, Eigen::VectorXd > fitnessMap    = problemInput_->getFitness();
    fitnessMap[ simulation_file_name_suffix ] = fitness_EigenFlagged ;
    problemInput_->setFitness( fitnessMap );

    //! Will print out depending on some input values. Each entry corresponds to
    //! a different type of output. Entries are binary, 0 or 1.
    if( ( int( problemInput_->getOutputSettings()[ 1 ] ) == 1 && saveTrajectoryOutput == true ) )
    {

        x_EigenFlagged( 0 ) = 1;
        fitness_EigenFlagged( 0 ) = 1;

        populationMap[ simulation_file_name_suffix ] = x_EigenFlagged ;
        problemInput_->setPopulation( populationMap );

        fitnessMap[ simulation_file_name_suffix ] = fitness_EigenFlagged ;
        problemInput_->setFitness( fitnessMap );

        if( debugInfo == 1 ){ std::cout << "Saving Propagation" << std::endl; }
        if( debugInfo == 1 ){ std::cout << "propagationTimeHistoryMap.size() = " << propagationTimeHistoryMap.size() << std::endl; }

        //! Write propagation history to file.
        tudat::input_output::writeDataMapToTextFile(
                    propagationTimeHistoryMap,
                    complete_file_name_Prop,
                    problemInput_->getOutputPath() + problemInput_->getOutputSubFolder() + "/" + "propagationHistory" + "/",
                    "",
                    std::numeric_limits< double >::digits10,
                    std::numeric_limits< double >::digits10,
                    "," );
        if( debugInfo == 1 ){ std::cout << "Propagation Saved" << std::endl; }
    }
    if( ( int( problemInput_->getOutputSettings()[ 2 ] ) == 1  && saveTrajectoryOutput == true ) )
    {
        if( debugInfo == 1 ){ std::cout << "Saving Dependent Variables" << std::endl; }
        if( debugInfo == 1 ){ std::cout << "depVarTimeHistoryMap.size() = " << depVarTimeHistoryMap.size() << std::endl; }


        //std::cout << "( propagationTimeHistoryMap.begin() )->first: " << ( propagationTimeHistoryMap.begin() )->first << std::endl;
        //std::cout << "( propagationTimeHistoryMap.rbegin() )->first: " << ( propagationTimeHistoryMap.rbegin() )->first << std::endl;
        //std::cout << "( propagationTimeHistoryMap.rbegin() )->second: " << ( propagationTimeHistoryMap.rbegin() )->second << std::endl;


        //std::cout << "( depVarTimeHistoryMap.begin() )->first: " << ( depVarTimeHistoryMap.begin() )->first << std::endl;
        //std::cout << "( depVarTimeHistoryMap.rbegin() )->first: " << ( depVarTimeHistoryMap.rbegin() )->first << std::endl;
        //std::cout << "( depVarTimeHistoryMap.rbegin() )->second: " << ( depVarTimeHistoryMap.rbegin() )->second << std::endl;

        //for ( unsigned int ii = 0; ii < timeVector.size(); ii++)
        //{
        //std::cout << "timeVector[ " << timeVector.size() - 1 << " ] =  " << timeVector[ timeVector.size() - 1 ]  << std::endl;
        //}


        /*  tudat::input_output::writeDataMapToTextFile(
                    depVarTimeHistoryMap,
                    complete_file_name_DepVar,
                    outputPath_ + outputSubFolder_,
                    "",
                    std::numeric_limits< double >::digits10,
                    std::numeric_limits< double >::digits10,
                    "," );
                    */

        Eigen::MatrixXd depVarTimeHistoryMatrixOutput( rowsTotal, columns + 1 );

        Eigen::MatrixXd timeVectorXd( rowsTotal, 1 );

        for ( int i = 0; i < rowsTotal; i++ )
        { timeVectorXd( i ) = timeVector[ i ]; }

        depVarTimeHistoryMatrixOutput << timeVectorXd, depVarTimeHistoryMatrix;

        tudat::input_output::writeMatrixToFile(
                    depVarTimeHistoryMatrixOutput,
                    complete_file_name_DepVar,
                    std::numeric_limits< double >::digits10,
                    problemInput_->getOutputPath() + problemInput_->getOutputSubFolder() + "/" + "dependentVariables" + "/",
                    ",",
                    "");


        if( debugInfo == 1 ){ std::cout << "Dependent Variables Saved" << std::endl; }
    }
    if( ( int( problemInput_->getOutputSettings()[ 3 ] ) == 1  && saveTrajectoryOutput == true ) )
    {
        if( debugInfo == 1 ){ std::cout << "Saving Evaluated Interpolators - Ascent" << std::endl; }

        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsAscent,
                                                     complete_file_name_interpolators_Ascent,
                                                     problemInput_->getOutputPath() + problemInput_->getOutputSubFolder() + "/" + "evaluatedInterpolatorsAscent" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );

        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsAscent_LB,
                                                     complete_file_name_interpolators_Ascent_LB,
                                                     problemInput_->getOutputPath() + problemInput_->getOutputSubFolder() + "/" + "evaluatedInterpolatorsAscent_LB" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsAscent_UB,
                                                     complete_file_name_interpolators_Ascent_UB,
                                                     problemInput_->getOutputPath() + problemInput_->getOutputSubFolder() + "/" + "evaluatedInterpolatorsAscent_UB" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );

        if( debugInfo == 1 ){ std::cout << "Evaluated Interpolators - Ascent Saved" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "Saving Evaluated Interpolators - Descent" << std::endl; }
        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsDescent,
                                                     complete_file_name_interpolators_Descent,
                                                     problemInput_->getOutputPath() + problemInput_->getOutputSubFolder() + "/" + "evaluatedInterpolatorsDescent" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsDescent_LB,
                                                     complete_file_name_interpolators_Descent_LB,
                                                     problemInput_->getOutputPath() + problemInput_->getOutputSubFolder() + "/" + "evaluatedInterpolatorsDescent_LB" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        tudat::input_output::writeDataMapToTextFile( evaluatedInterpolatorsDescent_UB,
                                                     complete_file_name_interpolators_Descent_UB,
                                                     problemInput_->getOutputPath() + problemInput_->getOutputSubFolder() + "/" + "evaluatedInterpolatorsDescent_UB" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );


        if( debugInfo == 1 ){ std::cout << "Evaluated Interpolators - Descent Saved" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "Saving Decision Vector - Ascent" << std::endl; }
        tudat::input_output::writeDataMapToTextFile( map_DV_mapped_Ascent,
                                                     complete_file_name_map_DV_mapped_Ascent,
                                                     problemInput_->getOutputPath() + problemInput_->getOutputSubFolder() + "/" + "map_DV_mapped_Ascent" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        if( debugInfo == 1 ){ std::cout << "Decision Vector - Ascent Saved" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "Saving Decision Vector - Descent" << std::endl; }
        tudat::input_output::writeDataMapToTextFile( map_DV_mapped_Descent,
                                                     complete_file_name_map_DV_mapped_Descent,
                                                     problemInput_->getOutputPath() + problemInput_->getOutputSubFolder() + "/" + "map_DV_mapped_Descent" + "/",
                                                     "",
                                                     std::numeric_limits< double >::digits10,
                                                     std::numeric_limits< double >::digits10,
                                                     "," );
        if( debugInfo == 1 ){ std::cout << "Decision Vector - Descent Saved" << std::endl; }
    }


    if( debugInfo == 1 ){ std::cout << "Printing to terminal" << std::endl; }
    for( unsigned int i = 0; i < fitness.size() - 1; i++)
    {
        std::cout  << fitness[ i ] << " | " ;
    }
    //std::cout  << fitness[ fitness.size() - 1 ] << "  ||  " << rowsAscent << "  ||  " << tof << "  ||  " << " Theoretical delV = " << theoreticaldelV_Ascent << "  ||  " << " Goal delV = " << goaldelV_Ascent << "  ||  " << " Actual delV = " << actualdelV_Ascent << std::endl;
    //    std::cout  << fitness[ fitness.size() - 1 ] << "  ||  " << rowsAscent << "  ||  " << tof <<  std::endl;
    std::cout  << fitness[ fitness.size() - 1 ] << " ||  " << timeOfFlight <<  std::endl;


    //  }



    //std::cout << "objectiveHeight_Ascent: " << objectiveHeight_Ascent << std::endl;
    //std::cout << "objectiveHeight_Ascent_calc: " << objectiveHeight_Ascent_calc << std::endl;



    // std::cout <<dynamicsSimulator.getPropagationTerminationReason( )->getPropagationTerminationReason( ) << std::endl;
    //std::this_thread::sleep_for( std::chrono::nanoseconds( 10 ) );
    //using namespace std::this_thread; // sleep_for, sleep_until
    //using namespace std::chrono; // nanoseconds, system_clock, seconds




    return fitness;

} // Fitness function.
