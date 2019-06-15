#include "updateGuidance.h"

namespace bislip {

void MyGuidance::updateGuidance( const double currentTime )
{

    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap_.at( vehicleName_ )->getVehicleSystems();
    double startingEpoch = bislipSystems->getStartingEpoch();
    double guidanceStepSize = bislipSystems->getGuidanceStepSize();

    //if( ( bislip::Variables::millis_since_midnight() - bislipSystems->getStartingMillis() ) > 120000 ) { bislipSystems->setDebugInfo( 1 ); }

    int debugInfo = bislipSystems->getDebugInfo();
    if( debugInfo == 1 ){ std::cout << "  " <<  std::endl; }
    if( debugInfo == 1 ){ std::cout << "Starting guidance for this evaluation" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "---------------------------------------------------------------------" <<  std::endl; }


    if( debugInfo == 4 ){std::cout << "Millis since start = " << ( bislip::Variables::millis_since_midnight() - bislipSystems->getStartingMillis() ) << std::endl; }

    //debugInfo = bislipSystems->getDebugInfo();

    std::string currentTrajectoryPhase = bislipSystems->getCurrentTrajectoryPhase();

    double timeOfFLight = 1.0;
    double updateConditionalPre = 1.0;
    double updateConditional = 1.0;
    double resolution = 100;

    //if ( std::isnan( currentTime ) == false )
    //{
    if( bislipSystems->getValidationFlag() == true ) { resolution = 100; }

    //timeOfFLight = ( std::trunc( ( currentTime - startingEpoch ) * resolution ) ) / resolution;
    timeOfFLight = ( std::round( ( currentTime - startingEpoch ) * resolution)  ) / resolution;
    updateConditionalPre = std::abs( std::remainder( timeOfFLight, guidanceStepSize ) );
    updateConditional = updateConditionalPre;
    if( debugInfo == 1 ){ std::cout << "  " <<  std::endl; }
    if( debugInfo == 1 ){ std::cout << "startingEpoch = " << startingEpoch << "  |  " << "currentTime = " << currentTime << "  |  " << "timeOfFLight = " <<  timeOfFLight << "  |  std::remainder( " << timeOfFLight << ", " << guidanceStepSize << " ) = " <<  updateConditionalPre << std::endl; }
    // if( debugInfo == 1 ){std::cout << "std::trunc( updateConditionalPre * 100 ) / 100 = " << updateConditional << std::endl; }

    //}


    if ( currentTime == currentTime )
    {
        if( FlightConditions_ == nullptr )
        {
            FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                        bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
        }
        if( debugInfo == 1 ){ std::cout << "Flight Conditions are Available  " <<  std::endl; }

        if( coefficientInterface_ == nullptr )
        {
            coefficientInterface_ = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                        bodyMap_.at( vehicleName_ )->getAerodynamicCoefficientInterface( ) );
        }
        if( debugInfo == 1 ){ std::cout << "Coefficient Interface is Available  " <<  std::endl; }
    }

    if ( std::isnan( currentTime ) == false )
    {
        if( debugInfo == 1 ){std::cout << "Calculate and Set Flight-Path Angle Rate" << std::endl; }
        // bislipSystems->setCurrentFlightPathAngleRate( bislip::Variables::computeFlightPathAngleRate( bodyMap_, vehicleName_, centralBodyName_ ) );

        if( debugInfo == 1 ){ std::cout << "Calculate and Set Current Local Gravity Vector" << std::endl; }
        bislipSystems->setCurrentLocalGravityVector( bislip::Variables::computeLocalGravity( bodyMap_, vehicleName_, centralBodyName_ ) );
    }


    if( debugInfo == 1 )
    {
        std::cout << "Choosing which guidance conditional to run" << std::endl;
        std::cout << "updateConditional: " << updateConditional << "  |  " << "Time of flight = " << timeOfFLight << std::endl;
    }

    if( updateConditional < 1e-3 )
    {
        if( debugInfo == 1 ){std::cout << "Evaluating Guidance Functions" << std::endl; }
        if( debugInfo == 1 ){std::cout << " " << std::endl; }
        if( debugInfo == 1 ){std::cout << " " << std::endl; }
        if( debugInfo == 1 ){std::cout << " " << std::endl; }

        if( bislipSystems->getValidationFlag() != true ) { evaluateGuidanceFunctions( bislipSystems, vehicleSystems, currentTrajectoryPhase, timeOfFLight ); }
        else { evaluateValidationGuidanceFunctions( bislipSystems, vehicleSystems, currentTrajectoryPhase, timeOfFLight ); }

        if( debugInfo == 1 ){ std::cout << "          Saving bank reversal trigger if applicable."<< std::endl; }
        if( debugInfo == 1 ){ std::cout << "            FlightConditions_->getCurrentTime() = "<< FlightConditions_->getCurrentTime() << std::endl; }
        if( debugInfo == 1 ){ std::cout << "            bislipSystems->getBankAngleReversalTimepoint() = "<< bislipSystems->getBankAngleReversalTimepoint() << std::endl; }

        if ( FlightConditions_->getCurrentTime() != bislipSystems->getBankAngleReversalTimepoint() )
        {
            if ( std::isnan( FlightConditions_->getCurrentTime() ) == false  )
            {
                if( debugInfo == 1 ){ std::cout << "               bislipSystems->getCurrentBankAngle() * bislipSystems->getTempBankAngle( ) = "<< bislipSystems->getCurrentBankAngle() * bislipSystems->getTempBankAngle( ) <<std::endl; }

                if ( bislipSystems->getCurrentBankAngle() * bislipSystems->getTempBankAngle( ) < 0.0 )
                {
                    bislipSystems->setBankAngleReversalTrigger( true );
                    if( debugInfo == 1 ){ std::cout << "               bislipSystems->getBankAngleReversalTrigger = "<< bislipSystems->getBankAngleReversalTrigger() <<std::endl; }

                    bislipSystems->setBankAngleReversalTimepoint( FlightConditions_->getCurrentTime() );
                    if( debugInfo == 1 ){ std::cout << "               NEW bislipSystems->getBankAngleReversalTimepoint = "<< bislipSystems->getBankAngleReversalTimepoint() <<std::endl; }
                }
                else { bislipSystems->setBankAngleReversalTrigger( false ); }
            }
        }


    }
    else
    {
        if( debugInfo == 1 ){std::cout << "No Guidance Conditional" << std::endl; }
        if( debugInfo == 1 ){std::cout << " " << std::endl; }
        if( debugInfo == 1 ){std::cout << " " << std::endl; }
        if( debugInfo == 1 ){std::cout << " " << std::endl; }

        currentAngleOfAttack_ = bislipSystems->getCurrentAngleOfAttack();

        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", bislipSystems->getCurrentBodyFlapAngle() );
        vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonLeft", bislipSystems->getCurrentElevonAngle() );
        vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonRight", bislipSystems->getCurrentElevonAngle() );

        currentBankAngle_ = bislipSystems->getCurrentBankAngle();

        currentAngleOfSideslip_ = 0.0;

    }
    // }

    if ( std::isnan( currentTime ) == false )
        //    if ( currentTime == currentTime )
    {

        //if( debugInfo == 1 ){ std::cout << "Calculate and Set Current Lift Force" << std::endl; }
        // bislipSystems->setCurrentLiftForce( bislip::Variables::computeCurrentLiftForce( bodyMap_, vehicleName_ ) );

        // if( debugInfo == 1 ){ std::cout << "Calculate and Set Current Drag Force" << std::endl; }
        // bislipSystems->setCurrentDragForce( bislip::Variables::computeCurrentDragForce( bodyMap_, vehicleName_ ) );



        if( debugInfo == 1 ){std::cout << "Calculating Heat Fluxes" << std::endl; }
        if( debugInfo == 1 ){std::cout << " " << std::endl; }
        if( debugInfo == 1 ){std::cout << " " << std::endl; }
        if( debugInfo == 1 ){std::cout << " " << std::endl; }

        if( debugInfo == 1 ){ std::cout << "     currentAirSpeed       = " << FlightConditions_->getCurrentAirspeed() << std::endl; }

        //************************************************** Heat Flux(es)
        if( debugInfo == 1 ){ std::cout << "Determine and set heat fluxes." << std::endl; }
        if( debugInfo == 1 ){ std::cout << "     Chapman Stagnation Heat Flux" << std::endl; }
        bislipSystems->setWorkingRadius( bislipSystems->getNoseRadius() );
        bislipSystems->setWallTemperature( bislipSystems->getChapmanWallTemp() );

        if( debugInfo == 1 ){ std::cout << "        Working Radius = " << bislipSystems->getWorkingRadius() << std::endl; }
        if( debugInfo == 1 ){ std::cout << "        Initial Wall Temperature = " <<  bislipSystems->getWallTemperature( ) << std::endl; }
        bislipSystems->setCurrentHeatFluxChapman( bislip::Variables::computeHeatingRateChapman( bodyMap_, vehicleName_) );
        bislipSystems->setChapmanWallTemp( bislipSystems->getWallTemperature() );

        if( debugInfo == 1 ){ std::cout << "        Current Chapman Stagnation Heat Flux = " <<  bislipSystems->getCurrentHeatFluxChapman( ) << std::endl; }
        if( debugInfo == 1 ){ std::cout << "        Current Chapman Stagnation Eq. Temp. = " <<  bislipSystems->getChapmanWallTemp( ) << std::endl; }

        if( debugInfo == 1 ){ std::cout << "     Tauber Heat Flux" << std::endl; }
        bislipSystems->setWorkingRadius( bislipSystems->getLeadingEdgeRadius() );

        if( debugInfo == 1 ){ std::cout << "        Working Radius = " << bislipSystems->getWorkingRadius() << std::endl; }
        bislipSystems->setCurrentHeatFluxTauber( bislip::Variables::computeHeatingRateTauber( bodyMap_, vehicleName_) );

        if( debugInfo == 1 ){ std::cout << "        Current Tauber Heat Flux = " <<  bislipSystems->getCurrentHeatFluxTauber( ) << std::endl; }
    }



    if( debugInfo == 1 ){ std::cout << "  " <<  std::endl; }
    if( debugInfo == 1 ){ std::cout << "  Various Relevant Values" <<  std::endl; }

    if( std::isnan( currentTime ) == false )
        //     if ( currentTime == currentTime )
    {
        if( debugInfo == 1 ){ std::cout << "  " <<  std::endl; }
        if( debugInfo == 1 ){ std::cout << "     Current Time     = " << currentTime - startingEpoch <<  std::endl; }
        if( debugInfo == 1 ){ std::cout << "     Current Height   = " << FlightConditions_->getCurrentAltitude() <<  std::endl; }
        if( debugInfo == 1 ){ std::cout << "     Current Airspeed = " << FlightConditions_->getCurrentAirspeed() <<  std::endl; }
        if( debugInfo == 1 ){ std::cout << "     Current density  = " << FlightConditions_->getCurrentDensity() <<  std::endl; }
    }


    if( debugInfo == 1 ){ std::cout << "  " <<  std::endl; }
    if( debugInfo == 1 ){ std::cout << "  " <<  std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Current Angle of Attack  = " << tudat::unit_conversions::convertRadiansToDegrees( currentAngleOfAttack_ ) <<  std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Current Bank Angle       = " << tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) <<  std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Current Sideslip Angle   = " << tudat::unit_conversions::convertRadiansToDegrees( currentAngleOfSideslip_ ) <<  std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Current Bodyflap angle   = " << tudat::unit_conversions::convertRadiansToDegrees( vehicleSystems->getCurrentControlSurfaceDeflection( "BodyFlap" ) ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Current Elevon angle     = " << tudat::unit_conversions::convertRadiansToDegrees( vehicleSystems->getCurrentControlSurfaceDeflection( "ElevonLeft" ) ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "     Current Thrust Magnitude = " << bislipSystems->getCurrentThrustMagnitude() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "  " <<  std::endl; }
    if( debugInfo == 1 ){ std::cout << "  " <<  std::endl; }


    if ( std::isnan( currentAngleOfAttack_ ) == true ) { throw std::runtime_error("Angle of attack not set."); }
    if ( std::isnan( currentBankAngle_ ) == true ) { throw std::runtime_error("Bank Angle not set."); }
    if ( std::isnan( currentAngleOfSideslip_ ) == true ) { throw std::runtime_error("Sideslip Angle not set."); }

    if( debugInfo == 1 ){ std::cout << "  " <<  std::endl; }
    if( debugInfo == 1 ){ std::cout << "Guidance for this evaluation is done" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "---------------------------------------------------------------------" <<  std::endl; }



    /*
    if( ( std::isnan( currentTime ) == false && bislipSystems->getCurrentThrottleSetting() == 0.0 &&  bislipSystems->getCurrentEngineStatus() == true ) )
    {
        std::cout << "  " << std::endl;
        std::cout << "  " << std::endl;
        std::cout << "HOLD UP!!!!!!!" << std::endl;
        std::cout << "  " << std::endl;
         std::cout << "  " << std::endl;
    }

    */

} // updateGuidance


void MyGuidance::evaluateGuidanceFunctions(
        // std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems,
        std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
        const std::string &currentTrajectoryPhase,
        const double &timeOfFLight )
{
    int debugInfo = bislipSystems->getDebugInfo();

    double currentMachNumber      = 0.0;
    double currentHeight          = 0.0;
    double currentAirspeed        = 0.0;
    double currentDynamicPressure = 0.0;
    double currentFlightPathAngle = 0.0;
    double currentHeadingAngle    = 0.0;
    double currentLatitude        = 0.0;
    double currentLongitude       = 0.0;
    double previousBankAngle      = 0.0;
    double timeOfFlight           = 0.0;

    if( debugInfo == 1 ){ std::cout << "      Selecting Source of Values" << std::endl; }

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentMachNumber      = bislipSystems->getInitialMachNumber();
        currentHeight          = bislipSystems->getInitialHeight();
        currentAirspeed        = bislipSystems->getInitialAirspeed();
        currentDynamicPressure = bislipSystems->getInitialDynamicPressure();
        currentFlightPathAngle = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialFlightPathAngle() );
        currentHeadingAngle    = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialHeading() );
        currentLatitude        = bislipSystems->getInitialLat();
        currentLongitude       = bislipSystems->getInitialLon();
        previousBankAngle      = 0.0;
        timeOfFlight           = 0;

        //bislipSystems->getCurrentMachNumber()
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentMachNumber      = FlightConditions_->getCurrentMachNumber();
        currentHeight          = FlightConditions_->getCurrentAltitude();
        currentAirspeed        = FlightConditions_->getCurrentAirspeed();
        currentDynamicPressure = FlightConditions_->getCurrentDynamicPressure();
        currentFlightPathAngle = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        currentHeadingAngle    = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
        currentLatitude        = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        currentLongitude       = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );
        previousBankAngle      = bislipSystems->getCurrentBankAngle();
        timeOfFlight           = FlightConditions_->getCurrentTime() - bislipSystems->getStartingEpoch();
    }

    if( debugInfo == 1 ){ std::cout << "    Determine and Set Engine Status" << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    -----------------------" << std::endl; }
    bislipSystems->setCurrentEngineStatus( bislip::Variables::determineEngineStatus( bodyMap_.at( vehicleName_ )->getBodyMass(), vehicleSystems->getDryMass() ) );
    if( debugInfo == 1 ){ std::cout << "    -----------------------" << std::endl; }


    if( debugInfo == 1 ){ std::cout << " "<< std::endl; }
    if( debugInfo == 1 ){ std::cout << "Validation case?          = " << bislipSystems->getValidationFlag() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Trajectory Phase          = " << bislipSystems->getCurrentTrajectoryPhase() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Time of Flight            = " << timeOfFlight << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Height            = " << currentHeight << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Airspeed          = " << currentAirspeed << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Dynamic Pressure  = " << currentDynamicPressure << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Mach Number       = " << currentMachNumber << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Flight-path Angle = " << tudat::unit_conversions::convertRadiansToDegrees( currentFlightPathAngle ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Mass              = " << bodyMap_.at( vehicleName_ )->getBodyMass() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Engine Status     = " << bislipSystems->getCurrentEngineStatus() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Starting Angle of Attack  = " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentAngleOfAttack() ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Starting Bank Angle       = " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentBankAngle() ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << " "<< std::endl; }


    //const double normalizedSpecificEnergy = bislip::Variables::computeNormalizedSpecificEnergy( flightConditions->getCurrentAltitude(), flightConditions->getCurrentAirspeed(), bislipSystems->getE_max() );
    // double resolution = 100;
    // if( bislipSystems->getValidationFlag() == true ) { resolution = 10; }

    //const double timeOfFLight = ( std::trunc( ( FlightConditions_->getCurrentTime() - bislipSystems->getStartingEpoch() ) * resolution ) ) / resolution;
    //const double timeOfFLight = std::trunc( ( FlightConditions_->getCurrentTime() - bislipSystems->getStartingEpoch() ) * 100 ) / 100;

    //! Determine and set values by evaluating interpolators.
    if( debugInfo == 1 ){ std::cout << "Determine NEW Values by Evaluating Interpolators." << std::endl; }

    //************************************************** Angle of Attack

    double evaluatedAngleOfAttack = 0.0;
    // if( normalizedSpecificEnergy < 1.0 )
    // {
    evaluatedAngleOfAttack = tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::AngleOfAttack, bodyMap_, vehicleName_, centralBodyName_ ) );
    //}
    bislipSystems->setCurrentAngleOfAttack( evaluatedAngleOfAttack );

    //double optimizedAngleofAttack = 0.0;
    //double bodyFixedTotal_g_Load_Magnitude = bislip::Variables::computeBodyFixedTotal_g_Load_Magnitude( bodyMap_, vehicleName_ );
    /*if( bislipSystems->getCurrentEngineStatus( ) == true )
        {
            const double throttleDifference1 = std::abs( bislipSystems->getCurrentThrottleSetting() - ( bislipSystems->getThrottleSettingLimits() ).first );
            const double throttleDifference2 = std::abs( bislipSystems->getCurrentThrottleSetting() - ( bislipSystems->getThrottleSettingLimits() ).second );

            if( ( throttleDifference1 < -1e-7 || throttleDifference2 < -1e-7 ) )
            {
                //if( bodyFixedTotal_g_Load_Magnitude < bislipSystems->getMechanicalLoadConstraint() )
                //{
                if( debugInfo == 1 ){ std::cout << "Verify conditions for Optimization of Angle of Attack" << std::endl; }
                if( debugInfo == 1 ){ std::cout << "   Engine is ON" << std::endl; }
                if( debugInfo == 1 ){ std::cout << "   Throttle Setting        = " << bislipSystems->getCurrentThrottleSetting() << std::endl; }
                if( debugInfo == 1 ){ std::cout << "   Throttle Difference 1   = " << throttleDifference1 << std::endl; }
                if( debugInfo == 1 ){ std::cout << "   Throttle Difference 1   = " << throttleDifference2 << std::endl; }
                if( debugInfo == 1 ){ std::cout << "   Current Mechanical Load = " << bodyFixedTotal_g_Load_Magnitude << std::endl; }

                optimizedAngleofAttack = bislip::Variables::optimizeAngleOfAttack( bodyMap_, vehicleName_ );
                bislipSystems->setCurrentAngleOfAttack( optimizedAngleofAttack );
                //}
            }

        }
        else if( bislipSystems->getCurrentEngineStatus( ) != true )
        {
            if( debugInfo == 1 ){ std::cout << "Verify conditions for Optimization of Angle of Attack" << std::endl; }
            if( debugInfo == 1 ){ std::cout << "   Engine is OFF" << std::endl; }
            if( debugInfo == 1 ){ std::cout << "   Current Mechanical Load = " << bodyFixedTotal_g_Load_Magnitude << std::endl; }

            optimizedAngleofAttack = bislip::Variables::optimizeAngleOfAttack( bodyMap_, vehicleName_ );
            bislipSystems->setCurrentAngleOfAttack( optimizedAngleofAttack );
        }
             */

    if( debugInfo == 1 ){ std::cout << "Extract bounds of angle of attack as a function of optimization parameter limits" << std::endl; }
    double normalizedSpecificEnergy = bislip::Variables::computeNormalizedSpecificEnergy( currentHeight, currentAirspeed, bislipSystems->getMaximumSpecificEnergy() );

    double angleOfAttackLowerBound_1 = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getParameterLowerBounds( bislip::Parameters::Interpolators::AngleOfAttack ) )->interpolate( normalizedSpecificEnergy ) );
    double angleOfAttackUpperBound_1 = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getParameterUpperBounds( bislip::Parameters::Interpolators::AngleOfAttack ) )->interpolate( normalizedSpecificEnergy ) );
    if( debugInfo == 1 ){ std::cout << "    Determined Angle of attack (interpolated)               = " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentAngleOfAttack() ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    Angle of attack Lower Bound due to parameter limits     = " << tudat::unit_conversions::convertRadiansToDegrees( angleOfAttackLowerBound_1 ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    Angle of attack Upper Bound due to parameter limits     = " << tudat::unit_conversions::convertRadiansToDegrees( angleOfAttackUpperBound_1 ) << std::endl; }

    if( debugInfo == 1 ){ std::cout << "Extract bounds of angle of attack as a function of Mach Number" << std::endl; }
    double angleOfAttackLowerBound_2 = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getAlphaMachEnvelopeLBInterpolator() )->interpolate( currentMachNumber ) );
    double angleOfAttackUpperBound_2 = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getAlphaMachEnvelopeUBInterpolator() )->interpolate( currentMachNumber ) );
    if( debugInfo == 1 ){ std::cout << "    Angle of attack Lower Bound due to Mach Number = " << tudat::unit_conversions::convertRadiansToDegrees( angleOfAttackLowerBound_2 ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    Angle of attack Upper Bound due to Mach Number = " << tudat::unit_conversions::convertRadiansToDegrees( angleOfAttackUpperBound_2 ) << std::endl; }

    double angleOfAttackLowerBound = std::max( angleOfAttackLowerBound_1, angleOfAttackLowerBound_2 );
    double angleOfAttackUpperBound = std::min( angleOfAttackUpperBound_1, angleOfAttackUpperBound_2 );

    if ( bislipSystems->getCurrentAngleOfAttack() >= angleOfAttackUpperBound ) { bislipSystems->setCurrentAngleOfAttack( angleOfAttackUpperBound ); }
    else if ( bislipSystems->getCurrentAngleOfAttack() <= angleOfAttackLowerBound ) { bislipSystems->setCurrentAngleOfAttack( angleOfAttackLowerBound ); }

    currentAngleOfAttack_ = bislipSystems->getCurrentAngleOfAttack();

    if( debugInfo == 1 ){ std::cout << "        Bounded angle of attack = " << bislipSystems->getCurrentAngleOfAttack() << std::endl; }


    //************************************************** Control Surface Deflections
    //! Compute and assign the current control surface deflection angles

    Eigen::Vector2d controlSurfaceDeflections;

    if( debugInfo == 1 ){ std::cout << "Compute and assign the current control surface deflections." << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Determine Control Surface Deflections for Validation Case" << std::endl; }

    //if( debugInfo == 1 ){ std::cout << "          Current Dynamic Pressure = " << currentDynamicPressure << std::endl; }
    //if( debugInfo == 1 ){ std::cout << "          Minimum Dynamic Pressure = " << bislipSystems->getMinimumDynamicPressureforControlSurface() << std::endl; }

    //if( currentDynamicPressure > bislipSystems->getMinimumDynamicPressureforControlSurface() )
    //{

        if( bislipSystems->getInitialValueFlag() == true )
        {
            if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values for Surface Deflections" << std::endl; }

            controlSurfaceDeflections << bislipSystems->getCurrentBodyFlapAngle(), bislipSystems->getCurrentElevonAngle(); }
        else
        {
            if( debugInfo == 1 ){ std::cout << "            Selecting Propagation Values for Surface Deflections" << std::endl; }

            controlSurfaceDeflections = bislip::Variables::computeControlSurfaceDeflection( bodyMap_, vehicleName_ );
        }

        bislipSystems->setCurrentBodyFlapAngle( controlSurfaceDeflections( 0 ) );
        bislipSystems->setCurrentElevonAngle( controlSurfaceDeflections( 1 ) );

        if( debugInfo == 1 ){ std::cout << "           Assigning Surface Deflections" << std::endl; }
        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", bislipSystems->getCurrentBodyFlapAngle() );
        vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonLeft", bislipSystems->getCurrentElevonAngle() );
        vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonRight", bislipSystems->getCurrentElevonAngle() );

    //}


    if( debugInfo == 1 ){ std::cout << "Control surface deflections have been assigned" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "Determine and set Full Coefficients" << std::endl; }

    bislipSystems->setFullCurrentCoefficients( bislip::Variables::computeFullCurrentCoefficients( bodyMap_, vehicleName_ ) );

    if( debugInfo == 1 ){ std::cout << "Determine and set Lift Force" << std::endl; }

    bislipSystems->setCurrentLiftForce( bislip::Variables::computeCurrentLiftForce( bodyMap_, vehicleName_ ) );
    if( debugInfo == 1 ){ std::cout << "    Calculated Lift Force = " << bislipSystems->getCurrentLiftForce() << std::endl; }

    if( debugInfo == 1 ){ std::cout << "Determine and set Drag Force" << std::endl; }

    bislipSystems->setCurrentDragForce( bislip::Variables::computeCurrentDragForce( bodyMap_, vehicleName_ ) );
    if( debugInfo == 1 ){ std::cout << "    Calculated Drag Force = " << bislipSystems->getCurrentDragForce() << std::endl; }



/*


    double bankAngleForSkipSuppression = 0.0;
    bislipSystems->setCurrentBankAngleLimit( bankAngleForSkipSuppression );

    if( debugInfo == 1 )
    {
        std::cout << "     Current Trajectory Phrase: " << currentTrajectoryPhase << std::endl;
        if( currentTrajectoryPhase == "Ascent" ){ std::cout << "          No bank angle computation." << std::endl; }
    }

    double bankAngle = 0.0;

    if ( currentTrajectoryPhase == "Descent" )
    {
        if( debugInfo == 1 ){ std::cout << "          Current Dynamic Pressure = " << currentDynamicPressure << std::endl; }
        if( debugInfo == 1 ){ std::cout << "          Minimum Dynamic Pressure = " << bislipSystems->getMinimumDynamicPressureforControlSurface() << std::endl; }

        if( currentDynamicPressure > bislipSystems->getMinimumDynamicPressureforControlSurface() )
        {

            bislipSystems->setEvaluatedBankAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::BankAngle, bodyMap_, vehicleName_, centralBodyName_ ) ) );
            if( debugInfo == 1 ){ std::cout << "    Evaluated Bank Angle = " << bislipSystems->getEvaluatedBankAngle() << std::endl; }

            bankAngle = bislip::Variables::computeBankAngle( bodyMap_, vehicleName_ );

            if( debugInfo == 1 ){ std::cout << "          Adjusted Bank Angle = " <<  tudat::unit_conversions::convertRadiansToDegrees( bankAngle ) << std::endl; }

            if( ( FlightConditions_->getCurrentTime() > startingEpoch_ + bislipSystems->getSkipSuppressionTimingTrigger( ) ) )
            {
                if( debugInfo == 1 ){ std::cout << "               Skip Suppression is Active" << std::endl; }
                if( debugInfo == 1 ){ std::cout << "                    Current Flight-Path Angle = " <<  currentFlightPathAngle << " rad" << std::endl; }

                //double currentFlightPathAngleRate = bislip::Variables::computeFlightPathAngleRate( bodyMap_, vehicleName_, "Earth" );

                //if ( currentFlightPathAngle >= 0 )
                //{
                if( debugInfo == 1 ){ std::cout << "                    Calculate Bank Angle Required to Suppress Skipping" << std::endl; }

                bankAngleForSkipSuppression = double( bislip::Variables::determineSignOfValue( bankAngle ) ) * std::abs( bislip::Variables::computeSkipSuppressionLimit( bodyMap_, vehicleName_, centralBodyName_ ) );

                bislipSystems->setCurrentBankAngleLimit( bankAngleForSkipSuppression );

                if( debugInfo == 1 ){ std::cout << "                    Current Bank Angle = " <<  tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentBankAngle( ) ) << " deg  |  " << "Calculated Bank Angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( bankAngle ) << " deg  |  " << "Bank Angle Required = " <<  tudat::unit_conversions::convertRadiansToDegrees( bankAngleForSkipSuppression ) << " deg" << std::endl; }

                //bankAngle = bankAngleForSkipSuppression;
                if ( std::abs( bankAngleForSkipSuppression ) > std::abs( bankAngle ) ) { bankAngle = bankAngleForSkipSuppression; }
                if( debugInfo == 1 ){ std::cout << "                    NEW Bank Angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( bankAngle ) << std::endl; }

                //}
            }
        }
    }
    bislipSystems->setCurrentBankAngle( bankAngle );
    currentBankAngle_ = bislipSystems->getCurrentBankAngle();


*/






    double evaluatedBankAngle       = 0.0;
    double skipSuppressionBankAngle = 0.0;
    double greaterBankAngle         = 0.0;
    double signedGreaterBankAngle   = 0.0;
    double bankAngle                = 0.0;

    if( debugInfo == 1 )
    {
        std::cout << "     Current Trajectory Phrase: " << currentTrajectoryPhase << std::endl;
        if( currentTrajectoryPhase == "Ascent" ){ std::cout << "          No bank angle computation." << std::endl; }
    }

    if ( currentTrajectoryPhase == "Descent" )
    {

        if( debugInfo == 1 ){ std::cout << "Compute and assign the Bank Angle" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "    Evaluating relevant interpolator" << std::endl; }
        evaluatedBankAngle = tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Interpolators::BankAngle, bodyMap_, vehicleName_, centralBodyName_ ) );

        if( ( FlightConditions_->getCurrentTime() > startingEpoch_ + bislipSystems->getSkipSuppressionTimingTrigger( ) ) )
        {
            if( debugInfo == 1 ){ std::cout << "               Skip Suppression is Active" << std::endl; }
            if( debugInfo == 1 ){ std::cout << "                    Calculate Bank Angle Required to Suppress Skipping" << std::endl; }

            skipSuppressionBankAngle = bislip::Variables::computeSkipSuppressionLimit( bodyMap_, vehicleName_, centralBodyName_ );
        }

        if( debugInfo == 1 ){ std::cout << "                    Choose the greater of these" << std::endl; }

        if( std::abs( evaluatedBankAngle ) > std::abs( skipSuppressionBankAngle ) ) { greaterBankAngle = evaluatedBankAngle; }
        else { greaterBankAngle = skipSuppressionBankAngle; }

        //greaterBankAngle = std::max( evaluatedBankAngle, skipSuppressionBankAngle );
        //std::cout << "                      Evaluated Bank Angle = " << evaluatedBankAngle << std::endl;
        //std::cout << "                      Limited Bank Angle   = " << skipSuppressionBankAngle << std::endl;
        //std::cout << "                      Greater Bank Angle   = " << greaterBankAngle << std::endl;
        if( std::abs( greaterBankAngle ) > 0.0 )
        {
            if( debugInfo == 1 ){ std::cout << "        Assigning the SIGN of the previous bank angle" << std::endl; }
            bankAngle = bislip::Variables::determineSignOfValue( previousBankAngle ) * greaterBankAngle;
        }

        if( debugInfo == 1 ){ std::cout << "    Calculate Heading Error" << std::endl; }
        const double currentHeadingAngleError = bislip::Variables::computeHeadingError( currentLatitude, currentLongitude, bislipSystems->getTargetLat(), bislipSystems->getTargetLon(), currentHeadingAngle );

        if( debugInfo == 1 ){ std::cout << "    Calculate Heading Error Deadband" << std::endl; }
        const double currentHeadingAngleErrorDeadBand = tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::computeHeadingErrorDeadBand( bodyMap_, vehicleName_ ) );

        if( debugInfo == 1 ){ std::cout << "    Calculate Angular Distance To Go" << std::endl; }
        const double angularDistanceToGo_deg = tudat::unit_conversions::convertRadiansToDegrees(
                    bislip::Variables::computeAngularDistance (
                        currentLatitude,
                        currentLongitude,
                        bislipSystems->getTargetLat(),
                        bislipSystems->getTargetLon() ) );

        //! Determine if vehicle is aiming away from target.
        //! Determine bank angle Reversal conditional. If result is lower than zero,
        //! bank reversal could happen. Should be less than zero if the vehicle is
        //! in fact aiming away from the target.

        //! This will happen if the bank angle is positive and the heading error is
        //! negative or if the bank angle is negative and the heading error is positive.

        if( debugInfo == 1 ){ std::cout << "    Calculate Reversal Conditional" << std::endl; }
        const double reversalConditional = bankAngle * currentHeadingAngleError;

        //! Determine if bank reversal should occur.
        bool reversal = false;
        if( reversalConditional < 0.0 )
        {
            if( debugInfo == 1 ){ std::cout << "        Reversal Conditional is NEGATIVE" << std::endl; }

            if( std::abs( currentHeadingAngleError ) >= currentHeadingAngleErrorDeadBand )
            {
                if( debugInfo == 1 ){ std::cout << "        Heading Error is GREATER than Deadband" << std::endl; }
                if( bislipSystems->getLowDistanceReversalCompleted( ) == false )
                {
                    if( debugInfo == 1 ){ std::cout << "        Low Distance Reversal has NOT been done yet" << std::endl; }
                    reversal = true;

                    if( angularDistanceToGo_deg < bislipSystems->getHeadingErrorDeadBandLowDistanceTrigger() )
                    {
                        if( debugInfo == 1 ){ std::cout << "        Low Distance Reversal has been completed" << std::endl; }
                        bislipSystems->setLowDistanceReversalCompleted( true );
                    }
                }
                else
                { if( debugInfo == 1 ){ std::cout << "        Low Distance Reversal has been satisfied" << std::endl; } }
            }
        }

        //! Impose bank angle reversal on new bank angle.
        if( reversal == true )
        {
            if( debugInfo == 1 ){ std::cout << "        Imposing Bank Angle Reversal " << std::endl; }
            bankAngle = -bankAngle;
        }

        //if( debugInfo == 1 ){ std::cout << "timeOfFLight = " << timeOfFLight << "  |  " << "Evaluated Bank Angle = " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getEvaluatedBankAngle() ) << "  ||  " << "Current Heading Error = " << currentHeadingAngleError << "  |  " << "Reversal Conditional = " << reversalConditional << "  |  " << " DeadBand = " << tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::computeHeadingErrorDeadBand( bodyMap_, vehicleName_ ) ) << "  |  " << "Reversal = " << reversal << std::endl; }

        //std::cout << "timeOfFLight = " << timeOfFLight << "  |  " << "Evaluated Bank Angle = " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getEvaluatedBankAngle() ) << "  |  " << "Current Heading Error = " << currentHeadingAngleError << "  |  " << "Reversal Conditional = " << reversalConditional << "  |  " << " DeadBand = " << tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::computeHeadingErrorDeadBand( bodyMap_, vehicleName_ ) ) << "  |  " << "Reversal = " << reversal << std::endl;
        //std::cout << timeOfFLight << "  " << tudat::unit_conversions::convertRadiansToDegrees( previousBankAngle ) << "  " << tudat::unit_conversions::convertRadiansToDegrees( evaluatedBankAngle ) << "  " << tudat::unit_conversions::convertRadiansToDegrees( signedEvaluatedBankAngle ) << "  " << tudat::unit_conversions::convertRadiansToDegrees( currentHeadingAngleError ) << " " << reversalConditional << " " <<  bislip::Variables::computeHeadingErrorDeadBand( bodyMap_, vehicleName_ ) << "  " << reversal << " " << tudat::unit_conversions::convertRadiansToDegrees( reversedBankAngle ) << std::endl;

        //! SIGN of bank angle is determined in here.
       // bankAngle = reversedBankAngle;//bislip::Variables::computeBankAngle( bodyMap_, vehicleName_ );

        if( debugInfo == 1 ){ std::cout << "          Calculated Bank Angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( bankAngle ) << std::endl; }
    }

    bislipSystems->setEvaluatedBankAngle( evaluatedBankAngle );
    bislipSystems->setCurrentBankAngleLimit( skipSuppressionBankAngle );
    bislipSystems->setCurrentBankAngle( bankAngle );
    currentBankAngle_ = bislipSystems->getCurrentBankAngle( );








































    currentAngleOfSideslip_ = 0.0;

    if( debugInfo == 1 ){ std::cout << "Guidance Angles are Set" << std::endl; }

    if( debugInfo == 1 ){ std::cout << "Determine and set engine related values by evaluating various functions." << std::endl; }

    if( bislipSystems->getCurrentEngineStatus( ) == true )
    {
        if( debugInfo == 1 ){ std::cout << "    Engine is ON" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "    Determine and Set Throttle Setting" << std::endl; }
        if( debugInfo == 1 ){ std::cout << "    --------------------------" << std::endl; }
        bislipSystems->setCurrentThrottleSetting( bislip::Variables::determineThrottleSetting( bodyMap_, vehicleName_, centralBodyName_ ) );
        if( debugInfo == 1 ){ std::cout << "    --------------------------" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "    Determine and Set Thrust Magnitude" << std::endl; }
        if( debugInfo == 1 ){ std::cout << "    --------------------------" << std::endl; }
        bislipSystems->setCurrentThrustMagnitude( bislip::Variables::computeThrustMagnitude( bodyMap_, vehicleName_ ) );
        if( debugInfo == 1 ){ std::cout << "    --------------------------" << std::endl; }

        /*
        double dif = bislipSystems->getCurrentThrustMagnitude()*std::cos(currentAngleOfAttack_) - bislipSystems->getCurrentDragForce();
        double weight = bodyMap_.at( vehicleName_ )->getBodyMass()*tudat::mathematical_constants::PI;
        double singamma = dif/weight;

        std::cout << "dif        = " << dif << "    |    weight     = " << weight<< std::endl;
        std::cout << "sin(gamma) = " << singamma << "    |    gamma      = " << std::asin(singamma )<< std::endl;
*/
        if( debugInfo == 1 ){ std::cout << "    Determine and Set Thrust Elevation Angle" << std::endl; }
        if( debugInfo == 1 ){ std::cout << "    --------------------------" << std::endl; }
        //if ( currentDynamicPressure < bislipSystems->getMinimumDynamicPressureforControlSurface() )
        // {
        //if( debugInfo == 1 ){ std::cout << "        Control Surfaces are not engaged, which requires alternate means of trim" << std::endl; }

        // bislipSystems->setCurrentThrustElevationAngle( bislip::Variables::determineThrustElevationAngle( bodyMap_, vehicleName_, centralBodyName_ ) );
        bislipSystems->setCurrentThrustElevationAngle( 0.0 );

        bislipSystems->setCurrentThrustAzimuthAngle( 0.0 );
        // }
        // else
        // {
        //   if( debugInfo == 1 ){ std::cout << "        Control Surfaces are engaged, nothing do to here" << std::endl; }
        //  bislipSystems->setCurrentThrustElevationAngle( bislip::Variables::determineThrustElevationAngle( bodyMap_, vehicleName_, centralBodyName_ ) );

        // bislipSystems->setCurrentThrustElevationAngle( 0.0 );
        //bislipSystems->setCurrentThrustAzimuthAngle( 0.0 );

        ///}
        if( debugInfo == 1 ){ std::cout << "    Determine and Set Thrust Direction" << std::endl; }
        if( debugInfo == 1 ){ std::cout << "    --------------------------" << std::endl; }
        bislipSystems->setCurrentBodyFixedThrustDirection( bislip::Variables::computeBodyFixedThrustDirection( bodyMap_, vehicleName_ ) );
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "    Engine if OFF" << std::endl; }
        bislipSystems->setCurrentThrottleSetting( 0.0 );
        bislipSystems->setCurrentThrustMagnitude( 0.0 );
        bislipSystems->setCurrentThrustElevationAngle( 0.0 );
        bislipSystems->setCurrentBodyFixedThrustDirection( Eigen::Vector3d::Zero( 3 ) );
    }


    //************************************************** Bank Angle
    if( debugInfo == 1 ){ std::cout << "Compute, assign, bound, and set the current bank angle." << std::endl; }






}

void MyGuidance::evaluateValidationGuidanceFunctions(
        // std::shared_ptr< tudat::aerodynamics::AtmosphericFlightConditions > &flightConditions,
        std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems,
        std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
        const std::string &currentTrajectoryPhase,
        const double &timeOfFLight )
{

    int debugInfo = bislipSystems->getDebugInfo();

    double currentMachNumber      = 0.0;
    double currentHeight          = 0.0;
    double currentAirspeed        = 0.0;
    double currentDynamicPressure = 0.0;
    double currentFlightPathAngle = 0.0;
    double currentHeadingAngle    = 0.0;
    double currentLatitude        = 0.0;
    double currentLongitude       = 0.0;
    double previousBankAngle      = 0.0;
    double timeOfFlight           = 0.0;

    if( debugInfo == 1 ){ std::cout << "      Selecting Source of Values" << std::endl; }

    if( bislipSystems->getInitialValueFlag() == true )
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values" << std::endl; }

        currentMachNumber      = bislipSystems->getInitialMachNumber();
        currentHeight          = bislipSystems->getInitialHeight();
        currentAirspeed        = bislipSystems->getInitialAirspeed();
        currentDynamicPressure = bislipSystems->getInitialDynamicPressure();
        currentFlightPathAngle = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialFlightPathAngle() );
        currentHeadingAngle    = tudat::unit_conversions::convertDegreesToRadians( bislipSystems->getInitialHeading() );
        currentLatitude        = bislipSystems->getInitialLat();
        currentLongitude       = bislipSystems->getInitialLon();
        previousBankAngle      = 0.0;
        timeOfFlight           = 0;

        //bislipSystems->getCurrentMachNumber()
    }
    else
    {
        if( debugInfo == 1 ){ std::cout << "            Selecting Propagated Values" << std::endl; }

        currentMachNumber      = FlightConditions_->getCurrentMachNumber();
        currentHeight          = FlightConditions_->getCurrentAltitude();
        currentAirspeed        = FlightConditions_->getCurrentAirspeed();
        currentDynamicPressure = FlightConditions_->getCurrentDynamicPressure();
        currentFlightPathAngle = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        currentHeadingAngle    = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );
        currentLatitude        = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        currentLongitude       = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );
        previousBankAngle      = bislipSystems->getCurrentBankAngle();
        timeOfFlight           = FlightConditions_->getCurrentTime() - bislipSystems->getStartingEpoch();

    }

    if( debugInfo == 1 ){ std::cout << " "<< std::endl; }
    if( debugInfo == 1 ){ std::cout << "Validation case?          = " << bislipSystems->getValidationFlag() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Trajectory Phase          = " << bislipSystems->getCurrentTrajectoryPhase() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Time of Flight            = " << timeOfFlight << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Height            = " << currentHeight << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Airspeed          = " << currentAirspeed << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Dynamic Pressure  = " << currentDynamicPressure << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Mach Number       = " << currentMachNumber << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Flight-path Angle = " << tudat::unit_conversions::convertRadiansToDegrees( currentFlightPathAngle ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Heading Angle     = " << tudat::unit_conversions::convertRadiansToDegrees( currentHeadingAngle ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Mass              = " << bodyMap_.at( vehicleName_ )->getBodyMass() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Current Engine Status     = " << bislipSystems->getCurrentEngineStatus() << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Starting Angle of Attack  = " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentAngleOfAttack() ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "Starting Bank Angle       = " << tudat::unit_conversions::convertRadiansToDegrees( previousBankAngle ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << " "<< std::endl; }


    double evaluatedAngleOfAttack = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getKourouAngleOfAttackInterpolator( ) )->interpolate( timeOfFLight ) );

    if( debugInfo == 1 ){ std::cout << "Extract bounds of angle of attack as a function of Mach Number" << std::endl; }
    double angleOfAttackLowerBound = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getAlphaMachEnvelopeLBInterpolator() )->interpolate( currentMachNumber ) );
    double angleOfAttackUpperBound = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getAlphaMachEnvelopeUBInterpolator() )->interpolate( currentMachNumber ) );
    if( debugInfo == 1 ){ std::cout << "    Angle of attack Lower Bound due to Mach Number = " << tudat::unit_conversions::convertRadiansToDegrees( angleOfAttackLowerBound ) << std::endl; }
    if( debugInfo == 1 ){ std::cout << "    Angle of attack Upper Bound due to Mach Number = " << tudat::unit_conversions::convertRadiansToDegrees( angleOfAttackUpperBound ) << std::endl; }

    if( evaluatedAngleOfAttack >= angleOfAttackUpperBound ) { evaluatedAngleOfAttack = angleOfAttackUpperBound; }
    else if( evaluatedAngleOfAttack <= angleOfAttackLowerBound ) { evaluatedAngleOfAttack = angleOfAttackLowerBound; }


    bislipSystems->setCurrentAngleOfAttack( evaluatedAngleOfAttack );
    currentAngleOfAttack_ = bislipSystems->getCurrentAngleOfAttack( );




    //************************************************** Control Surface Deflections
    //! Compute and assign the current control surface deflection angles

    Eigen::Vector2d controlSurfaceDeflections;

    if( debugInfo == 1 ){ std::cout << "Compute and assign the current control surface deflections." << std::endl; }

    if( debugInfo == 1 ){ std::cout << "    Determine Control Surface Deflections for Validation Case" << std::endl; }

    // if( timeOfFLight >= 194)
    // {
    if( debugInfo == 1 ){ std::cout << "          Current Dynamic Pressure = " << currentDynamicPressure << std::endl; }
    if( debugInfo == 1 ){ std::cout << "          Minimum Dynamic Pressure = " << bislipSystems->getMinimumDynamicPressureforControlSurface() << std::endl; }

    if( currentDynamicPressure > bislipSystems->getMinimumDynamicPressureforControlSurface() )
    {

        if( bislipSystems->getInitialValueFlag() == true )
        {
            if( debugInfo == 1 ){ std::cout << "            Selecting Initial Values for Surface Deflections" << std::endl; }

            controlSurfaceDeflections << bislipSystems->getCurrentBodyFlapAngle(), bislipSystems->getCurrentElevonAngle(); }
        else
        {
            if( debugInfo == 1 ){ std::cout << "            Selecting Propagation Values for Surface Deflections" << std::endl; }

            controlSurfaceDeflections = bislip::Variables::computeControlSurfaceDeflection( bodyMap_, vehicleName_ );
        }

        bislipSystems->setCurrentBodyFlapAngle( controlSurfaceDeflections( 0 ) );
        bislipSystems->setCurrentElevonAngle( controlSurfaceDeflections( 1 ) );

        if( debugInfo == 1 ){ std::cout << "           Assigning Surface Deflections" << std::endl; }
        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", bislipSystems->getCurrentBodyFlapAngle() );
        vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonLeft", bislipSystems->getCurrentElevonAngle() );
        vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonRight", bislipSystems->getCurrentElevonAngle() );

    }

    double evaluatedBankAngle = 0.0;
    double signedEvaluatedBankAngle = 0.0;
    double bankAngle = 0.0;

    if ( currentTrajectoryPhase == "Descent" )
    {

        if( debugInfo == 1 ){ std::cout << "Compute and assign the Bank Angle" << std::endl; }

        if( debugInfo == 1 ){ std::cout << "    Evaluating relevant interpolator" << std::endl; }
        evaluatedBankAngle = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getKourouBankAngleInterpolator( ) )->interpolate( timeOfFLight ) );

        if( std::abs( evaluatedBankAngle ) > 1e-9 )
        {
            if( debugInfo == 1 ){ std::cout << "        Assigning the SIGN of the previous bank angle" << std::endl; }
            signedEvaluatedBankAngle = bislip::Variables::determineSignOfValue( previousBankAngle ) * evaluatedBankAngle;
        }

        if( debugInfo == 1 ){ std::cout << "    Calculate Heading Error" << std::endl; }
        const double currentHeadingAngleError = bislip::Variables::computeHeadingError( currentLatitude, currentLongitude, bislipSystems->getTargetLat(), bislipSystems->getTargetLon(), currentHeadingAngle );

        if( debugInfo == 1 ){ std::cout << "    Calculate Heading Error Deadband" << std::endl; }
        const double currentHeadingAngleErrorDeadBand = tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::computeHeadingErrorDeadBand( bodyMap_, vehicleName_ ) );

        if( debugInfo == 1 ){ std::cout << "    Calculate Angular Distance To Go" << std::endl; }
        const double angularDistanceToGo_deg = tudat::unit_conversions::convertRadiansToDegrees(
                    bislip::Variables::computeAngularDistance (
                        currentLatitude,
                        currentLongitude,
                        bislipSystems->getTargetLat(),
                        bislipSystems->getTargetLon() ) );

        //! Determine if vehicle is aiming away from target.
        //! Determine bank angle Reversal conditional. If result is lower than zero,
        //! bank reversal could happen. Should be less than zero if the vehicle is
        //! in fact aiming away from the target.

        //! This will happen if the bank angle is positive and the heading error is
        //! negative or if the bank angle is negative and the heading error is positive.

        if( debugInfo == 1 ){ std::cout << "    Calculate Reversal Conditional" << std::endl; }
        const double reversalConditional = signedEvaluatedBankAngle * currentHeadingAngleError;

        //! Determine if bank reversal should occur.
        bool reversal = false;
        if( reversalConditional < 0.0 )
        {
            if( debugInfo == 1 ){ std::cout << "        Reversal Conditional is NEGATIVE" << std::endl; }

            if( std::abs( currentHeadingAngleError ) >= currentHeadingAngleErrorDeadBand )
            {
                if( debugInfo == 1 ){ std::cout << "        Heading Error is GREATER than Deadband" << std::endl; }
                if( bislipSystems->getLowDistanceReversalCompleted( ) == false )
                {
                    if( debugInfo == 1 ){ std::cout << "        Low Distance Reversal has NOT been done yet" << std::endl; }
                    reversal = true;

                    if( angularDistanceToGo_deg < bislipSystems->getHeadingErrorDeadBandLowDistanceTrigger() )
                    {
                        if( debugInfo == 1 ){ std::cout << "        Low Distance Reversal has been completed" << std::endl; }
                        bislipSystems->setLowDistanceReversalCompleted( true );
                    }
                    else
                    { if( debugInfo == 1 ){ std::cout << "        Low Distance Reversal has been satisfied" << std::endl; } }
                }
            }
        }

        double reversedBankAngle = signedEvaluatedBankAngle;
        //! Impose bank angle reversal on new bank angle.
        if( reversal == true )
        {
            if( debugInfo == 1 ){ std::cout << "        Imposing Bank Angle Reversal " << std::endl; }
            reversedBankAngle = -reversedBankAngle;
        }

        //if( debugInfo == 1 ){ std::cout << "timeOfFLight = " << timeOfFLight << "  |  " << "Evaluated Bank Angle = " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getEvaluatedBankAngle() ) << "  ||  " << "Current Heading Error = " << currentHeadingAngleError << "  |  " << "Reversal Conditional = " << reversalConditional << "  |  " << " DeadBand = " << tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::computeHeadingErrorDeadBand( bodyMap_, vehicleName_ ) ) << "  |  " << "Reversal = " << reversal << std::endl; }

        //std::cout << "timeOfFLight = " << timeOfFLight << "  |  " << "Evaluated Bank Angle = " << tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getEvaluatedBankAngle() ) << "  |  " << "Current Heading Error = " << currentHeadingAngleError << "  |  " << "Reversal Conditional = " << reversalConditional << "  |  " << " DeadBand = " << tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::computeHeadingErrorDeadBand( bodyMap_, vehicleName_ ) ) << "  |  " << "Reversal = " << reversal << std::endl;
        //std::cout << timeOfFLight << "  " << tudat::unit_conversions::convertRadiansToDegrees( previousBankAngle ) << "  " << tudat::unit_conversions::convertRadiansToDegrees( evaluatedBankAngle ) << "  " << tudat::unit_conversions::convertRadiansToDegrees( signedEvaluatedBankAngle ) << "  " << tudat::unit_conversions::convertRadiansToDegrees( currentHeadingAngleError ) << " " << reversalConditional << " " <<  bislip::Variables::computeHeadingErrorDeadBand( bodyMap_, vehicleName_ ) << "  " << reversal << " " << tudat::unit_conversions::convertRadiansToDegrees( reversedBankAngle ) << std::endl;

        //! SIGN of bank angle is determined in here.
        bankAngle = reversedBankAngle;//bislip::Variables::computeBankAngle( bodyMap_, vehicleName_ );

        if( debugInfo == 1 ){ std::cout << "          Calculated Bank Angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( bankAngle ) << std::endl; }
    }

    bislipSystems->setEvaluatedBankAngle( evaluatedBankAngle );
    bislipSystems->setCurrentBankAngle( bankAngle );
    currentBankAngle_ = bislipSystems->getCurrentBankAngle( );


    currentAngleOfSideslip_ = 0.0;

    bislipSystems->setCurrentEngineStatus( false );
    bislipSystems->setCurrentThrottleSetting( 0.0 );
    bislipSystems->setCurrentThrustMagnitude( 0.0 );
    bislipSystems->setCurrentThrustElevationAngle( 0.0 );
    bislipSystems->setCurrentThrustAzimuthAngle( 0.0 );


}





bool MyGuidance::getCurrentEngineStatus( )
{
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();
    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "Extracting the Current Engine Status" << std::endl; }

    return bislipSystems->getCurrentEngineStatus();
}

double MyGuidance::getCurrentThrustMagnitude( )
{
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();
    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "Extracting the Current Thrust Magnitude" << std::endl; }


    return bislipSystems->getCurrentThrustMagnitude();
}

Eigen::Vector3d MyGuidance::getCurrentBodyFixedThrustDirection( )
{
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();

    int debugInfo = bislipSystems->getDebugInfo();

    if( debugInfo == 1 ){ std::cout << "Extracting the Current Body-Fixed Thrust Direction" << std::endl; }

    return bislipSystems->getCurrentBodyFixedThrustDirection();
}
/*
Eigen::Vector6d MyGuidance::getPartialCurrentCoefficients( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    if( coefficientInterface_ == nullptr )
    {
        coefficientInterface_ = std::dynamic_pointer_cast< tudat::aerodynamics::AerodynamicCoefficientInterface >(
                    bodyMap_.at( vehicleName_ )->getAerodynamicCoefficientInterface( ) );
    }

    if( bislipSystems_ == nullptr )
    {
        bislipSystems_ = std::dynamic_pointer_cast< bislip::BislipVehicleSystems >(
                    bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
    }

    return bislip::Variables::computePartialCurrentCoefficients( bodyMap_, vehicleName_ );
}
*/



} // namespace bislip
