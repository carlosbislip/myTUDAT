#include "updateGuidance.h"

namespace bislip {

void MyGuidance::updateGuidance( const double currentTime )
{
    //std::cout << "Starting Aerodynamic guidance for this evaluation" << std::endl;
    std::shared_ptr< bislip::BislipVehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();
    std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap_.at( vehicleName_ )->getVehicleSystems();
    double startingEpoch = bislipSystems->getStartingEpoch();
    //double propagationStepSize = bislipSystems->getPropagationStepSize();
    double guidanceStepSize = bislipSystems->getGuidanceStepSize();
    int debugInfo = bislipSystems->getDebugInfo();
    //double samplingRatio = bislipSystems->getSamplingRatio();
    std::string currentTrajectoryPhase = bislipSystems->getCurrentTrajectoryPhase();

    if( debugInfo == 4 ){std::cout << "currentTime = " << currentTime << std::endl; }
    if( debugInfo == 4 ){std::cout << "currentTime - startingEpoch = " <<  currentTime - startingEpoch << std::endl; }
    if( debugInfo == 4 ){std::cout << "std::fmod( currentTime - startingEpoch, guidanceStepSize ) = " <<  std::fmod( currentTime - startingEpoch, guidanceStepSize )<< std::endl; }
    if( debugInfo == 4 ){std::cout << "(std::trunc((std::fmod( currentTime - startingEpoch, guidanceStepSize ) * 100)))/100 = " <<  (std::trunc((std::fmod( currentTime - startingEpoch, guidanceStepSize ) * 100)))/100 << std::endl; }

/*
    if ( bislipSystems->getValidationFlag( ) == true )
    {



        currentAngleOfAttack_ = ( bislipSystems->getKourouAngleOfAttackInterpolator( ) )->interpolate( timeOfFLight );

        Eigen::Vector2d controlSurfaceDeflections = bislip::Variables::computeControlSurfaceDeflection( bodyMap_, vehicleName_ );

        bislipSystems->setCurrentBodyFlapAngle( controlSurfaceDeflections( 0 ) );
        bislipSystems->setCurrentElevonAngle( controlSurfaceDeflections( 1 ) );
        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", bislipSystems->getCurrentBodyFlapAngle() );
        vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonLeft", bislipSystems->getCurrentElevonAngle() );
        vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonRight", bislipSystems->getCurrentElevonAngle() );

        currentBankAngle_ = ( bislipSystems->getKourouBankAngleInterpolator( ) )->interpolate( timeOfFLight );
        bislipSystems->setTempBankAngle( currentBankAngle_ );




        currentAngleOfSideslip_ = 0.0;

    }
    else
    {
*/

        if ( currentTime == currentTime )
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



        }

        const double timeOfFLight = std::trunc( ( currentTime - startingEpoch ) * 100 ) / 100;

        const double updateConditionalPre = std::fmod( timeOfFLight, guidanceStepSize );
        const double updateConditional = std::trunc( updateConditionalPre * 1000 ) / 1000;


        if( debugInfo == 2 ){
            std::cout << "Choosing which guidance conditional to run" << std::endl;
            std::cout << "updateConditional: " << updateConditional << "  |  " << "Time of flight = " << timeOfFLight << std::endl;
        }

        if( ( currentTime > startingEpoch ) && ( updateConditional == 0.0 ) )
        {
            if( debugInfo == 2 ){std::cout << "Actual Guidance Conditional" << std::endl; }
            if( debugInfo == 2 ){std::cout << " " << std::endl; }
            if( debugInfo == 2 ){std::cout << " " << std::endl; }
            if( debugInfo == 2 ){std::cout << " " << std::endl; }

            evaluateGuidanceFunctions( bislipSystems, vehicleSystems, currentTrajectoryPhase );

            if( debugInfo == 2 ){ std::cout << "          Saving bank reversal trigger if applicable."<< std::endl; }
            if( debugInfo == 2 ){ std::cout << "            FlightConditions_->getCurrentTime() = "<< FlightConditions_->getCurrentTime() << std::endl; }
            if( debugInfo == 2 ){ std::cout << "            bislipSystems->getBankAngleReversalTimepoint() = "<< bislipSystems->getBankAngleReversalTimepoint() << std::endl; }

            if ( FlightConditions_->getCurrentTime() != bislipSystems->getBankAngleReversalTimepoint() )
            {
                if ( std::isnan( FlightConditions_->getCurrentTime() ) == false  )
                {
                    if( debugInfo == 2 ){ std::cout << "               bislipSystems->getCurrentBankAngle() * bislipSystems->getTempBankAngle( ) = "<< bislipSystems->getCurrentBankAngle() * bislipSystems->getTempBankAngle( ) <<std::endl; }

                    if ( bislipSystems->getCurrentBankAngle() * bislipSystems->getTempBankAngle( ) < 0.0 )
                    {
                        bislipSystems->setBankAngleReversalTrigger( true );
                        if( debugInfo == 2 ){ std::cout << "               bislipSystems->getBankAngleReversalTrigger = "<< bislipSystems->getBankAngleReversalTrigger() <<std::endl; }

                        bislipSystems->setBankAngleReversalTimepoint( FlightConditions_->getCurrentTime() );
                        if( debugInfo == 2 ){ std::cout << "               NEW bislipSystems->getBankAngleReversalTimepoint = "<< bislipSystems->getBankAngleReversalTimepoint() <<std::endl; }
                    }
                    else { bislipSystems->setBankAngleReversalTrigger( false ); }
                }
            }

            if ( currentTime == currentTime )
            {
            //************************************************** Heat Flux(es)
            if( debugInfo == 1 ){ std::cout << "Determine and set heat fluxes." << std::endl; }
            if( debugInfo == 1 ){ std::cout << "     computeStagnationHeatFlux" << std::endl; }
            bislipSystems->setWorkingRadius( bislipSystems->getNoseRadius() );
            bislipSystems->setWallTemperature( bislipSystems->getChapmanWallTemp() );

            if( debugInfo == 40 ){ std::cout << "     Working Radius = " << bislipSystems->getWorkingRadius() << std::endl; }
            if( debugInfo == 40 ){ std::cout << "     Initial Wall Temperature = " <<  bislipSystems->getWallTemperature( ) << std::endl; }
            bislipSystems->setCurrentHeatFluxChapman( bislip::Variables::computeHeatingRateChapman( bodyMap_, vehicleName_) );
            bislipSystems->setChapmanWallTemp( bislipSystems->getWallTemperature() );


            if( debugInfo == 40 ){ std::cout << "     Current Stagnation Heat Flux = " <<  bislipSystems->getCurrentHeatFluxChapman( ) << std::endl; }
            if( debugInfo == 40 ){ std::cout << "     Current Stagnation Eq. Temp. = " <<  bislipSystems->getChapmanWallTemp( ) << std::endl; }

            if( debugInfo == 1 ){ std::cout << "     computeHeatingRateTauber" << std::endl; }
            bislipSystems->setWorkingRadius( bislipSystems->getLeadingEdgeRadius() );

            if( debugInfo == 1 ){ std::cout << "     Working Radius = " << bislipSystems->getWorkingRadius() << std::endl; }
            bislipSystems->setCurrentHeatFluxTauber( bislip::Variables::computeHeatingRateTauber( bodyMap_, vehicleName_) );

            if( debugInfo == 1 ){ std::cout << "     Current Tauber Heat Flux = " <<  bislipSystems->getCurrentHeatFluxTauber( ) << std::endl; }
}
        }
        else
        {
            if( debugInfo == 2 ){std::cout << "No Guidance Conditional" << std::endl; }

            currentAngleOfAttack_ = bislipSystems->getCurrentAngleOfAttack();

            vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", bislipSystems->getCurrentBodyFlapAngle() );
            vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonLeft", bislipSystems->getCurrentElevonAngle() );
            vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonRight", bislipSystems->getCurrentElevonAngle() );

            currentBankAngle_ = bislipSystems->getCurrentBankAngle();

            currentAngleOfSideslip_ = 0.0;

            if ( currentTime == currentTime )
            {
            //************************************************** Heat Flux(es)
            if( debugInfo == 1 ){ std::cout << "Determine and set heat fluxes." << std::endl; }
            if( debugInfo == 1 ){ std::cout << "     computeStagnationHeatFlux" << std::endl; }
            bislipSystems->setWorkingRadius( bislipSystems->getNoseRadius() );
            bislipSystems->setWallTemperature( bislipSystems->getChapmanWallTemp() );

            if( debugInfo == 40 ){ std::cout << "     Working Radius = " << bislipSystems->getWorkingRadius() << std::endl; }
            if( debugInfo == 40 ){ std::cout << "     Initial Wall Temperature = " <<  bislipSystems->getWallTemperature( ) << std::endl; }
            bislipSystems->setCurrentHeatFluxChapman( bislip::Variables::computeHeatingRateChapman( bodyMap_, vehicleName_) );
            bislipSystems->setChapmanWallTemp( bislipSystems->getWallTemperature() );
            if( debugInfo == 40 ){ std::cout << "     Current Stagnation Heat Flux = " <<  bislipSystems->getCurrentHeatFluxChapman( ) << std::endl; }
            if( debugInfo == 40 ){ std::cout << "     Current Stagnation Eq. Temp. = " <<  bislipSystems->getChapmanWallTemp( ) << std::endl; }



            if( debugInfo == 1 ){ std::cout << "     computeHeatingRateTauber" << std::endl; }
            bislipSystems->setWorkingRadius( bislipSystems->getLeadingEdgeRadius() );

            if( debugInfo == 1 ){ std::cout << "     Working Radius = " << bislipSystems->getWorkingRadius() << std::endl; }
            bislipSystems->setCurrentHeatFluxTauber( bislip::Variables::computeHeatingRateTauber( bodyMap_, vehicleName_) );

            if( debugInfo == 1 ){ std::cout << "     Current Tauber Heat Flux = " <<  bislipSystems->getCurrentHeatFluxTauber( ) << std::endl; }
}

        }
   // }






    if( debugInfo == 2 ){ std::cout << "  " <<  std::endl; }
    if( debugInfo == 2 ){ std::cout << "  " <<  std::endl; }
    if( debugInfo == 2 ){ std::cout << "     Current Angle of Attack = " << tudat::unit_conversions::convertRadiansToDegrees( currentAngleOfAttack_ ) <<  std::endl; }
    if( debugInfo == 2 ){ std::cout << "     Current Bank Angle      = " << tudat::unit_conversions::convertRadiansToDegrees( currentBankAngle_ ) <<  std::endl; }
    if( debugInfo == 2 ){ std::cout << "     Current Sideslip Angle  = " << tudat::unit_conversions::convertRadiansToDegrees( currentAngleOfSideslip_ ) <<  std::endl; }
    if( debugInfo == 2 ){ std::cout << "     Current bodyflap angle  = " << tudat::unit_conversions::convertRadiansToDegrees( vehicleSystems->getCurrentControlSurfaceDeflection( "BodyFlap" ) ) << std::endl; }
    if( debugInfo == 2 ){ std::cout << "     Current elevon angle    = " << tudat::unit_conversions::convertRadiansToDegrees( vehicleSystems->getCurrentControlSurfaceDeflection( "ElevonLeft" ) ) << std::endl; }
    if( debugInfo == 2 ){ std::cout << "  " <<  std::endl; }
    if( debugInfo == 2 ){ std::cout << "  " <<  std::endl; }


    if ( std::isnan( currentAngleOfAttack_ ) == true ) { throw std::runtime_error("Angle of attack not set."); }
    if ( std::isnan( currentBankAngle_ ) == true ) { throw std::runtime_error("Bank Angle not set."); }
    if ( std::isnan( currentAngleOfSideslip_ ) == true ) { throw std::runtime_error("Sideslip Angle not set."); }



} // updateGuidance


void MyGuidance::evaluateGuidanceFunctions(
        std::shared_ptr< bislip::BislipVehicleSystems > &bislipSystems,
        std::shared_ptr< tudat::system_models::VehicleSystems > &vehicleSystems,
        const std::string &currentTrajectoryPhase )
{

    int debugInfo = bislipSystems->getDebugInfo();
    const double timeOfFLight = std::trunc( ( FlightConditions_->getCurrentTime() - bislipSystems->getStartingEpoch() ) * 100 ) / 100;

    //! Determine and set values by evaluating interpolators.
    if( debugInfo == 1 ){ std::cout << "Determine and set values by evaluating interpolators." << std::endl; }

    if ( bislipSystems->getValidationFlag( ) == true )
    {
        currentAngleOfAttack_ = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getKourouAngleOfAttackInterpolator( ) )->interpolate( timeOfFLight ) );
        bislipSystems->setCurrentAngleOfAttack( currentAngleOfAttack_ );
    }
    else
    {
        bislipSystems->setCurrentAngleOfAttack( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::AngleOfAttack, bodyMap_, vehicleName_ ) ) );
    }

    if( debugInfo == 1 ){ std::cout << "bislipSystems->getCurrentAngleOfAttack() = " << bislipSystems->getCurrentAngleOfAttack() << std::endl; }

    if ( bislipSystems->getValidationFlag( ) == true )
    {
        //currentBankAngle_ = ( bislipSystems->getKourouAngleOfAttackInterpolator( ) )->interpolate( timeOfFLight );
        bislipSystems->setEvaluatedBankAngle( tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getKourouBankAngleInterpolator( ) )->interpolate( timeOfFLight ) ) );
    }
    else
    {
        bislipSystems->setEvaluatedBankAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::BankAngle, bodyMap_, vehicleName_ ) ) );
    }


    if( debugInfo == 1 ){ std::cout << "bislipSystems->getEvaluatedBankAngle() = " << bislipSystems->getEvaluatedBankAngle() << std::endl; }

    bislipSystems->setCurrentThrustElevationAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustElevationAngle, bodyMap_, vehicleName_ ) ) );

    if( debugInfo == 1 ){ std::cout << "bislipSystems->getCurrentThrustElevationAngle() = " << bislipSystems->getCurrentThrustElevationAngle() << std::endl; }

    bislipSystems->setCurrentThrustAzimuthAngle( tudat::unit_conversions::convertDegreesToRadians( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrustAzimuthAngle, bodyMap_, vehicleName_ ) ) );
    bislipSystems->setCurrentThrottleSetting( bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::ThrottleSetting, bodyMap_, vehicleName_ ) );

    //! Determine and set values by evaluating various functions.
    if( debugInfo == 1 ){ std::cout << "Determine and set values by evaluating various functions." << std::endl; }

    bislipSystems->setCurrentEngineStatus( bislip::Variables::determineEngineStatus( bodyMap_.at( vehicleName_ )->getBodyMass(), vehicleSystems->getDryMass() ) );
    bislipSystems->setCurrentBodyFixedThrustDirection( bislip::Variables::computeBodyFixedThrustDirection( bodyMap_, vehicleName_ ) );
    bislipSystems->setCurrentThrustMagnitude( bislip::Variables::computeThrustMagnitude( bodyMap_, vehicleName_ ) );




    //************************************************** Angle of Attack
    //! Bound and assign the current angle of attack.
    if( debugInfo == 1 ){ std::cout << "Bound and assign the current angle of attack." << std::endl; }
    double angleOfAttackUpperBound = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getAlphaMachEnvelopeUBInterpolator() )->interpolate( FlightConditions_->getCurrentMachNumber() ) );
    if( debugInfo == 1 ){ std::cout << "Angle of attack Upper Bound: " << angleOfAttackUpperBound << std::endl; }
    double angleOfAttackLowerBound = tudat::unit_conversions::convertDegreesToRadians( ( bislipSystems->getAlphaMachEnvelopeLBInterpolator() )->interpolate( FlightConditions_->getCurrentMachNumber() ) );
    if( debugInfo == 1 ){ std::cout << "Angle of attack Lower Bound: " << angleOfAttackLowerBound << std::endl; }

    if ( bislipSystems->getCurrentAngleOfAttack() > angleOfAttackUpperBound )
    {
        bislipSystems->setCurrentAngleOfAttack( angleOfAttackUpperBound );
        currentAngleOfAttack_ = angleOfAttackUpperBound;
    }
    else if ( bislipSystems->getCurrentAngleOfAttack() < angleOfAttackLowerBound )
    {
        bislipSystems->setCurrentAngleOfAttack( angleOfAttackLowerBound );
        currentAngleOfAttack_ = angleOfAttackLowerBound;
    }
    else
    { currentAngleOfAttack_ = bislipSystems->getCurrentAngleOfAttack(); }

    //************************************************** Control Surface Deflections
    //! Compute and assign the current control surface deflection angles
    if( debugInfo == 1 ){ std::cout << "Compute and assign the current control surface deflections." << std::endl; }

    if ( FlightConditions_->getCurrentDynamicPressure() > bislipSystems->getMinimumDynamicPressureforControlSurface() )
    {
        Eigen::Vector2d controlSurfaceDeflections = bislip::Variables::computeControlSurfaceDeflection( bodyMap_, vehicleName_ );

        bislipSystems->setCurrentBodyFlapAngle( controlSurfaceDeflections( 0 ) );
        bislipSystems->setCurrentElevonAngle( controlSurfaceDeflections( 1 ) );
        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", bislipSystems->getCurrentBodyFlapAngle() );
        vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonLeft", bislipSystems->getCurrentElevonAngle() );
        vehicleSystems->setCurrentControlSurfaceDeflection( "ElevonRight", bislipSystems->getCurrentElevonAngle() );
    }

    //************************************************** Bank Angle
    if( debugInfo == 2 ){ std::cout << "Compute, assign, bound, and set the current bank angle." << std::endl; }

    double bankAngleForSkipSuppression = 0.0;
    bislipSystems->setCurrentBankAngleLimit( bankAngleForSkipSuppression );

    if( debugInfo == 2 )
    {
        std::cout << "     Current Trajectory Phrase: " << currentTrajectoryPhase << std::endl;
        if( currentTrajectoryPhase == "Ascent" ){ std::cout << "          No bank angle computation." << std::endl; }
    }

    double bankAngle = 0.0;

    if ( currentTrajectoryPhase == "Descent" )
    {
        if( debugInfo == 2 ){ std::cout << "          Current Dynamic Pressure = " << FlightConditions_->getCurrentDynamicPressure() << std::endl; }
        if( debugInfo == 2 ){ std::cout << "          Minimum Dynamic Pressure = " << bislipSystems->getMinimumDynamicPressureforControlSurface() << std::endl; }

        if( FlightConditions_->getCurrentDynamicPressure() > bislipSystems->getMinimumDynamicPressureforControlSurface() )
        {

         bankAngle = bislip::Variables::computeBankAngle( bodyMap_, vehicleName_ );

                if( debugInfo == 2 ){ std::cout << "          Calculated Bank Angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( bankAngle ) << std::endl; }

                if( ( bislipSystems->getValidationFlag() == false ) && ( FlightConditions_->getCurrentTime() > startingEpoch_ + bislipSystems->getSkipSuppressionTimingTrigger( ) ) )
                {
                    if( debugInfo == 2 ){ std::cout << "               Skip Suppression is Active" << std::endl; }

                    double currentFlightPathAngle = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );

                    if( debugInfo == 2 ){ std::cout << "                    Current Flight-Path Angle = " <<  currentFlightPathAngle << " rad" << std::endl; }

                    //double currentFlightPathAngleRate = bislip::Variables::computeFlightPathAngleRate( bodyMap_, vehicleName_, "Earth" );

                    if ( currentFlightPathAngle > 0 )
                    {
                        if( debugInfo == 2 ){ std::cout << "                    Calculate Bank Angle Required to Suppress Skipping" << std::endl; }

                        bankAngleForSkipSuppression = double( bislip::Variables::determineBankAngleSign( bankAngle ) ) * std::abs( bislip::Variables::computeSkipSuppressionLimit( bodyMap_, vehicleName_, centralBodyName_ ) );

                        bislipSystems->setCurrentBankAngleLimit( bankAngleForSkipSuppression );

                        if( debugInfo == 2 ){ std::cout << "                    Current Bank Angle = " <<  tudat::unit_conversions::convertRadiansToDegrees( bislipSystems->getCurrentBankAngle( ) ) << " deg  |  " << "Calculated Bank Angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( bankAngle ) << " deg  |  " << "Bank Angle Required = " <<  tudat::unit_conversions::convertRadiansToDegrees( bankAngleForSkipSuppression ) << " deg" << std::endl; }

                        //bankAngle = bankAngleForSkipSuppression;
                        if ( std::abs( bankAngleForSkipSuppression ) > std::abs( bankAngle ) ) { bankAngle = bankAngleForSkipSuppression; }
                        if( debugInfo == 2 ){ std::cout << "                    NEW Bank Angle: " <<  tudat::unit_conversions::convertRadiansToDegrees( bankAngle ) << std::endl; }

                }
            }
        }
    }
    currentBankAngle_ = bankAngle;
    bislipSystems->setCurrentBankAngle( currentBankAngle_ );



    currentAngleOfSideslip_ = 0.0;

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
