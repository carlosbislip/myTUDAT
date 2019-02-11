#include "updateGuidance.h"

namespace bislip {

void MyGuidance::updateGuidance( const double currentTime )
{
    //std::cout << "Starting Aerodynamic guidance for this evaluation" << std::endl;

    if( currentTime == currentTime )
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

        std::shared_ptr< tudat::system_models::VehicleSystems > vehicleSystems = bodyMap_.at( vehicleName_ )->getVehicleSystems();
        std::shared_ptr< bislip::VehicleSystems > bislipSystems = bodyMap_.at( vehicleName_ )->getBislipSystems();


        //! Set of parameters that I am yet to figure out how to pass/extract them around.
        //const double b_ref = 13.0;
        //const double c_ref = 23.0;
        //const double del_x = 3.0;
        //const double del_x_T = 0.0;
        //const double del_z_T = -5.0;
        //const double R_E = 6.378137e6;
        //const double mu = 3.986004418e14;
        //const double J2 = 1082.626523e-6;
        //const double J3 = 2.532153e-7;
        //const double J4 = 1.6109876e-7;
        const double n = 1.2;


        //! Extract various parameters form current flight conditions
        double current_height = FlightConditions_->getCurrentAltitude( );
        //std::cout << "current_h:  " << current_h << std::endl;
        double currentDensity = FlightConditions_->getCurrentDensity( );
        //std::cout << "currentDensity:  " << currentDensity << std::endl;
        double currentDynamicPressure = FlightConditions_->getCurrentDynamicPressure();

        double current_M = FlightConditions_->getCurrentMachNumber( );
        double currentFlightPathAngle = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::flight_path_angle );
        double currentAirspeed = FlightConditions_->getCurrentAirspeed();
        double current_AoA = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::angle_of_attack );
        double surfaceArea = bodyMap_.at( vehicleName_ )->getAerodynamicCoefficientInterface( )->getReferenceArea( );
        double currentMass = bodyMap_.at( vehicleName_ )->getBodyMass( );

        const double currentLatitude = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::latitude_angle );
        const double currentLongitude = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::longitude_angle );
        const double currentHeading = FlightConditions_->getAerodynamicAngleCalculator( )->getAerodynamicAngle( tudat::reference_frames::heading_angle );

        //const double omega = 7.292115*1E-5;//bodyMap_.at( vehicleName_ )->getCentralBodyRotationRate( );
        const Eigen::Vector6d current_state = bodyMap_.at( vehicleName_ )->getFlightConditions( )->getCurrentBodyCenteredBodyFixedState( );
        const Eigen::Vector3d r = current_state.segment( 0, 3 );
        const Eigen::Vector3d V = current_state.segment( 3, 3 );
        const double r_norm = r.norm( );
        const double V_norm = currentAirspeed;


        double angleOfAttack = tudat::unit_conversions::convertDegreesToRadians(
                    bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::AngleOfAttack, bodyMap_, vehicleName_ ) );

        currentAngleOfAttack_ = angleOfAttack ;

        double bankAngle = tudat::unit_conversions::convertDegreesToRadians(
                    bislip::Variables::evaluateGuidanceInterpolator( bislip::Parameters::Optimization::BankAngle, bodyMap_, vehicleName_ ) );

        //! Update and retrieve current aerodynamic coefficients
        //Eigen::Vector6d currentCoefficients = getCurrentCoefficients( );


        //double  bodyFlapDeflection = getBodyFlapDeflection( );
        // if ( currentFlightPathAngle > 0 ){ bodyFlap = 0; }
        // if ( currentFlightPathAngle < 0 ){ bodyFlap = 20; }

        //        double bodyFlapDeflection = tudat::unit_conversions::convertDegreesToRadians( bodyFlap );

        vehicleSystems->setCurrentControlSurfaceDeflection( "BodyFlap", getBodyFlapDeflection( ) );


        //double maxBankAngle = tudat::unit_conversions::convertDegreesToRadians(
        //          bislip::Variables::computeEquilibriumGlideLimit( bodyMap_, vehicleName_, centralBodyName_ ) );

        //if ( std::isnan( maxBankAngle ) == false ){ if ( bankAngle > maxBankAngle ){ bankAngle = maxBankAngle; } }


        //! This if statement and following expression ensures that the new bank angle has the same sign as in previous timestep.
        int sign = 1;
        if ( currentBankAngle_ < 0.0 ){ sign = -1; }
        currentBankAngle_ = double(sign) * bankAngle;

        const double angularDistanceToGo_deg = tudat::unit_conversions::convertRadiansToDegrees(
                    bislip::Variables::computeAngularDistance (
                        currentLatitude,
                        currentLongitude,
                        bodyMap_.at( vehicleName_ )->getBislipSystems()->getTargetLat(),
                        bodyMap_.at( vehicleName_ )->getBislipSystems()->getTargetLon() ) );

        const double chi_err_deg = tudat::unit_conversions::convertRadiansToDegrees(
                    bislip::Variables::computeHeadingError (
                        currentLatitude,
                        currentLongitude,
                        bodyMap_.at( vehicleName_ )->getBislipSystems()->getTargetLat(),
                        bodyMap_.at( vehicleName_ )->getBislipSystems()->getTargetLon(),
                        currentHeading ) );

        // double chi_err_deg = tudat::unit_conversions::convertRadiansToDegrees( chi_err );
        const double abs_chi_err_deg = std::abs( chi_err_deg );

        //! Declare and determine bank bankReversal conditional. If result is lower than
        //! zero, bank bankReversal could happen. Should be less than zero if the vehicle is
        //! in fact aiming away from the target. (is it really?)
        const double reversal_conditional = chi_err_deg * sign;

        //! Determine if bank reversal should occur. Function of distance from target, absolute heading error, and heading error sign value.
        bool bankReversal = getReversal( angularDistanceToGo_deg, abs_chi_err_deg, reversal_conditional );
        if ( bankReversal == true ) { currentBankAngle_ = -currentBankAngle_; }

        currentAngleOfSideslip_ = 0;

    }

} // MyGuidance
} // namespace bislip
