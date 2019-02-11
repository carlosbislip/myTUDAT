#include "updateGuidance.h"

namespace bislip {

Eigen::Vector6d MyGuidance::getCurrentCoefficients( )
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
        bislipSystems_ = std::dynamic_pointer_cast< bislip::VehicleSystems >(
                    bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
    }

    return bislip::Variables::computeCurrentCoefficients( bodyMap_, vehicleName_ );
}

Eigen::Vector3d MyGuidance::getCurrentBodyFixedThrustDirection( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    return bislip::Variables::computeBodyFixedThrustDirection( bodyMap_, vehicleName_ );
}


bool MyGuidance::getCurrentEngineStatus( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    return bislip::Variables::determineEngineStatus( bodyMap_.at( vehicleName_ )->getBodyMass(), bodyMap_.at( vehicleName_ )->getVehicleSystems()->getDryMass() );
}

double MyGuidance::getCurrentThrustMagnitude( )
{
    if( FlightConditions_ == nullptr )
    {
        FlightConditions_ = std::dynamic_pointer_cast< tudat::aerodynamics::AtmosphericFlightConditions >(
                    bodyMap_.at( vehicleName_ )->getFlightConditions( ) );
    }

    return bislip::Variables::computeThrustMagnitude(  bodyMap_, vehicleName_ );
}

/*
double MyGuidance::getBodyflapCmIncrement( )
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
        bislipSystems_ = std::dynamic_pointer_cast< bislip::VehicleSystems >(
                    bodyMap_.at( vehicleName_ )->getBislipSystems( ) );
    }

    return bislip::Variables::computeBodyFlapCmIncrement( bodyMap_, vehicleName_ );
}
*/
double MyGuidance::getBodyFlapDeflection( )
{

    return bislip::Variables::computeBodyFlapDeflection( bodyMap_, vehicleName_ );
}

bool MyGuidance::getReversal( const double &angularDistanceToGo_deg, const double &abs_chi_err_deg, const double &reversal_conditional )
{

    bool bankReversal = false;

    //! Determine SIGN of bank angle. This overly complicated segment figures out
    //! the current conditions and if the bank angle sign indeed needs to change.
    //! It considers the distance from the target, the absolute heading error, and
    //! the bankReversalconditional.
    if ( ( angularDistanceToGo_deg >= 30 ) && ( abs_chi_err_deg > 15 ) && ( reversal_conditional < 0 ) )
    {
        bankReversal = true;
    }
    else if ( ( angularDistanceToGo_deg >= 30 ) && ( abs_chi_err_deg > 15 ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg >= 30 ) && ( abs_chi_err_deg < 15 ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg >= 30 ) && ( abs_chi_err_deg < 15 ) && ( reversal_conditional < 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 30 ) && ( angularDistanceToGo_deg >= 10 ) && ( abs_chi_err_deg > abs( 23 + (angularDistanceToGo_deg - 10) * ( 15 - 23 ) / ( 30 - 10 )) ) && ( reversal_conditional < 0 ) )
    {
        bankReversal = true;
    }
    else if ( ( angularDistanceToGo_deg < 30 ) && ( angularDistanceToGo_deg >= 10 ) && ( abs_chi_err_deg > abs( 23 + (angularDistanceToGo_deg - 10) * ( 15 - 23 ) / ( 30 - 10 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 30 ) && ( angularDistanceToGo_deg >= 10 ) && ( abs_chi_err_deg < abs( 23 + (angularDistanceToGo_deg - 10) * ( 15 - 23 ) / ( 30 - 10 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 10 ) && ( angularDistanceToGo_deg >= 0.94331 ) && ( abs_chi_err_deg > abs( 23 + (angularDistanceToGo_deg - 0.94331) * ( 23 - 23 ) / ( 10 - 0.94331 )) ) && ( reversal_conditional < 0 ) )
    {
        bankReversal = true;
    }
    else if ( ( angularDistanceToGo_deg < 10 ) && ( angularDistanceToGo_deg >= 0.94331 ) && ( abs_chi_err_deg > abs( 23 + (angularDistanceToGo_deg - 0.94331) * ( 23 - 23 ) / ( 10 - 0.94331 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 10 ) && ( angularDistanceToGo_deg >= 0.94331 ) && ( abs_chi_err_deg < abs( 23 + (angularDistanceToGo_deg - 0.94331) * ( 23 - 23 ) / ( 10 - 0.94331 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 0.94331 ) && ( angularDistanceToGo_deg >= 0.91452 ) && ( abs_chi_err_deg > abs( 20.36779 + (angularDistanceToGo_deg - 0.91452) * ( 23 - 20.36779 ) / ( 0.94331 - 0.91452 )) ) && ( reversal_conditional < 0 ) )
    {
        bankReversal = true;
    }
    else if ( ( angularDistanceToGo_deg < 0.94331 ) && ( angularDistanceToGo_deg >= 0.91452 ) && ( abs_chi_err_deg > abs( 20.36779 + (angularDistanceToGo_deg - 0.91452) * ( 23 - 20.36779 ) / ( 0.94331 - 0.91452 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 0.94331 ) && ( angularDistanceToGo_deg >= 0.91452 ) && ( abs_chi_err_deg < abs( 20.36779 + (angularDistanceToGo_deg - 0.91452) * ( 23 - 20.36779 ) / ( 0.94331 - 0.91452 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 0.91452 ) && ( angularDistanceToGo_deg >= 0.90739 ) && ( abs_chi_err_deg > abs( 19.62079 + (angularDistanceToGo_deg - 0.90739) * ( 20.36779 - 19.62079 ) / ( 0.91452 - 0.90739 )) ) && ( reversal_conditional < 0 ) )
    {
        bankReversal = true;
    }
    else if ( ( angularDistanceToGo_deg < 0.91452 ) && ( angularDistanceToGo_deg >= 0.90739 ) && ( abs_chi_err_deg > abs( 19.62079 + (angularDistanceToGo_deg - 0.90739) * ( 20.36779 - 19.62079 ) / ( 0.91452 - 0.90739 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 0.91452 ) && ( angularDistanceToGo_deg >= 0.90739 ) && ( abs_chi_err_deg < abs( 19.62079 + (angularDistanceToGo_deg - 0.90739) * ( 20.36779 - 19.62079 ) / ( 0.91452 - 0.90739 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 0.90739 ) && ( angularDistanceToGo_deg >= 0.86476 ) && ( abs_chi_err_deg > abs( 15.20725 + (angularDistanceToGo_deg - 0.86476) * ( 19.62079 - 15.20725 ) / ( 0.90739 - 0.86476 )) ) && ( reversal_conditional < 0 ) )
    {
        bankReversal = true;
    }
    else if ( ( angularDistanceToGo_deg < 0.90739 ) && ( angularDistanceToGo_deg >= 0.86476 ) && ( abs_chi_err_deg > abs( 15.20725 + (angularDistanceToGo_deg - 0.86476) * ( 19.62079 - 15.20725 ) / ( 0.90739 - 0.86476 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 0.90739 ) && ( angularDistanceToGo_deg >= 0.86476 ) && ( abs_chi_err_deg < abs( 15.20725 + (angularDistanceToGo_deg - 0.86476) * ( 19.62079 - 15.20725 ) / ( 0.90739 - 0.86476 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 0.86476 ) && ( angularDistanceToGo_deg >= 0.85776 ) && ( abs_chi_err_deg > abs( 14.35438 + (angularDistanceToGo_deg - 0.85776) * ( 15.20725 - 14.35438 ) / ( 0.86476 - 0.85776 )) ) && ( reversal_conditional < 0 ) )
    {
        bankReversal = true;
    }
    else if ( ( angularDistanceToGo_deg < 0.86476 ) && ( angularDistanceToGo_deg >= 0.85776 ) && ( abs_chi_err_deg > abs( 14.35438 + (angularDistanceToGo_deg - 0.85776) * ( 15.20725 - 14.35438 ) / ( 0.86476 - 0.85776 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 0.86476 ) && ( angularDistanceToGo_deg >= 0.85776 ) && ( abs_chi_err_deg < abs( 14.35438 + (angularDistanceToGo_deg - 0.85776) * ( 15.20725 - 14.35438 ) / ( 0.86476 - 0.85776 )) ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 0.85776 ) && ( abs_chi_err_deg > 14.35438 ) && ( reversal_conditional < 0 ) )
    {
        bankReversal = true;
    }
    else if ( ( angularDistanceToGo_deg < 0.85776 ) && ( abs_chi_err_deg > 14.35438 ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    else if ( ( angularDistanceToGo_deg < 0.85776 ) && ( abs_chi_err_deg < 14.35438 ) && ( reversal_conditional > 0 ) )
    {
        bankReversal = false;
    }
    return bankReversal;
}


} // namespace bislip
