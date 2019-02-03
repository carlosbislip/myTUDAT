#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>


//! Trying to implement a simple aerodynamic guidance. Initially taken from the
//! TUDAT website.
//! http://tudat.tudelft.nl/tutorials/tudatFeatures/accelerationSetup/aerodynamicGuidance.html#FlightConditions
namespace tudat
{
namespace aerodynamics
{
class ValidationAerodynamicGuidance: public AerodynamicGuidance
{
public:

    ValidationAerodynamicGuidance(
            const tudat::simulation_setup::NamedBodyMap& bodyMap,
            const std::string vehicleName)
//            const double chi_err_rad )
    {
        FlightConditions_ = std::dynamic_pointer_cast< AtmosphericFlightConditions >(
                    bodyMap.at( vehicleName )->getFlightConditions( ) );
        lat_f_rad_ = bodyMap.at( vehicleName )->getTargetLat( ) ;
        lon_f_rad_ = bodyMap.at( vehicleName )->getTargetLon( ) ;
        startingEpoch_ = bodyMap.at( vehicleName )->getStartingEpoch( ) ;
//     const double    lon_f_rad_ = lon_f_rad;
    }

    void updateGuidance(const double currentTime);

private:

    std::shared_ptr< AtmosphericFlightConditions > FlightConditions_;
    double lat_f_rad_;
    double lon_f_rad_;
    double startingEpoch_;


};

} // namespace aerodynamics

} // namespace tudat
