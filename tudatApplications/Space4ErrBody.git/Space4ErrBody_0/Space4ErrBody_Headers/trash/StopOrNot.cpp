#include "Space4ErrBodyProblem.h"
#include "StopOrNot.h"
#include "getStuff.h"
#include "bislipVariables.h"
//##include "updateGuidance.h"

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>

namespace bislip {

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
    const double lat_i_rad = bodyMap.at( vehicleName )->getVehicleSystems()->getInitialCoordinates().first;
    const double lon_i_rad = bodyMap.at( vehicleName )->getVehicleSystems()->getInitialCoordinates().second;

    //! Extract Final distance to target
    const double final_d_to_target_rad = bodyMap.at( vehicleName )->getVehicleSystems()->getFinalDistanceToTarget();
    const double h_UP = terminationConditionsValues[ 3 ];
    const double h_DN = terminationConditionsValues[ 4 ];
    const double V_max = terminationConditionsValues[ 5 ];
    const double n_max = terminationConditionsValues[ 6 ];
    const double q_dot_max = terminationConditionsValues[ 7 ];
    const double q_dyn_max = terminationConditionsValues[ 8 ];

    //! Extract and convert target coordinates
    const double lat_f_rad = tudat::unit_conversions::convertDegreesToRadians( bodyMap.at( vehicleName )->getVehicleSystems()->getTargetCoordinates().first );
    const double lon_f_rad = tudat::unit_conversions::convertDegreesToRadians( bodyMap.at( vehicleName )->getVehicleSystems()->getTargetCoordinates().second );

    //! Calculate current Distance to target.
    const double current_d_to_target_rad = bislip::variables::computeAngularDistance( lat_c_rad , lon_c_rad , lat_f_rad , lon_f_rad );

    //! Calculate total distance traveled
    const double total_d_traveled_rad = bislip::variables::computeAngularDistance( lat_i_rad , lon_i_rad , lat_c_rad , lon_c_rad );

    //! Extract Initial distance to target
    const double initial_d_to_target_rad = bodyMap.at( vehicleName )->getVehicleSystems()->getInitialDistanceToTarget();

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
    const double E_hat = bislip::variables::computeNormalizedSpecificEnergy( current_height, current_V, bodyMap.at( vehicleName )->getVehicleSystems()->getE_max() );

    //! Evaluate current throttle setting and thrust elevation angle.
    double throttle = bislip::variables::evaluateGuidanceInterpolator(
                    //current_gamma,
                    "Throttle Setting",
                    bodyMap.at( vehicleName )->getVehicleSystems(),
                    current_height,
                    current_V,
                    bodyMap.at( vehicleName )->getVehicleSystems()->getE_max() );

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
    }*/
    if ( ( ( throttle >= 0 ) && ( throttle <= 1.0 ) )
              &&  ( ( current_gamma > tudat::mathematical_constants::PI / 2.0 ) ) )
    {
        //! This condition terminates the simulation if flight-path angle is negative or beyond 90 deg while the engine is on.
        done = true;
        std::cout << " sixth " << std::endl;
        std::cout << "throttle:  " << throttle << std::endl;
        std::cout << "current_gamma:  " << tudat::unit_conversions::convertRadiansToDegrees( current_gamma ) << std::endl;
        std::cout << "Stopping propagation" << std::endl;
    }
   /* else if ( E_hat > 1.0 )
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

    //}*/
    else
    {
        done = false;
    }




    return done;
}
}
