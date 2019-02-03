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
        const std::vector< std::vector< double > > &bounds,
        const std::string &problem_name,
        const std::string &vehicleName,
        const std::vector< std::string > &parameterList_Ascent,
        const std::vector< double > &parameterBounds_Ascent,
        const std::vector< std::string > &parameterList_Descent,
        const std::vector< double > &parameterBounds_Descent,
        const std::vector< double > &vehicleParameterValues,
        const std::vector< std::string > &aeroCoeffFileList,
        const std::vector< double > &simulation_settingsValues,
        const std::vector< double > &initialConditionsValues,
        const std::vector< double > &terminationConditionsValues,
        const std::vector< double > &output_settingsValues,
        const std::string &outputSubFolder,
        const Eigen::Vector6d &initialState_spherical,
        const std::vector< std::string > &centralBodies,
        const std::vector< std::string > &bodiesToIntegrate,
        const std::map< bislip::Parameters::Optimization, std::pair < double, double > > &Bounds_Ascent,
        const std::map< bislip::Parameters::Optimization, std::pair < double, double > > &Bounds_Descent,
        const tudat::simulation_setup::NamedBodyMap& bodyMap,
        const tudat::basic_astrodynamics::AccelerationMap &accelerationsMap,
        const std::shared_ptr< tudat::ephemerides::RotationalEphemeris > &earthRotationalEphemeris,
        const std::shared_ptr< tudat::propagators::DependentVariableSaveSettings > &dependentVariablesToSave,
        const std::shared_ptr< tudat::propagators::PropagationTerminationSettings > &terminationSettings_Ascent,
        const std::shared_ptr< tudat::propagators::PropagationTerminationSettings > &terminationSettings_Descent ):
    problemBounds_( bounds ),
    problem_name_( problem_name ),
    vehicleName_( vehicleName ),
    parameterList_Ascent_( parameterList_Ascent ),
    parameterBounds_Ascent_( parameterBounds_Ascent ),
    parameterList_Descent_( parameterList_Descent ),
    parameterBounds_Descent_( parameterBounds_Descent ),
    vehicleParameterValues_( vehicleParameterValues ),
    aeroCoeffFileList_( aeroCoeffFileList ),
    simulation_settingsValues_( simulation_settingsValues ),
    initialConditionsValues_( initialConditionsValues ),
    terminationConditionsValues_( terminationConditionsValues ),
    output_settingsValues_( output_settingsValues ),
    outputSubFolder_( outputSubFolder ),
    initialState_spherical_( initialState_spherical) ,
   // bodiesWithMassToPropagate_( bodiesWithMassToPropagate ),
    centralBodies_( centralBodies ),
    bodiesToIntegrate_( bodiesToIntegrate ),
    Bounds_Ascent_( Bounds_Ascent ),
    Bounds_Descent_( Bounds_Descent ),
    bodyMap_( bodyMap ),
    accelerationsMap_( accelerationsMap ),
    earthRotationalEphemeris_( earthRotationalEphemeris ),
    dependentVariablesToSave_( dependentVariablesToSave ),
    terminationSettings_Ascent_( terminationSettings_Ascent ),
    terminationSettings_Descent_( terminationSettings_Descent )
    //massRateModels_( massRateModels ),
    //massPropagatorSettings_Ascent_( massPropagatorSettings_Ascent )
{ }

//! Descriptive name of the problem
std::string Space4ErrBodyProblem::get_name() const
{
    return problem_name_;
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > Space4ErrBodyProblem::get_bounds() const
{
    return { problemBounds_[ 0 ], problemBounds_[ 1 ] };
}

std::vector<double> Space4ErrBodyProblem::fitness( const std::vector< double > &x )  const
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
    using namespace tudat_applications;
    using namespace tudat::aerodynamics;
    using namespace bislip;
    //using namespace std::this_thread; // sleep_for, sleep_until
    //using namespace std::chrono; // nanoseconds, system_clock, seconds

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            UNPACK INPUT DATA             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //std::cout << "Unpacking data" << std::endl;

    //! Declare and initialize simulation start epoch.
    const double simulationStartEpoch = simulation_settingsValues_[ 0 ]; // 10/28/2018  11:00:00 AM  ---> 2458419.95833333000000000000000

    //! Declare and initialize simulation end epoch.
    const double simulationEndEpoch = simulationStartEpoch + simulation_settingsValues_[ 1 ];

    //! Declare and initialize numerical integration fixed step size.
    const double fixedStepSize = simulation_settingsValues_[ 2 ];

    //! Declare and initialize number of control nodes.
    const unsigned long nodesAscent = simulation_settingsValues_[ 3 ];
    const unsigned long nodesDescent = simulation_settingsValues_[ 4 ];

    //! Declare and initialize position vector of moment reference center
    const Eigen::Vector3d R_mrc( vehicleParameterValues_[ 3 ], vehicleParameterValues_[ 4 ], vehicleParameterValues_[ 5 ] ); // m

    //! Declare and initialize position vector of center of mass
    const Eigen::Vector3d R_com( vehicleParameterValues_[ 6 ], vehicleParameterValues_[ 7 ], vehicleParameterValues_[ 8 ] ); // m

    //! Declare and initialize position vector of center of thrust
    const Eigen::Vector3d R_cot( vehicleParameterValues_[ 9 ], vehicleParameterValues_[ 10 ], vehicleParameterValues_[ 11 ] ); // m

    //! Declare and initialize initial mass
    const double initialMass_Ascent = vehicleParameterValues_[ 12 ]; // kg

    //! Declare and initialize starting height
    const double h_i = initialConditionsValues_[ 2 ]; // m

    //! Declare and initialize starting position coordinates.
    const double initialLat_deg = initialConditionsValues_[ 0 ];
    const double initialLon_deg = initialConditionsValues_[ 1 ];

    //! Declare and initialize final position coordinates and additional termination conditions
    const double targetLat_deg = terminationConditionsValues_[ 0 ];
    const double targetLon_deg = terminationConditionsValues_[ 1 ];

    //! Declare and initialize various termination conditions
    const double finalDistanceToTarget_deg = terminationConditionsValues_[ 2 ];
    const double h_DN = terminationConditionsValues_[ 4 ];
    const double constraint_MechanicalLoad = terminationConditionsValues_[ 6 ];
    const double constraint_HeatingRate = terminationConditionsValues_[ 7 ];
    const double constraint_DynamicPressure = terminationConditionsValues_[ 8 ];


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


    const std::string centralBodyName = "Earth";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE NODAL STRUCTURE             /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Creating nodal structure" << std::endl;

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

    //! Various parameters
    const double R_E = spice_interface::getAverageRadius( centralBodyName );
    const double omega_E = 7.292115*1E-5;
    const double J2 = 1082.626523e-6;
    const double J3 = 2.532153e-7;
    const double J4 = 1.6109876e-7;

    //Eigen::Vector2d gravs = getGravs ( mu, J2, J3, J4, R_E, R_E, 0.0000001 );
    double g0 = tudat::physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION;

    //  std::cout << "Re-allocate Ascent DVs" << std::endl;
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
    const unsigned long  N = ( parameterList_Ascent_.size() - 4 ) * nodesAscent - 1;

    //! Declare and initialize various parameters common to the entire trajectory.
    //!     Initial velocity.
    const double initialAirspeed = x[ N ];

    //!     Velocity at the maximum height.
    const double upperLimitAirspeed = x[ N + 1 ];

    //!     Maximum height.
    const double h_UP = x[ N + 2 ];

    //!     Additional mass.
    const double additionalMass = x[ N + 3 ];

    const double goaldelV = upperLimitAirspeed - initialAirspeed;

    //! Set body Mass.
    bodyMap_.at( vehicleName_ )->setConstantBodyMass( initialMass_Ascent + additionalMass );
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setInitialMass( initialMass_Ascent + additionalMass );

    // std::cout << "Re-allocate Descent DVs" << std::endl;
    //! Re-allocate decision vector values into workable vectors for Descent phase.
    for( unsigned int i = 0; i < nodesDescent; i++ )
    {
        if ( i < ( nodesDescent - 1) ){ xn_interval_Descent( i ) = x[ ( N + 3 ) + i ]; }
        alpha_deg_Descent( i ) = x[ ( N + 3 ) + i + 1 * nodesDescent - 1 ];
        sigma_deg_Descent( i ) = x[ ( N + 3 ) + i + 2 * nodesDescent - 1 ];
        eps_T_deg_Descent( i ) = x[ ( N + 3 ) + i + 3 * nodesDescent - 1 ];
        phi_T_deg_Descent( i ) = x[ ( N + 3 ) + i + 4 * nodesDescent - 1 ];
        throttle_Descent( i )  = x[ ( N + 3 ) + i + 5 * nodesDescent - 1 ];
    }

    //! Declare and initialize number of parameters exclusive to Descent phase.
    const unsigned long NN = ( parameterList_Descent_.size() - 1 ) * nodesDescent - 1;

    //! Declare and initialize last parameter common to the entire trajectory.
    //!     Final velocity.
    const double lowerLimitAirspeed = x[ N + NN + 3 ];

    //std::cout << "Create vector of node locations for ascent" << std::endl;
    //! Create vector of node locations for Ascent phase.
    xn_Ascent( 0 ) = 0;
    for( unsigned int i = 0; i < nodesAscent - 1; i++ )
    {
        xn_Ascent( i + 1 )        = xn_Ascent( i ) + xn_interval_Ascent( i ) / xn_interval_Ascent.sum();
    }

    // std::cout << "Create vector of node locations for descent" << std::endl;
    //! Create vector of node locations for Descent phase.
    xn_Descent( 0 ) = 0;//xn_Ascent( xn_Ascent.size() - 1 ) ;//xn_interval_Descent.sum();
    for( unsigned int i = 0; i < nodesDescent - 1; i++ )
    {
        xn_Descent( i + 1 )        = xn_Descent( i ) + xn_interval_Descent( i ) / xn_interval_Descent.sum();
    }


    // xn_Descent = xn_Descent.reverse().eval();

    /*
 *
 *
    std::cout << "x =  " << std::endl;
    for( int i = 0; i < int( x.size( ) ); i++ )
    {
        std::cout << i << "   " << x[i] << std::endl;
    }
    std::cout << "-------" << std::endl;
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
    std::cout << "initialAirspeed =  " << initialAirspeed << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "upperLimitAirspeed =  " << upperLimitAirspeed << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "lowerLimitAirspeed =  " << lowerLimitAirspeed << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << "h_UP =  " << h_UP << std::endl;
    std::cout << "-------" << std::endl;
*/
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            COMPLETE KNOWN STATES            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::Vector6d initialStateAscent = initialState_spherical_;
    initialStateAscent( 3 ) = initialAirspeed;

    //! Two things are being done here
    //!     Converting state vector from spherical to Cartesian elements
    //!     Transforming state vector from Earth-Fixed frame to Inertial frame.
    const Eigen::Vector6d systemInitialState_Ascent = tudat::ephemerides::transformStateToGlobalFrame(
                tudat::orbital_element_conversions::convertSphericalOrbitalToCartesianState( initialStateAscent ),
                simulationStartEpoch,
                earthRotationalEphemeris_ );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE INTERPOLATORS             ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Creating interpolators" << std::endl;

    //! Calculate initial, maximum, and final specific energy levels.
    const double E_i = bislip::Variables::computeSpecificEnergy( h_i, initialAirspeed );
    const double E_max = bislip::Variables::computeSpecificEnergy( h_UP, upperLimitAirspeed );
    const double E_f = bislip::Variables::computeSpecificEnergy( h_DN, lowerLimitAirspeed );

    //! Set maximum Energy level.
    //vehicleSystems->setE_max( E_max );
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setE_max( E_max );

    //! Normalize initial, maximum, and final specific energy levels.
    //const double E_hat_i = bislip::Variables::computeNormalizedSpecificEnergy( h_i, initialAirspeed, E_max );
    const double E_hat_max = bislip::Variables::computeNormalizedSpecificEnergy( h_UP, upperLimitAirspeed, E_max );
    //const double E_hat_f = bislip::Variables::computeNormalizedSpecificEnergy( h_DN, lowerLimitAirspeed, E_max );

    //! Map normalized specific energy levels to control node locations.
    Eigen::VectorXd E_mapped_Ascent( xn_Ascent.size() );
    Eigen::VectorXd E_mapped_Descent( xn_Descent.size() );
    E_mapped_Ascent = ( ( E_max - E_i ) * xn_Ascent.array() + E_i ) / E_max;
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
    std::map< double, Eigen::VectorXd > map_DV_mapped_Ascent;
    Eigen::VectorXd DV_mapped_Ascent ( 6 );

    // std::cout << "Mapping Ascent DVs" << std::endl;
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
    }

    //! Declare data maps used for decision vector values related to Descent phase.
    std::map< double, double > map_alpha_deg_Descent, map_eps_T_deg_Descent, map_phi_T_deg_Descent, map_throttle_Descent, map_sigma_deg_Descent;
    std::map< double, Eigen::VectorXd > map_DV_mapped_Descent;
    Eigen::VectorXd DV_mapped_Descent ( 6 );

    // std::cout << "Mapping Descent DVs" << std::endl;
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
    }

    //! Calculate first derivatives for Hermite Spline Interpolation according to https://doi.org/10.1016/j.cam.2017.09.049.
    //std::vector< double > alpha_deg_Ascent_derivatives = bislip::Variables::HermiteDerivatives( E_mapped_Ascent, alpha_deg_Ascent, E_mapped_Ascent.size() );
    //std::vector< double > sigma_deg_Ascent_derivatives = bislip::Variables::HermiteDerivatives( E_mapped_Ascent, sigma_deg_Ascent );
    //   std::vector< double > eps_T_deg_Ascent_derivatives = bislip::Variables::HermiteDerivatives( E_mapped_Ascent, eps_T_deg_Ascent );
    //  std::vector< double > phi_T_deg_Ascent_derivatives = bislip::Variables::HermiteDerivatives( E_mapped_Ascent, phi_T_deg_Ascent );
    // std::vector< double > throttle_Ascent_derivatives = bislip::Variables::HermiteDerivatives( E_mapped_Ascent, throttle_Ascent );
    //   std::vector< double > alpha_deg_Descent_derivatives = bislip::Variables::HermiteDerivatives( E_mapped_Descent, alpha_deg_Descent );
    //  std::vector< double > sigma_deg_Descent_derivatives = bislip::Variables::HermiteDerivatives( E_mapped_Ascent, sigma_deg_Descent );
    // std::vector< double > eps_T_deg_Descent_derivatives = bislip::Variables::HermiteDerivatives( E_mapped_Ascent, eps_T_deg_Descent );
    // std::vector< double > phi_T_deg_Descent_derivatives = bislip::Variables::HermiteDerivatives( E_mapped_Ascent, phi_T_deg_Descent );
    // std::vector< double > throttle_Descent_derivatives = bislip::Variables::HermiteDerivatives( E_mapped_Ascent, throttle_Descent );

    //! Declare and initialize interpolator settings.
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = std::make_shared< interpolators::InterpolatorSettings >( hermite_spline_interpolator );

    //std::cout << "Creating Ascent Interpolators" << std::endl;
    //! Declare and initialize interpolators for Ascent phase.
    std::shared_ptr< OneDimensionalInterpolator< double, double > > interpolator_alpha_deg_Ascent, interpolator_sigma_deg_Ascent, interpolator_eps_T_deg_Ascent, interpolator_phi_T_deg_Ascent, interpolator_throttle_Ascent;
    interpolator_alpha_deg_Ascent = bislip::Variables::createOneDimensionalHermiteInterpolator( alpha_deg_Ascent, E_mapped_Ascent, map_alpha_deg_Ascent, interpolatorSettings );
    interpolator_sigma_deg_Ascent = bislip::Variables::createOneDimensionalHermiteInterpolator( sigma_deg_Ascent, E_mapped_Ascent, map_sigma_deg_Ascent, interpolatorSettings );
    interpolator_eps_T_deg_Ascent = bislip::Variables::createOneDimensionalHermiteInterpolator( eps_T_deg_Ascent, E_mapped_Ascent, map_eps_T_deg_Ascent, interpolatorSettings );
    interpolator_phi_T_deg_Ascent = bislip::Variables::createOneDimensionalHermiteInterpolator( phi_T_deg_Ascent, E_mapped_Ascent, map_phi_T_deg_Ascent, interpolatorSettings );
    interpolator_throttle_Ascent = bislip::Variables::createOneDimensionalHermiteInterpolator( throttle_Ascent, E_mapped_Ascent, map_throttle_Ascent, interpolatorSettings );


    //interpolator_alpha_deg_Ascent = interpolators::createOneDimensionalInterpolator< double, double >( map_alpha_deg_Ascent, interpolatorSettings, std::make_pair( alpha_deg_Ascent( 0 ), alpha_deg_Ascent( E_mapped_Ascent.size() - 1 ) ), alpha_deg_Ascent_derivatives );
    //interpolator_sigma_deg_Ascent = interpolators::createOneDimensionalInterpolator< double, double >( map_sigma_deg_Ascent, interpolatorSettings, std::make_pair( sigma_deg_Ascent( 0 ), sigma_deg_Ascent( E_mapped_Ascent.size() - 1 ) ), sigma_deg_Ascent_derivatives );
    //interpolator_eps_T_deg_Ascent = interpolators::createOneDimensionalInterpolator< double, double >( map_eps_T_deg_Ascent, interpolatorSettings, std::make_pair( eps_T_deg_Ascent( 0 ), eps_T_deg_Ascent( E_mapped_Ascent.size() - 1 ) ), eps_T_deg_Ascent_derivatives );
    //interpolator_phi_T_deg_Ascent = interpolators::createOneDimensionalInterpolator< double, double >( map_phi_T_deg_Ascent, interpolatorSettings, std::make_pair( phi_T_deg_Ascent( 0 ), phi_T_deg_Ascent( E_mapped_Ascent.size() - 1 ) ), phi_T_deg_Ascent_derivatives );
    //interpolator_throttle_Ascent = interpolators::createOneDimensionalInterpolator< double, double >( map_throttle_Ascent, interpolatorSettings, std::make_pair( throttle_Ascent( 0 ), throttle_Ascent( E_mapped_Ascent.size() - 1 ) ), throttle_Ascent_derivatives );

    //std::cout << "Creating Descent Interpolators" << std::endl;
    //! Declare and initialize interpolators for Descent phase.
    std::shared_ptr< OneDimensionalInterpolator< double, double > > interpolator_alpha_deg_Descent, interpolator_sigma_deg_Descent, interpolator_eps_T_deg_Descent, interpolator_phi_T_deg_Descent, interpolator_throttle_Descent;
    interpolator_alpha_deg_Descent = bislip::Variables::createOneDimensionalHermiteInterpolator( alpha_deg_Descent, E_mapped_Descent, map_alpha_deg_Descent, interpolatorSettings );
    interpolator_sigma_deg_Descent = bislip::Variables::createOneDimensionalHermiteInterpolator( sigma_deg_Descent, E_mapped_Descent, map_sigma_deg_Descent, interpolatorSettings );
    interpolator_eps_T_deg_Descent = bislip::Variables::createOneDimensionalHermiteInterpolator( eps_T_deg_Descent, E_mapped_Descent, map_eps_T_deg_Descent, interpolatorSettings );
    interpolator_phi_T_deg_Descent = bislip::Variables::createOneDimensionalHermiteInterpolator( phi_T_deg_Descent, E_mapped_Descent, map_phi_T_deg_Descent, interpolatorSettings );
    interpolator_throttle_Descent = bislip::Variables::createOneDimensionalHermiteInterpolator( throttle_Descent, E_mapped_Descent, map_throttle_Descent, interpolatorSettings );


    //! Declare vectors containing interpolated values.
    Eigen::VectorXd interpolated_values_Ascent( 5 ), interpolated_values_Descent( 5 );

    //! Declare data map to contain vectors of interpolated values.
    std::map< double, Eigen::VectorXd > interpolators_Ascent, interpolators_Descent;

    //! Declare evaluation variable.
    double eval;

    //std::cout << "Evaluate Interpolators and store to print out" << std::endl;
    //! Loop to populate vectors of interpolated values and then pass to data map.
    //!     Number of evaluations has been arbitrarily selected.
    double pp = 0;
    for ( unsigned int i = 0; i < 1001; ++i )
    {
        eval = pp * E_mapped_Ascent( E_mapped_Ascent.size() - 1 ) / 1000;
        interpolated_values_Ascent( 0 ) = interpolator_alpha_deg_Ascent->interpolate( eval );
        interpolated_values_Ascent( 1 ) = interpolator_sigma_deg_Ascent->interpolate( eval );
        interpolated_values_Ascent( 2 ) = interpolator_eps_T_deg_Ascent->interpolate( eval );
        interpolated_values_Ascent( 3 ) = interpolator_phi_T_deg_Ascent->interpolate( eval );
        interpolated_values_Ascent( 4 ) = interpolator_throttle_Ascent->interpolate( eval );

        interpolators_Ascent[ eval ] = interpolated_values_Ascent;

        eval = pp * E_mapped_Descent( 0 ) / 1000;
        interpolated_values_Descent( 0 ) = interpolator_alpha_deg_Descent->interpolate( eval );
        interpolated_values_Descent( 1 ) = interpolator_sigma_deg_Descent->interpolate( eval );
        interpolated_values_Descent( 2 ) = interpolator_eps_T_deg_Descent->interpolate( eval );
        interpolated_values_Descent( 3 ) = interpolator_phi_T_deg_Descent->interpolate( eval );
        interpolated_values_Descent( 4 ) = interpolator_throttle_Descent->interpolate( eval );

        interpolators_Descent[ eval ] = interpolated_values_Descent;

        pp += 1;
    }

    //std::cout << "Print out Interpolators and DVs" << std::endl;
    //! Print data maps containing vectors of evaluation of interpolators.
    writeDataMapToTextFile( interpolators_Ascent,
                            "interpolators_pre_Ascent",
                            tudat_applications::getOutputPath( ) + outputSubFolder_,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );
    writeDataMapToTextFile( map_DV_mapped_Ascent,
                            "DV_mapped_pre_Ascent",
                            tudat_applications::getOutputPath( ) + outputSubFolder_,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );
    writeDataMapToTextFile( interpolators_Descent,
                            "interpolators_pre_Descent",
                            tudat_applications::getOutputPath( ) + outputSubFolder_,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );
    writeDataMapToTextFile( map_DV_mapped_Descent,
                            "DV_mapped_pre_Descent",
                            tudat_applications::getOutputPath( ) + outputSubFolder_,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE MASS RATE SETTINGS: ASCENT            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Create mass rate models" << std::endl;

    //! Declare and initialize mass rate model settings.
    std::shared_ptr< MassRateModelSettings > massRateModelSettings =
            std::make_shared< FromThrustMassModelSettings >( true );

    //! Declare and initialize mass rate model settings.
    std::map< std::string, std::shared_ptr< MassRateModel > > massRateModels;
    massRateModels[ vehicleName_ ] = createMassRateModel(
                vehicleName_, massRateModelSettings, bodyMap_, accelerationsMap_ );

    //! Create settings for propagating the mass of the vehicle.
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back( vehicleName_ );

    //! Declare and initialize starting mass of vehicle.
    Eigen::VectorXd initialBodyMasses_Ascent( 1 );
    initialBodyMasses_Ascent( 0 ) = initialMass_Ascent + additionalMass;

    //! Declare and initialize mass propagation settings.
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings_Ascent =
            std::make_shared< MassPropagatorSettings< double > >(
                bodiesWithMassToPropagate, massRateModels, initialBodyMasses_Ascent, terminationSettings_Ascent_, dependentVariablesToSave_ );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE PROPAGATION SETTINGS: ASCENT              ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // std::cout << "Create Propagation Settings - Ascent" << std::endl;

    //! Create translational propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings_Ascent =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies_,
              accelerationsMap_,
              bodiesToIntegrate_,
              systemInitialState_Ascent,
              terminationSettings_Ascent_,
              cowell,
              dependentVariablesToSave_ );

    //! Declare and initialize list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector_Ascent;
    propagatorSettingsVector_Ascent.push_back( translationalPropagatorSettings_Ascent );
    propagatorSettingsVector_Ascent.push_back( massPropagatorSettings_Ascent );

    //! Declare and initialize propagation settings for both mass and translational dynamics.
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings_Ascent =
            std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector_Ascent, terminationSettings_Ascent_, dependentVariablesToSave_ );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE INTEGRATION SETTINGS: ASCENT           ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Create Integration Settings: Ascent" << std::endl;

    //! Declare and initialize integrator settings.
    std::shared_ptr< IntegratorSettings<  > > integratorSettings_Ascent =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationStartEpoch,
              fixedStepSize, 1, false );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           ASSIGN INTERPOLATORS & BOUNDS: ASCENT           /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Pass Ascent Interpolators to vehicle systems" << std::endl;
    //! Pass Ascent phase interpolators to vehicle systems.
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setAngleOfAttackInterpolator( interpolator_alpha_deg_Ascent );
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setBankAngleInterpolator( interpolator_sigma_deg_Ascent );
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setThrustElevationAngleInterpolator( interpolator_eps_T_deg_Ascent );
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setThrustAzimuthAngleInterpolator( interpolator_phi_T_deg_Ascent );
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setThrottleInterpolator( interpolator_throttle_Ascent );

    //    std::cout << "Pass Ascent Bounds to vehicle systems" << std::endl;

    //! Pass parameter bounds to vehicle systems.
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setParameterBounds( Bounds_Ascent_ );


    //map<int, int> m;
    //std::vector<std::string> v;
    //for(std::unordered_map< std::string, std::shared_ptr< tudat::simulation_setup::Body > >::iterator it = bodyMap_.begin(); it != bodyMap_.end(); ++it) {
      //v.push_back(it->first);
    // std::cout << it->first << "\n";
    //}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          PROPAGATE TRAJECTORY: ASCENT         /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Starting Ascent propagation" << std::endl;

    //! Propagate trajectory.
    SingleArcDynamicsSimulator< double > ascentSimulation(
                bodyMap_,
                integratorSettings_Ascent,
                propagatorSettings_Ascent );

    //std::cout << "Ascent Propagation done" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////          RETRIEVE ASCENT-DESCENT PROPAGATION LINKING DATA          //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Retrieve final Ascent epoch
    const double simulationEndEpoch_Ascent =  ( ascentSimulation.getEquationsOfMotionNumericalSolution().rbegin() )->first;

    //! Retrieve final Ascent state, equal to initial Descent state
    const Eigen::Vector6d systemInitialState_Descent =  ( ascentSimulation.getEquationsOfMotionNumericalSolution().rbegin() )->second;

    //! Retrieve dependent variables of final Ascent state
    const Eigen::VectorXd dependentVariableFinalState_Ascent = ( ascentSimulation.getDependentVariableHistory( ).rbegin() )->second;

    //! Retrieve final Ascent mass, equal to initial Descent mass
    const double initialMass_Descent = dependentVariableFinalState_Ascent[ 27 ];

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE PROPAGATION SETTINGS: DESCENT              ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Create Propagation Settings - Descent" << std::endl;

    //! Create translational propagation settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings_Descent =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies_,
              accelerationsMap_,
              bodiesToIntegrate_,
              systemInitialState_Descent,
              terminationSettings_Descent_,
              cowell,
              dependentVariablesToSave_ );

    //! Declare and initialize starting mass of vehicle.
    Eigen::VectorXd initialBodyMasses_Descent( 1 );
    initialBodyMasses_Descent( 0 ) = initialMass_Descent;

    //! Declare and initialize mass propagation settings.
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings_Descent =
            std::make_shared< MassPropagatorSettings< double > >(
                bodiesWithMassToPropagate, massRateModels, initialBodyMasses_Descent, terminationSettings_Descent_, dependentVariablesToSave_ );

    //! Declare and initialize list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector_Descent;
    propagatorSettingsVector_Descent.push_back( translationalPropagatorSettings_Descent );
    propagatorSettingsVector_Descent.push_back( massPropagatorSettings_Descent );

    //! Declare and initialize propagation settings for both mass and translational dynamics.
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings_Descent =
            std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector_Descent, terminationSettings_Descent_, dependentVariablesToSave_ );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           CREATE INTEGRATION SETTINGS: DESCENT           //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Create Integration Settings: Descent" << std::endl;

    std::shared_ptr< IntegratorSettings<  > > integratorSettings_Descent =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4,
              simulationEndEpoch_Ascent,
              fixedStepSize, 1, false );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           ASSIGN INTERPOLATORS & BOUNDS: DESCENT           ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Pass Descent Interpolators to vehicle systems" << std::endl;
    //! Pass Descent phase interpolators to vehicle systems.
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setAngleOfAttackInterpolator( interpolator_alpha_deg_Descent );
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setBankAngleInterpolator( interpolator_sigma_deg_Descent );
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setThrustElevationAngleInterpolator( interpolator_eps_T_deg_Descent );
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setThrustAzimuthAngleInterpolator( interpolator_phi_T_deg_Descent );
    bodyMap_.at( vehicleName_ )->getBislipSystems()->setThrottleInterpolator( interpolator_throttle_Descent );

    //std::cout << "Pass Descent Bounds to vehicle systems" << std::endl;

    bodyMap_.at( vehicleName_ )->getBislipSystems()->setParameterBounds( Bounds_Descent_ );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          PROPAGATE TRAJECTORY: DESCENT         ////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Starting Descent propagation" << std::endl;

    //! Propagate trajectory.
    SingleArcDynamicsSimulator< double > descentSimulation(
                bodyMap_,
                integratorSettings_Descent,
                propagatorSettings_Descent );

    //std::cout << "Descent Propagation done" << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          RETRIEVE DATA MAPS          //////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //std::cout << "Retrieve data maps" << std::endl;

    //! Retrieve propagation and dependent variable maps - Ascent
    const std::map< double, Eigen::VectorXd > propagationTimeHistoryMap_Ascent = ascentSimulation.getEquationsOfMotionNumericalSolution( );
    const std::map< double, Eigen::VectorXd > dependentVariableTimeHistoryMap_Ascent = ascentSimulation.getDependentVariableHistory( );

    //! Retrieve propagation and dependent variable maps - Descent
    const std::map< double, Eigen::VectorXd > propagationTimeHistoryMap_Descent = descentSimulation.getEquationsOfMotionNumericalSolution( );
    const std::map< double, Eigen::VectorXd > dependentVariableTimeHistoryMap_Descent = descentSimulation.getDependentVariableHistory( );

    //std::cout << "Merge Ascent & Descent data maps" << std::endl;

    //! Merge propagation data maps
    std::map< double, Eigen::VectorXd > propagationTimeHistoryMap;
    propagationTimeHistoryMap.insert( propagationTimeHistoryMap_Ascent.begin(), propagationTimeHistoryMap_Ascent.end() );
    propagationTimeHistoryMap.insert( propagationTimeHistoryMap_Descent.begin(), propagationTimeHistoryMap_Descent.end() );

    //! Merge dependent variable data maps
    std::map< double, Eigen::VectorXd > dependentVariableTimeHistoryMap;
    dependentVariableTimeHistoryMap.insert( dependentVariableTimeHistoryMap_Ascent.begin(), dependentVariableTimeHistoryMap_Ascent.end() );
    dependentVariableTimeHistoryMap.insert( dependentVariableTimeHistoryMap_Descent.begin(), dependentVariableTimeHistoryMap_Descent.end() );

    //! Extract final epoch
    const double simulationEndEpoch_Descent =  ( propagationTimeHistoryMap_Descent.rbegin() )->first;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        CLEAN UP AND CALCULATE FITNESS                //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Initialize variables used objectives functions.
    double finalDistanceToTarget_deg_calc = initialDistanceToTarget_deg;
    //std::cout << "Cleaning up and extracting Time Histories" << std::endl;
    double h_DN_calc = h_DN;
    double h_UP_calc = h_UP;
    double upperLimitAirspeed_calc = upperLimitAirspeed;
    double lowerLimitAirspeed_calc = lowerLimitAirspeed;
    double tof = simulation_settingsValues_[ 1 ];
    double maximum_NormalizedSpecificEnergy = E_hat_max;
    double maximumAirspeed = initialAirspeed;
    //double maximum_DynamicPressure = constraint_DynamicPressure;
    //double maximum_MechanicalLoad = constraint_MechanicalLoad;
    double finalMass_Descent = 0;
    double theoreticaldelV = g0 * bodyMap_.at( vehicleName_ )->getBislipSystems()->getSpecificImpulse() * log( bodyMap_.at( vehicleName_ )->getBislipSystems()->getInitialMass() / bodyMap_.at( vehicleName_ )->getVehicleSystems()->getDryMass() );

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
    tof = simulationEndEpoch_Descent - simulationStartEpoch;

    //! Transform final state from Inertial Frame to Earth-Fixed Frame
    //Eigen::Vector6d systemFinalState_EARTh_UPIXED = transformStateToTargetFrame( systemFinalState, simulationEndEpoch_calc, earthRotationalEphemeris );

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
    // const std::map< double, Eigen::VectorXd > dependentVariableTimeHistoryMap_Ascent = dynamicsSimulator.getDependentVariableHistory( );
    // std::cout << "dependentVariable_FINALSTATE: " << dependentVariable_FINALSTATE << std::endl;

    //! Declare number of elements in time history.
    long rowsAscent = dependentVariableTimeHistoryMap_Ascent.size();
    // long rowsDescent = dependentVariableTimeHistoryMap_Descent.size();
    long columns = ( ( ascentSimulation.getDependentVariableHistory( ).begin() )->second ).size();


    // std::cout << "Aqui" << std::endl;

    //std::cout << rowsAscent + rowsDescent << std::endl;

    Eigen::MatrixXd dependentVariableTimeHistoryMatrix( dependentVariableTimeHistoryMap.size() , columns );
    dependentVariableTimeHistoryMatrix = Eigen::MatrixXd::Zero( dependentVariableTimeHistoryMap.size() , columns );

    //    for ( std::map< double, Eigen::VectorXd >::iterator it = dependentVariableTimeHistoryMatrix.begin(); it != dependentVariableTimeHistoryMatrix.end(); ++it )
    for ( unsigned long i = 1; i < dependentVariableTimeHistoryMap.size(); i++ )
    {
        dependentVariableTimeHistoryMatrix.row( i - 1 ) = dependentVariableTimeHistoryMap.at( simulationStartEpoch + ( i - 1 ) * fixedStepSize );
    }

    //std::cout << "Aqui!!!!" << std::endl;

    //  rows = dependentVariableTimeHistoryMap_Descent.size();
    //  columns = ( ( descentSimulation.getDependentVariableHistory( ).begin() )->second ).size();

    Eigen::VectorXd dependentVariable_FlightPathAngle = dependentVariableTimeHistoryMatrix.col( 7 );
    Eigen::VectorXd dependentVariable_Height = dependentVariableTimeHistoryMatrix.col( 11 );
    Eigen::VectorXd dependentVariable_Airspeed = dependentVariableTimeHistoryMatrix.col( 13 );
    Eigen::VectorXd dependentVariable_MechanicalLoad = dependentVariableTimeHistoryMatrix.col( 14 );
    Eigen::VectorXd dependentVariable_DynamicPressure = dependentVariableTimeHistoryMatrix.col( 24 );
    Eigen::VectorXd dependentVariable_HeatingRate = dependentVariableTimeHistoryMatrix.col( 25 );
    Eigen::VectorXd dependentVariable_CurrentMass = dependentVariableTimeHistoryMatrix.col( 27 );
    Eigen::VectorXd dependentVariable_NormalizedSpecificEnergy = dependentVariableTimeHistoryMatrix.col( 29 );

    //! Extract Dependent variables of final state.
    const Eigen::VectorXd dependentVariableFinalState_Descent = ( descentSimulation.getDependentVariableHistory( ).rbegin() )->second;

    //const double altitude_f_calc = dependentVariable_FINALSTATE[3];
    const double targetLat_rad_calc  = dependentVariableTimeHistoryMatrix( dependentVariableTimeHistoryMatrix.rows() - 1, 4 );
    const double targetLon_rad_calc  = dependentVariableTimeHistoryMatrix( dependentVariableTimeHistoryMatrix.rows() - 1, 5 );
    const double E_hat_f_calc        = dependentVariableTimeHistoryMatrix( dependentVariableTimeHistoryMatrix.rows() - 1, 29 );
    std::ptrdiff_t ind_maximum_height;
    h_UP_calc                    = dependentVariable_Height.maxCoeff( &ind_maximum_height );
    upperLimitAirspeed_calc                    = dependentVariable_Airspeed( ind_maximum_height );
    h_DN_calc                    = dependentVariableTimeHistoryMatrix( dependentVariableTimeHistoryMatrix.rows() - 1, 11 );
    lowerLimitAirspeed_calc                    = dependentVariableTimeHistoryMatrix( dependentVariableTimeHistoryMatrix.rows() - 1, 13 );

    std::ptrdiff_t index_MaximumNormalizedSpecificEnergy;
    maximum_NormalizedSpecificEnergy  = dependentVariable_NormalizedSpecificEnergy.maxCoeff( &index_MaximumNormalizedSpecificEnergy );
    std::ptrdiff_t index_MaximumAirspeed;
    maximumAirspeed  = dependentVariable_Airspeed.maxCoeff( &index_MaximumAirspeed );
    const double actualdelV = maximumAirspeed - initialAirspeed;
    finalMass_Descent = dependentVariableFinalState_Descent[ 27 ];

    //std::ptrdiff_t index_MaximumDynamicPressure;
    //maximum_DynamicPressure               = dependentVariable_DynamicPressure.maxCoeff( &index_MaximumDynamicPressure );


    //std::cout << "lat_f_rad_calc: " << lat_f_rad_calc << std::endl;
    //std::cout << "lon_f_rad_calc: " << lon_f_rad_calc << std::endl;
    //std::cout << "h_UP_calc     : " << h_UP_calc << std::endl;
    //std::cout << "h_DN_calc     : " << h_DN_calc << std::endl;

    //const double altitude_f_calc = std::sqrt( pow(systemFinalState_EARTh_UPIXED[0],2) +
    //       pow(systemFinalState_EARTh_UPIXED[1],2) + pow(systemFinalState_EARTh_UPIXED[2],2) ) ;
    //const double lon_f_rad_calc = std::atan2(systemFinalState_EARTh_UPIXED[1] , systemFinalState_EARTh_UPIXED[0]);
    //const double lat_f_rad_calc = std::asin(systemFinalState_EARTh_UPIXED[2] / altitude_f_calc) ;

    //! Convert coordinates of final state to degrees: Earth-Fixed Frame
    targetLat_deg_calc = unit_conversions::convertRadiansToDegrees( targetLat_rad_calc );
    targetLon_deg_calc = unit_conversions::convertRadiansToDegrees( targetLon_rad_calc );

    //! Calculate angular distance of final state from target coordinates.
    // const double d_rad = bislip::Variables::computeAngularDistance( lat_f_rad_calc, lon_f_rad_calc, vehicleSystems->getTargetLat(), vehicleSystems->getTargetLon() );
    const double d_to_go_deg = unit_conversions::convertRadiansToDegrees( dependentVariableFinalState_Descent[ 35 ] );
    const double d_travelled_deg = unit_conversions::convertRadiansToDegrees( dependentVariableFinalState_Descent[ 34 ] );
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
    const double dif_h_UP = h_UP - h_UP_calc;

    const double dif_upperLimitAirspeed = abs( upperLimitAirspeed - upperLimitAirspeed_calc );

    //! Calculate offset from final elevation.
    const double dif_h_DN = h_DN - h_DN_calc;
    const double dif_lowerLimitAirspeed = abs( lowerLimitAirspeed - lowerLimitAirspeed_calc );


    // std::cout << "Penalty: Flight Path - Ascent" << std::endl;
    const double penaltyFlightPathAngleAscent = bislip::Variables::computePenalty( ( -1.0 ) * dependentVariable_FlightPathAngle, 1, rowsAscent, 0, fixedStepSize, tof, true ) * fixedStepSize / ( 2 * mathematical_constants::PI );

    // std::cout << "Penalty: Flight Path - Descent" << std::endl;
    const double penaltyFlightPathAngleDescent = bislip::Variables::computePenalty( dependentVariable_FlightPathAngle, rowsAscent, dependentVariable_FlightPathAngle.size(), 0, fixedStepSize, tof, true ) * fixedStepSize / ( 2 * mathematical_constants::PI );

    //std::cout << "Penalty: Monotonic - Ascent" << std::endl;
    const double penaltyMonotonicAscent = bislip::Variables::computePenalty( dependentVariable_NormalizedSpecificEnergy, 1, rowsAscent, 0, fixedStepSize, tof, true );

    //std::cout << "Penalty: Monotonic - Descent" << std::endl;
    const double penaltyMonotonicDescent = bislip::Variables::computePenalty( dependentVariable_NormalizedSpecificEnergy, rowsAscent, dependentVariable_NormalizedSpecificEnergy.size(), 0, fixedStepSize, tof, true );

    //std::cout << "Penalty: Dynamic Pressure - Ascent" << std::endl;
    const double penaltyDynamicPressureAscent = bislip::Variables::computePenalty( dependentVariable_DynamicPressure, 1, rowsAscent, constraint_DynamicPressure, fixedStepSize, tof, false );

    //std::cout << "Penalty: Dynamic Pressure - Descent" << std::endl;
    const double penaltyDynamicPressureDescent = bislip::Variables::computePenalty( dependentVariable_DynamicPressure, rowsAscent, dependentVariable_DynamicPressure.size(), constraint_DynamicPressure, fixedStepSize, tof, false );

    //std::cout << "Penalty: Mechanical Load - Ascent" << std::endl;
    const double penaltyMechanicalLoadAscent = bislip::Variables::computePenalty( dependentVariable_MechanicalLoad, 1, rowsAscent, constraint_MechanicalLoad, fixedStepSize, tof, false );

    //std::cout << "Penalty: Mechanical Load - Descent" << std::endl;
    const double penaltyMechanicalLoadDescent = bislip::Variables::computePenalty( dependentVariable_MechanicalLoad, rowsAscent, dependentVariable_MechanicalLoad.size(), constraint_MechanicalLoad, fixedStepSize, tof, false );

    //std::cout << "Penalty: Heating Rate - Ascent" << std::endl;
    const double penaltyHeatingRateAscent = bislip::Variables::computePenalty( dependentVariable_HeatingRate, 1, rowsAscent, constraint_HeatingRate, fixedStepSize, tof, false );

    //std::cout << "Penalty: Heating Rate - Descent" << std::endl;
    const double penaltyHeatingRateDescent = bislip::Variables::computePenalty( dependentVariable_HeatingRate, rowsAscent, dependentVariable_HeatingRate.size(), constraint_HeatingRate, fixedStepSize, tof, false );

    //! Calculate offset of final angular distance to termination condition distance.
    const double costDistanceTravelled                            = 10000.0 * ( bodyMap_.at( vehicleName_ )->getBislipSystems()->getInitialDistanceToTarget( ) - dependentVariableFinalState_Descent[ 36 ] ) / ( bodyMap_.at( vehicleName_ )->getBislipSystems()->getInitialDistanceToTarget( ) );

    //     = 1000 * ( d_to_go_deg - finalDistanceToTarget_deg) / finalDistanceToTarget_deg ;

    //std::cout << "Cost: Heating Rate - Ascent" << std::endl;
    const double costHeatingRateAscent                            = bislip::Variables::computePenalty( dependentVariable_HeatingRate, 1, rowsAscent, constraint_HeatingRate, fixedStepSize, tof, true ) * ( fixedStepSize / ( tof * constraint_HeatingRate ) );

    //std::cout << "Cost: Heating Rate - Descent" << std::endl;
    const double costHeatingRateDescent                           = bislip::Variables::computePenalty( dependentVariable_HeatingRate, rowsAscent, dependentVariable_HeatingRate.size(), constraint_HeatingRate, fixedStepSize, tof, true ) * ( fixedStepSize / ( tof * constraint_HeatingRate ) );

    const double costConsumedMassAscent                           = abs( initialMass_Ascent - dependentVariable_CurrentMass( rowsAscent - 1 ) ) / initialMass_Ascent;
    const double costConsumedMassDescent                          = abs( dependentVariable_CurrentMass( rowsAscent ) - finalMass_Descent ) / initialMass_Ascent;


    //const double penaltyLatitudeDifference = abs( targetLat_deg - targetLat_deg_calc );
    //const double penaltyLongitudeDifference = abs( targetLon_deg - targetLon_deg_calc );
    const double penaltyCoordinateOffsets                         = std::sqrt( ( dif_lat_deg * dif_lat_deg ) + ( dif_lon_deg * dif_lon_deg ) ) / initialDistanceToTarget_deg;

    //const double penaltyFinalNodeMagnitudeAscent                  = 1000 * ( xn_Ascent( xn_Ascent.size() - 1 ) - 1 );
    //const double penaltyInitialNodeMagnitudeDescent               = 1000 * ( xn_Descent( 0 ) - 1 );
    //const double penaltyFinalNodeMagnitudeDescent                 = 1000 * abs( xn_Descent( xn_Descent.size() - 1 ) );
    const double penaltyMaximumNormalizedSpecificEnergyAscent     = 1000.0 * abs( E_mapped_Ascent( E_mapped_Ascent.size() - 1 ) - maximum_NormalizedSpecificEnergy );
    const double penaltyMinimumNormalizedSpecificEnergyDescent    = abs( E_hat_f_calc - E_mapped_Descent( E_mapped_Descent.size() - 1 ) );

    //const double penaltyFinalMappedNormalizedSpecificEnergyAscent    = 1000 * ( E_mapped_Ascent( E_mapped_Ascent.size() - 1 ) - 1 );
    //const double penaltyInitialMappedNormalizedSpecificEnergyDescent = 1000 * ( E_mapped_Descent( 0 ) - 1 );
    //const double penaltyFinalMappedNormalizedSpecificEnergyDescent   = abs( E_mapped_Descent(  E_mapped_Descent.size() - 1 ) ) - E_mapped_Descent(  E_mapped_Descent.size() - 1 );
    //const double penaltyMappedNormalizedSpecificEnergyInterface      = abs( E_mapped_Descent( 0 ) - E_mapped_Ascent( E_mapped_Ascent.size() - 1 ) );

    const double penaltyNoFlight                                  = 1000.0 * ( 1.0 - maximum_NormalizedSpecificEnergy ) * ( 1.0 - maximum_NormalizedSpecificEnergy );
    const double penaltydelV                                      = 1000.0 * abs( goaldelV - actualdelV ) / theoreticaldelV;

    //! Assign values to Fitness vector! At the moment these are all 'objective
    //! functions'. To modify this I have change the header file and define how
    //! many are elements are objective function, equality contraints, and
    //! inequality constraints. This vector here must contain them is that exact order: nOF, nEC, nIC.
    std::vector< double > delta;

    /*
    delta.push_back( penaltyNoFlight );

    delta.push_back( costDistanceToGo + penaltyCoordinateOffsets );

    delta.push_back( costConsumedMassAscent );
    delta.push_back( penaltyMaximumNormalizedSpecificEnergyAscent + penaltyMonotonicAscent );
    delta.push_back( costHeatingRateAscent + penaltyHeatingRateAscent );
    delta.push_back( penaltyFlightPathAngleAscent + penaltyDynamicPressureAscent + penaltyMechanicalLoadAscent );

    delta.push_back( costConsumedMassDescent );
    delta.push_back( penaltyMinimumNormalizedSpecificEnergyDescent + penaltyMonotonicDescent );
    delta.push_back( costHeatingRateDescent + penaltyHeatingRateDescent );
    delta.push_back( penaltyFlightPathAngleDescent + penaltyDynamicPressureDescent + penaltyMechanicalLoadDescent );

*/

    delta.push_back( penaltyNoFlight + penaltydelV );

    delta.push_back( costDistanceTravelled + penaltyCoordinateOffsets );

    delta.push_back( costConsumedMassAscent +
                     penaltyMaximumNormalizedSpecificEnergyAscent +
                     penaltyMonotonicAscent + costHeatingRateAscent +
                     penaltyHeatingRateAscent + penaltyFlightPathAngleAscent +
                     penaltyDynamicPressureAscent + penaltyMechanicalLoadAscent );

    delta.push_back( costConsumedMassDescent +
                     penaltyMinimumNormalizedSpecificEnergyDescent +
                     penaltyMonotonicDescent + costHeatingRateDescent +
                     penaltyHeatingRateDescent + penaltyFlightPathAngleDescent +
                     penaltyDynamicPressureDescent + penaltyMechanicalLoadDescent );


    //delta.push_back( penaltyFinalNodeMagnitudeAscent + penaltyInitialNodeMagnitudeDescent + penaltyFinalNodeMagnitudeDescent + penaltyFinalMappedNormalizedSpecificEnergyAscent + penaltyInitialMappedNormalizedSpecificEnergyDescent + penaltyFinalMappedNormalizedSpecificEnergyDescent + penaltyMappedNormalizedSpecificEnergyInterface );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PRINT SIMULATION OUTPUT TO FILE               //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Get time stamp for this specific simulation. This avoids overwriting the
    //! file if another individual with the same properties shows up in other
    //! evolutions.
    std::string simulation_save_time = bislip::Variables::getCurrentDateTime( false );

    //! Create unique filename that cannot be overwritten due to the timestamp.
    std::string simulation_file_name_suffix =
            // std::to_string( Va_ordered( 0 ) ) + "_" +
            // std::to_string( gamma_rad_ordered( 0 ) ) + "_" +
            // std::to_string( chi_i_deg_ordered( 0 ) ) + "_" +
            std::to_string( tof ) + "_" +
            //std::to_string(x[3]) + "_" +
            std::to_string( targetLat_deg_calc ) + "_" +
            std::to_string( targetLon_deg_calc ) + "_" +
            simulation_save_time;


    for( unsigned int i = 0; i < delta.size() - 1; i++)
    {
        std::cout  << delta[ i ] << "  |  " ;
    }
    //std::cout  << delta[ delta.size() - 1 ] << "  ||  " << rowsAscent << "  ||  " << tof << "  ||  " << " Theoretical delV = " << theoreticaldelV << "  ||  " << " Goal delV = " << goaldelV << "  ||  " << " Actual delV = " << actualdelV << std::endl;
//    std::cout  << delta[ delta.size() - 1 ] << "  ||  " << rowsAscent << "  ||  " << tof <<  std::endl;
    std::cout  << delta[ delta.size() - 1 ] << "  ||  " << tof <<  std::endl;

    /*
    //! Print results to terminal. Used to gauge progress.
    std::cout << std::fixed << std::setprecision(10) <<
                 std::setw(15) << " initialAirspeed = " <<
                 std::setw(14) << initialAirspeed <<
                 std::setw(33) << " costFinalNodeMagnitudeAscent = " <<
                 std::setw(14) << costFinalNodeMagnitudeAscent<<
                 std::setw(40) << " costFinalMappedNormalizedSpecificEnergyAscent = " <<
                 std::setw(16) << costFinalMappedNormalizedSpecificEnergyAscent <<
                 std::setw(36) << " costInitialNodeMagnitudeDescent = " <<
                 std::setw(14) << costInitialNodeMagnitudeDescent<<
                 std::setw(53) << " costInitialMappedNormalizedSpecificEnergyDescent = " <<
                 std::setw(16) << costInitialMappedNormalizedSpecificEnergyDescent <<
                 std::setw(15) << " dif_mass = " <<
                 std::setw(16) << costConsumedMass <<
                 std::setw(15) << " costConsumedMass = " <<
                 std::setw(16) << dif_norm <<
                 std::setw(15) << " dif_lat_deg = " <<
                 std::setw(16) << dif_lat_deg <<
                 std::setw(15) << " dif_lon_deg = " <<
                 std::setw(16) << dif_lon_deg <<
                 std::setw(21) << " costDistanceToGo = " <<
                 std::setw(16) << costDistanceToGo <<
                 std::setw(9) << " dif_h_UP = " <<
                 std::setw(16) << dif_h_UP <<
                 std::setw(9) << " dif_h_DN = " <<
                 std::setw(16) << dif_h_DN <<
                 std::setw(7) << " tof = " <<
                 std::setw(16) << tof <<
                 std::setw(60) << simulation_file_name_suffix << std::endl;
                 */

    std::string complete_file_name_Prop = "HORUSPropHistory_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_DepVar = "HORUSDepVar_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_interpolators_Ascent = "interpolators_Ascent_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_map_DV_mapped_Ascent = "map_DV_mapped_Ascent_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_interpolators_Descent = "interpolators_Descent_" + simulation_file_name_suffix + ".dat";
    std::string complete_file_name_map_DV_mapped_Descent = "map_DV_mapped_Descent_" + simulation_file_name_suffix + ".dat";

    //! Will print out depending on some input values. Each entry corresponds to
    //! a different type of output. Entries are binary, 0 or 1.
    if ( int( output_settingsValues_[ 1 ] ) == 1 )
    {
        //   std::cout << "Saving EQM " << std::endl;

        //! Write propagation history to file.
        writeDataMapToTextFile( propagationTimeHistoryMap,
                                complete_file_name_Prop,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //   std::cout << "EQM Saved" << std::endl;

    }
    if ( int( output_settingsValues_[ 2 ] ) == 1 )
    {
        //   std::cout << "Saving DepVar " << std::endl;

        //! Write dependent variables' history to file.
        writeDataMapToTextFile( dependentVariableTimeHistoryMap,
                                complete_file_name_DepVar,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //    std::cout << "DepVar Saved" << std::endl;


    }
    if ( int( output_settingsValues_[ 3 ] ) == 1 )
    {
        //   std::cout << "Saving interpolators_Ascent " << std::endl;

        writeDataMapToTextFile( interpolators_Ascent,
                                complete_file_name_interpolators_Ascent,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //   std::cout << "interpolators_Ascent Saved" << std::endl;
        //   std::cout << "Saving map_DV_mapped_Ascent " << std::endl;

        writeDataMapToTextFile( map_DV_mapped_Ascent,
                                complete_file_name_map_DV_mapped_Ascent,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //   std::cout << "map_DV_mapped_Ascent Saved" << std::endl;

        //   std::cout << "Saving interpolators_Descent " << std::endl;

        writeDataMapToTextFile( interpolators_Descent,
                                complete_file_name_interpolators_Descent,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //  std::cout << "interpolators_Descent Saved" << std::endl;
        //   std::cout << "Saving map_DV_mapped_Descent " << std::endl;

        writeDataMapToTextFile( map_DV_mapped_Descent,
                                complete_file_name_map_DV_mapped_Descent,
                                tudat_applications::getOutputPath( ) + outputSubFolder_,
                                "",
                                std::numeric_limits< double >::digits10,
                                std::numeric_limits< double >::digits10,
                                "," );
        //std::cout << "map_DV_mapped_Descent Saved" << std::endl;


    }


    //  }



    //std::cout << "h_UP: " << h_UP << std::endl;
    //std::cout << "h_UP_calc: " << h_UP_calc << std::endl;



    // std::cout <<dynamicsSimulator.getPropagationTerminationReason( )->getPropagationTerminationReason( ) << std::endl;
    //std::this_thread::sleep_for( std::chrono::nanoseconds( 10 ) );
    //using namespace std::this_thread; // sleep_for, sleep_until
    //using namespace std::chrono; // nanoseconds, system_clock, seconds




    return delta;

} // Fitness function.
