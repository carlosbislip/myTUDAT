#include "getStuff.h"
#include "updateGuidance.h"

#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

namespace tudat { namespace aerodynamics {

double MyAerodynamicGuidance::getF (
        unsigned i,
        unsigned j,
        const Eigen::VectorXd &p,
        const double &eps,
        Eigen::Vector2d &x ) {

Eigen::Vector6d newCoefficients;
double del_x, del_x_T, del_z_T, S_ref, c_ref, Isp, m, q_d, current_M, n, g0, extra;
del_x     = p[0];
del_x_T   = p[1];
del_z_T   = p[2];
S_ref     = p[3];
c_ref     = p[4];
Isp       = p[5];
m         = p[6];
q_d       = p[7];
current_M = p[8];
n         = p[9];
g0        = p[10];
extra     = p[11];

double m_dot, f;
//std::cout << "x = "<< x[0] << "   " << x[1] << std::endl;









if ( i == 0 and j == 0 )
{
    x[0] = x[0] + eps;
    newCoefficients = MyAerodynamicGuidance::getNewCoefficients( x[0], current_M );
    m_dot = ( ( q_d * S_ref ) / ( g0 * Isp ) ) * ( ( del_x * ( newCoefficients[0] * std::sin( x[0] ) - newCoefficients[2] * std::cos ( x[0] ) ) - newCoefficients[4] * c_ref ) ) / ( del_x_T * std::sin( x[1] ) + del_z_T * std::cos( x[1] ) );
    f = n * m * g0  - sqrt( pow( m_dot * g0 * Isp * std::cos( x[0] + x[1] ) - q_d * S_ref * newCoefficients[0] , 2 ) + pow( m_dot * g0 * Isp * std::sin( x[0] + x[1] ) - q_d * S_ref * newCoefficients[2] , 2 ) );
}
if ( i == 0 and j == 1 )
{
    x[1] = x[1] + eps;
    newCoefficients = MyAerodynamicGuidance::getNewCoefficients( x[0], current_M );
    m_dot = ( ( q_d * S_ref ) / ( g0 * Isp ) ) * ( ( del_x * ( newCoefficients[0] * std::sin( x[0] ) - newCoefficients[2] * std::cos ( x[0] ) ) - newCoefficients[4] * c_ref ) ) / ( del_x_T * std::sin( x[1] ) + del_z_T * std::cos( x[1] ) );
    f = n * m * g0  - sqrt( pow( m_dot * g0 * Isp * std::cos( x[0] + x[1] ) - q_d * S_ref * newCoefficients[0] , 2 ) + pow( m_dot * g0 * Isp * std::sin( x[0] + x[1] ) - q_d * S_ref * newCoefficients[2] , 2 ) );
}
if ( i == 1 and j == 0 )
{
    x[0] = x[0] + eps;
    newCoefficients = MyAerodynamicGuidance::getNewCoefficients( x[0], current_M );
    m_dot = ( ( q_d * S_ref ) / ( g0 * Isp ) ) * ( ( del_x * ( newCoefficients[0] * std::sin( x[0] ) - newCoefficients[2] * std::cos ( x[0] ) ) - newCoefficients[4] * c_ref ) ) / ( del_x_T * std::sin( x[1] ) + del_z_T * std::cos( x[1] ) );
    f = n * m * g0 - ( m_dot * g0 * Isp * std::cos( x[0] + x[1] ) - q_d * S_ref * newCoefficients[0] ) + extra;
}
if ( i == 1 and j == 1 )
{
    x[1] = x[1] + eps;
    newCoefficients = MyAerodynamicGuidance::getNewCoefficients( x[0], current_M );
    m_dot = ( ( q_d * S_ref ) / ( g0 * Isp ) ) * ( ( del_x * ( newCoefficients[0] * std::sin( x[0] ) - newCoefficients[2] * std::cos ( x[0] ) ) - newCoefficients[4] * c_ref ) ) / ( del_x_T * std::sin( x[1] ) + del_z_T * std::cos( x[1] ) );
    f = n * m * g0 - ( m_dot * g0 * Isp * std::cos( x[0] + x[1] ) - q_d * S_ref * newCoefficients[0] ) + extra;
}

//std::cout << "Mach number = "<< current_M << std::endl;
std::cout << "new Aero Coefficients = "<< newCoefficients[0] << "   " << newCoefficients[2]  << "   " << newCoefficients[4]  << std::endl;
//std::cout << "x = "<< x[0] << "   " << x[1] << std::endl;
std::cout << "m_dot = "<< m_dot << std::endl;

/*
if ( eps != 0 )
{
std::cout << "J( "<<i<<" , "<<j<<" ) = "<< f << std::endl;
}
else
{
    std::cout << "f( "<<i<<" ) = "<< f << std::endl;
 }
 */
return f;

}

} // namespace aerodynamics
} // namespace tudat

