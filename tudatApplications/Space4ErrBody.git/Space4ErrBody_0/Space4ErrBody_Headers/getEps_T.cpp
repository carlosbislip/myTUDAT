#include "getStuff.h"

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace bislip {
double getEps_T ( const double &AoA,
                  const Eigen::VectorXd &p,
                  const Eigen::Vector6d &newCoefficients  ) {


    double s_x_new, c_x_new, s_AoA_x_new, s_2AoA_2x_new, c_2AoA_2x_new, s_AoA_x_new_phi, c_AoA_x_new_phi, s_AoA_x_new_sq, Thrust_terms_1, Thrust_terms_2;
    double del_x, del_x_T, del_z_T, S_ref, c_ref, Isp, m, q_d, current_M, n, g0, extra;
    double A, B, C, phi, f, fp_pre, fp_1, fp_2, fp, f2p_1, f2p_2, f2p;
    double dif;

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
    const double s_AoA = std::sin( AoA );
    const double c_AoA = std::cos( AoA );
    const double g0_Isp = g0 * Isp;
    const double extra_sq = extra * extra;
    const double q_d_S_ref = q_d * S_ref;
    const double q_d_S_ref_sq = q_d_S_ref * q_d_S_ref;

    double x_new = 0.1;
    double x_old = 0;
    int sgn = 1;

    do {

        A   = g0_Isp * ( del_x * ( newCoefficients[0] * s_AoA - newCoefficients[2] * c_AoA ) - newCoefficients[4] * c_ref );
        B   = sqrt( q_d_S_ref_sq * newCoefficients[2] * newCoefficients[2] +  extra_sq ) ;
        C   = q_d_S_ref_sq * newCoefficients[2] * newCoefficients[2] + q_d_S_ref * newCoefficients[0] * extra - extra_sq;
        phi = std::atan2( extra, q_d_S_ref *  newCoefficients[2] );

        s_x_new         = std::sin( x_new );
        c_x_new         = std::cos( x_new );
        s_AoA_x_new     = std::sin( AoA + x_new );
        s_AoA_x_new_sq  = s_AoA_x_new * s_AoA_x_new;
        s_2AoA_2x_new   = std::sin( 2 * ( AoA + x_new ) );
        c_2AoA_2x_new   = std::cos( 2 * ( AoA + x_new ) );
        Thrust_terms_1  = del_x_T * s_x_new + del_z_T * c_x_new;
        Thrust_terms_2  = del_x_T * c_x_new - del_z_T * s_x_new;
        s_AoA_x_new_phi = std::sin( AoA + x_new + phi );
        c_AoA_x_new_phi = std::cos( AoA + x_new + phi );

        f      = A * ( ( s_AoA_x_new_sq - 2 * B * s_AoA_x_new_phi ) / Thrust_terms_1 ) + C;

        fp_pre = A / ( Thrust_terms_1 * Thrust_terms_1 ) ;
        fp_1   = Thrust_terms_1 * ( s_2AoA_2x_new - 2 * B * c_AoA_x_new_phi );
        fp_2   = Thrust_terms_2 * ( s_AoA_x_new_sq - 2 * B * s_AoA_x_new_phi );
        fp     = fp_pre * ( fp_1 -fp_2 );

        f2p_1  = Thrust_terms_1 * Thrust_terms_1 * Thrust_terms_1 * A * ( 2 * c_2AoA_2x_new + s_AoA_x_new_sq );
        f2p_2  = 2 * A * ( fp_1 - fp_2 ) * Thrust_terms_1  * Thrust_terms_2;
        f2p    = ( f2p_1 - f2p_2 ) / ( Thrust_terms_1 * Thrust_terms_1 * Thrust_terms_1 * Thrust_terms_1);

        x_old = x_new;
        x_new = bislip::getRoot_Halley( x_old, f, fp, f2p );

        if ( std::abs( x_new ) > 90.0 * tudat::mathematical_constants::PI / 180 )
        {
            if ( x_new < 0 )
            {
                sgn = -1;
            }
            x_new = sgn * 90.0 * tudat::mathematical_constants::PI / 180;
        }
        //std::cout <<" old: " << x_old<< " new: " << x_new  << std::endl;
        dif = x_new - x_old;

    } while( abs( dif ) > 1e-10 );

return x_new ;

}
}
