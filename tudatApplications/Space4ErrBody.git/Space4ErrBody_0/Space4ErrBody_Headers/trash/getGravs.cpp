#include "getStuff.h"

namespace bislip {
Eigen::Vector2d getGravs (
        const double mu,
        const double J2,
        const double J3,
        const double J4,
        const double R_E,
        const double r,
        const double delta_rad ) {

    const double s_delta = std::sin( delta_rad );
    const double c_delta = std::cos( delta_rad );

    const double g_n_pre = -3 * ( mu / ( r * r ) ) * pow( R_E / r , 2 ) * s_delta * c_delta;
    const double g_n_sub_1 = J2;
    const double g_n_sub_2 = ( J3 / 2 ) * ( R_E / r ) * ( 1 / std::sin( delta_rad) ) * ( 5 * pow ( s_delta , 2 ) - 1 );
    const double g_n_sub_3 = ( 5 * J4 / 6 ) * pow( R_E / r , 2) * ( 7 * pow ( s_delta , 2 ) - 3 );
    const double g_n = g_n_pre * ( g_n_sub_1 + g_n_sub_2 + g_n_sub_3 );

    const double g_d_pre = ( mu / ( r * r ) );
    const double g_d_sub_1 = -( 3 * J2 / 2 ) * pow( R_E / r , 2 ) * ( 3 * pow  (s_delta , 2 ) - 1 );
    const double g_d_sub_2 = -2 * J3  * pow( R_E / r , 3 ) * ( s_delta ) * ( 5 * pow ( s_delta , 2 ) - 3 );
    const double g_d_sub_3 = -( 5 * J4 / 8 ) * pow( R_E / r , 4 ) * ( 35 * pow ( s_delta , 4 ) - 30 * pow ( s_delta , 2 ) + 3 );
    const double g_d = g_d_pre * ( 1 + g_d_sub_1 + g_d_sub_2 + g_d_sub_3 );

    Eigen::Vector2d gravs ( 2 );
    gravs << g_n, g_d;

    return gravs;

}
}
