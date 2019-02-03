#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

namespace bislip {
double getAngularDistance (
	const double &lat_i,
    const double &lon_i,
    const double &lat_f,
    const double &lon_f)
{
    const double angular_distance = std::acos( std::sin(lat_i) * std::sin(lat_f) + std::cos(lat_i) * std::cos(lat_f) * std::cos(lon_f - lon_i) );

    return angular_distance;
}
}
