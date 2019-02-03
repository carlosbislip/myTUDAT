#include "getStuff.h"

namespace bislip  {

double computeSpecificEnergy(
        const double height,
        const double airspeed)
{
    return 9.80665 * height + 0.5 * airspeed * airspeed;
}

} // namespace bislip
