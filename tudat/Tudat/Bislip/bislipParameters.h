#ifndef BISLIPPARAMETERS_H
#define BISLIPPARAMETERS_H

#include <Tudat/Bislip/bislipHeaders.h>


namespace bislip {

namespace Parameters {

enum Optimization
{
    AngleOfAttack,
    BankAngle,
    ThrustElevationAngle,
    ThrustAzimuthAngle,
    ThrottleSetting,
    NodeInterval,
    InitialVelocity,
    MaximumVelocity,
    MaximumHeight,
    AdditionalMass
};

}; // namespace Parameters

} // namespace bislip

#endif // BISLIPPARAMETERS_H
