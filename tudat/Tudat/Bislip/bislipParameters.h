#ifndef BISLIPPARAMETERS_H
#define BISLIPPARAMETERS_H

//#include <Tudat/Bislip/bislipHeaders.h>
#include "bislipHeaders.h"


namespace bislip {

namespace Parameters {

enum Bounds
{
    AngleOfAttack1,
    AngleOfAttack2,
    AngleOfAttack3,
    AngleOfAttack4,
    AngleOfAttack5,
    AngleOfAttack6,
    AngleOfAttack7,
    AngleOfAttack8,
    AngleOfAttack9,
    AngleOfAttack10,
    BankAngle1,
    BankAngle2,
    BankAngle3,
    BankAngle4,
    BankAngle5,
    BankAngle6,
    BankAngle7,
    BankAngle8,
    BankAngle9,
    BankAngle10,
    ThrustElevationAngle1,
    ThrustElevationAngle2,
    ThrustElevationAngle3,
    ThrustElevationAngle4,
    ThrustElevationAngle5,
    ThrustElevationAngle6,
    ThrustElevationAngle7,
    ThrustElevationAngle8,
    ThrustElevationAngle9,
    ThrustElevationAngle10,
    ThrustAzimuthAngle1,
    ThrustAzimuthAngle2,
    ThrustAzimuthAngle3,
    ThrustAzimuthAngle4,
    ThrustAzimuthAngle5,
    ThrustAzimuthAngle6,
    ThrustAzimuthAngle7,
    ThrustAzimuthAngle8,
    ThrustAzimuthAngle9,
    ThrustAzimuthAngle10,
    ThrottleSetting1,
    ThrottleSetting2,
    ThrottleSetting3,
    ThrottleSetting4,
    ThrottleSetting5,
    ThrottleSetting6,
    ThrottleSetting7,
    ThrottleSetting8,
    ThrottleSetting9,
    ThrottleSetting10,
    NodeInterval1,
    NodeInterval2,
    NodeInterval3,
    NodeInterval4,
    NodeInterval5,
    NodeInterval6,
    NodeInterval7,
    NodeInterval8,
    NodeInterval9,
    InitialVelocity,
    FinalVelocity,
    MaximumVelocity,
    MaximumHeight,
    AdditionalMass,

};

enum Interpolators
{
    AngleOfAttack,
    BankAngle,
    ThrustElevationAngle,
    ThrustAzimuthAngle,
    ThrottleSetting

};
}; // namespace Parameters

} // namespace bislip

#endif // BISLIPPARAMETERS_H
