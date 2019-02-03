#include "Space4ErrBody.h"
#include "updateGuidance.h"



namespace bislip {

double MyGuidance::getCurrentSpecificImpulse( )
{

    std::cout << "getCurrentSpecificImpulse"  << std::endl;

    currentSpecificImpulse_ = Isp_;

    return currentSpecificImpulse_;
}
} // namespace aerodynamics
//} // namespace tudat
