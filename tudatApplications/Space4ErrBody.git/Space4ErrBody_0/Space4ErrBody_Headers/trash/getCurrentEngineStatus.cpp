#include "Space4ErrBody.h"
#include "updateGuidance.h"



namespace bislip {

bool MyGuidance::getCurrentEngineStatus( )
{

    std::cout << "getCurrentEngineStatus"  << std::endl;

    if ( currentMass_ > finalMass_ )
    {
        currentEngineStatus_ = 1;
    }
    else
    {
        currentEngineStatus_ = 0;
    }
    return currentEngineStatus_;
}
} // namespace aerodynamics
//} // namespace tudat
