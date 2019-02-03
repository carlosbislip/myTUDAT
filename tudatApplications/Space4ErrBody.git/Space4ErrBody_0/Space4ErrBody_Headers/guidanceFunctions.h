#include "updateGuidance.h"


namespace bislip {

Eigen::Vector3d MyGuidance::getCurrentBodyFixedThrustDirection( );

bool MyGuidance::getCurrentEngineStatus( );

double MyGuidance::getCurrentThrustMagnitude( );

Eigen::Vector6d MyGuidance::getCurrentCoefficients( );

} // namespace bislip
