#ifndef GETGAMMA_H
#define GETGAMMA_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
namespace tudat
{

double getGamma ( const simulation_setup::NamedBodyMap& bodyMap,
                 const std::string vehicleName,
                 const std::vector< double > term_cond);

} // namespace tudat


#endif // GETGAMMA_H

