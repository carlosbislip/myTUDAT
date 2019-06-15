#ifndef BISLIPHEADERS_H
#define BISLIPHEADERS_H

#include <ctime>
#include <sstream>
#include <iomanip>
#include <utility>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>      // for sprintf()
#include <string>       // for std::string
#include <boost/date_time/posix_time/posix_time.hpp>

#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Mathematics/Interpolators/createInterpolator.h>
//#include <Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>
#include <Tudat/Astrodynamics/SystemModels/vehicleSystems.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h>
//#include <Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h>
#include <Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h>
#include <Tudat/InputOutput/solarActivityData.h>
#include <Tudat/Astrodynamics/Aerodynamics/nrlmsise00InputFunctions.h>
#include <Tudat/Mathematics/BasicMathematics/basicFunction.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
//#include <Tudat/External/SpiceInterface/spiceInterface.h>
#include <Tudat/Astrodynamics/ElectroMagnetism/basicElectroMagnetism.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h>
#include <Tudat/SimulationSetup/PropagationSetup/propagationTerminationSettings.h>

//#include <iostream>

#include <boost/filesystem.hpp>
//#include "applicationOutput_pagmo.h"

#include <Tudat/Basics/utilities.h>
#include <Tudat/InputOutput/basicInputOutput.h>
//#include <pagmo/problem.hpp>
//#include <pagmo/algorithms/nsga2.hpp>
//#include <pagmo/algorithms/moead.hpp>
//#include <pagmo/algorithms/ihs.hpp>

#include <limits>
//#include <iostream>
//#include <fstream>
//#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
//#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>// std::pair, std::get
#include <boost/filesystem/operations.hpp>
#include <chrono>
#include <thread>
#include <functional>

#include <Eigen/Core>



#endif // BISLIPHEADERS_H
