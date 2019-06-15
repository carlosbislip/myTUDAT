/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef SPACE4ERRBODY_HPP
#define SPACE4ERRBODY_HPP
/*
#include <utility>
#include <vector>
#include <limits>
#include <iostream>
//#include <iomanip>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <cstring>
//#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
//#include <sstream>
#include <utility>// std::pair, std::get
#include <boost/filesystem/operations.hpp>
#include <chrono>
#include <thread>

#include <Eigen/Core>

#include <pagmo/problem.hpp>
#include <pagmo/io.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/island.hpp>
#include <pagmo/serialization.hpp>
#include <pagmo/types.hpp>

//#include "getStuff.h"
//#include <Tudat/Bislip/Parameters/bislipOptimizationParameters.h>
//#include <Tudat/Bislip/Variables/bislipVariables.h>
//#include <Tudat/Bislip/SystemModels/bislipSystems.h>

#include <boost/format.hpp>


//#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <Tudat/Mathematics/RootFinders/secantRootFinder.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>
#include <Tudat/Astrodynamics/Aerodynamics/nrlmsise00Atmosphere.h>
#include <Tudat/Astrodynamics/Aerodynamics/flightConditions.h>
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamicGuidance.h>
#include <Tudat/Mathematics/Interpolators/createInterpolator.h>
#include <Tudat/Mathematics/Interpolators/interpolator.h>
#include <Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h>
#include <Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h>
#include <Tudat/Astrodynamics/SystemModels/vehicleSystems.h>
*/
#include "updateGuidance.h"
#include <Tudat/Bislip/bislipUtilities.h>
#include <Tudat/Bislip/bislipProblemInput.h>

#include <pagmo/pagmo.hpp>

//#include "getPagmoAlgorithm.h"
//#include "bislipUtilities.h"
//#include "applicationOutput_tudat.h"
//#include "applicationOutput_pagmo.h"


using namespace pagmo;
//using namespace pagmo::problem;

//! Define the problem PaGMO-style
struct Space4ErrBodyProblem {

    //! Empty constructor
    Space4ErrBodyProblem( ) { }

    //! Constructor.
    Space4ErrBodyProblem( const std::shared_ptr< bislip::ProblemInput > &problemInput,
                          const tudat::simulation_setup::NamedBodyMap& bodyMap );

    //! Calculate the fitness as a function of the parameter vector input_data
    std::vector< double > fitness( const std::vector< double >  &x ) const;

    //! Retrieve the allowable limits of the parameters to vary: pair
    //! containing minima and maxima of parameter values
    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    //! Retrieve the name of the problem
    std::string get_name( ) const;

    std::vector< double >::size_type get_nobj() const
    {
        //return 33u;
        return 3u;
    }
    std::vector< double >::size_type get_nec() const
    {
        return 0u;
    }

    //! Serialization function for Pagmo compatibility
    template < typename Archive >
    void serialize( Archive &ar )
    {
        ar( problemInput_->getDecisionVectorBounds() );
    }

private:

    const std::shared_ptr< bislip::ProblemInput > problemInput_;
    const tudat::simulation_setup::NamedBodyMap bodyMap_;

};

//PAGMO_REGISTER_PROBLEM(pagmo::Space4ErrBody)

#endif // SPACE4ERRBODY_OPTIMIZATION_HPP
