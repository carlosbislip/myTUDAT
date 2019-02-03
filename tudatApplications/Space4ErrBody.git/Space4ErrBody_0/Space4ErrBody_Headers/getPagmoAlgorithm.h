#include <iostream>

/*
#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/pso.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/bee_colony.hpp"
#include "pagmo/algorithms/sga.hpp"
#include "pagmo/algorithms/sea.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/bee_colony.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/nlopt.hpp"
#include "pagmo/algorithms/ihs.hpp"
#include "pagmo/algorithms/xnes.hpp"
*/
#include "pagmo/algorithms/nsga2.hpp"
#include "pagmo/algorithms/moead.hpp"
#include "pagmo/algorithms/ihs.hpp"

pagmo::algorithm getPagmoAlgorithm( const int index )
{
    //! Default multi-objective algorithm. (Arbitrary)
    pagmo::algorithm algo{ pagmo::nsga2( ) };

    switch( index )
    {
    case 0:
    {
        pagmo::algorithm algo{ pagmo::nsga2( ) };
        break;
    }
    case 1:
    {
        pagmo::algorithm algo{ pagmo::moead( ) };
        break;
    }
    case 2:
    {
        pagmo::algorithm algo{ pagmo::ihs( ) };
        break;
    }
    default:
    {
        throw std::runtime_error( "Error, multi-objective pagmo algorithm " + std::to_string( index ) + " was not found." );
    }
    }

    return algo;
}

