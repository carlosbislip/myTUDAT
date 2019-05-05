#ifndef BISLIPUTILITIES_H
#define BISLIPUTILITIES_H
/*
#include <iostream>


#include <boost/filesystem.hpp>
//#include "applicationOutput_pagmo.h"

#include "Tudat/Basics/utilities.h"
#include <Tudat/InputOutput/basicInputOutput.h>
//#include <pagmo/problem.hpp>
//#include <pagmo/algorithms/nsga2.hpp>
//#include <pagmo/algorithms/moead.hpp>
//#include <pagmo/algorithms/ihs.hpp>


#include <pagmo/pagmo.hpp>
*/


#include <Tudat/Bislip/bislipHeaders.h>

#include <pagmo/pagmo.hpp>

namespace bislip {

namespace Utilities {


//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    //    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
    //                                                std::string( "tudatBundle/tudatApplications/Space4ErrBody_Executables_testing/Space4ErrBody_Executables_and_Headers_testing/applicationOutput_tudat.h" ).length( ) );
    //    std::string outputPath = reducedPath + "Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/code/SimulationOutput/";
    std::string outputPath = "/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/tudatBundle.git/tudatApplications/Space4ErrBody.git/Space4ErrBody_0/matlab/SimulationOutput/";
    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}


/*
namespace tudat_pagmo_applications
{

//! Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    //    std::string reducedPath = filePath_.substr( 0, filePath_.length( ) -
    //                                std::string( "tudatBundle/tudatApplications/Space4ErrBody_Executables_testing/Space4ErrBody_Executables_and_Headers_testing/applicationOutput_pagmo.h" ).length( ) );
    //    std::string outputPath = reducedPath + "Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/code/SimulationOutput/";
    std::string outputPath = "/Users/bislip/Cloud Storage/OneDrive/School/TUDelft/Space Flight/Thesis/tudatBundle.git/tudatApplications/Space4ErrBody.git/Space4ErrBody_0/matlab/SimulationOutput/";
    if( extraDirectory != "" )
    {
        outputPath += extraDirectory;
    }

    if( outputPath.at( outputPath.size( ) - 1 ) != '/' )
    {
        outputPath += "/";
    }

    return outputPath;
}

}


*/


void printPopulationToFile( const std::vector< std::vector< double > >& population,
                            const std::string fileSuffix,
                            const std::string outputSubFolder,
                            const bool isFitness )
{

    Eigen::MatrixXd matrixToPrint( population.size( ), population.at( 0 ).size( ) );
    for( unsigned int i = 0; i < population.size( ); i++ )
    {
        for( unsigned int j = 0; j < population.at( 0 ).size( ); j++ )
        {
            matrixToPrint( i, j ) = population.at( i ).at( j );
        }
    }

    if( !isFitness )
    {
        tudat::input_output::writeMatrixToFile( matrixToPrint,
                                                "population_" + fileSuffix + ".dat",
                                                16,
                                                bislip::Utilities::getOutputPath( ) + outputSubFolder);
    }
    else
    {
        tudat::input_output::writeMatrixToFile( matrixToPrint,
                                                "fitness_" + fileSuffix + ".dat",
                                                16,
                                                bislip::Utilities::getOutputPath( ) + outputSubFolder);
    }
}



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


} // namespace Utilities
} // namespace bislip


#endif // BISLIPUTILITIES_H

