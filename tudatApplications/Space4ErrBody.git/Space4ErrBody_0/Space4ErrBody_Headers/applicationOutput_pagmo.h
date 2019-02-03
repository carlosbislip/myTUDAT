#ifndef TUDAT_PAGMO_APPLICATIONOUTPUT_H
#define TUDAT_PAGMO_APPLICATIONOUTPUT_H

#include <pagmo/problem.hpp>

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Basics/utilities.h"

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

#endif // TUDAT_PAGMO_APPLICATIONOUTPUT_H
