 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the source files.
set(BISLIP_SOURCES
  "${SRCROOT}${BISLIPDIR}/bislipVariables.cpp"
  "${SRCROOT}${BISLIPDIR}/bislipDecisionVectorEvaluation.cpp"
  "${SRCROOT}${BISLIPDIR}/updateGuidance.cpp"

)

# Set the header files.
set(BISLIP_HEADERS
    "${SRCROOT}${BISLIPDIR}/bislipHeaders.h"
    "${SRCROOT}${BISLIPDIR}/bislipUtilities.h"
    "${SRCROOT}${BISLIPDIR}/bislipVariables.h"
    "${SRCROOT}${BISLIPDIR}/bislipParameters.h"
    "${SRCROOT}${BISLIPDIR}/bislipVehicleSystems.h"
    "${SRCROOT}${BISLIPDIR}/bislipProblemInput.h"
    "${SRCROOT}${BISLIPDIR}/updateGuidance.h"

  )
# Find PaGMO library on local system.
if( USE_PAGMO )
    find_package(PaGMO2)

    # Include PaGMO directories.
    if(NOT APPLE OR APPLE_INCLUDE_FORCE)
        include_directories(SYSTEM AFTER "${PAGMO_INCLUDE_DIR}")
    else( )
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${PAGMO_INCLUDE_DIR}\"")
    endif( )

    list(APPEND TUDAT_APPLICATION_EXTERNAL_LIBRARIES pthread)

endif( )

# Add static libraries.
add_library(bislip STATIC ${BISLIP_SOURCES} ${BISLIP_HEADERS})
setup_tudat_library_target(bislip ${SRCROOT}${BISLIPDIR})
