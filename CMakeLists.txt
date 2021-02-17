cmake_minimum_required (VERSION 2.6)
project (BRiQ)

# The version number.
set (BRIQ_VERSION_MAJOR 1)
set (BRIQ_VERSION_MINOR 0)

SET( EXECUTABLE_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/bin" )
SET( LIBRARY_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/lib" )

include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include/")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -L ${LIBRARY_OUTPUT_PATH}" )


if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_subdirectory (math)
add_subdirectory (geometry)
add_subdirectory (tools)
add_subdirectory (model)
add_subdirectory (dataio)
add_subdirectory (forcefield)
add_subdirectory (para)
add_subdirectory (pred)

add_custom_target(combinelib 
        rm -f libbriq.a
	COMMAND ar rcT libbriq.a libnspmath.a libnspgeometry.a
libnspmodel.a libnsptools.a libnspdataio.a libnspforcefield.a libnsppara.a
libnsppred.a
	COMMAND ranlib libbriq.a
	DEPENDS nspmath nspgeometry nsptools nspmodel nspdataio
nspforcefield nsppara nsppred
        WORKING_DIRECTORY ${LIBRARY_OUTPUT_PATH})

add_subdirectory (exec)
