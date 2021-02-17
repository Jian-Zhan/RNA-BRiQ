cmake_minimum_required (VERSION 2.6)
project (migrated)

# The version number.
set (Dataio_VERSION_MAJOR 0)
set (Dataio_VERSION_MINOR 0)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include/")
SET(ENV{LD_LIBRARY_PATH} "${PROJECT_BINARY_DIR}/lib:/usr/local/lib:/usr/local/lib64:$ENV{LD_LIBRARY_PATH}")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -std=c++11 -L/usr/local/lib/ -pthread  -lloop" )

SET( EXECUTABLE_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/bin" )
SET( LIBRARY_OUTPUT_PATH  "${PROJECT_BINARY_DIR}/lib" )
add_library(loop1)
target_link_libraries(loop1 nspmath nspgeometry nsptools nspmodel
nspdataio nspforcefield nsppara nsppred)
add_dependencies(loop1 combinelib)
add_custom_target(combinelib 
        rm -f libloop.a
	COMMAND ar rcT libloop.a libnspmath.a libnspgeometry.a
libnspmodel.a libnsptools.a libnspdataio.a libnspforcefield.a libnsppara.a
libnsppred.a
	COMMAND ranlib libloop.a
	DEPENDS nspmath nspgeometry nsptools nspmodel nspdataio
nspforcefield nsppara nsppred
        WORKING_DIRECTORY ${LIBRARY_OUTPUT_PATH})
add_subdirectory (math)
add_subdirectory (geometry)
add_subdirectory (tools)
add_subdirectory (model)
add_subdirectory (dataio)
add_subdirectory (forcefield)
add_subdirectory (para)
add_subdirectory (pred)
