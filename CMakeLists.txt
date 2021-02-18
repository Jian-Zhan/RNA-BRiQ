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

add_library(briq
math/src/KDTree.cpp
geometry/src/Angles.cpp 
geometry/src/localframe.cpp
geometry/src/TransMatrix.cpp
geometry/src/CsMove.cpp
geometry/src/RMSD.cpp
tools/src/StringTool.cpp
tools/src/InputParser.cpp
tools/src/CmdArgs.cpp
model/src/AtomLib.cpp
model/src/AtomProperty.cpp
model/src/BaseDistanceMatrix.cpp
model/src/BasePairComposition.cpp
model/src/BaseRotamer.cpp
model/src/BaseRotamerLib.cpp
model/src/PhoDistanceMatrix.cpp
model/src/PhophateGroup.cpp
model/src/ProteinRep.cpp
model/src/ResName.cpp
model/src/RnaAtomLib.cpp
model/src/RNABaseName.cpp
model/src/RNARibose.cpp
model/src/Sub.cpp
model/src/RNABaseLib.cpp
model/src/BasePair.cpp
model/src/AssignRNASS.cpp
dataio/src/datapaths.cpp
forcefield/src/AtomicEnergyTable.cpp
forcefield/src/BasePairEnergyTable.cpp
forcefield/src/EnergyTable.cpp
forcefield/src/EnergyTable2.cpp
forcefield/src/NeighborConnectionGenerator.cpp
forcefield/src/DMSearcher.cpp
forcefield/src/RiboConnectHashMap.cpp
forcefield/src/RiboseOxygenEnergyTable.cpp
forcefield/src/RiboConnectToPO3.cpp
forcefield/src/SimpleRiboConnect.cpp
forcefield/src/PhoBasicEnergyTable.cpp
forcefield/src/BaseStackingEnergyTable.cpp
forcefield/src/BasePair6DEnergyTable.cpp
forcefield/src/PO3Builder.cpp
para/src/Parameter.cpp
pred/src/BRNode.cpp
pred/src/BRFoldingTree.cpp
pred/src/BRConnection.cpp
pred/src/MoveMutator.cpp
pred/src/BaseMoveLibrary.cpp
pred/src/MCRun.cpp
pred/src/PhoBasicLocal.cpp
pred/src/ThreeBaseMoveLibrary.cpp
pred/src/BackboneModelingTemplate.cpp
pred/src/FragmentLibrary.cpp
)

add_subdirectory (exec)
