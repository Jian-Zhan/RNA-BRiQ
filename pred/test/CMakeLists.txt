link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(BRiQ_Predict RunPredict.cpp)
add_executable(BRiQ_Refinement RunRefinement.cpp)
add_executable(BRiQ_assignSS TestAssignment.cpp)

target_link_libraries (BRiQ_Predict loop)
target_link_libraries (BRiQ_Refinement loop)
target_link_libraries (BRiQ_assignSS loop)
