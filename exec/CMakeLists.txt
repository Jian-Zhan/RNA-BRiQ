add_executable(BRiQ_Predict RunPredict.cpp)
add_executable(BRiQ_Refinement RunRefinement.cpp)
add_executable(BRiQ_assignSS TestAssignment.cpp)

add_dependencies(BRiQ_Predict combinelib)
add_dependencies(BRiQ_Refinement combinelib)
add_dependencies(BRiQ_assignSS combinelib)

target_link_libraries (BRiQ_Predict briq)
target_link_libraries (BRiQ_Refinement briq)
target_link_libraries (BRiQ_assignSS briq)
