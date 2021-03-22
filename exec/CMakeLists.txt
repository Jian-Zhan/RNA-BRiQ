add_executable(BRiQ_InitPDB GenerateInitModel.cpp)
add_executable(BRiQ_AssignSS TestAssignment.cpp)
add_executable(BRiQ_Predict RunPredict.cpp)
add_executable(BRiQ_Refinement RunRefinement.cpp)

add_dependencies(BRiQ_InitPDB combinelib)
add_dependencies(BRiQ_AssignSS combinelib)
add_dependencies(BRiQ_Predict combinelib)
add_dependencies(BRiQ_Refinement combinelib)

target_link_libraries (BRiQ_InitPDB briq)
target_link_libraries (BRiQ_AssignSS briq)
target_link_libraries (BRiQ_Predict briq)
target_link_libraries (BRiQ_Refinement briq)
