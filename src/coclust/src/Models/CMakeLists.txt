
FILE( GLOB models_files *.cpp)
#set(models_files ContinuousLBModel.cpp)
add_library(models ${models_files})

#target_link_libraries(models inparameters)