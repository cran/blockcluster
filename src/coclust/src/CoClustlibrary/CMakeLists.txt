FILE( GLOB library_files *.cpp )

add_library(CoCluster ${library_files})

target_link_libraries(CoCluster data algo models inparameters STKpp ${LIBS_ext})

