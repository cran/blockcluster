FILE( GLOB algo_files *.cpp )

add_library(algo ${algo_files})

#target_link_libraries(algo models stopcriteria)