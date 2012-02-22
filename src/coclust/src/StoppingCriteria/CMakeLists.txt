FILE( GLOB stopcriteria_files *.h )

add_library(stopcriteria ${stopcriteria_files})

#target_link_libraries(stopcriteria models)