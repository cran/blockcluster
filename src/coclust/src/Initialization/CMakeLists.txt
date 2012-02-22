FILE( GLOB init_files *.h )

add_library(init ${init_files})

#target_link_libraries(init models)