cmake_minimum_required (VERSION 3.18)

add_library (chopper_interface INTERFACE)
target_link_libraries (chopper_interface INTERFACE seqan3::seqan3)
target_link_libraries (chopper_interface INTERFACE sharg::sharg)
target_link_libraries (chopper_interface INTERFACE "seqan::hibf")
target_include_directories (chopper_interface INTERFACE ../include)
target_include_directories (chopper_interface INTERFACE ${CHOPPER_SUBMODULES_DIR}/robin-hood-hashing/src/include)
target_include_directories (chopper_interface INTERFACE ${CHOPPER_SUBMODULES_DIR}/simde/simde)
target_link_libraries (chopper_interface INTERFACE xxHash::xxhash)
target_compile_options (chopper_interface INTERFACE "-pedantic" "-Wall" "-Wextra")

add_library (chopper_layout_lib STATIC chopper_layout.cpp)
target_link_libraries (chopper_layout_lib "chopper_interface")

add_library (chopper_lib INTERFACE)
target_link_libraries (chopper_lib INTERFACE "chopper_interface" "chopper_layout_lib")

add_executable (chopper chopper.cpp)
target_link_libraries (chopper "chopper_lib")

add_executable (measure_hyperloglog EXCLUDE_FROM_ALL measure_hyperloglog.cpp)
target_link_libraries (measure_hyperloglog "chopper_interface")
target_compile_options (measure_hyperloglog PRIVATE "-Werror")

if (CHOPPER_INSTALL)
    install (TARGETS chopper RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif ()
