cmake_minimum_required(VERSION 3.10)
project(USSSProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -m64 -Wall -Wextra -DNDEBUG -g -w")

# Gurobi paths
set(GUROBI_INCLUDE_DIR /sw/external/gurobi/gurobi1003/linux64/include)
set(GUROBI_LIB_DIR /sw/external/gurobi/gurobi1003/linux64/lib)
set(GUROBI_LIBS
    ${GUROBI_LIB_DIR}/libgurobi_c++.a
    ${GUROBI_LIB_DIR}/libgurobi100.so
)
include_directories(${GUROBI_INCLUDE_DIR})

# Define a uniquely named object library for hi_pr to avoid conflict
add_library(usss_hi_pr_obj OBJECT
  ../external/exactDSP-cpp/hi_pr.c
)

target_compile_definitions(usss_hi_pr_obj PRIVATE
    NO_MAIN
    NDEBUG
    CUT_ONLY
    # PRINT_CUT  # Uncomment if needed
)

set_target_properties(usss_hi_pr_obj PROPERTIES
    COMPILE_FLAGS "-O4 -Wall"
)

# Add executable usss_driver and link
add_executable(usss_driver
    src/main.cpp
    src/frankwolfe.cpp
    src/supergreedy.cpp
    src/mnp.cpp
    src/stats.cpp
    src/flow.cpp
    $<TARGET_OBJECTS:usss_hi_pr_obj>
)

target_link_libraries(usss_driver PRIVATE ${GUROBI_LIBS} m)

if(USE_FETCHCONTENT_FOR_BOOST)
    target_link_libraries(usss_driver PRIVATE Boost::heap)
endif()

if (TARGET Eigen3::Eigen)
    target_link_libraries(usss_driver PRIVATE Eigen3::Eigen)
endif()
