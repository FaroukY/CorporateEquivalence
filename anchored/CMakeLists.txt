# Enable OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Add boost
# Boost (header-only)
set(BOOST_INCLUDE_DIR /sw/external/boost/boost_1_81_0)

# Include directories
include_directories(${BOOST_INCLUDE_DIR})

# include(FetchContent)
# FetchContent_Declare(
#   Boost
#   GIT_REPOSITORY https://github.com/boostorg/boost.git
#   GIT_TAG boost-1.80.0
#   GIT_SHALLOW TRUE
# )
# FetchContent_MakeAvailable(Boost)

# Try to find Eigen
find_package(Eigen3 3.3 NO_MODULE)

if (Eigen3_FOUND)
    message(STATUS "Found Eigen3 via find_package")
    include_directories(${EIGEN3_INCLUDE_DIRS})
else()
    message(WARNING "Eigen3 not found via find_package, falling back to manual include path")
    set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/external/eigen)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# Enable TBB (for parallel execution policies)
find_package(TBB REQUIRED)

add_library(anchored_hi_pr_obj OBJECT
    src/external/exactDSP-cpp/hi_pr.c
)

target_compile_definitions(anchored_hi_pr_obj PRIVATE
    NO_MAIN
    NDEBUG
    CUT_ONLY
)

set_target_properties(anchored_hi_pr_obj PROPERTIES
    COMPILE_FLAGS "-O4 -Wall"
)


# Add executable
add_executable(anchored_driver src/main.cpp src/supergreedy.h src/supergreedy.cpp src/mnp.h src/mnp.cpp src/frankwolfe.h src/frankwolfe.cpp src/graph.cpp src/parser.cpp src/stats.cpp src/flow.cpp src/flow.h src/stats.h src/logger.h $<TARGET_OBJECTS:anchored_hi_pr_obj>)
target_compile_definitions(anchored_driver PRIVATE ENABLE_LOGGING=1)

# Link OpenMP
target_link_libraries(anchored_driver
    PRIVATE
    OpenMP::OpenMP_CXX
    TBB::tbb
)

# And Eigen if found installed
if (TARGET Eigen3::Eigen)
    target_link_libraries(anchored_driver PRIVATE Eigen3::Eigen)
endif()

