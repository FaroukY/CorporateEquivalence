include_directories(${CMAKE_SOURCE_DIR}/src/external/exactDSP-cpp)
add_library(hi_pr_obj OBJECT
  src/external/exactDSP-cpp/hi_pr.c
)
target_compile_definitions(hi_pr_obj PRIVATE
  NO_MAIN       # don't include main()
  NDEBUG        # disable debug checks
  CUT_ONLY      # min-cut only
#   PRINT_CUT     # print sink size
)
# Match original optimization level and warnings
set_target_properties(hi_pr_obj PROPERTIES
  COMPILE_FLAGS "-O4 -Wall -march=native -flto"
)

add_library(incremental_lib src//incremental.cpp)
target_include_directories(incremental_lib PUBLIC src/external/incremental)
target_compile_definitions(incremental_lib PRIVATE FAST_IO DIRECTED_CASE)
target_compile_definitions(incremental_lib PRIVATE ENABLE_LOGGING=1)

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

# Enable OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Enable TBB (for parallel execution policies)
find_package(TBB REQUIRED)


# Add executable
add_executable(dsg_driver src/main.cpp src/fista.cpp src/graph.cpp src/parser.cpp src/mnp.cpp src/frankwolfe.cpp src/pushrelabel.cpp src/pushrelabel_contra.cpp src/rcdm_utilities.h src/rcdm_utilities.cpp src/rcdm_permutation.h src/rcdm_permutation.cpp src/pushrelabel_contra.h src/stats.cpp src/greedypp.cpp src/stats.h src/logger.h $<TARGET_OBJECTS:hi_pr_obj>)
target_compile_definitions(dsg_driver PRIVATE ENABLE_LOGGING=1)
target_compile_options(dsg_driver PRIVATE -O4 -Wall -march=native -flto -fopenmp)


# Link OpenMP
target_link_libraries(dsg_driver
    PRIVATE
    OpenMP::OpenMP_CXX
    TBB::tbb
    incremental_lib
)

# And Eigen if found installed
if (TARGET Eigen3::Eigen)
    target_link_libraries(dsg_driver PRIVATE Eigen3::Eigen)
endif()
