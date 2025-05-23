cmake_minimum_required(VERSION 3.10)
project(Equivalence)

set(CMAKE_POLICY_VERSION_MINIMUM 3.5) # Necessary flag to fetch Boost
set(CMAKE_CXX_STANDARD 17)

# Base optimization & CPU tuning flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -march=native -fopenmp")

# ABI compatibility
add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)

# Enable AddressSanitizer globally (compile & link flags)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -march=native -fopenmp")

# Enable OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Enable TBB (for parallel execution policies)
find_package(TBB REQUIRED)

# Try to find Eigen
find_package(Eigen3 3.3 NO_MODULE)
if (Eigen3_FOUND)
    message(STATUS "Found Eigen3 via find_package")
    include_directories(${EIGEN3_INCLUDE_DIRS})
else()
    message(WARNING "Eigen3 not found via find_package, falling back to manual include path")
    set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# Toggle this ON to fetch Boost using FetchContent, OFF to use system Boost
set(USE_FETCHCONTENT_FOR_BOOST ON)
if(USE_FETCHCONTENT_FOR_BOOST)
    # Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24+:
	if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
		cmake_policy(SET CMP0135 NEW)
	endif()

    # FetchContent for Boost
    message(STATUS "Fetching Boost using FetchContent. This can take a second; cmake did not hang! (see /build/dev/_deps)")
    set(BOOST_INCLUDE_LIBRARIES heap)
    set(BOOST_ENABLE_CMAKE ON)
    include(FetchContent)
    FetchContent_Declare(
        Boost
        URL https://github.com/boostorg/boost/releases/download/boost-1.81.0/boost-1.81.0.7z
        USES_TERMINAL_DOWNLOAD TRUE
        DOWNLOAD_NO_EXTRACT FALSE
    )
    FetchContent_MakeAvailable(Boost)
else()
  set(BOOST_INCLUDE_DIR /sw/external/boost/boost_1_81_0)
  include_directories(${BOOST_INCLUDE_DIR})
endif()

# Generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")

# Add subdirectories
add_subdirectory(dsg)
add_subdirectory(dss)
add_subdirectory(sfm)
add_subdirectory(anchored)
# add_subdirectory(usss)
