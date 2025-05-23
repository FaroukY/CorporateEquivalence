cmake_minimum_required(VERSION 3.10)
project(Equivalence)

set(CMAKE_CXX_STANDARD 17)

# Base optimization & CPU tuning flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -march=native -fopenmp")

# ABI compatibility
add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)

# Enable AddressSanitizer globally (compile & link flags)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -march=native -fopenmp")

# Generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")

# Add subdirectories
add_subdirectory(dsg)
add_subdirectory(dss)
add_subdirectory(sfm)
add_subdirectory(anchored)
add_subdirectory(usss)
