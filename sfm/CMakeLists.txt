# Add executable
add_executable(sfm_driver src/main.cpp src/supergreedy.h src/supergreedy.cpp src/mnp.cpp src/mnp.h src/frankwolfe.cpp src/frankwolfe.h src/graph.cpp src/parser.cpp src/stats.cpp src/stats.h src/logger.h)
target_compile_definitions(sfm_driver PRIVATE ENABLE_LOGGING=1)

# Link OpenMP
target_link_libraries(sfm_driver
    PRIVATE
    OpenMP::OpenMP_CXX
    TBB::tbb
)

if(USE_FETCHCONTENT_FOR_BOOST)
    target_link_libraries(sfm_driver PRIVATE Boost::heap)
endif()

# And Eigen if found installed
if (TARGET Eigen3::Eigen)
    target_link_libraries(sfm_driver PRIVATE Eigen3::Eigen)
endif()

