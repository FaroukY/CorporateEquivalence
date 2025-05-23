add_library(anchored_hi_pr_obj OBJECT
  ../external/exactDSP-cpp/hi_pr.c
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

if(USE_FETCHCONTENT_FOR_BOOST)
    target_link_libraries(anchored_driver PRIVATE Boost::heap)
endif()

# And Eigen if found installed
if (TARGET Eigen3::Eigen)
    target_link_libraries(anchored_driver PRIVATE Eigen3::Eigen)
endif()

