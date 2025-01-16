add_library(pot STATIC
  particle.cpp
  particle.h

  experiments.cpp
  experiments.h

  colormap.cpp
  colormap.h

  potential.cpp
  potential.h

  potentialparallel.cpp
  potentialparallel.h

  optparser.cpp
  optparser.hpp

)
target_include_directories(pot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pot PUBLIC TBB::tbb PNG::PNG pngpp Eigen3::Eigen)

add_executable(potential main.cpp)
target_link_libraries(potential PRIVATE pot)

add_executable(bench_potential bench_potential.cpp)
target_link_libraries(bench_potential PRIVATE pot)
