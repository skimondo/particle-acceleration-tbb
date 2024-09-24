add_library(pngpp INTERFACE)
target_include_directories(pngpp INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>")

