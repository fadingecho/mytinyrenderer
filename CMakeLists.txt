# https://blog.csdn.net/qq_38410730/article/details/102477162
# https://blog.csdn.net/machuanfei_c/article/details/90521997
# https://blog.csdn.net/qm5132/article/details/96116621

cmake_minimum_required(VERSION 2.8)
 
project(mytinyrenderer)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
aux_source_directory(. DIR_SRCS)
 
include_directories(./include)
 
add_subdirectory(src)
 
add_executable(test ${DIR_SRCS})
 
target_link_libraries(test source_lib)