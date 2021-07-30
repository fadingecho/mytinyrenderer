# CMakeLists.txt
aux_source_directory(. DIR_SOURCE_SRCS)
add_library(source_lib STATIC ${DIR_SOURCE_SRCS})