cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

project(BD_BWT_index)

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

file(GLOB LIBRARY_SOURCES dbwt/*.c src/io_tools.cpp)
file(GLOB TEST_SOURCES src/main.cpp)
#file(GLOB SOURCES src/*.cpp dbwt/*.c)

include_directories(include dbwt ../sdsl-lite/include ../sdsl-lite/build/external/libdivsufsort/include)

set(PROJECT_LINK_LIBS libsdsl.a libdivsufsort64.a libbdbwt.a)
link_directories(lib ../sdsl-lite/build/lib ../sdsl-lite/build/external/libdivsufsort/lib)

add_executable(tests ${TEST_SOURCES} ${COMMON_SOURCES})
target_link_libraries(tests ${PROJECT_LINK_LIBS})

add_library(bdbwt ${LIBRARY_SOURCES} ${COMMON_SOURCES})

# Copy headers to include
add_custom_command(TARGET bdbwt PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/../sdsl-lite/include ${CMAKE_SOURCE_DIR}/include)
add_custom_command(TARGET bdbwt PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/../sdsl-lite/build/external/libdivsufsort/include ${CMAKE_SOURCE_DIR}/include)
