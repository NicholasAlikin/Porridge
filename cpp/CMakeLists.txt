cmake_minimum_required(VERSION 3.20)

project(a)

set (VERSION_MAJOR 0)
set (VERSION_MINOR 0)
set (VERSION "${VERSION_MAJOR}.${VERSION_MINOR}")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++23 -Wall -O3")
# -fsanitize=address
file (GLOB
            CPP_SOURCES "math/source/*.cpp"
            CPP_SOURCES "fem/source/*.cpp"
            CPP_SOURCES "main.cpp")

include_directories("${PROJECT_SOURCE_DIR}/math/include")
include_directories("${PROJECT_SOURCE_DIR}/fem/include")

include_directories("${PROJECT_BINARY_DIR}")

add_executable(${PROJECT_NAME} ${CPP_SOURCES})