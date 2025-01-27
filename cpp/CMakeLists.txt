cmake_minimum_required(VERSION 3.20)

project(a)

set (VERSION_MAJOR 0)
set (VERSION_MINOR 0)
set (VERSION "${VERSION_MAJOR}.${VERSION_MINOR}")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++23 -Wall -O3 -D BEAM_RELEASES -D CORRECTOR_PRINTITER -g") # -D BEAM_RELEASES
# -fsanitize=address

file (GLOB
            CPP_SOURCES "math/source/*.cpp"
            CPP_SOURCES "fem/source/*.cpp"
            CPP_SOURCES "numerical_path_following/source/*.cpp"
            CPP_SOURCES "harmonic_balance/source/*.cpp"
            CPP_SOURCES "main.cpp")

include_directories("${PROJECT_SOURCE_DIR}/math/include")
include_directories("${PROJECT_SOURCE_DIR}/fem/include")
include_directories("${PROJECT_SOURCE_DIR}/numerical_path_following/include")
include_directories("${PROJECT_SOURCE_DIR}/harmonic_balance/include")

include_directories("${PROJECT_BINARY_DIR}")

add_executable(${PROJECT_NAME} ${CPP_SOURCES})