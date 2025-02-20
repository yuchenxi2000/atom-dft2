# Specify the minimum version of CMake
cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(atomdft VERSION 1.0 LANGUAGES CXX)

# release/debug
set(CMAKE_BUILD_TYPE "Release")

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Specify the source files
set(SOURCE_DIR "atom-dft/")

set(RAW_SOURCES
    ode.cpp
    linalg.cpp
    interp.cpp
    integrate.cpp
    radial_wave.cpp
    radial_grid.cpp
    poisson.cpp
)

string(REGEX REPLACE "([^;]+)" "${SOURCE_DIR}\\1" SOURCES "${RAW_SOURCES}")

set(RAW_HEADERS
    ode.hpp
    linalg.hpp
    interp.hpp
    integrate.hpp
    radial_wave.hpp
    radial_grid.hpp
    poisson.hpp
)

string(REGEX REPLACE "([^;]+)" "${SOURCE_DIR}\\1" HEADERS "${RAW_HEADERS}")

# Create a shared library
add_library(atomdft SHARED ${SOURCES} ${HEADERS})

# Set properties for the library
set_target_properties(atomdft PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
)
