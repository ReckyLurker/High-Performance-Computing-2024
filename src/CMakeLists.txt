# CMakeList to automate compliation and linking 

cmake_minimum_required(VERSION 3.10)
project(HPSC)

set(CMAKE_CXX_STANDARD 17) # Armadillo and OpenBLAS requirements 
set(CMAKE_CXX_STANDARD_REQUIRED ON) 

set(LAPACK_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/src/libraries/openblas/src/lapack-netlib/LAPACKE/include")
set(CBLAS_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/src/libraries/openblas/src/lapack-netlib/CBLAS/include")
set(ARMADILLO_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/src/libraries/armadillo/src/include")

include_directories("${CMAKE_SOURCE_DIR}/src/libraries/openblas/build")
include_directories("${CMAKE_SOURCE_DIR}/src/libraries/armadillo/src")

include_directories("${CMAKE_SOURCE_DIR}/matrix")
include_directories("${CMAKE_SOURCE_DIR}/solvers")

find_package(OpenBLAS REQUIRED) # Find libraries built from source

set(SOURCES
    "${CMAKE_SOURCE_DIR}/main.cpp" # Sources execution
)

# Create the executable
add_executable(HPSC ${SOURCES})

# Link the libraries
target_link_libraries(HPSC PRIVATE OpenBLAS::OpenBLAS "${CMAKE_SOURCE_DIR}/libraries/armadillo/src/libarmadillo.so")

# Include Header file directories
target_include_directories(HPSC PRIVATE
    ${ARMADILLO_INCLUDE_DIR}
    ${CBLAS_INCLUDE_DIR}
    ${LAPACK_INCLUDE_DIR}
    "${CMAKE_SOURCE_DIR}/matrix"
    "${CMAKE_SOURCE_DIR}/solvers"
) 

# Specify additional compiler flags if needed
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2") # Some Optimization 