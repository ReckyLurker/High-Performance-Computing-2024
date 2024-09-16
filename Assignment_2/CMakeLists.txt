cmake_minimum_required(VERSION 3.10)
project(HPSC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories("${CMAKE_SOURCE_DIR}/matrix/include")

# Main source file
set(EXTERN_SOURCES
    "${CMAKE_SOURCE_DIR}/main.cpp"
)

# Find all .cpp files in matrix/src directory
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/matrix/src/*.cpp")

# Add external sources to the list of sources
list(APPEND SOURCES ${EXTERN_SOURCES})

# Create the executable
add_executable(HPSC ${SOURCES})

# Include header directories
target_include_directories(HPSC PRIVATE
    "${CMAKE_SOURCE_DIR}/matrix/include"
)

# Specify additional compiler flags if needed
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
