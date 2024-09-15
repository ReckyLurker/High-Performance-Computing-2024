# Install script for directory: /home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lib/libopenblas.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openblas" TYPE FILE FILES "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/openblas_config.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openblas" TYPE FILE FILES "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/generated/f77blas.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openblas" TYPE FILE FILES "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/generated/cblas.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openblas" TYPE FILE FILES
    "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/LAPACKE/example/lapacke_example_aux.h"
    "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/LAPACKE/include/lapack.h"
    "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/LAPACKE/include/lapacke.h"
    "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/LAPACKE/include/lapacke_config.h"
    "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/LAPACKE/include/lapacke_mangling.h"
    "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/src/lapack-netlib/LAPACKE/include/lapacke_utils.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openblas" TYPE FILE FILES "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lapacke_mangling.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/openblas.pc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/OpenBLAS" TYPE FILE FILES "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/OpenBLASConfig.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/OpenBLAS" TYPE FILE RENAME "OpenBLASConfigVersion.cmake" FILES "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/OpenBLASConfigVersion.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/OpenBLAS/OpenBLASTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/OpenBLAS/OpenBLASTargets.cmake"
         "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles/Export/lib/cmake/OpenBLAS/OpenBLASTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/OpenBLAS/OpenBLASTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/OpenBLAS/OpenBLASTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/OpenBLAS" TYPE FILE FILES "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles/Export/lib/cmake/OpenBLAS/OpenBLASTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/OpenBLAS" TYPE FILE FILES "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/CMakeFiles/Export/lib/cmake/OpenBLAS/OpenBLASTargets-noconfig.cmake")
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/interface/cmake_install.cmake")
  include("/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/driver/level2/cmake_install.cmake")
  include("/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/driver/level3/cmake_install.cmake")
  include("/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/driver/others/cmake_install.cmake")
  include("/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/kernel/cmake_install.cmake")
  include("/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lapack/cmake_install.cmake")
  include("/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/utest/cmake_install.cmake")
  include("/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/test/cmake_install.cmake")
  include("/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/lapack-netlib/TESTING/cmake_install.cmake")
  include("/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/ctest/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/recklurker/home/recklurker/codes/High-Performance-Computing-2024/src/libraries/openblas/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
