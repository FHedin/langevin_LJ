# Copyright (c) 2011-2016, Florent Hédin, Markus Meuwly, and the University of Basel.
#All rights reserved.
#The 3-clause BSD license is applied to this software.
#See LICENSE.txt

cmake_minimum_required (VERSION 2.6)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(TGT "langevin_LJ" CACHE STRING "The executable name.")

set(USE_OMM "ON" CACHE STRING "Enable use of the OpenMM librqry")

project(${TGT} C)


# list of include folders with .h files
include_directories(
include
dSFMT
)

if (USE_OMM)
include_directories(
    /usr/local/openmm/include
)
endif()

# list all source files
set(
SRCS
src/io.c
src/logger.c
src/main.c
src/memory.c
src/parsing.c
src/rand.c
src/tools.c
dSFMT/dSFMT.c
src/ommInterface.c
)

# never remove -DHAVE_SSE2 -DDSFMT_MEXP=19937 as they are necessary for the dSFMT random numbers generator
add_definitions(-DHAVE_SSE2 -DDSFMT_MEXP=19937)
# if you want to disable dSFMT and use the standard C random numbers (NOT RECOMMENDED), comment previous line and uncomment the following :
#add_definitions(-DSTDRAND)

if (USE_OMM)
  link_directories(/usr/local/openmm/lib)
endif()

add_executable(${TGT} ${SRCS})

if ("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    # clang compiler : should perform well on all cpus
#     set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O3 -g")
#     set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O0 -g")

    # clang compiler : optimise for current machine
    set(CMAKE_C_FLAGS_RELEASE "-std=c11 -Wall -Wextra -msse2 -fno-strict-aliasing -O3 -g -march=native")
    set(CMAKE_C_FLAGS_DEBUG   "-std=c11 -Wall -Wextra -msse2 -fno-strict-aliasing -O0 -g")

    # clang
    target_link_libraries(${TGT} m)
    
    if (USE_OMM)
        target_link_libraries(${TGT} m OpenMM)
    endif()

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    # gnu compiler : should perform well on all cpus
#     set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O3 -g")
#     set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -fno-strict-aliasing -O0 -g")

    # gnu compiler : optimise for current machine
    set(CMAKE_C_FLAGS_RELEASE "-std=c11 -Wall -Wextra -pedantic -msse2 -fno-strict-aliasing -O3 -g -march=native")
    set(CMAKE_C_FLAGS_DEBUG   "-std=c11 -Wall -Wextra -pedantic -msse2 -fno-strict-aliasing -O0 -g")

    # gcc
    target_link_libraries(${TGT} m)
    
    if (USE_OMM)
        target_link_libraries(${TGT} m OpenMM)
    endif()
	
elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    # when using intel compiler : optimise for modern intel cpus
#     set(CMAKE_C_FLAGS_RELEASE "-std=c99 -Wall -Wextra -msse2 -ax,CORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSE4.1,SSSE3,SSE3 -fno-strict-aliasing -O3 -g")
#     set(CMAKE_C_FLAGS_DEBUG   "-std=c99 -Wall -Wextra -msse2 -ax,CORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSE4.1,SSSE3,SSE3 -fno-strict-aliasing -O0 -g")

    # when using intel compiler : optimise for current machine
    set(CMAKE_C_FLAGS_RELEASE "-std=c11 -Wall -Wextra -fno-strict-aliasing -O3 -g -xHost")
    set(CMAKE_C_FLAGS_DEBUG   "-std=c11 -Wall -Wextra -fno-strict-aliasing -O0 -g")
    
    # icc
    target_link_libraries(${TGT} m)

    
    if (USE_OMM)
        target_link_libraries(${TGT} m OpenMM)
    endif()
    
endif()
