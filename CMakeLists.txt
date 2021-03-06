#
# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their
# definition
# and dump it in the cache along with proper documentation, otherwise set
# CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()
#

# Project Definition
project(CPTimeFit)
cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0046 OLD)

# set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -pg" CACHE STRING "Debug options." FORCE)
# set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "Debug options." FORCE)

# SET(CMAKE_BUILD_TYPE Profiling)
# SET(CMAKE_CXX_FLAGS_PROFILING "-O3 -pg" CACHE STRING "Debug options." FORCE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{DOOMODULESYS}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{DOOFITSYS}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{DOOCORESYS}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{DOOSELECTIONSYS}/cmake/Modules/")

set(BASEPATH "${CMAKE_SOURCE_DIR}/src")
include_directories(${BASEPATH})

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


if(${APPLE})
MESSAGE( STATUS "Building for Mac OS X, switching on C++11 flags for Mac OS X/clang" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
endif(${APPLE})
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
MESSAGE( STATUS "Building for Linux, switching on C++11 flags for Linux/gcc" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -std=c++14 -pedantic -Wall -Wextra")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

FIND_PACKAGE(Boost 1.48 COMPONENTS program_options filesystem thread system random serialization REQUIRED)

FIND_PACKAGE(GSL REQUIRED)

find_package(ROOT REQUIRED)
find_package(DooCore REQUIRED)
find_package(DooSelection REQUIRED)
find_package(DooFit REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_INCLUDE_DIRS})

include_directories(${DooCore_INCLUDE_DIR})
link_directories(${DooCore_LIBRARY_DIR})

include_directories(${DooSelection_INCLUDE_DIR})
link_directories(${DooSelection_LIBRARY_DIR})

include_directories(${DooFit_INCLUDE_DIR})
link_directories(${DooFit_LIBRARY_DIR})

include_directories(${ROOT_INCLUDE_DIR})

file(GLOB EspressoSrcFiles "Espresso/src/*.cpp")
file(GLOB EspressoHeaderDir "Espresso/Espresso")

root_generate_dictionaries_local(espresso_dict Espresso/Espresso/RooGLMFunction.hh)
add_library(Espresso ${EspressoSrcFiles} RooGLMFunction_dict.h RooGLMFunction_dict.cc)

target_include_directories(Espresso PUBLIC ${EspressoHeaderDir})
target_link_libraries (Espresso ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${ROOT_LIBRARIES})

add_dependencies(Espresso espresso_dict)

set(ALL_LIBRARIES ${DooSelection_LIBRARIES} ${DooFit_LIBRARIES} ${DooCore_LIBRARIES} ${ROOT_LIBRARIES} ${ROOFIT_LIBRARIES} ${ADDITIONAL_LIBRARIES} ${Boost_LIBRARIES} "-lTreePlayer")

add_subdirectory(src)

include_directories($ENV{CPTOYMCSYS}/include)
link_directories($ENV{CPTOYMCSYS}/lib)

add_subdirectory(main)
