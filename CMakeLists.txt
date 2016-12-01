########################################################
# cmake file for building LRGRPC
# @author Lagarde François, IPNL
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
########################################################

add_compile_options(-std=c++1y)

project( Proto )
set( ${PROJECT_NAME}_VERSION_MAJOR 4 )
set( ${PROJECT_NAME}_VERSION_MINOR 8 )
set( ${PROJECT_NAME}_VERSION_PATCH 0 )

add_definitions( "-std=c++11 -O5 -Wall -ansi -Wno-long-long -Wno-sign-compare -fno-strict-aliasing -lpthread -ldl -lm" )
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/source/ROOT/)
list( APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake/)

include_directories(${CMAKE_SOURCE_DIR}/source/src/)
include_directories(${CMAKE_SOURCE_DIR}/source/include/)
include_directories(${CMAKE_SOURCE_DIR}/build/include/)

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})

ADD_EXECUTABLE( Proto ./source/Analysis.cpp)
TARGET_LINK_LIBRARIES( Proto ${ROOT_LIBRARIES})
INSTALL( TARGETS Proto DESTINATION ${CMAKE_SOURCE_DIR}/ )
