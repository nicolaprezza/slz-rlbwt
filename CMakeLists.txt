cmake_minimum_required(VERSION 2.6)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (slz-rlbwt)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/internal)
include_directories(~/include) #SDSL headers are here
include_directories(${PROJECT_SOURCE_DIR}/extern/BWTIL)
include_directories(${PROJECT_SOURCE_DIR}/extern/BWTIL/extern/bitvector/include)
include_directories(${PROJECT_SOURCE_DIR}/extern/DYNAMIC/include)
include_directories(${PROJECT_SOURCE_DIR}/extern/DYNAMIC/include/internal)

LINK_DIRECTORIES(~/lib) #SDSL lib are here

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++11")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-g -ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")


add_executable(slz-rlbwt-build slz-rlbwt-build.cpp)
TARGET_LINK_LIBRARIES(slz-rlbwt-build sdsl)
TARGET_LINK_LIBRARIES(slz-rlbwt-build divsufsort)
TARGET_LINK_LIBRARIES(slz-rlbwt-build divsufsort64)

add_executable(slz-rlbwt-locate slz-rlbwt-locate.cpp)
TARGET_LINK_LIBRARIES(slz-rlbwt-locate sdsl)
TARGET_LINK_LIBRARIES(slz-rlbwt-locate divsufsort)
TARGET_LINK_LIBRARIES(slz-rlbwt-locate divsufsort64)

add_executable(slz-rlbwt-count slz-rlbwt-count.cpp)
TARGET_LINK_LIBRARIES(slz-rlbwt-count sdsl)
TARGET_LINK_LIBRARIES(slz-rlbwt-count divsufsort)
TARGET_LINK_LIBRARIES(slz-rlbwt-count divsufsort64)
