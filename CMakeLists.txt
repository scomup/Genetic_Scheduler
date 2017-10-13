cmake_minimum_required(VERSION 2.8.12)  # Ships with Ubuntu 14.04 (Trusty)

project(Scheduler)

#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -Wall  -O0 -march=native ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall   -O0 -march=native")

set(CMAKE_BUILD_TYPE Relase)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11 -Wall  -O3 -march=native ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall   -O3 -march=native")

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

  
find_package(Boost REQUIRED COMPONENTS iostreams)



include_directories(${Boost_INCLUDE_DIRS})

file(GLOB_RECURSE ALL_SRCS "*.cc")

#Message(FATAL_ERROR ${ALL_SRCS})
add_executable(${PROJECT_NAME}  ${ALL_SRCS})


target_link_libraries(${PROJECT_NAME}
${Boost_LIBRARIES}
)
