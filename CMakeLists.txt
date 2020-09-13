cmake_minimum_required(VERSION 3.0)

project(AuctionAlgorithm)

add_subdirectory(test)

include_directories(./external/eigen/ ./include/)

add_executable(${PROJECT_NAME} LSAP.cpp)
target_link_libraries(${PROJECT_NAME} pthread)

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
endif()
