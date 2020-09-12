cmake_minimum_required(VERSION 3.0)

project(AuctionAlgorithm)

add_subdirectory(test)
aux_source_directory(. SRC_LIST)

include_directories(./external/)

add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} pthread)

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
endif()
