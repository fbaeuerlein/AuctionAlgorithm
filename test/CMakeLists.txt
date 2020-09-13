cmake_minimum_required(VERSION 3.0)

add_subdirectory(external/googletest)
enable_testing()

include_directories(external/googletest/googletest/include ../include ../external/eigen)

add_executable(AuctionAlgorithmTest AuctionAlgorithmTest.cpp)
target_link_libraries(AuctionAlgorithmTest gtest_main)