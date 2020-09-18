cmake_minimum_required(VERSION 3.0)

add_subdirectory(external/googletest)
enable_testing()

include_directories(external/googletest/googletest/include ../include ../external/eigen)

add_executable(AuctionAlgorithmTest AuctionAlgorithmTest.cpp)
target_link_libraries(AuctionAlgorithmTest gtest_main)

add_custom_command(TARGET AuctionAlgorithmTest PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                       ${CMAKE_SOURCE_DIR}/test/data $<TARGET_FILE_DIR:AuctionAlgorithmTest>)

add_executable(LocksTest LocksTest.cpp)
target_link_libraries(LocksTest gtest_main)