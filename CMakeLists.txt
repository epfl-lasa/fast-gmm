cmake_minimum_required(VERSION 3.9)

project(fast_gmm VERSION 0.0.0)

# Build options
option(BUILD_TESTING "Build all tests." OFF)

# Default to C99
if(NOT CMAKE_C_STANDARD)
 set(CMAKE_C_STANDARD 99)
endif()

# Default to C++17
if(NOT CMAKE_CXX_STANDARD)
 set(CMAKE_CXX_STANDARD 17)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
 add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_library(mathlib REQUIRED)

if(BUILD_TESTING)
 enable_testing()
 find_package(GTest REQUIRED)
 if (APPLE)
  add_definitions(-DGTEST_USE_OWN_TR1_TUPLE)
  add_definitions(-D__GLIBCXX__)
 endif (APPLE)
else()
 find_package(GTest QUIET)
endif()

include_directories(include)

FILE(GLOB CORE_SOURCES
  src/CDDynamics.cpp
  src/Gaussians.cpp
  src/GMRDynamics.cpp
  src/MJDynamics.cpp
  src/SVR.cpp
  src/ThirdPoly.cpp
)

add_library(${PROJECT_NAME} SHARED
  ${CORE_SOURCES}
)

target_link_libraries(${PROJECT_NAME}
  armadillo
)

install(DIRECTORY include/
  DESTINATION include
)

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

if (BUILD_TESTING)
 add_executable(${PROJECT_NAME}_test src/main.cpp)
 target_link_libraries(${PROJECT_NAME}_test
   ${PROJECT_NAME}
   ${GTEST_LIBRARIES}
   mathlib
   pthread
 )
 add_test(NAME ${PROJECT_NAME}_test COMMAND ${PROJECT_NAME}_test)
endif ()
