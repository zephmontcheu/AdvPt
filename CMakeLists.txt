cmake_minimum_required( VERSION 3.31 )
project( oktal LANGUAGES CXX )

set( CMAKE_CXX_STANDARD 23)
set( CMAKE_CXX_STANDARD_REQUIRED ON)
set( CMAKE_COMPILE_WARNING_AS_ERROR ON )
add_compile_options(-Wall -Wextra -Wconversion)

#   Create libraries
message(STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")

add_library( oktal )

target_include_directories( oktal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include )

add_subdirectory( src )

#   Pull Helper Library
option( ADVPT_ENABLE_HTGFILE "" ON )
message( STATUS "Downloading AdvPT helper libraries..." )
include( FetchContent )
FetchContent_Declare(
    advpt-helpers
    URL https://teaching.pages.i10git.cs.fau.de/advpt-exercise-ws2526/assets/advpt-helpers.tar.gz
)
FetchContent_MakeAvailable( advpt-helpers )

target_link_libraries(oktal
PUBLIC
    advpt::htgfile)

add_subdirectory( apps )

#   Add test suite
option( OKTAL_BUILD_TESTSUITE "Build the Oktal test suite" ON )

if( OKTAL_BUILD_TESTSUITE )
    include( CTest )
    add_subdirectory( tests )
endif()

