cmake_minimum_required(VERSION 3.15)

project(ewaldRG VERSION 1.0 LANGUAGES CXX)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find dependencies
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package (Boost COMPONENTS program_options  REQUIRED)
find_package(OpenMP REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBINT REQUIRED IMPORTED_TARGET libint2)

# Include directory that contains header/include files
include_directories(include ${Boost_INCLUDE_DIRS})

# Add all source files from src directory
file(GLOB SOURCES "src/*.cpp")

add_compile_options(-O3)

add_executable(ewaldRG ${SOURCES})

# Link against libraries
target_link_libraries (ewaldRG PUBLIC Eigen3::Eigen ${Boost_LIBRARIES} PkgConfig::LIBINT OpenMP::OpenMP_CXX)

install(TARGETS ewaldRG DESTINATION bin)

