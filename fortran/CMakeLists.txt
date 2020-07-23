cmake_minimum_required(VERSION 3.12)

project(maooam Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-O2 -Wall -cpp")
    set(debug "-g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow -fsanitize=address -cpp")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-assume byterecl -O2 -fpp")
    set(debug "-O0 -check all -traceback -fpe0 -check bounds -debug all -check uninit -ftrapuv -assume byterecl -fpp")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
find_package(LAPACK REQUIRED)

# Copy config files
file(GLOB nml_files  ./*.nml)
file(COPY ${nml_files} DESTINATION .)

#
# Compile.
#
file(GLOB sources  ./*.f90)
list(FILTER sources EXCLUDE REGEX test)
add_executable(maooam ${sources})
target_link_libraries(maooam ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

