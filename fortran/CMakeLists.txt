cmake_minimum_required(VERSION 3.13)

project(maooam Fortran)
#enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-O2 -Wall")
    set(debug "-g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow -fsanitize=address")
    set(links "-llapack -lblas" )
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-assume byterecl -O2")
    set(debug "-O0 -check all -traceback -fpe0 -check bounds -debug all -check uninit -ftrapuv -assume byterecl")
    set(links "-mkl")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
set(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} ${links}")

# Copy config files
file(GLOB nml_files  ./*.nml)
file(COPY ${nml_files} DESTINATION .)

#
# Compile.
#
file(GLOB sources  ./*.f90)
list(FILTER sources EXCLUDE REGEX test)
add_executable(maooam ${sources})


