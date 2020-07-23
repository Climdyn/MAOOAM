
ifndef COMPILER
  COMPILER=ifort
endif
# Define relevant compiler options for debug case and normal case for each compiler
# Normal Case:
IFORT_COMPILER_FLAGS = -assume byterecl -O2 -fpp
IFORT_LD_FLAGS = -mkl
GCC_COMPILER_FLAGS = -O2 -Wall -cpp
GCC_LD_FLAGS = -llapack -lblas
# Debug Case:
IFORT_COMPILER_FLAGS_DBG = -O0 -check all -traceback -fpe0 -check bounds -debug all -check uninit -ftrapuv -assume byterecl -fpp
IFORT_LD_FLAGS_DBG = -llapack -lblas
GCC_COMPILER_FLAGS_DBG = -g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow -fsanitize=address -cpp
GCC_LD_FLAGS_DBG = -llapack -lblas


SHELL = /bin/sh
DEBUG = false
#Adapt to host system
ifeq (,$(findstring Windows,$(OS)))
  RM = rm -f
else
  RM = del
endif


ifeq ($(COMPILER),gfortran)
  COMPILER_FLAGS = $(GCC_COMPILER_FLAGS)
  LD_FLAGS = $(GCC_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(GCC_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(GCC_LD_FLAGS_DBG)
else ifeq ($(COMPILER),ifort)
  COMPILER_FLAGS = $(IFORT_COMPILER_FLAGS)
  LD_FLAGS = $(IFORT_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(IFORT_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(IFORT_LD_FLAGS_DBG)
endif

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = maooam test_tl_ad

MODULES = params.mod util.mod inprod_analytic.mod tensor_def.mod aotensor_def.mod tl_ad_tensor.mod model_def.mod integrator_def.mod rk2_integrator.mod rk4_integrator.mod rk2_tl_integrator.mod rk2_ad_integrator.mod rk4_tl_integrator.mod rk4_ad_integrator.mod stat.mod
MODULE_OBJECTS = $(MODULES:.mod=.o)

TEST_PROGRAMS = test_inprod_analytic test_aotensor
TESTS = $(TEST_PROGRAMS:=.test)

all: $(PROGRAMS) 

debug: COMPILER_FLAGS = $(COMPILER_FLAGS_DBG) 
debug: LD_FLAGS = $(LD_FLAGS_DBG)
debug: all

test: COMPILER_FLAGS = $(COMPILER_FLAGS_DBG) 
test: LD_FLAGS = $(LD_FLAGS_DBG)
test: $(TESTS)

test_inprod_analytic.test: test_inprod_analytic
	./$< |sort  > tests/$<.out
	diff tests/$<.out tests/$<.ref && echo PASS: $@ || echo FAIL: $@

test_aotensor.test: test_aotensor 
	./$< |sort  > tests/$<.out
	diff tests/$<.out tests/$<.ref && echo PASS: $@ || echo FAIL: $@

%.mod: %.f90
	$(COMPILER) $(COMPILER_FLAGS) -c $< $(LD_FLAGS) 

%.o: %.f90 $(MODULES)
	$(COMPILER) $(COMPILER_FLAGS) -c $< $(LD_FLAGS) 

%: %.o $(MODULE_OBJECTS)
	$(COMPILER) $(COMPILER_FLAGS) $^ -o $@ $(LD_FLAGS) 

inprod_analytic.mod: params.mod util.mod
aotensor_def.mod: tensor_def.mod inprod_analytic.mod params.mod
integrator_def.mod: model_def.mod
rk2_integrator.mod: integrator_def.mod
tl_ad_tensor.mod: aotensor_def.mod tensor_def.mod
rk2_tl_integrator.mod: rk2_integrator.mod

clean:
	$(RM) *.o *.mod $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out

.PHONY: clean all test %.test
