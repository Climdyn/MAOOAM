
ifndef COMPILER
  COMPILER=gfortran
endif
ifdef INTMODE
  # Define integrator module
  INTEGRATOR = $(INTMODE)_integrator
  TL_AD_INTEGRATOR = $(INTMODE)_tl_ad_integrator
else
  # Define integrator module
  INTEGRATOR = rk2_integrator
  TL_AD_INTEGRATOR = rk2_tl_ad_integrator
endif
ifdef RES
  ST=_$(RES)
  $(warning Code will $(RES) the inner products in arrays! It may induce a huge memory usage!)
endif
# Define relevant compiler options for debug case and normal case for each compiler
# Normal Case:
IFORT_COMPILER_FLAGS = -assume byterecl -O2
IFORT_LD_FLAGS = -mkl
GCC_COMPILER_FLAGS = -O2 -Wall
GCC_LD_FLAGS = -llapack -lblas
# Debug Case:
IFORT_COMPILER_FLAGS_DBG = -O0 -check all -traceback -fpe0 -check bounds -debug all -check uninit -ftrapuv -assume byterecl
IFORT_LD_FLAGS_DBG = -llapack -lblas
GCC_COMPILER_FLAGS_DBG = -g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow 
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
  IFPORT = ifport.mod
else ifeq ($(COMPILER),ifort)
  $(warning This code was optimized for ifort >= 14.0.2 and has not been tested with other compilers)
  COMPILER_FLAGS = $(IFORT_COMPILER_FLAGS)
  LD_FLAGS = $(IFORT_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(IFORT_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(IFORT_LD_FLAGS_DBG)
endif

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = maooam test_tl_ad maooam_stoch  

PROGRAMS_MTV = maooam_MTV

PROGRAMS_WL = maooam_WL

MODULES = util.mod params.mod stoch_params.mod stoch_mod.mod inprod_analytic$(ST).mod aotensor_def.mod ic_def.mod $(INTEGRATOR).mod stat.mod tensor.mod tl_ad_tensor.mod sf_def.mod dec_tensor.mod WL_tensor.mod int_comp.mod corrmod.mod int_corr.mod MTV_int_tensor.mod MTV_sigma_tensor.mod corr_tensor.mod rk2_MTV_integrator.mod rk2_stoch_integrator.mod rk2_ss_integrator.mod sqrt_mod.mod MAR.mod rk2_WL_integrator.mod memory.mod $(TL_AD_INTEGRATOR).mod $(IFPORT)
MODULE_OBJECTS = $(MODULES:.mod=.o)

TEST_PROGRAMS = test_inprod_analytic test_aotensor
TESTS = $(TEST_PROGRAMS:=.test)

TEST_MTV_PROGRAMS = test_dec_tensor test_corr test_MTV_sigma_tensor test_MTV_int_tensor test_MTV_cubic_dampings
TESTS_MTV = $(TEST_MTV_PROGRAMS:=.test)

TEST_WL_PROGRAMS = test_dec_tensor test_corr test_MAR test_WL_tensor test_corr_tensor test_memory
TESTS_WL = $(TEST_WL_PROGRAMS:=.test)

all: $(PROGRAMS) stoch

corr: $(PROGRAMS_CORR)

mtv: corr $(PROGRAMS_MTV) 

wl: corr $(PROGRAMS_WL)

stoch: mtv wl

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

test_mtv: $(TESTS_MTV)

test_dec_tensor.test: test_dec_tensor
	./$< > tests/$<.out

test_MTV_int_tensor.test: test_MTV_int_tensor
	./$< > tests/$<.out

test_MTV_cubic_dampings.test: test_MTV_cubic_dampings
	./$< > tests/$<.out

test_corr.test: test_corr
	./$< > tests/$<.out

test_MTV_sigma_tensor.test: test_MTV_sigma_tensor
	./$< > tests/$<.out

test_wl: $(TESTS_WL)

test_MAR.test: test_MAR
	./$< > tests/$<.out

test_WL_tensor.test: test_WL_tensor
	./$< > tests/$<.out

test_corr_tensor.test: test_corr_tensor
	./$< > tests/$<.out

test_memory.test: test_memory
	./$< > tests/$<.out

%.mod: %.f90
	$(COMPILER) $(COMPILER_FLAGS) -c $< $(LD_FLAGS) 

%.o: %.f90 $(MODULES)
	$(COMPILER) $(COMPILER_FLAGS) -c $< $(LD_FLAGS) 

%: %.o $(MODULE_OBJECTS)
	$(COMPILER) $(COMPILER_FLAGS) $^ -o $@ $(LD_FLAGS) 

tensor.mod: util.mod
inprod_analytic$(ST).mod: params.mod util.mod
aotensor_def.mod: tensor.mod inprod_analytic$(ST).mod
$(INTEGRATOR).mod: aotensor_def.mod
ic_def.mod: util.mod aotensor_def.mod
stat.mod: params.mod
tl_ad_tensor.mod: aotensor_def.mod
dec_tensor.mod: aotensor_def.mod sf_def.mod params.mod stoch_params.mod
WL_tensor.mod: tensor.mod dec_tensor.mod corr_tensor.mod stoch_params.mod
sf_def.mod: util.mod aotensor_def.mod
corrmod.mod: params.mod sf_def.mod stoch_params.mod
stoch_params.mod: params.mod
int_corr.mod: params.mod sf_def.mod tensor.mod int_comp.mod corrmod.mod
MTV_int_tensor.mod: tensor.mod dec_tensor.mod int_corr.mod stoch_params.mod
corr_tensor.mod: tensor.mod params.mod stoch_params.mod corrmod.mod util.mod sf_def.mod
rk2_stoch_integrator.mod: params.mod tensor.mod dec_tensor.mod stoch_params.mod stoch_mod.mod
rk2_ss_integrator.mod: params.mod tensor.mod dec_tensor.mod stoch_params.mod stoch_mod.mod
rk2_MTV_integrator.mod: params.mod tensor.mod MTV_int_tensor.mod MTV_sigma_tensor.mod rk2_ss_integrator.mod stoch_mod.mod
rk2_WL_integrator.mod: params.mod tensor.mod WL_tensor.mod MAR.mod memory.mod rk2_ss_integrator.mod stoch_mod.mod
memory.mod: params.mod stoch_params.mod tensor.mod WL_tensor.mod rk2_ss_integrator.mod
MTV_sigma_tensor.mod: sqrt_mod.mod params.mod MTV_int_tensor.mod tensor.mod util.mod
sqrt_mod.mod: params.mod util.mod
$(TL_AD_INTEGRATOR).mod: tl_ad_tensor.mod
util.mod: $(IFPORT)
stoch_mod.mod: sf_def.mod

clean:
	$(RM) *.o *.mod $(PROGRAMS) $(PROGRAMS_MTV) $(PROGRAMS_CORR) $(PROGRAMS_WL) $(TEST_PROGRAMS) $(TEST_MTV_PROGRAMS) $(TEST_WL_PROGRAMS) tests/*.out tests/*.dat

.PHONY: clean all test test_mtv mtv test_wl wl %.test
