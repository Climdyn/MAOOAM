# Modular arbitrary-order ocean-atmosphere model: MAOOAM -- Fortran implementation #

## About ##

(c) 2013-2020 Lesley De Cruz and Jonathan Demaeyer

See [LICENSE.txt](LICENSE.txt) for license information.

This software is provided as supplementary material with:

* De Cruz, L., Demaeyer, J. and Vannitsem, S.: The Modular Arbitrary-Order
Ocean-Atmosphere Model: MAOOAM v1.0, Geosci. Model Dev., 9, 2793-2808,
[doi:10.5194/gmd-9-2793-2016](http://dx.doi.org/10.5194/gmd-9-2793-2016), 2016.

**Please cite this article if you use (a part of) this software for a
publication.**

The authors would appreciate it if you could also send a reprint of
your paper to <lesley.decruz@meteo.be>, <jonathan.demaeyer@meteo.be> and
<svn@meteo.be>. 

Consult the MAOOAM [code repository](http://www.github.com/Climdyn/MAOOAM)
for updates, and [our website](http://climdyn.meteo.be) for additional
resources.

------------------------------------------------------------------------

## Installation ##

### Using make alone

The program can be installed with Makefile. We provide configuration files for 
two compilers : gfortran and ifort.

By default, gfortran is selected. To select one or the other, simply modify the 
Makefile accordingly or pass the COMPILER flag to `make`.

To install, unpack the archive in a folder or clone with git:

```bash     
git clone https://github.com/Climdyn/MAOOAM.git
cd MAOOAM/fortran
```
and run:
```bash     
make
```
The program `maooam`  must now be available.
     
The command
```bash
make clean
```
removes the compiled files.

For Windows users, a minimalistic GNU development environment
 (including gfortran and make) is available at [www.mingw.org](http://www.mingw.org) .

### Using cmake and make

The program can also be installed with cmake and make. 

To install, first unpack the archive in a folder or clone with git:

```bash     
git clone https://github.com/Climdyn/MAOOAM.git
cd MAOOAM/fortran
```
and create the folder where you want to build the program sources (`build` in the example below):
```bash
mkdir build
cd build
```
and run cmake:
```bash
cmake ..
```
and run make:
```bash     
make
```
The program `maooam` and the configuration files must now be available in the `build` folder.
     
The command
```bash
make clean
```
inside this folder removes the compiled files.

------------------------------------------------------------------------

##  Description of the files ##

The model tendencies are represented through a tensor called atmoctensor which
includes all the coefficients. In the standard implementation using maooam.f90, this tensor is 
computed once at the program initialization.

* maooam.f90 : Main program.
* model_def.f90 : Main model class module.
* aotensor_def.f90 : Tensor class AtmOcTensor module.
* inprod_analytic.f90 : Inner products class module.
* integrator_def.f90 : A module holding the model's integrator base class definition.
* rk2_integrator.f90 : A module which contains the Heun integrator class for the model equations.
* rk2_tl_integrator.f90 : Heun Tangent Linear (TL) model integrator class module.
* rk2_ad_integrator.f90 : Heun Adjoint (AD) model integrator class module.
* rk4_integrator.f90 : A module which contains the RK4 integrator class for the model equations.
* rk4_tl_integrator.f90 : RK4 Tangent Linear (TL) model integrators module.
* rk4_ad_integrator.f90 : Adjoint (AD) model integrators module.
* Makefile : The Makefile.
* CMakeLists.txt : The CMake file.
* params.f90 : The model parameters classes module.
* tl_ad_tensor.f90 : Tangent Linear (TL) and Adjoint (AD) model tensors class definition module.
* test_tl_ad.f90 : Tests for the Tangent Linear (TL) and Adjoint (AD) model versions.
* README.md : A read me file.
* LICENSE.txt : The license text of the program.
* util.f90 : A module with various useful functions.
* tensor_def.f90 : Main tensor class utility module.
* stat.f90 : A module implementing a statistics accumulator class.
* params.nml : A namelist to specify the model parameters.
* int_params.nml : A namelist to specify the integration parameters.
* modeselection.nml : A namelist to specify which spectral decomposition will be used.

A documentation is available [here](./doc/html/index.html) (html) and [here](./doc/latex/Reference_manual.pdf) (pdf).
 
------------------------------------------------------------------------

## Usage ##

The user first has to fill the params.nml and int_params.nml namelist files according to their needs.
Indeed, model and integration parameters can be specified respectively in the params.nml and int_params.nml namelist files. Some examples related to already published article are available in the [params](./params/) folder.

The modeselection.nml namelist can then be filled : 
* NBOC and NBATM specify the number of blocks that will be used in respectively the ocean and
  the atmosphere. Each block corresponds to a given x and y wavenumber.
* The OMS and AMS arrays are integer arrays which specify which wavenumbers of
  the spectral decomposition will be used in respectively the ocean and the
  atmosphere. Their shapes are OMS(NBOC,2) and AMS(NBATM,2).
* The first dimension specifies the number attributed by the user to the block and the second
  dimension specifies the x and the y wavenumbers.
* The VDDG model is given as a default example. It is described in:
    - Vannitsem, S., Demaeyer, J., De Cruz, L., and Ghil, M.: Low-frequency
      variability and heat transport in a loworder nonlinear coupled ocean-atmosphere
      model, Physica D: Nonlinear Phenomena, 309, 71-85, [doi:10.1016/j.physd.2015.07.006](https://doi.org/10.1016/j.physd.2015.07.006), 2015.   
* Note that the variables of the model are numbered according to the chosen
  order of the blocks.

Finally, the IC.nml file specifying the initial condition should be defined. To
obtain an example of this configuration file corresponding to the model you
have previously defined, simply delete the current IC.nml file (if it exists)
and run the program :

    ./maooam

It will generate a new one and start with the 0 initial condition. If you want another 
initial condition, stop the program, fill the newly generated file and restart :

    ./maooam

It will generate two files :
 * evol_field.dat : the recorded time evolution of the variables.
 * mean_field.dat : the mean field (the climatology)
 
By default, the code uses the rk2_integrator class of integrator, which integrates the model with  
the [Heun algorithm](https://en.wikipedia.org/wiki/Heun%27s_method). However, by modifying the file maooam.f90, it is possible to use the 
rk4_integrator class which integrates the model with the [fourth-order Runge-Kutta algorithm (RK4)](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods). It is also possible to write an user-defined integrator by subclassing the base class integrator.

The tangent linear and adjoint models of MAOOAM are provided in the
tl_ad_tensor, with integrators provided in the rk2_tl_integrator, rk2_ad_integrator, rk4_tl_integrator and rk4_ad_integrator modules. 
It is documented [here](./doc/html/md_doc_tl_ad_doc.html).


------------------------------------------------------------------------

## Implementation notes ##

As the system of differential equations is at most bilinear in y[j] (j=1..n), y
being the array of variables, it can be expressed as a tensor contraction
(written using Einstein convention, i.e. indices that occur twice on one side
of an equation are summed over):

    dy  / dt =  T        y   y      (y  == 1)
      i          i,j,k    j   k       0

The tensor T that encodes the differential equations is composed so that:

* T[i][j][k] contains the contribution of dy[i]/dt proportional to y[j]*y[k].
* Furthermore, y[0] is always equal to 1, so that T[i][0][0] is the constant
contribution to var dy[i]/dt.
* T[i][j][0] + T[i][0][j] is the contribution to  dy[i]/dt which is linear in
y[j].

The tensor is composed as an upper triangular matrix 
(in the last two coordinates), and its computation uses the inner products defined in a inprod_analytic module.

The implementation is made using Fortran classes that are linked together.
It turns the model into an instanciated object that can be reused, allowing the usage of several different model versions in the same program.

------------------------------------------------------------------------

## Final Remarks ##

The authors would like to thank Kris for help with the lua2fortran project. It
has greatly reduced the amount of (error-prone) work.

  No animals were harmed during the coding process.
