
! model_def.f90
!
!>  Module to articulate the model classes and define a model version
!
!> @copyright
!> 2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.
!
!---------------------------------------------------------------------------!


MODULE model_def
  USE params
  USE inprod_analytic
  USE aotensor_def
  USE tl_ad_tensor
  USE tensor_def
  USE util, only: str,rstr,init_random_seed
  IMPLICIT NONE

  PRIVATE

  !> Class to hold the components of a model version.
  TYPE, PUBLIC :: Model
    TYPE(ModelConfiguration) :: model_configuration     !< Model configuration object of the model
    TYPE(InnerProducts) :: inner_products               !< Inner products object of the model
    TYPE(AtmOcTensor) :: aotensor                       !< Atmosphere-Ocean tendencies tensor of the model
    TYPE(TlTensor) :: tltensor                          !< Tangent linear model tendencies
    TYPE(AdTensor) :: adtensor                          !< Adjoint model tendencies
    INTEGER, POINTER :: ndim                            !< Dimension of the phase space of the model to integrate.
    LOGICAL :: initialized = .FALSE.
  CONTAINS
    PROCEDURE :: init => init_model
    PROCEDURE :: init_tl => init_tl_model
    PROCEDURE :: init_ad => init_ad_model
    PROCEDURE :: load_ic
    PROCEDURE :: tendencies
    PROCEDURE :: tl_tendencies
    PROCEDURE :: ad_tendencies
    PROCEDURE :: jacobian
    PROCEDURE :: jacobian_mat
    PROCEDURE :: clean => delete_model

  END TYPE Model

CONTAINS
  !-----------------------------------------------------!
  !                                                     !
  ! Jacobian functions                                  !
  !                                                     !
  !-----------------------------------------------------!

  !> Compute the Jacobian of MAOOAM in point \f$\boldsymbol{y}^\ast\f$ and return a tensor object.
  !> @param[in] imodel Model to return the Jacobian of.
  !> @param[in] ystar Vector \f$\boldsymbol{y}^\ast\f$ at which the jacobian should be evaluated.
  !> @return Jacobian in tensor form (table of tuples {i,j,0,value}).
  FUNCTION jacobian(imodel, ystar)
    CLASS(Model), INTENT(IN) :: imodel
    REAL(KIND=8), DIMENSION(0:imodel%ndim), INTENT(IN) :: ystar
    TYPE(Tensor) :: jacobian
    CALL jacobian%init(imodel%ndim)
    CALL imodel%aotensor%tensor%jsparse_mul(ystar,jacobian)
  END FUNCTION jacobian

  !> Compute the Jacobian of MAOOAM in point \f$\boldsymbol{y}^\ast\f$ and return a tensor object.
  !> @param[in] imodel Model to return the Jacobian of.
  !> @param[in] ystar Vector \f$\boldsymbol{y}^\ast\f$ at which the jacobian should be evaluated.
  !> @return Jacobian in matrix form.
  FUNCTION jacobian_mat(imodel, ystar)
    CLASS(Model), INTENT(IN) :: imodel
    REAL(KIND=8), DIMENSION(0:imodel%ndim), INTENT(IN) :: ystar
    REAL(KIND=8), DIMENSION(imodel%ndim,imodel%ndim) :: jacobian_mat
    CALL imodel%aotensor%tensor%jsparse_mul_mat(ystar,jacobian_mat)
  END FUNCTION jacobian_mat

  !> Routine computing the tendencies of the model
  !> @param[in] imodel Model to compute the tendencies of.
  !> @param[in] t Time at which the tendencies have to be computed. Actually not needed for autonomous systems.
  !> @param[in] y Point at which the tendencies have to be computed.
  !> @param[out] res Vector to store the tendencies.
  !> @remark Note that it is NOT safe to pass `y` as a result buffer,
  !> as this operation does multiple passes.
  SUBROUTINE tendencies(imodel, t, y,res)
    CLASS(Model), INTENT(IN) :: imodel
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:imodel%ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(0:imodel%ndim), INTENT(OUT) :: res
    CALL imodel%aotensor%tensor%sparse_mul3(y, y, res)
  END SUBROUTINE tendencies

  !> Tendencies for the AD model of MAOOAM in point \f$\boldsymbol{y}^\ast\f$ for a perturbation \f$\boldsymbol{\delta y}\f$.
  !> @param[in] imodel Model to compute the AD tendencies of.
  !> @param[in] t time
  !> @param[in] ystar Vector \f$\boldsymbol{y}^\ast\f$ (current point in model's trajectory).
  !> @param[in] deltay Vector \f$\boldsymbol{\delta y}\f$, i.e. the perturbation of the variables at time t.
  !> @param[out] res Vector to store the tendencies.
  SUBROUTINE ad_tendencies(imodel,t,ystar,deltay,res)
    CLASS(Model), INTENT(IN) :: imodel
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:imodel%ndim), INTENT(IN) :: ystar,deltay
    REAL(KIND=8), DIMENSION(0:imodel%ndim), INTENT(OUT) :: res
    CALL imodel%adtensor%tensor%sparse_mul3(deltay, ystar, res)
  END SUBROUTINE ad_tendencies

  !> Tendencies for the TL model of MAOOAM in point \f$\boldsymbol{y}^\ast\f$ for a perturbation \f$\boldsymbol{\delta y}\f$.
  !> @param[in] imodel Model to compute the TL tendencies of.
  !> @param[in] t time
  !> @param[in] ystar Vector \f$\boldsymbol{y}^\ast\f$ (current point in model's trajectory).
  !> @param[in] deltay Vector \f$\boldsymbol{\delta y}\f$, i.e. the perturbation of the variables at time t.
  !> @param[out] res Vector to store the tendencies.
  SUBROUTINE tl_tendencies(imodel,t,ystar,deltay,res)
    CLASS(Model), INTENT(IN) :: imodel
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:imodel%ndim), INTENT(IN) :: ystar,deltay
    REAL(KIND=8), DIMENSION(0:imodel%ndim), INTENT(OUT) :: res
    CALL imodel%tltensor%tensor%sparse_mul3(deltay,ystar,res)
  END SUBROUTINE tl_tendencies

  !> Subroutine to initialize the model object from NML files.
  !> @param[in,out] imodel Model object to initialize.
  !> @param[in] physics_nml Physical parameters namelist filename
  !> @param[in] mode_nml Modes configuration namelist filename
  !> @param[in] int_nml Numerical integration parameters namelist filename
  !> @remark If no NML filenames are provided, it will assume that the standard filenames of the model NML have to be used (e.g. "params.nml", "modeselection.nml" and "int_params.nml").
  SUBROUTINE init_model(imodel, physics_nml, mode_nml, int_nml)
    CLASS(Model), INTENT(INOUT), TARGET :: imodel
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: physics_nml
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: mode_nml
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: int_nml
    LOGICAL :: ok

    CALL imodel%model_configuration%clean
    CALL imodel%inner_products%clean
    CALL imodel%aotensor%clean

    CALL imodel%model_configuration%init(physics_nml=physics_nml, mode_nml=mode_nml, int_nml=int_nml)
    CALL imodel%inner_products%init(imodel%model_configuration)
    CALL imodel%aotensor%init(imodel%model_configuration, imodel%inner_products)

    imodel%ndim => imodel%model_configuration%modes%ndim

    ok = imodel%model_configuration%initialized .AND. imodel%inner_products%initialized .AND. imodel%aotensor%initialized

    if (ok) imodel%initialized = .TRUE.

  END SUBROUTINE init_model

  !> Subroutine to clean a model object.
  !> @param[in,out] imodel Model object to clean.
  SUBROUTINE delete_model(imodel)
    CLASS(Model), INTENT(INOUT) :: imodel

    NULLIFY(imodel%ndim)

    CALL imodel%model_configuration%clean
    CALL imodel%inner_products%clean
    CALL imodel%aotensor%clean

    CALL imodel%tltensor%clean
    CALL imodel%adtensor%clean

    imodel%initialized = .FALSE.

  END SUBROUTINE delete_model

  !> Subroutine to initialize the TL tendencies tensor of a model.
  !> @param[in,out] imodel Model object to initialize.
  SUBROUTINE init_tl_model(imodel)
    CLASS(Model), INTENT(INOUT), TARGET :: imodel

    CALL imodel%tltensor%clean

    IF (.NOT.imodel%initialized) THEN
      PRINT*, "*** init_tl_model: Trying to initialize TL model of an uninitialized model ! ***"
      PRINT*, "Please first initialize the model before trying to initialize the TL model."
      PRINT*, "Aborting operation."
      RETURN
    END IF
    CALL imodel%tltensor%init(imodel%aotensor)
  END SUBROUTINE init_tl_model

  !> Subroutine to initialize the AD tendencies tensor of a model.
  !> @param[in,out] imodel Model object to initialize.
  SUBROUTINE init_ad_model(imodel)
    CLASS(Model), INTENT(INOUT), TARGET :: imodel

    CALL imodel%adtensor%clean

    IF (.NOT.imodel%initialized) THEN
      PRINT*, "*** init_ad_model: Trying to initialize AD model of an uninitialized model ! ***"
      PRINT*, "Please first initialize the model before trying to initialize the AD model."
      PRINT*, "Aborting operation."
      RETURN
    END IF
    CALL imodel%adtensor%init(imodel%aotensor)
  END SUBROUTINE init_ad_model

  !> Subroutine to initialize the AD tendencies tensor of a model.
  !> @param[in] imodel Model object for wich to load the initial condition.
  !> @param[in] filename Filename of the initial condition NML file.
  !> @return A vector with the initial condition.
  FUNCTION load_ic(imodel, filename) RESULT(IC)
    CLASS(Model), INTENT(IN), TARGET :: imodel
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
    REAL(KIND=8), DIMENSION(0:imodel%ndim) :: IC

    INTEGER :: i,AllocStat,j
    INTEGER, POINTER :: ndim, natm, noc
    CHARACTER(len=20) :: fm
    REAL(KIND=8) :: size_of_random_noise
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    CHARACTER(LEN=4) :: init_type
    LOGICAL :: exists, std
    NAMELIST /IClist/ IC
    NAMELIST /RAND/ init_type,size_of_random_noise,seed

    fm(1:6)='(F3.1)'

    IF (.NOT.imodel%initialized) THEN
      PRINT*, 'Model not yet initialized, impossible to load any initial condition!'
      RETURN
    END IF

    CALL random_seed(size=j)

    ndim => imodel%model_configuration%modes%ndim
    natm => imodel%model_configuration%modes%natm
    noc => imodel%model_configuration%modes%noc

    ALLOCATE(seed(j), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** load_ic: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF


    IF (present(filename)) THEN
      INQUIRE(FILE=filename,EXIST=exists)
      std = .FALSE.
    ELSE
      PRINT*, "Warning: IC filename not provided."
      PRINT*, "Trying to load the standard file IC.nml instead ..."
      INQUIRE(FILE='./IC.nml',EXIST=exists)
      std = .TRUE.
    END IF

    IF (exists) THEN
      IF (std) THEN
        OPEN(8, file="IC.nml", status='OLD', recl=80, delim='APOSTROPHE')
      ELSE
        OPEN(8, file=filename, status='OLD', recl=80, delim='APOSTROPHE')
      END IF
      READ(8,nml=IClist)
      READ(8,nml=RAND)
      CLOSE(8)
      SELECT CASE (init_type)
        CASE ('seed')
          CALL random_seed(put=seed)
          CALL random_number(IC)
          IC=2*(IC-0.5)
          IC=IC*size_of_random_noise*10.D0
          IC(0)=1.0d0
          WRITE(6,*) "*** Namelist file written. Starting with 'seeded' random initial condition !***"
        CASE ('rand')
          CALL init_random_seed()
          CALL random_seed(get=seed)
          CALL random_number(IC)
          IC=2*(IC-0.5)
          IC=IC*size_of_random_noise*10.D0
          IC(0)=1.0d0
          WRITE(6,*) "*** Namelist file written. Starting with random initial condition !***"
        CASE ('zero')
          CALL init_random_seed()
          CALL random_seed(get=seed)
          IC=0
          IC(0)=1.0d0
          WRITE(6,*) "*** Namelist file written. Starting with initial condition in IC.nml !***"
        CASE ('read')
          CALL init_random_seed()
          CALL random_seed(get=seed)
          IC(0)=1.0d0
          ! except IC(0), nothing has to be done IC has already the right values
          WRITE(6,*) "*** Namelist file written. Starting with initial condition in IC.nml !***"
      END SELECT
    ELSE
      CALL init_random_seed()
      CALL random_seed(get=seed)
      IC=0
      IC(0)=1.0D0
      init_type="zero"
      size_of_random_noise=0.D0
      WRITE(6,*) "*** Namelist file written. Starting with 0 as initial condition !***"
    END IF
    IF (std) THEN
      OPEN(8, file="IC.nml", status='REPLACE')
    ELSE
      OPEN(8, file=filename, status='REPLACE')
    END IF
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! Namelist file :                                                              !"
    WRITE(8,'(a)') "! Initial condition.                                                           !"
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,*) ""
    WRITE(8,'(a)') "&ICLIST"
    WRITE(8,*) " ! psi variables"
    DO i=1,natm
      WRITE(8,*) " IC("//TRIM(str(i))//") = ",IC(i),"   ! typ= "&
        &//imodel%inner_products%awavenum(i)%typ//", Nx= "//TRIM(rstr(imodel%inner_products%awavenum(i)&
        &%Nx,fm))//", Ny= "//TRIM(rstr(imodel%inner_products%awavenum(i)%Ny,fm))
    END DO
    WRITE(8,*) " ! theta variables"
    DO i=1,natm
      WRITE(8,*) " IC("//TRIM(str(i+natm))//") = ",IC(i+natm),"   ! typ= "&
        &//imodel%inner_products%awavenum(i)%typ//", Nx= "//TRIM(rstr(imodel%inner_products%awavenum(i)&
        &%Nx,fm))//", Ny= "//TRIM(rstr(imodel%inner_products%awavenum(i)%Ny,fm))
    END DO

    WRITE(8,*) " ! A variables"
    DO i=1,noc
      WRITE(8,*) " IC("//TRIM(str(i+2*natm))//") = ",IC(i+2*natm),"   ! Nx&
            &= "      //TRIM(rstr(imodel%inner_products%owavenum(i)%Nx,fm))//", Ny= "&
        &//TRIM(rstr(imodel%inner_products%owavenum(i)%Ny,fm))
    END DO
    WRITE(8,*) " ! T variables"
    DO i=1,noc
      WRITE(8,*) " IC("//TRIM(str(i+noc+2*natm))//") = ",IC(i+2*natm+noc),"   &
            &! Nx= "      //TRIM(rstr(imodel%inner_products%owavenum(i)%Nx,fm))//", Ny= "&
        &//TRIM(rstr(imodel%inner_products%owavenum(i)%Ny,fm))
    END DO

    WRITE(8,'(a)') "&END"
    WRITE(8,*) ""
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! Initialisation type.                                                         !"
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! type = 'read': use IC above (will generate a new seed);"
    WRITE(8,'(a)') "!        'rand': random state (will generate a new seed);"
    WRITE(8,'(a)') "!        'zero': zero IC (will generate a new seed);"
    WRITE(8,'(a)') "!        'seed': use the seed below (generate the same IC)"
    WRITE(8,*) ""
    WRITE(8,'(a)') "&RAND"
    WRITE(8,'(a)') "  init_type= '"//init_type//"'"
    WRITE(8,'(a,d15.7)') "  size_of_random_noise = ",size_of_random_noise
    DO i=1,j
      WRITE(8,*) " seed("//TRIM(str(i))//") = ",seed(i)
    END DO
    WRITE(8,'(a)') "&END"
    WRITE(8,*) ""
    CLOSE(8)

  END FUNCTION load_IC
END MODULE model_def


