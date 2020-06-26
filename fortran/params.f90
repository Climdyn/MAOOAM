
! params.f90                                                                
!
!>  The model parameters module. 
!
!> @copyright                                                               
!> 2015-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------
!                                                                           
!                                                                           
!---------------------------------------------------------------------------

MODULE params

  IMPLICIT NONE

  PUBLIC

  !> The subclass containing the physical parameters of the model.
  TYPE, PUBLIC :: PhysicsConfiguration
    REAL(KIND=8) :: n         !< \f$n = 2 L_y / L_x\f$ - Aspect ratio
    REAL(KIND=8) :: phi0      !< Latitude in radian
    REAL(KIND=8) :: rra       !< Earth radius
    REAL(KIND=8) :: sig0      !< \f$\sigma_0\f$ - Non-dimensional static stability of the atmosphere.
    REAL(KIND=8) :: k         !< Bottom atmospheric friction coefficient.
    REAL(KIND=8) :: kp        !< \f$k'\f$ - Internal atmospheric friction coefficient.
    REAL(KIND=8) :: r         !< Frictional coefficient at the bottom of the ocean.
    REAL(KIND=8) :: d         !< Merchanical coupling parameter between the ocean and the atmosphere.
    REAL(KIND=8) :: f0        !< \f$f_0\f$ - Coriolis parameter
    REAL(KIND=8) :: gp        !< \f$g'\f$Reduced gravity
    REAL(KIND=8) :: H         !< Depth of the active water layer of the ocean.
    REAL(KIND=8) :: phi0_npi  !< Latitude exprimed in fraction of pi.

    REAL(KIND=8) :: lambda    !< \f$\lambda\f$ - Sensible + turbulent heat exchange between the ocean and the atmosphere.
    REAL(KIND=8) :: Co        !< \f$C_a\f$ - Constant short-wave radiation of the ocean.
    REAL(KIND=8) :: Go        !< \f$\gamma_o\f$ - Specific heat capacity of the ocean.
    REAL(KIND=8) :: Ca        !< \f$C_a\f$ - Constant short-wave radiation of the atmosphere.
    REAL(KIND=8) :: To0       !< \f$T_o^0\f$ -  Stationary solution for the 0-th order ocean temperature.
    REAL(KIND=8) :: Ta0       !< \f$T_a^0\f$ -  Stationary solution for the 0-th order atmospheric temperature.
    REAL(KIND=8) :: epsa      !< \f$\epsilon_a\f$ - Emissivity coefficient for the grey-body atmosphere.
    REAL(KIND=8) :: Ga        !< \f$\gamma_a\f$ - Specific heat capacity of the atmosphere.
    REAL(KIND=8) :: RR        !< \f$R\f$ - Gas constant of dry air

    REAL(KIND=8) :: scale     !< \f$L_y = L \, \pi\f$ - The characteristic space scale.
    REAL(KIND=8) :: pi        !< \f$\pi\f$
    REAL(KIND=8) :: LR        !< \f$L_R\f$ - Rossby deformation radius
    REAL(KIND=8) :: G         !< \f$\gamma\f$
    REAL(KIND=8) :: rp        !< \f$r'\f$ - Frictional coefficient at the bottom of the ocean.
    REAL(KIND=8) :: dp        !< \f$d'\f$ - Non-dimensional mechanical coupling parameter between the ocean and the atmosphere.
    REAL(KIND=8) :: kd        !< \f$k_d\f$ - Non-dimensional bottom atmospheric friction coefficient.
    REAL(KIND=8) :: kdp       !< \f$k'_d\f$ - Non-dimensional internal atmospheric friction coefficient.

    REAL(KIND=8) :: Cpo       !< \f$C'_a\f$ - Non-dimensional constant short-wave radiation of the ocean.
    REAL(KIND=8) :: Lpo       !< \f$\lambda'_o\f$ - Non-dimensional sensible + turbulent heat exchange from ocean to atmosphere.
    REAL(KIND=8) :: Cpa       !< \f$C'_a\f$ - Non-dimensional constant short-wave radiation of the atmosphere. @remark Cpa acts on psi1-psi3, not on theta.
    REAL(KIND=8) :: Lpa       !< \f$\lambda'_a\f$ - Non-dimensional sensible + turbulent heat exchange from atmosphere to ocean.
    REAL(KIND=8) :: sBpo      !< \f$\sigma'_{B,o}\f$ - Long wave radiation lost by ocean to atmosphere & space.
    REAL(KIND=8) :: sBpa      !< \f$\sigma'_{B,a}\f$ - Long wave radiation from atmosphere absorbed by ocean.
    REAL(KIND=8) :: LSBpo     !< \f$S'_{B,o}\f$ - Long wave radiation from ocean absorbed by atmosphere.
    REAL(KIND=8) :: LSBpa     !< \f$S'_{B,a}\f$ - Long wave radiation lost by atmosphere to space & ocean.
    REAL(KIND=8) :: L         !< \f$L\f$ - Domain length scale
    REAL(KIND=8) :: sc        !< Ratio of surface to atmosphere temperature.
    REAL(KIND=8) :: sB        !< Stefan–Boltzmann constant
    REAL(KIND=8) :: betp      !< \f$\beta'\f$ - Non-dimensional beta parameter

    REAL(KIND=8) :: nua=0.D0  !< Dissipation in the atmosphere
    REAL(KIND=8) :: nuo=0.D0  !< Dissipation in the ocean

    REAL(KIND=8) :: nuap      !< Non-dimensional dissipation in the atmosphere
    REAL(KIND=8) :: nuop      !< Non-dimensional dissipation in the ocean
    LOGICAL :: initialized = .FALSE.
  CONTAINS
    PROCEDURE :: init => init_physics_nml
  END TYPE PhysicsConfiguration

  !> The subclass containing the integration parameters.
  TYPE, PUBLIC :: IntegrationParameters
    REAL(KIND=8) :: t_trans   !< Transient time period
    REAL(KIND=8) :: t_run     !< Effective intergration time (length of the generated trajectory)
    REAL(KIND=8) :: dt        !< Integration time step
    REAL(KIND=8) :: tw        !< Write all variables every tw time units
    REAL(KIND=8) :: tw_snap   !< Write a snapshot every tw_snap time units
    LOGICAL :: writeout       !< Write to file boolean
    LOGICAL :: initialized = .FALSE.
  CONTAINS
    PROCEDURE :: init => init_int_params_nml
  END TYPE IntegrationParameters

  !> The subclass containing the modes parameters.
  TYPE, PUBLIC :: ModesConfiguration
    INTEGER :: nboc   !< Number of atmospheric blocks
    INTEGER :: nbatm  !< Number of oceanic blocks
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: oms   !< Ocean mode selection array
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ams   !< Atmospheric mode selection array
    INTEGER :: natm=0 !< Number of atmospheric basis functions
    INTEGER :: noc=0  !< Number of oceanic basis functions
    INTEGER :: ndim   !< Number of variables (dimension of the model)
    LOGICAL :: initialized = .FALSE.
  CONTAINS
    PROCEDURE :: init => init_modeselection_nml
    PROCEDURE :: clean => delete_modes_config
  END TYPE ModesConfiguration

  !> The general class holding the model configuration.
  TYPE, PUBLIC :: ModelConfiguration
    TYPE(PhysicsConfiguration) :: physics
    TYPE(ModesConfiguration) :: modes
    TYPE(IntegrationParameters) :: integration
    LOGICAL :: initialized = .FALSE.
  CONTAINS
    PROCEDURE :: init => init_model_config
    PROCEDURE :: clean => clean_model_config
  END TYPE ModelConfiguration


CONTAINS

  !> Subroutine to initialize the model configuration with NML files. Reads the physical parameters and mode selection from the namelist.
  !> @param[in,out] model_config Model configuration object
  !> @param[in] physics_nml Physical parameters namelist filename
  !> @param[in] mode_nml Modes configuration namelist filename
  !> @param[in] int_nml Numerical integration parameters namelist filename
  !> @remark If no NML filenames are provided, it will assume that the standard filenames of the model NML have to be used (e.g. "params.nml", "modeselection.nml" and "int_params.nml").
  SUBROUTINE init_model_config(model_config, physics_nml, mode_nml, int_nml)
    CLASS(ModelConfiguration), INTENT(INOUT) :: model_config
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: physics_nml
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: mode_nml
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: int_nml

    IF (present(physics_nml)) THEN
      CALL model_config%physics%init(physics_nml)
    ELSE
      CALL model_config%physics%init
    END IF

    IF (present(mode_nml)) THEN
      CALL model_config%modes%init(mode_nml)
    ELSE
      CALL model_config%modes%init
    END IF

    IF (present(int_nml)) THEN
      CALL model_config%integration%init(int_nml)
    ELSE
      CALL model_config%integration%init
    END IF

    IF ((model_config%physics%initialized).AND.(model_config%modes%initialized).AND.(model_config%integration%initialized)) THEN
      model_config%initialized = .TRUE.
    END IF

  END SUBROUTINE init_model_config

  !> Subroutine to clean the model configuraion object.
  !> @param[in,out] model_config Model configuration object
  SUBROUTINE clean_model_config(model_config)
    CLASS(ModelConfiguration), INTENT(INOUT) :: model_config

    CALL model_config%modes%clean
    model_config%initialized = .FALSE.
    model_config%integration%initialized = .FALSE.
    model_config%physics%initialized = .FALSE.

  END SUBROUTINE clean_model_config

  ! Subroutine to initialize the physical parameters subclass.
  SUBROUTINE init_physics_nml(physics_config, filename)
    CLASS(PhysicsConfiguration), INTENT(INOUT), TARGET :: physics_config
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
    CHARACTER(LEN=*), PARAMETER :: fs='params.nml'

    TYPE(PhysicsConfiguration), POINTER :: phys

    REAL(KIND=8) :: n
    REAL(KIND=8) :: rra
    REAL(KIND=8) :: sig0
    REAL(KIND=8) :: k
    REAL(KIND=8) :: kp
    REAL(KIND=8) :: r
    REAL(KIND=8) :: d
    REAL(KIND=8) :: f0
    REAL(KIND=8) :: gp
    REAL(KIND=8) :: H
    REAL(KIND=8) :: phi0_npi

    REAL(KIND=8) :: lambda
    REAL(KIND=8) :: Co
    REAL(KIND=8) :: Go
    REAL(KIND=8) :: Ca
    REAL(KIND=8) :: To0
    REAL(KIND=8) :: Ta0
    REAL(KIND=8) :: epsa
    REAL(KIND=8) :: Ga
    REAL(KIND=8) :: RR

    REAL(KIND=8) :: scale
    REAL(KIND=8) :: sc
    REAL(KIND=8) :: sB

    REAL(KIND=8) :: nua=0.D0
    REAL(KIND=8) :: nuo=0.D0

    NAMELIST /aoscale/  scale,f0,n,rra,phi0_npi
    NAMELIST /oparams/  gp,r,H,d,nuo
    NAMELIST /aparams/  k,kp,sig0,nua
    NAMELIST /toparams/ Go,Co,To0
    NAMELIST /taparams/ Ga,Ca,epsa,Ta0
    NAMELIST /otparams/ sc,lambda,RR,sB

    IF (present(filename)) THEN
      OPEN(8, file=filename, status='OLD', recl=80, delim='APOSTROPHE')
    ELSE
      OPEN(8, file=fs, status='OLD', recl=80, delim='APOSTROPHE')
    END IF

    READ(8,nml=aoscale)
    READ(8,nml=oparams)
    READ(8,nml=aparams)
    READ(8,nml=toparams)
    READ(8,nml=taparams)
    READ(8,nml=otparams)
    CLOSE(8)

    phys => physics_config

    phys%n = n
    phys%rra = rra
    phys%sig0 = sig0
    phys%k = k
    phys%kp = kp
    phys%r = r
    phys%d = d
    phys%f0 = f0
    phys%gp = gp
    phys%H = H
    phys%phi0_npi = phi0_npi

    phys%lambda = lambda
    phys%Co = Co
    phys%Go = Go
    phys%Ca = Ca
    phys%To0 = To0
    phys%Ta0 = Ta0
    phys%epsa = epsa
    phys%Ga = Ga
    phys%RR = RR

    phys%scale = scale

    phys%sc = sc
    phys%sB = sB
    
    phys%nua = nua
    phys%nuo = nuo

    ! Some general parameters (Domain, beta, gamma, coupling) !
    !---------------------------------------------------------!

    phys%pi=dacos(-1.D0)
    phys%L=phys%scale/phys%pi
    phys%phi0=phys%phi0_npi*phys%pi
    phys%LR=sqrt(phys%gp*phys%H)/phys%f0
    phys%G=-phys%L**2/phys%LR**2
    phys%betp=phys%L/phys%rra*cos(phys%phi0)/sin(phys%phi0)
    phys%rp=phys%r/phys%f0
    phys%dp=phys%d/phys%f0
    phys%kd=phys%k*2
    phys%kdp=phys%kp

    ! DERIVED QUANTITIES                                  !
    !-----------------------------------------------------!

    phys%Cpo=phys%Co/(phys%Go*phys%f0) * phys%RR/(phys%f0**2*phys%L**2)
    phys%Lpo=phys%lambda/(phys%Go*phys%f0)
    phys%Cpa=phys%Ca/(phys%Ga*phys%f0) * phys%RR/(phys%f0**2*phys%L**2)/2 ! Cpa acts on psi1-psi3, not on theta
    phys%Lpa=phys%lambda/(phys%Ga*phys%f0)
    phys%sBpo=4*phys%sB*phys%To0**3/(phys%Go*phys%f0) ! long wave radiation lost by ocean to atmosphere space
    phys%sBpa=8*phys%epsa*phys%sB*phys%Ta0**3/(phys%Go*phys%f0) ! long wave radiation from atmosphere absorbed by ocean
    phys%LSBpo=2*phys%epsa*phys%sB*phys%To0**3/(phys%Ga*phys%f0) ! long wave radiation from ocean absorbed by atmosphere
    phys%LSBpa=8*phys%epsa*phys%sB*phys%Ta0**3/(phys%Ga*phys%f0) ! long wave radiation lost by atmosphere to space & ocea
    phys%nuap=phys%nua/(phys%f0*phys%L**2)
    phys%nuop=phys%nuo/(phys%f0*phys%L**2)

    phys%initialized = .TRUE.

  END SUBROUTINE init_physics_nml

  ! Subroutine to initialize the modes parameters subclass.
  SUBROUTINE init_modeselection_nml(mode_selection, filename)
    CLASS(ModesConfiguration), INTENT(INOUT) :: mode_selection
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
    CHARACTER(LEN=*), PARAMETER :: fs='modeselection.nml'

    INTEGER :: nboc
    INTEGER :: nbatm
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: oms
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ams

    INTEGER :: i
    INTEGER :: AllocStat
    INTEGER, DIMENSION(2) :: s


    NAMELIST /numblocs/ nboc,nbatm
    NAMELIST /modeselection/ oms,ams

    IF (present(filename)) THEN
      OPEN(8, file=filename, status='OLD', recl=80, delim='APOSTROPHE')
    ELSE
      OPEN(8, file=fs, status='OLD', recl=80, delim='APOSTROPHE')
    END IF
    READ(8,nml=numblocs)

    ALLOCATE(oms(nboc,2),ams(nbatm,2), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** init_modeselection_nml: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF

    READ(8,nml=modeselection)
    CLOSE(8)

    mode_selection%nboc = nboc
    mode_selection%nbatm = nbatm

    IF (allocated(mode_selection%oms)) DEALLOCATE(mode_selection%oms)
    ALLOCATE(mode_selection%oms(nboc,2), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** init_modeselection_nml: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF

    IF (allocated(mode_selection%ams)) DEALLOCATE(mode_selection%ams)
    ALLOCATE(mode_selection%ams(nbatm,2), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** init_modeselection_nml: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF



    mode_selection%oms = oms
    mode_selection%ams = ams

    ! Computation of the dimension of the atmospheric         !
    ! and oceanic components                                  !
    !---------------------------------------------------------!

    mode_selection%natm=0
    DO i=1,nbatm
      IF (ams(i,1)==1) THEN
        mode_selection%natm=mode_selection%natm+3
      ELSE
        mode_selection%natm=mode_selection%natm+2
      ENDIF
    ENDDO
    s=shape(oms)
    mode_selection%noc=s(1)

    mode_selection%ndim=2*mode_selection%natm+2*mode_selection%noc

    DEALLOCATE(ams, oms)
    mode_selection%initialized = .TRUE.

  END SUBROUTINE init_modeselection_nml

  ! Subroutine to delete a modes configuration
  SUBROUTINE delete_modes_config(mode_selection)
    CLASS(ModesConfiguration), INTENT(INOUT) :: mode_selection

    IF (allocated(mode_selection%oms)) DEALLOCATE(mode_selection%oms)
    IF (allocated(mode_selection%ams)) DEALLOCATE(mode_selection%ams)

    mode_selection%natm=0
    mode_selection%noc=0
    mode_selection%ndim=0
    mode_selection%natm=0
    mode_selection%nboc=0
    mode_selection%nbatm=0

    mode_selection%initialized = .FALSE.

  END SUBROUTINE delete_modes_config

  ! Subroutine to initialize the integration parameters subclass.
  SUBROUTINE init_int_params_nml(i_params, filename)
    CLASS(IntegrationParameters), INTENT(INOUT) :: i_params
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename
    CHARACTER(LEN=*), PARAMETER :: fs='int_params.nml'

    REAL(KIND=8) :: t_trans
    REAL(KIND=8) :: t_run
    REAL(KIND=8) :: dt
    REAL(KIND=8) :: tw
    REAL(KIND=8) :: tw_snap
    LOGICAL :: writeout

    NAMELIST /int_params/ t_trans,t_run,dt,tw,tw_snap,writeout

    IF (present(filename)) THEN
      OPEN(8, file=filename, status='OLD', recl=80, delim='APOSTROPHE')
    ELSE
      OPEN(8, file=fs, status='OLD', recl=80, delim='APOSTROPHE')
    END IF

    READ(8,nml=int_params)
    CLOSE(8)

    i_params%t_trans = t_trans
    i_params%t_run = t_run
    i_params%dt = dt
    i_params%tw = tw
    i_params%tw_snap = tw_snap
    i_params%writeout = writeout

    i_params%initialized = .TRUE.

  END SUBROUTINE init_int_params_nml

END MODULE params
