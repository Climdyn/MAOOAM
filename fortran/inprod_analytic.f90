
! inprod_analytic.f90
!
!> Inner products between the truncated set of basis functions for the
!> ocean and atmosphere streamfunction fields.
!
!> @remark
!> These are partly calculated using the analytical expressions from
!> Cehelsky, P., & Tung, K. K. : Theories of multiple equilibria and
!> weather regimes-A critical reexamination. Part II: Baroclinic two-layer
!> models. Journal of the atmospheric sciences, 44(21), 3282-3303, 1987.
!
!
!> @copyright                                                               
!> 2015-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!                                                                           
!---------------------------------------------------------------------------!

MODULE inprod_analytic

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params
  USE util, only: isin,piksrt
  IMPLICIT NONE

  PRIVATE

  !> Atmospheric bloc specification object
  TYPE :: AtmosphericWavenumber
    CHARACTER :: typ=" "
    INTEGER :: M=0,P=0,H=0
    REAL(KIND=8) :: Nx=0.,Ny=0.
  END TYPE AtmosphericWavenumber

  !> Oceanic bloc specification object
  TYPE :: OceanicWavenumber
    INTEGER :: P=0,H=0
    REAL(KIND=8) :: Nx=0.,Ny=0.
  END TYPE OceanicWavenumber

  !> Class holding the atmospheric inner products functions
  TYPE :: AtmosphereInnerProducts
    TYPE(InnerProducts), POINTER :: inner_products   !< Pointer to a global inner products object
  CONTAINS
    PROCEDURE :: a => calculate_a         !< Eigenvalues of the Laplacian (atmospheric)
    PROCEDURE :: b => calculate_b         !< Streamfunction advection terms (atmospheric)
    PROCEDURE :: c => calculate_c_atm     !< Beta term for the atmosphere
    PROCEDURE :: d => calculate_d         !< Forcing of the ocean on the atmosphere
    PROCEDURE :: g => calculate_g         !< Temperature advection terms (atmospheric)
    PROCEDURE :: s => calculate_s         !< Forcing (thermal) of the ocean on the atmosphere
  END TYPE AtmosphereInnerProducts

  !> Class holding the oceanic inner products functions
  TYPE :: OceanInnerProducts
    TYPE(InnerProducts), POINTER :: inner_products   !< Pointer to a global inner products object
  CONTAINS
    PROCEDURE :: K => calculate_K         !< Forcing of the atmosphere on the ocean.
    PROCEDURE :: M => calculate_M         !< Forcing of the ocean fields on the ocean.
    PROCEDURE :: C => calculate_C_oc      !< Streamfunction advection terms (oceanic)
    PROCEDURE :: N => calculate_N         !< Beta term for the ocean
    PROCEDURE :: O => calculate_O         !< Temperature advection term (passive scalar)
    PROCEDURE :: W => calculate_W         !< Short-wave radiative forcing of the ocean
  END TYPE OceanInnerProducts

  !> Global class for the inner products. Contains also the modes informations.
  TYPE, PUBLIC :: InnerProducts

    TYPE(ModelConfiguration), POINTER :: model_config   !< Pointer to a model configuration object.

    !> Atmospheric blocs specification
    TYPE(AtmosphericWavenumber), DIMENSION(:), ALLOCATABLE, PUBLIC :: awavenum
    !> Oceanic blocs specification
    TYPE(OceanicWavenumber), DIMENSION(:), ALLOCATABLE, PUBLIC :: owavenum

    !> Atmospheric tensors
    TYPE(AtmosphereInnerProducts), PUBLIC :: atmos
    !> Oceanic tensors
    TYPE(OceanInnerProducts), PUBLIC :: ocean
    LOGICAL :: initialized = .FALSE.

  CONTAINS
    PROCEDURE :: init => init_inner_products        !< Procedure to initialize the inner products functions based on the modes configuration.
    PROCEDURE :: clean => delete_inner_products     !< Procedure to clean the inner products object.

  END TYPE InnerProducts



  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Definition of the Helper functions from Cehelsky    !
  ! & Tung                                              !
  !                                                     !
  !-----------------------------------------------------!

  ! Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION B1(Pi, Pj, Pk)
    INTEGER :: Pi,Pj,Pk
    B1 = (Pk + Pj) / REAL(Pi)
  END FUNCTION B1

  !  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION B2(Pi, Pj, Pk)
    INTEGER :: Pi,Pj,Pk
    B2 = (Pk - Pj) / REAL(Pi)
  END FUNCTION B2

  !  Integer Dirac delta function
  REAL(KIND=8) FUNCTION delta(r)
    INTEGER :: r
    IF (r==0) THEN
      delta = 1.D0
    ELSE
      delta = 0.D0
    ENDIF
  END FUNCTION delta

  !  "Odd or even" function
  REAL(KIND=8) FUNCTION flambda(r)
    INTEGER :: r
    IF (mod(r,2)==0) THEN
      flambda = 0.D0
    ELSE
      flambda = 1.D0
    ENDIF
  END FUNCTION flambda

  !  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S1(Pj, Pk, Mj, Hk)
    INTEGER :: Pk,Pj,Mj,Hk
    S1 = -((Pk * Mj + Pj * Hk)) / 2.D0
  END FUNCTION S1

  !  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S2(Pj, Pk, Mj, Hk)
    INTEGER :: Pk,Pj,Mj,Hk
    S2 = (Pk * Mj - Pj * Hk) / 2.D0
  END FUNCTION S2

  !  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S3(Pj, Pk, Hj, Hk)
    INTEGER :: Pj,Pk,Hj,Hk
    S3 = (Pk * Hj + Pj * Hk) / 2.D0
  END FUNCTION S3

  !  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S4(Pj, Pk, Hj, Hk)
    INTEGER :: Pj,Pk,Hj,Hk
    S4 = (Pk * Hj - Pj * Hk) / 2.D0
  END FUNCTION S4
 
  !-----------------------------------------------------!
  ! Inner products definition routines                  !
  !--------------------------------------------------------!
  ! 1. Inner products in the equations for the atmosphere  !
  !--------------------------------------------------------!
  
  !> Eigenvalues of the Laplacian (atmospheric)
  !> 
  !> \f$ a_{i,j} = (F_i, \nabla^2 F_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_a(self, i,j)
    CLASS(AtmosphereInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j
    TYPE(AtmosphericWavenumber) :: Ti
    
    calculate_a = 0.D0
    IF (i==j) THEN
      Ti = self%inner_products%awavenum(i)
      calculate_a = -(self%inner_products%model_config%physics%n**2) * Ti%Nx**2 - Ti%Ny**2
    END IF
  END FUNCTION calculate_a

  !> Streamfunction advection terms (atmospheric)
  !> 
  !> \f$ b_{i,j,k} = (F_i, J(F_j, \nabla^2 F_k))\f$ .
  REAL(KIND=8) FUNCTION calculate_b(self,i,j,k)
    CLASS(AtmosphereInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j,k

    calculate_b = self%a(k,k) * self%g(i,j,k)

  END FUNCTION calculate_b

  !> Beta term for the atmosphere
  !> 
  !> \f$ c_{i,j} = (F_i, \partial_x F_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_c_atm(self,i,j)
    CLASS(AtmosphereInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j
    TYPE(AtmosphericWavenumber) :: Ti, Tj

    Ti = self%inner_products%awavenum(i)
    Tj = self%inner_products%awavenum(j)
    calculate_c_atm = 0.D0
    IF ((Ti%typ == "K") .AND. (Tj%typ == "L")) THEN 
      calculate_c_atm = Ti%M * delta(Ti%M - Tj%H) * delta(Ti%P - Tj%P)
    ELSE IF ((Ti%typ == "L") .AND. (Tj%typ == "K")) THEN
      Ti = self%inner_products%awavenum(j)
      Tj = self%inner_products%awavenum(i)
      calculate_c_atm = - Ti%M * delta(Ti%M - Tj%H) * delta(Ti%P - Tj%P)
    END IF
    calculate_c_atm = self%inner_products%model_config%physics%n * calculate_c_atm
  END FUNCTION calculate_c_atm

  !> Forcing of the ocean on the atmosphere.
  !> 
  !> \f$ d_{i,j} = (F_i, \nabla^2 \eta_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_d(self,i,j)
    CLASS(AtmosphereInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j

    calculate_d=self%s(i,j) * self%inner_products%ocean%M(j,j)

  END FUNCTION calculate_d

  !> Temperature advection terms (atmospheric)
  !> 
  !> \f$ g_{i,j,k} = (F_i, J(F_j, F_k))\f$ .
  REAL(KIND=8) FUNCTION calculate_g(self,i,j,k)
    CLASS(AtmosphereInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j,k
    TYPE(AtmosphericWavenumber) :: Ti,Tj,Tk
    REAL(KIND=8) :: val,vb1, vb2, vs1, vs2, vs3, vs4
    INTEGER, DIMENSION(3) :: a,b
    INTEGER, DIMENSION(3,3) :: w
    CHARACTER, DIMENSION(3) :: s
    INTEGER :: par

    Ti = self%inner_products%awavenum(i)
    Tj = self%inner_products%awavenum(j)
    Tk = self%inner_products%awavenum(k)

    a(1)=i
    a(2)=j
    a(3)=k

    val=0.D0

    IF ((Ti%typ == "L") .AND. (Tj%typ == "L") .AND. (Tk%typ == "L")) THEN
       
      CALL piksrt(3,a,par)

      Ti = self%inner_products%awavenum(a(1))
      Tj = self%inner_products%awavenum(a(2))
      Tk = self%inner_products%awavenum(a(3))

      vs3 = S3(Tj%P,Tk%P,Tj%H,Tk%H)
      vs4 = S4(Tj%P,Tk%P,Tj%H,Tk%H)
      val = vs3 * ((delta(Tk%H - Tj%H - Ti%H) - delta(Tk%H &
        &- Tj%H + Ti%H)) * delta(Tk%P + Tj%P - Ti%P) +&
        & delta(Tk%H + Tj%H - Ti%H) * (delta(Tk%P - Tj%P&
        & + Ti%P) - delta(Tk%P - Tj%P - Ti%P))) + vs4 *&
        & ((delta(Tk%H + Tj%H - Ti%H) * delta(Tk%P - Tj&
        &%P - Ti%P)) + (delta(Tk%H - Tj%H + Ti%H) -&
        & delta(Tk%H - Tj%H - Ti%H)) * (delta(Tk%P - Tj&
        &%P - Ti%P) - delta(Tk%P - Tj%P + Ti%P)))
    ELSE

      s(1)=Ti%typ
      s(2)=Tj%typ
      s(3)=Tk%typ

      w(1,:)=isin("A",s)
      w(2,:)=isin("K",s)
      w(3,:)=isin("L",s)

      IF (ANY(w(1,:)/=0) .AND. ANY(w(2,:)/=0) .AND. ANY(w(3,:)/=0)) THEN
        b=w(:,1)
        Ti = self%inner_products%awavenum(a(b(1)))
        Tj = self%inner_products%awavenum(a(b(2)))
        Tk = self%inner_products%awavenum(a(b(3)))
        call piksrt(3,b,par)
        vb1 = B1(Ti%P,Tj%P,Tk%P)
        vb2 = B2(Ti%P,Tj%P,Tk%P)
        val = -2 * sqrt(2.) / self%inner_products%model_config%physics%pi * Tj%M * delta(Tj%M - Tk%H) * flambda(Ti%P + Tj%P + Tk%P)
        IF (val /= 0.D0) val = val * (vb1**2 / (vb1**2 - 1) - vb2**2 / (vb2**2 - 1))
      ELSEIF ((w(2,2)/=0) .AND. (w(2,3)==0) .AND. ANY(w(3,:)/=0)) THEN
        Ti = self%inner_products%awavenum(a(w(2,1)))
        Tj = self%inner_products%awavenum(a(w(2,2)))
        Tk = self%inner_products%awavenum(a(w(3,1)))
        b(1)=w(2,1)
        b(2)=w(2,2)
        b(3)=w(3,1)
        call piksrt(3,b,par)
        vs1 = S1(Tj%P,Tk%P,Tj%M,Tk%H)
        vs2 = S2(Tj%P,Tk%P,Tj%M,Tk%H)
        val = vs1 * (delta(Ti%M - Tk%H - Tj%M) * delta(Ti%P -&
          & Tk%P + Tj%P) - delta(Ti%M- Tk%H - Tj%M) *&
          & delta(Ti%P + Tk%P - Tj%P) + (delta(Tk%H - Tj%M&
          & + Ti%M) + delta(Tk%H - Tj%M - Ti%M)) *&
          & delta(Tk%P + Tj%P - Ti%P)) + vs2 * (delta(Ti%M&
          & - Tk%H - Tj%M) * delta(Ti%P - Tk%P - Tj%P) +&
          & (delta(Tk%H - Tj%M - Ti%M) + delta(Ti%M + Tk%H&
          & - Tj%M)) * (delta(Ti%P - Tk%P + Tj%P) -&
          & delta(Tk%P - Tj%P + Ti%P)))
      ENDIF
    ENDIF
    calculate_g=par*val*self%inner_products%model_config%physics%n
 
  END FUNCTION calculate_g

  !> Forcing (thermal) of the ocean on the atmosphere.
  !> 
  !> \f$ s_{i,j} = (F_i, \eta_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_s(self,i,j)
    CLASS(AtmosphereInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j
    TYPE(AtmosphericWavenumber) :: Ti
    TYPE(OceanicWavenumber) :: Dj
    REAL(KIND=8) :: val
    
    Ti = self%inner_products%awavenum(i)
    Dj = self%inner_products%owavenum(j)
    val=0.D0
    IF (Ti%typ == "A") THEN
      val = flambda(Dj%H) * flambda(Dj%P + Ti%P)
      IF (val /= 0.D0) THEN
        val = val*8*sqrt(2.)*Dj%P/(self%inner_products%model_config%physics%pi**2 * (Dj%P**2 - Ti%P**2) * Dj%H)
      END IF
    ELSEIF (Ti%typ == "K") THEN
      val = flambda(2 * Ti%M + Dj%H) * delta(Dj%P - Ti%P)
      IF (val /= 0.D0) THEN
        val = val*4*Dj%H/(self%inner_products%model_config%physics%pi * (-4 * Ti%M**2 + Dj%H**2))
      END IF
    ELSEIF (Ti%typ == "L") THEN
      val = delta(Dj%P - Ti%P) * delta(2 * Ti%H - Dj%H)
    END IF
    calculate_s=val
    
  END FUNCTION calculate_s

  !--------------------------------------------------------!
  ! 2. Inner products in the equations for the ocean       !
  !--------------------------------------------------------!
  
  !> Forcing of the atmosphere on the ocean.
  !> 
  !> \f$ K_{i,j} = (\eta_i, \nabla^2 F_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_K(self,i,j)
    CLASS(OceanInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j

    calculate_K = self%inner_products%atmos%s(j,i) * self%inner_products%atmos%a(j,j)
  END FUNCTION calculate_K

  !> Forcing of the ocean fields on the ocean.
  !> 
  !> \f$ M_{i,j} = (eta_i, \nabla^2 \eta_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_M(self,i,j)
    CLASS(OceanInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j
    TYPE(OceanicWavenumber) :: Di

    calculate_M=0.D0
    IF (i==j) THEN
      Di = self%inner_products%owavenum(i)
      calculate_M = -(self%inner_products%model_config%physics%n**2) * Di%Nx**2 - Di%Ny**2
    END IF
  END FUNCTION calculate_M

  !> Beta term for the ocean
  !> 
  !> \f$ N_{i,j} = (\eta_i, \partial_x \eta_j) \f$.
  REAL(KIND=8) FUNCTION calculate_N(self,i,j)
    CLASS(OceanInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j
    TYPE(OceanicWavenumber) :: Di,Dj
    REAL(KIND=8) :: val

    Di = self%inner_products%owavenum(i)
    Dj = self%inner_products%owavenum(j)
    calculate_N = 0.D0
    IF (Dj%H/=Di%H) THEN
      val = delta(Di%P - Dj%P) * flambda(Di%H + Dj%H)
      calculate_N = val * (-2) * Dj%H * Di%H * self%inner_products%model_config%physics%n
      calculate_N = calculate_N / ((Dj%H**2 - Di%H**2) * self%inner_products%model_config%physics%pi)
    ENDIF
        
  END FUNCTION calculate_N

  !> Temperature advection term (passive scalar)
  !> 
  !> \f$ O_{i,j,k} = (\eta_i, J(\eta_j, \eta_k))\f$ .
  REAL(KIND=8) FUNCTION calculate_O(self,i,j,k)
    CLASS(OceanInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j,k
    TYPE(OceanicWavenumber) :: Di,Dj,Dk
    REAL(KIND=8) :: vs3,vs4,val
    INTEGER, DIMENSION(3) :: a
    INTEGER :: par

    val=0.D0

    a(1)=i
    a(2)=j
    a(3)=k

    CALL piksrt(3,a,par)

    Di = self%inner_products%owavenum(a(1))
    Dj = self%inner_products%owavenum(a(2))
    Dk = self%inner_products%owavenum(a(3))

    vs3 = S3(Dj%P,Dk%P,Dj%H,Dk%H)
    vs4 = S4(Dj%P,Dk%P,Dj%H,Dk%H)
    val = vs3*((delta(Dk%H - Dj%H - Di%H) - delta(Dk%H - Dj&
      &%H + Di%H)) * delta(Dk%P + Dj%P - Di%P) + delta(Dk&
      &%H + Dj%H - Di%H) * (delta(Dk%P - Dj%P + Di%P) -&
      & delta(Dk%P - Dj%P - Di%P))) + vs4 * ((delta(Dk%H &
      &+ Dj%H - Di%H) * delta(Dk%P - Dj%P - Di%P)) +&
      & (delta(Dk%H - Dj%H + Di%H) - delta(Dk%H - Dj%H -&
      & Di%H)) * (delta(Dk%P - Dj%P - Di%P) - delta(Dk%P &
      &- Dj%P + Di%P)))
    calculate_O = par * val * self%inner_products%model_config%physics%n / 2
  END FUNCTION calculate_O

  !> Streamfunction advection terms (oceanic)
  !> 
  !> \f$ C_{i,j,k} = (\eta_i, J(\eta_j,\nabla^2 \eta_k))\f$ .
  REAL(KIND=8) FUNCTION calculate_C_oc(self,i,j,k)
    CLASS(OceanInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j,k

    calculate_C_oc = self%M(k,k) * self%O(i,j,k)

  END FUNCTION calculate_C_oc

  !> Short-wave radiative forcing of the ocean
  !> 
  !> \f$ W_{i,j} = (\eta_i, F_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_W(self,i,j)
    CLASS(OceanInnerProducts), INTENT(IN) :: self
    INTEGER, INTENT(IN) :: i,j
    
    calculate_W = self%inner_products%atmos%s(j,i)

  END FUNCTION calculate_W

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  !> Initialization routine for the inner products functions.
  !> @param[in,out] inner_products Inner products global object to initialize
  !> @param[in] model_config Global model configuration object to initialize the inner products with
  SUBROUTINE init_inner_products(inner_products, model_config)
    CLASS(InnerProducts), INTENT(INOUT), TARGET :: inner_products
    CLASS(ModelConfiguration), INTENT(IN), TARGET :: model_config

    TYPE(InnerProducts), POINTER :: ips

    INTEGER :: i,j
    INTEGER :: AllocStat

    ips => inner_products

    IF (.NOT.model_config%initialized) THEN
      PRINT*, "Warning: Model configuration not initialized."
      PRINT*, "Aborting inner products initialization."
      RETURN
    END IF

    ! Definition of the types and wave numbers tables

    IF (allocated(ips%owavenum)) DEALLOCATE(ips%owavenum)
    ALLOCATE(ips%owavenum(model_config%modes%noc), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** init_inner_products: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF


    IF (allocated(ips%awavenum)) DEALLOCATE(ips%awavenum)
    ALLOCATE(ips%awavenum(model_config%modes%natm), STAT=AllocStat)
    IF (AllocStat /= 0) THEN
      PRINT*, "*** init_inner_products: Problem with allocation! ***"
      STOP "Exiting ..."
    END IF

    j=0
    DO i=1,model_config%modes%nbatm
      IF (model_config%modes%ams(i,1)==1) THEN
        ips%awavenum(j+1)%typ='A'
        ips%awavenum(j+2)%typ='K'
        ips%awavenum(j+3)%typ='L'

        ips%awavenum(j+1)%P=model_config%modes%ams(i,2)
        ips%awavenum(j+2)%M=model_config%modes%ams(i,1)
        ips%awavenum(j+2)%P=model_config%modes%ams(i,2)
        ips%awavenum(j+3)%H=model_config%modes%ams(i,1)
        ips%awavenum(j+3)%P=model_config%modes%ams(i,2)

        ips%awavenum(j+1)%Ny=REAL(model_config%modes%ams(i,2))
        ips%awavenum(j+2)%Nx=REAL(model_config%modes%ams(i,1))
        ips%awavenum(j+2)%Ny=REAL(model_config%modes%ams(i,2))
        ips%awavenum(j+3)%Nx=REAL(model_config%modes%ams(i,1))
        ips%awavenum(j+3)%Ny=REAL(model_config%modes%ams(i,2))

        j=j+3
      ELSE
        ips%awavenum(j+1)%typ='K'
        ips%awavenum(j+2)%typ='L'

        ips%awavenum(j+1)%M=model_config%modes%ams(i,1)
        ips%awavenum(j+1)%P=model_config%modes%ams(i,2)
        ips%awavenum(j+2)%H=model_config%modes%ams(i,1)
        ips%awavenum(j+2)%P=model_config%modes%ams(i,2)

        ips%awavenum(j+1)%Nx=REAL(model_config%modes%ams(i,1))
        ips%awavenum(j+1)%Ny=REAL(model_config%modes%ams(i,2))
        ips%awavenum(j+2)%Nx=REAL(model_config%modes%ams(i,1))
        ips%awavenum(j+2)%Ny=REAL(model_config%modes%ams(i,2))

        j=j+2

      ENDIF
    ENDDO

    DO i=1,model_config%modes%noc
      ips%owavenum(i)%H=model_config%modes%oms(i,1)
      ips%owavenum(i)%P=model_config%modes%oms(i,2)

      ips%owavenum(i)%Nx=model_config%modes%oms(i,1)/2.D0
      ips%owavenum(i)%Ny=model_config%modes%oms(i,2)

    ENDDO

    inner_products%model_config => model_config
    inner_products%atmos%inner_products => inner_products
    inner_products%ocean%inner_products => inner_products

    inner_products%initialized = .TRUE.

  END SUBROUTINE init_inner_products

  !> Routine to clean a inner products global object
  !> @param[in,out] inner_products Inner products global object to initialize
  SUBROUTINE delete_inner_products(inner_products)
    CLASS(InnerProducts), INTENT(INOUT), TARGET :: inner_products

    IF (allocated(inner_products%owavenum)) DEALLOCATE(inner_products%owavenum)
    IF (allocated(inner_products%awavenum)) DEALLOCATE(inner_products%awavenum)

    inner_products%initialized = .FALSE.

  END SUBROUTINE delete_inner_products

END MODULE inprod_analytic

