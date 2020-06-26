
!  maooam.f90
!
!> Fortran 2003 implementation of the modular arbitrary-order ocean-atmosphere
!> model MAOOAM.
!
!> @copyright                                                               
!> 2015-2020 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam 
  USE model_def
  USE rk2_integrator
  USE stat
  IMPLICIT NONE

  TYPE(Model), TARGET :: maooam_model
  TYPE(RK2Integrator) :: integr
  TYPE(StatAccumulator) :: stat_acc

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: t_up
  REAL(KIND=8), POINTER :: tw, tw_snap, t_trans, t_run
  INTEGER :: i, next, IndexSnap, WRSTAT
  INTEGER, POINTER :: ndim
  CHARACTER(LEN=9) :: arg
  LOGICAL :: cont_evol    !< True if the initial state is to be read in snapshot_trans.dat (i.e. the previous evolution is to be continued)
  LOGICAL :: ex
  LOGICAL, POINTER :: writeout

  PRINT*, 'Model MAOOAM v1.4'
  PRINT*, 'Loading information...'

  ! Initializing cont_evol
  cont_evol=.FALSE.
  arg='arg'
  i=0
  DO WHILE ((.NOT. cont_evol) .AND. LEN_TRIM(arg) /= 0)
     CALL get_command_argument(i, arg)
     IF (TRIM(arg)=='continue') cont_evol=.TRUE.
     i=i+1
  END DO

  CALL maooam_model%init
  CALL integr%init(maooam_model)

  tw => maooam_model%model_configuration%integration%tw
  tw_snap => maooam_model%model_configuration%integration%tw_snap
  t_trans => maooam_model%model_configuration%integration%t_trans
  t_run => maooam_model%model_configuration%integration%t_run
  writeout => maooam_model%model_configuration%integration%writeout
  ndim => maooam_model%model_configuration%modes%ndim

  t_up=integr%dt/t_trans*100.D0

  IF (writeout) OPEN(10,file='evol_field.dat')

  ALLOCATE(X(0:ndim),Xnew(0:ndim))

  X = maooam_model%load_IC()

  IF (cont_evol) THEN
     INQUIRE(FILE='snapshots_trans.dat',EXIST=ex,NEXTREC=next)
     IF (ex) THEN
        OPEN(11,file='snapshots_trans.dat')
        READ(11,REC=next-1) X
        CLOSE(11)
     END IF
  END IF

  IF (writeout) OPEN(11,file='snapshots_trans.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)

  PRINT*, 'Starting the transient time evolution...'

  IndexSnap=0
  DO WHILE (t<t_trans)
     CALL integr%step(X,t,Xnew)
     X=Xnew
     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
     IF (writeout .AND. mod(t,tw_snap)<integr%dt) THEN
        IndexSnap=IndexSnap+1
        WRITE(11,rec=IndexSnap,iostat=WRSTAT) X
     END IF
  END DO

  IF (writeout) CLOSE(11)

  PRINT*, 'Starting the time evolution...'

  CALL stat_acc%init(ndim)
  
  t=0.D0
  t_up=integr%dt/t_run*100.D0

  IF (writeout) WRITE(10,*) t,X(1:ndim)

  DO WHILE (t<t_run)
     CALL integr%step(X,t,Xnew)
     X=Xnew
     IF (mod(t,tw)<integr%dt) THEN
        IF (writeout) WRITE(10,*) t, X(1:ndim)
        CALL stat_acc%accumulate(X(1:ndim))
     END IF
     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)

  IF (writeout) OPEN(10,file='mean_field.dat')

  IF (writeout) WRITE(10,*) stat_acc%mean()
  IF (writeout) CLOSE(10)

  DEALLOCATE(X,Xnew)
  CALL maooam_model%clean
  CALL integr%clean
  CALL stat_acc%clean

END PROGRAM maooam 
