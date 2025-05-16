MODULE modInit

  USE modTOOLS
  USE modOutputs
  USE modCST
  USE modANALYTIC
  USE modINTEGRATION

  PRIVATE
  PUBLIC :: INITIAL_CONDITION
  PUBLIC :: ALLOCATE_ALL
  PUBLIC :: OPEN_ALL
  PUBLIC :: CLOSE_ALL
  !PUBLIC :: Initialise

CONTAINS

  SUBROUTINE OPEN_ALL(DATA)

    TYPE(DATA_IN), INTENT(IN) :: DATA
    INTEGER :: statut
    CHARACTER (len=30) :: dir
    CHARACTER (LEN=60) ::  NAME

    IF(DATA%ERROR ) THEN

       dir = 'post/'

       OPEN(UNIT = 200, FILE = TRIM(dir)//"Error.Case."//TRIM(carac(DATA%NCI))//"."//TRIM(ADJUSTL(DATA%WBmeth))//"."//&
            TRIM(carac(DATA%Nbx))//"cells."//TRIM(ADJUSTL(DATA%space_order))//"."//&
            TRIM(ADJUSTL(DATA%time_scheme))//".txt", &
            STATUS = "REPLACE", ACTION = "WRITE", &
            IOSTAT = statut)
       ! IF(statut /= 0) THEN
       !    CALL read__stop( "open_err Error File", statut)
       ! END IF

       NAME = "Time_Error.CASE."//TRIM(carac(DATA%Nci))//"."//TRIM(ADJUSTL(DATA%WBmeth))//"."//TRIM(ADJUSTL(DATA%err_type))&
            //"."//TRIM(ADJUSTL(DATA%space_order))//&
            "."//TRIM(ADJUSTL(DATA%time_scheme))//"."//TRIM(carac(DATA%Nbx))//".txt"

       OPEN(UNIT = 239, FILE = TRIM(dir)//NAME, STATUS = "REPLACE", ACTION = "WRITE", IOSTAT = statut)

    END IF

  END SUBROUTINE OPEN_ALL

  SUBROUTINE CLOSE_ALL(DATA)

    IMPLICIT NONE

    TYPE(DATA_IN), INTENT(IN) :: DATA
    INTEGER :: ierr
    INTEGER :: ii
    INTEGER :: NUM
    !CHARACTER (len=30) :: dir

    IF(DATA%TIME_SERIE) THEN
       CLOSE ( UNIT = 100, IOSTAT=ierr)
       IF (ierr /= 0) THEN
          STOP "pb close 100"
       END IF
    END IF

    IF(DATA%ERROR) THEN
       CLOSE ( UNIT = 200, IOSTAT=ierr)
       IF (ierr /= 0) THEN
          STOP "pb close 200"
       END IF

       CLOSE ( UNIT = 239, IOSTAT=ierr)
       IF (ierr /= 0) THEN
          STOP "pb close 239"
       END IF
    END IF

    IF(DATA%TIME_SERIE) THEN
       DO ii = 1, DATA%point_number

          NUM = 170+ii
          CLOSE(UNIT = NUM, IOSTAT=ierr)
          IF (ierr /= 0) THEN
             STOP "pb close file time_series"
          END IF

       END DO
    END IF

    !WRITE(UNIT = 6, FMT = *) 'All files are closed now'

  END SUBROUTINE CLOSE_ALL

  SUBROUTINE ALLOCATE_ALL(V, FLU, MESH, DATA, RK)

    IMPLICIT NONE

    TYPE(MeshDef), INTENT(INOUT) :: MESH
    INTEGER :: istat
    TYPE(VFLU), INTENT(INOUT) :: FLU
    TYPE(DATA_IN), INTENT(IN) :: DATA
    TYPE(Variable), INTENT(INOUT) :: V
    TYPE(coeff_RK), DIMENSION(4), INTENT(INOUT) :: RK
    INTEGER :: Nbx
    INTEGER :: Nx
    !INTEGER :: Ndf
    INTEGER :: Nc
    INTEGER :: i
    !INTEGER :: j

    CHARACTER(LEN =*), PARAMETER :: sub__name = "ALLOCATE_ALL"

    Nbx = DATA%Nbx
    !NDf = DATA%NDf

    IF(debug) THEN
       !PRINT*, "Entering...", sub__name
    END IF

    ALLOCATE(MESH%X_VF(-2:DATA%Nbx+3), stat = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE X_VF"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating X_VF. STOP"
       STOP
    ENDIF

   ALLOCATE ( MESH%X_DF (-2: DATA%Nbx+4), STAT = istat )
    IF (debug) THEN
       !PRINT *, sub__name, " ALLOCATE X_DF"
    END IF
    IF (istat /= 0 ) THEN
       STOP "alloc_err L==3491"
    END IF

    ALLOCATE(MESH%array_DX(-2:DATA%Nbx+3), stat = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE array_DX"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating array_DX. STOP"
       STOP
    ENDIF


    ALLOCATE(V%t1(-2:Nbx+3)   , STAT = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE Tab"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating Tab. STOP"
       STOP
    ENDIF

    ALLOCATE(V%VF(-2:Nbx+3)  , STAT = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE TAB1"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating Tab1. STOP"
       STOP
    ENDIF


    ALLOCATE(V%t2(-2:Nbx+3)  , STAT = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE Tab2"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating Tab2. STOP"
       STOP
    ENDIF

    ALLOCATE(V%t3(-2:Nbx+3)  , STAT = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE Tab3"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating Tab3. STOP"
       STOP
    ENDIF

    ALLOCATE(V%t4(-2:Nbx+3)  , STAT = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE Tab3"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating Tab3. STOP"
       STOP
    ENDIF

    ALLOCATE(V%tmp(-2:Nbx+3), STAT = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE TabTemp"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating TabTemp. STOP"
       STOP
    ENDIF

    ALLOCATE(V%topoVF(-2:Nbx+3), STAT = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE TabTemp"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating TabTemp. STOP"
       STOP
    ENDIF


    ALLOCATE (FLU%Sx(0:Nbx), STAT = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE FLU_Sx"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating Flu_Sx. STOP"
       STOP
    ENDIF

    ALLOCATE (FLU%X1(0:Nbx), STAT = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE FLU_X1"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating FLU_X1. STOP"
       STOP
    ENDIF

    ALLOCATE (FLU%X2(0:Nbx), STAT = istat)
    IF(debug) THEN
       !PRINT*, sub__name, " ALLOCATE FLU_X2"
    END IF
    IF (istat /= 0) THEN
       !PRINT *,"Error allocating FLU_X2. STOP"
       STOP
    ENDIF


    IF(DATA%solex.OR.DATA%Sol_Ref) THEN
       ALLOCATE (V%Ex(-2:Nbx+3), STAT = istat )
       IF (debug ) THEN
          !PRINT *, sub__name, " ALLOCATE V%Ex"
       END IF
       IF (istat /= 0 ) THEN
          STOP "alloc_err V%Ex"
       END IF
    END IF


    IF(DATA%Sol_Ref) THEN
       ALLOCATE ( V%Ref(-2:DATA%NbxRef+3), STAT = istat )
       IF (debug ) THEN
          !PRINT *, sub__name, " ALLOCATE TabRef"
       END IF
       IF (istat /= 0 ) THEN
          STOP "alloc_err TabRef"
       END IF
    END IF

    IF(DATA%TIME_SERIE) THEN

       ALLOCATE ( V%t0 (-2:Nbx+3), STAT = istat)
       IF (debug) THEN
          !PRINT *, sub__name, " ALLOCATE V%t0"
       END IF
       IF (istat /= 0 ) THEN
          STOP "alloc_err L==2222"
       END IF

       DO i=-2, Nbx+3
          V%t0(i)%h = V%VF(i)%h
          V%t0(i)%qx = V%VF(i)%qx

       END DO


    END IF




    ALLOCATE(RK(1)%a(1:1,1:1))
    ALLOCATE(RK(1)%b(1:1))
    ALLOCATE(RK(1)%s_dt(1:1))

    ALLOCATE(RK(2)%a(1:2,1:2))
    ALLOCATE(RK(2)%b(1:2))
    ALLOCATE(RK(2)%s_dt(1:2))

    ALLOCATE(RK(3)%a(1:3,1:3))
    ALLOCATE(RK(3)%b(1:3))
    ALLOCATE(RK(3)%s_dt(1:3))


    ALLOCATE(RK(4)%a(1:4,1:4))
    ALLOCATE(RK(4)%b(1:4))
    ALLOCATE(RK(4)%s_dt(1:4))

    RK(1)%a(1,1) = 0.0D0
    RK(1)%b(1) = 1.0D0
    RK(1)%s_dt(1) = 1.0D0


    RK(2)%a(1,1) = 0.0D0
    RK(2)%a(1,2) = 0.0D0
    RK(2)%a(2,1) = 0.5D0
    RK(2)%a(2,2) = 0.0D0

    RK(2)%b(1) = 0.75D0
    RK(2)%b(2) = 0.25D0

    !RK(2)%b = (/0.75D0, 0.25D0/)
    RK(2)%s_dt = (/0.5D0, 1.0D0/)

    RK(3)%a(1,1) = 0.0D0
    RK(3)%a(1,2) = 0.0D0
    RK(3)%a(1,3) = 0.0D0

    RK(3)%a(2,1) = 0.5D0
    RK(3)%a(2,2) = 0.0D0
    RK(3)%a(2,3) = 0.0D0

    RK(3)%a(3,1) = -1.0D0
    RK(3)%a(3,2) = 2.0D0
    RK(3)%a(3,3) = 0.0D0

    RK(3)%b(1) = 1.0D0/6.0D0
    RK(3)%b(2) = 2.0D0/3.0D0
    RK(3)%b(3) = 1.0D0/6.0D0

    RK(4)%a(1,1) = 0.0D0
    RK(4)%a(1,2) = 0.0D0
    RK(4)%a(1,3) = 0.0D0
    RK(4)%a(1,4) = 0.0D0

    RK(4)%a(2,1) = 0.5D0
    RK(4)%a(2,2) = 0.0D0
    RK(4)%a(2,3) = 0.0D0
    RK(4)%a(2,4) = 0.0D0

    RK(4)%a(3,1) = 0.0D0
    RK(4)%a(3,2) = 0.5D0
    RK(4)%a(3,3) = 0.0D0
    RK(4)%a(3,4) = 0.0D0

    RK(4)%a(4,1) = 0.0D0
    RK(4)%a(4,2) = 0.0D0
    RK(4)%a(4,3) = 1.0D0
    RK(4)%a(4,4) = 0.0D0

    RK(4)%b(1) =  1.0D0/6.0D0
    RK(4)%b(2) = 1.0D0/3.0D0
    RK(4)%b(3) = 1.0D0/3.0D0
    RK(4)%b(4) = 1.0D0/6.0D0

    RK(4)%s_dt = (/0.5D0, 0.5D0, 1.0D0, 1.0D0/)


    !WRITE(6,*)'All allocations done !!! '


  END SUBROUTINE ALLOCATE_ALL




  SUBROUTINE INITIAL_CONDITION &
       (MESH, DATA, V, debug)

    IMPLICIT NONE

    TYPE(DATA_IN), INTENT(IN) :: DATA
    TYPE(MeshDef), INTENT(IN) :: MESH
    TYPE(Variable), INTENT(INOUT) :: V

    REAL(kind = 8) :: x0
    real*8 :: sigma , sd , sh , temp_x = 0.0
    REAL(kind = 8) :: x1
    REAL(kind = 8) :: cotbeta
    REAL(kind = 8) :: ksi
    REAL(kind = 8) :: U
    REAL(kind = 8) :: H
    REAL(kind = 8) :: L
    REAL(kind = 8) :: X
    REAL(kind = 8) :: D
    REAL(kind = 8) :: A
    REAL(kind = 8) :: h0
    REAL(kind = 8) :: f__4514
    REAL(kind = 8) :: k
    REAL(kind = 8) :: z0
    REAL(kind = 8) :: c
    INTEGER :: i
    INTEGER :: Nbx
    LOGICAL, INTENT(in) :: debug
    CHARACTER(LEN =*), PARAMETER :: sub__name = "Inicon"
    
    Nbx = DATA%Nbx

    IF(debug) THEN
       !PRINT*, "Entering...", sub__name
    END IF

    ! default initialisation
    DO i = -2, Nbx+3

       V%topoVF(i) = 0.0D0
    END DO

    SELECT CASE(DATA%NCI)


    CASE(0)

       !dam break
       

       DO i = -2, Nbx/2

     
          ! Initialisation partie VF
          V%topoVF(i) = -1.0D0
          V%VF(i)%h = 2.0d0
          V%VF(i)%qx = 0.0D0
       END DO

       DO i = Nbx/2+1, Nbx+3

          ! Initialisation partie VF
          V%topoVF(i) = -1.0D0
          V%VF(i)%h  = 1.0d0
          V%VF(i)%qx = 0.0D0
       END DO
       
    case(1964)

       ! Effondrement d une gaussienne
       x0 = 500.d0
       sigma = 50d0
       sd = 50.0d0/9.8d0 !50.0d0 / 9.8d0;
       sd = sd - 1.0d0
       sh = 0.1*sd;

       do i = -2 , Nbx/2
          V%topoVF(i) = -1.0d0 !2.0D0
       end do
       do i = Nbx/2 + 1 , Nbx + 3
          V%topoVF(i) = -1.0d0 !3.0D0
       end do


       do i = -2 , Nbx + 3
          !print*,"Voici Nbx = ",Nbx
          V%VF(i)%h = sd + sh*exp(-((MESH%X_DF(i)-x0)/sigma)**2) - V%topoVF(i) !- V%topoVF(i)
          !V%topoVF(i) =  0.0d0 !((0.02d0 + 0.03d0)/2.0d0) + (1.0d-02*tanh(10.0d0*(temp_x - 5.0d0))/2.0d0)
          V%VF(i)%qx = 0.0D0
          !temp_x = temp_x + 1.0d-02
       end do



       
    case(1996)

       ! Flow at rest
       do i = -2 , Nbx + 3
!!$          if( MESH%X_DF(i) < 12.0d0 .and. MESH%X_DF(i) > 8.0d0 ) then         
!!$             V%topoVF(i) = 0.2d0 - 0.05d0*(MESH%X_DF(i) - 10)**2
!!$          else
!!$             V%topoVF(i) = 0.0d0            
!!$          end if

          if ( MESH%X_DF(i) < 0.6d0 .and. MESH%X_DF(i) > 0.4d0 ) then
             V%topoVF(i) = 0.0d0 !-10.0d0*(MESH%X_DF(i) - 0.4d0)*(MESH%X_DF(i)-0.6d0) !0.2d0 - 0.1d0*((x - 5.0d-01)**2)
          else
             V%topoVF(i) = 0.0d0
          end if
          
          V%VF(i)%h = 1.0d0 !0.0d0 !2.0d0 !max(V%topoVF(i),0.15) 
          V%VF(i)%qx = 0.0D0
          !V%topoVF(i) = 0.0d0 !-1.0d0
          !temp_x = temp_x + 0.083 ! 
          !print*,"Voici x",temp_x
       end do
    CASE(1)

       !test de generation => soliton
       WRITE(6,*)'Test de génération d un soliton'
       H0 = 0.5D0

       DO i = -2, Nbx+3
          ! Initialisation partie VF
          V%topoVF(i) = 0.0D0
          V%VF(i)%h = H0
          V%VF(i)%qx = 0.0D0
       END DO


    CASE(2)

       !test de generation

       H0 = 4.20D0

       DO i = -2, Nbx+3

          ! Initialisation partie VF
          V%topoVF(i) = 0.0D0
          V%VF(i)%h = H0
          V%VF(i)%qx = 0.0D0

       END DO


    CASE(29)

       !Soliton sans topo 

       H0 = 0.50D0
       A = 0.2D0*H0
       k = sqrt(3.0D0*a/(4.0D0*H0**2*(H0+A)))
       c = sqrt(g_r__8*(H0+a))
       L = DATA%x_out - DATA%x_in
       d = L/2.D0
       !d=3.0D0

       DO i = -2, Nbx+3

          X = MESH%X_VF(i) !)(MESH%X_DF(i) + MESH%X_DF(i+1))/2.0D0
          !WRITE(6,*)'X IN INIT', X
          !Z0 = X - d

          Ksi = a/(cosh(k*(X-d))**2)

          U = c*(1.0D0-H0/(Ksi+H0))

          ! Initialisation partie VF
          V%topoVF(i) = 0.0D0
          V%VF(i)%h = H0 + ksi
          !WRITE(6,*)'V%VF(i)%h',V%VF(i)%h
          V%VF(i)%qx = U*(V%VF(i)%h)

       END DO


    CASE(290)

       !Soliton, configuration test Delis

       H0 = 10.0D0
       A = 0.2D0*H0
       k = sqrt(3.0D0*a/(4.0D0*H0**2*(H0+A)))
       c = sqrt(g_r__8*(H0+a))
       L = DATA%x_out - DATA%x_in
       d = 200.0D0
       !d=3.0D0

       DO i = -2, Nbx+3

          X = MESH%X_VF(i) !)(MESH%X_DF(i) + MESH%X_DF(i+1))/2.0D0
          !WRITE(6,*)'X IN INIT', X
          !Z0 = X - d

          Ksi = a/(cosh(k*(X-d))**2)

          U = c*(1.0D0-H0/(Ksi+H0))

          ! Initialisation partie VF
          V%topoVF(i) = 0.0D0
          V%VF(i)%h = H0 + ksi
          !WRITE(6,*)'V%VF(i)%h',V%VF(i)%h
          V%VF(i)%qx = U*(V%VF(i)%h)

       END DO



    CASE(40)

       !cas Cox + generation
       H0 = 0.4D0

       DO i = -2, Nbx+3

          X = MESH%X_VF(i) !)(MESH%X_DF(i) + MESH%X_DF(i+1))/2.0D0

          !ksi = H*1.0D0/(COSH(SQRT((3.0D0*H) & 
          !/(4.0D0*D**3))*(X-X1) )**2 )

          !U = KSI*SQRT(g_r__8/D)

          ! Initialisation partie VF
          IF(X<=2.0D0+4.5D0) THEN
             V%topoVF(i) = 0.0D0
          ELSE
             V%topoVF(i) = 0.0D0 !1.0D0/35.0D0*(X-6.5D0)
          END IF
          V%VF(i)%h = MAX(0.0D0,H0-V%topoVF(i))
          V%VF(i)%qx = 0.0D0 !U*(V%VF(i)%h)

       END DO

     

    CASE(41)

       !cas Dingemans + generation
       H0 = 0.4D0

       DO i = -2, Nbx+3

          X = MESH%X_VF(i) 

          ! Initialisation partie VF
          IF(X<=2.0D0+4.5D0) THEN
             V%topoVF(i) = 0.0D0
          ELSE
             V%topoVF(i) = 0.0D0 !1.0D0/35.0D0*(X-6.5D0)
          END IF
          V%VF(i)%h = MAX(0.0D0,H0-V%topoVF(i))
          V%VF(i)%qx = 0.0D0 !U*(V%VF(i)%h)

       END DO

       

    CASE(cte__30)

       !Soliton sans topo
       X0 = 60.15D0

       D = 1.0D0
       H = 0.19D0
       f__4514 = dacosh( SQRT (( 1.0D0 / 0.05D0 ) ))
       X1 = X0 - SQRT((4.0D0*D)/(3.0D0*H))*f__4514

       DO i = -2, Nbx+3

          X = MESH%X_VF(i) !)(MESH%X_DF(i) + MESH%X_DF(i+1))/2.0D0

          ksi = H*1.0D0/(COSH(SQRT((3.0D0*H) & 
               /(4.0D0*D**3))*(X-X1) )**2 )

          U = KSI*SQRT(g_r__8/D)

          ! Initialisation partie VF
          V%topoVF(i) = 0.0D0
          V%VF(i)%h = D + ksi
          V%VF(i)%qx = U*(V%VF(i)%h)

       END DO

       


    CASE(cte__31)

       !Synolakis

       X0 = 60.15D0
       cotBeta =  19.85D0
       !cotBeta = cte__9_d_85_r__8
       !cotbeta = 12.0D0

       D = 1.0D0
       H = 0.019D0
       f__4514 = acosh( SQRT (( 1.0D0 / 0.05D0 ) ))
       X1 = X0 - SQRT((4.0D0*D)/(3.0D0*H))*f__4514

       DO i = -2, Nbx+3

          X = (MESH%X_DF(i) + MESH%X_DF(i+1))/2.0D0

          ksi = H*1.0D0/(COSH(SQRT((3.0D0*H) & 
               /(4.0D0*D**3))*(X-X1) )**2 )

          U = KSI*SQRT(g_r__8/D)

          ! Initialisation partie VF
          IF (X <= X0) THEN
             V%topoVF(i) = 0.0D0
             V%VF(i)%h = D + ksi
             V%VF(i)%qx = U*(V%VF(i)%h)
          ELSE 
             V%topoVF(i) = 1.0D0/cotBeta*(X - X0)
             V%VF(i)%h = MAX(D + ksi - V%topoVF(i), 0.0D0)
             V%VF(i)%qx = U*V%VF(i)%h
          END IF

       END DO

       DO i = 1, DATA%NDf

          X = MESH%X_DF(i)

          ksi = H*1.0D0/(COSH(SQRT((3.0D0*H) & 
               /(4.0D0*D**3))*(X-X1) )**2 )

          U = KSI*SQRT(g_r__8/D)



       END DO


    END SELECT


    ! Update wet/dry + eta
    SELECT CASE(DATA%cond_lim)

    CASE(2,3,4)
       Var%TopoVf(0) = Var%TopoVf(1)
       Var%TopoVf(-1) = Var%TopoVf(2)
       Var%TopoVf(-2) = Var%TopoVf(3)

       Var%TopoVf(Nbx+1) = Var%TopoVf(Nbx)
       Var%TopoVf(Nbx+2) = Var%TopoVf(Nbx-1)
       Var%TopoVf(Nbx+3) = Var%TopoVf(Nbx-2)

    CASE(1)
       Var%TopoVf(0) = Var%TopoVf(Nbx)
       Var%TopoVf(-1) = Var%TopoVf(Nbx-1)
       Var%TopoVf(-2) = Var%TopoVf(Nbx-2)

       Var%TopoVf(Nbx+1) = Var%TopoVf(1)
       Var%TopoVf(Nbx+2) = Var%TopoVf(2)
       Var%TopoVf(Nbx+3) = Var%TopoVf(3)

    END SELECT


    WHERE (V%VF%h <= DATA%kappa)    
       V%VF%h  = DATA%kappa  
       V%VF%qx = 0.0D0
       V%VF%u = 0.0D0

    ELSEWHERE ((DATA%kappa <= V%VF%h) .AND. (V%VF%h <= DATA%kappabis))    
       V%VF%qx = 0.0D0
       V%VF%u = 0.0D0
    END WHERE

    V%VF%xi =  V%VF%h + V%topoVF


  END SUBROUTINE INITIAL_CONDITION

END MODULE ModInit
