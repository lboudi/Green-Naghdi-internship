MODULE modPOST

  !USE fortranf
  USE modCST
  USE modTypes
  USE modTOOLS
  USE modINTEGRATION
  USE modANALYTIC

  IMPLICIT NONE

  PUBLIC :: NORM_Error
  PUBLIC :: TIME_ERROR
  PUBLIC :: Solution_Exact
  PUBLIC :: Error
  PUBLIC :: Compute_errors

CONTAINS

  SUBROUTINE NORM_ERROR(Tab, TabEx, TabRef, DATA, Time, &
       sum_init, NN, ERR)

    IMPLICIT NONE

    TYPE(DATA_IN), INTENT(IN) :: DATA
    TYPE(T_ERR), INTENT(INOUT) :: ERR
    CHARACTER (LEN =*), PARAMETER :: c__format__8 = "(2ES16.8)"
    REAL (KIND = 8), INTENT ( IN ) :: Time
    REAL (KIND = 8), INTENT ( OUT ) :: sum_init
    INTEGER :: Nbx
    INTEGER :: i
    INTEGER :: j
    REAL (KIND = 8) :: sum
    REAL (KIND = 8) :: Error_L2
    REAL (KIND = 8) :: l
    REAL (KIND = 8) :: S
    REAL (KIND = 8) :: atime
    REAL (KIND = 8) :: Error_H
    REAL (KIND = 8) :: Error_U
    REAL (KIND = 8) :: XX
    TYPE(STATE), DIMENSION(:), POINTER :: Tab
    TYPE(STATE), DIMENSION( :), POINTER :: tabEx
    TYPE(STATE), DIMENSION( :), POINTER :: TabRef
    CHARACTER (LEN =*), PARAMETER :: sub__name = "func_error_l2"
    CHARACTER (len=30) :: dir
    INTEGER :: sta___tut38
    !INTEGER :: istat
    CHARACTER (LEN=* ), PARAMETER :: c__format__6 = "(6ES16.8)"
    INTEGER :: NbxRef
    INTEGER, INTENT(IN) :: NN

    IF (debug) THEN
       !PRINT *, "Entering ", sub__name
    END IF

    Nbx = DATA%Nbx

    NbxRef = 16384

    SELECT CASE(DATA%NCI) 

    CASE(cte__16)

       ! accuracy 1D CW Shu
       IF(NN == 1) THEN
          !read reference solution
          dir = 'DATA/'
          OPEN(UNIT = 120, FILE = TRIM(dir)//DATA%char__SolRef, STATUS = "OLD", ACTION = "READ", &
               IOSTAT = sta___tut38)
          ! IF(sta___tut38 /= 0) THEN
          !    CALL read__stop( "open_err reference solution - Cas 16", sta___tut38)
          ! END IF


          DO i = -2, NbxRef+3
             READ(UNIT = 120, FMT = c__format__6)XX, TabRef(i)%h, TabRef(i)%u, &
                  TabRef(i)%topo, TabRef(i)%xi, TabRef(i)%qx
          END DO

          !WRITE(6,*)'REFERENCE SOLUTION READ'

          !WRITE(6,*)'Computation of Exact solution with Nbx cells from'
          !WRITE(6,*)'     the reference solution with 16384 cells by integration'

          CALL Ex_from_Ref(TabEx, TabRef, DATA)

       END IF

       !WRITE(6,*)'Computation of L2-error'

       CALL Compute_Errors(Tab, TabEx, DATA, ERR)

    CASE(cte__44, cte__45)

       !        ! accuracy 1D, solution with friction+diffusion, sub/sub

       ! for this specific test, reference solution is read and sampled
       ! on current coarse mesh directly in the INITIALISE routine

       !WRITE(6,*)'Computation of L2-error'

       CALL Compute_Errors(Tab, TabEx, DATA, ERR)

    CASE(cte__47)

       !        ! accuracy 1D, solution with friction, sub/sub

       ! for this specific test, reference solution is read and sampled
       ! on current coarse mesh directly in the INITIALISE routine

       !WRITE(6,*)'Computation of L2-error'

       CALL Compute_Errors(Tab, TabEx, DATA, ERR)

    CASE(cte__37)

       !Cas Carrier & Greenspan periodic adimensionné

       l = cte__20_d_0_r__8
       S = 1.0D0 / cte__30_d_0_r__8
       j=1   

       Error_L2 =  cte__0_d_r__8
       Error_H =  cte__0_d_r__8
       Error_U =  cte__0_d_r__8
       sum =  cte__0_d_r__8
       sum_init =  cte__0_d_r__8

       !DO i= 1, Nbx
       !   sum = sum + (Tab(i)%h-TabEx(i)%h)**2 + (Tab(i)%qx-TabEx(i)%qx)**2
       !END DO

       !DO i= 1, Nbx
       !   sum_init = sum_init + TabEx(i)%h**2 + TabEx(i)%qx**2
       !END DO
       !Error_L2 = SQRT( sum/sum_init )

       DO i= 1, Nbx
          sum = sum + (Tab(i)%h-TabEx(i)%h)**2 
       END DO

       DO i= 1, Nbx
          sum_init = sum_init + TabEx(i)%h**2 + TabEx(i)%qx**2
       END DO

       Error_H = SQRT( sum/sum_init )
       aTime = Time/(SQRT(l/(S*g_r__8)))

       WRITE (UNIT = 210, FMT = c__format__8) atime, Error_H

       sum =  cte__0_d_r__8
       sum_init =  cte__0_d_r__8

       DO i= 1, Nbx
          sum = sum + (Tab(i)%qx-TabEx(i)%qx)**2 
       END DO

       DO i= 1, Nbx
          sum_init = sum_init + TabEx(i)%qx**2  + TabEx(i)%h**2
       END DO

       Error_U = SQRT( sum/sum_init )
       aTime = Time/(SQRT(l/(S*g_r__8)))

       WRITE (UNIT = 220, FMT = c__format__8) atime, Error_U

       !WRITE(200,'(2E16.8)')aTime, Error_L2


    CASE(12)

       !Sub-critical flow over a bump
       !WRITE(UNIT=6, FMT=*)'Computing error : Entering SolEx Case 12 ...'

       !PRINT*,'err_type',DATA%err_type

       CALL Compute_Errors(Tab, TabEx, DATA, ERR)

    CASE(1996)

       !Sub-critical flow over a bump
       !WRITE(UNIT=6, FMT=*)'Computing error : Entering SolEx Case 12 ...'

       !PRINT*,'err_type',DATA%err_type

       CALL Compute_Errors(Tab, TabEx, DATA, ERR)

    CASE(cte__26)

       !Sub-critical flow over a gaussian bump
       ! accuracy test for 2nd and 4th order schemes

       !WRITE(UNIT=6, FMT=*)'M_Tools : Entering SolEx Case 26 ...'
       !WRITE(UNIT=6, FMT=*)'err_type',DATA%err_type

       CALL Compute_Errors(Tab, TabEx, DATA, ERR)

    CASE(290)

       ! Soliton - cas Delis-Kazeola

       !WRITE(6,*)'Computation of L2-error'

       CALL Compute_Errors(Tab, TabEx, DATA, ERR)

    END SELECT

  END SUBROUTINE NORM_ERROR


  SUBROUTINE Compute_Errors(Tab, TabEx, DATA, ERR)

    IMPLICIT NONE

    TYPE(T_ERR), INTENT(INOUT) :: ERR
    TYPE(STATE), DIMENSION(:), POINTER :: TabEx
    TYPE(STATE), DIMENSION(:), POINTER :: Tab
    TYPE(DATA_IN), INTENT(IN) :: DATA
    INTEGER :: i
 
    INTEGER :: Nbx

    REAL (KIND = 8) :: sum_init
    REAL (KIND = 8) :: sum

    !REAL (KIND = 8) :: Error_L2
    REAL (KIND = 8) :: CError
    !REAL (KIND = 8) :: l
    !REAL (KIND = 8) :: S
    !REAL (KIND = 8) :: atime
    !REAL (KIND = 8) :: Error_H
    !REAL (KIND = 8) :: Error_U
    !REAL (KIND = 8) :: XX
    REAL (KIND = 8) :: dx

    dx = MESH%Dx

    Nbx = DATA%Nbx

    CError =  0.0D0
    sum =  0.0D0
    sum_init =  0.0D0

    SELECT CASE(DATA%err_type)

    CASE('L2')
       !WRITE(UNIT=6, FMT=*)'L2-Error computation ...'
       !! For h
       DO i= 1, Nbx
          sum_init = sum_init + TabEx(i)%h**2
       END DO

       DO i= 1, Nbx
          sum = sum + (Tab(i)%h-TabEx(i)%h)**2
       END DO

       ERR%L2_h =  SQRT( sum/sum_init )

       ! ===== ON SCREEN ======
       !ERR%L2_h = SQRT( sum*MESH%dx )
       !WRITE(UNIT=6, FMT=*)'===== L2-error ===='
       WRITE(UNIT=200, FMT=*)'===== L2-error ===='

       !WRITE (UNIT = 6, FMT = *)'h,', LOG(dx), LOG(ERR%L2_h)
       !WRITE(UNIT=6, FMT=*),'Error,', ERR%L2_h
       !WRITE(UNIT=6, FMT=*)'dx,',dx

       ! ====== IN FILE ========


       WRITE (UNIT = 200, FMT = *)'h,', LOG(dx), LOG(ERR%L2_h)
       WRITE(UNIT=200, FMT=*),'Error,', ERR%L2_h
       WRITE(UNIT=200, FMT=*)'dx,',dx


       !! For qx
       sum =  cte__0_d_r__8
       sum_init =  cte__0_d_r__8

       DO i= 1, Nbx
          sum = sum + (Tab(i)%qx-TabEx(i)%qx)**2 
       END DO

       DO i= 1, Nbx
          sum_init = sum_init + TabEx(i)%qx**2  
       END DO

       !ERR%L2_qx = SQRT(sum*MESH%dx )
       ERR%L2_qx  = SQRT( sum/sum_init )
       ! ======= ON SCREEN =====
       !WRITE (UNIT = 6, FMT = *)'qx,', LOG(dx), LOG(ERR%L2_qx)
       !WRITE(UNIT=6, FMT=*),'Error,', ERR%L2_qx
       !WRITE(UNIT=6, FMT=*)'dx,',dx

       ! ====== IN FILE =========
       WRITE(UNIT=200, FMT = *)'qx,', LOG(dx), LOG(ERR%L2_qx)
       WRITE(UNIT=200, FMT=*),'Error,', ERR%L2_qx
       WRITE(UNIT=200, FMT=*)'dx,',dx

    CASE('L1')

       DO i= 1, Nbx
          sum_init = sum_init + ABS(TabEx(i)%h)
       END DO

       DO i= 1, Nbx
          sum = sum + ABS(Tab(i)%h-TabEx(i)%h)
       END DO

       !CError = sum/sum_init
       ERR%L1_h = sum*MESH%dx

       !WRITE(UNIT=6, FMT=*)'==== L1-error ===='
       WRITE(UNIT=200, FMT=*)'==== L1-error ===='

       !WRITE (UNIT = 6, FMT = *)'h,', LOG(dx), LOG(ERR%L1_h)
       !WRITE(UNIT=6, FMT=*),'Error,', ERR%L1_h
       !WRITE(UNIT=6, FMT=*)'dx,',dx

       ! ====== IN FILE =========
       WRITE(UNIT=200, FMT = *)'h,', LOG(dx), LOG(ERR%L1_h)
       WRITE(UNIT=200, FMT=*),'Error,', ERR%L1_h
       WRITE(UNIT=200, FMT=*)'dx,',dx

       sum =  cte__0_d_r__8
       sum_init =  cte__0_d_r__8

       DO i= 1, Nbx
          sum_init = sum_init + ABS(TabEx(i)%qx)
       END DO

       DO i= 1, Nbx
          sum = sum + ABS(Tab(i)%qx-TabEx(i)%qx)
       END DO

       !CError = sum/sum_init
       ERR%L1_qx = sum*MESH%dx

       !WRITE (UNIT = 6, FMT = *)'qx,', LOG(dx), LOG(ERR%L1_qx)
       !WRITE(UNIT=6, FMT=*),'Error,',ERR%L1_qx
       !WRITE(UNIT=6, FMT=*)'dx,',dx

       ! ====== IN FILE =========
       WRITE(UNIT=200, FMT = *)'qx,', LOG(dx), LOG(ERR%L1_qx)
       WRITE(UNIT=200, FMT=*),'Error,',ERR%L1_qx
       WRITE(UNIT=200, FMT=*)'dx,',dx

    CASE('LINF')

       !DO i= 1, Nbx
       !   sum_init =  TabEx(i)%h**2
       !END DO

       ERR%LINF_h = MAXVAL(ABS(Tab(1:Nbx)%h-TabEx(1:Nbx)%h))

       !WRITE(UNIT=6, FMT=*)'==== LINF-error ===='
       WRITE(UNIT=200, FMT=*)'==== LINF-error ===='

       !WRITE (UNIT = 6, FMT = *)'h,', LOG(dx), LOG(ERR%LINF_h)
       !WRITE(UNIT=6, FMT=*),'Error,',ERR%LINF_h
       !WRITE(UNIT=6, FMT=*)'dx,',dx

       ! ====== IN FILE =========
       WRITE(UNIT=200, FMT = *)'h,', LOG(dx), LOG(ERR%LINF_h)
       WRITE(UNIT=200, FMT=*),'Error,',ERR%LINF_h
       WRITE(UNIT=200, FMT=*)'dx,',dx

       ERR%LINF_qx = MAXVAL(ABS(Tab(1:Nbx)%qx-TabEx(1:Nbx)%qx))

       !WRITE (UNIT = 6, FMT = *)'qx,', LOG(dx), LOG(ERR%LINF_qx)
       !WRITE(UNIT=6, FMT=*),'Error,',ERR%LINF_qx
       !WRITE(UNIT=6, FMT=*)'dx,',dx

       ! ====== IN FILE =========
       WRITE(UNIT=200, FMT = *)'qx,', LOG(dx), LOG(ERR%LINF_qx)
       WRITE(UNIT=200, FMT=*),'Error,',ERR%LINF_qx
       WRITE(UNIT=200, FMT=*)'dx,',dx

    CASE('ALL')

       !! L2
       !! For h
       sum = 0.0D0

       DO i= 1, Nbx
          sum_init = sum_init + TabEx(i)%h**2
       END DO

       DO i= 1, Nbx
          sum = sum + (Tab(i)%h-TabEx(i)%h)**2 !+ (Tab(i)%qx-TabEx(i)%qx)**2 &
          !+ (Tab(i)%qy-TabEx(i)%qy)**2
       END DO

       !CError = SQRT( sum/sum_init )
       ERR%L2_h = SQRT( sum*MESH%dx )

       !WRITE(UNIT=6, FMT = *)'====== L2-error ======'
       !WRITE (UNIT = 6, FMT = *)'Pour h,', 'dx=', dx, 'Error=', ERR%L2_h
       !WRITE (UNIT = 6, FMT = *)'Pour h,', LOG(dx), LOG(ERR%L2_h)

       ! ====== IN FILE =========
       WRITE(UNIT=200, FMT = *)'====== L2-error ======'
       WRITE (UNIT = 200, FMT = *)'Pour h,', 'dx=', dx, 'Error=', ERR%L2_h
       WRITE (UNIT = 200, FMT = *)'Pour h,', LOG(dx), LOG(ERR%L2_h)


       sum =  cte__0_d_r__8
       sum_init =  cte__0_d_r__8

       DO i= 1, Nbx
          sum = sum + (Tab(i)%qx-TabEx(i)%qx)**2 
       END DO

       DO i= 1, Nbx
          sum_init = sum_init + TabEx(i)%qx**2  
       END DO
       !print*,' For qx'
       !CError = SQRT( sum/sum_init )
       ERR%L2_qx = SQRT( sum*MESH%dx )

       !WRITE (UNIT = 6, FMT = *)'Pour qx,', 'dx=', dx, 'Error=', ERR%L2_qx
       !WRITE (UNIT = 6, FMT = *)'Pour qx,', LOG(dx), LOG(ERR%L2_qx)

       ! ====== IN FILE =========
       WRITE (UNIT = 200, FMT = *)'Pour qx,', 'dx=', dx, 'Error=', ERR%L2_qx
       WRITE (UNIT = 200, FMT = *)'Pour qx,', LOG(dx), LOG(ERR%L2_qx)


       !! L1

       !! For h
       sum =  cte__0_d_r__8
       sum_init =  cte__0_d_r__8

       DO i= 1, Nbx
          sum_init = sum_init + ABS(TabEx(i)%h)
       END DO

       DO i= 1, Nbx
          sum = sum + ABS(Tab(i)%h-TabEx(i)%h)
       END DO

       !CError = sum/sum_init
       ERR%L1_h = sum*MESH%dx

       !WRITE(UNIT=6, FMT=*)'====== L1-error ======'
       !WRITE (UNIT = 6, FMT = *)'Pour h,', 'dx=', dx, 'Error=', ERR%L1_h
       !WRITE (UNIT = 6, FMT = *)'Pour h,', LOG(dx), LOG(ERR%L1_h)

       ! ====== IN FILE =========
       WRITE(UNIT=200, FMT = *)'====== L1-error ======'
       WRITE (UNIT = 200, FMT = *)'Pour h,', 'dx=', dx, 'Error=', ERR%L1_h
       WRITE (UNIT = 200, FMT = *)'Pour h,', LOG(dx), LOG(ERR%L1_h)


       !! For qx
       sum =  cte__0_d_r__8
       sum_init =  cte__0_d_r__8

       DO i= 1, Nbx
          sum_init = sum_init + ABS(TabEx(i)%qx)
       END DO

       DO i= 1, Nbx
          sum = sum + ABS(Tab(i)%qx-TabEx(i)%qx)
       END DO

       !CError = sum/sum_init
       ERR%L1_qx = sum*MESH%dx

       !WRITE (UNIT = 6, FMT = *)'Pour qx,', 'dx=', dx, 'Error=', ERR%L1_qx
       !WRITE (UNIT = 6, FMT = *)'Pour qx,', LOG(dx), LOG(ERR%L1_qx)

       ! ====== IN FILE =========
       WRITE (UNIT = 200, FMT = *)'Pour qx,', 'dx=', dx, 'Error=', ERR%L1_qx
       WRITE (UNIT = 200, FMT = *)'Pour qx,', LOG(dx), LOG(ERR%L1_qx)


       !! LINF
       !! For h
       ERR%LINF_h = MAXVAL(ABS(Tab(1:Nbx)%h-TabEx(1:Nbx)%h))
       !WRITE(UNIT=6,FMT=*)'====== LINF-error ======'
       !WRITE (UNIT = 6, FMT = *)'Pour h,', 'dx=', dx, 'Error=', ERR%LINF_h
       !WRITE (UNIT = 6, FMT = *)'Pour h,', LOG(dx), LOG(ERR%LINF_h)

       ! ====== IN FILE =========
       WRITE(UNIT=200, FMT = *)'====== LINF-error ======'
       WRITE (UNIT = 200, FMT = *)'Pour h,', 'dx=', dx, 'Error=', ERR%LINF_h
       WRITE (UNIT = 200, FMT = *)'Pour h,', LOG(dx), LOG(ERR%LINF_h)


       !! For qx
       ERR%LINF_qx = MAXVAL(ABS(Tab(1:Nbx)%qx-TabEx(1:Nbx)%qx))
       !WRITE (UNIT = 6, FMT = *)'Pour qx,', 'dx=', dx, 'Error=', ERR%LINF_qx
       !WRITE (UNIT = 6, FMT = *)'Pour qx,', LOG(dx), LOG(ERR%LINF_qx)

       ! ====== IN FILE =========
       WRITE (UNIT = 200, FMT = *)'Pour qx,', 'dx=', dx, 'Error=', ERR%LINF_qx
       WRITE (UNIT = 200, FMT = *)'Pour qx,', LOG(dx), LOG(ERR%LINF_qx)


    END SELECT

  END SUBROUTINE Compute_Errors

  SUBROUTINE Solution_Exact(DATA, VEx, Time)

    !USE fortranf
    USE modCST

    IMPLICIT NONE
    TYPE(DATA_IN), INTENT(IN) :: DATA
    INTEGER :: Nbx
    !INTEGER, INTENT ( IN ) :: N
    real*8 :: g = 9.8d0
    REAL(KIND = 8 ), INTENT ( IN ) :: Time
    !REAL(KIND = 8 ) :: x0
     REAL(KIND = 8 ) :: Ksi
     REAL(KIND = 8 ) :: U
     !REAL(KIND = 8 ) :: H
     REAL(KIND = 8 ) :: l
     REAL(KIND = 8 ) :: X
     REAL(KIND = 8 ) :: D
    !REAL(KIND = 8 ) :: sigmat
    !REAL(KIND = 8 ) :: V
    !REAL(KIND = 8 ) :: Y
    INTEGER :: i
    !INTEGER :: j
    !INTEGER :: ind
    REAL(KIND = 8 ) :: A
    !REAL(KIND = 8 ) :: S
    ! REAL(KIND = 8 ) :: Hij
    ! REAL(KIND = 8 ) :: nu
     REAL(KIND = 8 ) :: h0
    ! REAL(KIND = 8 ) :: Om
    ! REAL(KIND = 8 ) :: Amp
    ! REAL(KIND = 8 ) :: r
    ! REAL(KIND = 8 ) :: r0
     REAL(KIND = 8 ) :: atime
     REAL(KIND = 8 ) :: c0
    REAL(KIND = 8 ) :: t0
    REAL(KIND = 8 ) :: Xadim
    !REAL(KIND = 8 ) :: Yadim
    !REAL(KIND = 8 ) :: epsilon
    REAL(KIND = 8 ), DIMENSION(:), POINTER :: ZF
    TYPE( STATE ), DIMENSION(:), POINTER :: VEx
    REAL(KIND = 8 ), DIMENSION(:), POINTER :: Cell
    CHARACTER(LEN =* ), PARAMETER :: sub__name = "solution_exact"    
     REAL(kind = 8) :: k
     REAL(kind = 8) :: z0
     REAL(kind = 8) :: c
     REAL(KIND = 8 ) :: m0
    REAL(KIND = 8 ) :: u0
     REAL(KIND = 8 ) :: href
     REAL(KIND = 8 ) ::hex , ex , temp_x

    IF(debug) THEN
       !PRINT*, "Entering ", sub__name
    END IF

    Nbx = DATA%Nbx
 
    SELECT CASE( DATA%NCI)


    CASE(29)

       ! Soliton sur fond plat

       !Soliton sans topo - Notations code scilab
       ! Test de Florent - Premier test de validation précise

     
       H0 = 0.50D0
       A = 0.2D0*H0
       k = sqrt(3.0D0*a/(4.0D0*H0**2*(H0+A)))
       c = sqrt(g_r__8*(H0+a))
       L = DATA%x_out - DATA%x_in
       d = L/2.0D0
       
       DO i = -2, Nbx+3

          X = MESH%X_VF(i) !)(MESH%X_DF(i) + MESH%X_DF(i+1))/2.0D0
          !Z0 = X - INT((x-d+L/3.0D0)/L)*L

          Ksi = a/(COSH(k*(X-c*Time - d))**2)

          U = c*(1.0D0 -H0/(Ksi+H0))

       
             VEx(i)%h = H0 + ksi
             VEx(i)%qx = U*(VEx(i)%h)

       END DO
      
    CASE(290)

       ! Soliton sur fond plat
       ! Config Delis

     
       H0 = 10.0D0
       A = 0.2D0*H0
       k = sqrt(3.0D0*a/(4.0D0*H0**2*(H0+A)))
       c = sqrt(g_r__8*(H0+a))
       L = DATA%x_out - DATA%x_in
       d = 200.0D0
       
       DO i = -2, Nbx+3

          X = MESH%X_VF(i) !)(MESH%X_DF(i) + MESH%X_DF(i+1))/2.0D0
          !Z0 = X - INT((x-d+L/3.0D0)/L)*L

          Ksi = a/(COSH(k*(X-c*Time - d))**2)

          U = c*(1.0D0 -H0/(Ksi+H0))

       
             VEx(i)%h = H0 + ksi
             VEx(i)%qx = U*(VEx(i)%h)

       END DO
       
    CASE( 0 ) ! cte_1 avant
       !Cas Dam-break sur fond sec en dimension 1 et adimensioné
       h0 = 1.0 !cte__5_d_0
       c0 = 1.0 !SQRT (g * h0 )
       t0 = 1.0 !c0 / g
       !PRINT*,cell(i)
       !WRITE (UNIT = 6, FMT =*) "Time", Time
       atime = Time / t0
       X =  (MESH%X_VF(i) +  MESH%X_VF(i + 1))/ cte__2_d_0_r__8 !( Cell(i) + Cell(i+1) ) / cte__2_d_0_r__8
       Xadim = X / h0


       DO i = 1, Nbx
          X =  (MESH%X_VF(i) +  MESH%X_VF(i + 1))/ cte__2_d_0_r__8    !( Cell(i) + Cell(i+1) ) / cte__2_d_0_r__8
          Xadim = X / h0
          IF (Xadim / atime<= cte__2_d_0_r__8) THEN
             Ksi =( 1.0D0 / cte__9_d_0_r__8) * ( cte__2_d_0_r__8 - Xadim / atime ) ** 2
          ELSE
             Ksi = 0.0D0
          END IF
          U =( cte__2_d_0_r__8 / cte__3_d_0_r__8 ) *( Xadim / atime + 1.0D0 )
          VEx(i)%h = h0 * Ksi
          VEx(i)%qx = VEx(i)%h * c0 * U

       END DO
  


    ! CASE(cte__18)

    !    ! tidal wave flow : very accurate asymptotically derived approximation
    !    ! Bermudez & Vazquez - Computer and fluids 1994

    !    L = 14000.0D0

    !    DO i = 1, Nbx
        
    !          X =( Cell(i) + Cell(i+1) ) / 2.0D0

    !          TabEx(i)%h = 64.5D0 - ZF(i) - 4.0D0*SIN(Pi_r__8*(4*Time/86400.0D0 + 1.0D0/2.0D0))
    !          TabEx(i)%qx = ((X-L)*Pi_r__8/5400.0D0)*COS(Pi_r__8*(4*Time/86400.0D0 + 1.0D0/2.0D0))
           

    !       END DO 


 
    ! CASE ( cte__40)
    !    !Onde linéaire 1D sur fond plat - Gaussienne
    !    x0 = cte__25_d_0
    !    sigmat = cte__5_d_0
    !    D = cte__10_d_0 / g
    !    H = cte__0_d_01_r__8 * D

    !    loop__19 : DO i = -2, Nbx+3
      
    !          X =( Cell ( i) + Cell (i+1) ) / cte__2_d_0_r__8
    !          Ksi = H * EXP (( -(( ( X - x0 ) - SQRT (( g * D) ) * Time) / sigmat) ** 2))
    !          U = Ksi * SQRT (( g / D ) )
    !          !print*,i,j,'Ksi', Ksi,'U',U
    !          !ZF(i) =  cte__0_d_r__8
    !          TabEx(i)%h= D + Ksi
    !          TabEx(i)%qx= U * TabEx(i)%h
            
         
    !    END DO loop__19
!!$
!!$    ! CASE(cte__12)
!!$    !    ! Ecoulement sub-critical
!!$
!!$    !    IF(N <= 1) THEN
!!$
!!$    !       ! Computation of the stationnary solution at the first iteration
!!$    !       !print*,'Je calcul la sol ex'
!!$    !       m0 = 4.42D0
!!$    !       h0 = cte__2_d_0_r__8
!!$    !       u0 = m0/h0
!!$
!!$       
!!$    !          DO i = -2, Nbx+3
!!$
!!$    !             href = 2.0D0
!!$    !             !ex = 1.0D0
!!$
!!$    !             hex = (0.5D0*u0**2 + g_r__8*h0 -g_r__8*ZF(i) - 0.5D0*m0**2/href**2)/g_r__8
!!$
!!$    !             DO WHILE((ABS(hex-href) >= 1.0D-17))
!!$    !                href = hex
!!$    !                hex = (0.5D0*u0**2 + g_r__8*h0 -g_r__8*ZF(i) - 0.5D0*m0**2/href**2)/g_r__8nnn
!!$    !             END DO
!!$
!!$    !             TabEx(i)%h = hex
!!$    !             TabEx(i)%qx = m0
!!$             
!!$    !             !PRINT*,'TabEx(i)%h',TabEx(i)%h
!!$           
!!$    !       END DO
!!$       !    END IF

    CASE(1996)
       !  Ecoulement sub-critical


       ! Computation of the stationnary solution at the first iteration
       !print*,'Je calcul la sol ex'
       temp_x = 0.0d0
       m0 = 4.42D0
       h0 = cte__2_d_0_r__8
       u0 = m0/h0
       do i = -2 , Nbx+3
          if( temp_x < 12 .and. temp_x > 8 ) then

             VEx(i)%topo =  0.2d0 - 0.05d0*(temp_x - 10)**2
          else
             VEx(i)%topo = 0.0d0
          end if
       end do

       temp_x = temp_x + 0.083

       DO i = -2, Nbx+3

          href = 2.0D0
          ex = 1.0D0
      
          hex = (0.5D0*u0**2 + 9.8d0*h0 - 9.8d0*VEx(i)%topo - 0.5D0*m0**2/href**2)/9.8d0

          DO WHILE((ABS(hex-href) >= 1.0D-17))
             href = hex
             hex = (0.5D0*u0**2 + 9.8d0*h0 - 9.8d0*VEx(i)%topo - 0.5D0*m0**2/href**2)/9.8d0
          END DO

          VEx(i)%h = hex
          VEx(i)%qx = m0

          !PRINT*,'TabEx(i)%h',VEx(i)%h

       END DO

    ! CASE(cte__26)
    !    ! sub-critical flow over a bump, accuracy version
    !    !print*,'ici'
    !    !read*
    !    IF(N <= 1) THEN

    !       ! Computation of the stationnary solution at the first iteration

         
    !          DO i = -2, Nbx+3

    !             CALL romberg(XX = MESH%X_VF(i), S_OUT = TabEx(i), FCT = ANALYTIC_SOL)
 
    !          END DO
         
    !    END IF


      
    END SELECT
  END SUBROUTINE Solution_Exact

 SUBROUTINE TIME_ERROR(Run, ERR)

!============================================================
!     Time evolution of the L2 error for a given mesh (time convergence study)
!============================================================
    !USE fortranf
    USE modCST

    IMPLICIT NONE
    TYPE(T_ERR), INTENT(IN) :: ERR
     TYPE(MyRun), INTENT(IN) :: Run

             ! ====== IN FILE =========
          WRITE(UNIT=239, FMT=*)log(Run%time), log(ERR%L2_h), log(ERR%L2_qx)

  END SUBROUTINE TIME_ERROR


END MODULE modPOST
