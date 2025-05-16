!=========================================================
  ! subroutines and functions for output and visualization
!=========================================================

MODULE modOutputs

  USE modTOOLS
  USE modPOST

  PRIVATE

  PUBLIC :: OUTPUT
  PUBLIC :: Out_1D
  PUBLIC :: OUT_TIME_SERIES

CONTAINS


  SUBROUTINE OUTPUT(N, MESH, V, DATA, Run, FILE_STATE)

    !USE fortranf
    USE modCST

    IMPLICIT NONE

    TYPE(Variable), INTENT(INOUT) :: V
    TYPE(MyRun), INTENT(IN) :: Run
    TYPE(MeshDef), INTENT(IN) :: MESH
    TYPE(DATA_IN), INTENT(IN) :: DATA
    CHARACTER (LEN=*), INTENT(IN) :: FILE_STATE
    !CHARACTER (LEN=* ), PARAMETER :: c__format__1 = "(5ES16.8)"
    !CHARACTER (LEN=* ), PARAMETER :: c__format__2 = "(9ES16.8)"
    !CHARACTER (LEN=* ), PARAMETER :: c__format__3 = "(7ES16.8)"
    !CHARACTER (LEN=* ), PARAMETER :: c__format__4 = "(6ES16.8)"
    !CHARACTER (LEN=* ), PARAMETER :: c__format__5 = "(10ES16.8)"
    INTEGER, INTENT ( IN ) :: N
    REAL(kind = 8), DIMENSION(-2:DATA%Nbx+3) :: tU

    INTEGER :: Nbx
    INTEGER :: i

    !Computation of the exact solution if $Solex$=TRUE

    Nbx = DATA%Nbx

    IF(DATA%solex) THEN
       !WRITE(UNIT = 6, FMT = *)'Computing analytical solution'
       CALL Solution_Exact(DATA, V%Ex,Run%time)
    END IF

    ! Preparing results for output...drying process
 
       DO i = -2, Nbx+3
          
          IF (V%VF(i)%h <= DATA%kappa) THEN
             V%VF(i)%h  =  DATA%kappa
             tU(i)  =0.0D0
             V%VF(i)%qx =0.0D0
          
          ELSEIF(V%VF(i)%h <= DATA%Kappabis) THEN
             tU(i)  = 0.0D0
             V%VF(i)%qx = 0.0D0
           
          ELSE
             tU(i)  = V%VF(i)%qx/ V%VF(i)%h 
          END IF

          IF (ABS(tU(i)) <= DATA%Utol) tU(i) = 0.0D0
        
       END DO

    CALL Out_1D(V%VF, MESH, tU, V%topoVF, Nbx, DATA%solex, V%Ex, DATA%NCI, N, FILE_STATE)

  END SUBROUTINE OUTPUT

  
  SUBROUTINE Out_1D(Tab, MESH, tU, ZF, Nbx,solex, TabEx, NCI, N, FILE_STATE)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: Nbx
    CHARACTER (LEN=*), INTENT(IN) :: FILE_STATE
    TYPE(STATE), DIMENSION(:), POINTER :: Tab
    TYPE(STATE), DIMENSION(:), POINTER :: TabEx
        TYPE(MeshDef), INTENT(IN) :: MESH
    REAL(kind = 8), DIMENSION(:), POINTER :: ZF
    LOGICAL, INTENT(in), OPTIONAL:: SolEx

    REAL(kind = 8), DIMENSION(-2:Nbx+3), INTENT(in) :: tU
    INTEGER, INTENT(in) :: NCI
    INTEGER, INTENT(in) :: N
    INTEGER :: i
    !INTEGER :: j
    !INTEGER :: k
    INTEGER :: ios
    !INTEGER :: m
    CHARACTER (len=30) :: dir
    !CHARACTER (len=30) :: str1
    REAL(kind = 8) :: c0
    REAL(kind = 8) :: h0
    REAL(kind = 8) :: sH
    REAL(kind = 8) :: sU
    REAL(kind = 8) :: sQ
    REAL(kind = 8) :: ZZ
    REAL(kind = 8) :: XX
    !REAL(kind = 8) :: YY
    !REAL(kind = 8) :: rr
    REAL(kind = 8) :: ll
    REAL(kind = 8) :: ss
    REAL(kind = 8) :: epsi
    REAL(kind = 8) :: sHEx
    REAL(kind = 8) :: sUEx
    REAL(kind = 8) :: sQEx


    CHARACTER (LEN=* ), PARAMETER :: c__format__9 = "(9ES16.8)"
    CHARACTER (LEN=* ), PARAMETER :: c__format__10 = "(10ES16.8)"
    CHARACTER (LEN=* ), PARAMETER :: c__format__6 = "(6ES16.8)"
    CHARACTER (LEN=* ), PARAMETER :: c__format__11 = "(11ES16.8)"

    dir = 'res/'

    !CALL Ouvrir(forout,TRIM(dir)//"SOL1D."// TRIM(carac(N))//".res", ios)
    CALL Ouvrir(forout,TRIM(dir)//TRIM(ADJUSTL(FILE_STATE))//TRIM(carac(N))//".res", ios)
    IF (ios /= 0) THEN
       WRITE(UNIT = 6, FMT =*) "File not found. Default values will be used."
       STOP "STOP OUT_1D"
    END IF


    DO i = -2, Nbx+3!, NxStep
       !print*,'i',i

       sH = Tab(i)%h
       sU = tU(i)  
       
       sQ = Tab(i)%qx
       
       ZZ = ZF(i)
       
       XX = MESH%X_VF(i)
       !print*,"invariant = ","   ",tU(20) + sqrt(9.8*Tab(20)%h)
   
       SELECT CASE(SolEx)

       CASE(.FALSE.)
          !print*,'Solex=FALSE'
          ! Output Format
          !#################################################
          !## 1:X  ## 2:H ## 3:U ## 4:Zf ## 5:H+Zf ## 6:Q  #
          !#################################################

          SELECT CASE(NCI)



          CASE(cte__36) ! C&G transient adimensionne
             ll = cte__20_d_0_r__8
             ss = 1.0D0 / cte__50_d_0_r__8
             epsi = cte__0_d_1_r__8
             WRITE (UNIT = forout, FMT = c__format__6)XX/ll, sH/(epsi*ss*ll), sU, ZZ/(epsi*ss*ll), (sH+ZZ)/(epsi*ss*ll), sQ


          CASE(cte__37) ! C&G periodic adim
             ll = cte__30_d_0_r__8
             SS = 1.0D0 / cte__30_d_0_r__8
             WRITE (UNIT = forout, FMT = c__format__6)XX/ll, sH/(ss*ll), sU, ZZ/(ss*ll), (sH+ZZ)/(ss*ll), sQ

          CASE DEFAULT ! Cas general

             !print*,'je passe ici'
              !print*,'xx',xx, i
             WRITE (UNIT = forout, FMT = c__format__6)XX, sH, sU, ZZ, sH+ZZ, sQ 
             !WRITE (UNIT = forout, FMT = c__format__6)XX,sU, sqrt(9.8*sH) , -sqrt(9.8*sH)
          

          END SELECT

       CASE(.TRUE.)
     
          IF (TabEX(i)%h.LE.DATA%kappa) THEN
             sHEx = kappa
             sUEx = cte__0_d_r__8  
             sQEx = cte__0_d_r__8
          ELSE
             sHEx = TabEx(i)%h
             sUEx = TabEx(i)%qx/ sHEx 
             sQEx = TabEx(i)%qx
          END IF

          ! Output format 1D with exact solution
          ! ## 1:X  ## 2:H ## 3:U ## 4:Zf ## 5:H+Zf ## 6:HEx+ZF ## 7:UEx ## 8:HEx

          SELECT CASE(NCI)

          CASE(cte__1)
             h0 = cte__5_d_0_r__8
             c0 = SQRT(g_r__8*h0)
             WRITE (UNIT = forout, FMT = c__format__11)XX/h0, sH/h0, sU/c0, ZZ/h0, (sH+ZZ)/h0, (sHEx+ZZ)/h0, &
                  sUEx/c0, sQ/(c0*h0), sQEx/(c0*h0)


          CASE(0) ! cte_0
             h0 = cte__5_d_0_r__8
             c0 = SQRT(g_r__8*h0)
             WRITE (UNIT = forout, FMT = c__format__11)XX/h0, sH/h0, sU/c0, ZZ/h0, (sH+ZZ)/h0, (sHEx+ZZ)/h0, &
                  sUEx/c0, sQ/(c0*h0), sQEx/(c0*h0)


          CASE(cte__36)
             ll = cte__20_d_0_r__8
             Ss = 1.0D0 / cte__50_d_0_r__8
             epsi = cte__0_d_1_r__8
             WRITE (UNIT = forout, FMT = c__format__9)XX/ll, sH/(epsi*ss*ll), sU, &
                  ZZ/(epsi*ss*ll), (sH+ZZ)/(epsi*ss*ll), (sHEx+ZZ)/(epsi*ss*ll), sUEx
 
          CASE(cte__37)
             ll = cte__20_d_0_r__8
             SS = 1.0D0 / cte__30_d_0_r__8
             WRITE (UNIT = forout, FMT = c__format__9)XX/ll, sH/(ss*ll), sU, ZZ/(ss*ll), (sH+ZZ)/(ss*ll), (sHEx+ZZ)/(ss*ll), sUEx


          CASE DEFAULT
             WRITE (UNIT = forout, FMT = c__format__10)XX, sH, sU, ZZ, sH+ZZ, sHEx+ZZ, sUEx, sHEx, sQ, sQEx

          END SELECT

       END SELECT

    END DO

    CALL Fermer(forout)

  END SUBROUTINE Out_1D


  SUBROUTINE Out_1D_DF(Tab, X_DF, tU, ZF, Ndf, N, FILE_STATE)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: Ndf
    CHARACTER (LEN=*), INTENT(IN) :: FILE_STATE
    TYPE(STATE), DIMENSION(:), POINTER :: Tab
    REAL(kind = 8), DIMENSION(:), POINTER :: X_DF
    REAL(kind = 8), DIMENSION(:), POINTER :: ZF
    REAL(kind = 8), DIMENSION(:), INTENT(in) :: tU

    INTEGER, INTENT(in) :: N
    INTEGER :: i
    !INTEGER :: j
    !INTEGER :: k
    INTEGER :: ios
    !INTEGER :: m
    CHARACTER (len=30) :: dir
    !CHARACTER (len=30) :: str1
    !REAL(kind = 8) :: c0
    !REAL(kind = 8) :: h0
    REAL(kind = 8) :: sH
    REAL(kind = 8) :: sU
    REAL(kind = 8) :: sQ
    REAL(kind = 8) :: ZZ
    REAL(kind = 8) :: XX
    !REAL(kind = 8) :: YY
    !REAL(kind = 8) :: rr
    !REAL(kind = 8) :: ll
    !REAL(kind = 8) :: ss
    !REAL(kind = 8) :: epsi
    !REAL(kind = 8) :: sHEx
    !REAL(kind = 8) :: sUEx
    !REAL(kind = 8) :: sQEx


    !CHARACTER (LEN=* ), PARAMETER :: c__format__9 = "(9ES16.8)"
    !CHARACTER (LEN=* ), PARAMETER :: c__format__10 = "(10ES16.8)"
    CHARACTER (LEN=* ), PARAMETER :: c__format__6 = "(6ES16.8)"
    !CHARACTER (LEN=* ), PARAMETER :: c__format__11 = "(11ES16.8)"
REAL(KIND = 8) :: eta

    dir = 'RESULTS/'

    !CALL Ouvrir(forout,TRIM(dir)//"SOL1D."// TRIM(carac(N))//".res", ios)
    CALL Ouvrir(forout,TRIM(dir)//TRIM(ADJUSTL(FILE_STATE))//TRIM(carac(N))//".res", ios)
    IF (ios /= 0) THEN
       WRITE(UNIT = 6, FMT =*) "File not found. Default values will be used."
       STOP "STOP OUT_1D"
    END IF

    DO i = 1, Ndf
 
       sH = Tab(i)%h
       sU = tU(i)  
       sQ = Tab(i)%qx
       ZZ = ZF(i)
       eta = Tab(i)%xi


       XX = X_DF(i)

             WRITE (UNIT = forout, FMT = c__format__6)XX, sH, sU, ZZ, eta, sQ


    END DO

    CALL Fermer(forout)

  END SUBROUTINE Out_1D_DF






END MODULE modOUTPUTS
