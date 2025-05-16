MODULE modTOOLS

  !USE fortranf
  USE modCST
  USE modTypes
  USE modCommon

  IMPLICIT NONE

  PUBLIC:: diag
  
  PUBLIC:: Initialise
  PUBLIC:: Ouvrir
  PUBLIC:: Fermer
  PUBLIC:: carac
  PUBLIC:: Solex_CG_Transient
  PUBLIC:: Solex_CG_Periodic
  PUBLIC:: TimeStep
  PUBLIC:: bessj1
  PUBLIC:: bessj0
  PUBLIC:: PSI
  PUBLIC:: Find_Last_Wet
  PUBLIC:: Find_Last_Wet_Reverse
  PUBLIC:: DEALLOCATE_ALL
  PUBLIC :: read__stop
  !PUBLIC :: inv

CONTAINS

  FUNCTION ROUND(r)

    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: r
    INTEGER :: ROUND
    INTEGER :: intm
    INTEGER :: intp


    intm = FLOOR(r)
    intp = FLOOR(r)+1

    IF(ABS(r-intm)>ABS(r-intp)) THEN
       ROUND = intp
    ELSE
       ROUND = intm
    END IF

  END FUNCTION ROUND


  SUBROUTINE read__stop(mssg, code) !-- Interface

    !CHARACTER(LEN=*), PARAMETER :: mod_name = "read__stop"

    !-- arguments

    CHARACTER(LEN=*), INTENT(IN) :: mssg
    INTEGER, INTENT(IN) :: code

    !-- locales

    WRITE(UNIT = 0, FMT=*) "Erreur ", mssg, ", code==", code
    STOP

  END SUBROUTINE read__stop

  PURE FUNCTION diag(diago, loc)

    ! ATTENTION, aucun garde-fou, le vecteur 'diago' doit etre Ã  la bonne taille pour loger a la diagonale 'loc'

    !real(kind=8), dimension(:,:), intent( INOUT ) :: Mat
    REAL(KIND = 8), DIMENSION(:), INTENT( IN ) :: diago
    INTEGER, INTENT(in) :: loc
    REAL(KIND = 8), DIMENSION(SIZE(diago)+ABS(loc), SIZE(diago)+ABS(loc)) :: diag
    INTEGER :: i
    !INTEGER :: istat

    diag = 0.0D0

    FORALL( i = 1: SIZE(diago))
       diag(i,i+ABS(loc)) = diago(i)
    END FORALL
    IF(loc < 0) diag = TRANSPOSE(diag)

  END FUNCTION diag

  ! PURE FUNCTION ACOSH( X) RESULT( acosh__res)
  !   !USE fortranf
  !   !USE modCST
  !   REAL(kind = 8) :: f__0
  !   INTERFACE
  !      PURE FUNCTION R1MACH( arg__0) RESULT( r1mach__res)
  !        !USE fortranf
  !        INTEGER, INTENT ( IN ) :: arg__0
  !        REAL(KIND=8) :: r1mach__res
  !      END FUNCTION R1MACH
  !   END INTERFACE
  !   REAL(KIND = 8) :: acosh__res
  !   REAL, PARAMETER :: ALN2 = 0.69314718055994530942e0
  !   REAL(KIND = 8) :: XMAX
  !   REAL(KIND = 8), INTENT ( IN ) :: X

  !   !***FIRST EXECUTABLE STATEMENT ACOSH
  !   f__0 = R1MACH( cte__3)
  !   XMAX = 1.0 / SQRT (f__0)

  !   IF (X < 1.0 ) THEN
  !      acosh__res = -9.99D99 !-- erreur
  !   ELSE

  !      IF (X < XMAX ) THEN
  !         acosh__res = LOG (X + SQRT(( X * X - 1.0)))
  !      ELSE
  !         acosh__res = ALN2 + LOG (X)
  !      END IF
  !   END IF

  ! END FUNCTION ACOSH

  ! SUBROUTINE DEALLOC_CSR(in_)


  !   TYPE(CSR), INTENT(INOUT) :: in_
  !   INTEGER :: istat
  !   IF(ALLOCATED(in_%a))  DEALLOCATE(in_%a, STAT = istat )
  !   IF(ALLOCATED(in_%ia))  DEALLOCATE(in_%ia, STAT = istat )
  !   IF(ALLOCATED(in_%ja))  DEALLOCATE(in_%ja, STAT = istat )



  ! END SUBROUTINE DEALLOC_CSR


  SUBROUTINE DEALLOCATE_ALL()

    IMPLICIT NONE
    !INTEGER :: NUM
    !INTEGER :: i__stat
    INTEGER :: istat
    !INTEGER :: ii

    CHARACTER(len=*), PARAMETER :: mod__name = 'DEALLOCATE_ALL'


    IF(ASSOCIATED(FLU%ALPHA))  DEALLOCATE(FLU%ALPHA, STAT = istat )
    IF (debug) THEN
       !PRINT *, mod__name, " DEALLOCATE FLU%ALPHA ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err FLU%ALPHA"
    END IF

    IF(ASSOCIATED(FLU%SR))  DEALLOCATE(FLU%SR, STAT = istat )
    IF (debug) THEN
       !PRINT *, mod__name, " DEALLOCATE FLU%SR ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err FLU%SR"
    END IF

    IF(ASSOCIATED(MESH%array_DX))  DEALLOCATE(MESH%array_DX, STAT = istat )
    IF (debug) THEN
       !PRINT *, mod__name, " DEALLOCATE array_DX ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err array_DX"
    END IF

    IF(ASSOCIATED(DATA%coor))  DEALLOCATE(DATA%coor, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE  DATA%coor", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err DATA%coor"
    END IF


    IF(ASSOCIATED(Var%t1)) DEALLOCATE(Var%t1, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE Tab1 ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err Tab1"
    END IF

    IF(ASSOCIATED(Var%t2)) DEALLOCATE(Var%t2, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE Var%t2 ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err Var%t2"
    END IF

    IF(ASSOCIATED(Var%t3)) DEALLOCATE(Var%t3, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE Var%t3 ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err Var%t3"
    END IF

    IF(ASSOCIATED(Var%t4)) DEALLOCATE(Var%t4, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE Var%t3 ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err Var%t3"
    END IF

    IF(ASSOCIATED(Var%t0)) DEALLOCATE(Var%t0, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE Var%t0 ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err Var%t0"
    END IF

    IF(ASSOCIATED(Var%tmp)) DEALLOCATE(Var%tmp, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE Tabtmp ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err Tamtmp"
    END IF

    IF(ASSOCIATED(Var%ref)) DEALLOCATE(Var%Ref, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE TabRef ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err TabRef"
    END IF

    IF(ASSOCIATED(FLU%X1)) DEALLOCATE(FLU%X1, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE FLU_X1 ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err FLU_X1"
    END IF

    IF(ASSOCIATED(FLU%X2)) DEALLOCATE( FLU%X2, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE FLU_X2 ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err FLU_X2"
    END IF

    IF(ASSOCIATED(FLU%Sx)) DEALLOCATE( FLU%Sx, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE FLU_Sx ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err FLU_Sx"
    END IF

    IF(ASSOCIATED(MESH%X_VF)) DEALLOCATE(MESH%X_VF, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE X_VF ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err X_VF"
    END IF

    IF(ASSOCIATED(MESH%X_DF)) DEALLOCATE(MESH%X_DF, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE X_DF ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "alloc_err X_DF"
    END IF

    IF(ASSOCIATED(Var%Ex)) DEALLOCATE(Var%Ex, STAT = istat )
    IF (debug ) THEN
       !PRINT *, mod__name, " DEALLOCATE TabEx ", istat
    END IF
    IF (istat /= 0 ) THEN
       STOP "De-alloc_err TabEx"
    END IF

    IF(ASSOCIATED(Var%Vf)) DEALLOCATE (Var%Vf, STAT = istat)
    IF (debug) THEN
       !PRINT *, mod__name, " DEALLOCATE Var%Vf"
    END IF
    IF (istat /= 0 ) THEN
       STOP "dealloc_err Var%Vf"
    END IF

    IF(ASSOCIATED(Var%topoVf)) DEALLOCATE (Var%topoVf, STAT = istat)
    IF (debug) THEN
       !PRINT *, mod__name, " DEALLOCATE Var%topoVf"
    END IF
    IF (istat /= 0 ) THEN
       STOP "dealloc_err Var%topoVf"
    END IF



    !PRINT*, 'All Tabs are de-allocated now, End of computation.'

  END SUBROUTINE DEALLOCATE_ALL

  SUBROUTINE Ouvrir( unit_, file_, iostat_)
    !USE fortranf
    USE modCST

    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: unit_
    CHARACTER(LEN =*), INTENT ( IN ) :: file_
    INTEGER, INTENT ( OUT ) :: iostat_
    CHARACTER(LEN = 5) :: A
    CHARACTER(LEN = 7) :: S
    CHARACTER(LEN = 11) :: F
    INTEGER :: simple_iostat
    !CHARACTER(LEN =*), PARAMETER :: mod__name = "ouvrir"

    ! IF(imprimer_la_trace ) THEN
    !    PRINT*, "Entering ", mod__name
    ! END IF

    SELECT CASE( unit_)
    CASE( forin, forin2, forin9 )
       F = "FORMATTED  "
       S = "OLD"
       A = "READ "

    CASE ( forout, forout2, forout9 )
       F = "FORMATTED  "
       S = "REPLACE"
       A = "WRITE"

    CASE ( unfin)
       F = "UNFORMATTED"
       S = "OLD"
       A = "READ "

    CASE ( unfout)
       F = "UNFORMATTED"
       S = "UNKNOWN"
       A = "WRITE"

    CASE DEFAULT
       iostat_ = cte__1000
       GO TO 99998
    END SELECT

    !-----------------------------------------------------------------------
    !                                                  -- Ouvrir.7
    !-----------------------------------------------------------------------
    !IF(env%verb >= 10) WRITE(forout9,*) "---- ", a, " ", TRIM(file_)
    ! certaines machines ne respectent pas la norme "any trailing blanks are ignored"

    OPEN (UNIT = unit_, FILE = TRIM( file_ ), STATUS = TRIM(S), FORM = TRIM (F), &
         ACTION = TRIM(A), IOSTAT = simple_iostat )

    iostat_ = simple_iostat

99998 CONTINUE

  END SUBROUTINE Ouvrir


  SUBROUTINE Fermer( std)
    !USE fortranf
    USE modCST

    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: std


    CLOSE (UNIT = std)

  END SUBROUTINE Fermer

  FUNCTION carac(i)
    CHARACTER*8 :: carac
    INTEGER :: i
    WRITE (carac, "(i8)") i
    carac = ADJUSTL(carac)
  END FUNCTION carac





  SUBROUTINE Initialise(i_case, DATA, Run)

    !USE fortranf
    USE modCST

    IMPLICIT NONE

    integer, intent( IN ) :: i_case
    TYPE(MyRun), INTENT(OUT) :: Run
    TYPE(DATA_IN), INTENT(OUT) :: DATA
    CHARACTER(LEN =* ), PARAMETER :: mod__name = "initialise"
    CHARACTER (len=30) :: dir

    INTEGER :: statut
    CHARACTER(len=30) :: bla
    INTEGER :: ii


    DATA%NCI = i_case

    select case(i_case)

    case(0) ! dam break on wet bottom

       DATA%NCLL = 2
       DATA%NCLR = 2
       DATA%cond_lim = 4

       DATA%Nbx = 7500
       DATA%x_in = 0.D0
       DATA%x_out = 15.D0
       DATA%NtMax = 5000000
       DATA%TimeInit = 0.0D0

       DATA%TimeOut = 5.0D0
       DATA%limiter = 1

       DATA%space_order = 'godunov'
       DATA%Riemann_solv = 'HLLC'
       DATA%time_order = 1

       DATA%kappa = 1.0D-4
       DATA%kappabis = 5.0D-3

       DATA%Utol = 1.0D-9

       DATA%CFL = 1.0d-02 !0.5D0
       DATA%ChoixDt = .FALSE.
       DATA%DtFix = 2.0D-2
       DATA%Solex = .FALSE.
       DATA%sol_ref = .FALSE.
       DATA%char__SolRef = 'solution_exact.res'
       DATA%NbxRef = 6400
       DATA%var_mat = .TRUE.
       DATA%NFreq = 10000
       DATA%ERROR = .FALSE.
       DATA%err_type = 'L2'
       DATA%TIME_SERIE = .FALSE.
       DATA%NFreqTime = 10
       DATA%point_number = 2
       DATA%file_time = 'time_series'
    case(2) ! Generation sinusoidale

       DATA%NCLL = 2
       DATA%NCLR = 2
       DATA%cond_lim = 4

       DATA%Nbx = 500
       DATA%x_in = 0.D0
       DATA%x_out = 1.D0
       DATA%NtMax = 90000000
       DATA%TimeInit = 0.0D0

       DATA%TimeOut = 900.0D0
       DATA%limiter = 1

       DATA%space_order = 'godunov'
       DATA%Riemann_solv = 'HLLC'
       DATA%time_order = 1

       DATA%kappa = 1.0D-4
       DATA%kappabis = 5.0D-3

       DATA%Utol = 1.0D-9

       DATA%CFL = 1.0d-02 !0.5D0
       DATA%ChoixDt = .FALSE.
       DATA%DtFix = 2.0D-2
       DATA%Solex = .FALSE.
       DATA%sol_ref = .FALSE.
       DATA%char__SolRef = 'solution_exact.res'
       DATA%NbxRef = 6400
       DATA%var_mat = .TRUE.
       DATA%NFreq = 10000
       DATA%ERROR = .FALSE.
       DATA%err_type = 'L2'
       DATA%TIME_SERIE = .FALSE.
       DATA%NFreqTime = 10
       DATA%point_number = 2
       DATA%file_time = 'time_series'
    case(1) ! Generation soliton

       DATA%NCLL = 2
       DATA%NCLR = 2
       DATA%cond_lim = 4

       DATA%Nbx = 7500
       DATA%x_in = 0.D0
       DATA%x_out = 15.D0
       DATA%NtMax = 5000000
       DATA%TimeInit = 0.0D0

       DATA%TimeOut = 5.0D0
       DATA%limiter = 1

       DATA%space_order = 'godunov'
       DATA%Riemann_solv = 'HLLC'
       DATA%time_order = 1

       DATA%kappa = 1.0D-4
       DATA%kappabis = 5.0D-3

       DATA%Utol = 1.0D-9

       DATA%CFL = 1.0d-02 !0.5D0
       DATA%ChoixDt = .FALSE.
       DATA%DtFix = 2.0D-2
       DATA%Solex = .FALSE.
       DATA%sol_ref = .FALSE.
       DATA%char__SolRef = 'solution_exact.res'
       DATA%NbxRef = 6400
       DATA%var_mat = .TRUE.
       DATA%NFreq = 10000
       DATA%ERROR = .FALSE.
       DATA%err_type = 'L2'
       DATA%TIME_SERIE = .FALSE.
       DATA%NFreqTime = 10
       DATA%point_number = 2
       DATA%file_time = 'time_series'
    case(1964)
       DATA%NCLL = 2
       DATA%NCLR = 2
       DATA%cond_lim = 3 ! 3 pour Neumann réfléxive

       DATA%Nbx = 5000 
       DATA%x_in = 0.0D0
       DATA%x_out = 1000.0D0 

             
       DATA%NtMax = 200000000 !300
       DATA%TimeInit = 0.0D0
 
       DATA%TimeOut = 200.0D0
       DATA%limiter = 1

       DATA%space_order = 'godunov'
       DATA%Riemann_solv = 'HLLC'
       DATA%time_order = 1

       DATA%kappa = 1.0D-4
       DATA%kappabis = 5.0D-3

       DATA%Utol = 1.0D-9

       DATA%CFL = 1.0d-02 !0.5D0
       DATA%ChoixDt = .FALSE.
       DATA%DtFix = 2.0D-2
       DATA%Solex = .FALSE.
       DATA%sol_ref = .FALSE.
       DATA%char__SolRef = 'no file'
       DATA%NbxRef = 6400
       DATA%var_mat = .TRUE.
       DATA%NFreq = 10000
       DATA%ERROR = .FALSE.
       DATA%err_type = 'L2'
       DATA%TIME_SERIE = .FALSE.
       DATA%NFreqTime = 10
       DATA%point_number = 2
       DATA%file_time = 'time_series'

    case(1996)
       DATA%NCLL = 2
       DATA%NCLR = 2
       DATA%cond_lim = 1996 ! 3 pour Neumann reflexive

       DATA%Nbx = 500  !500 !12500
       DATA%x_in = 0.0d0 !0.0D0
       DATA%x_out =  1.0d0 !25.0d0
       DATA%NtMax =  2000000 ! 90000000
       DATA%TimeInit = 0.0D0

       DATA%TimeOut = 200.0D0 ! 900.0d0
       DATA%limiter = 1

       DATA%space_order = 'godunov'
       DATA%Riemann_solv = 'VFROE'
       DATA%time_order = 1

       DATA%kappa = 1.0D-4
       DATA%kappabis = 5.0D-3

       DATA%Utol = 1.0D-9

       DATA%CFL = 1.0d-02 !0.4D0
       DATA%ChoixDt = .FALSE.
       DATA%DtFix = 2.0D-2
       DATA%Solex = .FALSE.
       DATA%sol_ref = .TRUE.
       DATA%char__SolRef = 'solution_exact.res'
       DATA%NbxRef = 6400
       DATA%var_mat = .TRUE.
       DATA%NFreq = 10000 !10
       DATA%ERROR = .FALSE.
       DATA%err_type = 'L2'
       DATA%TIME_SERIE = .FALSE.
       DATA%NFreqTime = 10
       DATA%point_number = 2
       DATA%file_time = 'time_series'
    end select



    ALLOCATE(DATA%coor(1:DATA%point_number))

    !READ(UNIT = 10, IOSTAT = sta___tut39, FMT =*)
    !READ(UNIT = 10, IOSTAT = sta___tut64, FMT =*)(DATA%coor(ii)%x, DATA%coor(ii)%y, bla, ii=1, DATA%point_number)




    !Ouverture des fichiers pour sortie graphique
    ! IF(iproc == 0) THEN
    dir = 'POST/'

    Run%UMax = 0.0D0
    Run%Time = DATA%TimeInit

    Run%parallelism = .TRUE.
    Run%debug = .FALSE.



  END SUBROUTINE Initialise



  SUBROUTINE TimeStep(DATA, Run, MESH)

    !USE fortranf
    USE modCST

    IMPLICIT NONE

    TYPE(MeshDef), INTENT(IN) :: MESH
    TYPE(MyRun), INTENT(INOUT) :: Run
    TYPE(DATA_IN), INTENT(IN) :: DATA
    REAL (KIND = 8 ) :: dr
    !REAL (KIND = 8 ) :: Dx1
    !REAL (KIND = 8 ) :: Dx2

    !REAL (KIND = 8 ), DIMENSION(1:DATA%Nbx-1) :: dele_i
    !INTEGER :: i
    CHARACTER (LEN =* ), PARAMETER :: sub__name = "timestep"
    LOGICAL :: debug
    LOGICAL :: parallelisme

    debug= .FALSE.
    parallelisme = .FALSE.

    IF (debug) THEN
       !PRINT *, "Entering ", sub__name!, " ", iproc
    END IF

    dr = MESH%dx

    !IF (debug) THEN
    !PRINT*,"dr",dr, " CFL==", DATA%CFL
    !PRINT*,"UMax", Run%UMax
    !END IF


    Run%Dt = MIN(dr*DATA%CFL/Run%UMAX, DATA%TimeOut - Run%Time)


    !WRITE (UNIT = 6, FMT = * )'DT = ',Run%Dt
    !WRITE (UNIT = 6, FMT = * )'DT = ',Run%Time
    !WRITE (UNIT = 6, FMT = * )'DT = ',Data%Timeout - Run%Time

  END SUBROUTINE TimeStep


  FUNCTION bessj1( X ) RESULT ( bessj1__res )
    !USE fortranf

    IMPLICIT NONE
    REAL (KIND = 8) :: bessj1__res
    REAL (KIND = 8), INTENT ( IN ) :: X
    REAL (KIND = 8) :: p1
    REAL (KIND = 8) :: p2
    REAL (KIND = 8) :: p3
    REAL (KIND = 8) :: p4
    REAL (KIND = 8) :: p5
    REAL (KIND = 8) :: q1
    REAL (KIND = 8) :: q2
    REAL (KIND = 8) :: q3
    REAL (KIND = 8) :: q4
    REAL (KIND = 8) :: q5
    REAL (KIND = 8) :: r1
    REAL (KIND = 8) :: r2
    REAL (KIND = 8) :: r3
    REAL (KIND = 8) :: r4
    REAL (KIND = 8) :: r5
    REAL (KIND = 8) :: r6
    REAL (KIND = 8) :: s1
    REAL (KIND = 8) :: s2
    REAL (KIND = 8) :: S3
    REAL (KIND = 8) :: s4
    REAL (KIND = 8) :: s5
    REAL (KIND = 8) :: s6
    REAL (KIND = 8) :: Y
    REAL (KIND = 8) :: Ax
    REAL (KIND = 8) :: xx
    REAL (KIND = 8) :: Z
    !CHARACTER (LEN = *), PARAMETER :: mod__name = "bessj1"

    ! IF (imprimer_la_trace ) THEN
    !    PRINT *, "Entering ", mod__name
    ! END IF
    p1 = 1.0D0
    p2 = cte__0_d_183105e_m_2_r__8
    p3 =( - cte__0_d_3516396496e_m_4_r__8)
    p4 = cte__0_d_2457520174e_m_5_r__8
    p5 =( - cte__0_d_240337019e_m_6_r__8)
    q1 = cte__0_d_04687499995_r__8
    q2 =( - cte__0_d_2002690873e_m_3_r__8)
    q3 = cte__0_d_8449199096 - cte__5
    q4 =( - cte__0_d_88228987e_m_6_r__8)
    q5 = cte__0_d_105787412e_m_6_r__8
    r1 = cte__72362614232_d_0_r__8
    r2 =( - cte__7895059235_d_0_r__8)
    r3 = cte__242396853_d_1_r__8
    r4 =( - cte__2972611_d_439_r__8)
    r5 = cte__15704_d_48260_r__8
    r6 =( - cte__30_d_16036606_r__8)
    s1 = cte__144725228442_d_0_r__8
    s2 = cte__2300535178_d_0_r__8
    S3 = cte__18583304_d_74_r__8
    s4 = cte__99447_d_43394_r__8
    s5 = cte__376_d_9991397_r__8
    s6 = 1.0D0

    IF (ABS (X)< cte__8_d_0_r__8 ) THEN
       Y = X ** 2
       bessj1__res =( X *( r1 + Y * ( r2 + Y * ( r3 + Y * ( r4 + Y * ( r5 + Y * r6 )) ) ))) &
            / ( s1 + Y *( s2 + Y *( S3 + Y *( s4 + Y *( s5 + Y * s6 ) )) ) )

    ELSE
       Ax = ABS ( X)
       Z = cte__8_d_0_r__8 / Ax
       Y = Z ** 2
       xx = Ax - cte__2_d_356194491_r__8
       bessj1__res = SQRT ( ( cte__0_d_636619772_r__8 / Ax ) ) &
            *( COS ( xx ) * ( p1 + Y * ( p2 + Y * ( p3 + Y *( p4 + Y * p5 ) ) )) &
            - Z * SIN( xx) * ( q1 + Y * ( q2 + Y * ( q3 + Y * ( q4 + Y * q5 )) )) ) &
            * SIGN(1.0D0, X)
    END IF
  END FUNCTION bessj1



  FUNCTION bessj0(x)
    REAL(kind=8) :: bessj0
    REAL(kind=8), INTENT(in) :: x
    REAL(kind=8) :: p1, p2, p3, p4, p5, q1, q2, q3, q4, q5, r1, r2, r3, r4, r5, r6
    REAL(kind=8) :: s1, s2, s3, s4, s5, s6, y, ax, xx, z

    p1 = 1.0D0
    p2 = -0.1098628627D-2
    p3 = 0.2734510407D-4
    p4 = -0.2073370639D-5
    p5 = 0.2093887211D-6
    q1 = -0.1562499995D-1
    q2 = 0.1430488765D-3
    q3 = -0.6911147651D-5
    q4 = 0.7621095161D-6
    q5 = -0.934945152D-7
    r1 = 57568490574.0D0
    r2 = -13362590354.0D0
    r3 = 651619640.7D0
    r4 = -11214424.18D0
    r5 = 77392.33017D0
    r6 = -184.9052456D0
    s1 = 57568490411.D0
    s2 = 1029532985.0D0
    s3 = 9494680.718D0
    s4 = 59272.64853D0
    s5 = 267.8532712D0
    s6 = 1.0D0

    IF(ABS(x).LT.8.0D0) THEN
       y = x**2
       bessj0 = (r1 + y*(r2 + y*(r3 + y*(r4 + y*(r5 + y*r6)))))/(s1 + y*(s2 + y*(s3 + y*(s4 + y*(s5 + y*s6)))))
    ELSE
       ax = ABS(x)
       z = 8.0D0/ax
       y = z**2
       xx = ax  - 0.785398164D0
       bessj0 = SQRT(0.636619772D0/ax)*(COS(xx)*(p1 + y*(p2 + y*(p3 + y*(p4 + y*p5)))) & 
            - z*SIN(xx)*(q1 + y*(q2 + y*(q3 + y*(q4 + y*q5)))))
    END IF

  END FUNCTION bessj0

  ELEMENTAL FUNCTION PSI( H, Q, ZF ) RESULT (psi__res)
    !USE fortranf
    USE modCST
    IMPLICIT NONE
    REAL (KIND = 8) :: psi__res
    REAL (KIND = 8), INTENT ( IN ) :: H
    REAL (KIND = 8), INTENT ( IN ) :: Q
    REAL (KIND = 8), INTENT ( IN ) :: ZF

    !-----------------------------------------------------------------------
    !                                                  -- PSI.1
    !-----------------------------------------------------------------------

    psi__res = Q ** 2 /( cte__2_d_0_r__8 * H ** 2) + g_r__8 * ( H + ZF )

  END FUNCTION PSI


  FUNCTION Find_Last_Wet(Tab, Nbx, FLW)

    IMPLICIT NONE
    TYPE(STATE), DIMENSION(:), POINTER :: Tab
    INTEGER, INTENT(in) :: Nbx
    REAL(kind=8), INTENT(in) :: FLW
    !INTEGER, INTENT(out) :: indlast
    INTEGER :: i
    INTEGER :: Find_Last_Wet


    FIND_LAST_WET=0     
    DO i=1, Nbx

       IF (Tab(i)%h.GT.FLW) THEN
          FIND_LAST_WET = FIND_LAST_WET + 1
       ELSE
          !PRINT*,'premiere cellule seche a i=', i
          !PRINT*,'derniere cellule mouillÃ©e a i=',FIND_LAST_WET
          EXIT
       END IF

    END DO

  END FUNCTION FIND_LAST_WET

  FUNCTION Find_Last_Wet_Reverse(Tab, Nbx, FLW)

    IMPLICIT NONE
    TYPE(STATE), DIMENSION(:), POINTER :: Tab
    INTEGER, INTENT(in) :: Nbx
    REAL(kind = 8), INTENT(in) :: FLW
    !INTEGER(kind=8), INTENT(out) :: indlast
    INTEGER :: i
    INTEGER :: Find_Last_Wet_Reverse


    FIND_LAST_WET_Reverse = Nbx     
    DO i=Nbx, 1,-1

       IF (Tab(i)%h.GT.FLW) THEN
          FIND_LAST_WET_Reverse = FIND_LAST_WET_Reverse - 1
       ELSE
          !PRINT*,'premiere cellule seche a i=', i
          !PRINT*,'derniere cellule mouillÃ©e a i=',FIND_LAST_WET
          EXIT
       END IF

    END DO

  END FUNCTION FIND_LAST_WET_REVERSE


  SUBROUTINE Solex_CG_Transient(i, X_DF, Time, lambda, sigma)

    IMPLICIT NONE

    INTEGER, INTENT ( IN ) :: i
    REAL (KIND = 8), INTENT ( IN ) :: Time
    REAL (KIND = 8), INTENT ( OUT ) :: lambda
    REAL (KIND = 8), INTENT ( OUT ) :: sigma
    REAL (KIND = 8) :: l
    REAL (KIND = 8) :: alpha
    REAL (KIND = 8) :: beta
    REAL (KIND = 8) :: rho
    REAL (KIND = 8) :: theta
    REAL (KIND = 8) :: epsilon
    REAL (KIND = 8) :: S
    REAL (KIND = 8) :: A
    REAL (KIND = 8) :: X
    REAL (KIND = 8) :: atime
    REAL (KIND = 8) :: lambdaref
    REAL (KIND = 8) :: sigmaref
    REAL (KIND = 8) :: V
    REAL (KIND = 8) :: Ksi
    REAL (KIND = 8) :: arg1
    REAL (KIND = 8) :: arg2
    REAL (KIND = 8) , DIMENSION(:), POINTER :: X_DF
    CHARACTER (LEN =*), PARAMETER :: sub__name = "solex_cg_transient"

    IF (debug) THEN
       PRINT *, "Entering ", sub__name
    END IF

    S = 1.0D0 / cte__50_d_0_r__8      ! pente
    epsilon = cte__0_d_1_r__8    ! amplitude
    A = cte__1_d_5_r__8 * SQRT (( 1.0D0 + cte__0_d_9_r__8 * epsilon ) )
    l = cte__20_d_0_r__8         ! longueur d'onde
    atime = Time / SQRT (( l /( S * g_r__8) ) )
    X =( X_DF(i) + X_DF(i+1)) / cte__2_d_0_r__8
    X = X / l


    ! on calcul sugma et lambda par iterations

    lambdaref =  cte__0_d_r__8
    sigmaref =  cte__0_d_r__8

    lambda = (cte__2_d_0_r__8/A)*aTIME
    sigma =  (cte__4_d_0_r__8/a) * SQRT(ABS(-X))

    DO WHILE((ABS(lambdaref-lambda) >= 1.0D-5).OR. &
         (ABS(sigmaref-sigma) >=  1.0D-5))

       lambdaref = lambda
       sigmaref = sigma

       alpha = cte__1_d_0_r__8 +sigmaref**2-lambdaref**2
       beta = cte__2_d_0_r__8 * lambdaref
       rho = SQRT(alpha**2 + beta**2)

       IF(alpha >= cte__0_d_r__8) THEN
          theta = ATAN(beta/alpha)
       ELSE
          theta = Pi_r__8 + ATAN(beta/alpha)
       END IF

       arg1 =( cte__3_d_0_r__8 * theta ) / cte__2_d_0_r__8
       arg2 =( cte__5_d_0_r__8 * theta ) / cte__2_d_0_r__8

       v = ( cte__8_d_0_r__8 *epsilon/A)*(rho**(-cte__3_d_0_r__8/cte__2_d_0_r__8) &
            *SIN(arg1) - cte__3_d_0_r__8 / cte__4_d_0_r__8*rho**(- cte__5_d_0_r__8/cte__2_d_0_r__8) &
            *(SIN(arg2) -lambda*COS(arg2)))

       lambda = (cte__2_d_0_r__8/A)*( aTIME + v )

       alpha = cte__1_d_0_r__8 + sigmaref**2-lambda**2

       beta = cte__2_d_0_r__8*lambda

       rho = SQRT(alpha**2 + beta**2)

       IF(alpha >= cte__0_d_r__8) THEN
          theta = ATAN(beta/alpha)
       ELSE
          theta = Pi_r__8 + ATAN(beta/alpha)
       END IF

       arg1 =  cte__3_d_0_r__8*theta/cte__2_d_0_r__8
       arg2 =  cte__5_d_0_r__8*theta/cte__2_d_0_r__8

       v = (cte__8_d_0_r__8*epsilon/a)*(rho**(-cte__3_d_0_r__8/cte__2_d_0_r__8)*SIN(arg1) - ( cte__3_d_0_r__8 &
            / cte__4_d_0_r__8)*rho**(-cte__5_d_0_r__8/cte__2_d_0_r__8)*(SIN(arg2) -lambda*COS(arg2)))

       ksi = -v**2/cte__2_d_0_r__8 + epsilon*(cte__1_d_0_r__8 - cte__2_d_0_r__8*rho**(-cte__3_d_0_r__8 /cte__2_d_0_r__8) &
            *((cte__5_d_0_r__8 / cte__4_d_0_r__8)*COS(arg1) + lambda*SIN(arg1)) &
            + cte__3_d_0_r__8/cte__2_d_0_r__8*rho**(-cte__5_d_0_r__8 /cte__2_d_0_r__8)*((1.0D0-lambda**2)*COS(arg2) &
            + cte__2_d_0_r__8*lambda*SIN(arg2)))

       sigma =  (cte__4_d_0_r__8/A)*SQRT(ABS(-X + ksi))   

    END DO

  END SUBROUTINE Solex_CG_Transient



  SUBROUTINE Solex_CG_Periodic(i, X_DF, Time, lambda, sigma)

    IMPLICIT NONE
    !REAL :: f__134
    !REAL :: f__136
    !REAL :: f__138
    INTEGER, INTENT ( IN ) :: i

    REAL(KIND = 8 ), INTENT ( IN ) :: Time
    REAL(KIND = 8 ), INTENT ( OUT ) :: lambda
    REAL(KIND = 8 ), INTENT ( OUT ) :: sigma
    REAL(KIND = 8 ) :: l
    REAL(KIND = 8 ) :: S
    REAL(KIND = 8 ) :: A
    REAL(KIND = 8 ) :: X
    REAL(KIND = 8 ) :: atime
    REAL(KIND = 8 ) :: lambdaref
    REAL(KIND = 8 ) :: sigmaref
    REAL(KIND = 8 ) :: U
    REAL(KIND = 8 ) :: Ksi
    REAL(KIND = 8 ), DIMENSION( :), POINTER :: X_DF
    CHARACTER(LEN =* ), PARAMETER :: sub__name = "solex_cg_periodic"
    LOGICAL :: debug
    debug = .FALSE.

    IF(debug) THEN
       PRINT*, "Entering ", sub__name
    END IF

    S = 1.0D0 / cte__30_d_0_r__8      ! pente
    A = cte__0_d_6_r__8
    l = cte__20_d_0_r__8         ! longueur d'onde

    aTime = Time/(SQRT(l/(S*g_r__8)))

    X = (X_DF(i) + X_DF(i+1))/cte__2_d_0_r__8
    x = x/l

    sigmaref = 0.0D0
    lambdaref = 0.0D0
    lambda = cte__2_d_0_r__8 * atime
    sigma = cte__4_d_0_r__8 * SQRT (ABS ( X) )

    DO WHILE((ABS(lambdaref-lambda)>= 1.0D-10).OR.(ABS(sigmaref-sigma) >= 1.0D-10))

       sigmaref = sigma
       lambdaref = lambda
       !f__134 = bessj1( sigmaref)
       !print*,'f__134',f__134
       u = -A*bessj1(sigmaref) * SIN(lambdaref)/sigmaref
       !u = -A*f__134 * SIN(lambdaref)/sigmaref
       lambda = cte__2_d_0_r__8*(u + aTime)

       u = -A*bessj1(sigmaref) * SIN(lambda)/sigmaref
       ksi = -u**2/cte__2_d_0_r__8 + A/cte__4_d_0_r__8*bessj0(sigmaref) * COS(lambda)


       sigma = cte__4_d_0_r__8*SQRT(ABS(ksi - X))

    END DO

  END SUBROUTINE Solex_CG_Periodic




END MODULE modTOOLS








