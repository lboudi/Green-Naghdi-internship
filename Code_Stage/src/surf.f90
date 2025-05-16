PROGRAM SURF1D

  USE modCOMMON
  USE modTOOLS        ! various routines 
  USE modBOUNDARY     ! Boundary conditions
  USE modINIT         ! initial datum definition
  USE modMESH         ! meshing routines
  USE modFLUX         ! Flux computation ! 
  USE modOUTPUTS      ! outputs routines
  USE modTYPES        ! structured types
  USE modTIMESCHEME   ! time marchng evolution schemes
  USE modPOST         ! post processing (error, time series, run-up ...))


  IMPLICIT NONE
  ! Various scalars or vectors

  REAL(kind = 8) :: Sum_Init

  ! Integer scalars
  INTEGER :: INFO
  INTEGER :: Niter
  INTEGER :: N0
  INTEGER :: Nbx
  ! Logical
  !LOGICAL :: SORTIE
  INTEGER :: Ndf
  INTEGER :: ordre_espace_df
  INTEGER :: cond_lim
  REAL(KIND=8)  :: dx
  INTEGER :: i_case

  CHARACTER (LEN =* ), PARAMETER :: mod__name = "SV1D"
  REAL(KIND=8) :: Timeloop__1
  REAL(KIND=8) :: Timeloop__2
  REAL(KIND=8) :: tol


     !IF (debug) THEN
     !WRITE (UNIT = 6, FMT =*)'Main : initialisation ...'
     !WRITE (UNIT = 6, FMT =*)"Entering ", mod__name
     !END IF

     ! choice of test case
     !i_case = 1964 !effondrement d une gaussienne
     !i_case = 0 ! Dam Break
     !i_case = 1 ! soliton
     i_case = 1996 ! subcritical flow over bump
     !i_case = 2 ! Test de generation
     !i_case = 1 ! Generation d'un soliton

     CALL Initialise(i_case, DATA, Run)

     Nbx = DATA%Nbx
     ordre_espace_df = DATA%ordre_espace_df
     cond_lim = DATA%cond_lim

     !PRINT*,'cond_lim',cond_lim
     !PRINT*,'Nbx', Nbx
     !PRINT*,'ordre_espace_df', ordre_espace_df
     !READ*


     !WRITE(UNIT = 6, FMT =* )'Main : opening files ...'
     CALL OPEN_ALL(DATA)

     !WRITE(UNIT = 6, FMT =* )'Main : allocation of structured types ...'
     CALL ALLOCATE_ALL(VAR, FLU, MESH, DATA, RK)

     !WRITE(UNIT = 6, FMT =* )'Main : creating mesh ...'
     ! En sortie de CREATE_MESH, les abscisses des noeuds du maillage (pour DF) dont stockés dans MESH%X_DF
     ! Les abscisses des centres des cellules (pour VF) sont stockés dans MESH%X_VF
     ! De plus, le premier noeud de maillage est dans X_DF(1), le dernier est dans X_DF(Nbx+1)
     ! Le premier abscisse de cellule est dans X_VF(1), le dernier dans X_VF(Nbx)
     ! Enfin, on a NDf = Nbx+1
     CALL CREATE_MESH(DATA, MESH)

     dx = MESH%dx
     !WRITE(6,*)'dx == ',dx

     !WRITE(UNIT = 6, FMT =* )'Starting time : ', Run%Time

     !WRITE(UNIT = 6, FMT =* )'Main : reading initial condition ...'
     CALL INITIAL_CONDITION &
          (MESH, DATA, Var,debug)





     !WRITE(6,*) 'Main : output of initial solution ...'

     N0 = 0
     CALL OUTPUT(N0, MESH, Var, DATA, Run, 'SOL1D.')


     ! IF(.NOT. DATA%SV_ONLY) THEN

     !    IF(debug) THEN
     !       WRITE(6,*)'avant creation des matrices'
     !    END IF

     !    CALL create_matdf(sD1H, sD2H, sD1Q, sD2Q, sD1B, sD2B, sD3B, ordre_espace_df, cond_lim, dx, Ndf)

     !    ! call csr2dns(sD1H, mat)
     !    ! WRITE(6,*)'mat'
     !    ! DO ii = 1, size(mat,1)
     !    !    WRITE(6,*)(mat(ii,jj),jj=1, size(mat,2))
     !    ! END DO
     !    ! READ*
     !    ! Auxiliary tabs for computations

     !    ALLOCATE(op1%ia(SIZE(sD2q%ia)))
     !    ALLOCATE(op1%a(SIZE(sD2q%a)))
     !    ALLOCATE(op1%ja(SIZE(sD2q%ja)))

     !    ALLOCATE(op2%ia(SIZE(sD1q%ia)))
     !    ALLOCATE(op2%a(SIZE(sD1q%a)))
     !    ALLOCATE(op2%ja(SIZE(sD1q%ja)))

     !    ALLOCATE(op3%ia(SIZE(sD1q%ia)))
     !    ALLOCATE(op3%a(SIZE(sD1q%a)))
     !    ALLOCATE(op3%ja(SIZE(sD1q%ja)))

     !    ALLOCATE(op4%ia(SIZE(op1%ia)))
     !    ALLOCATE(op4%a(SIZE(op1%a)))
     !    ALLOCATE(op4%ja(SIZE(op1%ja)))

     !    ALLOCATE(op5%ia(SIZE(op1%ia)))
     !    ALLOCATE(op5%a(SIZE(op1%a)+ndf))
     !    ALLOCATE(op5%ja(SIZE(op1%ja)+ndf))

     !    ALLOCATE(op6%ia(SIZE(op1%ia)))
     !    ALLOCATE(op6%a(SIZE(op1%a)+ndf))
     !    ALLOCATE(op6%ja(SIZE(op1%ja)+ndf))

     !    ALLOCATE(opGN%ia(SIZE(op1%ia)))
     !    ALLOCATE(opGN%a(SIZE(op1%a)+ndf))
     !    ALLOCATE(opGN%ja(SIZE(op1%ja)+ndf))

     !    !ALLOCATE(umf_newop%Ap(SIZE(op1%ia)+1))
     !    !ALLOCATE(umf_newop%Ax(SIZE(op1%a)+ndf))
     !    !ALLOCATE(umf_newop%Ai(SIZE(op1%ja)+ndf))


     !    !Vecteurs de bathymétrie
     !    !IF(debug) THEN
     !    WRITE(6,*)'Creation des D1_b, D2_b '
     !    !END IF
     !    ! WRITE(6,*)'sD1b%ncol',sD1b%ncol
  !    ! WRITE(6,*)'sD2b%ncol',sD2b%ncol
  !    ! WRITE(6,*)'sD3b%ncol',sD3b%ncol
  !    ALLOCATE (D1_b(1:NDf))

  !    ALLOCATE (D2_b(1:NDf))

  !    ALLOCATE (D3_b(1:NDf))


  !WRITE (UNIT = 6, FMT =* ) 'Main : starting evolution loop ...'
  !WRITE (UNIT = 6, FMT =* ) 'Press one touch to enter evolution loop'
  !READ*


  ! -- [MAIN TIME EVOLUTION LOOP -- BEGIN] -- 
  ! ----------------------------------------

  IF(DATA%ChoixDt) THEN
     Run%Dt = DATA%DtFix
     !WRITE (UNIT = 6, FMT =* )'DT', Run%Dt
  END IF
  !PRINT*,'before evolution loop'
  !READ*

  CALL CPU_TIME(Timeloop__1)
  loop__main : DO Niter = 1, DATA%NtMax !-- [MAIN TIME EVLUTION LOOP -- BEGIN] 

     !WRITE (UNIT = 6, FMT =* )'======   Iteration =====', Niter

     CALL time_loop(Var, FLU, MESH, DATA, Run)


     ! ### INCREMENTATION DU TEMPS TOTAL ###

     Run%time = Run%time + Run%Dt

     !PRINT*,'Time = ',Run%time

     IF(Run%time >= DATA%TimeOut) THEN

        !PRINT*,'TIME OUT =', DATA%TimeOut

        CALL OUTPUT(Niter, MESH, Var, DATA, Run, 'SOL1D.')

        IF(DATA%ERROR) THEN
           !WRITE(6,*)'COMPUTATION OF ERROR'
           CALL NORM_ERROR(Var%VF, Var%Ex, Var%Ref, DATA, Run%time, sum_init, NN=Niter, ERR=ERROR)
        END IF

        CALL DEALLOCATE_ALL()
        CALL CLOSE_ALL(DATA)
        STOP ! '-All Tabs are de-allocated now, End of computation.'

     END IF


     IF(0==MOD(Niter,DATA%NFreq))  THEN
        CALL OUTPUT(Niter, MESH, Var, DATA, Run,"SOL1D.")
     END IF

     ! L2-ERROR AND L2-TIME-ERROR
     IF((DATA%ERROR.AND.(0==MOD(Niter, DATA%NFreqTime))).OR.((Niter==1).AND.DATA%ERROR)) THEN

        !WRITE(6,*)'COMPUTATION OF L2 ERROR'
        CALL NORM_ERROR(Var%VF, Var%Ex, Var%Ref, DATA, Run%time, sum_init, Niter, ERR=ERROR)

        !WRITE(UNIT = 6, FMT = *)'L2 TIME ERROR'
        !CALL TIME_ERROR(RUN, ERROR)
     END IF

  END DO loop__main ! Time evolution loop -- END --

  CALL CPU_TIME(Timeloop__2)
  !WRITE(UNIT = 6, FMT = *)'CPU Time for evolution loop:', Timeloop__2-Timeloop__1, 'seconds'

  CALL OUTPUT(Niter, MESH, Var, DATA, Run, "SOL1D.")



  CALL DEALLOCATE_ALL()



  DEALLOCATE(RK(1)%a, RK(1)%b, RK(1)%s_dt, RK(2)%a, RK(2)%b,RK(2)%s_dt, RK(3)%a, RK(3)%b,RK(3)%s_dt, RK(4)%a, RK(4)%b,RK(4)%s_dt)

  CALL CLOSE_ALL(DATA)

END PROGRAM SURF1D








