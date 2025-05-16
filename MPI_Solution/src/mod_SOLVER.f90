MODULE modSolver

  ! Riemann Solvers

  USE modTOOLS

  PRIVATE
  PUBLIC :: Solver_VFRoe_ncv
  PUBLIC :: Solver_HLLC
  !PUBLIC :: Solver_Exact
  PUBLIC :: Solver_EDP
  !PUBLIC :: VFRoe_WB
  !PUBLIC :: VFRoe_WB_topo
  PUBLIC :: Solver_Godunov
CONTAINS

  ELEMENTAL SUBROUTINE Solver_EDP ( arg_sg, arg_sd, argS)         !-- [INLINE-DEBUT]


    !USE fortranf
    USE modCST
    IMPLICIT NONE

    !LOGICAL, INTENT ( IN ) :: arg_change_dir  !-- [INLINE-DECL]
    REAL(KIND = 8) :: acm !-- [INLINE-DECL]
    REAL(KIND = 8) :: aum !-- [INLINE-DECL]
    !REAL(KIND = 8) :: avm !-- [INLINE-DECL]
    REAL(KIND = 8) :: lambda1 !-- [INLINE-DECL]
    !REAL(KIND = 8) :: lambda2 !-- [INLINE-DECL]
    REAL(KIND = 8) :: lambda3 !-- [INLINE-DECL]
    !REAL(KIND = 8) :: xxx !-- [INLINE-DECL]
    TYPE( state ), INTENT ( IN OUT ) :: arg_sg  !-- [INLINE-DECL]
    TYPE( state ), INTENT ( IN OUT ) :: arg_sd  !-- [INLINE-DECL]
    TYPE( state ), INTENT ( OUT ) :: argS  !-- [INLINE-DECL]
    !CHARACTER(LEN =*), PARAMETER :: mod__name = "solver_EDP"


    acm = cte__0_d_5_r__8 *( arg_sg%C + arg_sd%C)
    aum = cte__0_d_5_r__8 *( arg_sg%U + arg_sd%U)

    ! ##### VALEURS PROPRES NUMERIQUES #####
    lambda1 = aum - acm
    lambda3 = aum + acm


    ! ##### CALCUL DE L'ETAT A L'INTERFACE #####
    IF (arg_sg%H<= DATA%kappa .AND. arg_sd%H<= DATA%kappa ) THEN

       ! 2 MAILLES SECHES
       !PRINT*,'2 MAILLES SECHES'
       argS%C = 0.0D0
       argS%U = 0.0D0

       argS%H = DATA%kappa

    ELSE
       ! CAS CLASSIQUE
       IF(lambda1 > 0.0D0 ) THEN
          !PRINT*,'CAS CLASSIQUE decentrement a droite'
          argS%C = arg_sg%C
          argS%U = arg_sg%U

          argS%H = argS%C ** 2 / g_r__8

       ELSE IF( lambda3< 0.0D0 ) THEN
          !PRINT*,'CAS CLASSIQUE decentrement a gauche'
          argS%C = arg_sd%C
          argS%U = arg_sd%U

          argS%H = argS%C ** 2 / g_r__8

       ELSE
          argS%h = 0.5D0*(arg_sg%h + arg_sd%h) &
               - 0.25D0*(arg_sd%u - arg_sg%u)*(arg_sg%h + arg_sd%h)/(arg_sg%c + arg_sd%c)
          argS%U = 0.5D0*(arg_sg%u + arg_sd%u) &
               - (arg_sd%h - arg_sg%h)*(arg_sg%c + arg_sd%c)/(arg_sg%h + arg_sd%h)

       END IF

       IF(arg_sg%U - arg_sg%C< 0.0D0 .AND. arg_sd%U - arg_sd%C > 0.0D0) THEN
          !PRINT*,'1-sonic rarefaction'
          argS%C = (1.0D0/3.0D0)*(arg_sg%U + 2.0D0*arg_sg%c)
          argS%U = (1.0D0/3.0D0)*(arg_sg%U + 2.0D0*arg_sg%c)

          argS%H = argS%C ** 2 / g_r__8
       END IF

       IF(arg_sg%U + arg_sg%C< 0.0D0 .AND. arg_sd%U + arg_sd%C > 0.0D0) THEN
          !PRINT*,'3-sonic rarefaction'
          argS%C = (1.0D0/3.0D0)*(-arg_sd%U + 2.0D0*arg_sd%c)
          argS%U = (1.0D0/3.0D0)*(arg_sd%U - 2.0D0*arg_sd%c)

          argS%H = argS%C ** 2 / g_r__8
       END IF


    END IF ! classical case

  END SUBROUTINE Solver_EDP

!   SUBROUTINE Solver_EDP(arg_sg, arg_sd, argS, arg_change_dir)

!     ! Approximate Riemann Solver based on Exact Depth Positivy property
!     LOGICAL, INTENT (IN) :: arg_change_dir  !-- [INLINE-DECL]
!     TYPE( state ), INTENT (INOUT) :: arg_sg  !-- [INLINE-DECL]
!     TYPE( state ), INTENT (INOUT) :: arg_sd  !-- [INLINE-DECL]
!     CHARACTER(LEN =*), PARAMETER :: mod__name = "solver_EDP"   
!     TYPE( state ), INTENT ( OUT ) :: argS  !-- [INLINE-DECL]

!     ! Permutation des variables si 'change_dir=.TRUE.'
!     ! on calculera dans la direction y
!     IF(arg_change_dir) THEN         !-- [INLINE-CORPS]
!        xxx = arg_sg%U
!        arg_sg%U = arg_sg%V
!        arg_sg%V = xxx
!        xxx = arg_sd%U
!        arg_sd%U = arg_sd%V
!        arg_sd%V = xxx
!     END IF

!     IF (arg_sg%H<= DATA%kappa .AND. arg_sd%H<= DATA%kappa) THEN

!        argS%C = 0.0D0
!        argS%U = 0.0D0
!        argS%V = 0.0D0
!        argS%H = 0.0D0

!     ELSE

!        ! CAS CLASSIQUE
!        IF(lambda1 > 0.0D0 ) THEN
!           !PRINT*,'CAS CLASSIQUE decentrement a droite'
!           argS%C = arg_sg%C
!           argS%U = arg_sg%U
!           argS%V = arg_sg%V

! !           IF(argS%C <= DATA%kappa) THEN
! !              argS%H = DATA%kappa
! !              argS%U = 0.0D0
! !              argS%V = 0.0D0
! !           ELSE
!              argS%H = argS%C ** 2 / g_r__8
!           !END IF

!        ELSE IF( lambda3< 0.0D0 ) THEN
!           !PRINT*,'upwind : right'
!           argS%C = arg_sd%C
!           argS%U = arg_sd%U
!           argS%V = arg_sd%V

! !           IF(argS%C <= DATA%kappa) THEN
! !              argS%H = DATA%kappa
! !              argS%U = 0.0D0
! !              argS%V = 0.0D0
! !           ELSE
!              argS%H = argS%C ** 2 / g_r__8
!           !END IF

!        ELSE
!           ! WRIRE(6,*) 'Star region'
!           argS%h = 0.5D0*(arg_sg%h + arg_sd%h) &
!                - 0.25D0*(arg_sd%u - arg_sg%u)*(arg_sg%h + arg_sd%h)/(arg_sg%c + arg_sd%c)
!           argS%U = 0.5D0*(arg_sg%u + arg_sd%u) &
!                - (arg_sd%h - arg_sg%h)*(arg_sg%c + arg_sd%c)/(arg_sg%h + arg_sd%h)

!           IF(aum <= 0.0D0 ) THEN
!              argS%V = arg_sd%V
!           ELSE
!              argS%V = arg_sg%V
!           END IF

!           IF(argS%h <= DATA%kappa) THEN
!              argS%h = DATA%kappa
!              argS%U = 0.0D0
!              argS%V = 0.0D0
!           END IF

!        END IF
!     END IF
!     ! Left or right sonic-rarefaction

! !     IF(arg_sg%U - arg_sg%C< 0.0D0 .AND. arg_sd%U - arg_sd%C > 0.0D0) THEN
! !        !PRINT*,'1-sonic rarefaction'
! !        argS%C = (1.0D0/3.0D0)*(arg_sg%U + 2.0D0*arg_sg%c)
! !        argS%U = (1.0D0/3.0D0)*(arg_sg%U + 2.0D0*arg_sg%c)
! !        IF(argS%C <= DATA%kappa) THEN
! !           argS%H = DATA%kappa
! !           argS%U = 0.0D0
! !           argS%V = 0.0D0
! !        ELSE
! !           argS%H = argS%C ** 2 / g_r__8
! !        END IF
! !     END IF

! !     IF(arg_sg%U + arg_sg%C< 0.0D0 .AND. arg_sd%U + arg_sd%C > 0.0D0) THEN
! !        !PRINT*,'3-sonic rarefaction'
! !        argS%C = (1.0D0/3.0D0)*(-arg_sd%U + 2.0D0*arg_sd%c)
! !        argS%U = (1.0D0/3.0D0)*(arg_sd%U - 2.0D0*arg_sd%c)
! !        IF(argS%C <= DATA%kappa) THEN
! !           argS%H = DATA%kappa
! !           argS%U = 0.0D0
! !           argS%V = 0.0D0
! !        ELSE
! !           argS%H = argS%C ** 2 / g_r__8
! !        END IF
! !     END IF

!     IF(arg_change_dir ) THEN
!        xxx = arg_sg%U
!        arg_sg%U = arg_sg%V
!        arg_sg%V = xxx
!        xxx = arg_sd%U
!        arg_sd%U = arg_sd%V
!        arg_sd%V = xxx
!        xxx = argS%U
!        argS%U = argS%V
!        argS%V = xxx
!     END IF       !-- [INLINE-FIN]

!   END SUBROUTINE Solver_EDP

  SUBROUTINE Solver_HLLC(arg_sg, arg_sd, FLU)

    !LOGICAL, INTENT (IN) :: arg_change_dir  !-- [INLINE-DECL]
    TYPE( state ), INTENT (INOUT) :: arg_sg  !-- [INLINE-DECL]
    TYPE( state ), INTENT (INOUT) :: arg_sd  !-- [INLINE-DECL]
    !CHARACTER(LEN =*), PARAMETER :: mod__name = "solver_HLLC"    
    !TYPE( state ), INTENT ( OUT ) :: argS  !-- [INLINE-DECL]
    !REAL(KIND=8), INTENT(out) :: S_L
    !REAL(KIND=8), INTENT(out) :: S_R
    REAL(KIND=8) :: S_L
    REAL(KIND=8) :: S_R
    REAL(KIND=8) :: h_star
    REAL(KIND=8) :: Q_L
    REAL(KIND=8) :: Q_R
    REAL(KIND=8) :: S_star
    REAL(KIND=8), DIMENSION(2) :: UL_star  
    REAL(KIND=8), DIMENSION(2) :: UR_star
    REAL(KIND=8), DIMENSION(2) :: UL
    REAL(KIND=8), DIMENSION(2) :: UR
    REAL(KIND=8), DIMENSION(2), INTENT(OUT) :: FLU


    !print*,'kappa',DATA%kappa
!     IF (arg_sg%H <= DATA%kappa .AND. arg_sd%H <= DATA%kappa ) THEN

!        !2 MAILLES SECHES
!        !PRINT*,'2 MAILLES SECHES'
!        FLU(2) = 0.0D0
!        FLU(2) = 0.0D0
!        FLU(3) = 0.0D0

!     ELSE

       !print*,'arg_sg%h',arg_sg%h
       !print*,'arg_sd%h',arg_sd%h

       ! H_star estimation from exact depth positivity Riemann Solver
               h_star = 0.5D0*(arg_sg%h + arg_sd%h) &
                    - 0.25D0*(arg_sd%u - arg_sg%u)*(arg_sg%h + arg_sd%h)/(arg_sg%c + arg_sd%c)

       ! H_star estimation from 2-rarefaction Riemann solver

       !h_star = (1.0D0/g_r__8)*(0.5D0*(arg_sg%c + arg_sd%c) + 0.25D0*(arg_sg%u - arg_sd%u))**2

       IF(h_star > arg_sg%h) THEN
          Q_L = SQRT(0.5D0*((h_star + arg_sg%h)*h_star)/arg_sg%h**2)
       ELSE
          Q_L = 1.0D0
       END IF

       IF(h_star > arg_sd%h) THEN
          Q_R = SQRT(0.5D0*((h_star + arg_sd%h)*h_star)/arg_sd%h**2)
       ELSE
          Q_R = 1.0D0
       END IF

       IF(arg_sg%h <= DATA%kappa) THEN
          S_L = arg_sd%u - 2.0D0*arg_sd%c
          !print*,'SLkappa',S_L
          !S_L = arg_sg%u + 2.0D0*arg_sg%c
       ELSE
          S_L = arg_sg%u - arg_sg%c*Q_L
          !print*,'SLpaskappa',S_L
       END IF

       IF(arg_sd%h <= DATA%kappa) THEN
          S_R = arg_sg%u + 2.0D0*arg_sg%c
          !print*,'SRkappa',S_R
          !S_R = arg_sd%u - 2.0D0*arg_sd%c
       ELSE
          S_R = arg_sd%u + arg_sd%c*Q_R
       END IF

       S_star = (S_L*arg_sd%h*(arg_sd%u - S_R) &
            - S_R*arg_sg%h*(arg_sg%u - S_L))/(arg_sd%h*(arg_sd%u - S_R) - arg_sg%h*(arg_sg%u - S_L))

!       IF(S_star <= 1.0D-10) S_star = 0.0D0

       !print*,'S_L',S_L
       !print*,'S_R',S_R
       !print*,'S_star',S_star
       !print*,
       !read*

       UL_star(1) = arg_sg%h*((S_L - arg_sg%u)/(S_L-S_star))
       UL_star(2) = arg_sg%h*((S_L - arg_sg%u)/(S_L-S_star))*S_star


       UR_star(1) = arg_sd%h*((S_R - arg_sd%u)/(S_R-S_star))
       UR_star(2) = arg_sd%h*((S_R - arg_sd%u)/(S_R-S_star))*S_star


       UL(1) = arg_sg%h
       UL(2) = arg_sg%h*arg_sg%u


       UR(1) = arg_sd%h
       UR(2) = arg_sd%h*arg_sd%u


       IF(0.0D0<= S_L) THEN
          !FLU(1) = arg_sg%h*arg_sg%u
          !FLU(2) = arg_sg%h*arg_sg%u**2 + g_r__8*arg_sg%h**2/2.0D0
          !FLU(3) = arg_sg%h*arg_sg%u*arg_sg%v
          FLU = FLUX(arg_sg)
       ELSEIF (S_L < 0.0D0 .AND. 0.0D0 <= S_star) THEN
          !FLU(1) = arg_sg%h*arg_sg%u + S_L*(UL_star(1) - )
          !FLU(2) = arg_sg%h*arg_sg%u**2 + g_r__8*arg_sg%h**2/2.0D0
          !FLU(3) = arg_sg%h*arg_sg%u*arg_sg%v
          FLU = FLUX(arg_sg) + S_L*(UL_star - UL)
       ELSEIF (S_star <= 0.0D0 .AND. 0.0D0 < S_R) THEN
          !FLU(1) = arg_sd%h*arg_sd%u
          !FLU(2) = arg_sd%h*arg_sd%u**2 + g_r__8*arg_sd%h**2/2.0D0
          !FLU(3) = arg_sd%h*arg_sd%u*arg_sd%v
          FLU = FLUX(arg_sd) + S_R*(UR_star - UR)
       ELSEIF(S_R <= 0.0D0) THEN
          !FLU(1) = arg_sd%h*arg_sd%u
          !FLU(2) = arg_sd%h*arg_sd%u**2 + g_r__8*arg_sd%h**2/2.0D0
          !FLU(3) = arg_sd%h*arg_sd%u*arg_sd%v
          FLU = FLUX(arg_sd)
       ELSE
          WRITE(6,*) 'ERROR IN HLLC solver'
          PRINT*,'S_L',S_L
          PRINT*,'S_R',S_R
          PRINT*,'S_star',S_star
          PRINT*,'arg_sg%h',arg_sg%h
          PRINT*,'arg_sd%h',arg_sd%h
          PRINT*,'arg_sg%u',arg_sg%u
          PRINT*,'arg_sd%u',arg_sd%u
          !print*,
       !read*
          READ*
       END IF


    !END IF

  END SUBROUTINE Solver_HLLC

  ! ELEMENTAL SUBROUTINE Solver_Exact(arg_sg, arg_sd, argS, arg_change_dir)

  !   LOGICAL, INTENT ( IN ) :: arg_change_dir  !-- [INLINE-DECL]
  !   TYPE( state ), INTENT ( IN OUT ) :: arg_sg  !-- [INLINE-DECL]
  !   TYPE( state ), INTENT ( IN OUT ) :: arg_sd  !-- [INLINE-DECL]
  !   !CHARACTER(LEN =*), PARAMETER :: mod__name = "solver_exact"    
  !   TYPE( state ), INTENT ( OUT ) :: argS  !-- [INLINE-DECL]




  ! END SUBROUTINE Solver_Exact

  ELEMENTAL SUBROUTINE Solver_Godunov( arg_g, arg_d, S)     


    !USE fortranf
    USE modCST
    IMPLICIT NONE

    TYPE(state), INTENT (IN OUT) :: arg_g  !--
    TYPE(state), INTENT (IN OUT) :: arg_d  !--
    TYPE(state), INTENT (OUT) :: S 

    !CHARACTER(LEN =*), PARAMETER :: mod__name = "solver_GODUNOV"

    REAL(KIND=8)::delta_h
    REAL(KIND=8)::h_s
    REAL(KIND=8)::u_s
    REAL(KIND=8)::c_s
    REAL(KIND=8)::h_ss
    REAL(KIND=8)::ges
    !REAL(KIND=8)::g
    REAL(KIND=8)::Q_R
    REAL(KIND=8)::S_R
    REAL(KIND=8)::Q_L
    REAL(KIND=8)::S_L
    REAL(KIND=8)::SHL
    REAL(KIND=8)::STL
    REAL(KIND=8)::SHR
    REAL(KIND=8)::STR
    REAL(KIND=8)::F_R
    REAL(KIND=8)::F_L
    REAL(KIND=8)::DF_R
    REAL(KIND=8)::DF_L
    REAL(KIND=8)::c
    REAL(KIND=8)::coord
    REAL(KIND=8)::TOL
    REAL(KIND=8)::h_min
    REAL(KIND=8)::GER
    REAL(KIND=8)::GEL
    REAL(KIND=8):: crit
    !REAL(KIND=8):: SHL
    REAL(KIND=8)::        SSL
    REAL(KIND=8)::        SSR
    !REAL(KIND=8)::       SHR
    !REAL(KIND=8):: STR
    !REAL(KIND=8):: STL

    ! Critère sec/mouillé
    crit = (arg_d%u - arg_g%u) - 2.0D0*(arg_g%c + arg_d%c)

    IF( arg_g%h == 0.0D0 .OR. arg_d%h == 0.0D0 .OR. crit >= 0.0D0) THEN

       !WRITE(6,*)'WARNING DRY BED'
       !STOP

       IF(arg_g%h == 0.0D0)  THEN
          ! Left state is initially dry

          SHR = arg_d%u + arg_d%c
          IF(0.0D0 > SHR) THEN
             ! the interface is at the right of the rarefaction wave

             s%h = arg_d%h
             s%u = arg_d%u
          ELSE
             STR = arg_d%u - 2.0D0*arg_d%c

             IF(0.0D0 > STR) THEN
                ! cells interface is IN the rarefaction wave
                s%u = (arg_d%u - 2.0D0*arg_d%c)/3.0D0
                s%c = (-arg_d%u +2.0D0*arg_d%c)/3.0D0
                s%h = s%c**2/g_r__8

             ELSE
                !cells interface is in the dry region
                s%h = arg_g%h
                s%u = arg_g%u
             END IF

          END IF


       ELSEIF (arg_d%h == 0.0D0) THEN
          ! Right state is initially dry

          SHL = arg_g%u - arg_g%c
          IF(0.0D0 < SHL) THEN
             ! the interface is at the left of the rarefaction wave

             s%h = arg_g%h
             s%u = arg_g%u
          ELSE
             STL = arg_g%u + 2.0D0*arg_g%c

             IF(0.0D0 < STL) THEN
                ! cells interface is IN the rarefaction wave
                s%u = (arg_g%u + 2.0D0*arg_g%c)/3.0D0
                s%c = (arg_g%u +2.0D0*arg_g%c)/3.0D0
                s%h = s%c**2/g_r__8

             ELSE
                !cells interface is in the dry region
                s%h = arg_d%h
                s%u = arg_d%u
             END IF

          END IF

       ELSE
          ! a dry state appears

          ! Waves speeds
          SHL = arg_g%u - arg_g%c
          SSL = arg_g%u + 2.0D0*arg_g%c
          SSR = arg_d%u - 2.0D0*arg_d%c
          SHR = arg_d%u + arg_d%c

          IF (0.0D0 <= SHL) THEN

             ! cells interface is at the left of the left rarefaction
             s%h = arg_g%h
             s%u = arg_g%u
          END IF

          IF(0.0D0>SHL .AND. 0.0D0<=SSL ) THEN
             ! interface is IN the left rarefaction
             s%u = (arg_g%u + 2.0D0*arg_g%c)/3.0D0
             s%c = (arg_g%u +2.0D0*arg_g%c)/3.0D0
             s%h = s%c**2/g_r__8
          END IF

          IF(0.0D0>SSL .AND. 0.0D0<SSR) THEN
             !interface is in the middle dry bed region
             s%h = 0.0D0
             s%u = 0.0D0
          END IF

          IF(0.0D0>SSR .AND. 0.0D0<=SHR) THEN
             ! interface is IN the right rarefaction
             s%u = (arg_d%u - 2.0D0*arg_d%c)/3.0D0
             s%c = (-arg_d%u +2.0D0*arg_d%c)/3.0D0
             s%h = s%c**2/g_r__8
          END IF

          IF(0.0D0 > SHR) THEN
             ! interface is at the right of the right rarefaction wave
             s%h = arg_d%h
             s%u = arg_d%u
          END IF

       END IF

    ELSE

       ! NO DRY BED

       delta_h = 1.0D0
       TOL = 1.0D-10

       ! Valeur d'initialisation

       h_s = (1.0D0/g_r__8)*(0.5D0*(arg_g%c+arg_d%c)-0.25D0*(arg_d%u-arg_g%u))**2
       h_min = MIN(arg_g%h,arg_d%h)

       IF(h_s <= h_min)THEN
          h_s = h_s
       ELSE
          GEL = SQRT(0.5D0*g_r__8*(arg_g%h + h_s)/h_s*arg_g%h)
          GER = SQRT(0.5D0*g_r__8*(arg_d%h + h_s)/h_s*arg_d%h) 
          h_s = (GEL*arg_g%h + GER*arg_d%h - (arg_d%h - arg_g%h))/(GEL+GER)
       ENDIF


       !  CALCUL DE H ET U DANS LA REGION ETOILE

       DO WHILE(delta_h > TOL)

          IF (h_s <= arg_g%h) THEN
             c = SQRT(g_r__8*h_s)
             F_L = 2.0D0*(c - arg_g%c)
             DF_L = g_r__8 / c
          ELSE
             ges = SQRT( 0.5D0*g_r__8*(h_s+arg_g%h)/(h_s*arg_g%h))     
             F_L = (h_s - arg_g%h)*ges
             DF_L = ges - 0.25D0*g_r__8*(h_s - arg_g%h)/(ges*h_s**2)
          ENDIF

          IF (h_s <= arg_d%h)THEN
             c = SQRT(g_r__8*h_s)
             F_R = 2.0D0*(c - arg_d%c)
             DF_R = g_r__8 / c
          ELSE
             ges = SQRT( 0.5D0*g_r__8*(h_s+arg_d%h)/(h_s*arg_d%h))     
             F_R = (h_s - arg_d%h)*ges
             DF_R = ges - 0.25D0*g_r__8*(h_s - arg_d%h)/(ges*h_s**2)    

          ENDIF


          h_ss = h_s - (F_L + F_R + arg_d%u - arg_g%u)/(DF_R+DF_L)


          !IF(h_s.LE.0.0D0)THEN
          !   h_s = TOL
          !ENDIF

          delta_h = 2.0D0*ABS(h_ss - h_s)/(h_ss + h_s)
          h_s = h_ss

       END DO

       u_s = (arg_g%u + arg_d%u )/2.0D0 + (F_R - F_L)/2.0D0
       c_s = SQRT(g_r__8*h_s)

       ! CALCUL DE LA SOLUTION A L INTERFACE

       IF(u_s > 0.0D0) THEN
          IF(h_s > arg_g%h) THEN
             Q_L = SQRT((h_s+arg_g%h)*h_s/(2*arg_g%h**2))
             S_L = arg_g%u-arg_g%c*Q_L

             IF(S_L > 0.0D0) THEN
                s%h=arg_g%h
                s%u=arg_g%u
             ELSE
                s%h=h_s
                s%u=u_s
             ENDIF
          ELSE
             SHL = arg_g%u - arg_g%c

             IF(SHL > 0.0D0) THEN
                s%h = arg_g%h
                s%u = arg_g%u
             ELSE
                STL = u_s - c_s

                IF(STL < 0.0D0) THEN
                   s%h=h_s
                   s%u=u_s
                ELSE
                   s%u = (arg_g%u + 2.0D0*arg_g%c + 2.0D0*coord)/3.0D0
                   s%c = (arg_g%u + 2.0D0*arg_g%c - coord)/3.0D0  
                   s%h = s%c**2 /g_r__8
                ENDIF

             ENDIF
          ENDIF


       ELSE

          IF (h_s > arg_d%h)THEN
             Q_R = SQRT((h_s+arg_d%h)*h_s/(2.0D0*arg_d%h**2))
             S_R = arg_d%u+arg_d%c*Q_R

             IF(S_R < 0.0D0) THEN
                s%h = arg_d%h
                s%u = arg_d%u
             ELSE
                s%h = h_s
                s%u = u_s
             ENDIF

          ELSE
             STR = u_s + c_s
             IF(STR > 0.0D0) THEN
                s%h = h_s
                s%u = u_s
             ELSE
                SHR = arg_d%u + arg_d%c 
                IF(SHR < 0.0D0)THEN
                   s%h = arg_d%h
                   s%u = arg_d%u
                ELSE
                   s%u = (arg_d%u - 2.0D0*arg_d%c + 2.0D0*coord)/3.0D0
                   s%c = (-arg_d%u + 2.0D0*arg_d%c + coord)/3.0D0  
                   s%h = s%c**2 /g_r__8
                ENDIF
             ENDIF
          ENDIF

       ENDIF

    END IF

  END SUBROUTINE Solver_Godunov



  ELEMENTAL SUBROUTINE Solver_VFRoe_ncv ( arg_sg, arg_sd, argS)         !-- [INLINE-DEBUT]
 
    !USE fortranf
    USE modCST
   IMPLICIT NONE
    !LOGICAL, INTENT ( IN ) :: arg_change_dir  !-- [INLINE-DECL]
    REAL(KIND = 8) :: acm !-- [INLINE-DECL]
    REAL(KIND = 8) :: aum !-- [INLINE-DECL]
    !REAL(KIND = 8) :: avm !-- [INLINE-DECL]
    REAL(KIND = 8) :: lambda1 !-- [INLINE-DECL]
    !REAL!(KIND = 8) :: lambda2 !-- [INLINE-DECL]
    REAL(KIND = 8) :: lambda3 !-- [INLINE-DECL]
    !REAL(KIND = r__8) :: xxx !-- [INLINE-DECL]
    TYPE( state ), INTENT ( IN OUT ) :: arg_sg  !-- [INLINE-DECL]
    TYPE( state ), INTENT ( IN OUT ) :: arg_sd  !-- [INLINE-DECL]
    TYPE( state ), INTENT ( OUT ) :: argS  !-- [INLINE-DECL]
    !CHARACTER(LEN =*), PARAMETER :: mod__name = "solver_vfroe_ncv"


    acm = cte__0_d_5_r__8 *( arg_sg%C + arg_sd%C)
    aum = cte__0_d_5_r__8 *( arg_sg%U + arg_sd%U)

    ! ##### VALEURS PROPRES NUMERIQUES #####
    lambda1 = aum - acm
    lambda3 = aum + acm


    ! ##### CALCUL DE L'ETAT A L'INTERFACE #####
    IF (arg_sg%H<= DATA%kappa .AND. arg_sd%H<= DATA%kappa ) THEN

       ! 2 MAILLES SECHES
       !PRINT*,'2 MAILLES SECHES'
       argS%C = 0.0D0
       argS%U = 0.0D0

    ELSE
       ! CAS CLASSIQUE
       IF(lambda1 > 0.0D0 ) THEN
          !PRINT*,'CAS CLASSIQUE decentrement a droite'
          argS%C = arg_sg%C
          argS%U = arg_sg%U

       ELSE IF( lambda3< 0.0D0 ) THEN
          !PRINT*,'CAS CLASSIQUE decentrement a gauche'
          argS%C = arg_sd%C
          argS%U = arg_sd%U

       ELSE
          !PRINT*,'CAS CLASSIQUE pas de decentrement'
          argS%C = acm - cte__0_d_25_r__8 * ( arg_sd%U - arg_sg%U )
          argS%U = aum -( g_r__8 /( cte__2_d_0_r__8 * acm )) * ( arg_sd%H - arg_sg%H )

       END IF

       ! CORRECTION ENTROPIQUE DU SCHEMA(2c, u, v)

       IF(arg_sg%U - arg_sg%C< 0.0D0 .AND. arg_sd%U - arg_sd%C > 0.0D0) THEN
          !PRINT*,'Entropic correction: 1-detente sonique'
          argS%C = cte__0_d_5_r__8 *( arg_sg%C + arg_sd%C)
          argS%U = cte__0_d_5_r__8 *( arg_sg%U + arg_sd%U)
       END IF

       IF(arg_sg%U + arg_sg%C< 0.0D0 .AND. arg_sd%U + arg_sd%C > 0.0D0) THEN
          !PRINT*,'Entropic correction : 3-detente sonique'
          argS%C = cte__0_d_5_r__8 *( arg_sg%C + arg_sd%C)
          argS%U = cte__0_d_5_r__8 *( arg_sg%U + arg_sd%U)
       END IF

    END IF ! classical case


    IF(argS%C <= kappa) THEN
       argS%H = DATA%kappa
       argS%U = 0.0D0
    ELSE
       argS%H = argS%C ** 2 / g_r__8
    END IF

  END SUBROUTINE Solver_VFRoe_ncv





  FUNCTION FLUX(S)
    !USE fortranf
    !USE les__ctes__z2f
    IMPLICIT NONE
    REAL (KIND=8), DIMENSION(2) :: FLUX
    TYPE( state ), INTENT (IN) :: S
    !CHARACTER (LEN =* ), PARAMETER :: mod__name = "FLUX"

    FLUX(1) = s%h*s%u
    FLUX(2) = s%h*s%u**2 + g_r__8*s%h**2/2.0D0


  END FUNCTION FLUX

END MODULE ModSolver
