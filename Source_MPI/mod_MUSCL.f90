MODULE modMUSCL

  !USE fortranf
  USE modCST
  USE modTOOLS

  PRIVATE
  PUBLIC :: MUSCL_o2a
  PUBLIC :: MUSCL_o2b
  PUBLIC :: MUSCL_o2c

  PUBLIC :: MUSCL_V4

  PUBLIC :: MUSCL_o4a
  PUBLIC :: MUSCL_o4b
  PUBLIC :: MUSCL_o4c
  PUBLIC :: MUSCL_o5

  PUBLIC :: MUSCL_3EQ

  PUBLIC :: ENO_o2
  PUBLIC :: ENOmod_o2

  PUBLIC :: minmod
  PUBLIC :: f_minmod
  PUBLIC :: van_albada
  PUBLIC :: van_leer
  PUBLIC :: Arora_Roe
  PUBLIC :: superbee
  PUBLIC :: Arora_Roe_nonlin
  PUBLIC :: Compute_Delta_star

CONTAINS
  

  ELEMENTAL SUBROUTINE MUSCL_o5(in_m2, in_m1, in, in_p1, in_p2, in_p3, out_L, out_R)

    IMPLICIT NONE

    ! Compute the four (third) order reconstruction TVD-MUSCL interface values q(i+1/2)- 
    ! and q(i+1/2)+ from q(i-1), q(i), q(i+1) and q(i+2) using Yamamoto - Computers Fluids (1993)

    ! WARNING = INTERFACE APPROACH
    ! out_L and out_R stands for the face states at left and right sides of the considred interface

    TYPE(state), INTENT(in) :: in_m2
    TYPE(state), INTENT(in) :: in_m1
    TYPE(STATE), INTENT(in) :: in
    TYPE(STATE), INTENT(in) :: in_p1
    TYPE(STATE), INTENT(in) :: in_p2
    TYPE(STATE), INTENT(in) :: in_p3

    TYPE(STATE), INTENT(out) :: out_L
    TYPE(STATE), INTENT(out) :: out_R

    !REAL(kind = 8) :: Prod_left
    !REAL(kind = 8) :: Prod_right
    !REAL!(kind = 8) :: QL
    !REAL(kind = 8) ::  QM
    !REAL(kind = 8) :: QR

    REAL(kind = 8) :: ri_L
    REAL(kind = 8) :: ri_R
    REAL(kind = 8) :: ri_m1_L
    REAL(kind = 8) :: ri_m1_R
    REAL(kind = 8) :: ri_p1_L
    REAL(kind = 8) :: ri_p1_R
    REAL(kind = 8) :: ri_p2_L
    REAL(kind = 8) :: ri_p2_R
    REAL(kind = 8) :: phi_star_L
    REAL(kind = 8) :: phi_star_R
    REAL(kind = 8) :: phi_lim_L
    REAL(kind = 8) :: phi_lim_R

    ! h reconstruction

    ri_L = (in_p1%h - in%h)/(in%h - in_m1%h)
    ri_R = (in%h - in_m1%h)/(in_p1%h - in%h)

    ri_m1_L = (in%h - in_m1%h)   / (in_m1%h - in_m2%h)
    ri_m1_R = (in_m1%h - in_m2%h)/ (in%h - in_m1%h)

    ri_p1_L = (in_p2%h - in_p1%h)/(in_p1%h - in%h)
    ri_p1_R = (in_p1%h - in%h)/(in_p2%h - in_p1%h)

    ri_p2_L = (in_p3%h - in_p2%h)/(in_p2%h - in_p1%h)
    ri_p2_R = (in_p2%h - in_p1%h)/(in_p3%h - in_p2%h)

    ! Fifth order limiter (Kim and Kim 2005))
    Phi_star_L = (-2.0D0/ri_m1_L + 11.0D0 + 24.0D0*ri_L - 3.0D0*ri_L*ri_p1_L)/30.0D0

    Phi_star_R = (-2.0D0/ri_p2_R + 11.0D0 + 24.0D0*ri_p1_R - 3.0D0*ri_p1_R*ri_R)/30.0D0

    ! Monotonicity preserving

    Phi_lim_L = MAX(0.0D0, MIN(2.0D0, 2.0D0*ri_L, Phi_star_L))

    Phi_lim_R = MAX(0.0D0, MIN(2.0D0, 2.0D0*ri_R, Phi_star_R))

    out_L%h = in%h + 0.5D0*Phi_lim_L*(in%h - in_m1%h)

    out_R%h = in_p1%h - 0.5D0*Phi_lim_R*(in_p2%h - in_p1%h)

    ! qx reconstruction

    ri_L = (in_p1%qx - in%qx)/(in%qx - in_m1%qx)
    ri_R = (in%qx - in_m1%qx)/(in_p1%qx - in%qx)

    ri_m1_L = (in%qx - in_m1%qx)   / (in_m1%qx - in_m2%qx)
    ri_m1_R = (in_m1%qx - in_m2%qx)/ (in%qx - in_m1%qx)

    ri_p1_L = (in_p2%qx - in_p1%qx)/(in_p1%qx - in%qx)
    ri_p1_R = (in_p1%qx - in%qx)/(in_p2%qx - in_p1%qx)

    ri_p2_L = (in_p3%qx - in_p2%qx)/(in_p2%qx - in_p1%qx)
    ri_p2_R = (in_p2%qx - in_p1%qx)/(in_p3%qx - in_p2%qx)

    ! Fifth order limiter (Kim and Kim 2005))
    Phi_star_L = (-2.0D0/ri_m1_L + 11.0D0 + 24.0D0*ri_L - 3.0D0*ri_L*ri_p1_L)/30.0D0

    Phi_star_R = (-2.0D0/ri_p2_R + 11.0D0 + 24.0D0*ri_p1_R - 3.0D0*ri_p1_R*ri_R)/30.0D0

    ! Monotonicity preserving

    Phi_lim_L = MAX(0.0D0, MIN(2.0D0, 2.0D0*ri_L, Phi_star_L))

    Phi_lim_R = MAX(0.0D0, MIN(2.0D0, 2.0D0*ri_R, Phi_star_R))

    out_L%qx = in%qx + 0.5D0*Phi_lim_L*(in%qx - in_m1%qx)

    out_R%qx = in_p1%qx - 0.5D0*Phi_lim_R*(in_p2%qx - in_p1%qx)


    ! h+z reconstruction

    ri_L = (in_p1%xi - in%xi)/(in%xi - in_m1%xi)
    ri_R = (in%xi - in_m1%xi)/(in_p1%xi - in%xi)

    ri_m1_L = (in%xi - in_m1%xi)   / (in_m1%xi - in_m2%xi)
    ri_m1_R = (in_m1%xi - in_m2%xi)/ (in%xi - in_m1%xi)

    ri_p1_L = (in_p2%xi - in_p1%xi)/(in_p1%xi - in%xi)
    ri_p1_R = (in_p1%xi - in%xi)/(in_p2%xi - in_p1%xi)

    ri_p2_L = (in_p3%xi - in_p2%xi)/(in_p2%xi - in_p1%xi)
    ri_p2_R = (in_p2%xi - in_p1%xi)/(in_p3%xi - in_p2%xi)

    ! Fifth order limiter (Kim and Kim 2005))
    Phi_star_L = (-2.0D0/ri_m1_L + 11.0D0 + 24.0D0*ri_L - 3.0D0*ri_L*ri_p1_L)/30.0D0

    Phi_star_R = (-2.0D0/ri_p2_R + 11.0D0 + 24.0D0*ri_p1_R - 3.0D0*ri_p1_R*ri_R)/30.0D0

    ! Monotonicity preserving

    Phi_lim_L = MAX(0.0D0, MIN(2.0D0, 2.0D0*ri_L, Phi_star_L))

    Phi_lim_R = MAX(0.0D0, MIN(2.0D0, 2.0D0*ri_R, Phi_star_R))

    out_L%xi = in%xi + 0.5D0*Phi_lim_L*(in%xi - in_m1%xi)

    out_R%xi = in_p1%xi - 0.5D0*Phi_lim_R*(in_p2%xi - in_p1%xi)

  END SUBROUTINE MUSCL_o5

  ELEMENTAL SUBROUTINE MUSCL_o2a(in_L, in_M, in_R, &
       GradH, GradHZ, GradQx, Dx)

    TYPE(STATE), INTENT(in) :: in_L
    TYPE(STATE), INTENT(in) :: in_M
    TYPE(STATE), INTENT(in) :: in_R

    REAL(kind = 8), INTENT(in) :: Dx

    REAL(kind = 8), INTENT(out) :: GradH
    REAL(kind = 8), INTENT(out) :: GradQx
    REAL(kind = 8), INTENT(out) :: GradHZ

    REAL(kind = 8) :: ZG
    REAL(kind = 8) :: ZM
    REAL(kind = 8) :: ZD

    REAL(kind = 8) :: PH
    REAL(kind = 8) :: PHZ
    REAL(kind = 8) :: PQx

    REAL(kind = 8), PARAMETER :: KAPPA = 0.0D0

    ! Initialisation de Z

    ZG = in_L%h + in_L%topo
    ZM = in_M%h + in_M%topo
    ZD = in_R%h + in_R%topo

    ! Produit pour la Monotonie

    PH = (in_M%h - in_L%h) * (in_R%h - in_M%h)
    PHZ  = (ZM - ZG) * (ZD - ZM)
    PQx = (in_M%qx - in_L%qx) * (in_R%qx - in_M%qx)

    ! Calcul du gradient du debit en hu

    ! Test sur la monotonie
    IF(PQx <= cte__0_d_r__8) THEN

       !QM est un extremum local
       GradQx = cte__0_d_r__8
    ELSE
       !Q est localement monotone
       GradQx = MIN(ABS(in_M%qx-in_L%qx), ABS(in_R%qx-in_M%qx))/Dx
       GradQx = SIGN(GradQx, in_R%qx-in_M%qx)
       !GradQx = ((in_R%qx - in_M%qx)*(in_M%qx-in_L%qx) + abs(in_M%qx-in_L%qx)*ABS(in_R%qx-in_M%qx))/(in_R%qx-in_L%qx)/Dx
    END IF
 
    ! Calcul du gradient pour la reconstruction de h
    IF ( (PHZ <= 0.0D0).OR.(in_L%h <= 2.D0*KAPPA) &
         .OR.(in_M%h <= 2.D0*KAPPA) &
         .OR.(in_R%h <= 2.D0*KAPPA) &
         .OR.(PH <= 0.0D0)) THEN

       !QM est un extremum local
       GradH = cte__0_d_r__8

    ELSE
       !Q est localement monotone

       GradH = MIN(ABS(in_M%h-in_L%h), ABS(in_R%h-in_M%h))/Dx
       GradH = SIGN(GradH, in_R%h-in_M%h)
       !GradH = (in_R%h - in_M%h)*(in_M%h-in_L%h) + abs(in_M%h-in_L%h)*ABS(in_R%h-in_M%h)/(in_R%h-in_L%h)/Dx

       !print*,'GradH',GradH
    END IF
    ! Calcul du gradient pour la reconstruction de h + z

    ! Test sur la monotonie

    IF ( (PHZ <= 0.0D0).OR.(in_L%h <= 2.D0*KAPPA) &
         .OR.(in_M%h <= 2.D0*KAPPA) &
         .OR.(in_R%h <= 2.D0*KAPPA)) THEN
       GradHZ = cte__0_d_r__8
       GradH = cte__0_d_r__8

    ELSE
       GradHZ = MIN(ABS(ZM-ZG), ABS(ZD-ZM))/ Dx
       GradHZ = SIGN(GradHZ, ZM-ZG)
       !GradHZ = (ZD - ZM)*(ZM-ZG) + abs(ZM-ZG)*ABS(ZD-ZM)/(ZD-ZG)/Dx
    END IF
    !PRINT*,'GradQx, GradQy et GradH', GradQx, GradQy, GradH

  END SUBROUTINE MUSCL_o2a


  ELEMENTAL SUBROUTINE MUSCL_o2b(Hi_m1, Hi, Hi_p1, Hi_p2, H_left, H_right, limiter, CFL)

    IMPLICIT NONE

    ! Non-conservative formulation, with variable limiters options
    ! Interface/vertex approach

    ! Compute the second order reconstruction TVD-MUSCL interface values q(i+1/2)- 
    ! and q(i+1/2)+ from q(i-1), q(i), q(i+1) and q(i+2) using classical TVD-limiters

    REAL(kind = 8), INTENT(in) :: Hi_m1
    REAL(kind = 8), INTENT(in) :: Hi
    REAL(kind = 8), INTENT(in) :: Hi_p1
    REAL(kind = 8), INTENT(in) :: Hi_p2
    !REAL(kind = 8), INTENT(in) :: Dxi
    !REAL(kind = 8), INTENT(in) :: Dxip
    !REAL(kind = 8), INTENT(in) :: Dx
    REAL(kind = 8), INTENT(in) :: CFL

    REAL(kind = 8), INTENT(out) :: H_left
    REAL(kind = 8), INTENT(out) :: H_right

    REAL(kind = 8) :: Prod_left
    REAL(kind = 8) :: Prod_right
    REAL(kind = 8) :: QL
    REAL(kind = 8) :: QM
    REAL(kind = 8) :: QR
    REAL(kind = 8) :: rr_1
    REAL(kind = 8) :: rr_2
    REAL(kind = 8) :: phi_1
    REAL(kind = 8) :: phi_2
    INTEGER, INTENT(in) :: limiter

    ! Produit pour la Monotonie

    Prod_left = (Hi - Hi_m1) * (Hi_p1 - Hi)
    Prod_right = (Hi_p1 - Hi) * (Hi_p2 - Hi_p1)

    QL = Hi - Hi_m1
    QM = Hi_p1 - Hi
    QR = Hi_p2 - Hi_p1

    ! #########

    rr_1 = QM/QL
    IF(Prod_left <= cte__0_d_r__8) THEN

       ! Hi is a local extremum
       phi_1 = cte__0_d_r__8
     
    ELSE

       ! slope limiter
       SELECT CASE(limiter)
       CASE(1)
          phi_1 = minmod(rr_1)
       CASE(2)
          phi_1 = van_albada(rr_1)
       CASE(3)
          phi_1 = van_leer(rr_1)
       CASE(4)
          phi_1 = Arora_Roe(rr_1, CFL)
       CASE(5)
          phi_1 = superbee(rr_1)
       CASE(6)
          phi_1 = Arora_Roe_nonlin(rr_1, CFL)
       CASE default
          !PRINT*,'BAD LIMITER CHOICE'

       END SELECT

    END IF

    H_left = Hi + cte__0_d_5_r__8*phi_1 * QL

 
    ! ##########
    rr_2 = QR/QM

    IF(Prod_right <= cte__0_d_r__8) THEN

       ! Hi+1 is a local extremum

       phi_2 = cte__0_d_r__8

    ELSE

       SELECT CASE(limiter)
       CASE(1)

          phi_2 = minmod(rr_2)
       CASE(2)
          phi_2 = van_albada(rr_2)
       CASE(3)
          phi_2 = van_leer(rr_2)
       CASE(4)
          phi_2 = Arora_Roe(rr_2, cte__0_d_8_r__8)
       CASE(5)
          phi_2 = superbee(rr_2)
       CASE(6)
          phi_2 = Arora_Roe_nonlin(rr_2, cte__0_d_8_r__8)
       CASE default
          !PRINT*,'BAD LIMITER CHOICE'
       END SELECT


    END IF

    H_right = Hi_p1 - cte__0_d_5_r__8*phi_2 * QM


  END SUBROUTINE MUSCL_O2b


  SUBROUTINE MUSCL_3EQ(in_L, in_M, in_R, Dx,&
       out_L, out_R)

    IMPLICIT NONE

    TYPE(STATE), INTENT(in) :: in_L
    TYPE(STATE), INTENT(in) :: in_M
    TYPE(STATE), INTENT(in) :: in_R

    REAL(kind = 8), INTENT(in) :: Dx

    TYPE(STATE), INTENT(out) :: out_L
    TYPE(STATE), INTENT(out) :: out_R

    REAL(kind = 8) :: GradH
    REAL(kind = 8) :: GradQx
    !REAL(kind = 8), :: Grad_topo

    REAL(kind = 8) :: ZG
    REAL(kind = 8) :: ZM
    REAL(kind = 8) :: ZD
    REAL(kind = 8) :: HMG
    REAL(kind = 8) :: HMD
    REAL(kind = 8) :: QMG
    REAL(kind = 8) :: QMD
    REAL(kind = 8) :: PsiG
    REAL(kind = 8) :: PsiM
    !REAL(kind = 8) :: ZPsiD

    REAL(kind = 8) :: PsiMG
    !REAL(kind = 8) :: PPsiMD
    !REAL!(kind = 8) :: PsiL
    !REAL(kind = 8) :: PsiR
    !REAL(kind = 8) :: PSI
    REAL(kind = 8) :: PsiD
    REAL(kind = 8) :: PsiMD
    REAL(kind = 8) :: PH
    REAL(kind = 8) :: PHZ
    REAL(kind = 8) :: PQx
    !REAL(kind = 8) :: Ptopo

    ! Initialisation de Z

    ZG = in_L%xi
    ZM = in_M%xi
    ZD = in_R%xi

    ! Produit pour la Monotonie

    PH = (in_M%h - in_L%h) * (in_R%h - in_M%h)
    PHZ  = (ZM - ZG) * (ZD - ZM)
    PQx = (in_M%qx - in_L%qx) * (in_R%qx - in_M%qx)

    !Ptopo = (ZFM - ZFG) * (ZFD - ZFM)

    ! Calcul du gradient du debit en hu

    ! Test sur la monotonie
    IF(PQx <= 0.0D0) THEN

       !QM est un extremum local
       GradQx = 0.0D0
    ELSE
       !Q est localement monotone
       !GradQx = MIN(ABS(in_M%qx-in_L%qx), ABS(in_R%qx-in_M%qx))/Dx
       !GradQx = SIGN(GradQx, in_R%qx-in_M%qx)
       GradQx = ((in_R%qx - in_M%qx)*(in_M%qx-in_L%qx) + ABS(in_M%qx-in_L%qx)*ABS(in_R%qx-in_M%qx))/(in_R%qx-in_L%qx)/Dx
    END IF

    ! Calcul du gradient pour la reconstruction de h

    IF((in_L%topo == in_M%topo).AND.(in_R%topo == in_M%topo)) THEN 

       ! On reconstruit h clasiquement
       IF (PH <= 0.0D0) THEN
          !QM est un extremum local
          GradH = cte__0_d_r__8
       ELSE
          GRADH = MIN( ABS(in_M%h-in_L%h), ABS(in_r%h-in_M%h))/DX
          GRADH = SIGN( GRADH, in_M%h-in_L%h)
       END IF

    ELSE
       ! on reconstruit h avec xi
       IF (PHZ <= 0.0D0) THEN
          !QM est un extremum local
          GradH = cte__0_d_r__8
       ELSE
          GRADH = MIN( MIN( ABS(ZM-ZG), ABS(ZD-ZM) ), in_M%h )/DX
          GRADH = SIGN( GRADH, ZM-ZG )
          !GradH = (ZD - ZM)*(ZM-ZG) + abs(ZM-ZG)*ABS(ZD-ZM)/(ZD-ZG)/Dx
       END IF
    END IF

    ! Valeurs aux interfaces

    HMG = in_M%h - (DX/2.D0) * GRADH
    HMD = in_M%h + (DX/2.D0) * GRADH

    QMG = in_M%qx - (DX/2.D0) * GRADQx
    QMD = in_M%qx + (DX/2.D0) * GRADQx

    !           CALCUL DES PSI
    PSIG  = fPSI( in_L%h, in_L%qx, in_L%topo)
    PSIM  = fPSI( in_M%h, in_M%qx, in_M%topo)
    PSID  = fPSI( in_R%h, in_R%qx, in_R%topo)

    PSIMG = fPSI(HMG, QMG, in_L%topo)
    PSIMD = fPSI(HMD, QMD, in_R%topo)

!     PSIL  = (PSIG+PSIM)/2.D0
!     PSIR  = (PSIM+PSID)/2.D0

!     PSIL  = ABS(PSIG-PSIM)/2.D0
!     PSIR  = ABS(PSIM-PSID)/2.D0

!     IF(.NOT.(ABS(PSIM-PSIMG) <= PSIL)) THEN
!        GRADH = 0.0D0
!     END IF
!     IF(.NOT.(ABS(PSIM-PSIMD) <= PSIR)) THEN
!        GRADH = 0.0D0
!     END IF

!               !TEST TVD+MONOTONIE SUR PSI
!     IF (PSIG < PSIM) THEN
!        IF ((PSIMG > PSIL).AND.(PSIMG < PSIM)) THEN
!           !                ON NE LIMITE PAS
!        ELSE
!           GRADH = 0.D0
!        END IF
!     ELSEIF (PSIG > PSIM) THEN
!        IF ((PSIMG > PSIM).AND.(PSIMG < PSIL)) THEN
!           !                 ON NE LIMITE PAS
!        ELSE
!           GRADH = 0.D0
!        END IF
!     ELSE
!        GRADH = 0.D0
!     END IF

!     IF (PSID < PSIM) THEN
!        IF ((PSIMD > PSIR).AND.(PSIMD < PSIM)) THEN
!           !                 ON NE LIMITE PAS
!        ELSE
!           GRADH = 0.D0
!        END IF
!     ELSEIF (PSID > PSIM) THEN
!        IF ((PSIMD > PSIM).AND.(PSIMD < PSIR)) THEN
!           !                 ON NE LIMITE PAS
!        ELSE
!           GRADH = 0.D0
!        END IF
!     ELSE
!        GRADH = 0.D0
!     END IF


    !           CALCUL DES VALEURS AUX INTERFACES
    out_L%h = in_M%h - (DX/2.D0) * GRADH
    out_R%h = in_M%h + (DX/2.D0) * GRADH

    out_L%qx = in_M%qx - (DX/2.D0) * GRADQx
    out_R%qx = in_M%qx + (DX/2.D0) * GRADQx

    out_L%xi = in_M%xi 
    out_R%xi = in_M%xi 




  END SUBROUTINE MUSCL_3EQ



  ELEMENTAL SUBROUTINE MUSCL_o2c(in_L, in_M, in_R, Dx, out_L, out_R)

    IMPLICIT NONE

    TYPE(STATE), INTENT(in) :: in_L
    TYPE(STATE), INTENT(in) :: in_M
    TYPE(STATE), INTENT(in) :: in_R

    REAL(kind = 8), INTENT(in) :: Dx

    TYPE(STATE), INTENT(out) :: out_L
    TYPE(STATE), INTENT(out) :: out_R

    REAL(kind = 8) :: xiG
    REAL(kind = 8) :: xiM
    REAL(kind = 8) :: xiD

    REAL(kind = 8) :: DH
    REAL(kind = 8) :: DQx
    REAL(kind = 8) :: Dxi

    ! Initialisation de xi

!     xiG = in_L%h + in_L%topo
!     xiM = in_M%h + in_M%topo
!     xiD = in_R%h + in_R%topo

    xiG = in_L%xi
    xiM = in_M%xi
    xiD = in_R%xi
    ! Calcul du gradient du debit en hu

    DQx = f_minmod((in_M%qx-in_L%qx)/Dx , (in_R%qx-in_M%qx)/Dx )
    out_L%qx = in_M%qx - (Dx/2.0D0)*DQx
    out_R%qx = in_M%qx + (Dx/2.0D0)*DQx

    ! Calcul du gradient du debit en h

    DH = f_minmod((in_M%h-in_L%h)/Dx , (in_R%h-in_M%h)/Dx )
    out_L%h = in_M%h - (Dx/2.0D0)*DH
    out_R%h = in_M%h + (Dx/2.0D0)*DH

    ! Calcul du gradient pour la reconstruction de h + z

    Dxi = f_minmod((xiM-xiG)/Dx , (xiD-xiM)/Dx )
    out_L%xi = xiM - (Dx/2.0D0)*Dxi
    out_R%xi = xiM + (Dx/2.0D0)*Dxi

  END SUBROUTINE MUSCL_o2c

  SUBROUTINE MUSCL_V4(in_L2, in_L1, in_M, in_R1, in_R2, out_L, out_R)

    IMPLICIT NONE
    TYPE(STATE), INTENT(in) :: in_L2
    TYPE(STATE), INTENT(in) :: in_L1
    TYPE(STATE), INTENT(in) :: in_M
    TYPE(STATE), INTENT(in) :: in_R1   
    TYPE(STATE), INTENT(in) :: in_R2

    !REAL(kind = 8), INTENT(in) :: Dx

    TYPE(STATE), INTENT(out) :: out_L
    TYPE(STATE), INTENT(out) :: out_R

    !REAL(kind = 8) :: xiG
    !REAL(kind = 8) :: xiM
    !REAL(kind = 8) :: xiD

    !REAL(kind = 8) :: PH
    !REAL(kind = 8) :: Pxi
    !REAL(kind = 8) :: PQx
    !REAL(kind = 8) :: DH
    !REAL(kind = 8) :: Dxi
    REAL(kind = 8) :: DHp
    REAL(kind = 8) :: DHm
    REAL(kind = 8) :: DHp_lim
    REAL(kind = 8) :: DHm_lim

    REAL(kind = 8) :: DQxp
    REAL(kind = 8) :: DQxm
    REAL(kind = 8) :: DQxp_lim
    REAL(kind = 8) :: DQxm_lim

    REAL(kind = 8) :: DXIp
    REAL(kind = 8) :: DXIm

    REAL(kind = 8) :: DXIp_lim
    REAL(kind = 8) :: DXIm_lim

    REAL(kind = 8), PARAMETER :: beta = 1.0D0/3.0D0
    REAL(kind = 8), PARAMETER :: xi_c = -1.0D0/10.0D0
    REAL(kind = 8), PARAMETER :: xi_d = -1.0D0/15.0D0

!     PH = (in_M%h - in_L1%h) * (in_R1%h - in_M%h)
!     Pxi  = (in_M%xi - in_L1%xi) * (in_R1%xi - in_M%xi)
!     PQx = (in_M%qx - in_L1%qx) * (in_R1%qx - in_M%qx)
!     PQy = (in_M%qy - in_L1%qy) * (in_R1%qy - in_M%qy)

    ! Calcul du gradient du debit en hu
    !    IF(PQx <= 0.0D0) THEN
    !       DQxp = 0.0D0
    !       DQxm = 0.0D0
    !    ELSE
    DQxp = (1.0D0 - beta)*(in_R1%qx - in_M%qx) + beta*(in_M%qx - in_L1%qx) &
         + xi_c*(-in_L1%qx + 3.0D0*in_M%qx - 3.0D0*in_R1%qx + in_R2%qx) &
         + xi_d*(-in_L2%qx + 3.0D0*in_L1%qx - 3.0D0*in_M%qx + in_R1%qx)

    DQxm = (1.0D0 - beta)*(in_M%qx - in_L1%qx) + beta*(in_R1%qx - in_M%qx) &
         + xi_c*(-in_L2%qx + 3.0D0*in_L1%qx - 3.0D0*in_M%qx + in_R1%qx) &
         + xi_d*(-in_L1%qx + 3.0D0*in_M%qx - 3.0D0*in_R1%qx + in_R2%qx)
    
    IF((in_M%qx - in_L1%qx)*(in_R1%qx-in_M%qx) <= 0.0D0) THEN
       DQXp_lim = 0.0D0
       DQXm_lim = 0.0D0
    ELSE
       DQXp_lim = SIGN( MIN(2.0D0*ABS(in_M%qx - in_L1%qx),&
            2.0D0*ABS(in_R1%qx-in_M%qx), ABS(DQXp)), in_M%qx - in_L1%qx)

       DQXm_lim = SIGN( MIN(2.0D0*ABS(in_R1%qx - in_M%qx),&
            2.0D0*ABS(in_M%qx-in_L1%qx), ABS(DQXm)), in_R1%qx - in_M%qx)

    END IF

    out_L%qx = in_M%qx - (1.0D0/2.0D0)*DQxm_lim
    out_R%qx = in_M%qx + (1.0D0/2.0D0)*DQxp_lim

    ! Calcul du gradient du debit en h


    !     IF(PH <= 0.0D0) THEN

    !        DHp = 0.0D0
    !        DHm = 0.0D0

    !     ELSE

    DHp = (1.0D0 - beta)*(in_R1%h - in_M%h) + beta*(in_M%h - in_L1%h) &
         + xi_c*(-in_L1%h + 3.0D0*in_M%h - 3.0D0*in_R1%h + in_R2%h) &
         + xi_d*(-in_L2%h + 3.0D0*in_L1%h - 3.0D0*in_M%h + in_R1%h)

    DHm = (1.0D0 - beta)*(in_M%h - in_L1%h) + beta*(in_R1%h - in_M%h) &
         + xi_c*(-in_L2%h + 3.0D0*in_L1%h - 3.0D0*in_M%h + in_R1%h) &
         + xi_d*(-in_L1%h + 3.0D0*in_M%h - 3.0D0*in_R1%h + in_R2%h)

    ! Limiter
    IF((in_M%h - in_L1%h)*(in_R1%h-in_M%h) <= 0.0D0) THEN
       DHp_lim = 0.0D0
       DHm_lim = 0.0D0
    ELSE
       DHp_lim = SIGN( MIN(2.0D0*ABS(in_M%h - in_L1%h),&
            2.0D0*ABS(in_R1%h-in_M%h), ABS(DHp)), in_M%h - in_L1%h)

       DHm_lim = SIGN( MIN(2.0D0*ABS(in_R1%h - in_M%h),&
            2.0D0*ABS(in_M%h-in_L1%h), ABS(DHm)), in_R1%h - in_M%h)

    END IF


    out_L%h = in_M%h - (1.0D0/2.0D0)*DHm_lim
    out_R%h = in_M%h + (1.0D0/2.0D0)*DHp_lim


    ! Calcul du gradient pour la reconstruction de h + z

    !     IF(Pxi <= 0.0D0) THEN

    !        DXIp = 0.0D0
    !        DXIm = 0.0D0

    !     ELSE


    DXIp = (1.0D0 - beta)*(in_R1%xi - in_M%xi) + beta*(in_M%xi - in_L1%xi) &
         + xi_c*(-in_L1%xi + 3.0D0*in_M%xi - 3.0D0*in_R1%xi + in_R2%xi) &
         + xi_d*(-in_L2%xi + 3.0D0*in_L1%xi - 3.0D0*in_M%xi + in_R1%xi)

    DXIm = (1.0D0 - beta)*(in_M%xi - in_L1%xi) + beta*(in_R1%xi - in_M%xi) &
         + xi_c*(-in_L2%xi + 3.0D0*in_L1%xi - 3.0D0*in_M%xi + in_R1%xi) &
         + xi_d*(-in_L1%xi + 3.0D0*in_M%xi - 3.0D0*in_R1%xi + in_R2%xi)

    ! Limiter
    IF((in_M%xi - in_L1%xi)*(in_R1%xi-in_M%xi)<=0.0D0) THEN
       DXIp_lim = 0.0D0
       DXIm_lim = 0.0D0
    ELSE
       DXIp_lim = SIGN(MIN(2.0D0*ABS(in_M%xi - in_L1%xi),&
            2.0D0*ABS(in_R1%xi-in_M%xi), ABS(DXIp)), in_M%xi - in_L1%xi)

       DXIm_lim = SIGN(MIN(2.0D0*ABS(in_R1%xi - in_M%xi),&
            2.0D0*ABS(in_M%xi-in_L1%xi), ABS(DXIm)), in_R1%xi - in_M%xi)

    end IF

    out_L%xi = in_M%xi - (1.0D0/2.0D0)*DXIm_lim
    out_R%xi = in_M%xi + (1.0D0/2.0D0)*DXIp_lim


  END SUBROUTINE MUSCL_V4

  ELEMENTAL SUBROUTINE ENO_o2(in_L2, in_L1, in_M, in_R1, in_R2, Dx, out_L, out_R)

    IMPLICIT NONE

    TYPE(STATE), INTENT(in) :: in_L2
    TYPE(STATE), INTENT(in) :: in_L1
    TYPE(STATE), INTENT(in) :: in_M
    TYPE(STATE), INTENT(in) :: in_R1   
    TYPE(STATE), INTENT(in) :: in_R2

    REAL(kind = 8), INTENT(in) :: Dx

    TYPE(STATE), INTENT(out) :: out_L
    TYPE(STATE), INTENT(out) :: out_R

    !REAL(kind = 8) :: xiG
    !REAL(kind = 8) :: xiM
    !REAL(kind = 8) :: xiD

    REAL(kind = 8) :: DH
    REAL(kind = 8) :: DQx
    REAL(kind = 8) :: Dxi

    REAL(kind = 8) :: D2m_H
    REAL(kind = 8) :: D2m_Qx
    REAL(kind = 8) :: D2m_xi

    REAL(kind = 8) :: D2p_H
    REAL(kind = 8) :: D2p_Qx
    REAL(kind = 8) :: D2p_xi

!     ! Initialisation de xi

!     xiG = in_L%h + in_L%topo
!     xiM = in_M%h + in_M%topo
!     xiD = in_R%h + in_R%topo

    ! Calcul du gradient du debit en hu

    D2m_Qx = f_minmod((in_M%qx-2.0D0*in_L1%qx + in_L2%qx)/Dx**2 &
         , (in_R1%qx-2.0D0*in_M%qx+in_L1%qx)/Dx**2)

    D2p_Qx = f_minmod((in_R1%qx-2.0D0*in_M%qx + in_L1%qx)/Dx**2 &
         , (in_R2%qx-2.0D0*in_R1%qx+in_M%qx)/Dx**2)

    DQx = f_minmod((in_M%qx-in_L1%qx)/Dx + (Dx/2.0D0)*D2m_Qx , &
         (in_R1%qx-in_M%qx)/Dx - (Dx/2.0D0)*D2p_Qx )

    out_L%qx = in_M%qx - (Dx/2.0D0)*DQx
    out_R%qx = in_M%qx + (Dx/2.0D0)*DQx

    ! Calcul du gradient du debit en h

    D2m_H = f_minmod((in_M%h-2.0D0*in_L1%h + in_L2%h)/Dx**2 &
         , (in_R1%h-2.0D0*in_M%h+in_L1%h)/Dx**2)

    D2p_H = f_minmod((in_R1%h-2.0D0*in_M%h + in_L1%h)/Dx**2 &
         , (in_R2%h-2.0D0*in_R1%h+in_M%h)/Dx**2)

    DH = f_minmod((in_M%h-in_L1%h)/Dx + (Dx/2.0D0)*D2m_H , &
         (in_R1%h-in_M%h)/Dx - (Dx/2.0D0)*D2p_H )

    out_L%h = in_M%h - (Dx/2.0D0)*DH
    out_R%h = in_M%h + (Dx/2.0D0)*DH

    ! Calcul du gradient pour la reconstruction de h + z

    D2m_xi = f_minmod((in_M%xi-2.0D0*in_L1%xi + in_L2%xi)/Dx**2 &
         , (in_R1%xi-2.0D0*in_M%xi+in_L1%xi)/Dx**2)

    D2p_xi = f_minmod((in_R1%xi-2.0D0*in_M%xi + in_L1%xi)/Dx**2 &
         , (in_R2%xi-2.0D0*in_R1%xi+in_M%xi)/Dx**2)

    Dxi = f_minmod((in_M%xi-in_L1%xi)/Dx + (Dx/2.0D0)*D2m_xi , &
         (in_R1%xi-in_M%xi)/Dx - (Dx/2.0D0)*D2p_xi )
    out_L%xi = in_M%xi - (Dx/2.0D0)*Dxi
    out_R%xi = in_M%xi + (Dx/2.0D0)*Dxi

  END SUBROUTINE ENO_o2

  ELEMENTAL SUBROUTINE ENOmod_o2(in_L2, in_L1, in_M, in_R1, in_R2, Dx, out_L, out_R)

    IMPLICIT NONE

    TYPE(STATE), INTENT(in) :: in_L2
    TYPE(STATE), INTENT(in) :: in_L1
    TYPE(STATE), INTENT(in) :: in_M
    TYPE(STATE), INTENT(in) :: in_R1   
    TYPE(STATE), INTENT(in) :: in_R2

    REAL(kind = 8), INTENT(in) :: Dx

    TYPE(STATE), INTENT(out) :: out_L
    TYPE(STATE), INTENT(out) :: out_R

    !REAL(kind = 8) :: xiG
    !REAL(kind = 8) :: xiM
    !REAL(kind = 8) :: xiD

    REAL(kind = 8) :: D_eno_H
    REAL(kind = 8) :: D_eno_Qx
    REAL(kind = 8) :: D_eno_xi

    REAL(kind = 8) :: D_mm_H
    REAL(kind = 8) :: D_mm_Qx
    REAL(kind = 8) :: D_mm_xi

    REAL(kind = 8) :: DH
    REAL(kind = 8) :: DQx
    REAL(kind = 8) :: Dxi

    REAL(kind = 8) :: D2m_H
    REAL(kind = 8) :: D2m_Qx
    REAL(kind = 8) :: D2m_xi

    REAL(kind = 8) :: D2p_H
    REAL(kind = 8) :: D2p_Qx
    REAL(kind = 8) :: D2p_xi

    ! Calcul du gradient du debit en hu

    D2m_Qx = f_minmod((in_M%qx-2.0D0*in_L1%qx + in_L2%qx)/(Dx**2) &
         , (in_R1%qx-2.0D0*in_M%qx+in_L1%qx)/(Dx**2))

    D2p_Qx = f_minmod((in_R1%qx-2.0D0*in_M%qx + in_L1%qx)/Dx**2 &
         , (in_R2%qx-2.0D0*in_R1%qx+in_M%qx)/Dx**2)

    D_eno_Qx = f_minmod((in_M%qx-in_L1%qx)/Dx + (Dx/2.0D0)*D2m_Qx , &
         (in_R1%qx-in_M%qx)/Dx - (Dx/2.0D0)*D2p_Qx )

!     D_eno_Qx = f_minmod((in_M%qx-in_L1%qx)/Dx , &
!          (in_R1%qx-in_M%qx)/Dx)

    D_mm_Qx = f_minmod((in_M%qx-in_L1%qx)/Dx , &
         (in_R1%qx-in_M%qx)/Dx)

    DQx = f_minmod(D_eno_Qx, 1.8D0*D_mm_Qx)
    !DQx = D_mm_Qx

    out_L%qx = in_M%qx - (Dx/2.0D0)*DQx
    out_R%qx = in_M%qx + (Dx/2.0D0)*DQx


    ! Calcul du gradient du debit en h

    D2m_H = f_minmod((in_M%h - 2.0D0*in_L1%h + in_L2%h)/Dx**2 &
         , (in_R1%h - 2.0D0*in_M%h + in_L1%h)/Dx**2)

    D2p_H = f_minmod((in_R1%h - 2.0D0*in_M%h + in_L1%h)/Dx**2 &
         , (in_R2%h - 2.0D0*in_R1%h + in_M%h)/Dx**2)

     D_eno_H = f_minmod((in_M%h - in_L1%h)/Dx + (Dx/2.0D0)*D2m_H , &
          (in_R1%h - in_M%h)/Dx - (Dx/2.0D0)*D2p_H )

!     D_eno_H = f_minmod((in_M%h - in_L1%h)/Dx , &
!          (in_R1%h - in_M%h)/Dx )

    D_mm_H = f_minmod((in_M%h - in_L1%h)/Dx , &
         (in_R1%h - in_M%h)/Dx )

    DH = f_minmod(D_eno_H, 1.8D0*D_mm_H)
    !DH = D_mm_H

    out_L%h = in_M%h - (Dx/2.0D0)*DH
    out_R%h = in_M%h + (Dx/2.0D0)*DH

    ! Calcul du gradient pour la reconstruction de h + z

    D2m_xi = f_minmod((in_M%xi - 2.0D0*in_L1%xi + in_L2%xi)/Dx**2 &
         , (in_R1%xi - 2.0D0*in_M%xi + in_L1%xi)/Dx**2)

    D2p_xi = f_minmod((in_R1%xi - 2.0D0*in_M%xi + in_L1%xi)/Dx**2 &
         , (in_R2%xi - 2.0D0*in_R1%xi + in_M%xi)/Dx**2)

    D_eno_xi = f_minmod((in_M%xi - in_L1%xi)/Dx + (Dx/2.0D0)*D2m_xi , &
         (in_R1%xi - in_M%xi)/Dx - (Dx/2.0D0)*D2p_xi )
!     D_eno_xi =  f_minmod((in_M%xi - in_L1%xi)/Dx  , &
!          (in_R1%xi - in_M%xi)/Dx )

    D_mm_xi = f_minmod((in_M%xi - in_L1%xi)/Dx  , &
         (in_R1%xi - in_M%xi)/Dx )

    Dxi = f_minmod(D_eno_xi, 1.8D0*D_mm_xi)
    !Dxi = D_mm_xi

    out_L%xi = in_M%xi - (Dx/2.0D0)*Dxi
    out_R%xi = in_M%xi + (Dx/2.0D0)*Dxi

  END SUBROUTINE ENOmod_o2





  ELEMENTAL SUBROUTINE MUSCL_o4a(H_im2, H_im1, H_i, H_ip1, H_ip2, Incr_H_L, Incr_H_R, r_b1, i_var)

    IMPLICIT NONE

    REAL(KIND = 8), INTENT(in) :: H_im2
    REAL(KIND = 8), INTENT(in) :: H_im1
    REAL(KIND = 8), INTENT(in) :: H_i
    REAL(KIND = 8), INTENT(in) :: H_ip1
    REAL(KIND = 8), INTENT(in) :: H_ip2

    REAL(KIND = 8), INTENT(in) :: r_b1
    REAL(KIND = 8), INTENT(out) :: Incr_H_L
    REAL(KIND = 8), INTENT(out) :: Incr_H_R

    !REAL(KIND = 8), INTENT(in) :: Dx
    REAL(KIND = 8) :: vv_L   
    REAL(KIND = 8) :: vv_R
    REAL(KIND = 8) :: Delta_star_L
    REAL(KIND = 8) :: Delta_star_R
    REAL(KIND = 8) :: ss

    REAL(KIND = 8) :: QL    
    REAL(KIND = 8) :: QM
    REAL(KIND = 8) :: QR
    REAL(KIND = 8) :: drt_top

    REAL(kind = 8) :: PH
    !INTEGER, INTENT(in) :: direction
    LOGICAL, INTENT(in) :: i_var

    ! Produit pour la Monotonie

    PH = (H_i - H_im1) * (H_ip1 - H_i)

    ! Evaluation of the left and right increments for the water height

    IF(PH <= 0.0D0) THEN
       ! H_i is a local extremum

       Incr_H_L = cte__0_d_r__8
       Incr_H_R = cte__0_d_r__8

    ELSE

       ! Evaluation of Delta_star_L (i.e. Delta_{i-1/2,j} - Delta**3_{i-1/2,j}/6)
       QL = H_im1 - H_im2
       QM = H_i - H_im1
       QR = H_ip1 - H_i

       CALL Compute_Delta_star(QL, QM, QR, Delta_star_L)

       ! Evaluation of Delta_star_R (i.e. Delta_{i+1/2,j} - Delta**3_{i+1/2,j}/6)
       QL = H_i - H_im1
       QM = H_ip1 - H_i
       QR = H_ip2 - H_ip1

       CALL Compute_Delta_star(QL, QM, QR, Delta_star_R)

       ss = SIGN(cte__1_d_0_r__8, Delta_star_L)
       vv_L = ss*MAX(cte__0_d_r__8, MIN(ABS(Delta_star_L), r_b1*ss*Delta_star_R))

       ss = SIGN(cte__1_d_0_r__8, Delta_star_R)
       vv_R = ss*MAX(cte__0_d_r__8, MIN(ABS(Delta_star_R), r_b1*ss*Delta_star_L))


       Incr_H_R =   (      vv_L + cte__2_d_0_r__8*vv_R )/cte__6_d_0_r__8
       Incr_H_L = - (cte__2_d_0_r__8*vv_L +       vv_R )/cte__6_d_0_r__8

       ! Projection step

       SELECT CASE(i_var)

       CASE(.TRUE.)
          drt_top = Incr_H_L + Incr_H_R - H_i

          IF(drt_top > cte__0_d_r__8) THEN
             !PRINT*,'limitation, case 1'
             Incr_H_L =  H_i /(cte__1_d_0_r__8 + Incr_H_R/Incr_H_L)
             Incr_H_R =  H_i - Incr_H_L
          END IF

          IF(drt_top < cte__0_d_r__8.AND.(Incr_H_L < -H_i).AND.(Incr_H_R >= Incr_H_L)) THEN
             !PRINT*,'limitation, case 2'

             Incr_H_R =  - (Incr_H_R/Incr_H_L)*H_i
             Incr_H_L =  -  H_i
          END IF

          IF(drt_top < cte__0_d_r__8.AND.(Incr_H_R < -H_i).AND.(Incr_H_L >= Incr_H_R)) THEN
             !PRINT*,'limitation, case 3'

             Incr_H_L =  - (Incr_H_L/Incr_H_R)*H_i
             Incr_H_R =  -  H_i
          END IF

       CASE default

          ! nothing is done, no further limitation

       END SELECT

    END IF

  END SUBROUTINE MUSCL_o4a

  ELEMENTAL SUBROUTINE MUSCL_o4b(Hi_m2, Hi_m1, Hi, Hi_p1, Hi_p2, Hi_p3, H_left, H_right, limiter, b1)

    IMPLICIT NONE

    ! Compute the four (third) order reconstruction TVD-MUSCL interface values q(i+1/2)- 
    ! and q(i+1/2)+ from q(i-1), q(i), q(i+1) and q(i+2) using Yamamoto - Computers Fluids (1993)

    REAL(kind = 8), INTENT(in) :: Hi_m2
    REAL(kind = 8), INTENT(in) :: Hi_m1
    REAL(kind = 8), INTENT(in) :: Hi
    REAL(kind = 8), INTENT(in) :: Hi_p1
    REAL(kind = 8), INTENT(in) :: Hi_p2
    REAL(kind = 8), INTENT(in) :: Hi_p3
    !REAL(kind = 8), INTENT(in) :: Dx
    REAL(KIND = 8), INTENT(in) :: b1
    REAL(kind = 8), INTENT(out) :: H_left
    REAL(kind = 8), INTENT(out) ::  H_right
    REAL(kind = 8) :: Prod_left
    REAL(kind = 8) :: Prod_right
    REAL(kind = 8) :: QL
    REAL(kind = 8) ::  QM
    REAL(kind = 8) :: QR
    REAL(kind = 8) :: rr_1
    REAL(kind = 8) :: rr_2
    REAL(kind = 8) :: inv_rr_1
    REAL(kind = 8) :: inv_rr_2
    REAL(kind = 8) :: phi_2
    REAL(kind = 8) :: inv_phi_1
    REAL(kind = 8) :: inv_phi_2
    REAL(kind = 8) :: phi_1
    REAL(kind = 8) :: Delta_star_L
    REAL(kind = 8) :: Delta_star_M
    REAL(kind = 8) :: Delta_star_R
    REAL(kind = 8) :: ss
    REAL(kind = 8) :: vv_1
    REAL(kind = 8) :: vv_2
    INTEGER, INTENT(in) :: limiter


    ! Produit pour la Monotonie
    !PRINT*,'MUSCL Ordre 4'

    !print*,'DANS MUSCL 4'
    Prod_left = (Hi - Hi_m1) * (Hi_p1 - Hi)
    Prod_right = (Hi_p1 - Hi) * (Hi_p2 - Hi_p1)

    QL = Hi - Hi_m1
    QM = Hi_p1 - Hi
    QR = Hi_p2 - Hi_p1



    CALL Compute_Delta_star(QL, QM, QR, Delta_star_M)


    ! ########
    QL = Hi_m1 - Hi_m2
    QM = Hi - Hi_m1
    QR = Hi_p1 - Hi


    CALL Compute_Delta_star(QL, QM, QR, Delta_star_L)

    ! #########

    QL = Hi_p1 - Hi
    QM = Hi_p2 - Hi_p1
    QR = Hi_p3 - Hi_p2



    CALL Compute_Delta_star(QL, QM, QR, Delta_star_R)

    ! #########

    rr_1 = Delta_star_M/Delta_star_L
    inv_rr_1 = Delta_star_L/Delta_star_M


    IF(Prod_left <= cte__0_d_r__8) THEN

       ! Hi is a local extremum
       phi_1 = cte__0_d_r__8
       inv_phi_1 = cte__0_d_r__8

    ELSE

       ! slope limiter
       SELECT CASE(limiter)

       CASE(0) ! Minmod for wet/dry


       CASE(cte__1) ! classical minmod

          phi_1 = minmod(b1*rr_1)
          inv_phi_1 = minmod(b1*inv_rr_1)

       CASE(cte__2)

          phi_1 = van_albada(b1*rr_1)
          inv_phi_1 = van_albada(b1*inv_rr_1)

       CASE(cte__3)
          phi_1 = van_leer(b1*rr_1)
          inv_phi_1 = van_leer(b1*inv_rr_1)

       CASE(cte__4)
          phi_1 = Arora_Roe(b1*rr_1, cte__0_d_8_r__8)
          inv_phi_1 = arora_roe(b1*inv_rr_1, cte__0_d_8_r__8)

       CASE(cte__5)
          phi_1 = superbee(b1*rr_1)
          inv_phi_1 = superbee(b1*inv_rr_1)

       CASE(cte__6)
          phi_1 = Arora_Roe_NonLin(b1*rr_1, cte__0_d_8_r__8)
          inv_phi_1 = Arora_Roe_NonLin(b1*inv_rr_1, cte__0_d_8_r__8)


       END SELECT

    END IF

    SELECT CASE(limiter)

    CASE(cte__1,cte__2,cte__3,cte__4,cte__5,cte__6)
       H_left = Hi + (phi_1 * Delta_star_L + cte__2_d_0_r__8 * inv_phi_1 * Delta_star_M)/cte__6_d_0_r__8
       

    CASE(cte__0)
       ss = SIGN(cte__1_d_0_r__8, Delta_star_L)
       vv_1 = ss*MAX(cte__0_d_r__8, MIN(ABS(Delta_star_L), b1*ss*Delta_star_M))

       ss = SIGN(cte__1_d_0_r__8, Delta_star_M)
       vv_2 = ss*MAX(cte__0_d_r__8, MIN(ABS(Delta_star_M), b1*ss*Delta_star_L))


       H_left =  Hi + (   vv_1 + cte__2_d_0_r__8* vv_2   )/cte__6_d_0_r__8
       !print*,'Hletf', H_left

    CASE default
       !PRINT*,'Error in the choice of limiter'

    END SELECT



    ! #########

    rr_2 = Delta_star_R/Delta_star_M
    !print*,'rr2',rr_2
    inv_rr_2 = cte__1_d_0_r__8/rr_2

    IF(Prod_right <= cte__0_d_r__8) THEN

       ! Hi+1 is a local extremum

       phi_2 = cte__0_d_r__8
       inv_phi_2 = cte__0_d_r__8

    ELSE
       SELECT CASE(limiter)
       CASE(1)

          phi_2 = minmod(b1*rr_2)
          inv_phi_2 = minmod(b1*inv_rr_2)

       CASE(2)

          phi_2 = van_albada(b1*rr_2)
          inv_phi_2 = van_albada(b1*inv_rr_2)

       CASE(3)
          phi_2 = van_leer(b1*rr_2)
          inv_phi_2 = van_leer(b1*inv_rr_2)

       CASE(4)
          phi_2 = Arora_Roe(b1*rr_2, cte__0_d_8_r__8)
          inv_phi_2 = arora_roe(b1*inv_rr_2, cte__0_d_8_r__8)

       CASE(5)
          phi_2 = superbee(b1*rr_2)
          inv_phi_2 = superbee(b1*inv_rr_2)

       CASE(6)
          phi_2 = Arora_Roe_NonLin(b1*rr_2, cte__0_d_8_r__8)
          inv_phi_2 = Arora_Roe_NonLin(b1*inv_rr_2, cte__0_d_8_r__8)

       END SELECT
    END IF


    SELECT CASE(limiter)

    CASE(1,2,3,4,5,6)
       H_right = Hi_p1 - (cte__2_d_0_r__8 * phi_2 * Delta_star_M + inv_phi_2 * Delta_star_R)/cte__6_d_0_r__8

    CASE(0)
       ss = SIGN(cte__1_d_0_r__8, Delta_star_M)
       vv_1 = ss*MAX(cte__0_d_r__8, MIN(ABS(Delta_star_M), b1*ss*Delta_star_R))

       ss = SIGN(cte__1_d_0_r__8, Delta_star_R)
       vv_2 = ss*MAX(cte__0_d_r__8, MIN(ABS(Delta_star_R), b1*ss*Delta_star_M))


       H_right =  Hi_p1 - (   cte__2_d_0_r__8*vv_1 + vv_2   )/cte__6_d_0_r__8
       !print*,'Hright', H_right
       

    CASE default
       !PRINT*,'Error in the choice of limiter'

    END SELECT

  END SUBROUTINE MUSCL_o4b

  ELEMENTAL SUBROUTINE Compute_Delta_star(QL, QM, QR, Delta_star)

    IMPLICIT NONE

    REAL(kind = 8), INTENT(in) :: QL, QM, QR
    REAL(kind = 8), INTENT(out) :: Delta_star
    REAL(kind = 8) :: ss, Delta_iL, Delta_iM, Delta_iR, Delta_O4 


    ss = SIGN(cte__1_d_0_r__8, QL)
    Delta_iL = ss*MAX(cte__0_d_r__8, MIN(ABS(QL), 2.0*ss*QM, 2.0*ss*QR))

    ss = SIGN(cte__1_d_0_r__8, QM)
    Delta_iM = ss*MAX(cte__0_d_r__8, MIN(ABS(QM), 2.0*ss*QL, 2.0*ss*QR))

    ss = SIGN(cte__1_d_0_r__8, QR)
    Delta_iR = ss*MAX(cte__0_d_r__8, MIN(ABS(QR), 2.0*ss*QL, 2.0*ss*QM))

    Delta_O4 = Delta_iR - 2.0*Delta_iM + Delta_iL

    Delta_star = QM - Delta_O4/cte__6_d_0_r__8

  END SUBROUTINE Compute_Delta_star


  ELEMENTAL FUNCTION van_leer(rr)

    IMPLICIT NONE


    REAL(KIND = 8), INTENT(in) :: rr
    REAL(KIND = 8) :: Van_leer

    van_leer = (rr + ABS(rr))/(cte__1_d_0_r__8 + rr)

  END FUNCTION van_leer

  ELEMENTAL FUNCTION minmod(rr)

    IMPLICIT NONE


    REAL(kind = 8), INTENT(in) :: rr
    REAL(kind = 8) :: minmod

    minmod = MAX(cte__0_d_r__8, MIN(cte__1_d_0_r__8, rr))

  END FUNCTION minmod

  ELEMENTAL FUNCTION f_minmod(xx, yy)

    IMPLICIT NONE

    REAL(kind = 8), INTENT(in) :: xx
    REAL(kind = 8), INTENT(in) :: yy
    REAL(kind = 8) :: f_minmod

    IF(xx >= 0.0D0 .AND. yy >= 0.0D0) THEN
       f_minmod = MIN(xx,yy)
    ELSEIF(xx <= 0.0D0 .AND. yy <= 0.0D0) THEN
       f_minmod = MAX(xx,yy)
    ELSE
       f_minmod = 0.0D0
    END IF

  END FUNCTION f_minmod

  ELEMENTAL FUNCTION superbee(rr)

    IMPLICIT NONE


    REAL(kind = 8), INTENT(in) :: rr
    REAL(kind = 8) :: superbee

    superbee = MAX(cte__0_d_r__8, MIN(cte__1_d_0_r__8, 2.0*rr), MIN(cte__2_d_0_r__8, rr))

  END FUNCTION superbee

  ELEMENTAL FUNCTION van_albada(rr)

    IMPLICIT NONE


    REAL(kind = 8), INTENT(in) :: rr
    REAL(kind = 8) :: van_albada

    van_albada = MAX(0.0, (rr + rr**2)/(cte__1_d_0_r__8 + rr**2))

  END FUNCTION van_albada

  ELEMENTAL FUNCTION Arora_Roe(rr, CFL)

    IMPLICIT NONE


    REAL(kind = 8), INTENT(in) :: rr
    REAL(kind = 8), INTENT(in) :: CFL
    REAL(kind = 8) :: arora_roe, s1, s2, phi_max

    s1 = cte__2_d_0_r__8/CFL
    phi_max = cte__2_d_0_r__8/(cte__1_d_0_r__8-CFL)
    s2 = (cte__1_d_0_r__8 + CFL)/3.0D0


    Arora_Roe = MAX(cte__0_d_r__8, MIN(s1*rr, cte__1_d_0_r__8 + s2*(rr-cte__1_d_0_r__8), phi_max))

  END FUNCTION Arora_Roe

  ELEMENTAL FUNCTION Arora_Roe_NonLin(rr, CFL)

    IMPLICIT NONE

    REAL(kind = 8), INTENT(in) :: rr
    REAL(kind = 8), INTENT(in) :: CFL
    REAL(kind = 8) :: s1, s2, phi_max
    REAL(kind = 8) :: Arora_Roe_NonLin

    s1 = cte__2_d_0_r__8
    phi_max = cte__2_d_0_r__8
    s2 = (cte__1_d_0_r__8 + CFL)/3.0D0


    Arora_Roe_Nonlin = MAX(cte__0_d_r__8, MIN(s1*rr, cte__1_d_0_r__8 + s2*(rr-cte__1_d_0_r__8), phi_max))

  END FUNCTION Arora_Roe_NonLin

  ELEMENTAL FUNCTION fPSI( H, Q, ZF ) 
    !USE fortranf
    !USE M_CST
    IMPLICIT NONE
    REAL (KIND = 8) :: fpsi
    REAL (KIND = 8), INTENT ( IN ) :: H
    REAL (KIND = 8), INTENT ( IN ) :: Q
    REAL (KIND = 8), INTENT ( IN ) :: ZF

    fpsi = Q**2/(2.0D0*H** 2)+g_r__8*(H+ZF)

  END FUNCTION fPSI

  ELEMENTAL SUBROUTINE MUSCL_o4c(H_im2, H_im1, H_i, H_ip1, H_ip2, Incr_H_L, Incr_H_R, r_b1)

    IMPLICIT NONE

    REAL(KIND = 8), INTENT(in) :: H_im2
    REAL(KIND = 8), INTENT(in) :: H_im1
    REAL(KIND = 8), INTENT(in) :: H_i
    REAL(KIND = 8), INTENT(in) :: H_ip1
    REAL(KIND = 8), INTENT(in) :: H_ip2

    REAL(KIND = 8), INTENT(in) :: r_b1
    REAL(KIND = 8), INTENT(out) :: Incr_H_L
    REAL(KIND = 8), INTENT(out) :: Incr_H_R

    !REAL(KIND = 8), INTENT(in) :: Dx
    REAL(KIND = 8) :: vv_L   
    REAL(KIND = 8) :: vv_R
    REAL(KIND = 8) :: Delta_star_L
    REAL(KIND = 8) :: Delta_star_R
    REAL(KIND = 8) :: ss

    REAL(KIND = 8) :: QL    
    REAL(KIND = 8) :: QM
    REAL(KIND = 8) :: QR
    !REAL(KIND = 8) :: drt_top

    !REAL(kind = 8) :: PH

    ! Evaluation of Delta_star_L (i.e. Delta_{i-1/2,j} - Delta**3_{i-1/2,j}/6)
    QL = H_im1 - H_im2
    QM = H_i - H_im1
    QR = H_ip1 - H_i

    CALL Compute_Delta_star(QL, QM, QR, Delta_star_L)

    ! Evaluation of Delta_star_R (i.e. Delta_{i+1/2,j} - Delta**3_{i+1/2,j}/6)
    QL = H_i - H_im1
    QM = H_ip1 - H_i
    QR = H_ip2 - H_ip1

    CALL Compute_Delta_star(QL, QM, QR, Delta_star_R)

    !PH = (H_i - H_im1)*(H_ip1 - H_i)

    !IF(PH <= 0.0D0) THEN
    !   vv_L = 0.0D0
    !ELSE 
       ss = SIGN(cte__1_d_0_r__8, Delta_star_L)
       vv_L = ss*MAX(cte__0_d_r__8, MIN(ABS(Delta_star_L), r_b1*ss*Delta_star_R))
    !END IF

    !IF(PH <= 0.0D0) THEN
    !   vv_R = 0.0D0
    !ELSE 
       ss = SIGN(cte__1_d_0_r__8, Delta_star_R)
       vv_R = ss*MAX(cte__0_d_r__8, MIN(ABS(Delta_star_R), r_b1*ss*Delta_star_L))

       Incr_H_R =   (      vv_L + cte__2_d_0_r__8*vv_R )/cte__6_d_0_r__8
       Incr_H_L = - (cte__2_d_0_r__8*vv_L +       vv_R )/cte__6_d_0_r__8
    !END IF

  END SUBROUTINE MUSCL_o4c

END MODULE modMUSCL
