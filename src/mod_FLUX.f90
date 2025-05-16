!=========================================================================================
!         M_FLUX : computation of interface well-balanced fuxes at order 1, 2, 4 or 5 
!                 with interface (non-conservative) or in-cell (conservative) approaches ...
!         F.Marche - I3M - 22/02/2008 version         
!=========================================================================================

MODULE modFlux

  USE modTOOLS
  USE modSolver
  USE modMUSCL
  USE modWENO
  USE modCST

  PRIVATE
  PUBLIC :: Flux_HR

CONTAINS

  !=== Numerical fluxes computation for Hydrostatic reconstruction based schemes ===

  SUBROUTINE Flux_HR &
       (Tab, MESH, ZF, DATA, FLU, UMax)

    USE modCST

    IMPLICIT NONE

    TYPE(MeshDef), INTENT(IN) :: MESH
    TYPE(DATA_IN), INTENT(IN) :: DATA
    INTEGER :: Nbx
    !INTEGER, INTENT(IN) :: step
    !REAL(KIND = 8) :: udele1
    !REAL(KIND = 8) :: vdele1
    !REAL(KIND = 8) :: udele2
    !REAL(KIND = 8) :: vdele2

    !REAL(KIND = 8), DIMENSION(2) :: twinv

    REAL(KIND = 8), DIMENSION(0:DATA%Nbx ) :: u_en_i

    !REAL(KIND = 8) :: u_en_j


    CHARACTER(LEN =*), PARAMETER :: sub__name = "flux_HR"
    TYPE( STATE ), DIMENSION(-2:DATA%Nbx+3), INTENT(INOUT):: Tab   
    TYPE(VFLU), INTENT(INOUT) :: FLU
    REAL(kind = 8), DIMENSION(:), POINTER :: ZF

    REAL(kind = 8), INTENT(out) :: UMax

    TYPE(state), DIMENSION(0:DATA%Nbx+1)  :: Grad
    !TYPE(state), DIMENSION(0:DATA%Nbx+1) :: S_G 
    !TYPE(state), DIMENSION(0:DATA%Nbx+1) :: S_D
    TYPE(state), DIMENSION(0:DATA%Nbx+1) :: S_L 
    TYPE(state), DIMENSION(0:DATA%Nbx+1) :: S_R
    TYPE(state), DIMENSION(0:DATA%Nbx+1) :: Incr_L
    TYPE(state), DIMENSION(0:DATA%Nbx+1) :: Incr_R
    TYPE(state), DIMENSION(0:DATA%Nbx+1) :: Val_C
    TYPE( STATE ), DIMENSION(0:DATA%Nbx) :: FLU_X
    REAL(KIND=8), DIMENSION(3) :: TFLU

    REAL(kind = 8) :: UM
    !REAL(kind = 8) :: VM
    REAL(kind = 8) :: CM
    REAL(kind = 8) :: Dx
    REAL(kind = 8) :: Dxi
    REAL(kind = 8) :: Dxip
    REAL(kind = 8) :: Zint
    REAL(kind = 8) :: Source1
    REAL(kind = 8) :: Source2
    REAL(kind = 8) :: hir
    REAL(kind = 8) :: hil
    REAL(kind = 8) :: hic
    REAL(kind = 8) :: zir
    REAL(kind = 8) :: zil
    REAL(kind = 8) :: zic
    REAL(kind = 8) :: hzil
    REAL(kind = 8) :: hzir
    REAL(kind = 8) :: hzic
    REAL(kind = 8) :: L12
    REAL(kind = 8) :: L11
    REAL(kind = 8) :: L21
    REAL(kind = 8) :: L22

    REAL(kind = 8) :: R12
    REAL(kind = 8) :: R11
    REAL(kind = 8) :: R21
    REAL(kind = 8) :: R22
    REAL(kind = 8), DIMENSION(-2:2) :: ubar_cv_1
    REAL(kind = 8), DIMENSION(-2:2) :: ubar_cv_2
    REAL(kind = 8) :: CV1_L
    REAL(kind = 8) :: CV1_C
    REAL(kind = 8) :: CV1_R
    REAL(kind = 8) :: CV2_L
    REAL(kind = 8) :: CV2_C
    REAL(kind = 8) :: CV2_R

    REAL(kind = 8) :: SourceCi
    !REAL(kind = r__8) :: Prod_tvd
    !REAL(kind = r__8) :: PHZ
    INTEGER :: i
    INTEGER :: j
    LOGICAL :: change_dir
    LOGICAL :: i_var
    !LOGICAL :: l_var
    LOGICAL :: debug 
    TYPE(state) :: sg
    TYPE(state) :: sd
    TYPE(state) :: s
    TYPE(state) :: s1
    TYPE(state) :: s2

    Nbx = DATA%Nbx

    Tab%topo = ZF
    Tab%xi = Tab%h + Tab%topo

    debug = .FALSE.
    !debug = .TRUE.

    IF (debug) THEN
       PRINT*,'Entering M_Flux unit ... '
       PRINT*,'M_Flux : variable initialisation ...'
       PRINT*,'M_Flux : computing gradients and increments for high order reconstructions ...'
    END IF

    !WRITE(6,*)'space_order ', DATA%space_order
    SELECT CASE(DATA%space_order)

    CASE('V4') 

       DO i = 0, Nbx+1

          Dx = MESH%dx

          CALL MUSCL_V4(Tab(i-2), Tab(i-1), Tab(i), Tab(i+1), Tab(i+2), &
               S_L(i), S_R(i))

       END DO


    CASE('eno2') !ENO  ORDRE 2

       DO i = 0, Nbx+1

          Dx = MESH%dx

          CALL ENO_o2(Tab(i-2), Tab(i-1), Tab(i), Tab(i+1), Tab(i+2), Dx, &
               S_L(i), S_R(i))
       END DO


    CASE('eno2mod') !ENO modifie ORDRE 2

       DO i = 0, Nbx+1

          Dx = MESH%dx

          CALL ENOmod_o2(Tab(i-2), Tab(i-1), Tab(i), Tab(i+1), Tab(i+2), Dx, &
               S_L(i), S_R(i))
       END DO

    CASE('muscl2c') !MUSCL ORDRE 2 (Bouchut)

       !DO i = 0, Nbx+1
       !  S_L(i) = Tab(i)
       !  S_R(i) = Tab(i)
       !END DO

       DO i = 0, Nbx+1

          Dx = MESH%dx

          CALL MUSCL_o2c(Tab(i-1), Tab(i), Tab(i+1), Dx, &
               S_L(i), S_R(i))
       END DO

    CASE('muscl2a') ! MUSCL initial (seguin))

       IF (debug) THEN
          PRINT*,'M_flux : computing "in-cell approach" second order gradients ...'
       END IF

       DO i = 0, Nbx+1

          Dx = MESH%dx

          CALL MUSCL_o2a(Tab(i-1), Tab(i), Tab(i+1), &
               Grad(i)%h, Grad(i)%xi, Grad(i)%qx, Dx)
       END DO

    CASE('muscl2b')

       IF (debug) THEN
          PRINT*,'M_flux : computing "interface approach" second order values ...'
       END IF

       S_L%h = cte__0_d_r__8
       S_R%h = cte__0_d_r__8
       S_L%xi = cte__0_d_r__8
       S_R%xi = cte__0_d_r__8
       S_L%qx = cte__0_d_r__8
       S_R%qx= cte__0_d_r__8


       Dx = MESH%dx

       DO i = 0, Nbx

          CALL MUSCL_o2b(Tab(i-1)%h, Tab(i)%h, Tab(i+1)%h, Tab(i+2)%h, S_L(i)%h, S_R(i)%h, &
               DATA%limiter, DATA%CFL)

          CALL MUSCL_o2b(Tab(i-1)%xi, Tab(i)%xi, Tab(i+1)%xi, Tab(i+2)%xi, S_L(i)%xi, S_R(i)%xi, &
               DATA%limiter, DATA%CFL)

          CALL MUSCL_o2b(Tab(i-1)%qx, Tab(i)%qx, Tab(i+1)%qx, Tab(i+2)%qx, S_L(i)%qx, S_R(i)%qx, &
               DATA%limiter, DATA%CFL)

       END DO

    CASE('muscl4a')

       IF (debug) THEN
          PRINT*,'M_flux : computing "in-cell approach" fourth order increments ...'
       END IF

       Dx = MESH%dx

       DO i = 0, Nbx+1

          i_var = .TRUE.

          CALL MUSCL_o4a(Tab(i-2)%xi, Tab(i-1)%xi, Tab(i)%xi, Tab(i+1)%xi, Tab(i+2)%xi, &
               Incr_L(i)%xi, Incr_R(i)%xi, DATA%b1, i_var)

          CALL MUSCL_o4a(Tab(i-2)%h, Tab(i-1)%h, Tab(i)%h, Tab(i+1)%h, Tab(i+2)%h, &
               Incr_L(i)%h, Incr_R(i)%h, DATA%b1, i_var)

          i_var = .FALSE.

          CALL MUSCL_o4a(Tab(i-2)%qx, Tab(i-1)%qx, Tab(i)%qx, Tab(i+1)%qx, Tab(i+2)%qx, &
               Incr_L(i)%qx, Incr_R(i)%qx, DATA%b1, i_var)

       END DO

    CASE('muscl4c')

       Dx = MESH%dx

       DO i = 0, Nbx+1

          CALL MUSCL_o4c(Tab(i-2)%xi, Tab(i-1)%xi, Tab(i)%xi, Tab(i+1)%xi, Tab(i+2)%xi, &
               Incr_L(i)%xi, Incr_R(i)%xi, DATA%b1)

          CALL MUSCL_o4c(Tab(i-2)%h, Tab(i-1)%h, Tab(i)%h, Tab(i+1)%h, Tab(i+2)%h, &
               Incr_L(i)%h, Incr_R(i)%h, DATA%b1)

          CALL MUSCL_o4c(Tab(i-2)%qx, Tab(i-1)%qx, Tab(i)%qx, Tab(i+1)%qx, Tab(i+2)%qx, &
               Incr_L(i)%qx, Incr_R(i)%qx, DATA%b1)

       END DO

    CASE('muscl4b')

       IF (debug) THEN
          PRINT*,'M_flux : computing "interface approach" fourth order interface values ...'
       END IF

       Dx = MESH%dx

       DO i = 0, Nbx

          CALL MUSCL_o4b(Tab(i-2)%h, Tab(i-1)%h, Tab(i)%h, Tab(i+1)%h, Tab(i+2)%h, &
               Tab(i+3)%h, S_L(i)%h, S_R(i)%h, DATA%limiter, DATA%b1)

          CALL MUSCL_o4b(Tab(i-2)%xi, Tab(i-1)%xi, Tab(i)%xi, Tab(i+1)%xi, Tab(i+2)%xi, Tab(i+3)%xi,  &
               S_L(i)%xi,  S_R(i)%xi, DATA%limiter, DATA%b1)

          CALL MUSCL_o4b(Tab(i-2)%qx, Tab(i-1)%qx, Tab(i)%qx, Tab(i+1)%qx, Tab(i+2)%qx, & 
               Tab(i+3)%qx, S_L(i)%qx, S_R(i)%qx, DATA%limiter, DATA%b1)

       END DO

    CASE('muscl5')


       Dx = MESH%dx

       DO i = 0, Nbx

          CALL MUSCL_o5(Tab(i-2), Tab(i-1), Tab(i), Tab(i+1), Tab(i+2), &
               Tab(i+3), S_L(i), S_R(i))

       END DO

    CASE('mixte5')


       Dx = MESH%dx

       DO i = 0, Nbx

          CALL MUSCL_o5(Tab(i-2), Tab(i-1), Tab(i), Tab(i+1), Tab(i+2), &
               Tab(i+3), S_L(i), S_R(i))

          CALL WENO_CO5(Tab(i-2)%xi, Tab(i-1)%xi, Tab(i)%xi, Tab(i+1)%xi, Tab(i+2)%xi, &
               Val_C(i)%xi)

          CALL WENO_CO5(Tab(i-2)%h, Tab(i-1)%h, Tab(i)%h, Tab(i+1)%h, Tab(i+2)%h, &
               Val_C(i)%h)

          CALL WENO_CO5(Tab(i-2)%qx, Tab(i-1)%qx, Tab(i)%qx, Tab(i+1)%qx, Tab(i+2)%qx, &
               Val_C(i)%qx)

       END DO

    CASE('weno3')

       DO i = 0, Nbx+1

          CALL WENO3(Tab(i-1)%xi, Tab(i)%xi, Tab(i+1)%xi, S_L(i)%xi, &
               S_R(i)%xi)

          CALL WENO3(Tab(i-1)%h, Tab(i)%h, Tab(i+1)%h, &
               S_L(i)%h, S_R(i)%h)

          CALL WENO3(Tab(i-1)%qx, Tab(i)%qx, Tab(i+1)%qx, &
               S_L(i)%qx, S_R(i)%qx)


       END DO


    CASE('weno5')

       IF (debug) THEN
          PRINT*,'M_Flux : computing WENO five order increments ...'
       END IF


       Dx = MESH%dx

       DO i = 0, Nbx+1

          CALL Built_WENO(Tab(i-2)%xi, Tab(i-1)%xi, Tab(i)%xi, Tab(i+1)%xi, Tab(i+2)%xi, S_L(i)%xi, &
               Val_C(i)%xi, S_R(i)%xi)

          CALL Built_WENO(Tab(i-2)%h, Tab(i-1)%h, Tab(i)%h, Tab(i+1)%h, Tab(i+2)%h, &
               S_L(i)%h, Val_C(i)%h, S_R(i)%h)

          CALL Built_WENO(Tab(i-2)%qx, Tab(i-1)%qx, Tab(i)%qx, Tab(i+1)%qx, Tab(i+2)%qx, &
               S_L(i)%qx, Val_C(i)%qx, S_R(i)%qx)
       END DO

    CASE('weno5cv')

       IF (debug) THEN
          WRITE(UNIT=6, FMT=*)'M_Flux : computing WENO fifth order...'
          !WRITE(UNIT=6, FMT=*)'Reconstruction in characteristic variables'
       END IF


       Dx = MESH%dx

       DO i = 0, Nbx+1

          CALL eigen_matrix(L11, L12, L21, L22, Tab(i)%h, Tab(i)%qx, &
               R11, R12, R21, R22)

          ubar_cv_1(-2) = L11*Tab(i-2)%xi + L12*Tab(i-2)%qx
          ubar_cv_1(-1) = L11*Tab(i-1)%xi + L12*Tab(i-1)%qx 
          ubar_cv_1(0) = L11*Tab(i)%xi   + L12*Tab(i)%qx
          ubar_cv_1(1) = L11*Tab(i+1)%xi + L12*Tab(i+1)%qx
          ubar_cv_1(2) = L11*Tab(i+2)%xi + L12*Tab(i+2)%qx

          ubar_cv_2(-2) = L21*Tab(i-2)%xi + L22*Tab(i-2)%qx
          ubar_cv_2(-1) = L21*Tab(i-1)%xi + L22*Tab(i-1)%qx 
          ubar_cv_2(0) = L21*Tab(i)%xi   + L22*Tab(i)%qx
          ubar_cv_2(1) = L21*Tab(i+1)%xi + L22*Tab(i+1)%qx 
          ubar_cv_2(2) = L21*Tab(i+2)%xi + L22*Tab(i+2)%qx

          CALL Built_WENO_CV(ubar_cv_1, &
               CV1_L, CV1_C, CV1_R)

          CALL Built_WENO_CV(ubar_cv_2, &
               CV2_L, CV2_C, CV2_R)

          S_L(i)%xi = R11*cv1_l + R12*cv2_l
          Val_C(i)%xi = R11*cv1_c + R12*cv2_c 
          S_R(i)%xi = R11*cv1_r + R12*cv2_r
          S_L(i)%qx = R21*cv1_l + R22*cv2_l
          Val_C(i)%qx=R21*cv1_c+R22*cv2_c
          S_R(i)%qx = R21*cv1_r + R22*cv2_r
          S_L(i)%h = Tab(i)%xi - ZF(i)
          Val_C(i)%h = Tab(i)%xi- ZF(i)
          S_R(i)%h = Tab(i)%xi- ZF(i)

       END DO


    CASE ('godunov') ! nothing is done
       !print*,"Je rentre ici"
    CASE default 

       PRINT*,'M_Flux : bad order of accuracy choice (space_order) ...'
       STOP

    END SELECT ! select case i_order : initialisation
    

    IF (debug) THEN
       PRINT*,'M_Flux : begining of x-fluxes computation ...'
       PRINT*,'M_Flux : left and right interface states initialisation ...'
    END IF

    change_dir = .FALSE.

    DO i = 0, Nbx

       SELECT CASE(DATA%space_order) 

       CASE('godunov') ! First order

          IF (debug) THEN
             PRINT*,'M_Flux : first order states initialisation ...'
          END IF

          sg%h = Tab(i)%h   
          sd%h = Tab(i+1)%h 

          sg%xi = ZF(i) + Tab(i)%h  
          sd%xi = ZF(i+1) + Tab(i+1)%h 

          sg%qx = Tab(i)%qx 
          sd%qx = Tab(i+1)%qx

       CASE('muscl2b') ! Second order, interface ("non-conservative") approach 

          IF (debug) THEN
             PRINT*,'M_Flux : second order "interface approach" states initialisation ...'
          END IF

          sg%h = S_L(i)%h   
          sd%h = S_R(i)%h

          sg%xi = S_L(i)%xi
          sd%xi = S_R(i)%xi

          sg%qx = S_L(i)%qx
          sd%qx = S_R(i)%qx

       CASE('muscl2c', 'eno2mod','eno2', 'V4') 

          IF (debug) THEN
             PRINT*,'M_Flux : second order "interface approach" states initialisation ...'
          END IF

          sg%h = S_R(i)%h   
          sd%h = S_L(i+1)%h

          sg%xi = S_R(i)%xi
          sd%xi = S_L(i+1)%xi

          sg%qx = S_R(i)%qx
          sd%qx = S_L(i+1)%qx

       CASE('muscl2a') ! Second order, in-cell ("conservative") approach

          IF (debug) THEN
             PRINT*,'M_Flux : second order "in-cell approach" states initialisation ...'
          END IF

          Dxi = MESH%dx
          Dxip = MESH%dx

          sg%h = Tab(i)%h     + Dxi/2.0D0*Grad(i)%h
          sd%h = Tab(i+1)%h   - Dxip/2.0D0*Grad(i+1)%h

          sg%xi = ZF(i) + Tab(i)%h + Dxi/2.0D0*Grad(i)%xi ; 
          sd%xi = ZF(i+1) + Tab(i+1)%h - Dxip/2.0D0*Grad(i+1)%xi

          sg%qx = Tab(i)%qx + Dxi/2.0D0*Grad(i)%qx
          sd%qx = Tab(i+1)%qx - Dxip/2.0D0*Grad(i+1)%qx

       CASE('muscl4b', 'muscl5','mixte5') ! Fourth order, interface ("non-conservative") approach

          IF (debug) THEN
             PRINT*,'M_Flux : fourth order "interface approach" states initialisation ...'
          END IF

          sg%h = S_L(i)%h   
          sd%h = S_R(i)%h

          sg%xi = S_L(i)%xi
          sd%xi = S_R(i)%xi

          sg%qx = S_L(i)%qx
          sd%qx = S_R(i)%qx

       CASE('muscl4a','muscl4c') ! Fourth order, in-cell ("conservative") approach

          IF(debug) THEN
             PRINT*,' M_Flux : fourth order "in-cell approach" states initialisation ...'
          END IF

          sg%h = Tab(i)%h     + Incr_R(i)%h
          sd%h = Tab(i+1)%h   + Incr_L(i+1)%h 

          sg%xi = ZF(i)   + Tab(i)%h + Incr_R(i)%xi 
          sd%xi = ZF(i+1) + Tab(i+1)%h +  Incr_L(i+1)%xi


          sg%qx = Tab(i)%qx   +  Incr_R(i)%qx      
          sd%qx = Tab(i+1)%qx +  Incr_L(i+1)%qx  

       CASE('weno3','weno5', 'weno5cv')

          IF(debug) THEN
             PRINT*,'M_Flux : WENO five order states initialisation ...'
          END IF

          sg%h = S_R(i)%h
          sd%h = S_L(i+1)%h 

          sg%xi = S_R(i)%xi 
          sd%xi = S_L(i+1)%xi

          sg%qx = S_R(i)%qx      
          sd%qx = S_L(i+1)%qx  

       CASE default         
          PRINT*,'M_Flux : bad order of accuracy choice (i_order) ...'

       END SELECT

       IF (debug) THEN
          PRINT*, 'M_Flux : z reconstruction ...'
       END IF


       sg%topo = sg%xi - sg%h  
       sd%topo = sd%xi - sd%h 

       IF (debug) THEN
          PRINT*, 'M_Flux : c construction ...'
          PRINT*, 'M_Flux : x-step -  (u, v) construction ...'
       END IF

       IF(sg%h <= DATA%kappa) THEN 
          sg%h = DATA%kappa
          sg%u = 0.0D0      
       ELSEIF ( sg%h <= DATA%kappabis) THEN
          sg%u = 0.0D0      
       ELSE
          sg%u =   sg%qx / sg%h
       END IF

       IF(sd%h <= DATA%kappa) THEN 
          sd%h = DATA%kappa
          sd%u = 0.0D0   
       ELSEIF(sd%h <= DATA%kappabis) THEN
          sd%u = 0.0D0    
       ELSE
          sd%u =   sd%qx / sd%h
       END IF

       IF (debug) THEN
          PRINT*, 'M_Flux : hydrostatic reconstruction ...'
       END IF

       Zint = MAX(sg%topo, sd%topo)
       s1%h = MAX(0.0D0, sg%h + sg%topo - Zint)
       s2%h = MAX(0.0D0, sd%h + sd%topo - Zint)

       s1%u = sg%u
       s2%u = sd%u

       s1%c = SQRT(g_r__8*s1%h) !- SQRT(g*kappa) 
       s2%c = SQRT(g_r__8*s2%h) !- SQRT(g*kappa)    


       CM = cte__0_d_5_r__8 * ( s1%C + s2%C)
       UM = cte__0_d_5_r__8 * ( s1%U + s2%U)

       u_en_i (i ) = ABS( UM ) + CM

       IF (debug) THEN
          PRINT*, 'M_Flux : computing interface states, calling the Riemann solver ...'
       END IF

       SELECT CASE(DATA%Riemann_Solv)

       CASE('VFROE')
          CALL Solver_VFRoe_ncv( s1, s2, S)

       CASE('EXACT')
          CALL Solver_Godunov( s1, s2, S)

       CASE('HLLC')
          CALL Solver_HLLC(s1, s2, TFLU)
          FLU_X(i)%h  = TFLU(1)
          FLU_X(i)%qx  = TFLU(2)


       END SELECT

       IF (debug) THEN
          PRINT*, 'M_Flux : computing interface fluxes ...'
       END IF

       SELECT CASE(DATA%Riemann_Solv)

       CASE('VFROE', 'EXACT')

          FLU_X(i)%h  = s%h*s%u 

          IF (debug) THEN
             PRINT *, "i,j==", i, j, " s1==", s1, " s2==", s2, " S==", S, " -> FLU_X==", FLU_X(i)%h
          END IF

          FLU_X(i)%qx  = s%h*s%u**2 + g_r__8*s%h**2/2.0D0

       END SELECT

       IF (debug) THEN
          PRINT*, 'M_Flux : computing first order source term ...'
       END IF
 
       Source1 = g_r__8/cte__2_d_0_r__8*sg%h**2 - g_r__8/cte__2_d_0_r__8*s1%h**2
       Source2 = g_r__8/cte__2_d_0_r__8*sd%h**2 - g_r__8/cte__2_d_0_r__8*s2%h**2

       IF (debug) THEN
          PRINT*, 'M_Flux : computing higher order centred source term ...'
       END IF

       SELECT CASE(DATA%space_order)

       CASE('godunov')

          SourceCi = cte__0_d_r__8

       CASE('muscl2b', 'muscl4b', 'muscl5')

          hir   =  sg%h
          zir   =  sg%topo
          hil   =  S_R(i-1)%h
          hzil  =  S_R(i-1)%xi
          zil   =  hzil - hil
          SourceCi = g_r__8*(hir + hil)/cte__2_d_0_r__8*(zil - zir)

       CASE('muscl2c', 'eno2mod', 'eno2', 'V4')

          hir   =  S_R(i)%h
          hzir  =  S_R(i)%xi
          zir   =  hzir - hir

          hil   =  S_L(i)%h
          hzil  =  S_L(i)%xi
          zil   =  hzil - hil

          SourceCi = g_r__8*(hir + hil)/cte__2_d_0_r__8*(zil - zir)

       CASE('muscl2a')

          Dxi = MESH%dx

          hir   =  sg%h
          zir   =  sg%topo

          hil   =  Tab(i)%h - Dxi/2.0D0*Grad(i)%h
          !hzil  =  Tab(i)%xi - Dxi/2.0D0*Grad(i)%xi 
          hzil  =  ZF(i) + Tab(i)%h - Dxi/2.0D0*Grad(i)%xi 
          zil   =  hzil - hil

          SourceCi = g_r__8*(hir + hil)/cte__2_d_0_r__8*(zil - zir)


       CASE('muscl4a','muscl4c')

          hir   =  sg%h
          zir   =  sg%topo
          hil   =  Tab(i)%h + Incr_L(i)%h
          !hzil  =  Tab(i)%xi + Incr_L(i)%xi
          hzil  =  ZF(i) + Tab(i)%h + Incr_L(i)%xi
          zil   =  hzil - hil

          SourceCi = g_r__8*(hir + hil)/cte__2_d_0_r__8*(zil - zir)

       CASE('weno3')

          hir   =  sg%h
          zir   =  sg%topo

          hil   =  S_L(i)%h
          hzil  =  S_L(i)%xi
          zil   =  hzil - hil


          !SourceCi = (g_r__8/6.0D0)*(4.0D0*((hil + hic)*(zil - zic) + &
          !     (hic + hir)*(zic - zir))- (hil + hir)*(zil - zir))

          SourceCi = g_r__8*(hir + hil)/2.0D0*(zil - zir)

       CASE('weno5','weno5cv')

          hir   =  sg%h
          zir   =  sg%topo

          hil   =  S_L(i)%h
          hzil  =  S_L(i)%xi
          zil   =  hzil - hil

          hic   = val_C(i)%h
          hzic  = val_C(i)%xi
          zic   = hzic - hic

          SourceCi = (g_r__8/6.0D0)*(4.0D0*((hil + hic)*(zil - zic) + &
               (hic + hir)*(zic - zir))- (hil + hir)*(zil - zir))

       CASE('mixte5')


          hir   =  sg%h
          zir   =  sg%topo

          hil   =  S_R(i-1)%h
          hzil  =  S_R(i-1)%xi
          zil   =  hzil - hil

          hic   = val_C(i)%h
          hzic  = val_C(i)%xi
          zic   = hzic - hic


          SourceCi = (g_r__8/6.0D0)*(4.0D0*((hil + hic)*(zil - zic) + &
               (hic + hir)*(zic - zir))- (hil + hir)*(zil - zir))

       CASE default
          IF (debug) THEN
             PRINT*,'M_Flux : bad order choice (i_order) ...'
             STOP
          END IF
       END SELECT


       FLU%Sx(i)%h = cte__0_d_r__8
       FLU%Sx(i)%qx = SourceCi

       IF (debug) THEN
          PRINT*,'M_Flux : computing interface fluxes including upwind source terms ...'
       END IF

       FLU%X1(i)%h  = FLU_X(i)%h
       FLU%X1(i)%qx = FLU_X(i)%qx + Source1

       FLU%X2(i)%h  = FLU_X(i)%h
       FLU%X2(i)%qx = FLU_X(i)%qx+ Source2


    END DO  !-- i = 0, Nbx

    IF (.NOT. DATA%ChoixDt) THEN
       UMAX  = MAXVAL( u_en_i )
    END IF


    IF (debug) THEN
       PRINT*,'M_Flux :  end of x-fluxes computing ...'
    END IF


    IF(debug ) THEN
       PRINT*, "Exiting ", sub__name !, iproc
    END IF

  END SUBROUTINE Flux_HR






END MODULE ModFlux



































