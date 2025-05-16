!###################################################################
!            COMPUTATION OF VARIOUS ANALYTICAL SOLUTIONS
!###################################################################

MODULE modANALYTIC

USE modTOOLS
USE modCST

PUBLIC :: SOLEX_Steady_Friction
PUBLIC :: ANALYTIC_Init

CONTAINS

  SUBROUTINE ANALYTIC_SOL(xx, RES)

    IMPLICIT NONE

    REAL(kIND = 8), INTENT(in) :: xx
    REAL(kIND = 8), DIMENSION(0:3), INTENT(out)  :: RES
    REAL(kIND = 8) :: c0
    REAL(kIND = 8) :: sigmaref
    REAL(kIND = 8) :: XM
    REAL(kIND = 8) :: m0
    REAL(kIND = 8) :: h0
    REAL(kIND = 8) :: u0
    REAL(kIND = 8) :: href
    REAL(KIND=8) :: hex
    !REAL(kIND = 8) :: chex
    REAL(kIND = 8) :: zz
    !INTEGER :: ii


    SELECT CASE(DATA%NCI)

    CASE(cte__26)
       ! Accuracy test, subcritical flow over a bump
       c0 = 1.0D0
       sigmaref = 2.0D0
       XM = 15.0D0
       m0 = 4.42D0
       h0 = 2.0D0
       u0 = m0/h0

       href=2.0D0
       zz = c0/(sigmaref*SQRT(2.0D0*Pi_r__8))*&
            EXP(-(0.5D0*(xx-XM)/sigmaref)**2)

       hex = (0.5D0*u0**2 + g_r__8*h0 -g_r__8*zz - &
            0.5D0*m0**2/href**2)/g_r__8

       DO WHILE((ABS(hex-href) >= 1.0D-17))
          href = hex
          hex = (0.5D0*u0**2 + g_r__8*h0 -g_r__8*zz &
               - 0.5D0*m0**2/href**2)/g_r__8
       END DO


       RES(0) = hex
       RES(1) = m0
       RES(2) = 0.0D0
       RES(3) = zz

    END SELECT

  END SUBROUTINE ANALYTIC_SOL

  SUBROUTINE ANALYTIC_Init(xx, RES)

    REAL(kIND = 8), INTENT(in) :: xx
    REAL(kIND = 8), DIMENSION(0:3), INTENT(out)  :: RES
    REAL(kIND = 8) :: c0
    REAL(kIND = 8) :: sigmaref
    REAL(kIND = 8) :: XM

    SELECT CASE(DATA%NCI)

    CASE(cte__16)
       ! Accuracy test Chi-Wang Shu

       RES(0) = 5.0D0 + EXP(COS(2.0D0*Pi_r__8*xx))
       RES(1) = SIN(COS(2.0D0*Pi_r__8*xx))
       RES(2) = 0.0D0
       !RES(3) = 2.0D0-(SIN(Pi_r__8*xx))**2
      RES(3) = (SIN(Pi_r__8*xx))**2

    CASE(cte__26)
       ! Accuracy test, subcritical flow over a bump
       c0 = 1.0D0
       sigmaref = 2.0D0
       XM = 15.0D0

       RES(3) = c0/(sigmaref*SQRT(2.0D0*Pi_r__8))*&
            EXP(-(0.5D0*(xx-XM)/sigmaref)**2)
       RES(0) = 2.0D0 - c0/(sigmaref*SQRT(2.0D0*Pi_r__8))*&
            EXP(-(0.5D0*(xx-XM)/sigmaref)**2)
       RES(1) = 0.0D0
       RES(2) = 0.0D0


    END SELECT
  END SUBROUTINE ANALYTIC_Init

  SUBROUTINE SOLEX_Steady_Source(VEX, MESH, DATA, ZEX)

    TYPE(STATE), DIMENSION(:), POINTER :: VEX
    TYPE(MeshDef), INTENT(IN) :: MESH
    TYPE(DATA_IN), INTENT(IN) :: DATA
    !INTEGER, INTENT(IN) :: Nbxref
    REAL(KIND=8), DIMENSION(-2:DATA%Nbx+3) :: DHEX
    REAL(KIND=8), DIMENSION(-2:DATA%Nbx+3) :: DDHEX
    REAL(KIND=8), DIMENSION(:), POINTER :: ZEX
    REAL(KIND=8) :: XX
    REAL(KIND=8) :: Dx
    REAL(KIND=8) :: varz
    REAL(KIND=8) :: discharge
    INTEGER :: i

    INTEGER :: Nbx


    Nbx = DATA%Nbx
    !Nbx = Nbxref



    SELECT CASE(DATA%NCI)

    CASE(cte__44)

       ! sub-critical/sub-critical with diffusion and friction
       ! Initial profile is specific for sub-critical/sub-critical configuration
       ! Diffusion is always switched "on"
       ! Friction coeffs (kl or kt) must be set to 0 in Init.txt
             !to switch off friction effects

!       discharge = 1.5D0
       discharge = 1.5D0
  

          VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
               *(1.0d0 + EXP(-16.d0*(1.d0/2.d0)**2.d0)/2.d0)

          dhex(0)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(-1.d0/2.d0) &
               *EXP(-16.d0*(1.d0/2.d0)**2.d0)/125.d0

          zex(0) = 0.d0
          VEX(0)%qx = discharge
  

       DX = MESH%DX
     
          DO i = 1, Nbx+3

             XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

             VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
                  *(1.d0+EXP(-16.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/2.d0)


             dhex(i)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(XX/1000.d0 -1.0D0/2.0D0) &
                  *EXP(-16.d0*(XX/1000.d0-1.d0/2.d0)**2.d0)/125.d0

             ddhex(i)=-(4.d0/g_r__8)**(1.d0/3.d0) &
                  *(1.d0-32.d0*(XX/1000.d0 - 1.d0/2.d0) &
                  *(XX/1000.0d0 - 1.d0/2.d0)) &
                  *EXP(-16.0D0*(XX/1000.0D0 - 1.d0/2.d0)*(XX/1000.0D0 - 1.d0/2.d0)) &
                  /62500.d0

             VEX(i)%qx = discharge
             
        
             CALL topo_kl_kt_diff(VEX(i)%qx, VEX(i)%h, dhex(i), ddhex(i), varz)
             
             zex(i) = -DX*varz+zex(i-1)
             !print*,'zex(i)',zex(i)
          END DO
    

    CASE(cte__45)

       ! friction + diffusion law in super-critical/super-critical
       ! Initial profile is specific for sub-critical/sub-critical configuration
       ! Diffusion is always switched "on"
       ! Friction coeffs (kl or kt) must be set to 0 in Init.txt
             !to switch off friction effects

   

          VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
               *(1.0d0 - EXP(-36.d0*(1.d0/2.d0)**2.d0)/5.d0)

          dhex(0)=((4.d0/g_r__8)**(1.0d0/3.0d0))*9.0D0 &
               *EXP(-36.0d0*(1.d0/2.d0)**2.d0)/(625.d0*(-2.0d0))

          zex(0) = 0.d0
          VEX(0)%qx = 2.5D0

   

       Dx = MESH%dx

     
          DO i = 1, Nbx

             XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

             VEX(i)%h = ((4.d0/g_r__8)**(1.d0/3.d0)) &
                  *(1.d0-EXP(-36.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/5.d0)

             dhex(i) = ((4.d0/g_r__8)**(1.d0/3.d0))*9.0D0 &
                  *EXP(-36.d0*(XX/1000.d0 - 1.d0/2.d0)**2.d0) &
                  *(XX/1000.d0 -1.0D0/2.0D0)/625.0D0 

             ddhex(i) = (4.d0/g_r__8)**(1.d0/3.d0) &
                  *9.d0*(1.d0-72.d0*(XX/1000.d0-1.d0/2.d0) &
                  *(XX/1000.d0-1.d0/2.d0)) &
                  *EXP(-36.0D0*(XX/1000-1.d0/2.d0)*(XX/1000-1.d0/2.d0)) &
                  /625000.d0

             VEX(i)%qx = 2.5D0

             CALL topo_kl_kt_diff(VEX(i)%qx, VEX(i)%h, dhex(i), ddhex(i), varz)

             zex(i) = -DX*varz+zex(i-1)
          END DO
    

!     CASE(cte__47)

!        ! laminar+turbulent laws, subcritical/subcritical with diffusion

!        DO j = -2, Nby+3

!           VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
!                *(1.0d0 + EXP(-16.d0*(-1.d0/2.d0)**2.d0)/2.d0)

!           dhex(0)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(-1.d0/2.d0) &
!                *EXP(-16.d0*(-1.d0/2.d0)**2.d0)/125.d0

!           zex(0) = 0.d0
!           VEX(0)%qx = 2.d0

!        END DO

!        Dx = MESH%dx

!        DO j = -2,Nby+3
!           DO i = 1, Nbx

!              XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

!              VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
!                   *(1.d0+EXP(-16.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/2.d0)


!              dhex(i)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(XX/1000.d0 -1.0D0/2.0D0) &
!                   *EXP(-16.d0*(XX/1000.d0-1.d0/2.d0)**2.d0)/125.d0

!              ddhex(i)=-(4.d0/g_r__8)**(1.d0/3.d0) &
!                   *(1.d0-32.d0*(XX/1000.d0 - 1.d0/2.d0) &
!                   *(XX/1000.0d0 - 1.d0/2.d0)) &
!                   *exp(-16.0D0*(XX/1000.0D0 - 1.d0/2.d0)*(XX/1000.0D0 - 1.d0/2.d0)) &
!                   /62500.d0

!              VEX(i)%qx = 2.0D0

!              CALL topo_kl_Kt_diff(VEX(i)%qx, VEX(i)%h, dhex(i), ddhex(i), varz)

!              zex(i) = -DX*varz+zex(i-1)
!           END DO
!        END DO

!     CASE(cte__48)

!        ! laminar+turbulent laws + diffusion law in supercritical/supercritical regime

!        DO j = -2, Nby+3

!           VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
!                *(1.0d0 - EXP(-36.d0*(-1.d0/2.d0)**2.d0)/5.d0)

!           dhex(0)=((4.d0/g_r__8)**(1.0d0/3.0d0))*9.0D0 &
!                *EXP(-36.0d0*(-1.d0/2.d0)**2.d0)/(625.d0*(-2.0d0))

!           zex(0) = 0.d0
!           VEX(0)%qx = 2.d0

!        END DO

!        Dx = MESH%dx

!        DO j = -2,Nby+3
!           DO i = 1, Nbx

!              XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

!              VEX(i)%h = ((4.d0/g_r__8)**(1.d0/3.d0)) &
!                   *(1.d0-EXP(-36.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/5.d0)

!              dhex(i) = ((4.d0/g_r__8)**(1.d0/3.d0))*9.0D0 &
!                   *EXP(-36.d0*(XX/1000.d0 - 1.d0/2.d0)**2.d0) &
!                   *(XX/1000.d0 -1.0D0/2.0D0)/625.0D0 

!              ddhex(i) = (4.d0/g_r__8)**(1.d0/3.d0) &
!                   *9.d0*(1.d0-72.d0*(XX/1000.d0-1.d0/2.d0) &
!                   *(XX/1000.d0-1.d0/2.d0)) &
!                   *exp(-36.0D0*(XX/1000-1.d0/2.d0)*(XX/1000-1.d0/2.d0)) &
!                   /625000.d0

!              VEX(i)%qx = 2.0D0

!              CALL topo_kl_kt_diff(VEX(i)%qx, VEX(i)%h, dhex(i), ddhex(i), varz)

!              zex(i) = -DX*varz+zex(i-1)
!           END DO
!        END DO

!     CASE(cte__49)

!        ! turbulent laws, subcritical/subcritical with diffusion

!        DO j = -2, Nby+3

!           VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
!                *(1.0d0 + EXP(-16.d0*(-1.d0/2.d0)**2.d0)/2.d0)

!           dhex(0)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(-1.d0/2.d0) &
!                *EXP(-16.d0*(-1.d0/2.d0)**2.d0)/125.d0

!           zex(0) = 0.d0
!           VEX(0)%qx = 2.d0

!        END DO

!        Dx = MESH%dx

!        DO j = -2,Nby+3
!           DO i = 1, Nbx

!              XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

!              VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
!                   *(1.d0+EXP(-16.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/2.d0)


!              dhex(i)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(XX/1000.d0 -1.0D0/2.0D0) &
!                   *EXP(-16.d0*(XX/1000.d0-1.d0/2.d0)**2.d0)/125.d0

!              ddhex(i)=-(4.d0/g_r__8)**(1.d0/3.d0) &
!                   *(1.d0-32.d0*(XX/1000.d0 - 1.d0/2.d0) &
!                   *(XX/1000.0d0 - 1.d0/2.d0)) &
!                   *exp(-16.0D0*(XX/1000.0D0 - 1.d0/2.d0)*(XX/1000.0D0 - 1.d0/2.d0)) &
!                   /62500.d0

!              VEX(i)%qx = 2.0D0

!              CALL topo_Kt_diff(VEX(i)%qx, VEX(i)%h, dhex(i), ddhex(i), varz)

!              zex(i) = -DX*varz+zex(i-1)
!           END DO
!        END DO

    CASE(cte__47)

       ! friction law from Marche (EJM/B-Fluid 2006) + Manning-Chezy + Darcy law, sub-critical/sub-critical
       ! Initial profile is specific for sub-critical/sub-critical configuration
       ! Diffusion is always switched "off"
       ! Friction coeffs (kl or kt, k_manning, k_darcy) must be set to 0 in Init.txt
             !to switch off particular friction law

          VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
               *(1.0d0 + EXP(-16.d0*(1.d0/2.d0)**2.d0)/2.d0)

          dhex(0)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(-1.d0/2.d0) &
               *EXP(-16.d0*(1.d0/2.d0)**2.d0)/125.d0

          zex(0) = 0.d0
          VEX(0)%qx = 2.d0


       Dx = MESH%dx

          DO i = 1, Nbx+3

             XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

             VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
                  *(1.d0+EXP(-16.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/2.d0)


             dhex(i)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(XX/1000.d0 -1.0D0/2.0D0) &
                  *EXP(-16.d0*(XX/1000.d0-1.d0/2.d0)**2.d0)/125.d0

             VEX(i)%qx = 2.0D0

             CALL topo_friction(VEX(i)%qx, VEX(i)%h, dhex(i), varz)

             zex(i) = -DX*varz+zex(i-1)
          END DO

!     CASE(cte__62)
!        ! turbulent law, subcritical/subcritical

!        DO j = -2, Nby+3

!           VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
!                *(1.0d0 + EXP(-16.d0*(-1.d0/2.d0)**2.d0)/2.d0)

!           dhex(0)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(-1.d0/2.d0) &
!                *EXP(-16.d0*(-1.d0/2.d0)**2.d0)/125.d0

!           zex(0) = 0.d0
!           VEX(0)%qx = 2.d0

!        END DO

!        Dx = MESH%dx

!        DO j = -2,Nby+3
!           DO i = 1, Nbx

!              XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

!              VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
!                   *(1.d0+EXP(-16.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/2.d0)


!              dhex(i)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(XX/1000.d0 -1.0D0/2.0D0) &
!                   *EXP(-16.d0*(XX/1000.d0-1.d0/2.d0)**2.d0)/125.d0

!              VEX(i)%qx = 2.0D0

!              CALL topo_kt(VEX(i)%qx, VEX(i)%h, dhex(i), varz)

!              zex(i) = -DX*varz+zex(i-1)
!           END DO
!        END DO


!     CASE(cte__63)
!        ! laminar & turbulent laws, subcritical/subcritical

!        DO j = -2, Nby+3

!           VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
!                *(1.0d0 + EXP(-16.d0*(-1.d0/2.d0)**2.d0)/2.d0)

!           dhex(0)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(-1.d0/2.d0) &
!                *EXP(-16.d0*(-1.d0/2.d0)**2.d0)/125.d0

!           zex(0) = 0.d0
!           VEX(0)%qx = 2.d0

!        END DO

!        Dx = MESH%dx

!        DO j = -2,Nby+3
!           DO i = 1, Nbx

!              XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

!              VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
!                   *(1.d0+EXP(-16.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/2.d0)


!              dhex(i)=-2.0D0*((4.d0/g_r__8)**(1.d0/3.d0))*(XX/1000.d0 -1.0D0/2.0D0) &
!                   *EXP(-16.d0*(XX/1000.d0-1.d0/2.d0)**2.d0)/125.d0

!              VEX(i)%qx = 2.0D0

!              CALL topo_kl_kt(VEX(i)%qx, VEX(i)%h, dhex(i), varz)

!              zex(i) = -DX*varz+zex(i-1)
!           END DO
!        END DO

    CASE(cte__48)

       ! friction laws in super-critical/super-critical
       ! Initial profile is specific for super-critical/super-critical configuration
       ! Diffusion is always switched "off"
       ! Friction coeffs (kl or kt, k_manning, k_darcy) must be set to 0 in Init.txt
             !to switch off particular friction law


          VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
               *(1.0d0 - EXP(-36.d0*(1.d0/2.d0)**2.d0)/5.d0)

          dhex(0)=((4.d0/g_r__8)**(1.0d0/3.0d0))*9.0D0 &
               *EXP(-36.0d0*(1.d0/2.d0)**2.d0)/(625.d0*(-2.0d0))

          zex(0) = 0.d0
          VEX(0)%qx = 2.d0


       Dx = MESH%dx

          DO i = 1, Nbx

             XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

             VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
                  *(1.d0-EXP(-36.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/5.d0)


             dhex(i)=((4.d0/g_r__8)**(1.d0/3.d0))*9.0D0 &
                  *EXP(-36.d0*(XX/1000.d0 - 1.d0/2.d0)**2.d0) &
                  *(XX/1000.d0 -1.0D0/2.0D0)/625.0D0 

             VEX(i)%qx = 2.0D0

             CALL topo_friction(VEX(i)%qx, VEX(i)%h, dhex(i), varz)

             zex(i) = -DX*varz+zex(i-1)
          END DO
 
!     CASE(cte__65)

!        ! quadratic law in supercritical/supercritical

!        DO j = -2, Nby+3

!           VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
!                *(1.0d0 - EXP(-36.d0*(-1.d0/2.d0)**2.d0)/5.d0)

!           dhex(0)=((4.d0/g_r__8)**(1.0d0/3.0d0))*9.0D0 &
!                *EXP(-36.0d0*(-1.d0/2.d0)**2.d0)/(625.d0*(-2.0d0))

!           zex(0) = 0.d0
!           VEX(0)%qx = 2.d0

!        END DO

!        Dx = MESH%dx

!        DO j = -2,Nby+3
!           DO i = 1, Nbx

!              XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

!              VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
!                   *(1.d0-EXP(-36.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/5.d0)


!              dhex(i)=((4.d0/g_r__8)**(1.d0/3.d0))*9.0D0 &
!                   *EXP(-36.d0*(XX/1000.d0 - 1.d0/2.d0)**2.d0) &
!                   *(XX/1000.d0 -1.0D0/2.0D0)/625.0D0 

!              VEX(i)%qx = 2.0D0

!              CALL topo_kt(VEX(i)%qx, VEX(i)%h, dhex(i), varz)

!              zex(i) = -DX*varz+zex(i-1)
!           END DO
!        END DO

!     CASE(cte__66)

!        ! laminar & quadratic law in supercritical/supercritical

!        DO j = -2, Nby+3

!           VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
!                *(1.0d0 - EXP(-36.d0*(-1.d0/2.d0)**2.d0)/5.d0)

!           dhex(0)=((4.d0/g_r__8)**(1.0d0/3.0d0))*9.0D0 &
!                *EXP(-36.0d0*(-1.d0/2.d0)**2.d0)/(625.d0*(-2.0d0))

!           zex(0) = 0.d0
!           VEX(0)%qx = 2.d0

!        END DO

!        Dx = MESH%dx

!        DO j = -2,Nby+3
!           DO i = 1, Nbx

!              XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

!              VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
!                   *(1.d0-EXP(-36.d0*(XX/1000.d0 - 1.d0 / 2.d0)**2.d0)/5.d0)


!              dhex(i)=((4.d0/g_r__8)**(1.d0/3.d0))*9.0D0 &
!                   *EXP(-36.d0*(XX/1000.d0 - 1.d0/2.d0)**2.d0) &
!                   *(XX/1000.d0 -1.0D0/2.0D0)/625.0D0 

!              VEX(i)%qx = 2.0D0

!              CALL topo_kl_kt(VEX(i)%qx, VEX(i)%h, dhex(i), varz)

!              zex(i) = -DX*varz+zex(i-1)
!           END DO
!        END DO

    CASE(cte__49)

       ! Friction law in sub-critical/super-critical
       ! Initial profile is specific for sub-critical/super-critical configuration
       ! Diffusion is always switched "off"
       ! Friction coeffs (kl or kt, k_manning, k_darcy) must be set to 0 in Init.txt
             !to switch off particular friction law

          VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
               *(1.0d0 - TANH(-3.d0/2.d0)/3.d0)

          dhex(0)=(-(4.d0/g_r__8)**(1.0d0/3.0d0)) &
               /(1000.d0*(COSH(-3.d0/2.d0))**2.d0)

          zex(0) = 0.d0
          VEX(0)%qx = 2.d0


       Dx = MESH%dx

          DO i = 1, Nbx

             XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

             IF(XX<=500.0D0) THEN

                VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
                     *(1.d0-(TANH(3.d0*(XX/1000.d0 -1.d0/2.d0)))/3.d0)


                dhex(i)=(-(4.d0/g_r__8)**(1.d0/3.d0)) &
                     /(1000.d0*(COSH(3.d0*(XX/1000.d0 - 1.d0/2.d0)))**2.d0)

             ELSE

                VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
                     *(1.d0-TANH(6.d0*(XX/1000.d0 -1.d0/2.d0))/6.d0)


                dhex(i)=(-(4.d0/g_r__8)**(1.d0/3.d0)) &
                     /(1000.d0*(COSH(6.d0*(XX/1000.d0 -1.d0/2.d0)))**2.d0)
             END IF

             VEX(i)%qx = 2.0D0

             CALL topo_friction(VEX(i)%qx, VEX(i)%h, dhex(i), varz)

             zex(i) = -DX*varz+zex(i-1)
          END DO

!     CASE(cte__68)

!        ! turbulent law in sub/sup

!        DO j = -2, Nby+3

!           VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
!                *(1.0d0 - TANH(-3.d0/2.d0)/3.d0)

!           dhex(0)=(-(4.d0/g_r__8)**(1.0d0/3.0d0)) &
!                /(1000.d0*(COSH(-3.d0/2.d0))**2.d0)

!           zex(0) = 0.d0
!           VEX(0)%qx = 2.d0

!        END DO

!        Dx = MESH%dx

!        DO j = -2,Nby+3
!           DO i = 1, Nbx

!              XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

!              IF(XX<=500.0D0) THEN

!                 VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
!                      *(1.d0-(TANH(3.d0*(XX/1000.d0-1.d0/2.d0)))/3.d0)


!                 dhex(i)=(-(4.d0/g_r__8)**(1.d0/3.d0)) &
!                      /(1000.d0*(COSH(3.d0*(XX/1000.d0 - 1.d0/2.d0)))**2.d0)

!              ELSE

!                 VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
!                      *(1.d0-TANH(6.d0*(XX/1000.d0-1.d0/2.d0))/6.d0)


!                 dhex(i)=(-(4.d0/g_r__8)**(1.d0/3.d0)) &
!                      /(1000.d0*(COSH(6.d0*(XX/1000.d0-1.d0/2.d0)))**2.d0)
!              END IF

!              VEX(i)%qx = 2.0D0

!              CALL topo_kt(VEX(i)%qx, VEX(i)%h, dhex(i), varz)

!              zex(i) = -DX*varz+zex(i-1)
!           END DO
!        END DO

    CASE(cte__46)

       ! Friction law in super-critical/sub-critical
       ! Initial profile is specific for super-critical/sub-critical configuration
       ! Diffusion is always switched "off"
       ! Friction coeffs (kl or kt, k_manning, k_darcy) must be set to 0 in Init.txt
             !to switch off particular friction law

          VEX(0)%h = ((4.0D0/g_r__8)**(1.0D0/3.d0)) &
               *(9.d0/10.d0-1.d0/6.d0)

          dhex(0)=((4.d0/g_r__8)**(1.0d0/3.0d0))/1500.d0 

          zex(0) = 0.d0
          VEX(0)%qx = 2.d0

       Dx = MESH%dx

          DO i = 1, Nbx

             XX = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

             IF(XX<=500.0D0) THEN

                VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
                     *(9.d0/10.d0-EXP(-XX/250.d0)/6.d0)

                dhex(i)=((4.d0/g_r__8)**(1.d0/3.d0)) &
                     *EXP(-XX/250.d0)/1500.d0
             ELSE
                VEX(i)%h=((4.d0/g_r__8)**(1.d0/3.d0)) &
                     *(1.d0 &
                     - 0.348427D0*EXP(-20.d0*(XX/1000.d0-1.d0/2.d0)) &
                     + 0.552264D0*EXP(-40.d0*(XX/1000.d0-1.d0/2.d0)) &
                     - 0.55558D0*EXP(-60.d0*(XX/1000.d0-1.d0/2.d0)) &
                     + 4.d0*EXP(XX/1000.d0-1.d0)/5.d0)


                dhex(i)=((4.d0/g_r__8)**(1.d0/3.d0)) &
                     *((-0.348427*EXP(-20.d0*(XX/1000.d0-1.d0/2.d0)) &
                     +2.d0*0.552264*EXP(-40.d0*(XX/1000.d0-1.d0/2.d0)) &
                     -3.d0*0.55558*EXP(-60.d0*(XX/1000.d0-1.d0/2.d0)))/(-50.d0) &
                     +EXP(XX/1000.d0-1.d0)/1250.d0)
             END IF

             VEX(i)%qx = 2.0D0

             CALL topo_friction(VEX(i)%qx, VEX(i)%h, dhex(i), varz)

             zex(i) = -DX*varz+zex(i-1)
          END DO

    END SELECT

  END SUBROUTINE SOLEX_Steady_Source
 
  SUBROUTINE topo_friction(q, h, dh, dz)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(in) :: q
    REAL(KIND=8), INTENT(in) :: h
    REAL(KIND=8), INTENT(in) :: dh
    REAL(KIND=8), INTENT(out) :: dz

    dz = (1.d0 - q**2.d0 / (g_r__8*h**3.d0))*dh &
     + DATA%kl*q/(g_r__8*h**2.d0*(1.d0+DATA%kl*h/(3.d0*DATA%mu_v))) &  ! laminar law
     + DATA%kt*q**2.0D0/(g_r__8*h**2.d0*(1.d0+DATA%kl*h/(3.d0*DATA%mu_v))**2.0D0) &  ! turbulent law
     + ((DATA%k_manning*q)**2.d0)/(h**(10.0d0/3.0d0)) & ! Manning law
     + DATA%k_darcy*q**2.d0/(8.d0*g_r__8*h**3.d0) !darcy-Weisbar law

    !print*,'dz',dz
  END SUBROUTINE topo_friction

!   SUBROUTINE topo_kl(q, h, dh, dz)

!     IMPLICIT NONE

!     REAL(KIND=8), INTENT(in) :: q
!     REAL(KIND=8), INTENT(in) :: h
!     REAL(KIND=8), INTENT(in) :: dh

!     REAL(KIND=8), INTENT(out) :: dz
!     !REAL(KIND=8), INTENT(out) :: kl
!     !REAL(KIND=8), INTENT(out) :: mu

!     dz = (1.d0 - q**2.d0 / (g_r__8*h**3.d0))*dh &
!      + DATA%kl*q/(g_r__8*h**2.d0*(1.d0+DATA%kl*h/(3.d0*DATA%mu)))


!   END SUBROUTINE topo_kl

!   SUBROUTINE topo_kt(q, h, dh, dz)

!     IMPLICIT NONE

!     REAL(KIND=8), INTENT(in) :: q
!     REAL(KIND=8), INTENT(in) :: h
!     REAL(KIND=8), INTENT(in) :: dh

!     REAL(KIND=8), INTENT(out) :: dz
!     !REAL(KIND=8), INTENT(out) :: kl
!     !REAL(KIND=8), INTENT(out) :: mu

!     dz = (1.d0 - q**2.d0 / (g_r__8*h**3.d0))*dh &
!      + DATA%kt*q**2.0D0/(g_r__8*h**2.d0*(1.d0+DATA%kl*h/(3.d0*DATA%mu))**2.0D0)
      
!   END SUBROUTINE topo_kt

!   SUBROUTINE topo_kl_kt(q, h, dh, dz)

!     IMPLICIT NONE

!     REAL(KIND=8), INTENT(in) :: q
!     REAL(KIND=8), INTENT(in) :: h
!     REAL(KIND=8), INTENT(in) :: dh

!     REAL(KIND=8), INTENT(out) :: dz
!     !REAL(KIND=8), INTENT(out) :: kl
!     !REAL(KIND=8), INTENT(out) :: mu

!     dz = (1.d0 - q**2.d0 / (g_r__8*h**3.d0))*dh &
!      + DATA%kt*q**2.0D0/(g_r__8*h**2.d0*(1.d0+DATA%kl*h/(3.d0*DATA%mu))**2.0D0) &
!      + DATA%kl*q/(g_r__8*h**2.d0*(1.d0+DATA%kl*h/(3.d0*DATA%mu)))

!   END SUBROUTINE topo_kl_kt

!   SUBROUTINE topo_Manning(q, h, dh, R, n, dz)

!     IMPLICIT NONE

!     REAL(KIND=8), INTENT(in) :: q
!     REAL(KIND=8), INTENT(in) :: h
!     REAL(KIND=8), INTENT(in) :: dh
!     REAL(KIND=8), INTENT(in) :: R
!     REAL(KIND=8), INTENT(in) :: n
!     REAL(KIND=8), INTENT(out) :: dz
!     !REAL(KIND=8), INTENT(out) :: kl
!     !REAL(KIND=8), INTENT(out) :: mu

!     dz = (1.d0 - q**2.d0 / (g_r__8*h**3.d0))*dh &
!          + 2.d0*q*R/(g_r__8*(h**2.0d0)) &
!          + ((n*q)**2.d0)/(h**(10.0d0/3.0d0))

!   END SUBROUTINE topo_Manning

!   SUBROUTINE topo_Darcy_Weisbach(q, h, dh, R, n, dz)

!     IMPLICIT NONE

!     REAL(KIND=8), INTENT(in) :: q
!     REAL(KIND=8), INTENT(in) :: h
!     REAL(KIND=8), INTENT(in) :: dh
!     REAL(KIND=8), INTENT(in) :: R
!     REAL(KIND=8), INTENT(in) :: n

!     REAL(KIND=8), INTENT(out) :: dz
!     !REAL(KIND=8), INTENT(out) :: kl
!     !REAL(KIND=8), INTENT(out) :: mu

!     dz = (1.d0 - q**2.d0 / (g_r__8*h**3.d0))*dh &
!          +2.d0*q*R/(g_r__8*h**2.d0) &
!          +n*q**2.d0/(8.d0*g_r__8*h**3.d0)

!   END SUBROUTINE topo_Darcy_Weisbach

!   SUBROUTINE topo_kl_diff(q, h, dh, ddh, dz)

!     IMPLICIT NONE

!     REAL(KIND=8), INTENT(in) :: q
!     REAL(KIND=8), INTENT(in) :: h
!     REAL(KIND=8), INTENT(in) :: dh
!     REAL(KIND=8), INTENT(in) :: ddh
!     REAL(KIND=8), INTENT(out) :: dz

!     dz= (1.d0-q**2.d0/(g_r__8*h**3.d0))*dh &
!          + DATA%kl*q/(g_r__8*h**2.d0*(1.d0+DATA%kl*h/(3.d0*DATA%mu))) &
!          + 4.d0*DATA%mu*(q*ddh+q*dh*dh/h)/(g_r__8*h*h)

!   END SUBROUTINE topo_kl_diff

!   SUBROUTINE topo_kt_diff(q, h, dh, ddh, dz)

!     IMPLICIT NONE

!     REAL(KIND=8), INTENT(in) :: q
!     REAL(KIND=8), INTENT(in) :: h
!     REAL(KIND=8), INTENT(in) :: dh
!     REAL(KIND=8), INTENT(in) :: ddh
!     REAL(KIND=8), INTENT(out) :: dz
!     !REAL(KIND=8), INTENT(out) :: kl
!     !REAL(KIND=8), INTENT(out) :: mu

!     dz = (1.d0 - q**2.d0 / (g_r__8*h**3.d0))*dh &
!      + DATA%kt*q**2.0D0/(g_r__8*h**2.d0*(1.d0+DATA%kl*h/(3.d0*DATA%mu))**2.0D0) &
!     + 4.d0*DATA%mu*(q*ddh+q*dh*dh/h)/(g_r__8*h*h)

!   END SUBROUTINE topo_kt_diff

  SUBROUTINE topo_kl_kt_diff(q, h, dh, ddh, dz)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(in) :: q
    REAL(KIND=8), INTENT(in) :: h
    REAL(KIND=8), INTENT(in) :: dh
    REAL(KIND=8), INTENT(in) :: ddh
    REAL(KIND=8), INTENT(out) :: dz
    !REAL(KIND=8), INTENT(out) :: kl
    !REAL(KIND=8), INTENT(out) :: mu

    dz = (1.d0 - q**2.d0 / (g_r__8*h**3.d0))*dh &
     + DATA%kt*q**2.0D0/(g_r__8*h**2.d0*(1.d0+DATA%kl*h/(3.d0*DATA%mu_v))**2.0D0) & ! laminar friction law
     + DATA%kl*q/(g_r__8*h**2.d0*(1.d0+DATA%kl*h/(3.d0*DATA%mu_v))) & ! quadratic friction law
     + 4.d0*DATA%mu_h*(q*ddh+q*dh*dh/h)/(g_r__8*h*h) & ! diffusion
     + ((DATA%k_manning*q)**2.d0)/(h**(10.0d0/3.0d0)) & ! Manning law
     + DATA%k_darcy*q**2.d0/(8.d0*g_r__8*h**3.d0) !darcy-Weisbar law

    !print*,'ligne diff', 4.d0*DATA%mu_h*(q*ddh+q*dh*dh/h)/(g_r__8n*h*h)
  END SUBROUTINE topo_kl_kt_diff

END MODULE modANALYTIC
