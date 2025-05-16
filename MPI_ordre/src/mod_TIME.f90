
MODULE modTimeScheme

  USE modCOMMON
  USE modTOOLS
  USE modBoundary
  USE modFLUX
  USE modOUTPUTS
  

  PUBLIC :: time_loop

CONTAINS


  SUBROUTINE time_loop(V, FLU, MESH, DATA, Run,Me,Np)

    IMPLICIT NONE

    TYPE(Variable), INTENT(INOUT) :: V
    TYPE(DATA_IN), INTENT(IN) :: DATA
    TYPE(MeshDef), INTENT(IN) :: MESH
    TYPE(MyRun), INTENT(INOUT) :: Run
    TYPE(VFLU), INTENT(INOUT) :: FLU
    REAL(KIND=8), DIMENSION(1:DATA%Ndf) :: tH
    REAL(KIND=8), DIMENSION(1:DATA%Ndf) :: tQ
    REAL(KIND=8), DIMENSION(1:DATA%Ndf) :: tU
    REAL(KIND=8), DIMENSION(1:DATA%Ndf) :: Eta
    TYPE(STATE), DIMENSION(-2:DATA%Nbx+3):: tmp
    TYPE(STATE), DIMENSION(-2:DATA%Nbx+3, DATA%time_order):: Kvf
    REAL(KIND=8), DIMENSION(-2:DATA%Nbx+3):: sum1
    REAL(KIND=8), DIMENSION(-2:DATA%Nbx+3):: sum2
    REAL(KIND=8), DIMENSION(1:DATA%Ndf):: sum3

    INTEGER :: Nbx
    INTEGER :: ii
    INTEGER :: jj
    !INTEGER :: i,j
    !INTEGER :: Ndf
    INTEGER :: in
    INTEGER :: out
    REAL(KIND=8) :: dx, dt
    !REAL(KIND=8),DIMENSION(1:2, 1:2) :: coeffs_a
    INTEGER :: to
    Integer,intent(in) :: Me , Np

    !WRITE(6,*)'entering time_loop'

    Nbx = DATA%nbx
    in = DATA%in
    out = DATA%out

    dx = MESH%dx
    dt = RUN%dt
    to = DATA%time_order
    !WRITE (UNIT = 6, FMT = * )'DT = ',Run%dt


    Run%step = 1
    Run%subDt = dt
    Run%subTime = Run%Time + Run%subDt
    !print*,Run%subTime

    Kvf(-2:Nbx+3,:)%h =  0.0D0
    Kvf(-2:Nbx+3,:)%qx = 0.0D0

    tmp(-2:Nbx+3)%h = 0.0D0
    tmp(-2:Nbx+3)%qx = 0.0D0

 

    DO ii =1, to

       sum1 = 0.0D0
       sum2 = 0.0D0

       DO jj = 1, to
          sum1 = sum1 + Kvf(-2:Nbx+3, jj)%h*RK(to)%a(ii, jj)
          sum2 = sum2 + Kvf(-2:Nbx+3, jj)%qx*RK(to)%a(ii, jj)
       END DO

       tmp(-2:Nbx+3)%h = V%VF(-2:Nbx+3)%h - (Run%Dt/2.0D0)*sum1 

       tmp(-2:Nbx+3)%qx = V%VF(-2:Nbx+3)%qx - (Run%Dt/2.0D0)*sum2 

       !WRITE(6,*)'1 stage of RK  scheme - VF1'

       CALL Boundary_VF(tmp, V%Ex, V%topoVF, MESH, DATA, Run, ii,Me,Np)

       CALL Flux_HR(tmp, MESH, V%topoVF, DATA, FLU, Run%UMax)

       IF((.NOT.DATA%ChoixDt).AND.(ii==1)) THEN
          CALL TimeStep( DATA, Run, MESH)
       END IF

       Kvf(1:Nbx,ii)%h = (1.0D0/dx)*(FLU%X1(1:Nbx)%h - FLU%X2(0:Nbx-1)%h) &
            
            - (1.0D0/dx)*FLU%Sx(1:Nbx)%h

       Kvf(1:Nbx,ii)%qx = (1.0D0/dx)*(FLU%X1(1:Nbx)%qx - FLU%X2(0:Nbx-1)%qx) &
            
            - (1.0D0/dx)*FLU%Sx(1:Nbx)%qx

    END DO

    !Etape finale schéma RK

    sum1 = 0.0D0
    sum2 = 0.0D0

    DO jj = 1, to
       sum1 = sum1 + Kvf(-2:Nbx+3, jj)%h*RK(to)%b(jj)
       sum2 = sum2 + Kvf(-2:Nbx+3, jj)%qx*RK(to)%b(jj)
    END DO

    V%Vf(-2:Nbx+3)%h = V%Vf(-2:Nbx+3)%h -Run%dt/2.0D0*sum1

    V%Vf(-2:Nbx+3)%qx = V%Vf(-2:Nbx+3)%qx -Run%dt/2.0D0*sum2


    V%VF%xi =  V%VF%h + V%topoVF



  END SUBROUTINE time_loop

END MODULE modTimeScheme
