!###################################################################
!#          FORCING BOUNDARY CONDITIONS                            #
!#                                                                 #
!###################################################################


MODULE modBoundary
  USE MPI
  USE modCOMMON
  USE modTOOLS
  USE modCST

  PRIVATE
  PUBLIC :: CL_Subcritique_Koba
  PUBLIC :: CL_Subcritique_2Invariants
  PUBLIC :: sub__Boundary_Cell
  PUBLIC :: F
  PUBLIC :: DEBIT_LEFT
  PUBLIC :: COTE_RIGHT
  PUBLIC :: sub__Boundary_FLU
  PUBLIC :: CL_Subcritique_Bonneton
  PUBLIC :: Sponge_layer
  PUBLIC :: Boundary_VF

CONTAINS

  SUBROUTINE Boundary_VF(Tab, TabEx, topoVF, MESH, DATA, Run, ii,Me,Np)      

    IMPLICIT NONE
    integer ::  statinfo
    integer,dimension(MPI_STATUS_SIZE) :: status

    TYPE(MyRun), INTENT(IN) :: Run
    TYPE(DATA_IN), INTENT(IN) :: DATA
    TYPE( STATE ), DIMENSION(-2:DATA%Nbx+3), INTENT(INOUT) :: Tab
    TYPE( STATE ), DIMENSION(-2:DATA%Nbx+3) :: Tab_tmp1 , Tab_tmp2 , Tab_tmp3
    TYPE( STATE ), DIMENSION(-2:DATA%Nbx+3) :: tmp_Tab1 , tmp_Tab2 , tmp_Tab3
    TYPE( STATE ), DIMENSION(:), POINTER :: TabEx 
    REAL(KIND=8), DIMENSION(:), POINTER :: topoVF 
    TYPE(MeshDef), INTENT(IN) :: MESH
    !REAL(KIND=8), DIMENSION(:), POINTER :: topo
    !REAL(kind = r__8) :: Time
    !INTEGER, INTENT(in) :: Nt
    INTEGER, INTENT(IN) :: ii ! iteration in RK scheme
    INTEGER :: in
    INTEGER :: out
    INTEGER :: Nbx
    REAL(KIND=8) :: dt

    REAL(KIND=8) :: Ksi
    REAL(KIND=8) :: H0
    REAL(KIND=8) :: eta0
    !REAL(KIND=8) :: H
    REAL(KIND=8) :: p
    REAL(KIND=8) :: dx
    REAL(KIND=8) :: u_0
    REAL(KIND=8) :: u_1
    REAL(KIND=8) :: c_0
    REAL(KIND=8) :: c_1
    REAL(KIND=8) :: period
    REAL(KIND=8) :: Time
    REAL(KIND=8) :: u
    INTEGER :: to
    REAL(KIND=8) :: a
    REAL(KIND=8) :: k
    REAL(KIND=8) :: c
    REAL(KIND=8) :: L
    REAL(KIND=8) :: d
    REAL(KIND=8) :: z0

    real*8 :: ha1 , ha2 , tmp_hg , tmp_hd , tmp_ug , tmp_ud
    real*8 :: x_R , x_L
    real*8 :: C1 , C2 , C3  , C4
    real*8,dimension(:),allocatable :: te ,  yy
    integer :: i,N
    Integer,intent(in) :: Me , Np
    real*8 :: u1 , u2 , u3 , h1 , h2 , h3 , q1 , q2 , q3
    real*8 :: uu1 , uu2 , uu3 , hh1 , hh2 , hh3 , qq1 , qq2 , qq3
    

    !print*,"ProId = " , "  ", Me , "  ","ProcNo = ","  ",Np

    Nbx = DATA%Nbx
    in = 1
    out = DATA%Nbx
    !TIME = Run%time
    !Time = Run%Time
    !Dt = Run%subDt
    to = DATA%time_order

    SELECT CASE(DATA%cond_lim)
    CASE( 1)
       ! Periodic
       ! Tab(in-3:in-1)%h = Tab(out-2:out)%h
       ! Tab(in-3:in-1)%qx = Tab(out-2:out)%qx
       ! Tab(out+1:out+3)%h = Tab(in:in+2)%h
       ! Tab(out+1:out+3)%qx = Tab(in:in+2)%qx


       Tab(0) = Tab(Nbx)
       Tab(-1) = Tab(Nbx-1)
       Tab(-2) = Tab(Nbx-2)

       Tab(Nbx+1) = Tab(1)
       Tab(Nbx+2) = Tab(2)
       Tab(Nbx+3) = Tab(3)

       TopoVf(Nbx+1) = TopoVf(1)
       TopoVf(Nbx+2) = TopoVf(2)
       TopoVf(Nbx+3) = TopoVf(3)

    CASE(2)
       ! Neumann homogène / transparent 
       Tab(0) = Tab(1)
       Tab(-1) = Tab(2)
       Tab(-2) = Tab(3)

       TopoVf(0) = TopoVf(1)
       TopoVf(-1) = TopoVf(2)
       TopoVf(-2) = TopoVf(3)

       Tab(Nbx+1) = Tab(Nbx)
       Tab(Nbx+2) = Tab(Nbx-1)
       Tab(Nbx+3) = Tab(Nbx-2)

       TopoVf(Nbx+1) = TopoVf(Nbx)
       TopoVf(Nbx+2) = TopoVf(Nbx-1)
       TopoVf(Nbx+3) = TopoVf(Nbx-2)

    case(1996)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       if(Np /= 1) then 

          if(Me == 0) then
             !CONDITION GAUCHE A METTRE
             Dx = MESH%dx
             dt = Run%subDt*RK(to)%s_dt(ii)
             time = Run%time + dt


             Tab(1)%qx = 0.0d0
             Tab(1)%u = Tab(1)%qx/Tab(1)%h
             !CONDITION GAUCHE


             N = 2

             allocate(te(N) , yy(N))

             do i = 1 , N
                te(i) = MESH%X_DF(i+1)
                yy(i) = (Tab(i+1)%u) - sqrt(g_r__8*Tab(i+1)%h)
             end do

             x_L = zero_fonction(te,yy,N, Data%x_in , Data%x_out, Dt)

             deallocate(te)
             deallocate(yy)

             allocate(te(N) , yy(N))

             do i = 1 , N
                te(i) = MESH%X_DF(i+1)
                yy(i) = (Tab(i+1)%u) - 2.0d0*sqrt(g_r__8*Tab(i+1)%h) 
             end do

             C1 = interpolation(x_L,te,yy,N)

             deallocate(te)
             deallocate(yy)

             !Tab(1)%u = C1 + 2.0d0*sqrt(g*Tab(1)%h) !!!!! a enlever

             Tab(1)%h = newton_droite(Tab(1)%h , Tab(1)%u , C1)
             Tab(0) = Tab(1)
             Tab(-1) = Tab(0)
             Tab(-2) = Tab(-1)
             
             u1 = Tab(Nbx)%u
             h1 = Tab(Nbx)%h
             q1 = Tab(Nbx)%qx

             u2 = Tab(Nbx-1)%u
             h2 = Tab(Nbx-1)%h
             q2 = Tab(Nbx-1)%qx
!!$
!!$             u3 = Tab(Nbx-2)%u
!!$             h3 = Tab(Nbx-2)%h
!!$             q3 = Tab(Nbx-2)%qx

             call MPI_SEND(u1 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 20 , MPI_COMM_WORLD, statinfo)
             call MPI_SEND(h1 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 21 , MPI_COMM_WORLD, statinfo)
             call MPI_SEND(q1 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 22 , MPI_COMM_WORLD, statinfo)


             call MPI_SEND(u2 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 23 , MPI_COMM_WORLD, statinfo)

             call MPI_SEND(h2 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 24 , MPI_COMM_WORLD, statinfo)
             call MPI_SEND(q2 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 25 , MPI_COMM_WORLD, statinfo)

!!$
!!$
!!$             call MPI_SEND(u3 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 26 , MPI_COMM_WORLD, statinfo)
!!$
!!$             call MPI_SEND(h3 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 27 , MPI_COMM_WORLD, statinfo)
!!$             call MPI_SEND(q3 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 28 , MPI_COMM_WORLD, statinfo)
!!$
!!$
             call MPI_RECV(uu1 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 29 , MPI_COMM_WORLD, status,statinfo)
             call MPI_RECV(hh1 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 30 , MPI_COMM_WORLD, status,statinfo)
             call MPI_RECV(qq1 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 31 , MPI_COMM_WORLD, status,statinfo)


             call MPI_RECV(uu2 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 32 , MPI_COMM_WORLD, status,statinfo)

             call MPI_RECV(hh2 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 33 , MPI_COMM_WORLD, status,statinfo)
             call MPI_RECV(qq2 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 34 , MPI_COMM_WORLD, status,statinfo)
!!$
!!$
!!$
!!$             call MPI_RECV(uu3 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 35 , MPI_COMM_WORLD, status,statinfo)
!!$
!!$             call MPI_RECV(hh3 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 36 , MPI_COMM_WORLD, status,statinfo)
!!$             call MPI_RECV(qq3 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 37 , MPI_COMM_WORLD, status,statinfo)
!!$
!!$
             Tab(Nbx+2)%u = uu1
             Tab(Nbx+2)%h = hh1
             Tab(Nbx+2)%qx = qq1

             Tab(Nbx+3)%u = uu2
             Tab(Nbx+3)%h = hh2
             Tab(Nbx+3)%qx = qq2
!!$
!!$             Tab(Nbx+4)%u = uu3
!!$             Tab(Nbx+4)%h = hh3
!!$             Tab(Nbx+4)%qx = qq3

             allocate(te(N) , yy(N))

             do i = 1 , N
                te(i) = MESH%X_DF(Nbx-i+1)
                yy(i) = (Tab(Nbx-i+1)%u) + sqrt(g_r__8*Tab(Nbx-i+1)%h)
             end do

             x_R = zero_fonction(te,yy,N, Data%x_out , Data%x_in, Dt)

             deallocate(te)
             deallocate(yy)


             allocate(te(N) , yy(N))

             do i = 1 , N
                te(i) = MESH%X_DF(Nbx-i+1)
                yy(i) = (Tab(Nbx-i+1)%u) + 2.0d0*sqrt(g_r__8*Tab(Nbx-i+1)%h) 
             end do
             C2 =  interpolation(x_R,te,yy,N)

             deallocate(te)
             deallocate(yy)

             call MPI_SEND(C2 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 10 , MPI_COMM_WORLD, statinfo)
             call MPI_RECV(C1 , 1 , MPI_DOUBLE_PRECISION , Me + 1 , 11 , MPI_COMM_WORLD,status, statinfo)
             Tab(Nbx+1)%u = ((C1 + C2)/(2.0d0))
             Tab(Nbx+1)%qx = Tab(Nbx+1)%h*Tab(Nbx+1)%u
             
             Tab(Nbx+1)%h = (((C2 - C1)**2)/(16.0d0*g_r__8))


          else
                N = 2
               Dx = MESH%dx
               dt = Run%subDt*RK(to)%s_dt(ii)
               time = Run%time + dt
                !CONDITION DROITE A METTRE

                !Tab(Nbx+1)%h = 1.0d0 !1.0d0 + 1.0d-01*(cos((8/(1.0d0*5.0d-02))*atan(1.0d0)*(time))) !1.0d0
               Tab(Nbx+1)%qx = 0.0d0
               Tab(Nbx + 1)%u = Tab(Nbx+1)%qx/Tab(Nbx+1)%h
                allocate(te(N) , yy(N))

                do i = 1 , N
                   te(i) = MESH%X_DF(Nbx-i+1)
                   yy(i) = (Tab(Nbx-i+1)%u) + sqrt(g_r__8*Tab(Nbx-i+1)%h)
                end do

                x_R = zero_fonction(te,yy,N, Data%x_out , Data%x_in, Dt)

                deallocate(te)
                deallocate(yy)

                allocate(te(N) , yy(N))

                do i = 1 , N
                   te(i) = MESH%X_DF(Nbx-i+1)
                   yy(i) = (Tab(Nbx-i+1)%u) + 2.0d0*sqrt(g_r__8*Tab(Nbx-i+1)%h) 
                end do

                Tab(Nbx+1)%h = newton_gauche(Tab(Nbx+1)%h , Tab(Nbx+1)%u , interpolation(x_R,te,yy,N))
                !U(Nbx+1)%zeta = U(Nbx+1)%h + topographie(X_DF(Nbx+1))
                !Tab(Nbx+1)%u = interpolation(x_R,te,yy,N) - 2.0d0*sqrt(g_r__8*Tab(Nbx+1)%h)
                !Tab(Nbx+1)%qx = Tab(Nbx+1)%h*Tab(Nbx+1)%u
                Tab(Nbx+2) = Tab(Nbx+1)
                Tab(Nbx + 3) = Tab(Nbx+2)

                deallocate(te)
                deallocate(yy)


                uu1 = Tab(2)%u
                hh1 = Tab(2)%h
                qq1 = Tab(2)%qx

                uu2 = Tab(3)%u
                hh2 = Tab(3)%h
                qq2 = Tab(3)%qx
!!$
!!$                uu3 =Tab(4)%u
!!$                hh3 = Tab(4)%h
!!$                qq3 = Tab(4)%qx
!!$
                call MPI_SEND(uu1 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 29 , MPI_COMM_WORLD, statinfo)
                call MPI_SEND(hh1 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 30 , MPI_COMM_WORLD, statinfo)
                call MPI_SEND(qq1 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 31 , MPI_COMM_WORLD, statinfo)


                call MPI_SEND(uu2 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 32 , MPI_COMM_WORLD, statinfo)

                call MPI_SEND(hh2 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 33 , MPI_COMM_WORLD, statinfo)
                call MPI_SEND(qq2 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 34 , MPI_COMM_WORLD, statinfo)

!!$
!!$
!!$                call MPI_SEND(uu3 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 35 , MPI_COMM_WORLD, statinfo)
!!$
!!$                call MPI_SEND(hh3 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 36 , MPI_COMM_WORLD, statinfo)
!!$                call MPI_SEND(qq3 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 37 , MPI_COMM_WORLD, statinfo)
!!$

                call MPI_RECV(u1 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 20 , MPI_COMM_WORLD, status,statinfo)
                call MPI_RECV(h1 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 21 , MPI_COMM_WORLD, status,statinfo)
                call MPI_RECV(q1 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 22 , MPI_COMM_WORLD, status,statinfo)


                call MPI_RECV(u2 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 23 , MPI_COMM_WORLD, status,statinfo)

                call MPI_RECV(h2 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 24 , MPI_COMM_WORLD, status,statinfo)
                call MPI_RECV(q2 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 25 , MPI_COMM_WORLD, status,statinfo)
!!$
!!$
!!$
!!$                call MPI_RECV(u3 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 26 , MPI_COMM_WORLD, status,statinfo)
!!$
!!$                call MPI_RECV(h3 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 27 , MPI_COMM_WORLD, status,statinfo)
!!$                call MPI_RECV(q3 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 28 , MPI_COMM_WORLD, status, statinfo)
!!$
                Tab(0)%u = u1
                Tab(0)%h = h1
                Tab(0)%qx = q1

                Tab(-1)%u = u2
                Tab(-1)%h = h2
                Tab(-1)%qx = q2
!!$
!!$                Tab(-2)%u = u3
!!$                Tab(-2)%h = h3
!!$                Tab(-2)%qx = q3

                allocate(te(N) , yy(N))

                do i = 1 , N
                   te(i) = MESH%X_DF(i+1)
                   yy(i) = (Tab(i+1)%u) - sqrt(g_r__8*Tab(i+1)%h)
                end do

                x_L = zero_fonction(te,yy,N, Data%x_in , Data%x_out, Dt)

                deallocate(te)
                deallocate(yy)

                allocate(te(N) , yy(N))

                do i = 1 , N
                   te(i) = MESH%X_DF(i+1)
                   yy(i) = (Tab(i+1)%u) - 2.0d0*sqrt(g_r__8*Tab(i+1)%h) 
                end do

                C1 = interpolation(x_L,te,yy,N)

                deallocate(te)
                deallocate(yy)
                call MPI_SEND(C1 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 11 , MPI_COMM_WORLD, statinfo)
                call MPI_RECV(C2 , 1 , MPI_DOUBLE_PRECISION , Me - 1 , 10 , MPI_COMM_WORLD,status, statinfo)

                Tab(1)%u = ((C1 + C2)/(2.0d0))
                Tab(1)%qx = Tab(1)%h*Tab(1)%u
                Tab(1)%h = (((C2 - C1)**2)/(16.0d0*g_r__8))


             end if
          end if
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



                

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CASE(3)
       !WRITE(6,*) 'Reflexive / solid wall'
       ! Tab(in-3:in-1)%h = Tab(in+2:-1:in)%h
       ! Tab(in-3:in-1)%qx = -Tab(in+2:-1:in)%qx
       ! Tab(out+1:out+3)%h = Tab(out:-1:out-2)%h
       ! Tab(out+1:out+3)%qx = -Tab(out:-1:out-2)%qx      
       Tab(0)%h = Tab(1)%h
       Tab(0)%qx = -Tab(1)%qx

       ! WRITE(6,*)'In Bondary'
       ! WRITE(6,*)'Tab(0)%h',Tab(0)%h
       ! WRITE(6,*)'Tab(1)%h',Tab(1)%h

       Tab(-1)%h = Tab(2)%h
       Tab(-1)%qx = -Tab(2)%qx

       Tab(-2)%h = Tab(3)%h
       Tab(-2)%qx = -Tab(3)%qx

       Tab(Nbx+1)%h = Tab(Nbx)%h
       Tab(Nbx+1)%qx = -Tab(Nbx)%qx

       ! WRITE(6,*)'Tab(Nbx)%h',Tab(Nbx)%h
       ! WRITE(6,*)'Tab(Nbx+1)%h',Tab(Nbx+1)%h
       ! WRITE(6,*)'Tab(Nbx)%qx',Tab(Nbx)%qx
       ! WRITE(6,*)'Tab(Nbx+1)%qx',Tab(Nbx+1)%qx

       Tab(Nbx+2)%h = Tab(Nbx-1)%h
       Tab(Nbx+2)%qx = - Tab(Nbx-1)%qx

       Tab(Nbx+3)%h = Tab(Nbx-2)%h
       Tab(Nbx+3)%qx = -Tab(Nbx-2)%qx

       TopoVf(0) = TopoVf(1)
       TopoVf(-1) = TopoVf(2)
       TopoVf(-2) = TopoVf(3)

       TopoVf(Nbx+1) = TopoVf(Nbx)
       TopoVf(Nbx+2) = TopoVf(Nbx-1)
       TopoVf(Nbx+3) = TopoVf(Nbx-2)

    CASE(4)

       ! Generation sur SV / Neumann sur dispersion

       SELECT CASE(DATA%NCI)
          
       CASE(0)
          !WRITE(6,*) 'Reflexive / solid wall'
          ! Tab(in-3:in-1)%h = Tab(in+2:-1:in)%h
          ! Tab(in-3:in-1)%qx = -Tab(in+2:-1:in)%qx
          ! Tab(out+1:out+3)%h = Tab(out:-1:out-2)%h
          ! Tab(out+1:out+3)%qx = -Tab(out:-1:out-2)%qx      
          Tab(0)%h = Tab(1)%h
          Tab(0)%qx = -Tab(1)%qx
          
          ! WRITE(6,*)'In Bondary'
          ! WRITE(6,*)'Tab(0)%h',Tab(0)%h
          ! WRITE(6,*)'Tab(1)%h',Tab(1)%h
          
          Tab(-1)%h = Tab(2)%h
          Tab(-1)%qx = -Tab(2)%qx

          Tab(-2)%h = Tab(3)%h
          Tab(-2)%qx = -Tab(3)%qx

          Tab(Nbx+1)%h = Tab(Nbx)%h
          Tab(Nbx+1)%qx = -Tab(Nbx)%qx

          ! WRITE(6,*)'Tab(Nbx)%h',Tab(Nbx)%h
          ! WRITE(6,*)'Tab(Nbx+1)%h',Tab(Nbx+1)%h
          ! WRITE(6,*)'Tab(Nbx)%qx',Tab(Nbx)%qx
          ! WRITE(6,*)'Tab(Nbx+1)%qx',Tab(Nbx+1)%qx

          Tab(Nbx+2)%h = Tab(Nbx-1)%h
          Tab(Nbx+2)%qx = - Tab(Nbx-1)%qx

          Tab(Nbx+3)%h = Tab(Nbx-2)%h
          Tab(Nbx+3)%qx = -Tab(Nbx-2)%qx

          TopoVf(0) = TopoVf(1)
          TopoVf(-1) = TopoVf(2)
          TopoVf(-2) = TopoVf(3)

          TopoVf(Nbx+1) = TopoVf(Nbx)
          TopoVf(Nbx+2) = TopoVf(Nbx-1)
          TopoVf(Nbx+3) = TopoVf(Nbx-2)

       CASE(1) ! génération d'un soliton

          Dx = MESH%dx
          dt = Run%subDt*RK(to)%s_dt(ii)
          time = Run%time + dt

          H0 = 0.5D0
          A = 0.2D0*H0
          k = sqrt(3.0D0*a/(4.0D0*H0**2*(H0+A)))
          c = sqrt(g_r__8*(H0+a))
          L = DATA%x_out - DATA%x_in
          d = -7.0D0
          z0 = 0.0D0

          Ksi = a/(cosh(k*(z0- c*Time-d))**2)

          ! u_0 = Tab(0)%qx/Tab(0)%h
          ! c_0 = SQRT(g_r__8*Tab(0)%h)

          ! u_1 = Tab(1)%qx/Tab(1)%h
          ! c_1 = SQRT(g_r__8*Tab(1)%h)

          ! !CALL CL_Subcritique_2Invariants(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u, D)

          ! !CALL CL_Subcritique_Koba(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u, D)
          ! call CL_Subcritique_Bonneton(ksi, H0, u_0, u_1, c_0, c_1, dt, Dx, u)

          ! Definition of the inlet boundary condition


          ! Tab(0)%h = H0 + ksi
          ! Tab(0)%qx = Tab(0)%h*u           
          ! Tab(-1) =  Tab(0)
          ! Tab(-2) =  Tab(-1)



          ! Dx = MESH%dx
          ! dt = Run%subDt*RK(to)%s_dt(ii)
          ! time = Run%time + dt

          ! H0 = 0.50D0
          ! A = 0.2D0*H0
          ! k = sqrt(3.0D0*a/(4.0D0*H0**2*(H0+A)))
          ! c = sqrt(g_r__8*(H0+a))
          ! L = DATA%x_out - DATA%x_in
          ! d = -7.0D0
          ! z0 = 0.0D0

          ! Ksi = a/(cosh(k*(z0- c*Time-d))**2)
          ! WRITE(6,*)'ksi',ksi
          ! !READ*

          ! !U = c*(1.0D0 -H0/(Ksi+H0))

          ! u_0 = Tab(0)%qx/Tab(0)%h
          ! c_0 = SQRT(g_r__8*Tab(0)%h)

          ! u_1 = Tab(1)%qx/Tab(1)%h
          ! c_1 = SQRT(g_r__8*Tab(1)%h)


          ! !CALL CL_Subcritique_2Invariants(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u, D)

          ! !CALL CL_Subcritique_Koba(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u, D)
          ! !call CL_Subcritique_Bonneton(ksi, H0, u_0, u_1, c_0, c_1, dt, Dx, u)

          ! ! Definition of the inlet boundary condition


          Tab(0)%h = H0 + ksi

          ! WRITE(6,*)'Tab(0)%h',Tab(0)%h
          ! !READ*

          Tab(0)%qx = c*Ksi    

          Tab(-1) =  Tab(0)
          Tab(-2) =  Tab(-1)



       CASE(2) !generation d'onde sinusoidale

          Dx = MESH%dx
          dt = Run%subDt*RK(to)%s_dt(ii)
          time = Run%time + dt

          H0 = 4.2D0
          eta0 = 0.1D0
          period = 2.5D0

          p = 0.0D0

          Ksi = eta0*SIN(2.0D0*Pi_r__8*Time/period + P)


          PRINT*,'RK(to)%s_dt(ii)',RK(to)%s_dt(ii)
          PRINT*,'ii',ii
          PRINT*,'dt',dt
          PRINT*,'time',time
          !READ*

          !H0 = 0.50D0
          A = 0.2D0*H0
          k = sqrt(3.0D0*a/(4.0D0*H0**2*(H0+A)))
          c = sqrt(g_r__8*(H0+a))
          L = DATA%x_out - DATA%x_in
          d = -7.0D0
          z0 = 0.0D0

          !Ksi = a/(cosh(k*(z0- c*Time-d))**2)

          u_0 = Tab(0)%qx/Tab(0)%h
          c_0 = SQRT(g_r__8*Tab(0)%h)

          u_1 = Tab(1)%qx/Tab(1)%h
          c_1 = SQRT(g_r__8*Tab(1)%h)

          print*,'u_0',u_0
          print*,'u_1',u_1
          print*,'c_0',c_0
          print*,'c_1',c_1
          print*,'Ksi',ksi
          !CALL CL_Subcritique_2Invariants(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u, D)

          !CALL CL_Subcritique_Koba(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u, D)
          call CL_Subcritique_Bonneton(ksi, H0, u_0, u_1, c_0, c_1, dt, Dx, u)

          ! Definition of the inlet boundary condition


          Tab(0)%h = H0 + ksi
          Tab(0)%qx = Tab(0)%h*u           
          print*,'Tab(0)%h',Tab(0)%h
          print*,'Tab(0)%qx',Tab(0)%qx


          Tab(-1) =  Tab(0)
          Tab(-2) =  Tab(-1)


       END SELECT

       ! TopoVf(0) = TopoVf(1)
       ! TopoVf(-1) = TopoVf(2)
       ! TopoVf(-2) = TopoVf(3)

       ! Tab(Nbx+1) = Tab(Nbx)
       ! Tab(Nbx+2) = Tab(Nbx-1)
       ! Tab(Nbx+3) = Tab(Nbx-2)

       ! TopoVf(Nbx+1) = TopoVf(Nbx)
       ! TopoVf(Nbx+2) = TopoVf(Nbx-1)
       ! TopoVf(Nbx+3) = TopoVf(Nbx-2)

    END SELECT

  END SUBROUTINE Boundary_VF


  SUBROUTINE sub__Boundary_Cell(Tab, TabEx, MESH, ZF, DATA, Run, Nt)      

    IMPLICIT NONE

    TYPE(MyRun), INTENT(IN) :: Run
    TYPE(DATA_IN), INTENT(IN) :: DATA
    TYPE( STATE ), DIMENSION(:), POINTER :: Tab
    TYPE( STATE ), DIMENSION(:), POINTER :: TabEx  

    TYPE(MeshDef), INTENT(IN) :: MESH
    REAL(kind = 8), DIMENSION(:), POINTER :: ZF
    REAL(kind = 8), DIMENSION(-2:DATA%Nbx+3) :: hh
    REAL(kind = 8), DIMENSION(:), SAVE, ALLOCATABLE :: HLimT
    REAL(kind = 8) :: Time
    !REAL(kind = 8) :: Dt
    INTEGER, INTENT(in) :: Nt
    !LOGICAL, INTENT(in) :: last_step
    INTEGER  :: Nbx
    REAL(kind = 8) :: H
    REAL(kind = 8) ::  U
    !REAL(kind = 8) ::  V
    REAL(kind = 8) ::  X
    REAL(kind = 8) ::  lambda
    REAL(kind = 8) ::  sigma
    REAL(kind = 8) ::  ksi
    !REAL(kind = 8) ::  discharge ! value of the incident discharge when using oblique incident boundary conditions
    REAL(kind = 8) :: HL
    REAL(kind = 8) :: UL
    !REAL(kind = 8) :: VL
    REAL(kind = 8) :: Hij
    !REAL(kind = 8) :: TT
    REAL(kind = 8) :: S
    REAL(kind = 8) :: A
    REAL(kind = 8) ::  h0
    REAL(kind = 8) :: c0
    REAL(kind = 8) :: eta0
    REAL(kind = 8) :: period
    REAL(kind = 8) :: tbf
    REAL(kind = 8) :: c_1
    REAL(kind = 8) :: u_1
    REAL(kind = 8) :: u_0
    REAL(kind = 8) :: c_0
    REAL(kind = 8) :: Dx
    REAL(kind = 8) :: p
    !REAL(kind = 8) :: ki
    INTEGER :: i
    !INTEGER :: j
    !INTEGER :: icox
    INTEGER :: N2

    !REAL(kind = 8), DIMENSION(:), ALLOCATABLE :: HCox
    !INTEGER :: it
    !REAL(kind = 8) :: TimeCox
    !REAL(kind = 8) :: alpha
    !REAL(kind = 8) :: beta
    REAL(kind = 8) :: epsilon
    REAL(kind = 8) :: atime
    REAL(kind = 8) :: l
    !REAL(kind = 8) :: Htot
    !REAL(kind = 8) :: Ldom
    REAL(kind = 8) :: D
    REAL(kind = 8) :: Dt
    CHARACTER(LEN =* ), PARAMETER :: sub__name = "boundary_cell"
    LOGICAL :: debug

    INTEGER :: sta___tut0
    INTEGER :: i__stat
    INTEGER :: sta___tut1

    Nbx = DATA%Nbx

    Time = Run%subTime
    Dt = Run%subDt

    !Time = Run%Time + Run%step*Run%Dt
    !Dt = Run%Dt

    debug=.FALSE.
    ! NCLL => Indice de conditions aux limites sur le bord gauche
    ! NCLR => Indice de conditions aux limites sur le bord droit

    IF(debug) THEN
       PRINT*, "Entering ", sub__name, " NCLL==", DATA%NCLL, "NCI==", DATA%NCI, " Dt==", Run%subDt, " Time==", Time
    END IF
    ! ##### CONDITIONS INITIALES SUR LES VALEURS #####

    ! Tab(-2)%h = Tab(1)%h
    ! Tab(0)%h = Tab(1)%h
    ! Tab(-1)%h = Tab(1)%h
    ! Tab(Nbx+1)%h = Tab(Nbx)%h
    ! Tab(Nbx+2)%h = Tab(Nbx)%h
    ! Tab(Nbx+3)%h = Tab(Nbx)%h

    ! Tab(-2)%qx = Tab(1)%qx
    ! Tab(0)%qx = Tab(1)%qx
    ! Tab(-1)%qx = Tab(1)%qx
    ! Tab(Nbx+1)%qx = Tab(Nbx)%qx
    ! Tab(Nbx+2)%qx = Tab(Nbx)%qx
    ! Tab(Nbx+3)%qx = Tab(Nbx)%qx

    SELECT CASE(DATA%NCLL)

    CASE(cte__1) 

       ! Neumann homogène

       Tab(0) = Tab(1)
       Tab(-1) = Tab(2)
       Tab(-2) = Tab(3)


    CASE(cte__2)

       !Entrante

       Tab(0) = Tab(1)
       Tab(0)%h = Tab(1)%h + cte__0_d_10_r__8

       Tab(-1) = Tab(0)
       Tab(-2) = Tab(0)


    CASE(cte__3)

       !Paroi

       Tab(0)%h = Tab(1)%h
       Tab(0)%qx = -Tab(1)%qx


       Tab(-1)%h = Tab(2)%h 
       Tab(-1)%qx = -Tab(2)%qx 


       Tab(-2)%h = Tab(3)%h 
       Tab(-2)%qx = -Tab(3)%qx   


    CASE(cte__6)

       !Dirichlet



       SELECT CASE(DATA%NCI)

       CASE default
          Tab(0)%h = cte__1_d_0_r__8
          Tab(0)%qx = cte__1_d_0_r__8


          Tab(-1)%h = cte__1_d_0_r__8
          Tab(-1)%qx = cte__1_d_0_r__8

          Tab(-2)%h = cte__1_d_0_r__8
          Tab(-2)%qx = cte__1_d_0_r__8


       CASE(cte__19)

          Tab(0)%h = Tab(1)%h
          Tab(0)%qx = 0.0D0


          Tab(-1) = Tab(0)
          Tab(-2) = Tab(0)

       CASE(26)

          !Tab(0)%h = Tab(1)%h
          Tab(0)%qx = 4.42D0
          Tab(-1) = Tab(0)
          Tab(-2) = Tab(0)


       END SELECT


    CASE(cte__9)

       !Periodic


       ZF(0) = ZF(Nbx)
       Tab(0)%h = Tab(Nbx)%h
       Tab(0)%qx = Tab(Nbx)%qx

       Tab(0)%topo = ZF(Nbx)

       ZF(-1) = ZF(Nbx-1)
       Tab(-1)%h = Tab(Nbx-1)%h
       Tab(-1)%qx = Tab(Nbx-1)%qx

       Tab(-1)%topo = ZF(Nbx-1)

       ZF(-2) = ZF(Nbx-2)
       Tab(-2)%h =  Tab(Nbx-2)%h
       Tab(-2)%qx = Tab(Nbx-2)%qx

       Tab(-2)%topo = ZF(Nbx-2)



    CASE(cte__7)

       SELECT CASE(DATA%NCI)

       CASE(cte__33)

          SELECT CASE(DATA%onde)

          CASE(cte__1)
             ! Cas overtoping  Forçage sinusoidal en hauteur d'eau

             h0 = cte__3_d_0_r__8
             eta0 = cte__0_d_3_r__8
             c0 = SQRT(g_r__8 * h0 )
             lambda = cte__200_d_0_r__8
             period = lambda / c0

             Ksi = eta0*COS(-cte__2_d_0_r__8*Pi_r__8*c0*Time/lambda &
                  + Pi_r__8/cte__2_d_0_r__8)

             ! Cas overtopping Forçage sinusoidal en hauteur d'eau + modulation basse frequence
          CASE(cte__2)

             h0 = cte__3_d_0_r__8

             eta0 = cte__0_d_2_r__8
             c0 = SQRT(g_r__8 * h0 )
             lambda = cte__20_d_0_r__8
             period = lambda / c0
             tbf = cte__10_d_0_r__8 * period

             Ksi = eta0*COS(-cte__2_d_0_r__8*Pi_r__8*Time/period + &
                  Pi_r__8/cte__2_d_0_r__8)*COS(cte__2_d_0_r__8*Pi_r__8*Time/tbf)

          CASE(cte__3)
             !Onde de Stokes 2eme ordre

             h0 = cte__3_d_0_r__8
             c0 = SQRT(g_r__8 * h0 )
             eta0 = cte__0_d_5
             period = cte__1_d_79_r__8

             Ksi = eta0*( COS(-cte__2_d_0_r__8*Pi_r__8*Time/period ) + &
                  cte__0_d_27_r__8*COS(-cte__4_d_0_r__8*Pi_r__8*Time/period))

          END SELECT

          Dx = (MESH%X_DF(2) - MESH%X_DF(0))/cte__2_d_0_r__8

          IF(Tab(0)%h <= DATA%kappa) THEN
             u_0 = cte__0_d_r__8 
             c_0 =cte__0_d_r__8
          ELSE
             u_0 = Tab(0)%qx/Tab(0)%h
             c_0 = SQRT(g_r__8*Tab(0)%h)
          END IF

          IF(Tab(1)%h <= DATA%kappa) THEN
             u_1 = cte__0_d_r__8 ; c_0 =cte__0_d_r__8
          ELSE
             u_1 = Tab(1)%qx/Tab(1)%h
             c_1 = SQRT(g_r__8*Tab(1)%h)
          END IF

          Tab(0)%h = Ksi + h0
          !print*,'Tab(0)%h',Tab(0)%h

          CALL CL_Subcritique_Koba(ksi, H0, u_0, u_1, c_0, c_1, Dt, Dx, u, D)

          Tab(0)%qx = Tab(0)%h*u



          Tab(-1) = Tab(0)
          Tab(-2) = Tab(0)



       CASE(cte__34)

          SELECT CASE(DATA%onde)

          CASE(cte__1)
             ! Cas overtoping  Forçage sinusoidal en hauteur d'eau

             h0 = cte__10_d_0_r__8
             eta0 = cte__1_d_5_r__8
             c0 = SQRT(g_r__8 * h0 )
             lambda = cte__200_d_0_r__8
             period = lambda / c0

             Ksi = eta0*COS(-cte__2_d_0_r__8*Pi_r__8*c0*Time/lambda + Pi_r__8/cte__2_d_0_r__8)

             ! Cas overtopping Forçage sinusoidal en hauteur d'eau + modulation basse frequence
          CASE(cte__2)

             h0 = cte__10_d_0_r__8

             eta0 = cte__3_d_0_r__8
             c0 = SQRT(g_r__8 * h0 )
             lambda = cte__20_d_0_r__8
             period = lambda / c0
             tbf = cte__10_d_0_r__8 * period

             Ksi = eta0*COS(-cte__2_d_0_r__8*Pi_r__8*Time/period + &
                  Pi_r__8/cte__2_d_0_r__8)*COS(cte__2_d_0_r__8*Pi_r__8*Time/tbf)

          CASE(cte__3)
             !Onde de Stokes 2eme ordre

             h0 = cte__10_d_0_r__8
             c0 = SQRT(g_r__8 * h0 )
             eta0 = cte__3_d_0_r__8
             period = cte__1_d_79_r__8

             Ksi = eta0*( COS(-cte__2_d_0_r__8*Pi_r__8*Time/period ) + &
                  cte__0_d_27_r__8*COS(-cte__4_d_0_r__8*Pi_r__8*Time/period))

          END SELECT



          Dx = (MESH%X_DF(2) - MESH%X_DF(0))/cte__2_d_0_r__8

          IF(Tab(0)%h.EQ.kappa) THEN
             u_0 = cte__0_d_r__8 
             c_0 =cte__0_d_r__8
          ELSE
             u_0 = Tab(0)%qx/Tab(0)%h
             c_0 = SQRT(g_r__8*Tab(0)%h)
          END IF

          IF(Tab(1)%h.EQ.kappa) THEN
             u_1 = cte__0_d_r__8 ; c_0 =cte__0_d_r__8
          ELSE
             u_1 = Tab(1)%qx/Tab(1)%h
             c_1 = SQRT(g_r__8*Tab(1)%h)
          END IF

          Tab(0)%h = Ksi + h0

          CALL CL_Subcritique_Koba(ksi, H0, u_0, u_1, c_0, c_1, Dt, Dx, u, D)

          Tab(0)%qx = Tab(0)%h*u



          Tab(-1) = Tab(0)
          Tab(-2) = Tab(0)



       CASE(cte__35)

          SELECT CASE(DATA%onde)

          CASE(cte__1)
             ! Cas overtoping  Forçage sinusoidal en hauteur d'eau

             h0 = cte__11_d_0_r__8
             eta0 = cte__0_d_5_r__8
             c0 = SQRT(g_r__8 * h0 )
             lambda = cte__200_d_0_r__8
             period = lambda / c0

             Ksi = eta0*COS(-cte__2_d_0_r__8*Pi_r__8*c0*Time/lambda + Pi_r__8/cte__2_d_0_r__8)

             ! Cas overtopping Forçage sinusoidal en hauteur d'eau + modulation basse frequence
          CASE(cte__2)

             h0 = cte__11_d_0_r__8

             eta0 = cte__0_d_5_r__8
             c0 = SQRT(g_r__8 * h0 )
             lambda = cte__20_d_0_r__8
             period = lambda / c0
             tbf = cte__10_d_0_r__8 * period

             Ksi = eta0*COS(-cte__2_d_0_r__8*Pi_r__8*Time/period + &
                  Pi_r__8/cte__2_d_0_r__8)*COS(cte__2_d_0_r__8*Pi_r__8*Time/tbf)

          CASE(cte__3)
             !Onde de Stokes 2eme ordre

             h0 = cte__11_d_0_r__8
             c0 = SQRT(g_r__8 * h0 )
             eta0 = cte__1_d_5_r__8
             period = cte__1_d_79_r__8

             Ksi = eta0*( COS(-cte__2_d_0_r__8*Pi_r__8*Time/period ) + &
                  cte__0_d_27_r__8*COS(-cte__4_d_0_r__8*Pi_r__8*Time/period))

          END SELECT



          Dx = (MESH%X_DF(2) - MESH%X_DF(0))/cte__2_d_0_r__8

          IF(Tab(0)%h == kappa) THEN
             u_0 = cte__0_d_r__8 
             c_0 =cte__0_d_r__8
          ELSE
             u_0 = Tab(0)%qx/Tab(0)%h
             c_0 = SQRT(g_r__8*Tab(0)%h)
          END IF

          IF(Tab(1)%h.EQ.kappa) THEN
             u_1 = cte__0_d_r__8 ; c_0 =cte__0_d_r__8
          ELSE
             u_1 = Tab(1)%qx/Tab(1)%h
             c_1 = SQRT(g_r__8*Tab(1)%h)
          END IF

          Tab(0)%h = Ksi + h0

          CALL CL_Subcritique_Koba(ksi, H0, u_0, u_1, c_0, c_1, Dt, Dx, u, D)

          Tab(0)%qx = Tab(0)%h*u


          Tab(-1) = Tab(0)
          Tab(-2) = Tab(0)



       CASE(cte__38)

          ! test Jonswap

          h0 = cte__3_d_0_r__8

          IF(NT == 1) THEN

             !  read Jonswap file
             OPEN (UNIT = 328, FILE = TRIM(dir_data)//"JSwap.txt", FORM = "FORMATTED", STATUS = "OLD", &
                  ACTION = "READ", POSITION = "REWIND", IOSTAT = sta___tut0 )

             ALLOCATE (HLimT(1: DATA%NtMax ), STAT = i__stat)
             IF (debug) THEN
                PRINT *, sub__name, " ALLOCATE HLimT"
             END IF
             IF (i__stat /= 0) THEN
                STOP "alloc_err HLimT"
             END IF

             !DO i=1,DATA%NtMax

             READ(UNIT = 328, IOSTAT = sta___tut1, FMT = *)(HLimT(i), i=1, DATA%NtMax)
             ! IF (sta___tut1 /= 0) THEN
             !    CALL read__stop( "read_err jSWap.txt", sta___tut1)
             ! END IF

             !END DO

          END IF

          h0 = cte__3_d_0_r__8
          eta0 = cte__0_d_3_r__8
          c0 = SQRT(g_r__8 * h0 )
          lambda = cte__200_d_0_r__8
          period = lambda / c0
          PRINT*,'period',period, 1.0/period

          !Ksi = eta0*COS(-cte__2_d_0_r__8*Pi_r__8*c0*Time/lambda &
          !     + Pi_r__8/cte__2_d_0_r__8)

          Ksi = HLimT(NT)
          PRINT*,'Ksi',Ksi



          Dx = (MESH%X_DF(2) - MESH%X_DF(0))/cte__2_d_0_r__8

          IF(Tab(0)%h <= DATA%kappa) THEN
             u_0 = cte__0_d_r__8 
             c_0 =cte__0_d_r__8
          ELSE
             u_0 = Tab(0)%qx/Tab(0)%h
             c_0 = SQRT(g_r__8*Tab(0)%h)
          END IF

          IF(Tab(1)%h <= DATA%kappa) THEN
             u_1 = cte__0_d_r__8 ; c_0 =cte__0_d_r__8
          ELSE
             u_1 = Tab(1)%qx/Tab(1)%h
             c_1 = SQRT(g_r__8*Tab(1)%h)
          END IF

          Tab(0)%h = Ksi + h0

          CALL CL_Subcritique_Koba(ksi, H0, u_0, u_1, c_0, c_1, Dt, Dx, u, D)

          Tab(0)%qx = Tab(0)%h*u



          Tab(-1) = Tab(0)
          Tab(-2) = Tab(0)



       CASE(cte__18)

          !tidal wave flow
          !print*,'je suis la'

          X = ( MESH%X_DF(i) +  MESH%X_DF(i+1))/cte__2_d_0_r__8

          Tab(0)%h = 64.5D0 - 4.0D0*SIN(Pi_r__8*(4*Time/86400.0D0 + 1.0D0/2.0D0))  
          !print*,'Tab(0)%h',Tab(0)%h
          Tab(0)%qx = 0.0D0

          Tab(-1) = Tab(0)
          Tab(-2) = Tab(0)

       CASE(cte__36)

          !Carrier & Greenspan transient case

          S = 1.0D0 / cte__50_d_0_r__8          ! slope
          epsilon = cte__0_d_1_r__8        ! amplitude
          A = cte__1_d_5_r__8 * SQRT (( 1.0D0 + cte__0_d_9_r__8 * epsilon ) )
          l = cte__20_d_0_r__8     ! wave length
          atime = Time / SQRT (( l /( S * g_r__8 ) ) ) ! non-dimensionalized time

          N2 = 1

          DO i = 0, -2, -1
             CALL Solex_CG_Transient(i, MESH%X_DF, Time, lambda, sigma)
             X = ( MESH%X_DF(i) +  MESH%X_DF(i+1))/cte__2_d_0_r__8
             X = X/l

             ! computation of non-dimensionalized ksi 
             Ksi = S*l*(X + a**2*sigma**2/cte__16_d_0_r__8)
             Tab(i)%h = MAX(Ksi-ZF(i),cte__0_d_r__8)
             !print*,'Tab(i)%h',i,Tab(i)%h

             ! computation of non-dimensionalized u
             Hij = Tab(i)%h
             Tab(i)%qx = SQRT((g_r__8*S*l))*Hij*((1.0D0/cte__2_d_0_r__8)*A*lambda-atime)

          END DO


       CASE(cte__37)

          !Onde non-linéaire 1D sur plage - Carrier et Greenspan cas periodic

          ! On adimensionne le temps
          !print*,'Time dans boundary',Time
          S = cte__1_d_0_r__8/cte__30_d_0_r__8            ! pente       
          l = cte__20_d_0_r__8  ! longueur d'onde
          !N1 = 0
          aTime = Time/(SQRT(l/(S*g_r__8)))

          N2 = 1
          DO i=0,-2,-1

             CALL Solex_CG_Periodic(i, MESH%X_DF, Time, lambda, sigma)

             X = ( MESH%X_DF(i) +  MESH%X_DF(i+1))/cte__2_d_0_r__8
             X = X/l

             ! calcul de ksi dimensionné
             Ksi = S*l*(X + sigma**2/cte__16_d_0_r__8)
             Tab(i)%h = MAX(Ksi-ZF(i),cte__0_d_r__8)

             ! calcul de u dimensionné
             Hij = Tab(i)%h
             Tab(i)%qx = SQRT(g_r__8*S*l)*Hij*((1.0D0/cte__2_d_0_r__8)*lambda - aTime)


          END DO



       CASE(cte__90) ! First Boussinesq validation case

          IF(debug) THEN
             PRINT*,'Dans M_Boundary'
             PRINT*,'Time',Time
          END IF

          H0 = cte__4_d_2_r__8
          eta0 = cte__0_d_1_r__8
          period = cte__2_d_5_r__8

          p = cte__0_d_r__8
          !omega =cte__2_d_0_r__8*Pi/period

          !k = 2.0_r__8*Pi_r__8/9.97

          Ksi = eta0*SIN(2.0D0*Pi_r__8*Time/period + P)

          IF (debug) THEN
             PRINT*,"ksi",ksi, " eta0==", eta0, " Time==", Time, " period==", period, " P==", P
          END IF

          Dx = (MESH%X_DF(2) - MESH%X_DF(1))/2.0D0 + (MESH%X_DF(1) - MESH%X_DF(0))/2.0D0

          IF(Tab(0)%h == DATA%kappa) THEN
             u_0 = cte__0_d_r__8 
             c_0 =cte__0_d_r__8
          ELSE
             u_0 = Tab(0)%qx/Tab(0)%h
             c_0 = SQRT(g_r__8*Tab(0)%h)
          END IF
          IF(Tab(1)%h == DATA%kappa) THEN
             u_1 = cte__0_d_r__8 
             c_1 =cte__0_d_r__8
          ELSE
             u_1 = Tab(1)%qx/Tab(1)%h
             c_1 = SQRT(g_r__8*Tab(1)%h)
          END IF

          CALL CL_Subcritique_2Invariants(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u, D)

          !CALL CL_Subcritique_Koba(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u, D)
          !call CL_Subcritique_Bonneton(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u)

          ! Definition of the inlet boundary condition


          Tab(0)%h = D
          Tab(0)%qx = Tab(0)%h*u           


          Tab(-1) =  Tab(0)
          Tab(-2) =  Tab(-1)



       CASE(cte__91) ! sinusoidal propagation in shallow water

          IF(debug) THEN
             PRINT*,'Dans M_Boundary'
             PRINT*,'Time',Time
          END IF

          H0 = cte__10_d_0_r__8
          eta0 = cte__0_d_04_r__8
          period = cte__20_d_193_r__8
          p = cte__0_d_r__8
          !omega = cte__2_d_0_r__8*Pi_r__8/period
          !k = 2.0*Pi_r__8/9.97
          !k = 2.0*Pi_r__8/9.68

          Ksi = eta0*SIN(cte__2_d_0_r__8*Pi_r__8*Time/period + P)

          Dx = (MESH%X_DF(2) - MESH%X_DF(1))/cte__2_d_0_r__8 + (MESH%X_DF(1) - MESH%X_DF(0))/cte__2_d_0_r__8

          IF(Tab(0)%h == kappa) THEN
             u_0 = cte__0_d_r__8 
             c_0 =cte__0_d_r__8
          ELSE
             u_0 = Tab(0)%qx/Tab(0)%h
             c_0 = SQRT(g_r__8*Tab(0)%h)
          END IF
          IF(Tab(1)%h == kappa) THEN
             u_1 = cte__0_d_r__8 ; c_1 =cte__0_d_r__8
          ELSE
             u_1 = Tab(1)%qx/Tab(1)%h
             c_1 = SQRT(g_r__8*Tab(1)%h)
          END IF

          CALL CL_Subcritique_2Invariants(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u, D)
          !CALL CL_Subcritique_Koba(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u, D)
          !CALL CL_Subcritique_Bonneton(ksi, H0, u_0, u_1, c_0, c_1, Run%Dt, Dx, u)
          !u = omega/(k*H0)*Ksi

          ! Definition of the inlet boundary condition

          !Tab(0)%h = H0 + Ksi
          Tab(0)%h = D
          Tab(0)%qx = Tab(0)%h*u

          !Tab(0)%qx = sqrt(g/H0)*Ksi*(H0+Ksi)
          !Tab(0)%qy = cte__0_d_r__8

          Tab(-1) =  Tab(0)
          Tab(-2) =  Tab(-1)

          ! ### COTE DROIT ###
       END SELECT
    END SELECT

    SELECT CASE(DATA%NCLR)

    CASE(cte__1)

       ! Neumann 

          Tab(Nbx+1) = Tab(Nbx)
          Tab(Nbx+2) = Tab(Nbx-1)
          Tab(Nbx+3) = Tab(Nbx-2)

    CASE(cte__3)

       !solid wall
 
          Tab(Nbx+1)%h = Tab(Nbx)%h
          Tab(Nbx+1)%qx = - Tab(Nbx)%qx 


          Tab(Nbx+2)%h = Tab(Nbx-1)%h
          Tab(Nbx+2)%qx = - Tab(Nbx-1)%qx 


          Tab(Nbx+3)%h = Tab(Nbx-2)%h
          Tab(Nbx+3)%qx = - Tab(Nbx-2)%qx 

    CASE(cte__6)

       !Dirichlet
   

          SELECT CASE(DATA%NCI)

          CASE default
             Tab(Nbx+1)%h = cte__1_d_0_r__8
             Tab(Nbx+1)%qx = cte__1_d_0_r__8


             Tab(Nbx+2) = Tab(Nbx+1)
             Tab(Nbx+3) = Tab(Nbx+1)


          CASE(cte__18)
             !print*,'je suis la 2'
             !Tab(Nbx+1)%h = cte__1_d_0_r__8
             Tab(Nbx+1)%qx = 0.0D0


             Tab(Nbx+2) = Tab(Nbx+1)
             Tab(Nbx+3) = Tab(Nbx+1)

          CASE(26)
             !print*,'je suis la 2'
             Tab(Nbx+1)%h = 2.0D0
             !Tab(Nbx+1)%qx = 0.0D0
             !Tab(Nbx+1)%qy = 0.0D0

             Tab(Nbx+2) = Tab(Nbx+1)
             Tab(Nbx+3) = Tab(Nbx+1)

          END SELECT

    CASE(cte__9)

       !Periodic

          ZF(Nbx+1) = ZF(1)
          Tab(Nbx+1)%topo = ZF(1)
          Tab(Nbx+1)%h = Tab(1)%h
          Tab(Nbx+1)%qx = Tab(1)%qx


          ZF(Nbx+2) = ZF(2)
          Tab(Nbx+2)%topo = ZF(2)
          Tab(Nbx+2)%h = Tab(2)%h
          Tab(Nbx+2)%qx = Tab(2)%qx


          ZF(Nbx+3) = ZF(3)
          Tab(Nbx+3)%topo = ZF(3)
          Tab(Nbx+3)%h = Tab(3)%h
          Tab(Nbx+3)%qx = Tab(3)%qx

    CASE(cte__8)


       !Condition de radiation (transparente)
       ! A utiliser avec une "sponge layer" éventuellement

       ! Filling of empty "temp" ghost cells
    
          Tab(Nbx+1) = Tab(Nbx)
          Tab(Nbx+2) = Tab(Nbx)
          Tab(Nbx+3) = Tab(Nbx)

          loop__504 : DO i = -2, Nbx+3
             hh(i) = Tab(i)%h
          END DO loop__504
    
          Dx = (MESH%X_DF(2) - MESH%X_DF(0))/2.0D0

          ! calcal de h par condition de radiation
          Tab(Nbx+1)%h = hh(Nbx+1) - Run%Dt*SQRT(g_r__8*hh(Nbx+1))*(hh(Nbx+1)-hh(Nbx))/Dx

          ! calcul de u par choc ou detente compatible
          H =  Tab(Nbx)%h
          U = Tab(Nbx)%qx/Tab(Nbx)%h

          HL = Tab(Nbx+1)%h
          CALL COTE_RIGHT(H, U, HL, UL)

          !   ! Calculé avec la théorie linéaire Saint-Venant
          !   Eta_s = Tab(Nbx+1)%h-d(Nbx+1)
          !   Tab(Nbx+1)%qx = (Eta_s*omega/(k*H0))*Tab(Nbx+1)%h

          Tab(Nbx+1)%qx = HL*UL


          Tab(Nbx+2) = Tab(Nbx+1)        
          Tab(Nbx+3) = Tab(Nbx+1)


    END SELECT

  END SUBROUTINE Sub__Boundary_Cell

  SUBROUTINE DEBIT_LEFT(QINT,HR,UR,H,U)

    IMPLICIT NONE

    REAL (KIND = 8 ), INTENT ( IN ) :: QINT
    REAL (KIND = 8 ), INTENT ( IN ) :: HR
    REAL (KIND = 8 ), INTENT ( IN ) :: UR

    REAL (KIND = 8 ), INTENT ( OUT ) :: H
    REAL (KIND = 8 ), INTENT ( OUT ) :: U

    REAL (KIND = 8 ) :: H1
    REAL (KIND = 8 ) :: H2
    REAL (KIND = 8 ) :: HINT
    REAL (KIND = 8 ) :: UINT
    REAL (KIND = 8 ) :: ERR
    INTEGER :: N
    CHARACTER (LEN =* ), PARAMETER :: sub__name = "debit_left"
    LOGICAL  :: debug

    debug = .FALSE.


    IF (debug) THEN
      PRINT *, "Entering ", sub__name
    END IF

    ERR = cte__1_d_0e_m_5_r__8 * HR
    N = 0
    H1 = cte__0_d_r__8
    H2 = cte__10_d_0_r__8 * HR

    loop__173 : DO WHILE ((ABS(H2-H1) >= ERR))

       HINT = (H1+H2)/cte__2_d_0_r__8
!       UINT = F(QINT,HR,UR,HINT)
      UINT = Fct ( HR, UR, HINT )

       IF (HINT*UINT > QINT) THEN
          H2 = HINT
       ELSE
          H1 = HINT
       END IF
       N = N+1

    END DO loop__173


    H = (H1+H2)/cte__2_d_0_r__8
    U = QINT/H

  END SUBROUTINE DEBIT_LEFT

  ELEMENTAL FUNCTION Fct(HR, UR, H) RESULT (f__res)
    !USE fortranf
    !USE les__ctes__z2f
    IMPLICIT NONE
    REAL (KIND = 8 ) :: f__res
    REAL (KIND = 8 ), INTENT ( IN ) :: HR
    REAL (KIND = 8 ), INTENT ( IN ) :: UR
    REAL (KIND = 8 ), INTENT ( IN ) :: H
    !CHARACTER (LEN =* ), PARAMETER :: mod__name = "fct"

    IF (H < HR ) THEN ! Detente
      f__res = UR + cte__2_d_0_r__8 *( SQRT(g_r__8* H ) - SQRT (g_r__8* HR) )
    ELSE ! Choc
      f__res = UR +( H - HR) * SQRT ( ( ( g_r__8 *( H + HR ))/(cte__2_d_0_r__8*H*HR)))
    END IF

  END FUNCTION Fct

!   FUNCTION F(QINT,HR,UR,H)

!     IMPLICIT NONE

!     REAL(kind = r__8):: F
!     REAL(kind = r__8), INTENT(in):: QINT, HR, UR, H

!     IF (H.LT.HR) THEN
!        ! Detente
!        F = UR + cte__2_d_0_r__8*(SQRT(g_r__8*H) - SQRT(g_r__8*HR))
!     ELSE
!        ! Choc
!        F = UR + (H-HR)*SQRT(g_r__8*(H+HR)/(cte__2_d_0_r__8*H*HR))
!     END IF

!   END FUNCTION F

  SUBROUTINE COTE_RIGHT(HL, UL, H, U)

    IMPLICIT NONE

    REAL (KIND = 8), INTENT ( IN ) :: HL
    REAL (KIND = 8), INTENT ( IN ) :: UL

    REAL (KIND = 8), INTENT ( OUT ) :: U
  
    REAL (KIND = 8), INTENT ( IN OUT ) :: H
    CHARACTER (LEN =*), PARAMETER :: sub__name = "cote_right"
    LOGICAL :: debug

    debug = .FALSE.

    IF (debug) THEN
       PRINT *, "Entering ", sub__name
    END IF

    IF (UL <= SQRT(g_r__8*HL)) THEN
       IF (H <= HL) THEN
          U = UL - cte__2_d_0_r__8*(SQRT(g_r__8*H) - SQRT(g_r__8*HL))
         
       ELSE 
          U = UL - (H-HL)*SQRT(g_r__8*(H+HL)/(cte__2_d_0_r__8*H*HL))
        
       END IF
    ELSE
       H = HL
       U = UL
     
    END IF
    !
  END SUBROUTINE COTE_RIGHT


  SUBROUTINE CL_Subcritique_Bonneton(ksi, H0, u0, u1, c0, c1, Dt, Dx, u)

! Computation of the limit velocity for a real and imposed water depth (issued from measurements)
!(i.e; the water depth takes into account both incident and reflected waves)
! Only the R- invariant is computed and used

!Ksi stands for real wave amplitude including the possible reflected waves 

    IMPLICIT NONE


    REAL(kind = 8), INTENT(in) :: Ksi
    REAL(kind = 8), INTENT(in) :: H0
    REAL(kind = 8), INTENT(in) :: u0
    REAL(kind = 8), INTENT(in) :: u1
    REAL(kind = 8), INTENT(in) :: c0
    REAL(kind = 8), INTENT(in) :: c1
    REAL(kind = 8), INTENT(in) :: Dt
    REAL(kind = 8), INTENT(in) :: Dx
    REAL(kind = 8), INTENT(out) :: u
    REAL(kind = 8) :: dp
    REAL(kind = 8) :: un
    REAL(kind = 8) :: cn
    REAL(kind = 8) :: dpt
    REAL(kind = 8) :: c
    REAL(kind = 8) :: Ip
    REAL(kind = 8) :: S
    REAL(kind = 8) :: D
    INTEGER :: i


    !I0 = u0 - 2.*c0
    !I1 = u1 - 2.*c1
    !c = sqrt(g*(H0+Ksi))

    !u = (I0 + 2.*c*(1. + Dt/(2.*dx)*(I1-I0)))/(1. + Dt/dx*(I1-I0))

    dp=cte__0_d_5_r__8
 
    !print*,'dt', dt
    !print*,'dx', dx
    i=0
    DO WHILE((ABS(dp-dpt).GT.1.D-10).AND.(i.LT.50))
       dpt = dp

       un = u0 + (dpt/Dx)*(u1 - u0)
       cn = c0 + (dpt/Dx)*(c1 - c0)

       dp= Dt*(cn-un)

       i=i+1
       !PRINT*,'i dans 1ere boucle',i
       IF(i.EQ.50) PRINT*,"ERROR -- Inlet_Left does not converge"
    END DO

    un = u0 + (dp/Dx)*(u1 - u0)
    cn = c0 + (dp/Dx)*(c1 - c0)

    Ip = un - 2.D0*cn 

    D = H0 + Ksi 
    c = SQRT(g_r__8*D)
    S = cte__1_d_0_r__8/35.0D0
    u = Ip + 2.D0*c + S*dt

  END SUBROUTINE CL_Subcritique_Bonneton


  SUBROUTINE CL_Subcritique_Koba(ksi, H0, u0, u1, c0, c1, Dt, Dx, u, D)

    ! Calcul la vitesse reelle et de la hauteur d'eau reelle a la limite pr une hauteur d'eau de l'onde incidente donnée.
    ! Seul l'invariant R- est calculé et utilisé mais on ajoute une contribution due 
    ! a l'onde reflechie (Kobayashi et al.) et une correction non-lineaire a l'ordre 2 (Hu et al.)

    ! ksi represente l'amplitude de l'onde incidente

    IMPLICIT NONE


    REAL(kind = 8), INTENT(in) :: Ksi
    REAL(kind = 8), INTENT(in) :: H0
    REAL(kind = 8), INTENT(in) :: u0
    REAL(kind = 8), INTENT(in) :: u1
    REAL(kind = 8), INTENT(in) :: c0
    REAL(kind = 8), INTENT(in) :: c1
    REAL(kind = 8), INTENT(in) :: Dt
    REAL(kind = 8), INTENT(in) :: Dx
    REAL(kind = 8), INTENT(out) :: u
    REAL(kind = 8), INTENT(out) :: D
    REAL(kind = 8) :: dp
    REAL(kind = 8) :: un
    REAL(kind = 8) :: cn
    REAL(kind = 8) :: dpt
    REAL(kind = 8) :: c
    REAL(kind = 8) :: Ip
    REAL(kind = 8) :: eta_r
    REAL(kind = 8) :: eta_rt
    !REAL(kind = 8) :: S
    INTEGER :: i

    dp=cte__0_d_5_r__8

    i=0
    DO WHILE((ABS(dp-dpt).GT.1.D-10).AND.(i.LT.50))
       dpt = dp

       un = u0 + (dpt/Dx)*(u1 - u0)
       cn = c0 + (dpt/Dx)*(c1 - c0)

       dp= Dt*(cn-un)

       i=i+1
       !PRINT*,'i dans 1ere boucle',i
       IF(i.EQ.50) PRINT*,"Inlet_Left n'as pas convergé...erreur"
    END DO

    un = u0 + (dp/Dx)*(u1 - u0)
    cn = c0 + (dp/Dx)*(c1 - c0)

    Ip = un - 2.D0*cn 

    ! Evaluation de l'amplitude de l'onde réfléchie (formule de Kobayashi)

    eta_r = -(1.0/2.0)*SQRT(h0/g_r__8)*Ip - h0 

    ! Correction non linéaire ordre 2 (methode de Hu et al.)
    i=0
    DO WHILE((ABS(eta_rt-eta_r).GT.1.D-10).AND.(i.LT.cte__50))

       eta_rt = eta_r

       D = H0 + Ksi + eta_rt
       eta_r = -cte__0_d_5_r__8*SQRT(H0/g_r__8)*Ip - H0 + (SQRT(D) - SQRT(H0))**2/cte__2_d_0_r__8

       i=i+1
       !PRINT*,'i dans 2ere boucle',i
       IF(i.EQ.50) PRINT*,"Inlet_Left n'as pas convergé...erreur"
    END DO


    D = H0 + Ksi + eta_r
    c = SQRT(g_r__8*D)
    !S = cte__1_d_0_r__8/35.0D0
    u = Ip + 2.D0*c !+ S*dt

  END SUBROUTINE CL_Subcritique_Koba

  SUBROUTINE CL_Subcritique_2Invariants(ksi, H0, u0, u1, c0, c1, Dt, Dx, u, D)

! Calcul la vitesse reelle et de la hauteur d'eau reelle a la limite pr une hauteur d'eau de l'onde incidente donnée.
! les'invariant R- et R+ sont calculés et utilisés pour tenir compte des amplitudes incidentes et reflehies

    IMPLICIT NONE


    REAL(kind = 8), INTENT(in) :: Ksi
    REAL(kind = 8), INTENT(in) :: H0
    REAL(kind = 8), INTENT(in) :: u0
    REAL(kind = 8), INTENT(in) :: u1
    REAL(kind = 8), INTENT(in) :: c0
    REAL(kind = 8), INTENT(in) :: c1
    REAL(kind = 8), INTENT(in) :: Dt
    REAL(kind = 8), INTENT(in) :: Dx
    REAL(kind = 8), INTENT(out) :: u
    REAL(kind = 8), INTENT(out) :: D
    REAL(kind = 8) :: dp
    REAL(kind = 8) :: un
    REAL(kind = 8) :: cn
    REAL(kind = 8) :: dpt
    !REAL(kind = 8) :: c
    REAL(kind = 8) :: Ipm
    REAL(kind = 8) :: Ipp
    !REAL(kind = 8) :: I1
    !REAL(kind = 8) :: I0
    INTEGER :: i


! Calcul de l'invariant sortant R-
    dp=cte__0_d_5_r__8

    i=0
    DO WHILE((ABS(dp-dpt).GT.1.D-10).AND.(i.LT.50))
       dpt = dp

       un = u0 + (dpt/Dx)*(u1 - u0)
       cn = c0 + (dpt/Dx)*(c1 - c0)

       dp= Dt*(cn-un)

       i=i+1
       !PRINT*,'i dans 1ere boucle',i
       IF(i.EQ.50) PRINT*,"Inlet_Left n'as pas convergé...erreur"
    END DO

    un = u0 + (dp/Dx)*(u1 - u0)
    cn = c0 + (dp/Dx)*(c1 - c0)

    Ipm = un - 2.D0*cn 

    ! Calcul de l'invariant entrant R+ par la théorie linéaire en utilisant 
    ! l'amplitude de l'onde incidente et la profondeur au repos
    Ipp = 2.0*SQRT(g_r__8*(H0+Ksi)) + SQRT(g_r__8/h0)*Ksi

    !S = cte__1_d_0_r__8/35.0D0
    u = (1.0/2.0) * (Ipm + Ipp)
    D = 1.0/(16.0*g_r__8) * (Ipp - Ipm)**2
    

  END SUBROUTINE CL_Subcritique_2Invariants


  ! SUBROUTINE Sponge_layer(Tab, X_DF, DATA)

  !   IMPLICIT NONE

  !   TYPE( STATE ), DIMENSION(:), POINTER :: Tab
  !    REAL(kind = r__8), DIMENSION(:), POINTER :: X_DF
  !   REAL(kind = r__8), DIMENSION(-2:DATA%Nbx+3) :: mu
  !   REAL(kind = r__8), DIMENSION(-2:DATA%Nbx+3) :: U
  !   REAL(kind = r__8), DIMENSION(-2:DATA%Nbx+3) :: SEta
  !   TYPE(DATA_IN), INTENT(IN) :: DATA
  !   REAL(kind = r__8) :: D_Layer
  !   REAL(kind = r__8) :: alpha
  !   REAL(kind = r__8) :: X
  !   REAL(kind = r__8) :: dist
  !   REAL(kind = r__8) :: Dx
  !   REAL(kind = r__8) ::  H0
  !   INTEGER :: i


  !   PRINT*,'M_Boundary : entering sub SPONGE_LAYER ..'
  !   SELECT CASE(DATA%NCI)

  !   CASE(90)

  !      H0 = 4.2D0
  !      D_Layer = 15.D0

  !   CASE(91)

  !      H0 = cte__10_d_0_r__8
  !      D_Layer = 200.D0

  !   CASE default

  !      PRINT*,'M_Boundary : sub SPONGE_LAYER : do not use a sponge layer for this test ...'

  !   END SELECT

  !   alpha  = 4.D0


  !      DO i =-2, DATA%Nbx+3

  !         X = (X_DF(i) + X_DF(i+1))/2.0
  !         !Y = (X_DF(i)%Y + X_DF(i,j+1)%Y)/2.0

  !         Dx = X_DF(i+1) - X_DF(i)
  !         dist = DATA%DomlenX - X

  !         IF(dist >= D_layer) THEN
  !            mu(i) = cte__1_d_0_r__8
  !         ELSE
  !            mu(i) = EXP((2.**(-dist/Dx) - 2.**(-D_layer/Dx))*LOG(alpha))
  !         END IF
  
  !      END DO
  
  
  !   PRINT*,' M_Boundary : sub SPONGE LAYER : absorbing layer application ...'
  
  
  !      DO i=-2, DATA%Nbx+3
  
  !         SEta(i) =  Tab(i)%h - H0  ! Eta = variations de la surface libre
  !         U(i) = Tab(i)%qx/Tab(i)%h
  
  !      END DO
  
  
  
  
  !      DO i =-2, DATA%Nbx+3
  
  !         SEta(i) = SEta(i)/mu(i)
  !         Tab(i)%h = H0+SEta(i)
  !         U(i) = U(i)/mu(i)
  !         Tab(i)%qx = Tab(i)%h*U(i)
  
  
  !      END DO
  
  
  ! END SUBROUTINE Sponge_layer
  
  subroutine second_degre(x1 , x2 , a , b , c)
    implicit none
    real*8 , intent(inout) :: x1 , x2
    real*8 , intent(in) :: a , b , c
    real*8 :: delta
    delta = (b**2) - (4.0d0*a*c)
    x1 = (-b - sqrt(delta))/(2*a)
    x2 = (-b + sqrt(delta))/(2*a)
  end subroutine second_degre


  function zero_fonction(te,y,N,x_in,x_out,Dt)
    implicit none
    real*8 :: zero_fonction
    real*8,dimension(:),allocatable,intent(in) :: te , y
    integer , intent(in) :: N
    integer :: cpt = 0
    real*8 :: x0 , x1 , xx_in , xx_out
    real*8,parameter :: eps = 1.0d-06
    real*8,intent(in) :: x_in , Dt , x_out
    xx_in = min(x_in,x_out)
    xx_out = max(x_in,x_out)
    x0 = x_in
    !print*,"iteration =   ",cpt
10  x1 = x0 - (fonction(x0,x_in,Dt,te,y,N)/derivee_fonction(x0,Dt,te,y,N))

    do while(abs((fonction(x1,x_in,Dt,te,y,N))) >= eps)
       cpt = cpt + 1
       x0 = x1
       x1 = x0 - (fonction(x0,x_in,Dt,te,y,N)/derivee_fonction(x0,Dt,te,y,N))
       !print*,"iteration =   ",cpt,"    ","valeur = ","   ", fonction(x0,x_in,Dt,te,y,N)
    end do
    if(x1 < xx_in .or. x1 > xx_out) then
       !print*,"je suis ici"
       !goto 10
    end if

    zero_fonction = x1
  end function zero_fonction


  function interpolation(x,te,y,N)
    implicit none
    real*8,dimension(:),allocatable,intent(in) :: te , y
    integer , intent(in) :: N
    real*8 ::  produit , somme
    integer :: i , j
    real*8 :: interpolation
    real*8 , intent(in) :: x

    produit = 1.0d0
    somme = 0.0d0
    do i = 1 , N
       do j = 1 , N
          if(j /= i) then
             produit = produit*((x - te(j))/(te(i) - te(j)))
          end if
       end do
       somme = somme +  y(i)*produit
       produit = 1.0d0
    end do


    interpolation = somme

  end function interpolation

  function derivee_interpolation(x , te , y , N)
    implicit none
    real*8,dimension(:),allocatable,intent(in) :: te , y
    integer , intent(in) :: N
    integer :: i , j , k
    real*8 :: produit , product, somme , sum , prod
    real*8 :: derivee_interpolation
    real*8 , intent(in) :: x

    somme = 0.0d0
    sum = 0.0d0
    product = 1.0d0
    produit = 1.0d0
    prod = 1.0d0


    do i = 1 , N
       do j = 1 , N
          if( j /= i ) then
             do k = 1 , N
                if(k /= j .and. k /= i) then
                   product = product*(x - te(k))

                   produit = (te(i) - te(j))
                   prod = prod*produit
                   sum = sum + product
                end if
             end do
             product = 1.0d0
          end if
       end do
       somme = somme +  y(i)*sum/prod
       prod = 1.0d0
       sum = 0.0d0
    end do


    derivee_interpolation = somme
  end function derivee_interpolation



  function fonction(x,x_in,Dt,te,y,N)
    implicit none
    real*8,dimension(:),allocatable,intent(in) :: te , y
    integer,intent(in) :: N
    real*8 :: fonction
    real*8 , intent(in) :: x , x_in , Dt
    real*8 ::  produit , somme
    integer :: i , j
    produit = 1.0d0
    somme = 0.0d0
    do i = 1 , N
       do j = 1 , N
          if(j /= i) then
             produit = produit*((x - te(j))/(te(i) - te(j)))
          end if
       end do
       somme = somme +  y(i)*produit
       produit = 1.0d0
    end do

    fonction = x - x_in + Dt*(somme)

  end function fonction

  function derivee_fonction(x ,Dt, te , y , N)
    implicit none
    real*8,dimension(:),allocatable,intent(in) :: te , y
    integer,intent(in) :: N
    real*8 :: derivee_fonction
    real*8 , intent(in) :: x , Dt
    integer :: i , j , k
    real*8 :: produit , product, somme , sum , prod
    somme = 0.0d0
    sum = 0.0d0
    product = 1.0d0
    produit = 1.0d0
    prod = 1.0d0


    do i = 1 , N
       do j = 1 , N
          if( j /= i ) then
             do k = 1 , N
                if(k /= j .and. k /= i) then
                   product = product*(x - te(k))

                   produit = (te(i) - te(j))
                   prod = prod*produit
                   sum = sum + product
                end if
             end do
             product = 1.0d0
          end if
       end do
       somme = somme +  y(i)*sum/prod
       prod = 1.0d0
       sum = 0.0d0
    end do


    derivee_fonction = 1.0d0 + Dt*(somme)
  end function derivee_fonction

  function fg(x ,  q , I)
    implicit none
    real*8,intent(in) :: x , q , I
    real*8,parameter :: g = 9.81d0
    real*8 :: fg
    fg = q/x + 2*sqrt(g*x) - I
  end function fg

  function fg_prime(x ,  q )
    implicit none
    real*8,intent(in) :: x , q 
    real*8,parameter :: g = 9.81d0
    real*8 :: fg_prime
    fg_prime = -q/(x**2) + sqrt(g/x)

  end function fg_prime


  function fd(x ,  q , I)
    implicit none
    real*8,intent(in) :: x , q , I
    real*8,parameter :: g = 9.81d0
    real*8 :: fd
    fd = q/x - 2*sqrt(g*x) - I
  end function fd

  function fd_prime(x ,  q )
    implicit none
    real*8,intent(in) :: x , q 
    real*8,parameter :: g = 9.81d0
    real*8 :: fd_prime
    fd_prime = -q/(x**2) - sqrt(g/x)

  end function fd_prime
  
  
  
  function newton_gauche(x , q , I)
    implicit none
    real*8,intent(in) :: x , q , I
    real*8 :: newton_gauche
    real*8 :: x1 , x0
    real*8,parameter :: eps = 1.0d-06
    x0 = x
    x1 = x0 - fg(x0,q,I)/fg_prime(x,q)
    do while (abs(x1 - x0) >= eps)
       x0 = x1
       x1 = x0 - fg(x0,q,I)/fg_prime(x,q)
    end do
    newton_gauche = x1
  end function newton_gauche

  function newton_droite(x,q,I)
    implicit none
    real*8,intent(in) :: x , q ,I
    real*8 :: newton_droite
    real*8 :: x1 , x0
    real*8,parameter :: eps = 1.0d-06
    x0 = x
    x1 = x0 - fd(x0,q,I)/fd_prime(x,q)
    do while (abs(x1 - x0) >= eps)
       x0 = x1
       x1 = x0 - fd(x0,q,I)/fd_prime(x,q)
    end do
    newton_droite = x1
  end function newton_droite


  function zero(te,y,N)
    implicit none
    real*8 :: zero
    real*8,dimension(:),allocatable,intent(in) :: te , y
    integer , intent(in) :: N
    integer :: cpt = 0
    real*8 :: x0 , x1
    real*8,parameter :: eps = 1.0d-06

    x0 = 10.0d0
    !print*,"iteration =   ",cpt
    x1 = x0 - (interpolation(x0,te,y,N)/derivee_interpolation(x0,te,y,N))

    do while(abs(x1 - x0) >= eps)
       cpt = cpt + 1
       x0 = x1
       x1 = x0 - (interpolation(x0,te,y,N)/derivee_interpolation(x0,te,y,N))
       !print*,"iteration =   ",cpt
    end do
    zero = x1
  end function zero



END MODULE ModBoundary










