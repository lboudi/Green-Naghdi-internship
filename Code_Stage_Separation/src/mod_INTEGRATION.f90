MODULE modINTEGRATION

  USE modTOOLS
  USE modCST
  USE modTYPES
  USE modANALYTIC

  PRIVATE
  PUBLIC :: romberg
  PUBLIC :: Ex_from_Ref
  !PUBLIC :: romberg_sol


CONTAINS


  ! romberg integration 

  SUBROUTINE romberg(XX, S_OUT, R_OUT, FCT)

    IMPLICIT NONE

    REAL(KIND = 8),  INTENT(in) :: xx
    !TYPE( STATE ), INTENT(in) :: s_in
    TYPE( STATE ), INTENT(out) :: S_OUT
    REAL(KIND = 8),  INTENT(out), OPTIONAL :: R_OUT
    REAL(KIND = 8), DIMENSION(0:4) :: weights
    !REAL(KIND = 8), DIMENSION(1:4), INTENT(out) :: integ 
    REAL(KIND = 8), DIMENSION(0:4) :: position_x
    REAL(KIND = 8), DIMENSION(0:3) :: RES
    REAL(KIND = 8), DIMENSION(0:3) :: integ

    INTERFACE
       SUBROUTINE FCT(xx, RES)
         REAL(kIND = 8), INTENT(in) :: xx
         REAL(kIND = 8), DIMENSION(0:3), INTENT(out)  :: RES
       END SUBROUTINE FCT
    END INTERFACE

    INTEGER :: ii
    INTEGER :: kk

    weights(0) =   7.0D0/90.0D0
    weights(1) =  32.0D0/90.0D0
    weights(2) =  12.0D0/90.0D0
    weights(3) =  32.0D0/90.0D0
    weights(4) =   7.0D0/90.0D0

    position_x(0) = xx - MESH%DX/2.0D0
    position_x(1) = xx - MESH%DX/4.0D0
    position_x(2) = xx
    position_x(3) = xx + MESH%DX/4.0D0
    position_x(4) = xx + MESH%DX/2.0D0 

    integ = 0.0D0

    ! computation of function values (data) at sampling points 

    DO ii=0,4
       CALL FCT(position_x(ii), RES) 
       DO kk=0,3
          integ(kk) = integ(kk) + weights(ii) * RES(kk)
       END DO
    END DO

    S_OUT%h  = integ(0)
    S_OUT%qx = integ(1)

    IF(PRESENT(R_OUT)) THEN
       R_OUT = integ(3)
    END IF

  END SUBROUTINE romberg

  SUBROUTINE Ex_from_Ref(TabEx, TabRef, DATA)

    IMPLICIT NONE

    TYPE(STATE), DIMENSION(:), POINTER :: TabREF
    TYPE(STATE), DIMENSION(:), POINTER :: TabEx
    TYPE(DATA_IN), INTENT(IN) :: DATA
    REAL(KIND = 8) :: dx_ref
    INTEGER :: NbxRef
    INTEGER :: iref
    INTEGER :: ii
    INTEGER :: NbCell

    SELECT CASE(DATA%NCI)

    CASE(cte__16) 
       ! Accuracy test of Shu with 16384 points for the
       !   reference solution

       Nbxref = 16384
       Dx_ref = 0.00006103515625D0

       WRITE(UNIT = 6, FMT = *)'NbxRef/Nbx', NbxRef/DATA%Nbx


 
          iref = 1
          DO ii = 1, DATA%Nbx
             DO NbCell = 1, NbxRef/DATA%Nbx
                !print*,'ii',ii
                !print*,'jj',jj
                !print*,'iref',iref
                TabEx(ii)%h = TabEx(ii)%h + TabRef(iref)%h*Dx_ref
                TabEx(ii)%qx = TabEx(ii)%qx + TabRef(iref)%qx*Dx_ref
            

                iref = iref + 1

             END DO
          TabEx(ii)%h = (1.0D0/MESH%Dx)*TabEx(ii)%h
          TabEx(ii)%qx = (1.0D0/MESH%Dx)*TabEx(ii)%qx
          END DO
    

       
    END SELECT
  END SUBROUTINE Ex_from_Ref

  !   SUBROUTINE romberg_sol(xx, s_out)

  !   IMPLICIT NONE

  !     REAL(KIND = 8),  INTENT(in) :: xx
  !     !TYPE( STATE ), INTENT(in) :: s_in
  !     TYPE( STATE ), INTENT(out) :: s_out
  !     !REAL(KIND = 8),  INTENT(out) :: topo_out
  !     REAL(KIND = 8), DIMENSION(0:4) :: weights
  !     !REAL(KIND = 8), DIMENSION(1:4), INTENT(out) :: integ 
  !     REAL(KIND = 8), DIMENSION(0:4) :: position_x
  !     REAL(KIND = 8), DIMENSION(0:3) :: RES
  !     REAL(KIND = 8y), DIMENSION(0:3) :: integ

  !     INTEGER :: ii
  !     INTEGER :: kk

  !     weights(0) =   7.0D0/90.0D0
  !     weights(1) =  32.0D0/90.0D0
  !     weights(2) =  12.0D0/90.0D0
  !     weights(3) =  32.0D0/90.0D0
  !     weights(4) =   7.0D0/90.0D0

  !     position_x(0) = xx - MESH%DX/2.0D0
  !     position_x(1) = xx - MESH%DX/4.0D0	   
  !     position_x(2) = xx
  !     position_x(3) = xx + MESH%DX/4.0D0
  !     position_x(4) = xx + MESH%DX/2.0D0 

  !     integ = 0.0D0

  !     ! computation of function values (data) at sampling points 

  !     DO ii=0,4
  !        CALL ANALYTIC_SOL(position_x(ii), RES) 
  !        DO kk=0,2
  !           integ(kk) = integ(kk) + weights(ii) * RES(kk)
  !        END DO
  !     END DO

  !     s_out%h  = integ(0)
  !     s_out%qx = integ(1)
  !     s_out%qy = integ(2)
  !     !topo_out = integ(3)

  !   END SUBROUTINE romberg_sol

END MODULE modINTEGRATION
