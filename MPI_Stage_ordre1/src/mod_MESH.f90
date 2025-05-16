MODULE modMesh

  USE modTOOLS

  PRIVATE

  PUBLIC :: coordinates
  PUBLIC :: CREATE_MESH
  PUBLIC :: MESH_Uniform
  PUBLIC :: COMPUTE_SPACIAL_STEPS

CONTAINS

  SUBROUTINE COMPUTE_SPACIAL_STEPS(MESH, DATA)

    IMPLICIT NONE
    TYPE(MeshDef), INTENT(INOUT) :: MESH
    TYPE(DATA_IN), INTENT(IN) :: DATA
    INTEGER :: i
    !INTEGER :: j
    !INTEGER :: istat
    !CHARACTER(LEN =*), PARAMETER :: sub__name = "COMPUTE_SPACIAL_STEPS"

       DO i = -2, DATA%Nbx+3

          MESH%array_DX(i) = MESH%X_DF(i+1) - MESH%X_DF(i)
    END DO

  END SUBROUTINE COMPUTE_SPACIAL_STEPS

  SUBROUTINE COORDINATES(MESH, DATA)

    ! compute the cell-centred coordinates for output purpose

    IMPLICIT NONE
    TYPE(MeshDef), INTENT(INOUT) :: MESH
    TYPE(DATA_IN), INTENT(IN) :: DATA
    INTEGER :: i
    !INTEGER :: istat
    !CHARACTER(LEN =*), PARAMETER :: sub__name = "COORDINATES"

   DO i = -2, DATA%Nbx+3

          MESH%X_VF(i) = 0.5D0*(MESH%X_DF(i) + MESH%X_DF(i+1))

       END DO

  END SUBROUTINE COORDINATES


  SUBROUTINE MESH_UNIFORM(MESH, DATA)

    IMPLICIT NONE
    TYPE(MeshDef), INTENT(INOUT) :: MESH
    !TYPE(Bipoint), DIMENSION(:,:), POINTER :: Cell
    TYPE(DATA_IN), INTENT(IN) :: DATA
    !INTEGER :: i__stat
    REAL(kind = 8) :: Dx
    INTEGER :: i

    CHARACTER (LEN =* ), PARAMETER :: sub__name = "MESH_Uniform"
    LOGICAL :: debug

    debug = .FALSE.

    IF (debug) THEN
       !PRINT *, "Entering ", sub__name
    END IF

    Dx = ABS(DATA%X_out-DATA%X_in) / REAL(DATA%Nbx) ; !PRINT*,'Dx =',Dx 

    DO i = -2, DATA%Nbx+4

       MESH%X_DF(i) = cte__0_d_r__8

    END DO

    ! ### Maillage du domaine d'etude ###

    DO i = 1, DATA%Nbx+1

       MESH%X_DF(i) = DATA%X_in + (i-1)*Dx

    END DO


    ! ### AJOUT DE CELLULES POUR IMPOSER LES C.L. ###



    ! Ordre 1
    MESH%X_DF(0) = MESH%X_DF(1) - Dx


    ! Ordre 2
    MESH%X_DF(-1) = MESH%X_DF(0) - Dx


    ! Ordre 4/5
    MESH%X_DF(-2) = MESH%X_DF(-1) - Dx


    ! First order
    MESH%X_DF(DATA%Nbx+2) = MESH%X_DF(DATA%Nbx+1) + Dx

    ! Second order
    MESH%X_DF(DATA%Nbx+3) = MESH%X_DF(DATA%Nbx+2) + Dx


    ! Third/Fourth Order
    MESH%X_DF(DATA%Nbx+4) = MESH%X_DF(DATA%Nbx+3) + Dx


    ! DO i=-2,DATA%Nbx+4
    !    PRINT*,'MESH%X_DF dans MESH',MESH%X_DF(i),i
    ! end do
    ! read*
  END SUBROUTINE MESH_UNIFORM

  SUBROUTINE CREATE_MESH(DATA, MESH)

    IMPLICIT NONE

    TYPE(MeshDef), INTENT(OUT):: MESH
    TYPE(DATA_IN), INTENT(IN) :: DATA
    !INTEGER :: i__stat
    !REAL(kind = r__8) :: Dx

    !INTEGER :: i

    CHARACTER (LEN =* ), PARAMETER :: sub__name = "CREATE_MESH"
    LOGICAL :: debug
    !INTEGER :: istat
    debug = .FALSE.

    MESH%Nbx = DATA%Nbx
    !Nbx = DATA%Nbx


    IF (debug) THEN
       !PRINT *, "Entering ", sub__name
    END IF

       !WRITE(UNIT = 6, FMT = *)'Creating uniform cartesian mesh ...'
       CALL MESH_UNIFORM(MESH, DATA)
       CALL COMPUTE_SPACIAL_STEPS(MESH, DATA)
       !WRITE(UNIT = 6, FMT = *)'Uniform cartesian mesh created !'

       MESH%Dx = MESH%X_DF(1) - MESH%X_DF(0)
 

       IF(MESH%Dx == 0.0D0) STOP 'ERROR : Dx==0'

       MESH%array_DX = MESH%Dx

    CALL COORDINATES(MESH, DATA)

    ! DO i=-2,DATA%Nbx+4
    !    PRINT*,'MESH%X_DF dans create_mesh',MESH%X_DF(i),i
    ! end do
    ! read*

  END SUBROUTINE CREATE_MESH



END MODULE ModMesh





