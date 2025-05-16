MODULE modCommon

  !USE fortranf
  USE modCST
  USE modTypes

  IMPLICIT NONE

  !REAL(KIND=8), PARAMETER, PUBLIC :: nr = 7.5D0

  INTEGER, PARAMETER, PUBLIC :: forin = 11
  INTEGER, PARAMETER, PUBLIC :: forin2 = 12
  INTEGER, PARAMETER, PUBLIC :: forin9 = 13
  INTEGER, PARAMETER, PUBLIC :: forout = 14
  INTEGER, PARAMETER, PUBLIC :: forout2 = 15
  INTEGER, PARAMETER, PUBLIC :: forout9 = 16
  INTEGER, PARAMETER, PUBLIC :: unfin = 17
  INTEGER, PARAMETER, PUBLIC :: unfout = 18
  REAL(KIND = 8 ), PUBLIC :: r0
  REAL(KIND = 8 ), PUBLIC :: r1
  REAL(KIND = 8 ), PUBLIC :: kappa
  REAL(KIND = 8 ), PUBLIC :: Kappabis
  REAL(KIND = 8 ), PUBLIC :: Utol
  REAL(KIND = 8 ), PUBLIC :: CFL
  REAL(KIND = 8 ), PUBLIC :: Reynolds
  REAL(KIND = 8 ), PUBLIC :: Rossby
  CHARACTER(len=8),   PUBLIC :: def_bathy
  CHARACTER(LEN = 30),PUBLIC :: name_bathy
  TYPE(Data_In), PUBLIC :: DATA
  TYPE(MeshDef), PUBLIC :: MESH
  TYPE(MyRun), PUBLIC :: Run
  !TYPE(Diffusion), PUBLIC :: visc
  TYPE(VFLU), PUBLIC :: FLU
  TYPE(Variable), PUBLIC :: Var
  TYPE(T_ERR), PUBLIC :: ERROR
  CHARACTER (len=30), PARAMETER, PUBLIC :: dir_tmp = 'TMP/'
  CHARACTER (len=30), PARAMETER, PUBLIC :: dir_data = 'DATA/'
  CHARACTER (len=30), PARAMETER, PUBLIC :: dir_res = 'RESULTS/'
  CHARACTER (len=30), PARAMETER, PUBLIC :: dir_bath = 'BATHY/'

  LOGICAL, PARAMETER, PUBLIC :: debug = .FALSE.

  TYPE(coeff_RK), DIMENSION(4) :: RK

  ! INTEGER, DIMENSION(2) :: Mh_numeric
  ! INTEGER, DIMENSION(2) :: Mq_numeric
  ! TYPE(csr), SAVE :: sP
  ! TYPE(csr), SAVE :: sN

  ! TYPE(csr), SAVE :: sPh
  ! TYPE(csr), SAVE :: sPq

  ! TYPE(csr), SAVE :: sNh
  ! TYPE(csr), SAVE :: sNq
  ! TYPE(csr), SAVE :: sD1b
  ! TYPE(csr), SAVE :: sD2b
  ! TYPE(csr), SAVE :: sD3b

  ! TYPE(csr), SAVE :: op1
  ! TYPE(csr), SAVE :: op2
  ! TYPE(csr), SAVE :: op3
  ! TYPE(csr), SAVE :: op4
  ! TYPE(csr), SAVE :: op5
  ! TYPE(csr), SAVE :: op6
  ! TYPE(csr), SAVE :: opGN

  REAL(kind = 8) :: Time__1
  REAL(kind = 8) :: Time__2
  !INTEGER :: compteur
  !INTEGER :: iproc

END MODULE modCOMMON
