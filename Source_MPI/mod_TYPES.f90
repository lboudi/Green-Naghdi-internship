MODULE modTypes

  !USE fortranf

  IMPLICIT NONE

  TYPE, PUBLIC :: coeff_RK
     REAL(KIND = 8), DIMENSION(:,:), POINTER :: a
     REAL(KIND = 8), DIMENSION(:), POINTER :: b
     REAL(KIND = 8), DIMENSION(:), POINTER :: s_dt
  END TYPE coeff_RK


  TYPE, PUBLIC ::state
     REAL(KIND = 8 ) :: H    = 0.0D0
     REAL(KIND = 8 ) :: C    = 0.0D0
     REAL(KIND = 8 ) :: U    = 0.0D0
     REAL(KIND = 8 ) :: Qx   = 0.0D0
     REAL(KIND = 8 ) :: xi   = 0.0D0
     REAL(KIND = 8 ) :: topo = 0.0D0
     REAL(KIND = 8 ) :: vorti = 0.0D0
  END TYPE state
  ! ===================================================================
  ! ===================================================================
  TYPE, PUBLIC :: Variable

     !TYPE( STATE ), DIMENSION(:), POINTER :: t
     !TYPE( STATE ), DIMENSION(:), POINTER :: Df
     TYPE( STATE ), DIMENSION(:), POINTER :: Vf
     TYPE( STATE ), DIMENSION(:), POINTER :: t1
     TYPE( STATE ), DIMENSION(:), POINTER :: t2
     TYPE( STATE ), DIMENSION(:), POINTER :: t3
     TYPE( STATE ), DIMENSION(:), POINTER :: t4
     TYPE( STATE ), DIMENSION(:), POINTER :: tmp
     TYPE( STATE ), DIMENSION(:), POINTER :: t0
     !REAL(KIND=8), DIMENSION(:), POINTER :: topo
     !REAL(KIND=8), DIMENSION(:), POINTER :: topoDF
     REAL(KIND=8), DIMENSION(:), POINTER :: topoVF
     !REAL(KIND=8), DIMENSION(:), POINTER :: topoR
     ! Exact main variables : hex, hu_ex, hv_ex
     TYPE( STATE ), DIMENSION(:), POINTER :: Ex
     TYPE( STATE ), DIMENSION(:), POINTER :: Ref
     
  END TYPE Variable

  ! ===================================================================
  ! ===================================================================
  TYPE Vecteur 
     REAL(kind=8), DIMENSION(3) :: vect
  END TYPE Vecteur
  ! ===================================================================
  ! ===================================================================
  TYPE cellule
     INTEGER                :: val
     TYPE(cellule), POINTER :: suiv
  END TYPE cellule
  ! ===================================================================
  ! ===================================================================
  TYPE cellule2d !-- pour les facettes 2D (triangles) des éléments 3D
     INTEGER                  :: is2, is3 !!-- seront rangés suivant l'ordre lexicographique
     INTEGER                  :: jt1, jt2 !!-- numéro de l'élément contenant la facette; une facette appartient à 2 elts au plus
     TYPE(cellule2d), POINTER :: suiv
  END TYPE cellule2d
  ! ===================================================================
  ! ===================================================================
  TYPE rpointer
     REAL, DIMENSION(:) ,    POINTER :: Vp !!-- [[1:data%nvarphy]] Variables phys. au bord
  END TYPE rpointer
  ! ===================================================================
  ! ===================================================================
  TYPE, PUBLIC ::Bipoint
     REAL(KIND = 8 ) :: X
     REAL(KIND = 8 ) :: Y
  END TYPE Bipoint

  TYPE, PUBLIC ::BiInt
     INTEGER :: i
     INTEGER :: j
  END TYPE BiInt
  ! ===================================================================
  ! ===================================================================
  TYPE, PUBLIC :: MyRun
     INTEGER :: Niter
     REAL(KIND=8) :: Time
     REAL(KIND=8) :: UMax
     REAL(KIND=8) :: Dt
     LOGICAL :: debug
     INTEGER :: NProc
     LOGICAL :: parallelism
     INTEGER :: step
     REAL(KIND=8) :: subDt
     REAL(KIND=8) :: subTime
  END TYPE MyRun

  TYPE, PUBLIC :: T_ERR

     REAL(KIND=8) :: Niter
     REAL(KIND=8) :: time
     REAL(KIND=8) :: L2_h
     REAL(KIND=8) :: L2_qx

     REAL(KIND=8) :: L1_h
     REAL(KIND=8) :: L1_qx

     REAL(KIND=8) :: LINF_h
     REAL(KIND=8) :: LINF_qx

  END TYPE T_ERR
  ! ===================================================================
  ! ===================================================================




  TYPE, PUBLIC :: VFLU
     TYPE( STATE ), DIMENSION(:), POINTER :: X1 
     TYPE( STATE ), DIMENSION(:), POINTER :: X2

     TYPE( STATE ), DIMENSION(:), POINTER :: Sx

     REAL(KIND=8), DIMENSION(:), POINTER :: SR
     REAL(KIND=8), DIMENSION(:), POINTER :: ALPHA

  END TYPE VFLU


  ! ===================================================================
  ! ===================================================================
  TYPE, PUBLIC:: DATA_In
     ! pilote.data
     CHARACTER(LEN=80)        :: RootName 

     !!-- valeurs possibles pour Data%GeoType
     !!-- "1D"
     !!-- "2Dplan"
     !!-- "2Dper"
     !!-- "2Daxisym"
     !!-- "2DaxisymX"
     !!-- "2DaxisymY"
     !!-- "3D"
     REAL(KIND=8)             :: H0
     LOGICAL                  :: debug
     CHARACTER(LEN=10)        :: Geotype !!-- géométrie du problème
     INTEGER                  :: sel_GeoType
     CHARACTER(LEN=10)        :: MeshType !!-- type du maillage (les valeurs autorisées sont données dans donmesh)
     INTEGER                  :: NCI
     INTEGER                  :: iproc
     INTEGER                  :: limiter
     INTEGER                  :: NCLL
     INTEGER                  :: NCLR
     INTEGER                  :: NCLT
     INTEGER                  :: NCLB
     INTEGER                  :: Nbx ! Nbre de mailles, partie VF
     INTEGER                  :: Nby
     INTEGER                  :: NDf ! Nbre de noeuds, partie DF
     INTEGER                  :: in
     INTEGER                  :: out
     INTEGER                  :: ordre_tps_df
     INTEGER                  :: NtMax
     INTEGER                  :: NFreq
     INTEGER                  :: NStart
     INTEGER                  :: ChoixDom
     INTEGER                  :: Dim
     INTEGER                  :: onde
     INTEGER                  :: N1
     INTEGER                  :: N3
     INTEGER                  :: Raf
     REAL(KIND = 8 )       :: b1
     REAL(KIND = 8 )       :: DomlenX
     REAL(KIND = 8 )       :: DomlenY
     REAL(KIND = 8 )       :: TimeOut
     REAL(KIND = 8 )       :: XG
     REAL(KIND = 8 )       :: XD
     REAL(KIND = 8 )       :: X_in
     REAL(KIND = 8 )       :: X_out
     REAL(KIND = 8 )       :: YB
     REAL(KIND = 8 )       :: YH
     REAL(KIND = 8 )       :: FLW
     REAL(KIND = 8 )       :: TimeInit
     REAL(KIND = 8 )       :: x1
     REAL(KIND = 8 )       :: X2
     REAL(KIND = 8 )       :: alpha
     REAL(KIND = 8 )       :: theta
     REAL(KIND = 8 )       :: gamma
     INTEGER                  :: cond_lim
     LOGICAL                  :: SV_ONLY
     INTEGER                  :: ordre_espace_df
     INTEGER                  :: time_order

     INTEGER                  :: NbxRef
     REAL(KIND = 8 )       :: kl
     REAL(KIND = 8 )       :: kt
     REAL(KIND = 8)        :: k_manning
     REAL(KIND = 8)        :: k_darcy

     REAL(KIND = 8 )       :: kappa
     REAL(KIND = 8 )       :: Kappabis
     REAL(KIND = 8 )       :: Utol
     REAL(KIND = 8 )       :: CFL
     REAL(KIND = 8 )       :: DtFix
     INTEGER                  :: NFreqTime
     INTEGER                  :: NFreqVort
     CHARACTER(LEN=30)        :: space_order
     LOGICAL                  :: source
     CHARACTER(LEN=30)        :: source_term
     CHARACTER(LEN=30)        :: source_fric

     LOGICAL                  :: Solex
     LOGICAL                  :: ERROR
     LOGICAL                  :: ChoixDt
     LOGICAL                  :: TIME_SERIE
     LOGICAL                  :: VORTI
     LOGICAL                  :: var_plt
     LOGICAL                  :: var_vtk
     LOGICAL                  :: var_surf
     LOGICAL                  :: var_mat
     LOGICAL                  :: var_ensight
     CHARACTER(LEN=30)        :: friction


     LOGICAL                  :: sol_ref
     REAL(KIND = 8 )       :: theta_in
     REAL(KIND = 8 )       :: reynolds
     REAL(KIND = 8 )       :: mu_h
     REAL(KIND = 8 )       :: mu_v
     REAL(KIND = 8 )       :: rossby
     CHARACTER(len=8)         :: def_bathy
     CHARACTER(LEN = 30)      :: name_bathy
     CHARACTER(LEN = 30)      :: char__SolRef
     CHARACTER(LEN = 30)      :: Riemann_Solv
     CHARACTER(LEN = 30)      :: time_scheme
     CHARACTER(LEN = 30)      :: err_type
     REAL(KIND = 8 )       :: kwave
     INTEGER                  :: point_number
     TYPE(BiPoint), DIMENSION(:), POINTER :: coor
     CHARACTER(LEN = 30)      :: file_time
     CHARACTER(LEN = 30)      :: WBmeth

     LOGICAL                  :: Reprise !!-- reprise d'un job, lu sur fichier, appelle l'interface reprise_in
     CHARACTER(LEN=50)        :: Modele !!-- choix du sous-modèle
     INTEGER                  :: sel_Modele
     CHARACTER(LEN=20)        :: Approche !!-- approche mathématique utilisée
     INTEGER                  :: sel_Approche
     !!-- valeurs possibles pour Data%Methode
     !!-- "EXP-RK1"

     !!-- "EXP-RK2"

     !!-- "EXP-RK3"
     !!-- "EXP-RK4"
     !!-- "EXP-PC2"

     !!-- "IMPLIN-GEAR"
     !!-- "IMPLIN-GEAR-RD"
     !!-- "IMP-GEAR"
     !!-- "IMP-DCE"
     !!-- "IMP-GEAR-RD"
     !!-- "IMP-DCE-RD"
     !!-- "IMPLIN-RK1-RD"
     !      CHARACTER(LEN=20)        :: Methode !!-- méthode temporelle (EXPlicite ou IMPlicite)
     !      INTEGER                  :: sel_Methode
     !      !!-- valeurs possibles pour Data%LoiVisq
     !                               !!--"Constant" : EOSvisco=MUREF
     !                               !!--"Sutherland" : EOSvisco = visco_sutherland(MUREF,T)
     !                               !!--"Pant" : EOSvisco = visco_pant(T)
     !                               !!--"No" : EOSvisco = 0.0
     !      CHARACTER(len=20)        :: LoiVisq !!-- loi pour la viscosité
     !      INTEGER                  :: sel_LoiVisq
     !      INTEGER                  :: Impre !!-- niveau de mise au point, pour contrôler les impressions
     !      !!-- valeurs possibles pour Data%Iflux : euler         | eulnsvf              | multif | incomp : indifférencié

     !      INTEGER                  :: Iflux !!-- sélecteur de méthode pour calculer le flux
     !      !!-- valeurs possibles pour Data%LoiCFL
     !                               !!-- 0=clfmax ! data%cflmax
     !                               !!-- 2=minmult ! MIN(data%mult*var*kt,data%cflmax)
     !      INTEGER                  :: LoiCFL !!-- mode de calcul de la CFL variable
     !      INTEGER                  :: Ifre !!-- frequence d'affichage des resultats à l'écran
     !      CHARACTER(len=10),DIMENSION(:),POINTER  :: VarNames !!-- [[1:Data%NVarPhy]] noms des variables du modele traite : rho, u, v, w, ei ou epsi, c, t, a1, a1*rho1, a2*rho2, p, t1, t2, phi; initialisés dans <modele>/ModelInterfaces.f90
     !      INTEGER                  :: Nvar !!-- nombre de variables principales du problème : vitesse (2 ou 3) + d'autres qui dépendent du modèle utilisé
     !                               !!-- 1 ="rho"
     !                               !!-- 2 =" u "
     !                               !!-- 3 =" v "
     !                               !!-- 1+data%ndim =" w "
     !                               !!-- Data%Nvar   =" ei"
     !                               !!-- 1+Data%Nvar =" p "
     !                               !!-- 2+Data%Nvar =" C " !-- vitesse du son non normalisée (voir gp_son)
     !                               !!-- 3+Data%Nvar =" T "
     !      INTEGER                  :: NvarT !!-- nombre de variables supplémentaires, dépend du modèle
     !      INTEGER                  :: NVarPhy !!-- somme de Nvar et NvarT
     !      CHARACTER(LEN=80)        :: PbName !!-- nom du problème à résoudre (%%%%-pm : documenter les emplois)
     !      CHARACTER(LEN=80)        :: OutputName !!-- === RootName

  END TYPE DATA_In
  ! ===================================================================
  ! ===================================================================
  TYPE MeshDef
     !      !! Maillage composante volumique

     CHARACTER(LEN=30) :: mesh_type
     INTEGER :: Nbx
     REAL(KIND = 8), DIMENSION(:), POINTER :: X_DF
     LOGICAL :: refined
     REAL(KIND=8) :: dx

     REAL(KIND = 8), DIMENSION(:), POINTER :: X_VF

     REAL(KIND = 8), DIMENSION(:), POINTER :: array_DX

  END TYPE MeshDef


END MODULE ModTypes
