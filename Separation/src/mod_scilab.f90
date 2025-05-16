MODULE modSCILAB
IMPLICIT NONE

CONTAINS

  PURE FUNCTION diag(diago, loc)

    ! ATTENTION, aucun garde-fou, le vecteur 'diago' doit etre à la bonne taille pour loger a la diagonale 'loc'

    !real(kind=8), dimension(:,:), intent( INOUT ) :: Mat
    REAL(KIND = 8), DIMENSION(:), INTENT( IN ) :: diago
    INTEGER, INTENT(in) :: loc
    REAL(KIND = 8), DIMENSION(SIZE(diago)+ABS(loc), SIZE(diago)+ABS(loc)) :: diag
    INTEGER :: i
    INTEGER :: istat

    diag = 0.0D0

    FORALL( i = 1: SIZE(diago))
       diag(i,i+ABS(loc)) = diago(i)
    END FORALL
    IF(loc < 0) diag = TRANSPOSE(diag)

  END FUNCTION diag

END MODULE MODSCILAB
