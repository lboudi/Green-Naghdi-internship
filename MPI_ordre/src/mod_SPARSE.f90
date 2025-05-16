MODULE modSPARSE

  USE modTOOLS
  USE modTYPES

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pre_ilut
  PUBLIC :: inv_bicgstab
  PUBLIC :: bicgstab_wp

  PUBLIC :: apldia
  PUBLIC :: diamua
  PUBLIC :: amudia
  PUBLIC :: filter
  PUBLIC :: csort
  PUBLIC :: clncsr
  PUBLIC :: transp
  PUBLIC :: chngelm
  PUBLIC :: getelm
  PUBLIC :: ivperm
  PUBLIC :: dvperm
  PUBLIC :: diapos
  PUBLIC :: aplsca
  PUBLIC :: aplsb1
  PUBLIC :: aplb1
  PUBLIC :: aplb
  PUBLIC :: csr2dns
  PUBLIC :: amub
  PUBLIC :: csr2umf
  PUBLIC :: dns2csr
  PUBLIC :: full_2_HB
  PUBLIC :: amux_
  PUBLIC :: LU_UMF
  PUBLIC :: LU_CSR
  PUBLIC :: UMF__SOLVE_FROM_NUM
  PUBLIC :: UMF__FULL_INVERSION
  PUBLIC :: UMF_INV

CONTAINS

  SUBROUTINE bicgstab_wp(in_, rhs, out_, nmax, nzmax, au, jau, ju, vv)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhs
    TYPE(csr), INTENT(IN) :: in_
    INTEGER, INTENT(IN) :: nmax
    INTEGER, INTENT(IN) :: nzmax
    REAL(KIND=8), DIMENSION(nmax), INTENT(OUT) :: out_
    INTEGER, DIMENSION(:), INTENT(IN) :: jau
    INTEGER, DIMENSION(:), INTENT(IN) :: ju
    !INTEGER, DIMENSION(:) :: iw
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: au
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: vv

    REAL(KIND=8), DIMENSION(nmax) :: sol
    REAL(kind=8), DIMENSION(16) :: fpar
    INTEGER, DIMENSION(16) :: ipar
    INTEGER :: maxits
    INTEGER :: im
    REAL(KIND=8) :: eps
    REAL(KIND=8) :: permtol
    REAL(KIND=8) :: tol
    INTEGER :: lfil
    INTEGER :: ierr
    INTEGER :: i


    ipar = 0
    fpar = 0
    !vv = 0.0D0

    im   = 10
    eps  = 1.0D-06
    maxits = 200 
    !iout = 6
    permtol = 1.0
    ipar(1) = 0
    ipar(2) = 2
    ipar(3) = 1
    ipar(4) = 20*nmax
    ipar(5) = im
    ipar(6) = maxits
    fpar(1) = eps
    fpar(2) = 2.22D-16
    tol = 0.0001D0 
    lfil = 10

    !CALL ilut (in_%nrow, in_%a, in_%ja, in_%ia, lfil, tol, au, jau, ju, 2*nzmax, vv, iw, ierr)

    DO i = 1, nmax
       sol(i) = 0.0D0 !K_tmp_2x(i)
    ENDDO

    ! ITERATIVE RESOLUTION WITH BICGSTAB, AND REVERSE COMMUNICATION (SPARSKIT)

11  CALL bcgstab(in_%nrow, rhs, sol, ipar, fpar, vv)


    IF (ipar(1) == 1) THEN
       CALL amux(in_%nrow, vv(ipar(8)), vv(ipar(9)), in_%a, in_%ja, in_%ia)
       GOTO 11
    ELSE IF (ipar(1) == 2) THEN
       CALL atmux(in_%nrow, vv(ipar(8)), vv(ipar(9)), in_%a, in_%ja, in_%ia)
       GOTO 11
    ELSE IF (ipar(1)==3 .OR. ipar(1)==5) THEN
       CALL lusol(in_%nrow, vv(ipar(8)),vv(ipar(9)),au,jau,ju)
       GOTO 11
    ELSE IF (ipar(1)==4 .OR. ipar(1)==6) THEN
       CALL lutsol(in_%nrow, vv(ipar(8)),vv(ipar(9)),au,jau,ju)
       GOTO 11
    ELSE IF (ipar(1) <= 0) THEN
       IF (ipar(1) == 0) THEN
          PRINT *, 'Iterative solver has satisfied convergence test.'
       ELSE IF (ipar(1)==-1) THEN
          PRINT *, 'Iterative solver has iterated too many times.'
          STOP
       ELSE IF (ipar(1)==-2) THEN
          PRINT *, 'Iterative solver was not given enough work space.'
          PRINT *, 'The work space should at least have ', ipar(4),' elements.'
          STOP
       ELSE IF (ipar(1)==-3) THEN
          PRINT *, 'Iterative solver is facing a break-down.'
          STOP
       ELSE
          PRINT *, 'Iterative solver terminated. code =', ipar(1)
          !STOP
       ENDIF
    ENDIF

    out_ = sol

  END SUBROUTINE bicgstab_wp

  SUBROUTINE pre_ilut(in_, rhs, out_, nmax, nzmax)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhs
    TYPE(csr), INTENT(IN) :: in_
    INTEGER, INTENT(IN) :: nmax
    INTEGER, INTENT(IN) :: nzmax
    REAL(KIND=8), DIMENSION(nmax), INTENT(OUT) :: out_
    INTEGER, DIMENSION(2*nzmax) :: jau
    INTEGER, DIMENSION(2*nzmax) :: ju
    INTEGER, DIMENSION(3*nmax) :: iw
    REAL(KIND=8), DIMENSION(2*nzmax) :: au
    REAL(KIND=8), DIMENSION(nmax*8) :: vv

    REAL(KIND=8), DIMENSION(nmax) :: sol
    REAL(kind=8), DIMENSION(16) :: fpar
    INTEGER, DIMENSION(16) :: ipar
    INTEGER :: maxits
    INTEGER :: im
    REAL(KIND=8) :: eps
    REAL(KIND=8) :: permtol
    REAL(KIND=8) :: tol
    INTEGER :: lfil
    INTEGER :: ierr
    INTEGER :: i


    ipar = 0
    fpar = 0
    vv = 0.0D0

    im   = 10
    eps  = 1.0D-06
    maxits = 200 
    !iout = 6
    permtol = 1.0
    ipar(1) = 0
    ipar(2) = 2
    ipar(3) = 1
    ipar(4) = 20*nmax
    ipar(5) = im
    ipar(6) = maxits
    fpar(1) = eps
    fpar(2) = 2.22D-16
    tol = 0.0001D0 
    lfil = 10

    CALL ilut (in_%nrow, in_%a, in_%ja, in_%ia, lfil, tol, au, jau, ju, 2*nzmax, vv, iw, ierr)


  END SUBROUTINE pre_ilut

  SUBROUTINE LU_UMF(A, numeric)

    REAL(KIND=8), DIMENSION(:,:), POINTER :: A
    TYPE(UMF) :: UMF_VAR
    INTEGER, DIMENSION(2), INTENT(OUT) :: numeric
    INTEGER :: Nb
    !INTEGER :: ierr
    !REAL(KIND = 8), DIMENSION(:), POINTER :: Ax
    !INTEGER, DIMENSION(:), POINTER :: Ai
    !INTEGER, DIMENSION(:), POINTER :: Ap

    !REAL(KIND = 8), DIMENSION(:), POINTER :: Axo
    !INTEGER, DIMENSION(:), POINTER :: Aio
    !INTEGER, DIMENSION(:), POINTER :: Apo
    TYPE(CSR) :: tt1
    Nb = SIZE(A,1)


    !CALL full_2_HB(A, UMF_VAR%Ax, UMF_VAR%Ap, UMF_VAR%Ai)


    CALL dns2csr(A, tt1)
    CALL csr2umf (tt1, UMF_VAR)

    !WRITE(6,*)'UMF_VAR%Ax',UMF_VAR%Ax
    !WRITE(6,*)'UMF_VAR%Ap',UMF_VAR%Ap
    !WRITE(6,*)'UMF_VAR%Ai',UMF_VAR%Ai


    !WRITE(6,*)'ici'
    !READ*

    CALL umf4def (UMF_VAR%control)

    UMF_VAR%control(1) = 2

    CALL umf4sym (Nb, Nb, UMF_VAR%Ap, UMF_VAR%Ai, UMF_VAR%Ax, &
         UMF_VAR%symbolic, UMF_VAR%control, UMF_VAR%info)
    IF (UMF_VAR%info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4sym: ', UMF_VAR%info(1)
       STOP
    ENDIF

    CALL umf4num (UMF_VAR%Ap, UMF_VAR%Ai, UMF_VAR%Ax, &
         UMF_VAR%symbolic, UMF_VAR%numeric, UMF_VAR%control, UMF_VAR%info)
    IF (UMF_VAR%info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4num: ', UMF_VAR%info(1)
       STOP
    ENDIF

    CALL umf4fsym (UMF_VAR%symbolic)

    numeric = UMF_VAR%numeric

    DEALLOCATE(A)
    DEALLOCATE(UMF_VAR%Ax)
    DEALLOCATE(UMF_VAR%Ai)
    DEALLOCATE(UMF_VAR%Ap)

    !DEALLOCATE(Ax)
    !DEALLOCATE(Ai)
    !DEALLOCATE(Ap)

    DEALLOCATE(tt1%a)
    DEALLOCATE(tt1%ja)
    DEALLOCATE(tt1%ia)




    IF(debug) THEN
       WRITE(6,*)'Matrice factorisée dans'!, UMF_Var
    END IF

  END SUBROUTINE LU_UMF

  SUBROUTINE LU_CSR(in_, numeric)

    TYPE(CSR), INTENT(INOUT) :: in_
    TYPE(UMF) :: UMF_VAR
    INTEGER, DIMENSION(2), INTENT(OUT) :: numeric
    INTEGER :: Nb
    !INTEGER :: ierr
    !REAL(KIND = 8), DIMENSION(:), POINTER :: Ax
    !INTEGER, DIMENSION(:), POINTER :: Ai
    !INTEGER, DIMENSION(:), POINTER :: Ap

    !REAL(KIND = 8), DIMENSION(:), POINTER :: Axo
    !INTEGER, DIMENSION(:), POINTER :: Aio
    !INTEGER, DIMENSION(:), POINTER :: Apo
    !TYPE(CSR) :: tt1

    Nb = in_%nrow


    !CALL full_2_HB(A, UMF_VAR%Ax, UMF_VAR%Ap, UMF_VAR%Ai)


    !CALL dns2csr(A, tt1)
    CALL csr2umf (in_, UMF_VAR)

    !WRITE(6,*)'UMF_VAR%Ax',UMF_VAR%Ax
    !WRITE(6,*)'UMF_VAR%Ap',UMF_VAR%Ap
    !WRITE(6,*)'UMF_VAR%Ai',UMF_VAR%Ai


    !WRITE(6,*)'ici'
    !READ*

    CALL umf4def (UMF_VAR%control)

    UMF_VAR%control(1) = 2

    CALL umf4sym (Nb, Nb, UMF_VAR%Ap, UMF_VAR%Ai, UMF_VAR%Ax, &
         UMF_VAR%symbolic, UMF_VAR%control, UMF_VAR%info)
    IF (UMF_VAR%info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4sym: ', UMF_VAR%info(1)
       STOP
    ENDIF

    CALL umf4num (UMF_VAR%Ap, UMF_VAR%Ai, UMF_VAR%Ax, &
         UMF_VAR%symbolic, UMF_VAR%numeric, UMF_VAR%control, UMF_VAR%info)
    IF (UMF_VAR%info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4num: ', UMF_VAR%info(1)
       STOP
    ENDIF

    CALL umf4fsym (UMF_VAR%symbolic)

    numeric = UMF_VAR%numeric

    !DEALLOCATE(A)
    DEALLOCATE(UMF_VAR%Ax)
    DEALLOCATE(UMF_VAR%Ai)
    DEALLOCATE(UMF_VAR%Ap)

    !DEALLOCATE(Ax)
    !DEALLOCATE(Ai)
    !DEALLOCATE(Ap)

    ! DEALLOCATE(in_%a)
    ! DEALLOCATE(in_%ja)
    ! DEALLOCATE(in_%ia)

    IF(debug) THEN
       WRITE(6,*)'Matrice factorisée dans'!, UMF_Var
    END IF

  END SUBROUTINE LU_CSR


  SUBROUTINE UMF__SOLVE_FROM_NUM(solu_st, BB, numeric)

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(INOUT):: solu_st
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: BB
    INTEGER, DIMENSION(2), INTENT(IN) :: numeric
    REAL(KIND=8),DIMENSION(20) :: control
    REAL(KIND=8),DIMENSION(90) :: info
    INTEGER :: sys

    sys = 0
    CALL umf4sol (sys, solu_st, BB, numeric, control, info)

    IF (info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4sol: ', info(1)
       STOP
    ENDIF

  END SUBROUTINE UMF__SOLVE_FROM_NUM

  SUBROUTINE UMF__FULL_INVERSION(A, BB, solu_st)

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    TYPE(UMF) :: UMF_VAR
    INTEGER :: Nb
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT):: solu_st
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: BB

    !REAL(KIND = 8), DIMENSION(:), POINTER :: Ax
    !INTEGER, DIMENSION(:), POINTER :: Ai
    !INTEGER, DIMENSION(:), POINTER :: Ap

    !REAL(KIND = 8), DIMENSION(:), POINTER :: Axo
    !INTEGER, DIMENSION(:), POINTER :: Aio
    !INTEGER, DIMENSION(:), POINTER :: Apo

    Nb = SIZE(A,1)

    CALL full_2_HB(A, UMF_VAR%Ax, UMF_VAR%Ap, UMF_VAR%Ai)



    !CALL dns2csr(A, Ax, Ai, Ap, ierr)
    !CALL csrcsc2 (nb, nb, Ax, Ai, Ap, UMF_VAR%Ax, UMF_VAR%Ai, UMF_VAR%Ap)


    ! WRITE(6,*)'Ax',Ax
    !   WRITE(6,*)'Ap',Ap
    !   WRITE(6,*)'Ai',Ai


    !  WRITE(6,*)'Axo',Axo
    !  WRITE(6,*)'Apo',Apo
    !  WRITE(6,*)'Aio',Aio

      !WRITE(6,*)'UMF_VAR%Ax in UMF_FULL',UMF_VAR%Ax
      !WRITE(6,*)'UMF_VAR%Ap in UMF_FULL',UMF_VAR%Ap
      !WRITE(6,*)'UMF_VAR%Ai in UMF_FULL',UMF_VAR%Ai
      !READ*
    !  READ*
    CALL umf4def (UMF_VAR%control)

    UMF_VAR%control(1) = 2

    CALL umf4sym (Nb, Nb, UMF_VAR%Ap, UMF_VAR%Ai, UMF_VAR%Ax, &
         UMF_VAR%symbolic, UMF_VAR%control, UMF_VAR%info)
    IF (UMF_VAR%info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4sym: ', UMF_VAR%info(1)
       STOP
    ENDIF

    CALL umf4num (UMF_VAR%Ap, UMF_VAR%Ai, UMF_VAR%Ax, &
         UMF_VAR%symbolic, UMF_VAR%numeric, UMF_VAR%control, UMF_VAR%info)
    IF (UMF_VAR%info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4num: ', UMF_VAR%info(1)
       STOP
    ENDIF

    CALL umf4fsym (UMF_VAR%symbolic)

    !DEALLOCATE(A)
    DEALLOCATE(UMF_VAR%Ax)
    DEALLOCATE(UMF_VAR%Ai)
    DEALLOCATE(UMF_VAR%Ap)

    !DEALLOCATE(Ax)
    !DEALLOCATE(Ai)
    !DEALLOCATE(Ap)

    !DEALLOCATE(Axo)
    ! DEALLOCATE(Aio)
    ! DEALLOCATE(Apo)

    IF(debug) THEN
       WRITE(6,*)'Matrice factorisée dans'
    END IF

    UMF_VAR%sys = 0
    CALL umf4sol (UMF_VAR%sys, solu_st, BB, UMF_VAR%numeric, &
         UMF_VAR%control, UMF_VAR%info)

    IF (UMF_VAR%info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4sol: ', UMF_VAR%info(1)
       STOP
    ENDIF


  END SUBROUTINE UMF__FULL_INVERSION

  SUBROUTINE UMF_INV(UMF_VAR, BB, solu_st)

    IMPLICIT NONE

    TYPE(UMF), INTENT(INOUT) :: UMF_VAR
    !REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: A
    !TYPE(UMF) :: UMF_VAR
    INTEGER :: Nb
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT):: solu_st
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: BB

    !REAL(KIND = 8), DIMENSION(:), POINTER :: Ax
    !INTEGER, DIMENSION(:), POINTER :: Ai
    !INTEGER, DIMENSION(:), POINTER :: Ap

    !REAL(KIND = 8), DIMENSION(:), POINTER :: Axo
    !INTEGER, DIMENSION(:), POINTER :: Aio
    !INTEGER, DIMENSION(:), POINTER :: Apo

    !Nb = SIZE(A,1)

    !CALL full_2_HB(A, UMF_VAR%Ax, UMF_VAR%Ap, UMF_VAR%Ai)

    Nb = UMF_VAR%nrow


    !CALL dns2csr(A, Ax, Ai, Ap, ierr)
    !CALL csrcsc2 (nb, nb, Ax, Ai, Ap, UMF_VAR%Ax, UMF_VAR%Ai, UMF_VAR%Ap)


    ! WRITE(6,*)'Ax',Ax
    !   WRITE(6,*)'Ap',Ap
    !   WRITE(6,*)'Ai',Ai


    !  WRITE(6,*)'Axo',Axo
    !  WRITE(6,*)'Apo',Apo
    !  WRITE(6,*)'Aio',Aio

    !  WRITE(6,*)'UMF_VAR%Ax',UMF_VAR%Ax
    !  WRITE(6,*)'UMF_VAR%Ap',UMF_VAR%Ap
    !  WRITE(6,*)'UMF_VAR%Ai',UMF_VAR%Ai

    !  READ*
    CALL umf4def (UMF_VAR%control)
    !WRITE(6,*)'In UMF_INV 1'
    UMF_VAR%control(1) = 2

    CALL umf4sym (Nb, Nb, UMF_VAR%Ap, UMF_VAR%Ai, UMF_VAR%Ax, &
         UMF_VAR%symbolic, UMF_VAR%control, UMF_VAR%info)
    IF (UMF_VAR%info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4sym: ', UMF_VAR%info(1)
       STOP
    ENDIF
 !WRITE(6,*)'In UMF_INV 2'
    CALL umf4num (UMF_VAR%Ap, UMF_VAR%Ai, UMF_VAR%Ax, &
         UMF_VAR%symbolic, UMF_VAR%numeric, UMF_VAR%control, UMF_VAR%info)
    IF (UMF_VAR%info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4num: ', UMF_VAR%info(1)
       STOP
    ENDIF
 !WRITE(6,*)'In UMF_INV 3'
    CALL umf4fsym (UMF_VAR%symbolic)
 !WRITE(6,*)'In UMF_INV 3'
    !DEALLOCATE(A)
    !DEALLOCATE(UMF_VAR%Ax)
    !DEALLOCATE(UMF_VAR%Ai)
    !DEALLOCATE(UMF_VAR%Ap)

    !DEALLOCATE(Ax)
    !DEALLOCATE(Ai)
    !DEALLOCATE(Ap)

    !DEALLOCATE(Axo)
    ! DEALLOCATE(Aio)
    ! DEALLOCATE(Apo)

    IF(debug) THEN
       WRITE(6,*)'Matrice factorisée dans'
    END IF
     !WRITE(6,*)'In UMF_INV 4'

    UMF_VAR%sys = 0
    CALL umf4sol (UMF_VAR%sys, solu_st, BB, UMF_VAR%numeric, &
         UMF_VAR%control, UMF_VAR%info)

    IF (UMF_VAR%info(1) < 0) THEN
       PRINT *, 'Error occurred in umf4sol: ', UMF_VAR%info(1)
       STOP
    ENDIF
 !WRITE(6,*)'In UMF_INV 5'

    CALL umf4fnum (UMF_VAR%numeric)

  END SUBROUTINE UMF_INV

  SUBROUTINE inv_bicgstab(in_, rhs, out_, nmax, nzmax)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhs
    TYPE(csr), INTENT(IN) :: in_
    INTEGER, INTENT(IN) :: nmax
    INTEGER, INTENT(IN) :: nzmax
    REAL(KIND=8), DIMENSION(nmax), INTENT(OUT) :: out_
    INTEGER, DIMENSION(2*nzmax) :: jau
    INTEGER, DIMENSION(2*nzmax) :: ju
    INTEGER, DIMENSION(3*nmax) :: iw
    REAL(KIND=8), DIMENSION(2*nzmax) :: au
    REAL(KIND=8), DIMENSION(nmax*8) :: vv

    REAL(KIND=8), DIMENSION(nmax) :: sol
    REAL(kind=8), DIMENSION(16) :: fpar
    INTEGER, DIMENSION(16) :: ipar
    INTEGER :: maxits
    INTEGER :: im
    REAL(KIND=8) :: eps
    REAL(KIND=8) :: permtol
    REAL(KIND=8) :: tol
    INTEGER :: lfil
    INTEGER :: ierr
    INTEGER :: i


    ipar = 0
    fpar = 0
    vv = 0.0D0

    im   = 10
    eps  = 1.0D-06
    maxits = 200 
    !iout = 6
    permtol = 1.0
    ipar(1) = 0
    ipar(2) = 2
    ipar(3) = 1
    ipar(4) = 20*nmax
    ipar(5) = im
    ipar(6) = maxits
    fpar(1) = eps
    fpar(2) = 2.22D-16
    tol = 0.0001D0 
    lfil = 10

    CALL ilut (in_%nrow, in_%a, in_%ja, in_%ia, lfil, tol, au, jau, ju, 2*nzmax, vv, iw, ierr)

    DO i = 1, nmax
       sol(i) = 0.0D0 !K_tmp_2x(i)
    ENDDO

    ! ITERATIVE RESOLUTION WITH BICGSTAB, AND REVERSE COMMUNICATION (SPARSKIT)

11  CALL bcgstab(in_%nrow, rhs, sol, ipar, fpar, vv)


    IF (ipar(1) == 1) THEN
       CALL amux(in_%nrow, vv(ipar(8)), vv(ipar(9)), in_%a, in_%ja, in_%ia)
       GOTO 11
    ELSE IF (ipar(1) == 2) THEN
       CALL atmux(in_%nrow, vv(ipar(8)), vv(ipar(9)), in_%a, in_%ja, in_%ia)
       GOTO 11
    ELSE IF (ipar(1)==3 .OR. ipar(1)==5) THEN
       CALL lusol(in_%nrow, vv(ipar(8)),vv(ipar(9)),au,jau,ju)
       GOTO 11
    ELSE IF (ipar(1)==4 .OR. ipar(1)==6) THEN
       CALL lutsol(in_%nrow, vv(ipar(8)),vv(ipar(9)),au,jau,ju)
       GOTO 11
    ELSE IF (ipar(1) <= 0) THEN
       IF (ipar(1) == 0) THEN
          PRINT *, 'Iterative solver has satisfied convergence test.'
       ELSE IF (ipar(1)==-1) THEN
          PRINT *, 'Iterative solver has iterated too many times.'
          STOP
       ELSE IF (ipar(1)==-2) THEN
          PRINT *, 'Iterative solver was not given enough work space.'
          PRINT *, 'The work space should at least have ', ipar(4),' elements.'
          STOP
       ELSE IF (ipar(1)==-3) THEN
          PRINT *, 'Iterative solver is facing a break-down.'
          STOP
       ELSE
          PRINT *, 'Iterative solver terminated. code =', ipar(1)
          !STOP
       ENDIF
    ENDIF

    out_ = sol

  END SUBROUTINE inv_bicgstab

  FUNCTION amux_ (x, in_) 

    INTEGER :: n
    TYPE(CSR), INTENT(IN) :: in_
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
    REAL(KIND=8), DIMENSION(in_%nrow) :: amux_
    REAL(KIND=8) :: t
    INTEGER :: i, k

    !-----------------------------------------------------------------------
    !         A times a vector
    !----------------------------------------------------------------------- 
    ! multiplies a matrix by a vector using the dot product form
    ! Matrix A is stored in compressed sparse row storage.

    ! on entry:
    !----------
    ! n     = row dimension of A
    !  x     = real array of length equal to the column dimension of
    !          the A matrix.
    !  a, ja,
    !     ia = input matrix in compressed sparse row format.
    ! 
    !  on return:
    ! -----------
    !  y     = real array of length n, containing the product y=Ax
    ! 
    ! -----------------------------------------------------------------------
    !  local variables

    n = in_%nrow

    IF(in_%ncol /= SIZE(X)) STOP 'amux: ERROR in matrix profile'

    DO i = 1,n
       !WRITE(6,*)'i in amux',i
       !     compute the inner product of row i with vector x

       t = 0.0d0
       DO  k=in_%ia(i), in_%ia(i+1)-1 
          !WRITE(6,*)'k in amux',k
          t = t + in_%a(k)*x(in_%ja(k))
       END DO

       !     store result in y(i) 
       !IF(ABS(t) > 1.0D-14 THEN
       amux_(i) = t
       !ELSE
       !   amux(i) = 0.0D0
       !END IF

    END DO
    WHERE(ABS(amux_) <=1.0D-14)
       amux_ = 0.0D0
    END WHERE

  END FUNCTION amux_

  subroutine chngelm (i, j, in_, alpha) 

    !     purpose:
    !     -------- 
    !     this function replace the element a(i,j) of a matrix a, 
    !     for any pair (i,j), by the new value alpha.  the matrix is assumed to be stored 
    !     in compressed sparse row (csr) format. it performs a
    !     binary search in the case where it is known that the elements 
    !     are sorted so that the column indices are in increasing order. 


    ! on entry: 
    !---------- 
    !     i      = the row index of the element sought (input).
    !     j      = the column index of the element sought (input).
    !     a      = the matrix a in compressed sparse row format (input).
    !     ja     = the array of column indices (input).
    !     ia     = the array of pointers to the rows' data (input).
    !     sorted = logical indicating whether the matrix is knonw to 
    !              have its column indices sorted in increasing order 
    !              (sorted=.true.) or not (sorted=.false.).
    !              (input). 
    !     alpha = the value to be replaced with
    ! on return:
    !----------- 
    !     getelm = value of a(i,j). 
    !     iadd   = address of element a(i,j) in arrays a, ja if found,
    !              zero if not found. (output) 
    !
    !     note: the inputs i and j are not checked for validity. 

    INTEGER, INTENT(IN) :: i,j
    !LOGICAL, INTENT(IN) :: sorted
    REAL(KIND=8), INTENT(IN) :: alpha
    !REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: a
    !INTEGER, DIMENSION(:), INTENT(IN) :: ia
    !INTEGER, DIMENSION(:), INTENT(IN) :: ja

    TYPE(CSR), INTENT(INOUT) :: in_

    !INTEGER, DIMENSION(SIZE(in_%ja)) :: iwk

    !     local variables.

    INTEGER :: ibeg, iend, imid, k, iadd
    REAL(KIND=8) :: getelm

    !     initialization 

    iadd = 0 
    getelm = 0.0
    ibeg = in_%ia(i)
    iend = in_%ia(i+1)-1

    !!    case where matrix is not necessarily sorted

    IF (.NOT.in_%sorted) THEN 

       ! scan the row - exit as soon as a(i,j) is found

       DO  k=ibeg, iend
          IF (in_%ja(k) ==  j) THEN
             iadd = k 
             GOTO 20 
          ENDIF
       END DO

       !    end unsorted case. begin sorted case

    ELSE

       !     begin binary search.   compute the middle index.

10     imid = ( ibeg + iend ) / 2

    !    test if  found

       IF (in_%ja(imid) == j) THEN
          iadd = imid 
          GOTO 20
       ENDIF
       IF (ibeg >= iend) GOTO 20

       !     else     update the interval bounds. 

       IF (in_%ja(imid) > j) THEN
          iend = imid -1
       ELSE 
          ibeg = imid +1
       ENDIF
       GOTO 10  

       !     end both cases

    ENDIF

20  IF (iadd /= 0) THEN 
       in_%a(iadd) = alpha
    ELSE
       STOP 'ERROR IN CHNGELM'
    END IF

  END subroutine chngelm


  FUNCTION getelm (i, j, in_) 

    !     purpose:
    !     -------- 
    !     this function returns the element a(i,j) of a matrix a, 
    !     for any pair (i,j).  the matrix is assumed to be stored 
    !     in compressed sparse row (csr) format. getelm performs a
    !     binary search in the case where it is known that the elements 
    !     are sorted so that the column indices are in increasing order. 
    !     also returns (in iadd) the address of the element a(i,j) in 
    !     arrays a and ja when the search is successsful (zero if not).

    ! on entry: 
    !---------- 
    !     i      = the row index of the element sought (input).
    !     j      = the column index of the element sought (input).
    !     a      = the matrix a in compressed sparse row format (input).
    !     ja     = the array of column indices (input).
    !     ia     = the array of pointers to the rows' data (input).
    !     sorted = logical indicating whether the matrix is knonw to 
    !              have its column indices sorted in increasing order 
    !              (sorted=.true.) or not (sorted=.false.).
    !              (input). 
    ! on return:
    !----------- 
    !     getelm = value of a(i,j). 
    !     iadd   = address of element a(i,j) in arrays a, ja if found,
    !              zero if not found. (output) 
    !
    !     note: the inputs i and j are not checked for validity. 

    INTEGER, INTENT(IN) :: i,j
    !LOGICAL, INTENT(IN) :: sorted

    !REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: a
    !INTEGER, DIMENSION(:), INTENT(IN) :: ia
    !INTEGER, DIMENSION(:), INTENT(IN) :: ja

    TYPE(CSR), INTENT(IN) :: in_

    REAL(KIND=8) :: getelm

    !INTEGER, DIMENSION(SIZE(in_%ja)) :: iwk

    !     local variables.

    INTEGER :: ibeg, iend, imid, k, iadd

    !     initialization 

    iadd = 0 
    getelm = 0.0
    ibeg = in_%ia(i)
    iend = in_%ia(i+1)-1

    !!    case where matrix is not necessarily sorted

    IF (.NOT.in_%sorted) THEN 

       ! scan the row - exit as soon as a(i,j) is found

       DO  k=ibeg, iend
          IF (in_%ja(k) ==  j) THEN
             iadd = k 
             GOTO 20 
          ENDIF
       END DO

       !    end unsorted case. begin sorted case

    ELSE

       !     begin binary search.   compute the middle index.

10     imid = ( ibeg + iend ) / 2

    !    test if  found

       IF (in_%ja(imid) == j) THEN
          iadd = imid 
          GOTO 20
       ENDIF
       IF (ibeg >= iend) GOTO 20

       !     else     update the interval bounds. 

       IF (in_%ja(imid) > j) THEN
          iend = imid -1
       ELSE 
          ibeg = imid +1
       ENDIF
       GOTO 10  

       !     end both cases

    ENDIF

20  IF (iadd /= 0) getelm = in_%a(iadd) 

  END FUNCTION getelm


  SUBROUTINE transp (nrow, ncol, in_)

    INTEGER, INTENT(IN) :: nrow
    INTEGER, INTENT(INOUT) :: ncol

    !REAL(KIND = 8), DIMENSION(:), INTENT(INOUT) :: a
    !INTEGER, DIMENSION(:), INTENT(INOUT) :: ia
    !INTEGER, DIMENSION(:), INTENT(INOUT) :: ja

    TYPE(CSR), INTENT(INOUT) :: in_

    INTEGER, DIMENSION(SIZE(in_%ja)) :: iwk

    ! local variables
    INTEGER:: i, k, j, jcol, init, nnz, ierr,l, inext
    REAL(KIND=8) :: t, t1


    ! In-place transposition routine.

    ! this subroutine transposes a matrix stored in compressed sparse row 
    ! format. the transposition is done in place in that the arrays a,ja,ia
    ! of the transpose are overwritten onto the original arrays.
    !------------------------------------------------------------------------
    ! on entry:
    !--------- 
    ! nrow	= integer. The row dimension of A.
    ! ncol	= integer. The column dimension of A.
    ! a	= real array of size nnz (number of nonzero elements in A).
    !         containing the nonzero elements 
    ! ja	= integer array of length nnz containing the column positions
    ! 	  of the corresponding elements in a.
    ! ia	= integer of size n+1, where n = max(nrow,ncol). On entry

    !
    ! on return:
    !----------
    !
    ! ncol	= actual row dimension of the transpose of the input matrix.
    !         Note that this may be .le. the input value for ncol, in
    !         case some of the last columns of the input matrix are zero
    !         columns. In the case where the actual number of rows found
    !         in transp(A) exceeds the input value of ncol, transp will
    !         return without completing the transposition. see ierr.
    ! a,
    ! ja,
    ! ia	= contains the transposed matrix in compressed sparse
    !         row format. The row dimension of a, ja, ia is now ncol.

    ! Note: 
    !----- 1) If you do not need the transposition to be done in place
    !         it is preferrable to use the conversion routine csrcsc 
    !         (see conversion routines in formats).
    !      2) the entries of the output matrix are not sorted (the column
    !         indices in each are not in increasing order) use csrcsc
    !         if you want them sorted.


    nnz = in_%ia(nrow+1)-1

    !     determine column dimension

    jcol = 0
    DO  k=1, nnz
       jcol = MAX(jcol,in_%ja(k))
    END DO
    IF (jcol > ncol) THEN
       ierr = jcol
       WRITE(6,*)'ERROR IN TRANSP'
       RETURN
    ENDIF

    !     convert to coordinate format. use iwk for row indices.

    ncol = jcol

    DO  i=1,nrow
       DO  k=in_%ia(i),in_%ia(i+1)-1
          iwk(k) = i
       END DO
    END DO
    !     find pointer array for transpose. 
    DO  i=1,ncol+1
       in_%ia(i) = 0
    END DO
    DO  k=1,nnz
       i = in_%ja(k)
       in_%ia(i+1) = in_%ia(i+1)+1
    END DO
    in_%ia(1) = 1 

    DO  i=1,ncol
       in_%ia(i+1) = in_%ia(i) + in_%ia(i+1)
    END DO

    !     loop for a cycle in chasing process. 

    init = 1
    k = 0
5   t = in_%a(init)
    i = in_%ja(init)
    j = iwk(init)
    iwk(init) = -1

6   k=k+1
    !current row number is i.  determine  where to go. 
    l = in_%ia(i)
    !     save the chased element. 
    t1 = in_%a(l)
    inext = in_%ja(l)
    !     then occupy its location.
    in_%a(l)  = t
    in_%ja(l) = j
    !     update pointer information for next element to be put in row i. 
    in_%ia(i) = l+1
    !     determine  next element to be chased
    IF (iwk(l) < 0) GOTO 65
    t = t1
    i = inext
    j = iwk(l)
    iwk(l) = -1
    IF (k < nnz) GOTO 6
    GOTO 70
65  init = init+1
    IF (init .GT. nnz) GOTO 70
    IF (iwk(init) < 0) GOTO 65
    !     restart chasing --	
    GOTO 5
70  CONTINUE
    DO  i=ncol,1,-1 
       in_%ia(i+1) = in_%ia(i)
    END DO
    in_%ia(1) = 1


    in_%sorted = .FALSE.

  END SUBROUTINE transp

  SUBROUTINE clncsr(in_)

    INTEGER :: nrow

    TYPE(CSR), INTENT(INOUT) :: in_
    INTEGER, DIMENSION(in_%nrow+1) :: iwk
    INTEGER, DIMENSION(in_%nrow) :: indu
    ! local variables
    INTEGER:: i, k, j, klast, ipos

    nrow = in_%nrow

    !     This routine cleans up a CSR matrix
    !     -- remove duplicate/zero entries,

    !on entry:
    !     nrow    -- row dimension of the matrix
    !     a,ja,ia -- input matrix in CSR format
    !

    !     On return:
    !     a,ja,ia -- cleaned matrix
    !     Work space:
    !     iwk     -- integer work space of size nrow+1
    !
    !
    !     .. eliminate duplicate entries --
    !     array INDU is used as marker for existing indices, it is also the
    !     location of the entry.
    !     IWK is used to stored the old IA array.
    !     matrix is copied to squeeze out the space taken by the duplicated
    !     entries.
    !
    DO i = 1, nrow
       indu(i) = 0
       iwk(i) = in_%ia(i)
    END DO
    iwk(nrow+1) = in_%ia(nrow+1)
    k = 1
    DO i = 1, nrow
       in_%ia(i) = k
       ipos = iwk(i)
       klast = iwk(i+1)

       DO WHILE (ipos < klast)

          j = in_%ja(ipos)
          IF (indu(j) == 0) THEN
             !     .. new entry ..
             IF (in_%a(ipos) /= 0.0D0) THEN
                indu(j) = k
                in_%ja(k) = in_%ja(ipos)
                in_%a(k) = in_%a(ipos)
                k = k + 1
             ENDIF

          ELSE
             !     .. duplicate entry ..
             in_%a(indu(j)) = in_%a(indu(j)) + in_%a(ipos)
          ENDIF
          ipos = ipos + 1

       END DO
       !     .. remove marks before working on the next row ..
       DO  ipos = in_%ia(i), k - 1
          indu(in_%ja(ipos)) = 0
       END DO
    END DO
    !WRITE(6,*)'ipos',ipos

    !    On nettoie la place gagnée en mettant des zéro (pour les vérifs de visu))
    ! 
    !    A vérifier
    in_%a(ipos:) = 0.0D0
    in_%ja(ipos:) = 0
    in_%ia(nrow+1) = k


  END SUBROUTINE clncsr

  ! SUBROUTINE clncsr(in_)

  !   !INTEGER, INTENT(IN) :: nrow

  !   TYPE(CSR), INTENT(INOUT) :: in_
  !   INTEGER, DIMENSION(in_%nrow+1) :: iwk
  !   INTEGER, DIMENSION(in_%nrow) :: indu
  !   ! local variables
  !   INTEGER:: i, k, j, kfirst, klast, ko, ipos, nrow
  !   REAL(KIND=8) :: tmp

  !   nrow = in_%nrow
  !   !     This routine cleans up a CSR matrix
  !   !     -- remove duplicate/zero entries,

  !   !on entry:
  !   !     nrow    -- row dimension of the matrix
  !   !     a,ja,ia -- input matrix in CSR format
  !   !

  !   !     On return:
  !   !     a,ja,ia -- cleaned matrix
  !   !     Work space:
  !   !     iwk     -- integer work space of size nrow+1
  !   !
  !   !
  !   !     .. eliminate duplicate entries --
  !   !     array INDU is used as marker for existing indices, it is also the
  !   !     location of the entry.
  !   !     IWK is used to stored the old IA array.
  !   !     matrix is copied to squeeze out the space taken by the duplicated
  !   !     entries.
  !   !

  !   DO i = 1, nrow
  !      indu(i) = 0
  !      iwk(i) = in_%ia(i)
  !   END DO
  !   iwk(nrow+1) = in_%ia(nrow+1)
  !   k = 1
  !   DO i = 1, nrow
  !      in_%ia(i) = k
  !      ipos = iwk(i)
  !      klast = iwk(i+1)

  !      DO WHILE (ipos < klast)

  !         j = in_%ja(ipos)
  !         IF (indu(j) == 0) THEN
  !            !     .. new entry ..
  !            IF (in_%a(ipos) /= 0.0D0) THEN
  !               indu(j) = k
  !               in_%ja(k) = in_%ja(ipos)
  !               in_%a(k) = in_%a(ipos)
  !               k = k + 1
  !            ENDIF

  !         ELSE
  !            !     .. duplicate entry ..
  !            in_%a(indu(j)) = in_%a(indu(j)) + in_%a(ipos)
  !         ENDIF
  !         ipos = ipos + 1

  !      END DO
  !      !     .. remove marks before working on the next row ..
  !      DO  ipos = in_%ia(i), k - 1
  !         indu(in_%ja(ipos)) = 0
  !      END DO
  !   END DO
  !   !WRITE(6,*)'ipos',ipos

  !   !    On nettoie la place gagnée en mettant des zéro (pour les vérifs de visu))
  !   ! 
  !   !    A vérifier
  !   in_%a(ipos:) = 0.0D0
  !   in_%ia(nrow+1) = k


  ! END SUBROUTINE clncsr

  SUBROUTINE ivperm (n, ix, perm) 

    INTEGER, INTENT(IN) :: n
    INTEGER, DIMENSION(n), INTENT(INOUT) :: perm
    INTEGER, DIMENSION(n), INTENT(INOUT) :: ix
    INTEGER :: init, ii, k, next, j
    ! this subroutine performs an in-place permutation of an integer vector 
    ! ix according to the permutation array perm(*), i.e., on return, 
    ! the vector x satisfies,
    !
    !	ix(perm(j)) :== ix(j), j=1,2,.., n
    !
    !-----------------------------------------------------------------------
    ! on entry:
    !---------
    ! n 	= length of vector x.
    ! perm 	= integer array of length n containing the permutation  array.
    ! ix	= input vector
    !
    ! on return:
    !---------- 
    ! ix	= vector x permuted according to ix(perm(*)) :=  ix(*)

    ! local variables
    INTEGER tmp, tmp1

    init      = 1
    tmp=ix(init)
    ii= perm(init)
    perm(init)= -perm(init)
    k         = 0

    ! loop

6   k = k+1

    ! save the chased element --

    tmp1= ix(ii) 
    ix(ii)= tmp
    next= perm(ii) 
    IF (next< 0 ) GOTO 65

    ! test for end 

    IF (k .GT. n) GOTO 101
    tmp       = tmp1
    perm(ii)  = - perm(ii)
    ii        = next 

    ! end loop 

    GOTO 6

    ! reinitilaize cycle --

65  init      = init+1
    IF (init > n) GOTO 101
    IF (perm(init) < 0) GOTO 65
    tmp=ix(init)
    ii= perm(init)
    perm(init)=-perm(init)
    GOTO 6

101 CONTINUE
    DO j=1, n
       perm(j) = -perm(j)
    END DO

  END SUBROUTINE ivperm


  SUBROUTINE dvperm (n, x, perm) 


    INTEGER, INTENT(IN) :: n
    INTEGER, DIMENSION(n), INTENT(INOUT) :: perm
    REAL(KIND=8), DIMENSION(n), INTENT(INOUT) :: x
    INTEGER :: init, ii, k, next, j
    REAL(KIND=8) :: tmp, tmp1
    !-----------------------------------------------------------------------
    ! this subroutine performs an in-place permutation of a real vector x 
    ! according to the permutation array perm(*), i.e., on return, 
    ! the vector x satisfies,
    !
    !	x(perm(j)) :== x(j), j=1,2,.., n
    !
    !-----------------------------------------------------------------------
    ! on entry:
    !---------
    ! n 	= length of vector x.
    ! perm 	= integer array of length n containing the permutation  array.
    ! x	= input vector
    !
    ! on return:
    !---------- 
    ! x	= vector x permuted according to x(perm(*)) :=  x(*)
    !

    ! local variables 


    init      = 1
    tmp=x(init)
    ii= perm(init)
    perm(init)= -perm(init)
    k         = 0

    ! loop

6   k = k+1

    ! save the chased element --

    tmp1= x(ii) 
    x(ii)     = tmp
    next= perm(ii) 
    IF (next < 0 ) GOTO 65

    ! test for end 

    IF (k > n) GOTO 101
    tmp       = tmp1
    perm(ii)  = - perm(ii)
    ii        = next 

    ! end loop 

    GOTO 6

    ! reinitilaize cycle --

65  init      = init+1
    IF (init > n) GOTO 101
    IF (perm(init) .LT. 0) GOTO 65
    tmp= x(init)
    ii= perm(init)
    perm(init)=-perm(init)
    GOTO 6

101 CONTINUE
    DO j=1, n
       perm(j) = -perm(j)
    END DO

  END SUBROUTINE dvperm

  SUBROUTINE csort (in_)  

    TYPE(CSR), INTENT(INOUT) :: in_
    INTEGER :: nrow, ncol
    !REAL(KIND = 8), DIMENSION(:), INTENT(INOUT) :: a
    !INTEGER, DIMENSION(:), INTENT(INOUT) :: ia
    !INTEGER, DIMENSION(:), INTENT(INOUT) :: ja

    INTEGER, DIMENSION(:), POINTER :: iwork
    ! local variables
    INTEGER:: i, k, j, m, ifirst, nnz, next, size, ko, irow

    nrow = in_%nrow
    ncol = in_%ncol

    !-----------------------------------------------------------------------
    ! This routine sorts the elements of  a matrix (stored in Compressed
    ! Sparse Row Format) in increasing order of their column indices within 
    ! each row. It uses a form of bucket sort with a cost of O(nnz) where
    ! nnz = number of nonzero elements. 
    ! requires an integer work array of length 2*nnz.  
    !-----------------------------------------------------------------------
    ! on entry:
    !--------- 
    ! n     = the row dimension of the matrix
    ! a     = the matrix A in compressed sparse row format.
    ! ja    = the array of column indices of the elements in array a.
    ! ia    = the array of pointers to the rows. 
    ! iwork = integer work array of length max ( m+1, 2*nnz ) 
    !         where m   = column dimension of the matrix
    !               nnz = (ia(n+1)-ia(1))

    ! 
    ! on return:
    !----------
    ! the matrix stored in the structure a, ja, ia is permuted in such a
    ! way that the column indices are in increasing order within each row.
    ! iwork(1:nnz) contains the permutation used  to rearrange the elements.

    !-----------------------------------------------------------------------

    ! count the number of elements in each column
    nnz = in_%ia(nrow+1)-in_%ia(1)

    size = MAX(ncol+1, 2*nnz)
    !WRITE(6,*)'size',size

    ALLOCATE(iwork(size))

    m = 0
    DO i=1, nrow
       DO k=in_%ia(i), in_%ia(i+1)-1 
          m = MAX( m, in_%ja(k) )
       END DO
    END DO

    DO  j=1,m
       iwork(j+1) = 0
    END DO
    DO i=1, nrow
       DO k=in_%ia(i), in_%ia(i+1)-1 
          j = in_%ja(k)
          iwork(j+1) = iwork(j+1)+1
       END DO
    END DO
    !
    ! compute pointers from lengths. 

    iwork(1) = 1
    DO i=1,m
       iwork(i+1) = iwork(i) + iwork(i+1)
    END DO

    ! get the positions of the nonzero elements in order of columns.

    ifirst = in_%ia(1) 
    nnz = in_%ia(nrow+1)-ifirst
    DO i=1,nrow
       DO k=in_%ia(i),in_%ia(i+1)-1 
          j = in_%ja(k) 
          next = iwork(j)
          iwork(nnz+next) = k
          iwork(j) = next+1
       END DO
    END DO

    ! convert to coordinate format

    DO  i=1, nrow
       DO  k=in_%ia(i), in_%ia(i+1)-1 
          iwork(k) = i
       END DO
    END DO
    ! loop to find permutation: for each element find the correct 
    ! position in (sorted) arrays a, ja. Record this in iwork. 

    DO  k=1, nnz
       ko = iwork(nnz+k) 
       irow = iwork(ko)
       next = in_%ia(irow)

       ! the current element should go in next position in row. iwork  records this position. 

       iwork(ko) = next
       in_%ia(irow)  = next+1
    END DO

    ! perform an in-place permutation of the  arrays 
    !========================
    ! ATTENTION BUG A REVOIR
    !CALL ivperm (nnz, in_%ja(ifirst), iwork) 
    !CALL dvperm (nnz, in_%a(ifirst), iwork) 

    CALL ivperm (nnz, in_%ja, iwork) 
    CALL dvperm (nnz, in_%a, iwork) 

    !========================

    ! reshift the pointers of the original matrix back.

    DO  i=nrow,1,-1
       in_%ia(i+1) = in_%ia(i)
    END DO
    in_%ia(1) = ifirst 

    DEALLOCATE(iwork)
    in_%sorted = .TRUE.

  END SUBROUTINE csort

  SUBROUTINE filter(job, drptol, in_, out_)

    INTEGER, INTENT(IN) :: job

    TYPE(CSR), INTENT(IN) :: in_
    TYPE(CSR), INTENT(OUT) :: out_
    REAL(KIND=8), INTENT(IN) :: drptol
    INTEGER ::  n
    REAL(KIND=8) :: norm,loctol
    INTEGER:: index,row,k,k1,k2 

    n = in_%nrow
    !-----------------------------------------------------------------------
    !      This module removes any elements whose absolute value
    !      is small from an input matrix A and puts the resulting
    !      matrix in B.  The input parameter job selects a definition
    !      of small.
    !-----------------------------------------------------------------------
    ! on entry:
    !---------
    !  n	 = integer. row dimension of matrix
    !  job   = integer. used to determine strategy chosen by caller to
    !         drop elements from matrix A. 
    !          job = 1  
    !              Elements whose absolute value is less than the
    !              drop tolerance are removed.
    !          job = 2
    !              Elements whose absolute value is less than the 
    !              product of the drop tolerance and the Euclidean
    !              norm of the row are removed. 
    !          job = 3
    !              Elements whose absolute value is less that the
    !              product of the drop tolerance and the largest
    !              element in the row are removed.
    ! 
    ! drptol = real. drop tolerance used for dropping strategy.
    ! a	
    ! ja
    ! ia     = input matrix in compressed sparse format
    ! len	 = integer. the amount of space available in arrays b and jb.
    !
    ! on return:
    !---------- 
    ! b	
    ! jb
    ! ib    = resulting matrix in compressed sparse format. 

    !------ This module is in place. (b,jb,ib can ne the same as 
    !       a, ja, ia in which case the result will be overwritten).

    ! local variables


    ALLOCATE(out_%ia(SIZE(in_%ia)))
    ALLOCATE(out_%a(SIZE(in_%a)))
    ALLOCATE(out_%ja(SIZE(in_%ja)))

    index = 1
    DO row= 1,n
       k1 = in_%ia(row)
       k2 = in_%ia(row+1) - 1
       out_%ia(row) = index
       GOTO (100,200,300) job
100    norm = 1.0d0
       GOTO 400
200    norm = 0.0d0
       DO k = k1,k2
          norm = norm + in_%a(k) * in_%a(k)
       END DO
       norm = SQRT(norm)
       GOTO 400
300    norm = 0.0d0
       DO k = k1,k2
          IF( ABS(in_%a(k)) > norm) THEN
             norm = ABS(in_%a(k))
          ENDIF
       END DO
400    loctol = drptol * norm
       DO k = k1,k2
          IF( ABS(in_%a(k)) > loctol)THEN 
             out_%a(index) =  in_%a(k)
             out_%ja(index) = in_%ja(k)
             index = index + 1
          ENDIF
       END DO
    END DO
    out_%ia(n+1) = index

    out_%sorted = in_%sorted
    out_%nrow = in_%nrow
    out_%ncol = in_%ncol

  END SUBROUTINE filter

  ! FUNCTION amux (x, in_) 

  !   INTEGER :: n
  !   TYPE(CSR), INTENT(IN) :: in_
  !   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
  !   REAL(KIND=8), DIMENSION(in_%nrow) :: amux
  !   REAL(KIND=8) :: t
  !   INTEGER :: i, k

  !   !-----------------------------------------------------------------------
  !   !         A times a vector
  !   !----------------------------------------------------------------------- 
  !   ! multiplies a matrix by a vector using the dot product form
  !   ! Matrix A is stored in compressed sparse row storage.

  !   ! on entry:
  !   !----------
  !   ! n     = row dimension of A
  !   !  x     = real array of length equal to the column dimension of
  !   !          the A matrix.
  !   !  a, ja,
  !   !     ia = input matrix in compressed sparse row format.
  !   ! 
  !   !  on return:
  !   ! -----------
  !   !  y     = real array of length n, containing the product y=Ax
  !   ! 
  !   ! -----------------------------------------------------------------------
  !   !  local variables

  !   n = in_%nrow

  !   IF(in_%ncol /= SIZE(X)) STOP 'amux: ERROR in matrix profile'

  !   DO i = 1,n

  !      !     compute the inner product of row i with vector x

  !      t = 0.0d0
  !      DO  k=in_%ia(i), in_%ia(i+1)-1 
  !         t = t + in_%a(k)*x(in_%ja(k))
  !      END DO

  !      !     store result in y(i) 
  !      amux(i) = t
  !   END DO

  ! END FUNCTION amux


  SUBROUTINE amudia (in_, diag, out_)

    INTEGER :: nrow
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: diag

    TYPE(CSR), INTENT(IN) :: in_
    TYPE(CSR), INTENT(OUT) :: out_
    !INTEGER, DIMENSION(in_%nrow) :: iw
    INTEGER :: k, ii, k1, k2


    !-----------------------------------------------------------------------
    ! performs the matrix by matrix product B = A * Diag  (in place) 
    !-----------------------------------------------------------------------
    ! on entry:
    ! ---------
    ! nrow	= integer. The row dimension of A
    !
    ! job   = integer. job indicator. Job=0 means get array b only
    !         job = 1 means get b, and the integer arrays ib, jb.
    !
    ! a,
    ! ja,
    ! ia   = Matrix A in compressed sparse row format.
    ! 
    ! diag = diagonal matrix stored as a vector dig(1:n)
    !
    ! on return:
    !----------
    !
    ! b, 
    ! jb, 
    ! ib	= resulting matrix B in compressed sparse row sparse format.

    nrow = in_%nrow

    IF(in_%ncol /= size(diag)) STOP "amudia : error in matrix profiles "

    !WRITE(6,*)'nrow in amudia', nrow

    ALLOCATE(out_%ia(SIZE(in_%ia)))
    ALLOCATE(out_%a(SIZE(in_%a)))
    ALLOCATE(out_%ja(SIZE(in_%ja)))

    DO  ii=1, nrow

       !     scale each element 

       k1 = in_%ia(ii)
       k2 = in_%ia(ii+1)-1
       !WRITE(6,*)'k1',k1,'k2',k2
       DO  k=k1, k2
          out_%a(k) = in_%a(k)*diag(in_%ja(k)) 
       END DO
    END DO

    DO  ii=1, nrow+1
       out_%ia(ii) = in_%ia(ii)
    END DO

    DO k=in_%ia(1), in_%ia(nrow+1) -1 
       out_%ja(k) = in_%ja(k)
    END DO

    out_%nrow = in_%nrow
    out_%ncol = in_%ncol
    out_%sorted = in_%sorted

  END SUBROUTINE amudia

  SUBROUTINE diamua (in_, diag, out_)

    TYPE(CSR), INTENT(IN) :: in_
    !TYPE(CSR), INTENT(INOUT) :: out_
    TYPE(CSR), INTENT(INOUT) :: out_
    INTEGER :: nrow
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: diag

    !INTEGER, DIMENSION(in_%nrow) :: iw
    INTEGER :: k, ii, k1, k2
    REAL(KIND=8) :: scal

    nrow = in_%nrow

    IF(in_%nrow /= size(diag)) STOP "diamua : error in matrix profiles "

    ! performs the matrix by matrix product B = Diag * A  (in place) 

    ! on entry:
    ! ---------
    ! nrow == integer. The row dimension of A

    ! a,
    ! ja,
    ! ia   = Matrix A in compressed sparse row format.
    ! 
    ! diag = diagonal matrix stored as a vector dig(1:n)
    !
    ! on return:
    !----------
    !
    ! b, 
    ! jb, 
    ! ib	= resulting matrix B in compressed sparse row sparse format.
    !	    
    ! Notes:
    !-------
    ! 1)        The column dimension of A is not needed. 
    ! 2)        algorithm in place (B can take the place of A).
    !           in this case use job=0.
    !-----------------------------------------------------------------
    !out_%a =0.0D0
    !out_%ia = 0
    !out_%ja = 0

    !ALLOCATE(out_%ia(SIZE(in_%ia)))
    !ALLOCATE(out_%a(SIZE(in_%a)))
    !ALLOCATE(out_%ja(SIZE(in_%ja)))
    !WRITE(6,*)'size(in_%a)',size(in_%a)
    !WRITE(6,*)'size(out_%a)',size(out_%a)
    !read*
    !WRITE(6,*)'size stmp1%a', size(out_%a)
    !WRITE(6,*)'size in%a', size(in_%a)

    DO  ii=1,nrow

       !     normalize each row 

       k1 = in_%ia(ii)
       k2 = in_%ia(ii+1)-1
       scal = diag(ii) 
       DO  k=k1, k2
          out_%a(k) = in_%a(k)*scal
       END DO
    END DO

    DO  ii=1, nrow+1
       out_%ia(ii) = in_%ia(ii)
    END DO
    DO  k=in_%ia(1), in_%ia(nrow+1) -1 
       out_%ja(k) = in_%ja(k)
    END DO

    out_%nrow = in_%nrow
    out_%ncol = in_%ncol
    out_%sorted = in_%sorted

  END SUBROUTINE diamua

  SUBROUTINE apldia (in_, diag, out_) 

    INTEGER :: nrow

    TYPE(CSR), INTENT(IN) :: in_
    TYPE(CSR), INTENT(INOUT) :: out_

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: diag
    !REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: a
    !INTEGER, DIMENSION(:), INTENT(IN) :: ia
    !INTEGER, DIMENSION(:), INTENT(IN) :: ja

    !REAL(KIND = 8), DIMENSION(SIZE(in_%a)+in_%nrow) :: tmp_c
    !INTEGER, DIMENSION(SIZE(in_%a)+in_%nrow) :: tmp_jc

    !REAL(KIND = 8), DIMENSION(:), POINTER :: c
    !INTEGER, DIMENSION(:), POINTER :: jc

    !INTEGER, DIMENSION(:), POINTER :: ic


    INTEGER, DIMENSION(in_%nrow) :: iw
    INTEGER :: j, k, ii,icount, ko, k1, k2, nnz
    LOGICAL :: test

    nrow = in_%nrow

    !WRITE(6,*)'nrow in apldia',nrow

    IF(in_%nrow /= size(diag)) STOP 'apldia: ERROR IN MATRIX PROFILE'
    IF(in_%ncol /= size(diag)) STOP 'apldia: ERROR IN MATRIX PROFILE'

    !-----------------------------------------------------------------------
    ! Adds a diagonal matrix to a general sparse matrix:  B = A + Diag 
    !-----------------------------------------------------------------------
    ! on entry:
    ! ---------
    ! nrow	= integer. The row dimension of A

    ! a,
    ! ja,
    ! ia   = Matrix A in compressed sparse row format.
    !     
    ! diag = diagonal matrix stored as a vector dig(1:n)
    !
    ! on return:
    !----------
    !
    ! b, 
    ! jb, 
    ! ib	= resulting matrix B in compressed sparse row sparse format.
    !
    !
    ! iw    = integer work array of length n. On return iw will
    !         contain  the positions of the diagonal entries in the 
    !         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
    !         are the values/column indices of the diagonal elements 
    !         of the output matrix. ). 



    !     copy integer arrays into b's data structure if required

    !ALLOCATE(out_%ia(SIZE(in_%ia)))
    !ALLOCATE(out_%a(SIZE(in_%a)+nrow))
    !ALLOCATE(out_%ja(SIZE(in_%ja)+nrow))

    !WRITE(6,*)'out_%ia',SIZE(out_%ia)
    !WRITE(6,*)'out_%a',SIZE(out_%a)
    !WRITE(6,*)'out_%ja',SIZE(out_%ja)

    nnz = in_%ia(nrow+1)-1
    DO   k=1, nnz
       out_%ja(k) = in_%ja(k)
       out_%a(k)  = in_%a(k) 
    END DO

    out_%ia = in_%ia


    !     get positions of diagonal elements in data structure.

    CALL diapos (nrow,in_%ja,in_%ia,iw)

    !     count number of holes in diagonal and add diag(*) elements to
    !     valid diagonal entries.

    icount = 0
    DO j=1, nrow
       IF (iw(j)== 0) THEN
          icount = icount+1
       ELSE
          out_%a(iw(j)) = in_%a(iw(j)) + diag(j) 
       ENDIF
    END DO

    !     if no diagonal elements to insert return

    IF (icount == 0) THEN
       out_%nrow = in_%nrow
       out_%ncol = in_%ncol
       out_%sorted = in_%sorted
       RETURN
    END IF

    !     shift the nonzero elements if needed, to allow for created 
    !     diagonal elements. 

    ko = out_%ia(nrow+1)+icount

    !     copy rows backward

    DO  ii=nrow, 1, -1 

       !     go through  row ii

       k1 = out_%ia(ii)
       k2 = out_%ia(ii+1)-1 
       out_%ia(ii+1) = ko
       test = (iw(ii) == 0) 
       DO  k = k2,k1,-1 
          j = out_%ja(k)
          IF (test .AND. (j < ii)) THEN 
             test = .FALSE. 
             ko = ko - 1
             out_%a(ko) = diag(ii) 
             out_%ja(ko) = ii
             iw(ii) = ko
          ENDIF
          ko = ko-1
          out_%a(ko) = out_%a(k) 
          out_%ja(ko) = j
       END DO
       !     diagonal element has not been added yet.
       IF (test) THEN
          ko = ko-1
          out_%a(ko) =  diag(ii) 
          out_%ja(ko) = ii
          iw(ii) = ko
       ENDIF
    END DO
    out_%ia(1) = ko 

    out_%nrow = in_%nrow
    out_%ncol = in_%ncol
    out_%sorted = in_%sorted

    !WRITE(6,*)'out_%nrow in apldia',out_%nrow
    !WRITE(6,*)'out_%ncol in apldia',out_%ncol


  END SUBROUTINE apldia

  SUBROUTINE diapos(n,ja,ia,idiag) 

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER, DIMENSION(:), INTENT(in) ::ia
    INTEGER, DIMENSION(:), INTENT(in) ::ja
    INTEGER, DIMENSION(n), INTENT(OUT) :: idiag
    INTEGER :: i,k
    !  this subroutine returns the positions of the diagonal elements of a
    !  sparse matrix a, ja, ia, in the array idiag.

    ! on entry:
    !  -----------
    !  n	= integer. row dimension of the matrix a.
    !  a,ja,ia = matrix stored compressed sparse row format. a array skipped.

    !   on return:
    !  -----------
    !    idiag  = integer array of length n. The i-th entry of idiag 
    !             points to the diagonal element a(i,i) in the arrays
    !             a, ja. (i.e., a(idiag(i)) = element A(i,i) of matrix A)
    !             if no diagonal element is found the entry is set to 0.

    idiag = 0

    !     sweep through data structure. 

    DO   i=1,n
       DO  k= ia(i),ia(i+1) -1
          IF (ja(k)== i) idiag(i) = k
       END DO
    END DO

  END SUBROUTINE diapos



  SUBROUTINE aplsca (in_, scal, out_)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: scal
    INTEGER :: nrow
    TYPE(CSR), INTENT(IN) :: in_
    TYPE(CSR), INTENT(INOUT) :: out_

    INTEGER, DIMENSION(in_%nrow) :: iw
    INTEGER :: j, k, ii, icount, ko, k1, k2
    LOGICAL :: test

    nrow = in_%nrow

    ! Adds a scalar to the diagonal entries of a sparse matrix A :=A + s I 

    ! on entry:
    ! ---------
    ! nrow The row dimension of A
    !
    ! a,
    ! ja,
    ! ia    = Matrix A in compressed sparse row format.
    ! 
    ! scal  = real. scalar to add to the diagonal entries. 
    !
    ! on return:
    !----------
    !
    ! c, 
    ! jc, 
    ! ic	= matrix A with diagonal elements shifted (or created).
    !	    
    ! iw    = integer work array of length n. On return iw will
    !         contain  the positions of the diagonal entries in the 
    !         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
    !         are the values/column indices of the diagonal elements 
    !         of the output matrix. ). 

    !ALLOCATE(out_%ia(SIZE(in_%ia)))
    !ALLOCATE(out_%a(SIZE(in_%a)+nrow))
    !ALLOCATE(out_%ja(SIZE(in_%ja)+nrow))

    out_%a = 0.0D0
    out_%ja = 0
    out_%a(1:SIZE(in_%a)) = in_%a
    out_%ja(1:SIZE(in_%ja)) = in_%ja
    out_%ia = in_%ia

    iw = 0

    CALL diapos (nrow,in_%ja,in_%ia,iw)
    !WRITE(6,*)'iw dans iplsca',iw
    icount = 0
    DO  j=1, nrow
       IF (iw(j)== 0) THEN
          icount = icount+1
       ELSE
          out_%a(iw(j)) = out_%a(iw(j)) + scal 
       ENDIF
    END DO

    !     if no diagonal elements to insert in data structure return.

    IF (icount== 0) THEN
       out_%nrow = in_%nrow
       out_%ncol = in_%ncol
       out_%sorted = in_%sorted

       RETURN
    END IF
    ! shift the nonzero elements if needed, to allow for created 
    ! diagonal elements. 

    ko = in_%ia(nrow+1)+icount

    !     copy rows backward

    DO  ii=nrow, 1, -1 

       !     go through  row ii

       k1 = out_%ia(ii)


       k2 = out_%ia(ii+1)-1 

       out_%ia(ii+1) = ko
       test = (iw(ii) == 0) 

       DO  k = k2,k1,-1 
          j = out_%ja(k)
          IF (test .AND. (j < ii)) THEN 
             test = .FALSE. 
             ko = ko - 1
             out_%a(ko) = scal 
             out_%ja(ko) = ii
             iw(ii) = ko
          ENDIF
          ko = ko-1
          out_%a(ko) = out_%a(k) 
          out_%ja(ko) = j
       END DO
       !     diagonal element has not been added yet.
       IF (test) THEN
          ko = ko-1
          out_%a(ko) = scal 
          out_%ja(ko) = ii
          iw(ii) = ko
          !WRITE(6,*)'c(',ko,')',c(ko)

       ENDIF
    END DO
    out_%ia(1) = ko 
    !write(6,*)'c end of aplsca',c
    out_%nrow = in_%nrow
    out_%ncol = in_%ncol
    out_%sorted = in_%sorted

  END SUBROUTINE aplsca


  SUBROUTINE aplsb1 (in1, s, in2, out_)

    REAL(KIND=8), INTENT(IN) :: s
    TYPE(CSR), INTENT(IN) :: in1
    TYPE(CSR), INTENT(IN) :: in2
    TYPE(CSR), INTENT(INOUT) :: out_

    !REAL(KIND = 8), DIMENSION(SIZE(in1%a)+SIZE(in2%a)) :: tmp_c
    !INTEGER, DIMENSION(SIZE(in1%a)+SIZE(in2%a)) :: tmp_jc

    !INTEGER, DIMENSION(in1%ncol) :: iw
    INTEGER :: i, kb, ka, len, kc,kamax, kbmax, j1, j2, nrow, ncol



    ! performs the operation C = A+s B for matrices in sorted CSR format.
    ! the difference with aplsb is that the resulting matrix is such that
    ! the elements of each row are sorted with increasing column indices in
    ! each row, provided the original matrices are sorted in the same way. 
    !-----------------------------------------------------------------------
    ! on entry:
    ! ---------
    ! nrow	= integer. The row dimension of A and B
    ! ncol  = integer. The column dimension of A and B.
    !
    ! a,
    ! ja,
    ! ia   = Matrix A in compressed sparse row format with entries sorted
    !
    ! s	= real. scalar factor for B.
    ! 
    ! b, 
    ! jb, 
    ! ib	=  Matrix B in compressed sparse row format with entries sorted
    !        ascendly in each row   

    ! on return:
    !----------
    ! c, 
    ! jc, 
    ! ic	= resulting matrix C in compressed sparse row sparse format
    !         with entries sorted ascendly in each row

    ! Notes: 
    !-------
    !     this will not work if any of the two input matrices is not sorted
    !-----------------------------------------------------------------------


    !ierr = 0

    IF(in1%nrow /= in2%nrow) STOP "aplb : error in matrix profiles "
    IF(in1%ncol /= in2%ncol) STOP "aplb : error in matrix profiles "

    nrow = in1%nrow
    ncol = in1%ncol

    !ALLOCATE(out_%ia(SIZE(in1%ia)))

    len = 0
    kc = 1
    out_%ia(1) = kc 

    !     the following loop does a merge of two sparse rows + adds  them.

    DO  i=1, nrow
       ka = in1%ia(i)
       kb = in2%ia(i)
       kamax = in1%ia(i+1)-1
       kbmax = in2%ia(i+1)-1 

       DO WHILE(ka <= kamax .OR. kb <= kbmax)

          !5      CONTINUE 

          !     this is a while  -- do loop -- 

          !IF (ka <= kamax .OR. kb <= kbmax) THEN      
          IF (ka <= kamax) THEN
             j1 = in1%ja(ka)
          ELSE
             !     take j1 large enough  that always j2 .lt. j1
             j1 = ncol+1
          ENDIF
          IF (kb <= kbmax) THEN 
             j2 = in2%ja(kb)         
          ELSE 
             !     similarly take j2 large enough  that always j1 .lt. j2 
             j2 = ncol+1
          ENDIF

          !     three cases

          IF (j1== j2) THEN 
             out_%a(kc) = in1%a(ka)+s*in2%a(kb)
             out_%ja(kc) = j1
             ka = ka+1
             kb = kb+1
             kc = kc+1
             len = len+1
          ELSE IF (j1 < j2) THEN
             out_%ja(kc) = j1
             out_%a(kc) = in1%a(ka)
             ka = ka+1
             kc = kc+1
             len = len+1
          ELSE IF (j1 > j2) THEN
             out_%ja(kc) = j2
             out_%a(kc) = s*in2%a(kb)
             kb = kb+1
             kc = kc+1
             len = len+1
          ENDIF

          !GOTO 5

          !     end while loop
          !ENDIF
       END DO

       out_%ia(i+1) = kc
    END DO

    out_%sorted = .TRUE.

    out_%nrow = in1%nrow
    out_%ncol = in1%ncol

  END SUBROUTINE aplsb1


  SUBROUTINE aplb1(in1, in2, out_)

    IMPLICIT NONE

    !INTEGER, INTENT(IN) :: nrow
    !INTEGER, INTENT(IN) :: ncol

    TYPE(CSR), INTENT(IN) :: in1
    TYPE(CSR), INTENT(IN) :: in2
    TYPE(CSR), INTENT(OUT) :: out_
    INTEGER :: nrow
    INTEGER :: ncol
    REAL(KIND = 8), DIMENSION(SIZE(in1%a)+SIZE(in2%a)) :: tmp_c
    INTEGER, DIMENSION(SIZE(in1%a)+SIZE(in2%a)) :: tmp_jc

    !INTEGER, DIMENSION(in1%ncol) :: iw
    INTEGER :: i, kb, ka, len, kc,kamax, kbmax, j1, j2


    IF(in1%nrow /= in2%nrow) STOP "aplb : error in matrix profiles "
    IF(in1%ncol /= in2%ncol) STOP "aplb : error in matrix profiles "

    ALLOCATE(out_%ia(SIZE(in1%ia))) ! nrow +1

    !-----------------------------------------------------------------------
    ! performs the matrix sum  C = A+B for matrices in sorted CSR format.
    ! the difference with aplb  is that the resulting matrix is such that
    ! the elements of each row are sorted with increasing column indices in
    ! each row, provided the original matrices are sorted in the same way. 
    !-----------------------------------------------------------------------
    ! on entry:
    ! ---------
    ! nrow	= integer. The row dimension of A and B
    ! ncol  = integer. The column dimension of A and B.

    !
    ! a,
    ! in1%ja,
    ! in1%ia   = Matrix A in compressed sparse row format with entries sorted
    ! 
    ! b, 
    ! in2%ja, 
    ! in2%ia	=  Matrix B in compressed sparse row format with entries sorted
    !        ascendly in each row   
    ! 
    ! on return:
    !----------
    ! c, 
    ! jc, 
    ! ic	= resulting matrix C in compressed sparse row sparse format
    !         with entries sorted ascendly in each row. 

    ! Notes: 
    !-------
    !     this will not work if any of the two input matrices is not sorted
    !-----------------------------------------------------------------------

    nrow = in1%nrow
    ncol = in1%ncol

    kc = 1
    out_%ia(1) = kc 
    len = 0

    DO i=1, nrow
       ka = in1%ia(i)
       kb = in2%ia(i)
       kamax = in1%ia(i+1)-1
       kbmax = in2%ia(i+1)-1 

       DO WHILE(ka .LE. kamax .OR. kb .LE. kbmax)

          !5      CONTINUE 
          IF (ka .LE. kamax) THEN
             j1 = in1%ja(ka)
          ELSE
             j1 = ncol+1
          ENDIF
          IF (kb .LE. kbmax) THEN 
             j2 = in2%ja(kb)         
          ELSE 
             j2 = ncol+1
          ENDIF

          !     three cases

          IF (j1 .EQ. j2) THEN 
             tmp_c(kc) = in1%a(ka)+in2%a(kb)
             tmp_jc(kc) = j1
             ka = ka+1
             kb = kb+1
             kc = kc+1
             len = len+1
          ELSE IF (j1 .LT. j2) THEN
             tmp_jc(kc) = j1
             tmp_c(kc) = in1%a(ka)
             ka = ka+1
             kc = kc+1
             len = len+1
          ELSE IF (j1 .GT. j2) THEN
             tmp_jc(kc) = j2
             tmp_c(kc) = in2%a(kb)
             kb = kb+1
             kc = kc+1
             len = len+1
          ENDIF

          !IF (ka .LE. kamax .OR. kb .LE. kbmax) GOTO 5
          out_%ia(i+1) = kc
       END DO
    END DO
    !WRITE(6,*)'tmp_c in aplb1',(tmp_c(i),i=1,SIZE(tmp_c))
    !WRITE(6,*)'len en fin d aplb1',len
    !READ*

    ALLOCATE(out_%a(1:len))
    ALLOCATE(out_%ja(1:len))
    out_%a(1:len) = tmp_c(1:len)
    out_%ja(1:len) = tmp_jc(1:len)

    out_%sorted = .TRUE.

    out_%nrow = in1%nrow
    out_%ncol = in1%ncol

  END SUBROUTINE aplb1


  SUBROUTINE aplb (in1, in2, out_)

    IMPLICIT NONE

    TYPE(CSR), INTENT(IN) :: in1
    TYPE(CSR), INTENT(IN) :: in2
    TYPE(CSR), INTENT(OUT) :: out_


    REAL(KIND = 8), DIMENSION(SIZE(in1%a)+SIZE(in2%a)) :: tmp_c
    INTEGER, DIMENSION(SIZE(in1%a)+SIZE(in2%a)) :: tmp_jc


    INTEGER, DIMENSION(in1%ncol) :: iw
    INTEGER :: k, ii, kb, ka, len, jcol, jpos

    ALLOCATE(out_%ia(SIZE(in1%ia))) ! nrow +1

    !-----------------------------------------------------------------------
    ! performs the matrix sum  C = A+B. 
    !-----------------------------------------------------------------------
    ! on entry:
    ! ---------
    ! nrow	= integer. The row dimension of A and B
    ! ncol  = integer. The column dimension of A and B.

    !
    ! a,
    ! ja,
    ! ia   = Matrix A in compressed sparse row format.
    ! 
    ! b, 
    ! jb, 
    ! ib	=  Matrix B in compressed sparse row format.

    ! 
    ! on return:
    !----------
    ! c, 
    ! jc, 
    ! ic	= resulting matrix C in compressed sparse row sparse format.
    !	    
    !
    ! work arrays:
    !------------
    ! iw	= integer work array of length equal to the number of
    !         columns in A.
    !
    !-----------------------------------------------------------------------


    ! CHECK MATRIX PROFILES
    IF(in1%nrow /= in2%nrow) STOP "aplb : error in matrix profiles "
    IF(in1%ncol /= in2%ncol) STOP "aplb : error in matrix profiles "

    len = 0
    out_%ia(1) = 1 

    iw = 0

    DO ii=1, in1%nrow
       !     row i 
       !WRITE(6,*)'ii',ii

       DO ka=in1%ia(ii), in1%ia(ii+1)-1 

          !WRITE(6,*)'ka',ka
          len = len+1
          jcol    = in1%ja(ka)

          tmp_jc(len) = jcol 
          tmp_c(len)  = in1%a(ka) 
          iw(jcol)= len
       END DO

       DO kb=in2%ia(ii),in2%ia(ii+1)-1
          !WRITE(6,*)'kb',kb
          jcol = in2%ja(kb)
          jpos = iw(jcol)
          IF (jpos == 0) THEN
             len = len+1

             tmp_jc(len) = jcol
             tmp_c(len)  = in2%a(kb)
             iw(jcol)= len
          ELSE
             tmp_c(jpos) = tmp_c(jpos) + in2%a(kb)
          ENDIF
       END DO

       DO  k=out_%ia(ii), len
          !WRITE(6,*)'k',k
          iw(tmp_jc(k)) = 0
       END DO
       !WRITE(6,*)'here'
       out_%ia(ii+1) = len+1
    END DO

    ALLOCATE(out_%a(1:len))
    ALLOCATE(out_%ja(1:len))

    out_%a(1:len) = tmp_c(1:len)
    out_%ja(1:len) = tmp_jc(1:len)

    out_%sorted = .FALSE.
    out_%nrow = in1%nrow
    out_%ncol = in1%ncol


  END SUBROUTINE aplb



  SUBROUTINE csr2dns(in_, dns) 

    IMPLICIT NONE


    TYPE(CSR), INTENT(IN) :: in_
    REAL(KIND=8), DIMENSION(in_%nrow,in_%ncol), INTENT(out) :: dns
    INTEGER :: ierr
    INTEGER :: i,j,k
    INTEGER :: nrow, ncol

    nrow = in_%nrow
    ncol = in_%ncol

    ! Compressed Sparse Row    to    Dense 

    !converts a row-stored sparse matrix into a densely stored one

    ! On entry:
    !---------- 
    !
    ! nrow	= row-dimension of a
    ! ncol	= column dimension of a
    ! a, 
    ! ja, 
    ! ia    = input matrix in compressed sparse row format. 
    !         (a=value array, ja=column array, ia=pointer array)
    ! dns   = array where to store dense matrix

    !
    ! on return: 
    !----------- 
    ! dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
    ! 
    ! ierr  = integer error indicator. 
    !         ierr .eq. 0  means normal return
    !         ierr .eq. i  means that the code has stopped when processing
    !         row number i, because it found a column number .gt. ncol.

    ierr = 0

    dns = 0.0d0

    DO  i=1,nrow
       DO  k=in_%ia(i),in_%ia(i+1)-1
          j = in_%ja(k) 
          IF (j > ncol) THEN
             ierr = i
             STOP 'ERROR in CSRDNS'
          ENDIF
          dns(i,j) = in_%a(k)
       END DO
    END DO

  END SUBROUTINE csr2dns


  SUBROUTINE amub (in1, in2, out_) 

    IMPLICIT NONE

    INTEGER :: nrow
    INTEGER :: ncol

    TYPE(CSR), INTENT(IN) :: in1
    TYPE(CSR), INTENT(IN) :: in2
    TYPE(CSR), INTENT(OUT) :: out_

    REAL(KIND = 8), DIMENSION(in1%nrow*in2%ncol) :: tmp_c
    INTEGER, DIMENSION(in1%nrow*in2%ncol) :: tmp_jc

    INTEGER, DIMENSION(in2%ncol) :: iw
    INTEGER :: k, ii, jj, kb, ka, len, jcol, jpos
    REAL(KIND=8) :: scal
    ! performs the matrix by matrix product C = A B 


    ! on entry:
    ! ---------
    ! nrow  = integer. The row dimension of A = row dimension of C
    ! ncol  = integer. The column dimension of B = column dimension of C

    !
    ! a,
    ! ja,
    ! ia   = Matrix A in compressed sparse row format.
    ! 
    ! b, 
    ! jb, 
    ! ib    =  Matrix B in compressed sparse row format
    ! 
    ! on return:
    !----------
    ! c, 
    ! jc, 
    ! ic    = resulting matrix C in compressed sparse row sparse format.
    !           
    !
    ! work arrays:
    !------------
    ! iw    = integer work array of length equal to the number of
    !         columns in A.
    ! Note: 
    !-------
    !   The row dimension of B is not needed. However there is no checking 
    !   on the condition that ncol(A) = nrow(B). 
    !
    !----------------------------------------------------------------------- 

    IF(in1%ncol /= in2%nrow) STOP "amub : error in matrix profiles for product "

    ALLOCATE(out_%ia(SIZE(in1%ia))) ! nrow +1

    !WRITE(6,*)SIZE(tmp_c)

    nrow = in1%nrow
    ncol = in2%ncol

    len = 0
    out_%ia(1) = 1 

    !     initialize array iw.

    iw = 0

    DO ii=1, nrow 
       !     row i 
       DO ka=in1%ia(ii), in1%ia(ii+1)-1 
          scal = in1%a(ka)
          jj   = in1%ja(ka)
          DO kb=in2%ia(jj),in2%ia(jj+1)-1
             jcol = in2%ja(kb)
             jpos = iw(jcol)
             IF (jpos == 0) THEN
                len = len+1

                tmp_jc(len) = jcol
                iw(jcol)= len
                tmp_c(len)  = scal*in2%a(kb)
             ELSE
                tmp_c(jpos) = tmp_c(jpos) + scal*in2%a(kb)
             ENDIF
          END DO
       END DO
       DO  k=out_%ia(ii), len
          iw(tmp_jc(k)) = 0
       END DO
       out_%ia(ii+1) = len+1
    END DO

    ALLOCATE(out_%a(1:len))
    ALLOCATE(out_%ja(1:len))
    out_%a(1:len) = tmp_c(1:len)
    out_%ja(1:len) = tmp_jc(1:len)

    out_%nrow = in1%nrow
    out_%ncol = in2%ncol
    out_%sorted = .FALSE.

  END SUBROUTINE amub

  SUBROUTINE csr2umf (in_, umf_)

    IMPLICIT NONE

    INTEGER :: nrow, ncol
    TYPE(CSR), INTENT(IN) :: in_
    TYPE(UMF), INTENT(OUT) :: umf_
    INTEGER :: ipos, i, j, k, next


    nrow = in_%nrow
    ncol = in_%ncol

    !TYPE(CSR), INTENT(OUT) :: csc

    ! ATTENTION !
    ! On renvoie un type UMF car la finalité de cette conversion est de pouvoir utiliser UMFPACK en sortie
    ! Dans le cas général il suffit de sortir un tupe CSR


    ! Compressed Sparse Row     to      Compressed Sparse Column

    ! (transposition operation)   Not in place. 

    ! Rectangular version.  n is number of rows of CSR matrix,
    !                       n2 (input) is number of columns of CSC matrix.

    !   on entry:

    !   nrow	= number of rows of CSR matrix.
    !   ncol    = number of columns of CSC matrix.

    !   a	= real array of length nnz (nnz=number of nonzero elements in input 
    !           matrix) containing the nonzero elements.
    !   ja	= integer array of length nnz containing the column positions
    !   	  of the corresponding elements in a.
    !   ia	= integer of size nrow+1. ia(k) contains the position in a, ja of
    !  	  the beginning of the k-th row.

    !   on return:

    !   output arguments:
    !   ao	= real array of size nzz containing the "a" part of the transpose
    !   jao	= integer array of size nnz containing the column indices.
    !   iao	= integer array of size ncol+1 containing the "ia" index array of
    !  	  the transpose. 

    ipos = 1

    ALLOCATE(umf_%Ap(ncol+1))
    ALLOCATE(umf_%Ai(SIZE(in_%ja)))
    ALLOCATE(umf_%Ax(SIZE(in_%a)))

    umf_%Ai = 0
    umf_%Ax = 0.0D0
    umf_%Ap = 0

    DO i=1,ncol+1
       umf_%Ap(i) = 0
       !1      continue
    END DO
    DO i=1, nrow
       DO k=in_%ia(i), in_%ia(i+1)-1 
          j = in_%ja(k)+1
          umf_%Ap(j) = umf_%Ap(j)+1
          !2            continue 
          !3            continue
       END DO
    END DO
    !---------- compute pointers from lengths ------------------------------
    umf_%Ap(1) = ipos 
    DO  i=1,ncol
       umf_%Ap(i+1) = umf_%Ap(i) + umf_%Ap(i+1)
    END DO
    !4               continue
    !--------------- now do the actual copying ----------------------------- 
    DO  i=1,nrow
       DO  k=in_%ia(i),in_%ia(i+1)-1 
          j = in_%ja(k) 
          next = umf_%Ap(j)
          umf_%Ax(next) = in_%a(k)
          umf_%Ai(next) = i
          umf_%Ap(j) = next+1
       END DO
    END DO

    !-------------------------- reshift iao and leave ---------------------- 
    DO  i=ncol,1,-1
       umf_%Ap(i+1) = umf_%Ap(i)
    END DO

    umf_%Ap(1) = ipos

    !----------------------------------------------------------------
    ! convert from 1-based to 0-based f(for UMFPACK ONLY !!!!))
    !----------------------------------------------------------------
    umf_%Ap = umf_%Ap - 1
    umf_%Ai = umf_%Ai - 1

    umf_%nrow = in_%nrow
    umf_%ncol = in_%ncol

  END SUBROUTINE csr2umf

  SUBROUTINE dns2csr(dns, csr_)

    ! -----------------------------------------------------------------------
    !  Dense		to    Compressed Row Sparse 
    !----------------------------------------------------------------------- 
    !
    !onverts a densely stored matrix into a row orientied
    !ompactly sparse matrix. ( reverse of csrdns )
    ! Note: this routine does not check whether an element 
    ! is small. It considers that a(i,j) is zero if it is exactly
    ! equal to zero: see test below.
    !-----------------------------------------------------------------------
    ! on entry:
    !---------
    !
    ! nrow	= row-dimension of a
    ! ncol	= column dimension of a
    ! nzmax = maximum number of nonzero elements allowed. This
    !         should be set to be the lengths of the arrays a and ja.
    ! dns   = input nrow x ncol (dense) matrix.
    ! ndns	= first dimension of dns. 
    !
    ! on return:
    !---------- 
    ! 
    ! a, ja, ia = value, column, pointer  arrays for output matrix 
    !
    ! ierr	= integer error indicator: 
    !         ierr .eq. 0 means normal retur
    !         ierr .eq. i means that the the code stopped while
    !         processing row number i, because there was no space left in
    !         a, and ja (as defined by parameter nzmax).
    !-----------------------------------------------------------------------



    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: dns
    TYPE(CSR), INTENT(OUT) :: csr_

    INTEGER :: ii,i
    INTEGER :: jj,j
    INTEGER :: iv
    INTEGER :: Ni, Nj
    INTEGER :: count, next

    !INTEGER, INTENT(OUT) :: ierr

    TYPE :: CELL
       REAL(KIND=8) ::  val
       INTEGER :: i
       INTEGER :: j
       INTEGER :: number
       TYPE(CELL), POINTER :: suiv
    END TYPE CELL

    TYPE( CELL ), POINTER :: The_Cell
    !TYPE( CELL ), POINTER :: Prev_cell
    TYPE( CELL ), POINTER :: Next_cell
    TYPE( CELL ), POINTER :: Debut
    TYPE( CELL ), POINTER :: New_Cell

    Ni = SIZE(dns,1)
    Nj = SIZE(dns,2)
    !print*,'Nb',Nb

    !=====================
    ! Count the number of non-zero element

    ALLOCATE(The_cell)
    Debut => The_cell
    The_Cell%suiv => NULL()
    count = 0

    DO jj = 1, Nj
       DO ii = 1, Ni

          IF(dns(ii,jj) /= 0.0D0) THEN
             count = count + 1

             The_Cell%val = dns(ii,jj)
             The_Cell%i = ii
             The_Cell%j = jj
             The_Cell%number = count

             !IF((ii==Nb).AND.(jj==Nb)) EXIT

             ALLOCATE(New_cell)
             New_Cell%suiv => NULL()
             The_Cell%suiv => New_Cell
             The_Cell => New_Cell

          END IF

       END DO
    END DO
    !PRINT*,'count',count
    ! ======================================

    ! Creation of Ax and Ai

    ALLOCATE(csr_%ia(Ni+1))
    ALLOCATE(csr_%ja(count))
    ALLOCATE(csr_%a(count))

    next = 1
    csr_%ia(1) = 1
    DO i=1, SIZE(dns,1)
       DO j=1, SIZE(dns,2) 
          IF(dns(i,j) /= 0.0D0) THEN

             csr_%ja(next) = j
             csr_%a(next) = dns(i,j)
             next = next+1
          END IF
       END DO
       csr_%ia(i+1) = next
    END DO

    ! ADDITIONNAL DATA
    csr_%sorted = .TRUE.
    csr_%nrow   = ni
    csr_%ncol   = nj

    !WRITE(6,*)'A',A
    !===== CLEAN THE HASHTABLE
    The_Cell => Debut

    DO iv = 1, count
       Next_Cell => The_Cell%suiv
       DEALLOCATE(The_Cell)
       The_Cell => Next_Cell
    END DO
    DEALLOCATE(The_Cell)

  END SUBROUTINE dns2csr


  SUBROUTINE full_2_HB(Mat, Ax, Ap, Ai)

    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: Mat
    REAL(KIND = 8), DIMENSION(:), POINTER :: Ax
    INTEGER, DIMENSION(:), POINTER :: Ai
    INTEGER, DIMENSION(:), POINTER :: Ap

    INTEGER :: ii
    INTEGER :: jj
    INTEGER :: iv
    INTEGER :: Nb
    INTEGER :: N2
    INTEGER :: count

    TYPE :: CELL
       REAL(KIND=8) ::  val
       INTEGER :: i
       INTEGER :: j
       INTEGER :: number
       TYPE(CELL), POINTER :: suiv
    END TYPE CELL

    TYPE( CELL ), POINTER :: The_Cell
    !TYPE( CELL ), POINTER :: Prev_cell
    TYPE( CELL ), POINTER :: Next_cell
    TYPE( CELL ), POINTER :: Debut
    TYPE( CELL ), POINTER :: New_Cell

    Nb = SIZE(Mat,1)
    N2 = SIZE(Mat,2)
    !print*,'Nb',Nb

    ALLOCATE(The_cell)
    Debut => The_cell
    The_Cell%suiv => NULL()
    count = 0

    DO jj = 1, Nb
       DO ii = 1, Nb

          IF(Mat(ii,jj) /= 0.0D0) THEN
             count = count + 1

             The_Cell%val = Mat(ii,jj)
             The_Cell%i = ii
             The_Cell%j = jj
             The_Cell%number = count

             !IF((ii==Nb).AND.(jj==Nb)) EXIT

             ALLOCATE(New_cell)
             New_Cell%suiv => NULL()
             The_Cell%suiv => New_Cell
             The_Cell => New_Cell

          END IF

       END DO
    END DO
    PRINT*,'count',count
    ! Creation of Ax and Ai
    ALLOCATE(Ax(count))
    ALLOCATE(Ai(count))
    ALLOCATE(Ap(Nb+1))

    The_Cell => Debut
    iv = 1

    DO WHILE (ASSOCIATED(The_Cell%suiv))
       !Do iv = 1, count
       Ax(iv) = The_Cell%val
       Ai(iv) = The_Cell%i
       The_Cell => The_Cell%suiv
       iv  =iv+1
    END DO


    ! Creation of Ap

    The_Cell => Debut
    jj=1
    iv = 1
    Ap(iv)= The_Cell%number
    The_Cell => The_Cell%suiv
    iv = iv+1
    DO WHILE (ASSOCIATED(The_Cell%suiv))

       IF(The_Cell%j /= jj) THEN
          jj = The_Cell%j
          Ap(iv) = The_Cell%number
          The_Cell => The_Cell%suiv
          iv=iv+1
       ELSE
          The_Cell => The_Cell%suiv
       END IF

    END DO
    Ap(Nb+1)= count+1

    !----------------------------------------------------------------
    !onvert from 1-based to 0-based
    !----------------------------------------------------------------
    Ap = Ap - 1
    Ai = Ai - 1

    !WRITE(UNIT = 6, FMT=*)'Ax',Ax
    !WRITE(UNIT = 6, FMT=*)'Ai',Ai
    !WRITE(UNIT = 6, FMT=*)'Ap',Ap

    ! De-allocation of the POINTER-CHAIN
    The_Cell => Debut

    DO iv = 1, count
       Next_Cell => The_Cell%suiv
       DEALLOCATE(The_Cell)
       The_Cell => Next_Cell
    END DO
    DEALLOCATE(The_Cell)


  END SUBROUTINE full_2_HB


END MODULE MODSPARSE
