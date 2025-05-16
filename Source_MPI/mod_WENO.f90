MODULE modWENO

  USE modTOOLS
  USE modCST

  IMPLICIT NONE

CONTAINS
  SUBROUTINE WENO3(h_im1, h_i, h_ip1, val_l, val_r)


    ! This function computes the WENO reconstruction from cell-averages ubar
    ! in the x-direction in the three points at center - dx/2, center and center + dx/2.
    ! It uses the si_0,1,2 indicators computed by si() and negative weight splitting methode of Shi and Shu

    IMPLICIT NONE


    REAL(KIND = 8),INTENT( IN ) :: h_im1
    REAL(KIND = 8),INTENT( IN ) :: h_i
    REAL(KIND = 8),INTENT( IN ) :: h_ip1


    REAL(KIND = 8), PARAMETER :: epsilon = 1.0D-6

    REAL(KIND = 8),INTENT( OUT ) :: val_L
    REAL(KIND = 8),INTENT( OUT ) :: val_R


    REAL(KIND = 8) :: h_ip12_0
    REAL(KIND = 8) :: h_ip12_1
    REAL(KIND = 8) :: h_im12_0
    REAL(KIND = 8) :: h_im12_1
    REAL(KIND = 8) :: beta0
    REAL(KIND = 8) :: beta1
    REAL(KIND = 8) :: d0
    REAL(KIND = 8) :: d1
    REAL(KIND = 8) :: d0_tilde 
    REAL(KIND = 8) :: d1_tilde
    REAL(KIND = 8) :: omega0
    REAL(KIND = 8) :: omega1
    REAL(KIND = 8) :: omega0_tilde
    REAL(KIND = 8) :: omega1_tilde
    REAL(KIND = 8) :: alpha0
    REAL(KIND = 8) :: alpha1
    REAL(KIND = 8) :: alpha0_tilde
    REAL(KIND = 8) :: alpha1_tilde

    !INTEGER  :: k

    ! smooth indicators

    beta0 = (h_ip1 - h_i)**2
    beta1 = (h_i - h_im1)**2

    !weno3 weights

    d0 = 2.0D0/3.0D0
    d1 = 1.0D0/3.0D0
    d0_tilde = d1
    d1_tilde = d0

    alpha0 = d0/(epsilon+beta0)
    alpha1 = d1/(epsilon+beta1)

    alpha0_tilde = d0_tilde/(epsilon+beta0)
    alpha1_tilde = d1_tilde/(epsilon+beta1)

    omega0 = alpha0/(alpha0+alpha1)
    omega1 = alpha1/(alpha0+alpha1)
    omega0_tilde = alpha0_tilde/(alpha0_tilde+alpha1_tilde)
    omega1_tilde = alpha1_tilde/(alpha0_tilde+alpha1_tilde)


    ! 3rd order reconstructed values

    h_ip12_0 = 0.5D0*(h_i+h_ip1)
    h_ip12_1 = 0.5D0*(-h_im1 + 3.0D0*h_i)

    h_im12_0 = 0.5D0*(3.0D0*h_i - h_ip1)
    h_im12_1 = 0.5D0*(h_im1 + h_i)




    !//Finally, compute the reconstruction

    val_l = omega0_tilde*h_im12_0 + omega1_tilde*h_im12_1


    val_r = omega0*h_ip12_0 + omega1*h_ip12_1


  END SUBROUTINE WENO3

  SUBROUTINE Built_WENO(h_im2, h_im1, h_i, h_ip1, h_ip2, val_l, val_c, val_r)


! This function computes the WENO reconstruction from cell-averages ubar
! in the x-direction in the three points at center - dx/2, center and center + dx/2.
! It uses the si_0,1,2 indicators computed by si() and negative weight splitting methode of Shi and Shu

    IMPLICIT NONE

    REAL(KIND = 8),INTENT( IN ) :: h_im2
    REAL(KIND = 8),INTENT( IN ) :: h_im1
    REAL(KIND = 8),INTENT( IN ) :: h_i
    REAL(KIND = 8),INTENT( IN ) :: h_ip1
    REAL(KIND = 8),INTENT( IN ) :: h_ip2
    !REAL(KIND = 8) :: val_l
    !REAL(KIND = 8) :: val_c
    !REAL(KIND = 8) :: val_r
    !LOGICAL, INTENT( IN ) :: l_var
    REAL(KIND = 8), PARAMETER :: theta = 3.0D0
    REAL(KIND = 8), PARAMETER :: c13 = 1./3
    REAL(KIND = 8), PARAMETER :: c56 = 5.D0/6.D0
    REAL(KIND = 8), PARAMETER :: c16 = 1.D0/6.D0
    REAL(KIND = 8), PARAMETER :: c76 = 7.D0/6.D0
    REAL(KIND = 8), PARAMETER :: c116 = 11.D0/6.D0
    !REAL(KIND = 8), PARAMETER :: c112 = 11.D0/2.D0
    !REAL(KIND = 8), PARAMETER :: c124 = 1.D0/24.D0
    REAL(KIND = 8), PARAMETER :: eps = 1.0D-6
    !REAL(KIND = 8), PARAMETER :: c1312 =13.D0/12.D0
    !REAL(KIND = 8), PARAMETER :: c2324 =23.D0/24.D0

    REAL(KIND = 8), DIMENSION(0:2)  :: dc_p(0:2)
    REAL(KIND = 8), DIMENSION(0:2)  :: dc_m
    REAL(KIND = 8), DIMENSION(0:2)  :: dr
    REAL(KIND = 8), DIMENSION(0:2)  :: dl
    REAL(KIND = 8), DIMENSION(0:2)  :: dc
    REAL(KIND = 8), DIMENSION(0:2)  :: cl
    REAL(KIND = 8), DIMENSION(0:2)  :: cc_p
    REAL(KIND = 8), DIMENSION(0:2)  :: cc_m
    REAL(KIND = 8), DIMENSION(0:2)  :: cr
    !REAL(KIND = 8), DIMENSION(0:2)  :: cc

    REAL(KIND = 8), DIMENSION(0:2)  :: rl
    REAL(KIND = 8), DIMENSION(0:2)  :: rc
    REAL(KIND = 8), DIMENSION(0:2)  :: rr

    REAL(KIND = 8),INTENT( OUT ) :: val_L
    REAL(KIND = 8),INTENT( OUT ) :: val_R
    REAL(KIND = 8), INTENT( OUT ) :: val_C

    REAL(KIND = 8) :: sigmac_p
    REAL(KIND = 8) :: sigmac_m
    REAL(KIND = 8) :: q0
    REAL(KIND = 8) :: q1
    REAL(KIND = 8) :: q2
    REAL(KIND = 8) :: clave
    REAL(KIND = 8) :: ccave_p
    REAL(KIND = 8) :: ccave_m
    REAL(KIND = 8) :: crave
    !REAL(KIND = 8) :: ccave

    REAL(KIND = 8) :: si_0
    REAL(KIND = 8) :: si_1
    REAL(KIND = 8) :: si_2
    !REAL(KIND = 8) :: drt_top

    INTEGER  :: k

    CALL smooth_indicator(h_im2, h_im1, h_i, h_ip1, h_ip2, si_0, si_1, si_2)

    dl(0) = 0.3D0
    dl(1) = 0.6D0
    dl(2) = 0.1D0

    ! Fith order WENO weight
    ! warning !! negative weight can lead to oscilations near discontinuities
    ! the cure are 1) the negative weight splitting technique of Shu 
    !              2) ""decrease"" accuracy to use fourth order positive WENO weights

    ! First method : fifth order WENO weights
    dc(0)= -0.11250D0
    dc(1)=  1.2250D0
    dc(2)= -0.11250D0

    ! Second method : fourth order WENO weights
!    dc(0)= 0.250D0
!    dc(1)=  0.5D0
!    dc(2)= 0.250D0

    ! Third method : negative weights splitting (for central reconstruction)

    dc_p(0) = 0.50D0*(dc(0)+theta*ABS(dc(0)))
    dc_p(1) = 0.50D0*(dc(1)+theta*ABS(dc(1)))
    dc_p(2) = 0.50D0*(dc(2)+theta*ABS(dc(2)))
    dc_m(0) = dc_p(0)-dc(0)
    dc_m(1) = dc_p(1)-dc(1)
    dc_m(2) = dc_p(2)-dc(2)

    sigmac_p = dc_p(0) + dc_p(1) + dc_p(2)
    sigmac_m = dc_m(0) + dc_m(1) + dc_m(2)


    dc_p(0) = dc_p(0)/sigmac_p
    dc_p(1) = dc_p(1)/sigmac_p
    dc_p(2) = dc_p(2)/sigmac_p

    dc_m(0) = dc_m(0)/sigmac_m
    dc_m(1) = dc_m(1)/sigmac_m
    dc_m(2) = dc_m(2)/sigmac_m


    dr(0)=0.10D0
    dr(1)=0.60D0
    dr(2)=0.30D0


    !Compute a bunch of interpolating coefficients

    !for(int i=2 i<N+4 i++){
    q0 = (eps + si_0)*(eps + si_0)
    q1 = (eps + si_1)*(eps + si_1)
    q2 = (eps + si_2)*(eps + si_2)
    !Compute the weights:
    cl(0) = dl(0)/q0
    cl(1) = dl(1)/q1
    cl(2) = dl(2)/q2

! === With negative WENO weights splitting
    cc_p(0) = dc_p(0)/q0
    cc_p(1) = dc_p(1)/q1
    cc_p(2) = dc_p(2)/q2

    cc_m(0) = dc_m(0)/q0
    cc_m(1) = dc_m(1)/q1
    cc_m(2) = dc_m(2)/q2
! ===

! === Without negative WENO weights splitting

!       cc(0) = dc(0)/q0
!       cc(1) = dc(1)/q1
!       cc(2) = dc(2)/q2
! ===

    cr(0) = dr(0)/q0
    cr(1) = dr(1)/q1
    cr(2) = dr(2)/q2

    !// Normalize the weights    
    clave = cl(0) + cl(1) + cl(2)     

!  === With negative weights splitting
   ccave_p = cc_p(0) + cc_p(1) + cc_p(2)  
   ccave_m = cc_m(0) + cc_m(1) + cc_m(2)  
! ===

!  === Without negative weights splitting
!     ccave = cc(0) + cc(1) + cc(2)
! ===

    crave = cr(0) + cr(1) + cr(2)

    DO k = 0, 2


       cl(k) = cl(k)/clave;

!  === With negative weights splitting
      cc_p(k) = cc_p(k)/ccave_p
      cc_m(k) = cc_m(k)/ccave_m
! ===

!  === Without negative weights splitting
!       cc(k) = cc(k)/ccave	
! ===

       cr(k) = cr(k)/crave
    END DO

    !//Compute the three right and left boundary values

    rl(0) = -c16*h_im2 + c56*h_im1 + c13*h_i 
    rl(1) =  c13*h_im1 + c56*h_i   - c16*h_ip1
    rl(2) = c116*h_i   - c76*h_ip1 + c13*h_ip2

!  === With negative weights splitting
   rc(0) = h_i - (h_im2 - 2.0D0*h_im1 + h_i)/24.0D0 
   rc(1) = h_i - (h_im1 - 2.0D0*h_i   + h_ip1)/24.0D0
   rc(2) = h_i - (h_i   - 2.0D0*h_ip1 + h_ip2)/24.0D0
! ===

!  === Without negative weights splitting
!     rc(0) = -c124*h_im2 +  c112*h_im1 + c2324*h_i 
!      rc(1) = -c124*h_im1 + c1312*h_i   -  c124*h_ip1
!      rc(2) = c2324*h_i   +  c112*h_ip1 -  c124*h_ip2
     
! ===

    rr(0) =  c13*h_im2 - c76*h_im1 + c116*h_i 
    rr(1) = -c16*h_im1 + c56*h_i   +  c13*h_ip1
    rr(2) =  c13*h_i   + c56*h_ip1 -  c16*h_ip2

    !//Finally, compute the reconstruction

    val_l = cl(0)*rl(0) + cl(1)*rl(1) + cl(2)*rl(2)
!  === With negative weights splitting
    val_c = sigmac_p*(cc_p(0)*rc(0)+cc_p(1)*rc(1)+cc_p(2)*rc(2)) - sigmac_m*(cc_m(0)*rc(0)+cc_m(1)*rc(1)+cc_m(2)*rc(2))
! ===
!  === Without negative weights splitting
!    val_c = cc(0)*rc(0) +cc(1)*rc(1) +cc(2)*rc(2)
! ===
    val_r = cr(0)*rr(0) +cr(1)*rr(1) +cr(2)*rr(2)  

  END SUBROUTINE Built_WENO

  SUBROUTINE Built_WENO_CV(u_bar, val_l, val_c, val_r)


! This function computes the WENO reconstruction from cell-averages ubar
! in the x-direction in the three points at center - dx/2, center and center + dx/2.
! It uses the si_0,1,2 indicators computed by si() and negative weight splitting methode of Shi and Shu

    IMPLICIT NONE

    REAL(KIND = 8), DIMENSION(-2:2), INTENT( IN ) :: u_bar
    REAL(KIND = 8), INTENT( OUT ) :: val_l
    REAL(KIND = 8), INTENT( OUT ) :: val_c
    REAL(KIND = 8), INTENT( OUT ) :: val_r

    REAL(KIND = 8), PARAMETER :: theta = 3.0D0
    REAL(KIND = 8), PARAMETER :: c13 = 1./3
    REAL(KIND = 8), PARAMETER :: c56 = 5.D0/6.D0
    REAL(KIND = 8), PARAMETER :: c16 = 1.D0/6.D0
    REAL(KIND = 8), PARAMETER :: c76 = 7.D0/6.D0
    REAL(KIND = 8), PARAMETER :: c116 = 11.D0/6.D0
    !REAL(KIND = 8), PARAMETER :: c112 = 11.D0/2.D0
    !REAL(KIND = 8), PARAMETER :: c124 = 1.D0/24.D0
    REAL(KIND = 8), PARAMETER :: eps = 1.0D-6
    !REAL(KIND = 8), PARAMETER :: c1312 =13.D0/12.D0
    !REAL(KIND = 8), PARAMETER :: c2324 =23.D0/24.D0

    REAL(KIND = 8), DIMENSION(0:2)  :: dc_p(0:2)
    REAL(KIND = 8), DIMENSION(0:2)  :: dc_m
    REAL(KIND = 8), DIMENSION(0:2)  :: dr
    REAL(KIND = 8), DIMENSION(0:2)  :: dl
    REAL(KIND = 8), DIMENSION(0:2)  :: dc
    REAL(KIND = 8), DIMENSION(0:2)  :: cl
    REAL(KIND = 8), DIMENSION(0:2)  :: cc_p
    REAL(KIND = 8), DIMENSION(0:2)  :: cc_m
    REAL(KIND = 8), DIMENSION(0:2)  :: cr
    !REAL(KIND = 8), DIMENSION(0:2)  :: cc

    REAL(KIND = 8), DIMENSION(0:2)  :: rl
    REAL(KIND = 8), DIMENSION(0:2)  :: rc
    REAL(KIND = 8), DIMENSION(0:2)  :: rr

    REAL(KIND = 8) :: sigmac_p
    REAL(KIND = 8) :: sigmac_m
    REAL(KIND = 8) :: q0
    REAL(KIND = 8) :: q1
    REAL(KIND = 8) :: q2
    REAL(KIND = 8) :: clave
    REAL(KIND = 8) :: ccave_p
    REAL(KIND = 8) :: ccave_m
    REAL(KIND = 8) :: crave
    !REAL(KIND = 8) :: ccave

    REAL(KIND = 8) :: si_0
    REAL(KIND = 8) :: si_1
    REAL(KIND = 8) :: si_2
    !REAL(KIND = 8) :: drt_top


    INTEGER  :: k

    !Call eigen_matrix(L11L12, L21, L22, Tab%h, h_i, , R11, R12, R21, R22)


    CALL smooth_indicator_cv(u_bar, si_0, si_1, si_2)

    dl(0) = 0.3D0
    dl(1) = 0.6D0
    dl(2) = 0.1D0

    ! Fith order WENO weight
    ! warning !! negative weight can lead to oscilations near discontinuities
    ! the cure are 1) the negative weight splitting technique of Shu 
    !              2) ""decrease"" accuracy to use fourth order positive WENO weights

    ! First method : fifth order WENO weights
    dc(0)= -0.11250D0
    dc(1)=  1.2250D0
    dc(2)= -0.11250D0

    ! Second method : fourth order WENO weights
!    dc(0)= 0.250D0
!    dc(1)=  0.5D0
!    dc(2)= 0.250D0

    ! Third method : negative weights splitting (for central reconstruction)

    dc_p(0) = 0.50D0*(dc(0)+theta*ABS(dc(0)))
    dc_p(1) = 0.50D0*(dc(1)+theta*ABS(dc(1)))
    dc_p(2) = 0.50D0*(dc(2)+theta*ABS(dc(2)))
    dc_m(0) = dc_p(0)-dc(0)
    dc_m(1) = dc_p(1)-dc(1)
    dc_m(2) = dc_p(2)-dc(2)

    sigmac_p = dc_p(0) + dc_p(1) + dc_p(2)
    sigmac_m = dc_m(0) + dc_m(1) + dc_m(2)


    dc_p(0) = dc_p(0)/sigmac_p
    dc_p(1) = dc_p(1)/sigmac_p
    dc_p(2) = dc_p(2)/sigmac_p

    dc_m(0) = dc_m(0)/sigmac_m
    dc_m(1) = dc_m(1)/sigmac_m
    dc_m(2) = dc_m(2)/sigmac_m


    dr(0)=0.10D0
    dr(1)=0.60D0
    dr(2)=0.30D0


    !Compute a bunch of interpolating coefficients

    !for(int i=2 i<N+4 i++){
    q0 = (eps + si_0)*(eps + si_0)
    q1 = (eps + si_1)*(eps + si_1)
    q2 = (eps + si_2)*(eps + si_2)

    !Compute the weights:
    cl(0) = dl(0)/q0
    cl(1) = dl(1)/q1
    cl(2) = dl(2)/q2

! === With negative WENO weights splitting
    cc_p(0) = dc_p(0)/q0
    cc_p(1) = dc_p(1)/q1
    cc_p(2) = dc_p(2)/q2

    cc_m(0) = dc_m(0)/q0
    cc_m(1) = dc_m(1)/q1
    cc_m(2) = dc_m(2)/q2
! ===

! === Without negative WENO weights splitting

!       cc(0) = dc(0)/q0
!       cc(1) = dc(1)/q1
!       cc(2) = dc(2)/q2
! ===

    cr(0) = dr(0)/q0
    cr(1) = dr(1)/q1
    cr(2) = dr(2)/q2

    !// Normalize the weights    
    clave = cl(0) + cl(1) + cl(2)     

!  === With negative weights splitting
   ccave_p = cc_p(0) + cc_p(1) + cc_p(2)  
   ccave_m = cc_m(0) + cc_m(1) + cc_m(2)  
! ===
    crave = cr(0) + cr(1) + cr(2)

    DO k = 0, 2


       cl(k) = cl(k)/clave;

!  === With negative weights splitting
      cc_p(k) = cc_p(k)/ccave_p
      cc_m(k) = cc_m(k)/ccave_m
! ===
       cr(k) = cr(k)/crave
    END DO

    !//Compute the three right and left boundary values

    rl(0) = -c16*u_bar(-2) + c56*u_bar(-1) + c13*u_bar(0) 
    rl(1) =  c13*u_bar(-1) + c56*u_bar(0)   - c16*u_bar(1)
    rl(2) = c116*u_bar(0)   - c76*u_bar(1) + c13*u_bar(2)

!  === With negative weights splitting
   rc(0) = u_bar(0) - (u_bar(-2) - 2.0D0*u_bar(-1) + u_bar(0))/24.0D0 
   rc(1) = u_bar(0) - (u_bar(-1) - 2.0D0*u_bar(0)   + u_bar(1))/24.0D0
   rc(2) = u_bar(0) - (u_bar(0)   - 2.0D0*u_bar(1) + u_bar(2))/24.0D0
! ===
    rr(0) =  c13*u_bar(-2) - c76*u_bar(-1) + c116*u_bar(0) 
    rr(1) = -c16*u_bar(-1) + c56*u_bar(0)   +  c13*u_bar(1)
    rr(2) =  c13*u_bar(0)   + c56*u_bar(1) -  c16*u_bar(2)

    !//Finally, compute the reconstruction

    val_l = cl(0)*rl(0) + cl(1)*rl(1) + cl(2)*rl(2)

    val_c = sigmac_p*(cc_p(0)*rc(0)+cc_p(1)*rc(1)+cc_p(2)*rc(2)) - sigmac_m*(cc_m(0)*rc(0)+cc_m(1)*rc(1)+cc_m(2)*rc(2))

    val_r = cr(0)*rr(0) +cr(1)*rr(1) +cr(2)*rr(2)  

  END SUBROUTINE Built_WENO_CV

  SUBROUTINE WENO_CO5(h_im2, h_im1, h_i, h_ip1, h_ip2, val_c)

! This function computes the WENO reconstruction from cell-averages ubar
! in the x-direction in the three points at center - dx/2, center and center + dx/2.
! It uses the si_0,1,2 indicators computed by si() and negative weight splitting methode of Shi and Shu

    IMPLICIT NONE

    REAL(KIND = 8),INTENT( IN ) :: h_im2
    REAL(KIND = 8),INTENT( IN ) :: h_im1
    REAL(KIND = 8),INTENT( IN ) :: h_i
    REAL(KIND = 8),INTENT( IN ) :: h_ip1
    REAL(KIND = 8),INTENT( IN ) :: h_ip2
    !REAL(KIND = 8) :: val_l
    !REAL(KIND = 8) :: val_c
    !REAL(KIND = 8) :: val_r
    !LOGICAL, INTENT( IN ) :: l_var
    REAL(KIND = 8), PARAMETER :: theta = 3.0D0
    !REAL(KIND = 8), PARAMETER :: c13 = 1./3
    !REAL(KIND = 8), PARAMETER :: c56 = 5.D0/6.D0
    !REAL(KIND = 8), PARAMETER :: c16 = 1.D0/6.D0
    !REAL(KIND = 8), PARAMETER :: c76 = 7.D0/6.D0
    !REAL(KIND = 8), PARAMETER :: c116 = 11.D0/6.D0
    REAL(KIND = 8), PARAMETER :: eps = 1.0D-12
    !REAL(KIND = 8), PARAMETER :: c1312 =13.D0/12.D0

    REAL(KIND = 8), DIMENSION(0:2)  :: dc_p(0:2)
    REAL(KIND = 8), DIMENSION(0:2)  :: dc_m
    REAL(KIND = 8), DIMENSION(0:2)  :: dr
    REAL(KIND = 8), DIMENSION(0:2)  :: dl
    REAL(KIND = 8), DIMENSION(0:2)  :: dc
    REAL(KIND = 8), DIMENSION(0:2)  :: cl
    REAL(KIND = 8), DIMENSION(0:2)  :: cc_p
    REAL(KIND = 8), DIMENSION(0:2)  :: cc_m
    REAL(KIND = 8), DIMENSION(0:2)  :: cr

    !REAL(KIND = 8), DIMENSION(0:2)  :: rl
    REAL(KIND = 8), DIMENSION(0:2)  :: rc
    !REAL(KIND = 8), DIMENSION(0:2)  :: rr

    REAL(KIND = 8), INTENT( OUT ) :: val_C

    REAL(KIND = 8) :: sigmac_p
    REAL(KIND = 8) :: sigmac_m
    REAL(KIND = 8) :: q0
    REAL(KIND = 8) :: q1
    REAL(KIND = 8) :: q2
    REAL(KIND = 8) :: clave
    REAL(KIND = 8) :: ccave_p
    REAL(KIND = 8) :: ccave_m
    REAL(KIND = 8) :: crave

    REAL(KIND = 8) :: si_0
    REAL(KIND = 8) :: si_1
    REAL(KIND = 8) :: si_2
    !REAL(KIND = 8) :: drt_top


    INTEGER  :: k


    CALL smooth_indicator(h_im2, h_im1, h_i, h_ip1, h_ip2, si_0, si_1, si_2)

    dl(0) = 0.3D0
    dl(1) = 0.6D0
    dl(2) = 0.1D0


    dc(0)= -0.11250D0
    dc(1)=  1.2250D0
    dc(2)= -0.11250D0

    ! weight splitting
    dc_p(0) = 0.50D0*(dc(0)+theta*ABS(dc(0)))
    dc_p(1) = 0.50D0*(dc(1)+theta*ABS(dc(1)))
    dc_p(2) = 0.50D0*(dc(2)+theta*ABS(dc(2)))
    dc_m(0) = dc_p(0)-dc(0)
    dc_m(1) = dc_p(1)-dc(1)
    dc_m(2) = dc_p(2)-dc(2)

    sigmac_p = dc_p(0) + dc_p(1) + dc_p(2)
    sigmac_m = dc_m(0) + dc_m(1) + dc_m(2)


    dc_p(0) = dc_p(0)/sigmac_p
    dc_p(1) = dc_p(1)/sigmac_p
    dc_p(2) = dc_p(2)/sigmac_p

    dc_m(0) = dc_m(0)/sigmac_m
    dc_m(1) = dc_m(1)/sigmac_m
    dc_m(2) = dc_m(2)/sigmac_m


    dr(0)=0.10D0
    dr(1)=0.60D0
    dr(2)=0.30D0


    !Compute a bunch of interpolating coefficients

    !for(int i=2 i<N+4 i++){
    q0 = (eps + si_0)*(eps + si_0)
    q1 = (eps + si_1)*(eps + si_1)
    q2 = (eps + si_2)*(eps + si_2)

    !Compute the weights:
    cl(0) = dl(0)/q0
    cl(1) = dl(1)/q1
    cl(2) = dl(2)/q2

    cc_p(0) = dc_p(0)/q0
    cc_p(1) = dc_p(1)/q1
    cc_p(2) = dc_p(2)/q2

    cc_m(0) = dc_m(0)/q0
    cc_m(1) = dc_m(1)/q1
    cc_m(2) = dc_m(2)/q2
    !/*
    !  cc(0) = dc(0)/q0
    !  cc(1) = dc(1)/q1
    !  cc(2) = dc(2)/q2
    !*/
    cr(0) = dr(0)/q0
    cr(1) = dr(1)/q1
    cr(2) = dr(2)/q2

    !// Normalize the weights    
    clave = cl(0) + cl(1) + cl(2)     

    ccave_p = cc_p(0) + cc_p(1) + cc_p(2)  
    ccave_m = cc_m(0) + cc_m(1) + cc_m(2)   
    !//double ccave = cc(0) + cc(1) + cc(2)

    crave = cr(0) + cr(1) + cr(2)

    DO k = 0, 2


       cl(k) = cl(k)/clave;

       cc_p(k) = cc_p(k)/ccave_p
       cc_m(k) = cc_m(k)/ccave_m

       !//cc(k) = cc(k)/ccave	

       cr(k) = cr(k)/crave
    END DO

    !//Compute the three right and left boundary values

!     rl(0) = -c16*h_im2 + c56*h_im1 + c13*h_i 
!     rl(1) =  c13*h_im1 + c56*h_i   - c16*h_ip1
!     rl(2) = c116*h_i   - c76*h_ip1 + c13*h_ip2


    rc(0) = h_i - (h_im2 - 2.0D0*h_im1 + h_i)/24.0D0 
    rc(1) = h_i - (h_im1 - 2.0D0*h_i   + h_ip1)/24.0D0
    rc(2) = h_i - (h_i   - 2.0D0*h_ip1 + h_ip2)/24.0D0

    !rc(0) = -c124*ubar[i-2] +  c112*ubar[i-1] + c2324*ubar[i]; 
    ! rc(1) = -c124*ubar[i-1] + c1312*ubar[i]   -  c124*ubar[i+1];
    ! rc(2) = c2324*ubar[i]   +  c112*ubar[i+1] -  c124*ubar[i+2]; 


    !rr(0) =  c13*h_im2 - c76*h_im1 + c116*h_i 
    !rr(1) = -c16*h_im1 + c56*h_i   +  c13*h_ip1
    !rr(2) =  c13*h_i   + c56*h_ip1 -  c16*h_ip2

    !//Finally, compute the reconstruction

    val_c = sigmac_p*(cc_p(0)*rc(0)+cc_p(1)*rc(1)+cc_p(2)*rc(2)) - sigmac_m*(cc_m(0)*rc(0)+cc_m(1)*rc(1)+cc_m(2)*rc(2))


  END SUBROUTINE WENO_CO5

  SUBROUTINE smooth_indicator(h_im2, h_im1, h_i, h_ip1, h_ip2, si_0, si_1, si_2)

    IMPLICIT NONE

    REAL(KIND = 8), INTENT( IN ) :: h_im2
    REAL(KIND = 8), INTENT( IN ) :: h_im1
    REAL(KIND = 8), INTENT( IN ) :: h_i
    REAL(KIND = 8), INTENT( IN ) :: h_ip1
    REAL(KIND = 8), INTENT( IN ) :: h_ip2
    REAL(KIND = 8), INTENT( OUT ) :: si_0
    REAL(KIND = 8), INTENT( OUT ) :: si_1
    REAL(KIND = 8), INTENT( OUT ) :: si_2

    REAL(KIND = 8), PARAMETER :: c1312 = 13.0D0/12.0D0

    !//Compute the Smoothness Indicators


    si_0 = c1312*((h_im2 - 2.0D0*h_im1 + h_i)*(h_im2 - 2.0D0*h_im1 + h_i)) & 
         + ((h_im2 - 4.0D0*h_im1 + 3.0D0*h_i)*(h_im2 - 4.0D0*h_im1 + 3.0D0*h_i))/4.0D0

    si_1 = c1312*((h_im1 - 2.0D0*h_i + h_ip1)*(h_im1 - 2.0D0*h_i   + h_ip1)) & 
         + ((h_im1 - h_ip1)*(h_im1 - h_ip1))/4.0D0

    si_2 = c1312*((h_i   - 2.0D0*h_ip1 + h_ip2)*(h_i - 2.0D0*h_ip1 + h_ip2))  &
         + ((3.0D0*h_i - 4.0D0*h_ip1 + h_ip2)*(3.0D0*h_i - 4.0D0*h_ip1 + h_ip2))/4.0D0  

  END SUBROUTINE smooth_indicator

  SUBROUTINE smooth_indicator_cv(u_bar, si_0, si_1, si_2)

    IMPLICIT NONE

    REAL(KIND = 8), DIMENSION(-2:2), INTENT( IN ) :: u_bar
    REAL(KIND = 8), INTENT( OUT ) :: si_0
    REAL(KIND = 8), INTENT( OUT ) :: si_1
    REAL(KIND = 8), INTENT( OUT ) :: si_2

    REAL(KIND = 8), PARAMETER :: c1312 = 13.0D0/12.0D0

    !//Compute the Smoothness Indicators

    si_0 = c1312*((u_bar(-2) - 2.0D0*u_bar(-1) + u_bar(0))*(u_bar(-2) - 2.0D0*u_bar(-1) + u_bar(0))) & 
         + ((u_bar(-2) - 4.0D0*u_bar(-1) + 3.0D0*u_bar(0))*(u_bar(-2) - 4.0D0*u_bar(-1) + 3.0D0*u_bar(0)))/4.0D0

    si_1 = c1312*((u_bar(-1) - 2.0D0*u_bar(0) + u_bar(1))*(u_bar(-1) - 2.0D0*u_bar(0)   + u_bar(1))) & 
         + ((u_bar(-1) - u_bar(1))*(u_bar(-1) - u_bar(1)))/4.0D0

    si_2 = c1312*((u_bar(0)   - 2.0D0*u_bar(1) + u_bar(2))*(u_bar(0) - 2.0D0*u_bar(1) + u_bar(2)))  &
         + ((3.0D0*u_bar(0) - 4.0D0*u_bar(1) + u_bar(2))*(3.0D0*u_bar(0) - 4.0D0*u_bar(1) + u_bar(2)))/4.0D0  

  END SUBROUTINE smooth_indicator_cv

subroutine eigen_matrix(L11, L12, L21, L22, h, hu,&
     R11, R12, R21, R22)
 
real(KIND = 8), intent(out) :: L11
real(KIND = 8), intent(out) :: L12
real(KIND = 8), intent(out) :: L21
real(KIND = 8), intent(out) :: L22
real(KIND = 8), intent(out) :: R11
real(KIND = 8), intent(out) :: R12
real(KIND = 8), intent(out) :: R21
real(KIND = 8), intent(out) :: R22
real(KIND = 8), intent(in) :: h
real(KIND = 8), intent(in) :: hu
real(KIND = 8) :: u
real(KIND = 8) :: c

! left eigenvectors
  u = hu/h
  c = sqrt(g_r__8*h)

 L11 = (u + c)/(2.0D0*c)
 L12 = -1.0D0/(2.0D0*c)
 L21 = (-u + c)/(2.0D0*c)  
 L22 = 1.0D0/(2.0D0*c)

 ! right eigenvectors  
 R11 = 1.0D0
 R21 = u - c
 R12 = 1.0D0
 R22 = u + c

end subroutine eigen_matrix

END MODULE modWENO
