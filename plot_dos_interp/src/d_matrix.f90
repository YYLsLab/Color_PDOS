subroutine d_matrix (dy1, dy2, dy3,s,nsym,at,bg)
   !---------------------------------------------------------------
   !
   !USE symm_base, ONLY:  nsym, sr
   implicit none
   integer,parameter::max_num_sym=48
   real*8  :: dy1 (3, 3, max_num_sym), dy2 (5, 5, max_num_sym), dy3 (7, 7, max_num_sym)
   integer :: s(3,3,max_num_sym)
   real*8  :: sa(3,3)
   real*8  :: sb(3,3)
   real*8  :: sr(3,3,max_num_sym)
   integer :: nsym
   real*8  :: at(3,3)
   real*8  :: bg(3,3)
   !
   integer, parameter :: maxl = 3, maxm = 2*maxl+1, &
      maxlm = (maxl+1)*(maxl+1)
   ! maxl = max value of l allowed
   ! maxm = number of m components for l=maxl
   ! maxlm= number of l,m spherical harmonics for l <= maxl
   integer :: m, n, isym
   real*8 :: ylm(maxm, maxlm),  yl1 (3, 3), yl2(5, 5), yl3(7,7), &
      yl1_inv (3, 3), yl2_inv(5, 5),  yl3_inv(7, 7), ylms(maxm, maxlm), &
      rl(3,maxm), rrl (maxm), srl(3,maxm), delta(7,7), capel
   real*8, parameter :: eps = 1.0d-9
   !
   !  randomly distributed points on a sphere
   !
   DO isym = 1,nsym
      sa (:,:) = dble ( s(:,:,isym) )
      sb = matmul ( bg, sa )
      sr (:,:, isym) = matmul ( at, transpose (sb) )
   ENDDO

   do m = 1, maxm
      rl (1, m) = randl () - 0.5d0
      rl (2, m) = randl () - 0.5d0
      rl (3, m) = randl () - 0.5d0
      rrl (m) = rl (1,m)**2 + rl (2,m)**2 + rl (3,m)**2
   enddo
   call ylmrl ( maxlm, 2*maxl+1, rl, rrl, ylm )
   !
   !  invert Yl for each block of definite l (note the transpose operation)
   !
   !  l = 1 block
   !
   do m = 1, 3
      do n = 1, 3
         yl1 (m, n) = ylm (n, 1+m)
      end do
   end do
   call invmatl (3, yl1, yl1_inv, capel)
   !
   !  l = 2 block
   !
   do m = 1, 5
      do n = 1, 5
         yl2 (m, n) = ylm (n, 4+m)
      end do
   end do
   call invmatl (5, yl2, yl2_inv, capel)
   !
   !  l = 3 block
   !
   do m = 1, 7
      do n = 1, 7
         yl3 (m, n) = ylm (n, 9+m)
      end do
   end do
   call invmatl (7, yl3, yl3_inv, capel)
   !
   ! now for each symmetry operation of the point-group ...
   !
   do isym = 1, nsym
      !
      ! srl(:,m) = rotated rl(:,m) vectors
      !
      srl = matmul (sr(:,:,isym), rl)
      !
      call ylmrl ( maxlm, maxm, srl, rrl, ylms )
      !
      !  find  D_S = Yl_S * Yl_inv (again, beware the transpose)
      !
      !  l = 1
      !
      do m = 1, 3
         do n = 1, 3
            yl1 (m, n) = ylms (n, 1+m)
         end do
      end do
      dy1 (:, :, isym) = matmul (yl1(:,:), yl1_inv(:,:))
      !
      !  l = 2 block
      !
      do m = 1, 5
         do n = 1, 5
            yl2 (m, n) = ylms (n, 4+m)
         end do
      end do
      dy2 (:, :, isym) = matmul (yl2(:,:), yl2_inv(:,:))
      !
      !  l = 3 block
      !
      do m = 1, 7
         do n = 1, 7
            yl3 (m, n) = ylms (n, 9+m)
         end do
      end do
      dy3 (:, :, isym) = matmul (yl3(:,:), yl3_inv(:,:))
      !
   enddo

   return
contains

   function randl ()
      ! x=randl() : generate uniform real(DP) numbers x in [0,1]
      !
      REAL*8 :: randl
      !
      INTEGER , PARAMETER  :: m    = 714025, &
         ia   = 1366, &
         ic   = 150889, &
         ntab = 97
      REAL*8, PARAMETER  :: rm = 1.d0 / m
      INTEGER              :: j
      INTEGER, SAVE        :: ir(ntab), iy, idum=0
      LOGICAL, SAVE        :: first=.true.
      !

      IF ( first ) THEN
         !
         first = .false.
         idum = MOD( ic - idum, m )
         !
         DO j=1,ntab
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
         END DO
         idum=mod(ia*idum+ic,m)
         iy=idum
      END IF
      j=1+(ntab*iy)/m
      iy=ir(j)
      randl=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      !
      RETURN
      !
   end function randl

   subroutine ylmrl (lmax2, ng, g, gg, ylm)
      !-----------------------------------------------------------------------
      !
      !     Real spherical harmonics ylm(G) up to l=lmax
      !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
      !     Numerical recursive algorithm based on the one given in Numerical
      !     Recipes but avoiding the calculation of factorials that generate
      !     overflow for lmax > 11
      !
      implicit none
      !
      integer, intent(in) :: lmax2, ng
      real*8, intent(in) :: g (3, ng), gg (ng)
      !
      ! BEWARE: gg = g(1)^2 + g(2)^2 +g(3)^2  is not checked on input
      !         incorrect results will ensue if the above does not hold
      !
      real*8, intent(out) :: ylm (ng,lmax2)
      !
      ! local variables
      !
      real*8, parameter :: pi=3.1415926535898
      real*8, parameter :: fpi=12.566370614359
      real*8, parameter :: eps = 1.0d-9
      real*8, allocatable :: cost (:), sent(:), phi (:), Q(:,:,:)
      real*8 :: c, gmod
      integer :: lmax, ig, l, m, lm
      !
      if (ng < 1 .or. lmax2 < 1) return
      do lmax = 0, 25
         if ((lmax+1)**2 == lmax2) go to 10
      end do
10    continue
      !
      if (lmax == 0) then
         ylm(:,1) =  sqrt (1.d0 / fpi)
         return
      end if
      !
      !  theta and phi are polar angles, cost = cos(theta)
      !
      allocate(cost(ng))
      allocate(sent(ng))
      allocate(phi(ng))
      allocate(Q(ng,0:lmax,0:lmax) )
      !
      !$omp parallel default(shared), private(ig,gmod,lm,l,c,m)

      !$omp do
      do ig = 1, ng
         gmod = sqrt (gg (ig) )
         if (gmod < eps) then
            cost(ig) = 0.d0
         else
            cost(ig) = g(3,ig)/gmod
         endif
         !
         !  beware the arc tan, it is defined modulo pi
         !
         if (g(1,ig) > eps) then
            phi (ig) = atan( g(2,ig)/g(1,ig) )
         else if (g(1,ig) < -eps) then
            phi (ig) = atan( g(2,ig)/g(1,ig) ) + pi
         else
            phi (ig) = sign( pi/2.d0,g(2,ig) )
         end if
         sent(ig) = sqrt(max(0d0,1.d0-cost(ig)**2))
         ! write(*,*) "phi=",phi(1),"cos=",cost(1),"sin=",sent(1)
      enddo
      !
      !  Q(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m) where
      !  P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
      !
      lm = 0
      do l = 0, lmax
         c = sqrt (DBLE(2*l+1) / fpi)
         if ( l == 0 ) then
            !$omp do
            do ig = 1, ng
               Q (ig,0,0) = 1.d0
            end do
         else if ( l == 1 ) then
            !$omp do
            do ig = 1, ng
               Q (ig,1,0) = cost(ig)
               Q (ig,1,1) =-sent(ig)/sqrt(2.d0)
            end do
         else
            !
            !  recursion on l for Q(:,l,m)
            !
            do m = 0, l - 2
               !$omp do
               do ig = 1, ng
                  Q(ig,l,m) = cost(ig)*(2*l-1)/sqrt(DBLE(l*l-m*m))*Q(ig,l-1,m) &
                     - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*Q(ig,l-2,m)
               end do
            end do
            !$omp do
            do ig = 1, ng
               Q(ig,l,l-1) = cost(ig) * sqrt(DBLE(2*l-1)) * Q(ig,l-1,l-1)
            end do
            !$omp do
            do ig = 1, ng
               Q(ig,l,l)   = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent(ig)*Q(ig,l-1,l-1)
            end do
         end if
         !
         ! Y_lm, m = 0
         !
         lm = lm + 1
         !$omp do
         do ig = 1, ng
            ylm(ig, lm) = c * Q(ig,l,0)
         end do
         !
         do m = 1, l
            !
            ! Y_lm, m > 0
            !
            lm = lm + 1
            !$omp do
            do ig = 1, ng
               ylm(ig, lm) = c * sqrt(2.d0) * Q(ig,l,m) * cos (m*phi(ig))
            end do
            !
            ! Y_lm, m < 0
            !
            lm = lm + 1
            !$omp do
            do ig = 1, ng
               ylm(ig, lm) = c * sqrt(2.d0) * Q(ig,l,m) * sin (m*phi(ig))
            end do
         end do
      end do
      !
      !$omp end parallel
      !
      deallocate(cost, sent, phi, Q)
      !
      return
   end subroutine ylmrl

   subroutine invmatl (n, a, a_inv, da)
      !-----------------------------------------------------------------------
      ! computes the inverse "a_inv" of matrix "a", both dimensioned (n,n)
      ! if the matrix is dimensioned 3x3, it also computes determinant "da"
      ! matrix "a" is unchanged on output - LAPACK
      !
      implicit none
      integer :: n
      real*8, DIMENSION (n,n) :: a, a_inv
      real*8 :: da
      !
      integer :: info, lda, lwork, ipiv (n)
      ! info=0: inversion was successful
      ! lda   : leading dimension (the same as n)
      ! ipiv  : work space for pivoting (assumed of length lwork=n)
      real*8 :: work (n)
      ! more work space
      !
      lda = n
      lwork=n
      !
      a_inv(:,:) = a(:,:)
      !
      call dgetrf (n, n, a_inv, lda, ipiv, info)
      call dgetri (n, a_inv, lda, ipiv, work, lwork, info)
      !
      if (n == 3) then
         da = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
            a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
            a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
      else
         da = 0.d0
      end if

      return
   end subroutine invmatl

end subroutine d_matrix
