*                              D   I   S   C   L   A   I   M   E   R
*
*       You are using a BETA version of the program fourier_lib_module.f, which is currently
*       under development by D. A. Spong of the Fusion Energy Division,
*       Oak Ridge National Laboratory.  Please report any problems or comments
*       to him.  As a BETA version, this program is subject to periodic change
*       and improvement without notice.
*
      module kind_spec
c
c     Double precision
c
      integer, parameter :: rprec = selected_real_kind(12,100)
      integer, parameter :: iprec = selected_int_kind(8)
      end module kind_spec
c
c
c Variables that need to be set for each time the configuration is changed:
c   ith, izt, nfp, mpol, ntor, mp_col, nt_col, mode_family, tokamak
c
      module fourier_lib
      use kind_spec
      implicit none
      integer :: ith, izt, nfp, mode_family
      integer :: mpol, ntor
      integer :: nznt, mnmx
      integer :: i, j, lg, nu, nl, m, n, mn, istat
c   Equilibrium coefficient arrays
      real(kind=rprec), allocatable :: thtgrd(:),ztgrd(:),rn(:),rm(:)
      real(kind=rprec), allocatable :: fnm(:), f(:), anm(:)
      real(kind=rprec), allocatable :: cos_ar(:,:),sin_ar(:,:)
      real(kind=rprec), allocatable :: cos_toF(:,:),sin_toF(:,:)
      integer, allocatable :: in(:), im(:)
      real(kind=rprec) twopi, arg
c   Eigenfunction arrays,variables
      integer :: count, mn_col, ith_col, izt_col, ntors
      real(kind=rprec), allocatable :: rn_col(:),rm_col(:),
     >     rho(:), y_phi(:)
      real(kind=rprec), allocatable :: rn2_col(:),rm2_col(:),rnm_col(:)
      integer, allocatable :: in_col(:), im_col(:)
      integer, allocatable :: ntora(:), mpol_lwr(:), mpol_upr(:)
c   Fourier mode block sub-matrices
      real(kind=rprec), allocatable :: a(:,:,:), b(:,:,:), d(:,:,:),
     1  aa(:,:,:), bb(:,:,:), dd(:,:,:), dd1(:,:,:), dd2(:,:,:),
     2  cc(:,:,:)

      contains

      subroutine readin
       open(unit=20,file="fourier.dat",status="old")
       read(20,*) nfp, ith, izt, mode_family
c       mpol = ith/5; ntor = izt/5
       mpol = ith*2/5; ntor = izt*2/5
       nznt=ith*izt; mnmx=(2*ntor+1)*mpol-ntor
       read(20,*) ntors
       allocate(ntora(ntors), stat=istat)
       allocate(mpol_lwr(ntors), stat=istat)
       allocate(mpol_upr(ntors), stat=istat)
       do i=1,ntors
        read(20,*) ntora(i), mpol_lwr(i), mpol_upr(i)
       end do       
       close(unit=20)
      end subroutine readin
c            
      subroutine trig_array
      use kind_spec
c      use fourier_internal_vars
      implicit none
      real(kind=rprec) :: dum,dnorm
      allocate(thtgrd(nznt), stat=istat)
      allocate(ztgrd(nznt), stat=istat)
      allocate(im(mnmx), stat=istat)
      allocate(in(mnmx), stat=istat)
      allocate(rm(mnmx), stat=istat)
      allocate(rn(mnmx), stat=istat)
      allocate(fnm(mnmx), stat=istat)
      allocate(f(nznt), stat=istat)
      allocate(anm(mnmx), stat=istat)
      allocate(cos_ar(nznt,mnmx), stat=istat)
      allocate(sin_ar(nznt,mnmx), stat=istat)
      allocate(cos_toF(nznt,mnmx), stat=istat)
      allocate(sin_toF(nznt,mnmx), stat=istat)
      twopi = 8.d0*atan(1.d0)
c    Generate theta, zeta grid
      lg = 0
      do i=1,izt
      do j=1,ith
       lg = lg + 1
       ztgrd(lg) = twopi*dble(i-1)/(dble(nfp*izt))
       thtgrd(lg) = twopi*dble(j-1)/dble(ith)
      end do
      end do
c    Generate Fourier mode distribution
      mn = 0
      nu = ntor
      do m=0,mpol-1
       nl = -ntor
       if(m .eq. 0) nl = 0
        do n = nl,nu
         mn = mn + 1
         rm(mn) = dble(m)
         rn(mn) = dble(n*nfp)
         im(mn) = m
         in(mn) = n*nfp
        enddo
      enddo
      do i=1,nznt
       do mn=1,mnmx
        arg = -rn(mn)*ztgrd(i) + rm(mn)*thtgrd(i)
        cos_ar(i,mn) = cos(arg)
        sin_ar(i,mn) = sin(arg)
        dnorm = 2.d0/dble(nznt)
        dum = abs(rn(mn)) + abs(rm(mn))
        if(nint(dum) .eq. 0) dnorm = .5d0*dnorm
        cos_toF(i,mn) = cos_ar(i,mn)*dnorm
        sin_toF(i,mn) = sin_ar(i,mn)*dnorm
       end do
      end do
      end subroutine trig_array
      
      subroutine convolution_array
      use kind_spec
      implicit none
c     First, count the number of modes to be used
      count = 0
      do i=1,ntors
       count = count + mpol_upr(i) - mpol_lwr(i) + 1
      end do
      mn_col = count
      write(*,'("mn_col = ",i5)') mn_col
c
c     Allocate memory for arrays needed in eigenfunction calculations
      allocate(rm_col(mn_col), stat=istat)
      allocate(rn_col(mn_col), stat=istat)
      allocate(im_col(mn_col), stat=istat)
      allocate(in_col(mn_col), stat=istat)
      allocate(rm2_col(mn_col), stat=istat)
      allocate(rn2_col(mn_col), stat=istat)
      allocate(rnm_col(mn_col), stat=istat)
c     Create mode number arrays needed for eigenfunction
       count = 0
       do n = 1,ntors
        do m = mpol_lwr(n), mpol_upr(n)
         count = count + 1
         rm_col(count) = real(m)
         rn_col(count) = real(ntora(n))
         im_col(count) = m
         in_col(count) = ntora(n)
	 write(*,*) count, im_col(count),in_col(count)
         rm2_col(count) = rm_col(count)*rm_col(count)
         rn2_col(count) = rn_col(count)*rn_col(count)
         rnm_col(count) = rm_col(count)*rn_col(count)
        end do
      end do
      mn_col = count   
c      write(*,*) mn_col, count
      end subroutine convolution_array
c
c
      subroutine scs_convolve(ans,m1,n1,m2,n2,meq,neq)
      use kind_spec
      implicit none
      real(kind=rprec) :: tht_int1, tht_int2, tht_int3, tht_int4,
     1  zeta_int1, zeta_int2, zeta_int3, zeta_int4, ans
      integer :: m1, m2, n1, n2, meq, neq
      integer :: sm1, sm2, sn1, sn2, smeq, sneq
c
c     This subroutine calculates the 2D integral for tht and zeta running
c     from 0 to 2*PI of:
c     sin(m1*tht - n1*zeta)*cos(meq*tht - neq*zeta)*sin(m2*tht - n2*zeta)
c     The value returned is actually (2/PI) times this integral.
c
      sm1=1;sm2=1;sn1=1;sn2=1;smeq=1;sneq=1
      if(m1 .ne. 0) sm1 = m1/abs(m1)
      if(m2 .ne. 0) sm2 = m2/abs(m2)
      if(n1 .ne. 0) sn1 = n1/abs(n1)
      if(n2 .ne. 0) sn2 = n2/abs(n2)
      if(meq .ne. 0) smeq = meq/abs(meq)
      if(neq .ne. 0) sneq = neq/abs(neq)
      m1 = sm1*m1; n1 = sn1*n1
      m2 = sm2*m2; n2 = sn2*n2
      meq = smeq*meq; neq = sneq*neq
      call css(tht_int1, meq, m1, m2)
      call ccc(zeta_int1, n1, n2, neq)
      call ccc(tht_int2, m1, m2, meq)
      call css(zeta_int2, neq, n1, n2)
      call css(tht_int3, m2, m1, meq)
      call css(zeta_int3, n1, n2, neq)
      call css(tht_int4, m1, m2, meq)
      call css(zeta_int4, n2, n1, neq)
      ans = tht_int1*zeta_int1*sm1*sm2
     1    + tht_int2*zeta_int2*sn1*sn2
     2    - tht_int3*zeta_int3*sm1*smeq*sn2*sneq
     3    - tht_int4*zeta_int4*sm2*smeq*sn1*sneq
      m1 = sm1*m1; n1 = sn1*n1
      m2 = sm2*m2; n2 = sn2*n2
      meq = smeq*meq; neq = sneq*neq
      end subroutine scs_convolve
c
c
c
      subroutine ccc_convolve(ans,m1,n1,m2,n2,meq,neq)
      use kind_spec
      implicit none
      real(kind=rprec) :: tht_int1, tht_int2, tht_int3, tht_int4,
     1  zeta_int1, zeta_int2, zeta_int3, zeta_int4, ans
      integer :: m1, m2, n1, n2, meq, neq
      integer :: sm1, sm2, sn1, sn2, smeq, sneq
c
c     This subroutine calculates the 2D integral for tht and zeta running
c     from 0 to 2*PI of:
c     cos(m1*tht - n1*zeta)*cos(meq*tht - neq*zeta)*cos(m2*tht - n2*zeta)
c     The value returned is actually (2/PI) times this integral.
c
      sm1=1;sm2=1;sn1=1;sn2=1;smeq=1;sneq=1
      if(m1 .ne. 0) sm1 = m1/abs(m1)
      if(m2 .ne. 0) sm2 = m2/abs(m2)
      if(n1 .ne. 0) sn1 = n1/abs(n1)
      if(n2 .ne. 0) sn2 = n2/abs(n2)
      if(meq .ne. 0) smeq = meq/abs(meq)
      if(neq .ne. 0) sneq = neq/abs(neq)
      m1 = sm1*m1; n1 = sn1*n1
      m2 = sm2*m2; n2 = sn2*n2
      meq = smeq*meq; neq = sneq*neq
      call ccc(tht_int1, m1, m2, meq)
      call ccc(zeta_int1, n1, n2, neq)
      call css(tht_int2, meq, m1, m2)
      call css(zeta_int2, neq, n1, n2)
      call css(tht_int3, m1, m2, meq)
      call css(zeta_int3, n1, n2, neq)
      call css(tht_int4, m2, m1, meq)
      call css(zeta_int4, n2, n1, neq)
      ans = tht_int1*zeta_int1
     1    + tht_int2*zeta_int2*sm1*sm2*sn1*sn2
     2    + tht_int3*zeta_int3*sm2*smeq*sn2*sneq
     3    + tht_int4*zeta_int4*sm1*smeq*sn1*sneq
      m1 = sm1*m1; n1 = sn1*n1
      m2 = sm2*m2; n2 = sn2*n2
      meq = smeq*meq; neq = sneq*neq
      end subroutine ccc_convolve
c
c
c
      subroutine ccc(result, i, j, k)
      use kind_spec
      real(kind=rprec), parameter :: zero = 0, one = 1,
     1    two = 2, four = 4
      real(kind=rprec) :: result
      integer :: i, j, k, izeros
c
c     This subroutine calculates the 1D integral of
c     cos(i*x)*cos(j*x)*cos(k*x)
c     for x running from 0 to 2*PI.
c     The value returned is actually (2/PI) times this integral.
c
      result = zero
      izeros = 0
      if(i .eq. 0) izeros = izeros + 1
      if(j .eq. 0) izeros = izeros + 1
      if(k .eq. 0) izeros = izeros + 1
      if(izeros .ne. 0) then
       if(izeros .eq. 3) then
        result = four
        return
       else if(izeros .eq. 2) then
        result = zero
        return
       else if(izeros .eq. 1) then
        if(i .eq. 0) then
         result = zero
         if(j*j .eq. k*k) then
          result = two
         endif
         return
        else if(j .eq. 0) then
          result = zero
          if(i*i .eq. k*k) then
          result = two
         endif
         return
        else if(k .eq. 0) then
          result = zero
          if(i*i .eq. j*j) then
          result = two
         endif
         return
        endif
       endif
      else if (izeros .eq. 0) then
       if(k .eq. (i-j)) result = one
       if(k .eq. (i+j)) result = one
       if(k .eq. (j-i)) result = one
       if(k .eq. -(i+j)) result = one
      endif
      end subroutine ccc
c
c
c
      subroutine css(result, k, i, j)
      use kind_spec
      real(kind=rprec), parameter :: one = 1.d0, neg_one = -1.d0,
     1   zero = 0.d0, two = 2.d0, neg_two = -2.d0
      integer :: i, j, k
      real(kind=rprec) :: result
c
c     This subroutine calculates the 1D integral of
c     cos(k*x)*sin(i*x)*sin(j*x)
c     for x running from 0 to 2*PI.
c     The value returned is actually (2/PI) times this integral.
c
      result = zero
      if(i .eq. 0 .or. j .eq. 0) return
      if(k .eq. 0) then
        if(i*i .ne. j*j) then
         return
        else if(i*i .eq. j*j) then
         result = two
         return
        endif
      endif
      if(k .eq. (i-j)) result = one
      if(k .eq. (i+j)) result = neg_one
      if(k .eq. (j-i)) result = one
      if(k .eq. -(i+j)) result = neg_one
      end subroutine css
c
c
      subroutine trig_deallocate
      integer :: istat
c      deallocate(thtgrd,ztgrd,rm,rn,cos_ar,sin_ar,
c     1   rm_col,rn_col,im_col,in_col,
c     2   rm2_col,rn2_col,rnm_col,stat=istat)
      deallocate(thtgrd,ztgrd,cos_ar,sin_ar,stat=istat)
      end subroutine trig_deallocate
      
      subroutine toFourier(ftype)
      use kind_spec
      implicit none
      character*1 ftype
c
c   Do Fourier transform integrations needed to convert data on a
c   theta, zeta grid [stored in array f(i=1,nznt)] to a set
c   of Fourier amplitudes [stored in array fnm(mn=1,mnmx)].
c   Typically, the number of grid points in each direction needs
c   to be > 3*number of modes used in each direction to avoid
c   aliasing errors(implies nznt > 9*mnmx).
c   
c       if(ftype .eq. 's') then
c         fnm = matmul(f,sin_toF)
c       else if(ftype .eq. 'c') then
c         fnm = matmul(f,cos_toF)
c       endif
      do mn=1,mnmx     !loop over Fourier modes
       fnm(mn) = 0.d0
      do i=1,nznt      !loop over theta,zeta grid
       if(ftype .eq. 's') then
         fnm(mn) = fnm(mn) + f(i)*sin_toF(i,mn)
       else if(ftype .eq. 'c') then
         fnm(mn) = fnm(mn) + f(i)*cos_toF(i,mn)
       endif
      end do
      end do
c
      return
      end subroutine toFourier
c
c
      subroutine toReal(ftype)
c
c    Convert Fourier mode representation [stored in array anm(mn=1,mnmx)]
c    to values of function on a regularly spaced 2D grid
c    [stored in array f(i=1,nznt)].
c
      use kind_spec
      implicit none
      character*1 ftype
      do i=1,nznt
      f(i) = 0.d0
        do mn=1,mnmx
         if(ftype .eq. 's') then
          f(i) = f(i) + anm(mn)*sin_ar(i,mn)
         else if(ftype .eq. 'c') then
          f(i) = f(i) + anm(mn)*cos_ar(i,mn)
         endif
        end do
       end do
      return
      end subroutine toReal
c
c
c
      subroutine dbydth(ftype)
c
c    Take the theta derivative of the input Fourier amplitude array, fnm
c    and place the result in the output Fourier amplitude array, anm.
c    Changes to the sin/cos parity are reflected through the sin_type and
c    cos_type variables.
c
      use kind_spec
      implicit none
      character*1 ftype
      do i=1,mnmx
         if(ftype .eq. 's') then
          anm(i) = rm(i)*fnm(i)
         else if(ftype .eq. 'c') then
          anm(i) = -rm(i)*fnm(i)
         endif
        end do
c
c   following lines for resetting ftype cannot be used if
c    dbydth is called with constant arg (i.e., dbydth('c')
c    as it is in ae_solve.f
c
c        if(ftype .eq. 's') then
c          ftype = 'c'
c         else if(ftype .eq. 'c') then
c          ftype = 's'
c        endif
       return
       end subroutine dbydth
c
c
c
      subroutine dbydzt(ftype)
c
c    Take the zeta derivative of the input Fourier amplitude array, fnm
c    and place the result in the output Fourier amplitude array, anm.
c    Changes to the sin/cos parity are reflected through the sin_type and
c    cos_type variables.
c
      use kind_spec
      implicit none
      character*1 ftype
      do i=1,mnmx
         if(ftype .eq. 's') then
          anm(i) = -rn(i)*fnm(i)
          else if(ftype .eq. 'c') then
          anm(i) = rn(i)*fnm(i)
         endif
        end do
c
c   following lines for resetting ftype cannot be used if
c    dbydzt is called with constant arg (i.e., dbydzt('c')
c    as it is in ae_solve.f
c
c        if(ftype .eq. 's') then
c          ftype = 'c'
c         else if(ftype .eq. 'c') then
c          ftype = 's'
c        endif
       return
       end subroutine dbydzt
c
c
      end module fourier_lib
