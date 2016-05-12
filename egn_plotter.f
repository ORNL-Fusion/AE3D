*                              D   I   S   C   L   A   I   M   E   R
*
*       You are using a BETA version of the program egn_plotter.f, which is currently
*       under development by D. A. Spong of the Fusion Energy Division,
*       Oak Ridge National Laboratory.  Please report any problems or comments
*       to him.  As a BETA version, this program is subject to periodic change
*       and improvement without notice.
*
*      To compile (Intel compiler):
*
*      ifort -c -O egn_plotter.f
*      ifort -O -o xplot_egn_nw egn_plotter.o `PKG_CONFIG_PATH=/usr/local/plplot5.9.7/lib/pkgconfig pkg-config --cflags --libs plplotd-f77`
*
*   3/24/2010 Set up to plot vs. sqrt(toroidal flux) as independent variable
*
*       To change, see lines 274, 367
*   7/22/2009 changes: Modified to plot eigenmodes with negative values of
*      omega**2 (if they are present). To keep these separate, plots for modes with
*      omega**2 .ge. 0. are labeled with omega while those with omega**2 .lt. 0.
*      are labeled with the omega**2 value.
*
*    Command line arguments (these are optional, if no arguments are given,
*     then the following default values are used):
*    1 - inum ( = 0 by default) - format i10
*        if inum = 0, then all eigenvectors are plotted, provided
*        their associated frequency (eigenvalue) < freq_max. A single
*        quantity [phi(rho,theta=0,zeta=0)] is plotted along with the
*        frequency and dominant m,n pair
*        if inum > 0, then only the inum-th eigenvector is plotted and the
*        the individual m,n components of phi are plotted. Different colored
*        curves are used to distinguish them, but at this time, these are
*        not labeled with m,n values
*
*    2 - fscale ( = 1.0 by default) - format f10.3
*        fscale is a scaling factor for the frequencies. If the density
*        normalization and magnetic field level are set to what one wants
*        in AE3D and if the frequency in kHz is desired, then fscale = 1
*        is appropriate. If not, then other values of fscale can be used
*
*    3 - freq_max (= 500 by default) - format f10.3
*        as indicated above, when inum = 0, freq_max limits the eigenfunctions
*        displayed to only those with f < freq_max
*
      module egn_dat
       real*8, DIMENSION(:,:,:), ALLOCATABLE :: egn_vectors
       real*8, DIMENSION(:), ALLOCATABLE :: dm, rho, y_phi, dm_rd
      end module egn_dat
      
      program egn_plot
      USE PLPLOT
      use egn_dat
      implicit none
      integer :: mn_col, ns, i, j, k, istat, nmat, ns1, icol, inum,
     >  kz, kt, itheta, izeta, nfp, m, ic, nmat_rd, mpol
      real*8, DIMENSION(:,:,:), ALLOCATABLE :: egn_vec_extnd
      real*8, DIMENSION(:,:), ALLOCATABLE :: y_plot, egn_vectors0,
     >  egn_vectors_sngl, egn_vec_extnd_sngl, phi_00plane
      real*8, DIMENSION(:), ALLOCATABLE :: rr, yy, egn_max_local,
     >  ymn, ymx, ax
      real*8, DIMENSION(:), ALLOCATABLE :: egn_max,es_energy,em_energy
      integer, DIMENSION(:), ALLOCATABLE :: im_col, in_col,
     >  mmm, nnn, iacept, i_orig
      real*8 :: rho_eval, phi, ymin, ymax, arg, thetang, zetang, twopi,
     >  potential, fscale, freq_max, dm_test, dum1, egn_max_local0,
     >  cutoff, viz_flux, xcol, moment_0, moment_r, moment_r2
      logical, dimension(:), allocatable :: egn_keep
      logical :: kp, b_and_w
     
      integer :: clrs(8)
      character*5 case_num
      character*35 title
      character*15 eigenvalue
      character*4 :: mmm_ch, nnn_ch
      character*18 :: title1
      character*10 :: arg1, arg2, arg3
      character*4 :: input_form
      logical :: visit_out, displacement, visit_vs_radius
      displacement = .false.    !choice to plot displacement rather than potential
      cutoff = 1.e+5  !only put eigenmodes in curve file that are within cutoff factor of max.
      input_form = "asci"
      b_and_w = .false.
      visit_out = .true.      
      visit_vs_radius = .true.
      write(*,*)input_form
      viz_flux = 0.8   !selects surface for AVS data - sync with metric_element_create.f
      freq_max = 500.
      inum = 0
      fscale = 1.d0
      clrs(1) = 1; clrs(2) = 3; clrs(3) = 4; clrs(4) = 8
      clrs(5) = 9; clrs(6) =10; clrs(7) =13; clrs(8) =15
      call getarg(1, arg1)
      call getarg(2, arg2)
      call getarg(3, arg3)
      read(arg1,'(i10)') inum
      read(arg2,'(f10.3)') fscale
      read(arg3,'(f10.3)') freq_max
      if(fscale .eq. 0.d0) fscale = 1.d0
      if(freq_max .eq. 0.d0) freq_max = 500.
      write(*,'("Command line arguments:")')
      write(*,'("inum = ",i5,3x,"fscale = ",
     >  f10.3,3x,"freq_max = ",f10.3)')
     >    inum, fscale, freq_max
      twopi = 8.d0*atan(1.d0)
      itheta = 200; izeta = 200
      if(input_form .eq. "binr") then
       open(unit=33,file="egn_mode_bin.dat",form="unformatted",
     >       status="old")
      else if(input_form .eq. "asci") then
       open(unit=33,file="egn_mode_asci.dat",
     >       status="old")
      end if
      open(unit=14,file="potential_3d",status="unknown")
      open(unit=15,file="potential_2d",status="unknown")
      if(inum .eq. 0) open(unit=78,file="cntm_extras.dat",
     >        status="unknown")
      if(input_form .eq. "binr") then
       read(33) nmat
       read(33) mn_col
       read(33) ns
       ns1 = ns     !make radial points consistent with AE3D calculation
c       ns1 = 300   !enhanced resolution for plots
       
      else if(input_form .eq. "asci") then

       read(33,*) nmat
       read(33,*) mn_col
       read(33,*) ns
       ns1 = ns     !make radial points consistent with AE3D calculation
c       ns1 = 300   !enhanced resolution for plots

       write(*,'("# of eigenmodes available = ",i10,/,
     >  "# of m,n pairs included = ",i6,/,
     >  "# of flux surfaces = ",i4)') nmat, mn_col, ns
      endif
      open(unit=25,file="egn_values.dat",status="old")
      
       allocate(rho(ns), stat=istat)
       allocate(yy(ns1), stat=istat)
       allocate(ax(ns1), stat=istat)
       allocate(im_col(mn_col), stat=istat)
       allocate(in_col(mn_col), stat=istat)
       allocate(mmm(nmat), stat=istat)
       allocate(nnn(nmat), stat=istat)
       allocate(egn_max(nmat), stat=istat)
       allocate(dm(nmat), stat=istat)
       allocate(dm_rd(nmat), stat=istat)
       allocate(ymn(nmat), stat=istat)
       allocate(ymx(nmat), stat=istat)
       allocate(y_phi(ns), stat=istat)
       allocate(y_plot(nmat,ns1), stat=istat)
       allocate(rr(ns1), stat=istat)
       allocate(egn_vectors0(mn_col,ns1), stat=istat)
       allocate(egn_vectors_sngl(mn_col,ns1), stat=istat)
       allocate(egn_vec_extnd_sngl(mn_col,ns1), stat=istat)
       allocate(iacept(nmat), stat=istat)
       allocate(i_orig(nmat), stat=istat)
       allocate(es_energy(nmat), stat=istat)
       allocate(em_energy(nmat), stat=istat)
       allocate(egn_max_local(mn_col), stat=istat)
       allocate(egn_keep(mn_col), stat=istat)
       if(inum .ne. 0) then
        do i=1,nmat
          read(25,*) dum1, es_energy(i), em_energy(i)  !new format
c          read(25,*) dum1      !old format
        end do
       end if
       do j=1,mn_col
       	egn_keep(j) = .true.
      end do
      if(input_form .eq. "binr") then
      
       read(33) im_col
       read(33) in_col
       read(33) dm
       read(33) rho
       allocate(egn_vectors(nmat,mn_col,ns), stat=istat)
       read(33) egn_vectors
       
      else if(input_form .eq. "asci") then
      
       do i=1,mn_col
        read(33,*) im_col(i)
        read(33,*) in_col(i)
       end do
       ic = 0
       do i=1,nmat
        read(33,*) dm(i)
c	dm_test = fscale*sqrt(abs(dm(i)))
	dm_test = (fscale**2)*dm(i)
c	if(dm_test .lt. freq_max) then
	if(dm_test .gt. 0. .and. dm_test .lt. freq_max**2) then
	 ic = ic + 1
	 iacept(i) = ic
	 dm_rd(ic) = (fscale**2)*dm(i)
	 i_orig(ic) = i
	else
	 iacept(i) = 0
	endif
       end do
       nmat_rd = ic
       allocate (egn_vectors(nmat_rd,mn_col,ns), stat=istat)
c       write(*,*) dm(nmat)
       do i=1,ns
        read(33,*) rho(i)
       end do
c       write(*,*) rho(ns)
       
       do j=1,nmat
        do m = 1,ns
         do i=1,mn_col
          read(33,*,END=51) egn_vectors0(i,m)
          end do
         end do
         if (inum .eq. 0) then
         if(iacept(j) .ne. 0) then
         do m = 1,ns
          do i=1,mn_col
           egn_vectors(iacept(j),i,m) = egn_vectors0(i,m)
          end do
         end do
	 endif    !if(iacept(j) .ne. 0)
        if((j/1000)*1000 .eq. j) write(*,'(i12,3x,i12)') j, nmat
	 else if(inum .ne. 0 .and. inum .eq. j) then
	  do m = 1,ns
          do i=1,mn_col
           egn_vectors_sngl(i,m) = egn_vectors0(i,m)
          end do
         end do
	end if   !if (inum .eq. 0 or 1
       end do
   51  continue
       nmat = nmat_rd
       end if    !if(input_form .eq. "asci"
       write(*,'(i12)') nmat_rd
       
       allocate(egn_vec_extnd(nmat,mn_col,ns1), stat=istat)
       allocate(phi_00plane(nmat,ns1), stat=istat)
       
       dm(:) = (fscale**2)*dm(:)
c
c       do j=1,mn_col
c        write(*,*) im_col(j), in_col(j)
c       end do
       do i=1,nmat
	egn_max(i) = -1.e+30
        do j=1,mn_col
	  do k=1,ns
	   if(abs(egn_vectors(i,j,k)) .gt. egn_max(i)) then
            egn_max(i) = abs(egn_vectors(i,j,k))
	    mmm(i) = im_col(j); nnn(i) = in_col(j)
	   end if
	   end do
          end do
	 end do
	 
        if(inum .ne. 0) then
	 do j = 1,mn_col
	  egn_max_local(j) = -1.e+30
	  do k = 1,ns
	   if(abs(egn_vectors_sngl(j,k)) .gt. egn_max_local(j)) then
	     egn_max_local(j) = abs(egn_vectors_sngl(j,k))
	   end if
	  end do
	 end do
	 egn_max_local0 = maxval(egn_max_local)
	  do j = 1,mn_col
	   if(egn_max_local(j) .gt. egn_max_local0/cutoff) then
	    egn_keep(j) = .true.
	   else
	    egn_keep(j) = .false.
	   end if
	  end do
	endif
c
      if(inum .eq. 0) then
      call plscol0(0, 255, 255, 255)
      call plscol0(15, 0, 0, 0)
      call plssub(2,2)
c      call plscolbg(255,255,255)
      call plinit
      call plscolor(1)
      call plcol0(15)
      call plfont(2)
       do i=1,nmat    !loop over eigenfunctions
	ymax = -1.e+30; ymin = 1.e+30
        do j=1,mn_col
	 mpol = im_col(j)
	 do k=1,ns
	  y_phi(k) = egn_vectors(i,j,k)
	 end do
	 do k=1,ns1
	  rr(k) = dble(k-1)/dble(ns1-1)
c	linear_fe_eval must be called using toroidal flux
	  rho_eval = rr(k)**2         !plot vs. <r>/<a>
c	  rho_eval = rr(k)            !plot vs. psi/psi_edge
	  call linear_fe_eval(rho_eval,ns,phi,mpol)
	  egn_vec_extnd(i,j,k) = phi
	  ymin = min(ymin,phi)
	  ymax = max(ymax,phi)
	 end do
	end do
	ymx(i) = ymax
	ymn(i) = ymin
       end do

      do i=1,nmat
       call plcol0(15)
c       if(dm_rd(i) .lt. 0.) cycle
       if(input_form.eq."asci")write(case_num,'(i5)') i_orig(i)
       if(input_form .eq. "binr")write(case_num,'(i5)') i
       write(mmm_ch,'(i4)') mmm(i)
       write(nnn_ch,'(i4)') nnn(i)
       if(dm_rd(i) .ge. 0.) then
        if(input_form.eq."asci")
     >    write(eigenvalue,'(e12.5)')sqrt(dm_rd(i))
        if(input_form .eq. "binr")
     >    write(eigenvalue,'(e12.5)') sqrt(dm(i))
        title = case_num // ", omega =" // eigenvalue
       else if(dm_rd(i) .lt. 0.) then
        if(input_form.eq."asci")
     >    write(eigenvalue,'(e12.5)') dm_rd(i)
        if(input_form .eq. "binr")
     >    write(eigenvalue,'(e12.5)') dm(i)
        title = case_num // ", omgsq =" // eigenvalue
       end if
       title1 = "m =" // mmm_ch // ", n =" // nnn_ch
       ymin = ymn(i); ymax = ymx(i)
       call plenv(0.d0, 1.d0, ymin, ymax, 0, 0)
       call pllab('<r>/<a>','#gf(#gh=0,#gz=0)',title)
       call plmtex("t",-2.0d0,0.5d0,0.5d0,title1)
       do j = 1,mn_col
        do k=1,ns1
         yy(k) = egn_vec_extnd(i,j,k)
        end do
	icol = mod(j,15) + 1
	if(icol .eq. 2) icol = 12
c	if(icol .eq. 5) icol = 9
c	if(icol .eq. 6) icol = 1
	call plcol0(icol)
c	xcol = dble(clrs(icol))/15.
c	call plcol1(xcol)
	call plwidth(2.d0)
        call plline(rr,yy)
       end do
       phi_00plane(i,:) = 0.
       do k = 1,ns1
        do j = 1,mn_col
	 phi_00plane(i,k) = phi_00plane(i,k) + egn_vec_extnd(i,j,k)
        end do
c	yy(k) = phi_00plane(i,k)/2.
       end do
c       call plpoin(ns1,rr,yy,0)
      end do
c
      do i=1,nmat
       do k=1,ns1
        yy(k) = abs(phi_00plane(i,k))
       end do
       call simpun(rr,yy,ns1,1,ax)
       moment_0 = ax(ns1)
       do k=1,ns1
        yy(k) = rr(k)*abs(phi_00plane(i,k))
       end do
       call simpun(rr,yy,ns1,1,ax)
       moment_r = ax(ns1)/moment_0       
       do k=1,ns1
        yy(k) = (abs(rr(k)-moment_r))*abs(phi_00plane(i,k))
       end do
       call simpun(rr,yy,ns1,1,ax)
       moment_r2 = 2.*ax(ns1)/moment_0
c       write(78,'(i3,2x,4(e12.4,2(3x,e12.4)))') i,sqrt(dm_rd(i)),
c     >       moment_0,moment_r, moment_r2       
       write(78,'(e12.4,2(3x,e12.4))') moment_r,sqrt(dm_rd(i)),
     >       moment_r2
      end do
c
      else if(inum .ne. 0) then
      if(visit_out) open(unit=21,file="egn_out.curve",status="unknown")
      if(.not.visit_out) open(unit=21,file="egn_out.dat",
     >   status="unknown")
      if(.not.visit_out) then
       if(dm(inum) .ge. 0.)
     >   write(21,'(i4,2x,i4,3(2x,e15.7))') ns1, mn_col, sqrt(dm(inum)),
     >             es_energy(inum), em_energy(inum)
       if(dm(inum) .lt. 0.)
     >   write(21,'(i4,2x,i4,3(2x,e15.7))') ns1, mn_col, dm(inum),
     >             es_energy(inum), em_energy(inum)
       do k=1,ns1
        write(21,*) real(k-1)/real(ns1-1)
       end do
      end if
      call plscol0(0, 255, 255, 255)
      call plscol0(15, 0, 0, 0)
      call plssub(1,1)
      call plinit
      call plcol0(15)
      call plfont(2)
      ymin = 1.d+30; ymax = -1.d+30
        do j=1,mn_col
	 kp = egn_keep(j)
c	 write(*,*) j, kp, egn_keep(j)
	 mpol = im_col(j)
	 if(visit_out .and. kp)
     >    write(21,'("#",i2,",",i3)') im_col(j),in_col(j)
	 if(.not.visit_out) write(21,*) im_col(j),in_col(j)
	 do k=1,ns
	  if(displacement) then
	    y_phi(k) = egn_vectors_sngl(j,k)*dble(mpol)
	  else
	    y_phi(k) = egn_vectors_sngl(j,k)
	  endif
	 end do
	 do k=1,ns1
	  rr(k) = dble(k-1)/dble(ns1-1)
c	linear_fe_eval must be called using toroidal flux
	  if(visit_vs_radius) rho_eval = rr(k)**2         !plot vs. <r>/<a>
	  if(.not.visit_vs_radius) rho_eval = rr(k)            !plot vs. psi/psi_edge
	  call linear_fe_eval(rho_eval,ns,phi,mpol)
	  egn_vec_extnd_sngl(j,k) = phi
c	  if(visit_out .and. kp) write(21,*) sqrt(rr(k)), phi
	  if(visit_out .and. kp) write(21,*) rr(k), phi
	  if(.not.visit_out) write(21,*) phi
	  ymax = max(ymax,phi); ymin = min(ymin,phi)
	 end do
	end do
	write(*,*) ymin, ymax
       if(ymin .lt. 0.d0) ymin = 1.1*ymin
       if(ymax .gt. 0.d0) ymax = 1.1*ymax
       call plenv(0.d0, 1.d0, ymin, ymax, 0, 0)
       if(dm(inum) .ge. 0.) write(eigenvalue,'(e12.5)') sqrt(dm(inum))
       if(dm(inum) .lt. 0.) write(eigenvalue,'(e12.5)') dm(inum)
       write(case_num,'(i5)') inum
       if(dm(inum) .ge. 0.) 
     >  title = case_num // ", omega =" // eigenvalue
       if(dm(inum) .lt. 0.) 
     >  title = case_num // ", omgsq =" // eigenvalue
c       call pllab('<r>/<a>','#gf#dmn#u(#gr)',title)
       call pllab('#gr','#gf#dmn#u(#gr)',title)
       do j = 1,mn_col
        do k=1,ns1
	 yy(k) = egn_vec_extnd_sngl(j,k)
	end do
	icol = mod(j,15) + 1
	if(icol .eq. 2) icol = 12
c	if(icol .eq. 5) icol = 9
c	if(icol .eq. 6) icol = 1
	call plcol0(icol)
	if(b_and_w) call plcol0(15)
c	xcol = dble(clrs(icol))/15.
c	call plcol1(xcol)
	call plwidth(2.d0)
        call plline(rr,yy)
       end do
       itheta = 300; izeta = 300
       k = (ns-1)*viz_flux
       toroidal: do kz = 1,izeta
       poloidal: do kt = 1,itheta
       thetang = twopi*dble(kt-1)/dble(itheta-1)
c       zetang = .9*twopi*dble(kz-1)/dble(izeta-1)
       zetang = twopi*dble(kz-1)/dble(izeta-1)      !for silo file generation
       potential = 0.d0
       do j = 1,mn_col
        arg = dble(im_col(j))*thetang - dble(in_col(j))*zetang
        potential = potential + egn_vectors_sngl(j,k)*cos(arg)
       end do
       write(14,'(e15.7)') potential
       end do poloidal
       end do toroidal
       
c       zetang = 0.9*twopi
       zetang = 0.
c       zetang = 1.5*(twopi/2.)/4.
       write(*,*) ns, itheta
       write(*,'("num. flux surfaces for 2d plot = ",i5)')
     >       int((ns-1)*viz_flux-1)
       flux: do k = 2,(ns-1)*viz_flux
c       flux: do k = 1,(ns-1)
       poloidal_1: do kt = 1,itheta
       thetang = twopi*dble(kt-1)/dble(itheta-1)       
       potential = 0.d0
       do j = 1,mn_col
        arg = dble(im_col(j))*thetang - dble(in_col(j))*zetang
	if(displacement) then
          potential = potential
     >            - real(im_col(j))*egn_vectors_sngl(j,k)*sin(arg)
        else
          potential = potential
     >            + egn_vectors_sngl(j,k)*cos(arg)
	endif
       end do
       write(15,'(e15.7)') potential
       end do poloidal_1
       end do flux
      endif
      write(*,'("j = 1 to ",i3)') ns
       call plend
       end program egn_plot
c
c
      subroutine linear_fe_eval(rho_eval,ns,phi_eval,m)
      use egn_dat
      implicit none
      real*8 :: rho_eval, phi_eval
      integer :: l, ns, m
      if(rho_eval .ge. 0.d0 .and. rho_eval .le. rho(1)) then
       if(m .ne. 0) then
        phi_eval = rho_eval*y_phi(1)/rho(1)
       else if(m .eq. 0) then
c        phi_eval = y_phi(1)
	phi_eval = (y_phi(2)*(rho(1)**2) - y_phi(1)*(rho(2)**2)
     >     + (y_phi(1) - y_phi(2))*(rho_eval**2))/
     >     ((rho(1)**2) - (rho(2)**2))
       end if
       return
      else if(rho_eval .ge. rho(ns) .and. rho_eval .le. 1.d0) then
       phi_eval = (rho_eval - 1.d0)*y_phi(ns)/(rho(ns) - 1.d0)
       return
      else
       l = 1
       do while (.not.(rho_eval .gt. rho(l) .and. rho_eval
     >           .le. rho(l+1)))
        l = l + 1
       end do
       phi_eval = y_phi(l)*(rho_eval-rho(l+1))/(rho(l)-rho(l+1))
     >     + y_phi(l+1)*(rho_eval-rho(l))/(rho(l+1)-rho(l))
      endif
      return
      end
c
c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine simpun(xx,fx,nx,i,ax)
c     ================================================================
c
c     program author      j. barish,
c     computing technology center, union carbide corp., nuclear div.,
c     oak ridge, tenn.
c
      implicit   none
c
      integer    nx
      integer    i
      integer    ix
      integer    ic
c
      real*8     xx(nx),fx(nx),ax(nx)
      real*8     d1,d2,d3,a2,a3
c
      if (i.gt.0) then
         ax(1)=0.0d+00
         do ix=2,nx,2
            d1=xx(ix)-xx(ix-1)
            ax(ix)=ax(ix-1)+d1/2.0d+00*(fx(ix)+fx(ix-1))
            if (nx.eq.ix) return
            d2=xx(ix+1)-xx(ix-1)
            d3=d2/d1
            a2=d3/6.0d+00*d2**2/(xx(ix+1)-xx(ix))
            a3=d2/2.0d+00-a2/d3
            ax(ix+1)=ax(ix-1)+(d2-a2-a3)*fx(ix-1)+a2*fx(ix)+a3*fx(ix+1)
	 enddo
      else
         ax(nx)=0.0d+00
         do ix=2,nx,2
            ic=nx+1-ix
            d1=xx(ic+1)-xx(ic)
            ax(ic)=ax(ic+1)+d1/2.0d+00*(fx(ic+1)+fx(ic))
            if (nx.eq.ix) return
            d2=xx(ic+1)-xx(ic-1)
            d3=d2/(xx(ic)-xx(ic-1))
            a2=d3/6.0d+00*d2**2/d1
            a3=d2/2.0d+00-a2/d3
            ax(ic-1)=ax(ic+1)+(d2-a2-a3)*fx(ic-1)+a2*fx(ic)+a3*fx(ic+1)
	 enddo
      endif
c
      return
      end
