*                              D   I   S   C   L   A   I   M   E   R
*
*       You are using a BETA version of the program ae_solve.f, which is currently
*       under development by D. A. Spong of the Fusion Energy Division,
*       Oak Ridge National Laboratory.  Please report any problems or comments
*       to him.  As a BETA version, this program is subject to periodic change
*       and improvement without notice.
*
      program ae_eigensolver
      use kind_spec
      use fourier_lib
      implicit none
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: iotaf,
     1   iotapf, jpolf, jpolpf, jtorf, jtorpf, phipf, phippf,
     2   work, dm, em, tau, ion_density, mu0_rho_ion, jprl_coef0,
     3   jprl_coef1, jprl_coef2, prespf,
     4   q1, q2, q3, wtht, wzet, d_wtht_ds, d_wzet_ds,
     5   wtht_half, wzet_half
c
c
c  Note: this code uses normalized toroidal flux (variables: rho, r_pt) as
c        an independent variable. If a radius-like variable is desired
c        (i.e., sqrt(tor_flux) or radius based on sqrt(area/pi)), this must be
c        constructed in subsequent codes that use the AE3D eigenfunction data
c
c  Version 31 includes options for the RFP (activated if lrfp = T)
c   (i.e., toroidal flux => poloidal flux as a radial coordinate) so that
c   both non-field-reversed and field-reversed cases can be treated without
c   singularities entering in through iota. This assumes that VMEC and xmetric
c   have both been run with the RFP option turned on.
c
c  Version 30 will be based on symmetric (self-adjoint) versions of the parallel
c   current and pressure terms. Previous non-self-adjoint versions will be removed.
c   Also, references to eigenvalue solvers other than DGEEV will be removed and
c   options for the RFP will be included (i.e., toroidal flux => poloidal flux
c   as a radial coordinate) so that both non-field-reversed and field-reversed cases
c   can be treated without singularities entering in through iota.
c   Version 30 has been derived from Versio 22 (7/26/2013)
c
c  Put in option to calculate matrix elements only and write them out to files
c   (without allocating storage for them in the code) and then stop (stop 25).
c   This is activated by setting build_matrix_only = .true. in the plasma.dat input file.
c   One can also do a full-matrix solve with output of matrix elements by
c   setting build_matrix_only = .false. and jdqz_data = .true.  1/26/2011.
c
c  Added Jprl related terms 1/5/2009
c
c  Modified 12/2008 so that ion density profile, ion mass and several other
c   parameters are read in through the plasma.dat file rather than being
c   hardwired into the code. (see explanation below for variables contained
c   in plasma.dat file).
c
c   6/6/2008 - corrected magnetic axis boundary conditions for m = 0 modes (d/dr = 0)
c   5/29/2008 - Added calculation of inertial and bending energies for each eigenmode
c   5/2008 - ion density, mass density, rad/sec -> kHz conversion factored in so that
c     eigenvalues come out in units of kHz
c   5/2008 - extraplation introduced for outer radial point. Previously this tended to
c     be erroneous for some equilibria and resulted in the inertia matrix not having
c     the expected positive-definiteness. The extrapolation ensures that the metrics,
c     profiles, etc. are continuous in going to this point.
c   This version deallocates memory as it goes along. Also, uneeded spaces following
c   the unary minus signs for rho_int_2drv terms in subroutine element were removed.
c   Finally, an option to use the more general eigensolver dggev is introduced.
c
c  This version of ae_solve.f has been modified to be compatible with stellgap.f.
c  The frequencies are output in kHz. The mode selection is done using the fourier.dat
c  file. The equilibrium data comes from the ae_metric.dat file (created by
c  running the metric_element_create.f code). Following the
c  ae3d run, a post-processing code, egn_plotter.f can be run. Egn_plotter.f can
c  either plot all of the eigenfunctions and frequencies, or plot individual
c  eigenfunctions.
c  The central ion density, ion density profile, and ion mass are set from the
c  plasma.dat input file.
c  In AE3D the number of radial (flux) mesh points is set by the original equilibrium
c  run. AE3D does not at this time allow for a higher resolution radial mesh as
c  stellgap.f does.
c
c
c  Input parameters:   (fourier.dat file).
c
c   ith, izt = number of theta and zeta grid points in the original mesh data
c     (in file tae_data_vmec) calculated by xcobra_vmec_tae.
c    This file contains data
c    for iota, phipc, |B|, gssup (the psi-psi contra-variant metric element,
c    the Jacobian, and the contravariant poloidal and toroidal components of
c    the magnetic field)
c   nfp = number of field periods (for a tokamak run, set nfp = 1)
c   mpol, ntor = number of poloidal and toroidal modes used in the Fourier
c    expansion of the above gridded quantities suplied by the tae_data_vmec
c    data file.  To avoid anti-aliasing errors, mpol < 0.5*ith and
c    ntor < 0.5*izt.
c   mp_col = number of m,n mode pairs used in the representation of the
c    eigebfunction
c   nt_col = number of toroidal modes used in the representation of the
c    eigenfunction (should be an integral multiple of nfp)
c   mode_family = toroidal mode number about which the toroidal eigenfunction
c     spectrum is built. i.e., for m = 0,
c     n = mode_family, mode_family + nfp, mode_family + 2*nfp,
c         ..., mode_family + nt_col
c      while for m .ne. 0,
c      n = -nt_col + mode_family, -nt_col + nfp + mode_family, -nt_col
c            + 2*nfp + mode_family, ..., nt_col + mode_family
c
c   Input parameters:   (plasma.dat file).
c
c   ion_to_proton_mass = m_ion/m_proton
c   ion_density_0 = n_ion(0) = ion density (in m**-3) at magnetic axis
c   ion_profile = integer parameter that determines form of n_ion(r)/n_ion(0) fit
c      for ion_profile = 0 ion_density = [iota(rho)/iota(0)]**2
c      for ion_profile = 1 ion_density = polynomial fit
c                          = nion(1) + nion(2)*rho + nion(3)*(rho**2)
c                           + nion(4)*(rho**3) + ... + nion(9)*(rho**8) + nion(10)*(rho**9)
c               (here rho = normalized toroidal flux)
c      for ion_profile = 2 ion_density = ion_density = constant = 1.
c      for ion_profile = 3 ion_density = [1 - aion*(rho**bion)]**cion
c
c   nion(10) = array of polynomial fit coefficients for ion_profile = 1
c   aion, bion, cion = parameters for ion_profile = 3 fit
c   jdqz_data = logical variable that used to request output of the matrices
c               so that the Jacobi-Davidson solver code can be used to find
c               the cluster of 40 eigenmodes centered around a pre-specified
c               frequency
c  build_matrix_only = logical variable that controls whether the code
c               stops after writing out the matrices needed for the 
c               Jacobi-Davidson solver (run separately) or goes on to
c               solve for all of the eigenvalues using the complete solver.
c               Normally, for large problems where calculating all of the
c               eigenvalues becomes slow, one will probably want to set
c               build_matrix_only = .true. in order to avoid the lengthy
c               process of calculating all eigenvalues and go directly to
c               using the clustered solver.
c
c   egnout_form = "binr" or "asci" - specifies whether the output file
c                  (unit=33, egn_mode_bin.dat or egn_mode_asc.dat)
c                 is written in binary or text format. This is the output
c                 file that is read into the post processing egn_plotter.f
c                 code. A similar parmeter should be set in the egn_plotter.f
c                 code that is consistent with the form of the ae3d output
c                 file.
c
c   To compile (on Mac using ifort compiler with ifort math libraries):
c    MKLPATH="/Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/32"
c    MKLINCLUDE="/Library/Frameworks/Intel_MKL.framework/Versions/Current/include"
c    OPT="-O -vec-report0"
c    ifort -c $OPT fourier_lib_module.f
c    ifort -c $OPT ae_solve.f
c    ifort $OPT -o xae3d fourier_lib_module.o ae_solve.o \
c      -I. -L$MKLPATH -I$MKLINCLUDE -lmkl_intel -lmkl_lapack \
c      -lmkl_core -lguide -lmkl_intel_thread
c
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: bfield,
     1  rjacob,gsssup, gttsup, gzzsup, gstsup,
     2  gszsup, gtzsup, bfields, bfieldth, bfieldze,
     3  djds, djdt, djdz, jfcn, brho, jprl_over_b_re,
     4  jprl_over_b12_re, d_brho_dth_re, d_brho_dzt_re
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE ::
     1  gsssup_ft, gttsup_ft, gzzsup_ft, gstsup_ft,
     2  gszsup_ft, gtzsup_ft, qq, brho_ft, d_brho_dth_ft,
     3  d_brho_dzt_ft, jprl_over_b_ft, jprl_over_b12_ft,
     4  bbinv2_ft, brho_binv2_ft,
     5  jprl_ovr_b_brho_ft, jprl_ovr_b_dbrho_dtht_ft,
     6  jprl_ovr_b_dbrho_dzet_ft
      REAL(kind=rprec) :: kprl_i_nw, kprl_j_nw, d_kprl_i_nw_dtht,
     1  d_kprl_i_nw_dzet
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: alf_tht,
     1  alf_zet, iota_inv_pf
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE ::
     1  pres_a_tht, pres_a_zet, pres_c_tht, pres_c_zet,
     2  pres_d_tht, pres_d_zet
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE ::
     1  pres_a_tht_ft, pres_a_zet_ft, pres_c_tht_ft, pres_c_zet_ft,
     2  pres_d_tht_ft, pres_d_zet_ft, grad_to_prp_ft,
     3  bgradb2_ft, b2gradb_ft, b3_ft
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE ::
     1  gss_ji_djds2,gss_j,gss_djds,gst_djdt,     
     2  gst_j,gst_djds,gst_ji_djds_djdt,     
     3  gtt_djdt,gtt_j,gtt_ji_djdt2,gzz_djdz,     
     4  gzz_j,gsz_djdz,gsz_ji_djds_djdz,     
     5  gsz_djds,gsz_j,gzt_j,     
     6  gzt_djdt,gzt_djdz,gzz_ji_djdz2,
     7  gzt_ji_djdz_djdt
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: amat, bmat,
     1  xx, zz
      REAL(kind=rprec), DIMENSION(:,:,:), ALLOCATABLE :: bb1, bb2, cc2
      real(kind=rprec) :: tht, zt, two_pi, ans_ccc, ans_scs,
     1  ans_ssc, ans_css, kprl_ij, kprl_i, kprl_j, rme, rne,
     2  rni, rmi, rnj, rmj, mat_inertia_ij, mat_bending_ij, sum,
     3  dds, denom, denom2, res, norm0, res_a, alpha, beta, r_pt
      REAL(kind=rprec), DIMENSION(:,:,:), ALLOCATABLE :: egn_vectors
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: mmx, nnx, ix, jx
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
      REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE ::
     1  mat_inertia, mat_bending, lambda_i, temp1, temp2, vr, vl,
     2  ctmp1, ctmp2
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: egn_test,
     1  egn_test1, alphar, alphai, betar
c   diagnostic arrays
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: a_mag, b_mag,
     1  d_mag, aa_mag, bb_mag, dd1_mag, dd2_mag
      REAL(kind=rprec), DIMENSION(:), ALLOCATABLE ::
     >   energy_inertia, energy_bending,
     >   energy_inertia_1, energy_bending_1
      REAL(kind=rprec) :: ion_to_proton_mass
      REAL(kind=rprec) :: stell_rfp_specific_term1,
     > stell_rfp_specific_term2, stell_rfp_specific_term3,
     > stell_rfp_specific_term4, stell_rfp_specific_term5,
     > stell_rfp_specific_term6, stell_rfp_specific_term7
      integer :: nm_c
      real(kind=rprec) :: nion(10)
      real :: tsec, t0, t00, ddot, stemp1, stemp2
      real(kind=rprec) :: mu0, scale_khz, one, zero, tsec0,
     1   aion, bion, cion
      real(kind=rprec) mass_ion, ion_density_0        !ion mass, central ion density, SI units
      real(kind=rprec), parameter :: mass_proton = 1.67d-27   !proton mass, SI units
      real(kind=rprec), parameter :: z0 = 0.d0
      integer :: k, ic, ig, ns, izeta, itheta, nmat, ku, kl,
     1  info, lwork, ku1, liwork, ion_profile
      integer :: mi, mj, ni, nj, meq, neq, icount,
     1  mm, nn, ii, jj, jup, ie, je, incx, incy, itype
      character*1 :: uplo
      character*1 :: tb
      character*4 :: egnout_form
      logical :: residue_test,band_matrix_test,full_matrix_test
      logical :: jdqz_data, build_matrix_only, lrfp
      namelist /plasma_input/ ion_to_proton_mass, ion_density_0,
     >    ion_profile, jdqz_data, build_matrix_only, egnout_form,
     >    nion, aion, bion, cion
      uplo = 'U'; tb = char(9)
      build_matrix_only = .false.
      lrfp = .false.
      open(unit=4,file="plasma.dat",status="old")
      read(4,plasma_input)
      close(unit=4)
      mass_ion = mass_proton*ion_to_proton_mass
c      
c      The following solver is used:
c
c         DGGEV - Generalized non-symmetric solver (Ax = lambda*Bx)
c
      residue_test = .false.
      full_matrix_test = .true.
      band_matrix_test = .false.
      t0 = secnds(0.0)
      t00 = secnds(0.0)
      two_pi = 8.d0*atan(1.d0)
      mu0 = 2.d-7*two_pi
      zero = 0.d0
c
c     Read mode selections from fourier.dat file,
c     set-up arrays for Fourier conversions
c
      call readin
      call trig_array
      call convolution_array
      tsec = secnds(t0)
      t0 = secnds(0.0)
      write(*,'("time setup(sec) = ",e12.5)') tsec
c
c     Allocate arrays and read metric elements
c     on s, theta, zeta grid from ae_metric.dat file
c
      open(unit=15,file="ae_metric.dat",
     >       status="old")
       read(15,*) ns, izeta, itheta, ig
       write(*,*) ns, izeta, itheta, ig
       write(*,*) mnmx, ith, izt, mpol, ntor
       if(ig .ne. nznt) then
        write(*,'("Discrepancy between poloidal-toroidal",/,
     1         " grid size from fourier.dat file and that",/,
     2         " from ae_metric.dat file")')
        stop 10
       end if
      open(unit=18,file="egn_values.dat",
     >       status="unknown")
      open(unit=47,file="omega2.dat",
     >       status="unknown")
      if(egnout_form .eq. "binr") then
       open(unit=33,file="egn_mode_bin.dat",form="unformatted",
     >       status="unknown")
      else if(egnout_form .eq. "asci") then
        open(unit=33,file="egn_mode_asci.dat",
     >       status="unknown")
      endif
      if(jdqz_data .or. build_matrix_only) then
       open(unit=35,file="a_matrix.dat",status="unknown")
       open(unit=36,file="b_matrix.dat",status="unknown")
       open(unit=32,file="jdqz.dat",status="unknown")
      end if
c
       allocate(iotaf(ns),iotapf(ns), jpolf(ns), jpolpf(ns),
     1  jtorf(ns), jtorpf(ns), phipf(ns), phippf(ns),
     2  ion_density(ns), mu0_rho_ion(ns), jprl_coef0(ns),
     3  jprl_coef1(ns), jprl_coef2(ns), prespf(ns),stat=istat)
       allocate(q1(ns), q2(ns), q3(ns), wtht(ns), wzet(ns),
     1  d_wtht_ds(ns), d_wzet_ds(ns),stat=istat)
       allocate(wtht_half(ns+1), wzet_half(ns+1),stat=istat)
       allocate(bfield(ig,ns), rjacob(ig,ns),
     1  gsssup(ig,ns), gttsup(ig,ns), gzzsup(ig,ns),
     2  gstsup(ig,ns), gszsup(ig,ns),
     3  gtzsup(ig,ns), bfields(ig,ns), brho(ig,ns),
     4  bfieldth(ig,ns), bfieldze(ig,ns), djds(ig,ns),
     5  djdt(ig,ns), djdz(ig,ns), jfcn(ig,ns), stat=istat)
       do k = 1,ns
        read(15,19) iotaf(k), iotapf(k), jpolf(k), jpolpf(k),
     >    jtorf(k), jtorpf(k), phipf(k), phippf(k)
       end do
       do k = 1,ns
        do i=1,nznt
      read(15,18) rjacob(i,k), bfield(i,k), gsssup(i,k),
     >               gttsup(i,k), gzzsup(i,k), gstsup(i,k),
     >               gszsup(i,k), gtzsup(i,k), bfields(i,k),
     >               bfieldth(i,k), bfieldze(i,k)
        end do
       end do
       do k = 1,ns
        do i=1,nznt
	 rjacob(i,k) = -rjacob(i,k)
        end do
       end do
       open(unit=77,file="test.dat",status="unknown")
       do k = 1,ns
        write(77,*) k,
     1   (rjacob(nznt/2,k)*(bfield(nznt/2,k)**2)),
     2   (jtorf(k) - jpolf(k)/iotaf(k))
       end do
       close(unit=77)
c
c  Write out coefficients and B_rho that will be used subsequently
c  in AE3D to form J_prl/B:
c
       do k = 1,ns
        read(15,48) jprl_coef0(k),jprl_coef1(k),
     1              jprl_coef2(k),prespf(k)
c        if((k/4)*4 .eq. k) write(*,*) k,prespf(k)
       end do
c
      do k = 1,ns
c       prespf(k) = -prespf(k)
c       prespf(k) = 0.d0
      end do
c
       do k = 1,ns
        do i=1,nznt
	 read(15,49) brho(i,k)
        end do
       end do
   48  format(e15.7,3(2x,e15.7))
   49  format(e15.7)
c
c      Do linear extrapolation to get k = ns point since the
c      accuracy of this point coming from the metric_element_create.f
c      is often not good (to check averaged blocks vs. radius, see
c      ae_diag.dat file).
c
       iotaf(ns) = 2.*iotaf(ns-1) - iotaf(ns-2)
       iotapf(ns) = 2.*iotapf(ns-1) - iotapf(ns-2)
       jpolf(ns) = 2.*jpolf(ns-1) - jpolf(ns-2)
       jpolpf(ns) = 2.*jpolpf(ns-1) - jpolpf(ns-2)
       jtorf(ns) = 2.*jtorf(ns-1) - jtorf(ns-2)
       jtorpf(ns) = 2.*jtorpf(ns-1) - jtorpf(ns-2)
       prespf(ns) = 2.*prespf(ns-1) - prespf(ns-2)
       phipf(ns) = 2.*phipf(ns-1) - phipf(ns-2)
       phippf(ns) = 2.*phippf(ns-1) - phippf(ns-2)
       jprl_coef0(ns) = 2.*jprl_coef0(ns-1) - jprl_coef0(ns-2)
       jprl_coef1(ns) = 2.*jprl_coef1(ns-1) - jprl_coef1(ns-2)
       jprl_coef2(ns) = 2.*jprl_coef2(ns-1) - jprl_coef2(ns-2)
       do i=1,nznt
        rjacob(i,ns) = 2.*rjacob(i,ns-1) - rjacob(i,ns-2)
        bfield(i,ns) = 2.*bfield(i,ns-1) - bfield(i,ns-2)
        gsssup(i,ns) = 2.*gsssup(i,ns-1) - gsssup(i,ns-2)
        gttsup(i,ns) = 2.*gttsup(i,ns-1) - gttsup(i,ns-2)
        gzzsup(i,ns) = 2.*gzzsup(i,ns-1) - gzzsup(i,ns-2)
        gstsup(i,ns) = 2.*gstsup(i,ns-1) - gstsup(i,ns-2)
        gszsup(i,ns) = 2.*gszsup(i,ns-1) - gszsup(i,ns-2)
        gtzsup(i,ns) = 2.*gtzsup(i,ns-1) - gtzsup(i,ns-2)
        bfields(i,ns) = 2.*bfields(i,ns-1) - bfields(i,ns-2)
        bfieldth(i,ns) = 2.*bfieldth(i,ns-1) - bfieldth(i,ns-2)
        bfieldze(i,ns) = 2.*bfieldze(i,ns-1) - bfieldze(i,ns-2)
        djds(i,ns) = 2.*djds(i,ns-1) - djds(i,ns-2)
        djdt(i,ns) = 2.*djdt(i,ns-1) - djdt(i,ns-2)
        djdz(i,ns) = 2.*djdz(i,ns-1) - djdz(i,ns-2)
        jfcn(i,ns) = 2.*jfcn(i,ns-1) - jfcn(i,ns-2)
        brho(i,ns) = 2.*brho(i,ns-1) - brho(i,ns-2)
       end do
       
   18  format(e15.7,10(2x,e15.7))
   19  format(e15.7,7(2x,e15.7))
c
c
       allocate(rho(ns), stat=istat)
       allocate(iota_inv_pf(ns), stat=istat)
       scale_khz = (1.e+3*two_pi)**2
       do k=1,ns
        rho(k) = dble(k)/dble(ns+1)
	r_pt = rho(k)
        iota_inv_pf(m) = iotapf(m)/(iotaf(m)**2)   !neg. sign factored in futher on
	if(ion_profile .eq. 0) then
	  ion_density(k) = (iotaf(k)/iotaf(1))**2   !profile that lines up gaps
	else if(ion_profile .eq. 1) then
	  ion_density(k) = nion(1) + r_pt*nion(2) + nion(3)*(r_pt**2)
     1     + nion(4)*(r_pt**3) + nion(5)*(r_pt**4) + nion(6)*(r_pt**5)
     2     + nion(7)*(r_pt**6) + nion(8)*(r_pt**7) + nion(9)*(r_pt**8)
     3     + nion(10)*(r_pt**9)
        else if(ion_profile .eq. 2) then
	  ion_density(k) = 1.d0
	else if(ion_profile .eq. 3) then
	  ion_density(k) = (1. - aion*(r_pt**bion))**cion
	end if
        mu0_rho_ion(k) = mass_ion*mu0*ion_density_0
     >                   *ion_density(k)*scale_khz
       end do
       
       open(unit=97,file="profiles.dat",status="unknown")
       write(97,'("rho",a1,"den_ion",a1,"iota",a1,"iotap",a1,"jpol",
     1   a1,"jpolp",a1,"jtor",a1,"jtorp",a1,"presp",a1,"phip",a1,
     2   "phipp",a1,"jprl0",a1,"jprl1",a1,"jprl2")') tb,tb,tb,tb,tb,
     3   tb,tb,tb,tb,tb,tb,tb,tb
       do k=1,ns
        write(97,'(e15.7,13(a1,e15.7))') rho(k),tb,ion_density(k),tb,
     1  iotaf(k),tb,iotapf(k),tb,jpolf(k),tb,jpolpf(k),tb,jtorf(k),tb,
     2  jtorpf(k),tb,prespf(k),tb,phipf(k),tb,phippf(k),tb,
     3  jprl_coef0(k),tb,jprl_coef1(k),tb,jprl_coef2(k)
       end do
       close(unit=97)
c
c      Transform equilibrium metric elements to
c      Fourier space and store in arrays.
c
c    (1) Metric FT arrays for inertia term:
c        (build the required inertia coefficients
c         in real space and then transform to Fourier
c         space using the toFourier subroutine) 
c
       allocate(gsssup_ft(mnmx,ns), gttsup_ft(mnmx,ns),
     1  gzzsup_ft(mnmx,ns), gstsup_ft(mnmx,ns), 
     2  gszsup_ft(mnmx,ns), gtzsup_ft(mnmx,ns), stat=istat)
       allocate(grad_to_prp_ft(mnmx,ns), bgradb2_ft(mnmx,ns),
     >  b2gradb_ft(mnmx,ns), b3_ft(mnmx,ns), stat=istat)
      do m = 1,ns
      
       f(:)=gsssup(:,m)*mu0_rho_ion(m)*rjacob(:,m)/(bfield(:,m)**2)
       call toFourier('c')
       gsssup_ft(:,m) = fnm(:)

       f(:)=gttsup(:,m)*mu0_rho_ion(m)*rjacob(:,m)/(bfield(:,m)**2)
       call toFourier('c')
       gttsup_ft(:,m) = fnm(:)

       f(:)=gzzsup(:,m)*mu0_rho_ion(m)*rjacob(:,m)/(bfield(:,m)**2)
       call toFourier('c')
       gzzsup_ft(:,m) = fnm(:)

       f(:)=gstsup(:,m)*mu0_rho_ion(m)*rjacob(:,m)/(bfield(:,m)**2)
       call toFourier('s')
       gstsup_ft(:,m) = fnm(:)

       f(:)=gszsup(:,m)*mu0_rho_ion(m)*rjacob(:,m)/(bfield(:,m)**2)
       call toFourier('s')
       gszsup_ft(:,m) = fnm(:)

       f(:)=gtzsup(:,m)*mu0_rho_ion(m)*rjacob(:,m)/(bfield(:,m)**2)
       call toFourier('c')
       gtzsup_ft(:,m) = fnm(:)
       
       if(lrfp) then
        f(:)=-mu0_rho_ion(m)/(rjacob(:,m)*bfield(:,m)**4)
       else
        f(:)=-mu0_rho_ion(m)*(phipf(m)**2)/(rjacob(:,m)*bfield(:,m)**4)
       endif
       call toFourier('c')
       grad_to_prp_ft(:,m) = fnm(:)
       
       if(lrfp) then
        f(:) = rjacob(:,m)*bfield(:,m)*((bfieldze(:,m)/iotaf(m)
     >  + bfieldth(:,m))**2)/((jtorf(m)-jpolf(m)/iotaf(m))**4)
       else
        f(:) = rjacob(:,m)*bfield(:,m)*((bfieldze(:,m)
     >  + iotaf(m)*bfieldth(:,m))**2)/((iotaf(m)*jtorf(m)-jpolf(m))**4)
       endif
       call toFourier('c')
       bgradb2_ft(:,m) = fnm(:)
       
       if(lrfp) then
        f(:) = rjacob(:,m)*(bfield(:,m)**2)*(bfieldze(:,m)/iotaf(m)
     >    + bfieldth(:,m))/((jtorf(m)-jpolf(m)/iotaf(m))**4)
       else
        f(:) = rjacob(:,m)*(bfield(:,m)**2)*(bfieldze(:,m)
     >    + iotaf(m)*bfieldth(:,m))/((iotaf(m)*jtorf(m)-jpolf(m))**4)
       endif
       call toFourier('s')
       b2gradb_ft(:,m) = fnm(:)

       if(lrfp) then
        f(:) = rjacob(:,m)*(bfield(:,m)**3)
     >       /((jtorf(m)-jpolf(m)/iotaf(m))**4)
       else
        f(:) = rjacob(:,m)*(bfield(:,m)**3)
     >       /((iotaf(m)*jtorf(m)-jpolf(m))**4)
       endif
       call toFourier('c')
       b3_ft(:,m) = fnm(:)

      end do   ! m = 1,ns   
      tsec = secnds(t0)
      t0 = secnds(0.0)
      write(*,'("time inertia FT(sec) = ",e12.5)') tsec
c
c
c    Fill-in arrays for pressure-gradient terms
c
c  New symmetric approach
      allocate(pres_a_tht_ft(mnmx,ns), pres_a_zet_ft(mnmx,ns),
     >  pres_c_tht_ft(mnmx,ns), pres_c_zet_ft(mnmx,ns),
     >  pres_d_tht_ft(mnmx,ns), pres_d_zet_ft(mnmx,ns),
     >  stat=istat)
c
      allocate(pres_a_tht(ig,ns), pres_a_zet(ig,ns),
     >  pres_c_tht(ig,ns), pres_c_zet(ig,ns),
     >  pres_d_tht(ig,ns), pres_d_zet(ig,ns),
     >  stat=istat)
      do m = 1,ns
       f(:) = -(prespf(m)/(bfield(:,m)**5))*jtorf(m)
     >        *(jtorf(m)*bfieldze(:,m) + jpolf(m)*bfieldth(:,m))
       call toFourier('s')
       pres_a_tht_ft(:,m) = fnm(:)
c       
       f(:) = (prespf(m)/(bfield(:,m)**5))*jpolf(m)
     >        *(jtorf(m)*bfieldze(:,m) + jpolf(m)*bfieldth(:,m))
       call toFourier('s')
       pres_a_zet_ft(:,m) = fnm(:)
c       
       f(:) = -(prespf(m)/(bfield(:,m)**5))*jtorf(m)
     >        *(-brho(:,m)*bfieldze(:,m) - jpolf(m)*bfields(:,m))
       call toFourier('c')
       pres_c_tht_ft(:,m) = fnm(:)
c
       f(:) = (prespf(m)/(bfield(:,m)**5))*jpolf(m)
     >        *(-brho(:,m)*bfieldze(:,m) - jpolf(m)*bfields(:,m))
       call toFourier('c')
       pres_c_zet_ft(:,m) = fnm(:)
c
       f(:) = -(prespf(m)/(bfield(:,m)**5))*jtorf(m)
     >        *(brho(:,m)*bfieldth(:,m) - jtorf(m)*bfields(:,m))
       call toFourier('c')
       pres_d_tht_ft(:,m) = fnm(:)
c
       f(:) = (prespf(m)/(bfield(:,m)**5))*jpolf(m)
     >        *(brho(:,m)*bfieldth(:,m) - jtorf(m)*bfields(:,m))
       call toFourier('c')
       pres_d_zet_ft(:,m) = fnm(:)
c
      end do   ! m = 1,ns   
c
c     (2) Metric FT arrays for bending energy terms:
c        (build the required bending energy coefficients
c         in real space and then transform to Fourier
c         space using the toFourier subroutine) 
c
c       Calculate Jacobian derivatives:
c
       do m = 1,ns
        if(lrfp) then
         denom = jtorf(m) - jpolf(m)/iotaf(m)
         denom2 = denom**2
         dds = jtorpf(m) - jpolpf(m)/iotaf(m)
     1         + jpolf(m)*iotapf(m)/(iotaf(m)**2)
	else
         denom = iotaf(m)*jtorf(m) - jpolf(m)
         denom2 = denom**2
         dds = iotaf(m)*jtorpf(m) + jtorf(m)*iotapf(m) - jpolpf(m)
	endif
        do i=1,nznt
        if(lrfp) then
         jfcn(i,m) = 1.d0/(bfield(i,m)*rjacob(i,m))
	else
         jfcn(i,m) = phipf(m)/(bfield(i,m)*rjacob(i,m))
	endif
         djds(i,m) = bfields(i,m)/denom - bfield(i,m)*dds/denom2
         djdt(i,m) = bfieldth(i,m)/denom
         djdz(i,m) = bfieldze(i,m)/denom
        end do
       end do
      allocate(gss_ji_djds2(mnmx,ns),gss_j(mnmx,ns),      
     1  gss_djds(mnmx,ns),gst_djdt(mnmx,ns),     
     2  gst_j(mnmx,ns),gst_djds(mnmx,ns),gst_ji_djds_djdt(mnmx,ns),     
     3  gtt_djdt(mnmx,ns),gtt_j(mnmx,ns),gtt_ji_djdt2(mnmx,ns),
     4  gzz_djdz(mnmx,ns), gzz_ji_djdz2(mnmx,ns),   
     5  gzz_j(mnmx,ns),gsz_djdz(mnmx,ns),gsz_ji_djds_djdz(mnmx,ns),     
     6  gsz_djds(mnmx,ns),gsz_j(mnmx,ns),gzt_j(mnmx,ns),     
     7  gzt_djdt(mnmx,ns),gzt_djdz(mnmx,ns),
     8  gzt_ji_djdz_djdt(mnmx,ns), stat=istat)
     
      do m = 1,ns
        f(:) = gsssup(:,m)*(djds(:,m)**2)*rjacob(:,m)
        call toFourier('c')
        gss_ji_djds2(:,m) = fnm(:)
    
        f(:) = gsssup(:,m)*rjacob(:,m)*(jfcn(:,m)**2)
        call toFourier('c')
        gss_j(:,m) = fnm(:)

        f(:) = gsssup(:,m)*djds(:,m)*jfcn(:,m)*rjacob(:,m)
        call toFourier('c')
        gss_djds(:,m) = fnm(:)

        f(:) = gstsup(:,m)*djdt(:,m)*jfcn(:,m)*rjacob(:,m)
        call toFourier('c')
        gst_djdt(:,m) = fnm(:)

        f(:) = gstsup(:,m)*rjacob(:,m)*(jfcn(:,m)**2)
        call toFourier('s')
        gst_j(:,m) = fnm(:)

        f(:) = gstsup(:,m)*djds(:,m)*jfcn(:,m)*rjacob(:,m)
        call toFourier('s')
        gst_djds(:,m) = fnm(:)

        f(:) = gstsup(:,m)*djdt(:,m)*djds(:,m)*rjacob(:,m)
        call toFourier('c')
        gst_ji_djds_djdt(:,m) = fnm(:)

        f(:) = gttsup(:,m)*djdt(:,m)*jfcn(:,m)*rjacob(:,m)
        call toFourier('s')
        gtt_djdt(:,m) = fnm(:)

        f(:) = gttsup(:,m)*rjacob(:,m)*(jfcn(:,m)**2)
        call toFourier('c')
        gtt_j(:,m) = fnm(:)

        f(:) = gttsup(:,m)*(djdt(:,m)**2)*rjacob(:,m)
        call toFourier('c')
        gtt_ji_djdt2(:,m) = fnm(:)

        f(:) = gzzsup(:,m)*djdz(:,m)*jfcn(:,m)*rjacob(:,m)
        call toFourier('s')
        gzz_djdz(:,m) = fnm(:)

        f(:) = gzzsup(:,m)*rjacob(:,m)*(jfcn(:,m)**2)
        call toFourier('c')
        gzz_j(:,m) = fnm(:)

        f(:) = gzzsup(:,m)*(djdz(:,m)**2)*rjacob(:,m)
        call toFourier('c')
        gzz_ji_djdz2(:,m) = fnm(:)

        f(:) = gszsup(:,m)*djdz(:,m)*jfcn(:,m)*rjacob(:,m)
        call toFourier('c')
        gsz_djdz(:,m) = fnm(:)

        f(:) = gszsup(:,m)*djds(:,m)*djdz(:,m)*rjacob(:,m)
        call toFourier('c')
        gsz_ji_djds_djdz(:,m) = fnm(:)

        f(:) = gszsup(:,m)*djds(:,m)*jfcn(:,m)*rjacob(:,m)
        call toFourier('s')
        gsz_djds(:,m) = fnm(:)

        f(:) = gszsup(:,m)*rjacob(:,m)*(jfcn(:,m)**2)
        call toFourier('s')
        gsz_j(:,m) = fnm(:)

        f(:) = gtzsup(:,m)*rjacob(:,m)*(jfcn(:,m)**2)
        call toFourier('c')
        gzt_j(:,m) = fnm(:)

        f(:) = gtzsup(:,m)*djdt(:,m)*jfcn(:,m)*rjacob(:,m)
        call toFourier('s')
        gzt_djdt(:,m) = fnm(:)

        f(:) = gtzsup(:,m)*djdz(:,m)*jfcn(:,m)*rjacob(:,m)
        call toFourier('s')
        gzt_djdz(:,m) = fnm(:)

        f(:) = gtzsup(:,m)*djdt(:,m)*djdz(:,m)*rjacob(:,m)
        call toFourier('c')
        gzt_ji_djdz_djdt(:,m) = fnm(:)

      end do   ! m = 1,ns   
     
      tsec = secnds(t0)
      t0 = secnds(0.0)
      write(*,'("time bending FT(sec) = ",e12.5)') tsec
c
c
c    (3) Metric FT arrays for J_parallel term:
c        (build the required J|| coefficients
c         in real space and then transform to Fourier
c         space using the toFourier subroutine) 
c
       allocate(brho_ft(mnmx,ns), d_brho_dth_ft(mnmx,ns),
     1  d_brho_dzt_ft(mnmx,ns), stat=istat)
       allocate(jprl_over_b_ft(mnmx,ns),
     1  jprl_over_b12_ft(mnmx,ns), stat=istat)
       allocate(jprl_over_b_re(nznt,ns),
     1  jprl_over_b12_re(nznt,ns), stat=istat)
       allocate(jprl_ovr_b_brho_ft(mnmx,ns),
     1  jprl_ovr_b_dbrho_dtht_ft(mnmx,ns),
     2  jprl_ovr_b_dbrho_dzet_ft(mnmx,ns), stat=istat)
       allocate(d_brho_dth_re(nznt,ns),
     1  d_brho_dzt_re(nznt,ns), stat=istat)
       allocate(alf_tht(ns), alf_zet(ns), stat=istat)
      do m = 1,ns
       f(:) = brho(:,m)
       call toFourier('s')
       brho_ft(:,m) = fnm(:)
       
       call dbydth('s')
       
       d_brho_dth_ft(:,m) = anm(:)
       call dbydzt('s')
       d_brho_dzt_ft(:,m) = anm(:)
      end do   ! m = 1,ns
c
c   Form the part of Jprl/B that can be done
c    in Fourier space
c      
      do m = 1,ns
       jprl_over_b12_ft(:,m) =
     1      jprl_coef1(m)*d_brho_dzt_ft(:,m)
     2    + jprl_coef2(m)*d_brho_dth_ft(:,m)
      end do
c
c   Take Jprl/B to real space
c
      do m = 1,ns
       anm(:) = jprl_over_b12_ft(:,m)
       call toReal('c')
       jprl_over_b12_re(:,m) = f(:)
      end do
c
c   Add flux surface term that could be
c   directly added to Fourier space version
c
      do m = 1,ns
       jprl_over_b_re(:,m) = jprl_coef0(m)
     1   + jprl_over_b12_re(:,m)
      end do
c
c   Take total Jprl/B back to Fourier space
c
      do m = 1,ns
       f(:) = jprl_over_b_re(:,m)
       call toFourier('c')
       jprl_over_b_ft(:,m) = fnm(:)
      end do
c
c   Calculate various combinations of Jprl/B and Brho that are needed
c
      do m = 1,ns
       f(:) = jprl_over_b_re(:,m)*brho(:,m)
       call toFourier('s')
       jprl_ovr_b_brho_ft(:,m) = fnm(:)
      end do
c
      do m = 1,ns
       anm(:) = d_brho_dth_ft(:,m)
       call toReal('c')
       d_brho_dth_re(:,m) = f(:)
      end do
      do m = 1,ns
       f(:) = jprl_over_b_re(:,m)*d_brho_dth_re(:,m)
       call toFourier('c')
       jprl_ovr_b_dbrho_dtht_ft(:,m) = fnm(:)
      end do
c
      do m = 1,ns
       anm(:) = d_brho_dzt_ft(:,m)
       call toReal('c')
       d_brho_dzt_re(:,m) = f(:)
      end do
      do m = 1,ns
       f(:) = jprl_over_b_re(:,m)*d_brho_dzt_re(:,m)
       call toFourier('c')
       jprl_ovr_b_dbrho_dzet_ft(:,m) = fnm(:)
      end do
c
      do m = 1,ns
       if(lrfp) then
        alf_zet(m) = (1./iotaf(m))/(-jpolf(m)/iotaf(m) + jtorf(m))
        alf_tht(m) = 1./(-jpolf(m)/iotaf(m) + jtorf(m))
       else
        alf_zet(m) = 1./(-jpolf(m) + iotaf(m)*jtorf(m))
        alf_tht(m) = iotaf(m)/(-jpolf(m) + iotaf(m)*jtorf(m))
       endif
      end do     

      tsec = secnds(t0)
      t0 = secnds(0.0)
      write(*,'("time Jparallel FT(sec) = ",e12.5)') tsec
      
      
      deallocate(bfield, rjacob,gsssup, gttsup,gzzsup,
     1  gstsup, gszsup,gtzsup,bfields,bfieldth, bfieldze,
     2  djds,djdt,djdz,stat=istat)
      deallocate(d_brho_dth_re,
     1  d_brho_dzt_re,stat=istat)
      call trig_deallocate
     
c
c     Form Fourier mode sub-block matrices for the
c      inertia term
c
      allocate(a(mn_col,mn_col,ns), stat=istat)
      allocate(b(mn_col,mn_col,ns), stat=istat)
      allocate(d(mn_col,mn_col,ns), stat=istat)
      allocate(y_phi(ns), stat=istat)
      allocate(a_mag(ns), stat=istat)
      allocate(b_mag(ns), stat=istat)
      allocate(d_mag(ns), stat=istat)
      do m = 1,ns
       do i = 1,mn_col
        mi = im_col(i); ni = in_col(i)
        do j = 1,mn_col
         mj = im_col(j); nj = in_col(j)
         a(i,j,m) = 0.d0; b(i,j,m) = 0.d0; d(i,j,m) = 0.d0
         do k = 1,mnmx
          meq = im(k); neq = in(k)
          call ccc_convolve(ans_ccc,mi,ni,mj,nj,meq,neq)
          a(i,j,m) = a(i,j,m) + ans_ccc*gsssup_ft(k,m)
          call scs_convolve(ans_ssc,meq,neq,mi,ni,mj,nj)
          b(i,j,m) = b(i,j,m) - dble(mi)*ans_ssc*gstsup_ft(k,m)
     >               + dble(ni)*ans_ssc*gszsup_ft(k,m)
          call scs_convolve(ans_scs,mi,ni,mj,nj,meq,neq)
	  if(lrfp) then
           stell_rfp_specific_term1 = ans_scs*grad_to_prp_ft(k,m)    !grad to grad_perp term
     >               *(dble(nj)/iotaf(m) - dble(mj))    !grad to grad_perp term
     >               *(dble(ni)/iotaf(m) - dble(mi))    !grad to grad_perp term
	  else
           stell_rfp_specific_term1 = ans_scs*grad_to_prp_ft(k,m)       !grad to grad_perp term
     >               *(dble(ni)*dble(nj)-iotaf(m)*(dble(ni)*dble(mj)    !grad to grad_perp term
     >               + dble(nj)*dble(mi) - iotaf(m)*dble(mi)*dble(mj))) !grad to grad_perp term
	  endif
          d(i,j,m) = d(i,j,m) + (dble(mi)*dble(mj)*gttsup_ft(k,m)
     >               - dble(ni)*dble(mj)*gtzsup_ft(k,m)
     >               - dble(nj)*dble(mi)*gtzsup_ft(k,m)
     >               + dble(ni)*dble(nj)*gzzsup_ft(k,m))*ans_scs
c     >               + stell_rfp_specific_term1
         end do
c     if(m .eq. 1) write(*,99) mi, ni, mj, nj,
c     >         a(i,j,m),b(i,j,m),d(i,j,m)
        end do
       end do
      end do
 99    format(i4,3(1x,i4),3(2x,e12.4))
      tsec = secnds(t0)
      t0 = secnds(0.0)
      write(*,'("time inertia blocks(sec) = ",e12.5)') tsec
      deallocate(gsssup_ft, gttsup_ft,
     1  gzzsup_ft, gstsup_ft,gszsup_ft, gtzsup_ft,
     2  stat=istat)
c
c     Form Fourier mode sub-block matrices for the
c      bending energy term
c
      allocate(aa(mn_col,mn_col,ns), stat=istat)
      allocate(bb(mn_col,mn_col,ns), stat=istat)
      allocate(bb1(mn_col,mn_col,ns), stat=istat)
      allocate(bb2(mn_col,mn_col,ns), stat=istat)
      allocate(cc(mn_col,mn_col,ns), stat=istat)
      allocate(cc2(mn_col,mn_col,ns), stat=istat)
      allocate(dd(mn_col,mn_col,ns), stat=istat)
      allocate(dd1(mn_col,mn_col,ns), stat=istat)
      allocate(dd2(mn_col,mn_col,ns), stat=istat)
      allocate(aa_mag(ns), stat=istat)
      allocate(bb_mag(ns), stat=istat)
      allocate(dd1_mag(ns), stat=istat)
      allocate(dd2_mag(ns), stat=istat)
      do m = 1,ns
       do i = 1,mn_col
        mi = im_col(i); ni = in_col(i)
        do j = 1,mn_col
         mj = im_col(j); nj = in_col(j)
         aa(i,j,m) = 0.d0; bb(i,j,m) = 0.d0
         dd1(i,j,m) = 0.d0; dd2(i,j,m) = 0.d0
         do k = 1,mnmx
           meq = im(k); neq = in(k)
           rmi = dble(mi); rni = dble(ni)
           rmj = dble(mj); rnj = dble(nj)
           rme = dble(meq); rne = dble(neq)
	   if(lrfp) then
            kprl_i = rni/iotaf(m) - rmi
	    kprl_j = rnj/iotaf(m) - rmj
	   else
            kprl_i = rni - iotaf(m)*rmi
	    kprl_j = rnj - iotaf(m)*rmj
	   endif
	   kprl_i_nw = alf_zet(m)*rni - alf_tht(m)*rmi
	   kprl_j_nw = alf_zet(m)*rnj - alf_tht(m)*rmj
	   d_kprl_i_nw_dtht = rmi*rni*alf_zet(m) - rmi*rmi*alf_tht(m)
	   d_kprl_i_nw_dzet = rmi*rni*alf_tht(m) - rni*rni*alf_zet(m)
           kprl_ij = kprl_i*kprl_j
           call scs_convolve(ans_ssc,meq,neq,mi,ni,mj,nj)  !eq=sin, i=cos, j=sin
           call scs_convolve(ans_css,meq,neq,mj,nj,mi,ni)  !eq=sin, i=sin, j=cos
           call scs_convolve(ans_scs,mi,ni,mj,nj,meq,neq)  !eq=cos, i=sin, j=sin
           call ccc_convolve(ans_ccc,mi,ni,mj,nj,meq,neq)  !eq=cos, i=cos, j=cos
	  if(lrfp) then
           stell_rfp_specific_term2 = 
     >         -ans_scs*gss_j(k,m)*rnj*kprl_i*iota_inv_pf(m)
           stell_rfp_specific_term3 = 
     >          ans_scs*(gss_ji_djds2(k,m)*kprl_ij
     >          - gss_djds(k,m)*iota_inv_pf(m)*(kprl_i*rnj - kprl_j*rni)
     >          + rni*rnj*gss_j(k,m)*(iota_inv_pf(m)**2))
           stell_rfp_specific_term4 = 
     >          -ans_scs*gst_djdt(k,m)*kprl_j*rni*iota_inv_pf(m)
           stell_rfp_specific_term5 = 
     >          -ans_css*gst_j(k,m)*rni*rmj*kprl_j*iota_inv_pf(m)
           stell_rfp_specific_term6 = 
     >          -ans_scs*gsz_djdz(k,m)*kprl_j*rni*iota_inv_pf(m)
           stell_rfp_specific_term7 = 
     >           ans_css*gsz_j(k,m)*kprl_j*rnj*rni*iota_inv_pf(m)
	  else
           stell_rfp_specific_term2 =
     >         -ans_scs*gss_j(k,m)*rmj*iotapf(m)*kprl_i
           stell_rfp_specific_term3 =
     >          ans_scs*(gss_ji_djds2(k,m)*kprl_ij
     >          - gss_djds(k,m)*iotapf(m)*(kprl_i*rmj - kprl_j*rmi)
     >          + rmi*rmj*gss_j(k,m)*iotapf(m)*iotapf(m))
           stell_rfp_specific_term4 =
     >          -ans_scs*gst_djdt(k,m)*kprl_j*rmi*iotapf(m)
           stell_rfp_specific_term5 =
     >          -ans_css*gst_j(k,m)*rmi*rmj*kprl_j*iotapf(m)
           stell_rfp_specific_term6 =
     >          -ans_scs*gsz_djdz(k,m)*kprl_j*rmi*iotapf(m)
           stell_rfp_specific_term7 =
     >           ans_css*gsz_j(k,m)*kprl_j*rnj*rmi*iotapf(m)
	  endif
           aa(i,j,m) = aa(i,j,m) + ans_scs*gss_j(k,m)*kprl_ij

           bb1(i,j,m) = bb1(i,j,m)
     1           + ans_scs*gss_djds(k,m)*kprl_ij
     3           + ans_scs*gst_djdt(k,m)*kprl_ij
     4           + ans_css*gst_j(k,m)*rmj*kprl_ij
     5           + ans_scs*gsz_djdz(k,m)*kprl_ij
     1           
           bb2(i,j,m) = bb2(i,j,m) + stell_rfp_specific_term2
     6           - ans_css*gsz_j(k,m)*rnj*kprl_j
     8           - ans_css*(rmi*pres_a_zet_ft(k,m)               !new sym pressure term
     9           + rni*pres_a_tht_ft(k,m))                       !new sym pressure term
     1           - 0.5*ans_ccc*jprl_over_b_ft(k,m)
     2             *(jpolf(m)*d_kprl_i_nw_dtht
     3             + jtorf(m)*d_kprl_i_nw_dzet)
           cc2(i,j,m) = cc2(i,j,m)
     1            - 0.5*ans_scs*jprl_over_b_ft(k,m)*jpolf(m)*kprl_i_nw*
     2               (jtorf(m)*kprl_j_nw + rmj)
     3            + 0.5*ans_scs*jprl_over_b_ft(k,m)*jtorf(m)*kprl_i_nw*
     4               (rnj + jpolf(m)*kprl_j_nw)
           dd1(i,j,m) = dd1(i,j,m) + stell_rfp_specific_term3
     4          + ans_scs*gtt_ji_djdt2(k,m)*kprl_ij
     5          + ans_css*gtt_djdt(k,m)*kprl_ij*rmj
     6          + ans_ssc*gtt_djdt(k,m)*kprl_ij*rmi
     7          + ans_ccc*gtt_j(k,m)*kprl_ij*rmi*rmj
     8          + ans_scs*gzz_ji_djdz2(k,m)*kprl_ij
     9          - ans_css*gzz_djdz(k,m)*rnj*kprl_ij
     1          - ans_ssc*gzz_djdz(k,m)*rni*kprl_ij
     2          + ans_ccc*gzz_j(k,m)*rni*rnj*kprl_ij
c     4          - ans_scs*bgradb2_ft(k,m)*kprl_i*kprl_j           !grad to grad_perp term
c     5          + ans_ssc*b2gradb_ft(k,m)*kprl_j*kprl_i*kprl_i    !grad to grad_perp term
c     6          + ans_css*b2gradb_ft(k,m)*kprl_i*kprl_j*kprl_j    !grad to grad_perp term
c     7          - ans_ccc*b3_ft(k,m)*((kprl_i*kprl_j)**2)         !grad to grad_perp term
           dd2(i,j,m) = dd2(i,j,m)
     1          + ans_scs*gst_ji_djds_djdt(k,m)*kprl_ij
     2          + ans_css*gst_djds(k,m)*rmj*kprl_ij
     >          + stell_rfp_specific_term4
     3          + stell_rfp_specific_term5          
     5          + ans_scs*gsz_ji_djds_djdz(k,m)*kprl_ij
     6          - ans_css*gsz_djds(k,m)*kprl_ij*rnj
c     >          + stell_rfp_specific_term6
c     7          + stell_rfp_specific_term7
     9          + ans_scs*gzt_ji_djdz_djdt(k,m)*kprl_ij
     1          - ans_css*gzt_djdt(k,m)*rnj*kprl_ij
     2          + ans_ssc*gzt_djdz(k,m)*rmi*kprl_ij
     3          - ans_ccc*gzt_j(k,m)*rmi*rnj*kprl_ij
     4          + ans_scs*(rmi*rmj*pres_c_zet_ft(k,m)         !new sym pressure term
     5          + rni*rmj*pres_c_tht_ft(k,m)                  !new sym pressure term
     6          - rmi*rnj*pres_d_zet_ft(k,m)                  !new sym pressure term
     7          + rni*rnj*pres_d_tht_ft(k,m))                 !new sym pressure term
     4           - 0.5*ans_ssc*jprl_ovr_b_brho_ft(k,m)
     5             *(jpolf(m)*d_kprl_i_nw_dtht
     6             + jtorf(m)*d_kprl_i_nw_dzet)*kprl_j_nw
     7           - 0.5*ans_ssc*jprl_ovr_b_brho_ft(k,m)
     8             *((jpolf(m)*d_kprl_i_nw_dtht
     9             + jtorf(m)*d_kprl_i_nw_dzet)*kprl_j_nw
     1             + d_kprl_i_nw_dzet*(jtorf(m)*kprl_j_nw + rmj)
     2             + d_kprl_i_nw_dtht*(rnj + jpolf(m)*kprl_j_nw))
     3           - 0.5*ans_scs*jprl_ovr_b_dbrho_dzet_ft(k,m)*kprl_i_nw
     4             *(jtorf(m)*kprl_j_nw + rmj)
     5           - 0.5*ans_scs*jprl_ovr_b_dbrho_dtht_ft(k,m)*kprl_i_nw
     6             *(rnj + jpolf(m)*kprl_j_nw)
     7           - 0.5*ans_scs*jprl_over_b_ft(k,m)*jpolpf(m)*kprl_i_nw
     8              *(jtorf(m)*kprl_j_nw + rmj)
     9           + 0.5*ans_scs*jprl_over_b_ft(k,m)*jtorpf(m)*kprl_i_nw
     1              *(rnj + jpolf(m)*kprl_j_nw)
       
          end do
         end do
        end do
       end do
      do m = 1,ns
       do i = 1,mn_col
        do j = 1,mn_col
         bb(i,j,m) = bb1(i,j,m) + bb2(i,j,m)
         cc(i,j,m) = cc2(i,j,m)
         dd(i,j,m) = dd1(i,j,m) + dd2(i,j,m) + dd2(j,i,m)
        end do
       end do
      end do
      tsec = secnds(t0)
      t0 = secnds(0.0)
      write(*,'("time bending blocks(sec) = ",e12.5)') tsec
c
c   collect diagnostic info. on the block sizes
c 
      open(unit=44,file="ae_diag.dat",status="unknown")     
      do m = 1,ns
       a_mag(m) = 0.;b_mag(m) = 0.;d_mag(m) = 0.
       aa_mag(m) = 0.;bb_mag(m) = 0.;dd1_mag(m) = 0.;dd2_mag(m) = 0.
       do i = 1,mn_col
        do j = 1,mn_col
	 a_mag(m) = a_mag(m) + a(i,j,m)**2
	 b_mag(m) = b_mag(m) + b(i,j,m)**2
	 d_mag(m) = d_mag(m) + d(i,j,m)**2
	 aa_mag(m) = aa_mag(m) + aa(i,j,m)**2
	 bb_mag(m) = bb_mag(m) + bb(i,j,m)**2
	 dd1_mag(m) = dd1_mag(m) + dd1(i,j,m)**2
	 dd2_mag(m) = dd2_mag(m) + dd2(i,j,m)**2
        end do
       end do
       write(44,'(i3,7(3x,e15.7))') m,a_mag(m),b_mag(m),d_mag(m),
     1       aa_mag(m),bb_mag(m),dd1_mag(m),dd2_mag(m)
      end do
      close(unit=44)
      
      deallocate(gss_ji_djds2,gss_j,gss_djds,gst_djdt,     
     1  gst_j,gst_djds,gst_ji_djds_djdt,     
     2  gtt_djdt,gtt_j,gtt_ji_djdt2,
     3  gzz_djdz, gzz_ji_djdz2,   
     4  gzz_j,gsz_djdz,gsz_ji_djds_djdz,     
     5  gsz_djds,gsz_j,gzt_j,gzt_djdt,gzt_djdz,
     6  gzt_ji_djdz_djdt, stat=istat)
      deallocate(brho_ft, d_brho_dth_ft, d_brho_dzt_ft,
     1  jprl_over_b_ft, jprl_over_b_re,
     3  jprl_over_b12_ft, jprl_over_b12_re, stat=istat)
      deallocate(jprl_ovr_b_brho_ft,
     1  jprl_ovr_b_dbrho_dtht_ft,
     2  jprl_ovr_b_dbrho_dzet_ft, stat=istat)
      deallocate(bbinv2_ft, brho_binv2_ft, stat=istat)
c
c      setup indirect addressing arrays:
c
      allocate(mmx(mn_col*ns,mn_col*ns), stat=istat)
      allocate(nnx(mn_col*ns,mn_col*ns), stat=istat)
      allocate(ix(mn_col*ns,mn_col*ns), stat=istat)
      allocate(jx(mn_col*ns,mn_col*ns), stat=istat)
      if(.not.build_matrix_only) then
       allocate(mat_inertia(mn_col*ns,mn_col*ns), stat=istat)
       allocate(mat_bending(mn_col*ns,mn_col*ns), stat=istat)
       mat_inertia(:,:) = 0.d0; mat_bending(:,:) = 0.d0
      end if
      ii = 0
      do mm = 1,ns
       do i = 1,mn_col
        mi = im_col(i)
        ii = ii + 1
        jj = 0
        do nn = 1,ns
         do j = 1,mn_col
          mj = im_col(j)
          jj = jj + 1
	  mmx(ii,jj) = mm; nnx(ii,jj) = nn
	  ix(ii,jj) = i; jx(ii,jj) = j
          call element(mat_inertia_ij,mat_bending_ij,mm,i,nn,j,
     >           ns, mi, mj)
c   Direct incorporation of boundary conditions:       NIFS test
c             if(mm .eq. 1 .and. nn .eq. 1) then
c	       if(i .eq. j) mat_inertia_ij = 1.
c	       if(i .ne. j) mat_inertia_ij = 0.
c	       mat_bending_ij = 0.
c             else if(mm .eq. 1 .and. nn .eq. 2) then
c	        mat_inertia_ij = 0.
c	        mat_bending_ij = 0.
c		if(i .eq. j .and. mi .eq. 0) mat_inertia_ij = -1.
c             else if(mm .eq. ns .and. nn .eq. ns) then
c	       if(i .eq. j) mat_inertia_ij = 1.
c	       if(i .ne. j) mat_inertia_ij = 0.
c	       mat_bending_ij = 0.
c             else if(mm .eq. ns .and. nn .eq. ns-1) then
c	        mat_inertia_ij = 0.
c	        mat_bending_ij = 0.
c	     end if
c       NIFS test
             if(.not.build_matrix_only) then
              mat_inertia(ii,jj) = -mat_inertia_ij
              mat_bending(ii,jj) = -mat_bending_ij
	     else if(build_matrix_only) then
	      if(mat_inertia_ij .ne. 0.d0) then
	       write(36,'(i5,2x,i5,2x,e15.7)') ii,jj,-mat_inertia_ij
	      end if
	      if(mat_bending_ij .ne. 0.d0) then
	       write(35,'(i5,2x,i5,2x,e15.7)') ii,jj,-mat_bending_ij
	      end if
	    end if
         end do
        end do
       end do
      end do
      tsec = secnds(t0)
      t0 = secnds(0.0)
      write(*,'("time full matrix assembly(sec) = ",e12.5)') tsec
c
      if(build_matrix_only) then
       write(32,*) ns, mn_col
       do i=1,mn_col
        write(32,*) im_col(i),in_col(i)
       end do
       do j=1,ns
        write(32,*) rho(j)
       end do
       close(unit=32)
       close(unit=35)
       close(unit=36)
       stop 25
      end if
c      
      nmat = mn_col*ns
      if(jdqz_data) then
       write(32,*) ns, mn_col
       do i=1,mn_col
        write(32,*) im_col(i),in_col(i)
       end do
       do j=1,ns
        write(32,*) rho(j)
       end do
        do i=1,nmat
	 do j=1,nmat
	  if(mat_inertia(i,j) .ne. 0.d0) then
	   write(36,'(i5,2x,i5,2x,e15.7)') i, j, mat_inertia(i,j)
	  end if
	  if(mat_bending(i,j) .ne. 0.d0) then
	   write(35,'(i5,2x,i5,2x,e15.7)') i, j, mat_bending(i,j)
	  end if
	 end do
	end do
	close(unit=32)
	close(unit=35)
	close(unit=36)
      end if
c
       allocate(ctmp1(nmat,nmat), stat=istat)
       allocate(ctmp2(nmat,nmat), stat=istat)
c     
c
c
c    Solve A*z = lambda*B*z for eigenvalues (lambda) and eigenvectors (z)
c     where A = amat = mat_bending, B = bmat = mat_inertia.
c
      ku1 = ku+1
c    Allocate work arrays for Lapack routines
      lwork = 80*nmat*mn_col
      allocate(work(lwork),stat=istat)
      if(istat .ne. 0) write(*,'("work allocation failed")')
      allocate(dm(nmat),stat=istat)
      if(istat .ne. 0) write(*,'("dm allocation failed")')
      allocate(zz(nmat,nmat), stat=istat)
      allocate(alphar(nmat), stat=istat)
      allocate(alphai(nmat), stat=istat)
      allocate(betar(nmat), stat=istat)
      allocate(vr(nmat,nmat), stat=istat)
      allocate(vl(nmat,1), stat=istat)
c     Symmetry unconstrained:
c       ctmp1(:,:) = mat_bending(:,:)
c       ctmp2(:,:) = mat_inertia(:,:)
c     Force symmetry:
       do i = 1,nmat
        do j = 1,nmat
	  ctmp1(i,j) = 0.5d0*(mat_bending(i,j) + mat_bending(j,i))
	  ctmp2(i,j) = 0.5d0*(mat_inertia(i,j) + mat_inertia(i,j))
	end do
       end do
       tsec0 = secnds(t00)
       write(*,'("just before dggev:",e12.4," sec")') tsec0
       t00 = secnds(0.0)
       call dggev('N','V',nmat,ctmp1,nmat,ctmp2,nmat,
     >   alphar,alphai,betar,vl,nmat,vr,nmat,work,lwork,info)
       tsec0 = secnds(t00)
       write(*,'("just after dggev:",e12.4," sec")') tsec0
       if(info .ne. 0) write(*,'("info dggev = ",i9)') info
c       write(18,*) nmat
c       do i=1,nmat
c        write(18,'(e15.7,2(2x,e15.7))') alphar(i), alphai(i), betar(i)
c       end do
       
      
c
c    Calculate inertial (electrostatic) and bending (electromagnetic)
c     energies for each eigenmode:
c
       allocate(energy_inertia(nmat), stat=istat)
       allocate(energy_bending(nmat), stat=istat)
       allocate(energy_inertia_1(nmat), stat=istat)
       allocate(energy_bending_1(nmat), stat=istat)
       one = 1.d0; zero = 0.d0
       ctmp1(:,:) = 0.d0; ctmp2(:,:) = 0.d0
        call dgemm('N','N',nmat,nmat,nmat,one,mat_inertia,
     1   nmat,vr,nmat,zero,ctmp1,nmat)
        call dgemm('T','N',nmat,nmat,nmat,one,vr,
     1   nmat,ctmp1,nmat,zero,ctmp2,nmat)
       do i=1,nmat
        energy_inertia(i) = ctmp2(i,i)
       end do
       one = 1.d0; zero = 0.d0
       ctmp1(:,:) = 0.d0; ctmp2(:,:) = 0.d0
        call dgemm('N','N',nmat,nmat,nmat,one,mat_bending,
     1    nmat,vr,nmat,zero,ctmp1,nmat)
        call dgemm('T','N',nmat,nmat,nmat,one,vr,
     1    nmat,ctmp1,nmat,zero,ctmp2,nmat)
       do i=1,nmat
        energy_bending(i) = ctmp2(i,i)
       end do
       deallocate(ctmp1, ctmp2, stat=istat)
c
       icount = 0
       do i=1,nmat
c        if(alphai(i) .ne. 0.d0) cycle   !these lines commented out so egn_value file
c	if(betar(i) .ne. 0.d0) then      !has nmat lines instead of .le. nmat
	 icount = icount + 1
c	 if(betar(i) .ne. 0.d0 .and. alphar(i) .ge. 0.) then
	 if(betar(i) .ne. 0.d0) then
	  dm(icount) = alphar(i)/betar(i)
         energy_inertia_1(icount) =
     >     energy_inertia(i)/(mu0*1.d+6*(two_pi**2))
         energy_bending_1(icount) =
     >     energy_bending(i)/(mu0*dm(icount)*1.d+6*(two_pi**2))
	 else
	  dm(icount) = 1.e+10
         energy_inertia_1(icount) =  1.e+10
         energy_bending_1(icount) =  1.e+10
	 end if
	 write(47,'(i8,4(2x,e15.7))') i,alphar(i),alphai(i),
     >       betar(i),dm(i)
	 zz(:,icount) = vr(:,i)
c	end if
       end do
       nmat = icount
c
c
      do i=1,nmat
       if(dm(i) .ge. 0.) then
        write(18,87)sqrt(dm(i)),energy_inertia_1(i),energy_bending_1(i)
       else if(dm(i) .lt. 0.) then
        write(18,87) dm(i),energy_inertia_1(i),energy_bending_1(i)
       end if
      end do
  87  format(e15.7,2(2x,e15.7))      
c
c      
      allocate(egn_vectors(nmat,mn_col,ns), stat=istat)
      do j = 1,nmat
       jj = 0
       do m = 1,ns
        do i=1,mn_col
         jj = jj + 1
          egn_vectors(j,i,m) = zz(jj,j)
        end do
       end do
      end do
      write(*,*) mn_col, ns
c      jj = 64
c      do m = 1,ns
c       write(20,23) rho(m), (egn_vectors(jj,i,m), i=1,mn_col)
c      end do
  23  format(37(1x,e12.5))
      if(egnout_form .eq. "binr") then
       write(33) nmat
       write(33) mn_col
       write(33) ns
       write(33) im_col
       write(33) in_col
       write(33) dm
       write(33) rho
       write(33) egn_vectors
      else if(egnout_form .eq. "asci") then
       write(33,*) nmat
       write(33,*) mn_col
       write(33,*) ns
       do i=1,mn_col
        write(33,*) im_col(i)
        write(33,*) in_col(i)
       end do
       write(*,*) im_col(mn_col),im_col(mn_col)
       do i=1,nmat
        write(33,'(e15.7)') dm(i)
       end do
       write(*,'(e15.7)') dm(nmat)
       do i=1,ns
        write(33,'(e15.7)') rho(i)
       end do
       write(*,'(e15.7)') rho(ns)
       do j=1,nmat
        do m = 1,ns
         do i=1,mn_col
          write(33,'(e15.7)') egn_vectors(j,i,m)
          end do
         end do
       end do
      endif
      
      if(residue_test .and. full_matrix_test) then
      
       allocate(lambda_i(nmat,nmat), stat=istat)
       lambda_i(:,:) = 0.d0
       do i=1,nmat
        lambda_i(i,i) = dm(i)
       end do
       allocate(temp1(nmat,nmat), stat=istat)
       temp1(:,:) = 0.d0
       temp1 = matmul(lambda_i,mat_inertia)
       deallocate(lambda_i, mat_inertia, stat=istat)
       allocate(temp2(nmat,nmat), stat=istat)
       temp2(:,:) = 0.d0
       temp2(:,:) = mat_bending(:,:) - temp1(:,:)
       deallocate(mat_bending, stat=istat)
       temp1 = matmul(temp2,zz)
       res = 0.d0; norm0 = 0.d0
       do i=1,nmat
        do j=1,nmat
	 norm0 = norm0 + abs(temp2(i,j))/dble(nmat*nmat)
	 res = res + abs(temp1(i,j))/dble(nmat*nmat)
	end do
       end do
       write(*,'("norm(A - lambda*B) = ",e15.7,/,
     1   "norm(Ay - lambda*By) = ",e15.7)') norm0, res

      else if(residue_test .and. band_matrix_test) then
 
       allocate(egn_test(nmat), stat=istat)
       allocate(egn_test1(nmat), stat=istat)
       do i=1,nmat
        do j=1,nmat
	 egn_test(j) = zz(j,i)
	end do
	alpha = 1.d0; beta = 0.d0; incx = 1; incy = 1
	call dsbmv(uplo,nmat,ku,alpha,mat_bending,ku1,
     1         egn_test,incx,beta,egn_test1,incy)
	incx = 1; incy = 1
        res_a = ddot(nmat,egn_test1,incx,egn_test1,incy)
        res_a = sqrt(res_a)/dble(nmat)
	alpha = -dm(i); beta = 1.d0; incx = 1; incy = 1
	call dsbmv(uplo,nmat,ku,alpha,mat_inertia,ku1,egn_test,
     1         incx,beta,egn_test1,incy)
	incx = 1; incy = 1
	res = ddot(nmat,egn_test1,incx,egn_test1,incy)
        res = sqrt(res)/dble(nmat)
        write(22,*) res_a, res
       end do
      end if   !if(residue_test .and. full(or)band_matrix_test)
      
      end program ae_eigensolver
c
c
c
      subroutine element(xmat1,xmat2,mm,i1,nn,j1,ns,mr,mc)
      use kind_spec
      use fourier_lib
      implicit none
      real(kind=rprec) :: xmat1, xmat2, rho_int_2drv,
     >  rho_int_1drv_a, rho_int_1drv_b, rho_int0,
     >  rho_lwr_mm, rho_upr_mm, rho_lwr_nn, rho_upr_nn
      integer :: mm, nn, i1, j1, ns, mr, mc
      if(mm .eq. 1) rho_lwr_mm = 0.d0
      if(nn .eq. 1) rho_lwr_nn = 0.d0
      if(mm .eq. ns) rho_upr_mm = 1.d0
      if(nn .eq. ns) rho_upr_nn = 1.d0
      if(mm .gt. 1) rho_lwr_mm = rho(mm-1)
      if(mm .lt. ns) rho_upr_mm = rho(mm+1)
      if(nn .gt. 1) rho_lwr_nn = rho(nn-1)
      if(nn .lt. ns) rho_upr_nn = rho(nn+1)
      if(nn .eq. mm) then
       rho_int_2drv =
     >     1.d0/(rho(mm) - rho_lwr_mm)
     >   + 1.d0/(rho_upr_mm - rho(mm))
       rho_int_1drv_a = 0.d0
       rho_int_1drv_b = 0.d0
       rho_int0 = (rho_upr_mm - rho_lwr_mm)/3.d0
      else if(nn .eq. (mm-1)) then
       rho_int_2drv =
     >   -1.d0/(rho(mm) - rho_lwr_mm)
       rho_int_1drv_a = -0.5d0
       rho_int_1drv_b = 0.5d0
       rho_int0 = (rho(mm) - rho_lwr_mm)/6.d0
      else if(nn .eq. (mm+1)) then
       rho_int_2drv =
     >   -1.d0/(rho_upr_mm - rho(mm))
       rho_int_1drv_a = 0.5d0
       rho_int_1drv_b = -0.5d0
       rho_int0 = (rho_upr_mm - rho(mm))/6.d0     

   
c      else if(mm .eq. (nn-1)) then
c       rho_int_2drv =
c     >   -1.d0/(rho(nn) - rho_lwr_nn)
c       rho_int_1drv_a = 0.5d0
c       rho_int_1drv_b = -0.5d0
c       rho_int0 = (rho(nn) - rho_lwr_nn)/6.d0
c      else if(mm .eq. (nn+1)) then
c       rho_int_2drv =
c     >   -1.d0/(rho_upr_nn - rho(nn))
c       rho_int_1drv_a = -0.5d0
c       rho_int_1drv_b = 0.5d0
c       rho_int0 = (rho_upr_nn - rho(nn))/6.d0       
      else
       rho_int_2drv = 0.d0; rho_int0 = 0.d0
       rho_int_1drv_a = 0.d0; rho_int_1drv_b = 0.d0
      end if
c
c   B.c.'s for m = 0 modes at magnetic axis:
c
      if(mr .eq. 0 .or. mc .eq. 0) then
       if(mm .eq. 1 .and. nn .eq. 1) then
        rho_int0 = rho(1) + (rho(2)-rho(1))/3.d0
	rho_int_1drv_a = 0.5d0*(rho(1)+rho(2))/(rho(1)-rho(2))
	rho_int_1drv_b = rho_int_1drv_a
	rho_int_2drv = 1.d0/(rho(2) - rho(1))
       end if
      end if
c
      xmat1 = ((a(i1,j1,mm)+a(i1,j1,nn))*rho_int_2drv
     >      + (b(i1,j1,mm)+b(i1,j1,nn))*rho_int_1drv_a
     >      + (b(j1,i1,mm)+b(j1,i1,nn))*rho_int_1drv_b
     >      + (d(i1,j1,mm)+d(i1,j1,nn))*rho_int0)/2.d0
      xmat2 = ((aa(i1,j1,mm)+aa(i1,j1,nn))*rho_int_2drv
     >      + (bb(i1,j1,mm)+bb(i1,j1,nn))*rho_int_1drv_a
     >      + (bb(j1,i1,mm)+bb(j1,i1,nn))*rho_int_1drv_b
     >      + (cc(i1,j1,mm)+cc(i1,j1,nn))*rho_int_1drv_b
     >      + (cc(j1,i1,mm)+cc(j1,i1,nn))*rho_int_1drv_a
     >      + (dd(i1,j1,mm)+dd(i1,j1,nn))*rho_int0)/2.d0
     

c     
c      xmat1 = a(i1,j1,mm)*rho_int_2drv
c     >      + b(i1,j1,mm)*rho_int_1drv_a
c     >      + b(j1,i1,mm)*rho_int_1drv_b
c     >      + d(i1,j1,mm)*rho_int0
c      xmat2 = aa(i1,j1,mm)*rho_int_2drv
c     >      + bb(i1,j1,mm)*rho_int_1drv_a
c     >      + bb(j1,i1,mm)*rho_int_1drv_b
c     >      + dd(i1,j1,mm)*rho_int0
      return
      end
      
      
      subroutine linear_fe_eval(rho_eval,ns,phi_eval)
      use kind_spec
      use fourier_lib
      implicit none
      real(kind=rprec) :: rho_eval, phi_eval
      integer :: l, ns
      if(rho_eval .ge. 0.d0 .and. rho_eval .le. rho(1)) then
       phi_eval = rho_eval*y_phi(1)/rho(1)
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
      
