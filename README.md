# AE3D
AE3D calculates shear Alfvén eigenmodes for 3D configurations (stellarators, RFPs, 3D tokamaks)

These codes are used to calculate shear Alfven eigenmodes for 2D and 3D configurations, without sound wave coupling effects. The relevant publication for this code is: D. A. Spong, E. D’Azevedo, and Y. Todo, Phys. Plasmas 17, 022106 (2010).

The workflow is as follows:

(a) prepare a VMEC equilibrium for the case of interest

(b) Run xbooz_xform (from Stellopt code suite) to convert from VMEC coordinates to Boozer coordinates. This run needs the in_booz.* file, which tells it what range of m/n modes to use and the selected surfaces to use from the VMEC run. Typically, the first and last of the VMEC surfaces are removed, because they can sometimes have noisy data.

(c) Using the boozmn.* file produced in (b), run xmetric_ver* to extract needed data and calculate the metric elements used in the continuum calculation. This produces a file called ae_metric.dat which is used as input for the AE3D executable (xae3d_ver*).

(d) The eigenmode solver code is comprised of ae_solve_version*.f and fourier_lib_module.f. These are compiled and loaded using the scirpt bld_xae_verX, which one runs with the verison number of ae_solve.f added as a command line argument. The input files for xae3d are ae_metric.dat, fourier.dat and plasma.dat. fourier.dat and plasma.dat are the same as used for the continuum calculation by xstgap. However, plasma.dat will differ from the version used in the xstgap_snd (continuum with sound coupling) in that parameters relating to the electron temperature (last two lines of plasma.dat) should be removed for the xae3d run.

(e) xae3d solves for all of the eigenvalues and writes this data to the text files egn_mode_ascii.dat and egn_values.dat. The eigenmodes in these files can be plotted (labeled with their eigenfrequencies) by running the code xplt_egn_nw. This is based on the source code egn_plotter.f, which is compiled using bld_plt_egn_nw; this plotter requires the Plplot graphics library (http://plplot.sourceforge.net). egn_plotter.f is set up with some filters (based on the eigenfrequency) on which eigenmodes are displayed to limit the number of plots.

(f) For large problems the direct solution of all eigenmodes by xae3d can become quite time-consuming, so another option is provided. Xae3d writes out the non-zero elements of the matrices it generates in the text files a_matrix.dat and b_matrix.dat. The file jdqz.dat also provides information about the Fourier modes used and the radial grid. This data can be used as input to a Jacobi-Davidson solver, xjdqz, which solves for a more limited number (40) of eigenmodes near a specified target frequency (given as a command line argument to xjdqz). Xjdqz is compiled using the script bld_jdqz and requires a library libjdqz.a. The source to build this libary can be obtained from http://www.staff.science.uu.nl/~sleij101/JD_software/jd.html. Another option is to use the two matrix files with Jacobi-Davidson solvers in the SLEPC library (http://slepc.upv.es). Some limited testing has been done with the Slepc solvers, indicating agreement with results from xjdqz and xae3d solvers. SLEPC requires the installation of PETSC (must be the same versions), and has the advantgage over xjdqz that it is parallelized and can be used for very large eigenvalue problems.
