!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%% MODULES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module data_grid
 integer, parameter::      Nt =  2500     ! number of time steps
 integer, parameter::      Nx =  1024     ! grid points for x1 & x2
 integer, parameter::      NR =    64     ! grid points for R
 integer, parameter::      NE =   500     ! grid points for energy
 integer, parameter:: Nstates =     4     ! number of electronic states
 integer, parameter:: Vstates =     2     ! number of vibrational states
 integer         :: sym                   ! symmetry of spatial eigenfunctions in ITP
 integer		 :: init_state            ! initial properties: state
 double precision:: R_init, beta_R        ! initial properties: wave packet
 double precision:: beta_nucl, Req        ! initial properties: nuclear  distribution
 double precision:: beta_elec, xeq1, xeq2 ! initial properties: electron distribution
 double precision:: dx, dR                ! spacing of real-space
 double precision:: dpx, dpR              ! spacing of momenta
 double precision:: dE                    ! spacing of energy
 double precision:: dt, dt_adia, dt_vib   ! spacing of time (propagation, elec. ITP, nuclear ITP)
 double precision:: thresh                ! threshold value for ITPs
 double precision:: fwhm_ir, fwhm_xuv     ! FWHM of envelopes (IR & XUV pulse)
 double precision:: beta_ir, beta_xuv     ! handy variables for describing envelopes (derived from fwhm)
 double precision:: lambda_ir, lambda_xuv ! central pulse wavelength
 double precision:: omega_ir, omega_xuv   ! central pulse frequecies
 double precision:: t0_ir, t0_xuv         ! center of pulse envelopes
 double precision:: E0_ir, E0_xuv         ! field strength of IR and XUV pulses 
 double precision:: cep                   ! Carrier-Envelope-Phase of IR pulse (in units of PI)
end module

module data_au
 double precision, parameter::  au2A = 0.52917706d0        ! conversion of length (a.u. to Angstrom)
 double precision, parameter:: au2fs = 0.024d0             ! conversion of time (a.u. to fs)
 double precision, parameter:: cm2au = 4.5554927d-6        ! conversion of energy (wavenumber to a.u.)
 double precision, parameter::  j2eV = 6.242D18            ! conversion of energy (J to eV)
 double precision, parameter:: au2eV = 27.2116d0           ! conversion of energy (a.u. to eV)
 double precision, parameter::  i2au = 2.0997496D-9        ! ???
 double precision, parameter::    pi = 3.141592653589793d0 ! pi
 double precision, parameter::  mass = 1836.1526738917d0   ! proton mass
 double precision, parameter::    me = 1.d0                ! electron mass
 complex*16, parameter::          im = (0.d0,1.d0)         ! imaginary unit
end module

module pot_param
use data_au
 double precision::   rr    ! Position of right nucleus
 double precision::   rl    ! Position of left nucleus
 double precision::   rc    ! screening parameter: electron-(moving)-nucleus
 double precision::   rf    ! screening parameter: electron-(fixed)-nuclei
 double precision::   re    ! screening parameter: electron-electron
 double precision::   z1    ! charge of left nucleus (in a.u.)
 double precision::   z2    ! charge of right nucleus (in a.u.)
 double precision::   z3    ! charge of center nucleus (in a.u.)
 double precision::   R0    ! R-grid start
 double precision:: Rend    ! R-grid end
 double precision::   x0    ! x-grid start
 double precision:: xend    ! x-grid end
 double precision::   E0    ! Energy-grid start
 double precision:: Eend    ! Energy-grid end
end module pot_param

module global_arrays
use data_grid
 double precision:: PR(NR), Px(Nx)          ! nuclear & electronic momentum grid
 double precision:: pot(NR,Nx,Nx)           ! potential with electrons
 double precision:: ionicpot(NR)            ! potential without electrons
 double precision:: mu(NR,Nstates-1)        ! dipole moment
 complex*16      :: mu_ionic(NR,Nx)         ! ionic dipole moment (NR,Nx)
 complex*16      :: mu_ionic_e(NR,NE)       ! ionic dipole moment (NR,NE)
 double precision:: adb(NR,Nstates)         ! electronic Eigenenergies (potential energy surface)
 double precision:: ewf(NR,Nx,Nx,Nstates)   ! electronic wave function
 double precision:: nwf(NR,Vstates)         ! nuclear wave function
 double precision:: cofR(NR), cofx(Nx)      ! nuclear & electronic cutoff-functions
end module

module FFTW3
  use, intrinsic :: iso_c_binding
  ! include '/usr/include/fftw3.f03'                                        ! Desktop packet
  ! include '/cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/include/fftw3.f03' ! ARA cluster
  include '/home/ki96zib/apps/fftw/include/fftw3.f'
end module

module options
 logical::  option_adiabatic    ! true: do imaginary time propagation, false: use existing data
 logical::  xuv						! status of xuv pulse. false for off, true for on
end module