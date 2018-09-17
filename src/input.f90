module MODinput
  contains

  subroutine input

   use pot_param
   use data_grid
   use options

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%% Reading input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  open(1,file='input',status='old')

      read(1,*) dt      ! length of time step during propagation (in fs)
      read(1,*) dt_adia ! length of time step during calc. for elec. wf   (in fs)
      read(1,*) dt_vib  ! length of time step during calc. for nuclear wf (in fs)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*) rc      ! screening parameter: electron-(moving)-nucleus (in A)
      read(1,*) rf      ! screening parameter: electron-(fixed)-nuclei   (in A)
      read(1,*) re      ! screening parameter: electron-electron         (in A)
      read(1,*) rl      ! Position of left nucleus  (in A)
      read(1,*) rr      ! Position of right nucleus (in A)
      read(1,*) z1      ! charge of left nucleus   (in a.u.)
      read(1,*) z2      ! charge of right nucleus  (in a.u.)
      read(1,*) z3      ! charge of center nucleus (in a.u.)
      read(1,*) R0      ! R-grid start (in A)
      read(1,*) Rend    ! R-grid end   (in A)
      read(1,*) x0      ! x-grid start (in A)
      read(1,*) xend    ! x-grid end   (in A)
      read(1,*) E0      ! Energy-grid start (in eV, used for ionic wf with pert. theory)
      read(1,*) Eend    ! Energy-grid end (in eV)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*) init_state ! initial state for wave packet (1 = groundstate)
      read(1,*) R_init     ! center of initial wave packet used for propagation (in A)
      read(1,*) beta_R     ! width of initial wave packet (exp(-beta_R *(R-R_init)^2) ) (in A^-2)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*) lambda_ir   ! wavelength of IR pulse (in nm)
      read(1,*) lambda_xuv  ! wavelength of XUV pulse (in nm)
      read(1,*) E0_ir       ! field strength of IR pulse (in a.u.)
      read(1,*) E0_xuv      ! field strength of XUV pulse (in a.u.)
      read(1,*) xuv         ! status of xuv pulse. false for off, true for on
      read(1,*) fwhm_ir     ! fwhm of IR pulse (in fs)
      read(1,*) fwhm_xuv    ! fwhm of XUV pulse (in fs)
      read(1,*) t0_ir       ! center of IR pulse (in fs)
      read(1,*) t0_xuv      ! center of XUV pulse (in fs)
      read(1,*) cep         ! Carrier-Envelope-Phase of IR pulse (in units of PI)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*) Req       ! center of initial wavefunction for nuclear imaginary time propagation
      read(1,*) xeq1      ! center of initial wavefunction (x1) for electronic imaginary time propagation
      read(1,*) xeq2      ! center of initial wavefunction (x2) for electronic imaginary time propagation
      read(1,*) beta_elec ! width of initial wavefunction for electronic ITP (exp(-beta...))
      read(1,*) beta_nucl ! width of initial wavefunction for nuclear ITP (exp(-beta...))
      read(1,*) thresh    ! threshold value for ITPs
      read(1,*) sym       ! spatial symmetry of eigenfunctions (1 even, -1 uneven)
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*) option_adiabatic  ! true: do imaginary time propagation, false: use existing data

  close(1,status='keep')

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%% Conversion to atomic units %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt      = dt      / au2fs
    dt_adia = dt_adia / au2fs
    dt_vib  = dt_vib  / au2fs

    rc   = rc   / au2A
    rf   = rf   / au2A
    re   = re   / au2A
    rl   = rl   / au2A
    rr   = rr   / au2A

    R0   = R0   / au2A
    Rend = Rend / au2A
    x0   = x0   / au2A
    xend = xend / au2A

    E0   = E0   / au2eV
    Eend = Eend / au2eV
    R_init = R_init / au2A
    beta_R  = beta_R  * au2A**2  ! conversion from A^-2 to a.u.^-2

    t0_ir  = t0_ir  / au2fs
    t0_xuv = t0_xuv / au2fs
    omega_ir  = (1d0 / (lambda_ir  *1d-7)) *cm2au ! first conversion to wavenumbers then to a.u.
    omega_xuv = (1d0 / (lambda_xuv *1d-7)) *cm2au ! first conversion to wavenumbers then to a.u.

    Req  = Req   / au2A
    xeq  = xeq   / au2A
    
    beta_elec = beta_elec  !* au2A**2   ! conversion from A^-2 to a.u.^-2
    beta_nucl = beta_nucl

    beta_ir  = 1d0 / ( 2*(fwhm_ir  / (2d0 * sqrt(2d0 * log(2d0))) /au2fs)**2 )   ! e^(-beta*(...)^2)
    beta_xuv = 1d0 / ( 2*(fwhm_xuv / (2d0 * sqrt(2d0 * log(2d0))) /au2fs)**2 )

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%% Defining grid parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dR = (Rend - R0) / (NR - 1)   ! spacing of R grid
    dx = (xend - x0) / (Nx - 1)   ! spacing of x grid
    dE = (Eend - E0) / (NE - 1)   ! spacing of energy grid

    dpR = (2d0 * pi) / (NR * dR)  ! spacing of R momentum grid
    dpx = (2d0 * pi) / (Nx * dx)  ! spacing of x momentum grid

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%% Print to Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    print*
    print*, '================================================================================'
    print*, '%%%%%%%%%%%%%%%%%%%%%%%%% Shin-Metiu 2e Model System %%%%%%%%%%%%%%%%%%%%%%%%%%%'  
    print*, '================================================================================'
    print*
    print*, '   Parameters'
    print*, '______________________________________'
    print*
    write(6,"(A,F4.1,A,F4.1,A,F5.1,A)") '  Re = ', sngl(re *au2A), ' A    Rc = ', sngl(rc *au2A), &
    &                                   ' A    R_init = ', sngl(R_init*au2A), ' A'
    print*
    write(6,"(A,F5.1,A,F5.3,A)") '   T = ', sngl(Nt*dt *au2fs), ' fs     @ dt = ', sngl(dt *au2fs), ' fs'
    print*
    write(6,"(A,F4.1,A)")          ' cep = ', sngl(cep), ' pi'
    print*, '______________________________________'
    print*
    print*, 'Progress:'

    return
  end subroutine

end module