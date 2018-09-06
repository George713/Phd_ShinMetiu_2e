module MODprop3d
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%% Propagation in 2D: Nuclear & electronic coordinate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine prop3D

  use global_arrays
  use pot_param
  use MODsubroutines
  use FFTW3
  use omp_lib

  implicit none

  integer I,J1,J2,V,T

  integer*8 planF, planB, planFR, planBR, planFx, planBx

  double precision x1, x2, R, E_grid(2)
  double precision norm, norm1e, norm2e
  double precision time
  double precision A, E_tot    ! vector potential & total electric field (= IR + XUV)
  double precision A_ir, A_xuv ! Electric fields of IR & XUV pulses

  double precision erR, erx     ! expectation value of R & x in real space
  double precision epR, epx     ! expectation value of R & x in momentum space
  double precision epot         ! expectation value of potential energy
  double precision pes_1e       ! pes for one energy
  double precision asy_opt_1e, asy_opt_2e ! condensed asymmetry (fwd/bwd asymmetry summed up)
  double precision sum1_1e, sum2_1e   ! sum of pes for sector I & III
  double precision sum1_2e, sum2_2e   ! sum of pes for sector I & III

  double precision, allocatable, dimension(:)  :: densR   ! nuclear density
  double precision, allocatable, dimension(:)  :: densx   ! total electronic density
  double precision, allocatable, dimension(:,:):: pes1e,pes2e     ! photoelectron spectrum
  double precision, allocatable, dimension(:,:):: asy2d_1e, asy2d_2e  ! photoelectron spectrum

  complex*16, allocatable, dimension(:,:,:):: psi                       ! inside wave function 
  complex*16, allocatable, dimension(:)    :: psiR                      ! auxillary wf for R FT
  complex*16, allocatable, dimension(:,:)  :: psix                      ! auxillary wf x FT
  complex*16, allocatable, dimension(:,:,:):: psi_out1e, psi_out1e_new  ! outside wf & outrunning wf for single ionization
  complex*16, allocatable, dimension(:,:,:):: psi_out2e, psi_out2e_new  ! outside wf & outrunning wf for double ionization
  complex*16, allocatable, dimension(:,:,:):: psix_out                  ! auxillary wf x FT of outrunning part
  complex*16, allocatable, dimension(:,:,:):: kprop                     ! kinetic propagator

  allocate( densR(NR) , densx(Nx) , pes1e(Nx,Nx), pes2e(Nx,Nx), asy2d_1e(Nx/2,Nx/2) , asy2d_2e(Nx/2,Nx/2) )
  allocate( psi(NR,Nx,Nx) , psiR(NR) , psix(Nx,Nx) , kprop(NR,Nx,Nx) )
  allocate( psi_out1e(NR,Nx,Nx) , psi_out1e_new(NR,Nx,Nx) , psi_out2e(NR,Nx,Nx) , psi_out2e_new(NR,Nx,Nx) , psix_out(NR,Nx,Nx) )

  open( 1,file='out/3d_initial_wavefunction.out', status='replace')
  open( 2,file='out/vecP_and_Efields.out',        status='replace')
  open( 3,file='out/3d_norm.out',                 status='replace')
  open( 4,file='out/3d_erR.out',                  status='replace')
  open( 5,file='out/3d_erx.out',                  status='replace')
  open( 7,file='out/3d_epR.out',                  status='replace')
  open( 8,file='out/3d_epx.out',                  status='replace')
  ! open( 9,file='out/3d_epot.out',                 status='replace')
  open(10,file='out/3d_densR.out',                status='replace')
  open(11,file='out/3d_densx.out',                status='replace')
  open(12,file='out/3d_pop.out',                  status='replace')
  open(13,file='out/3d_pop_in_superposition.out', status='replace')
  open(14,file='out/3d_nwf_in_elec_states.out',   status='replace')
  open(15,file='out/3d_nwf_in_sup_states.out',    status='replace')
  open(16,file='out/3d_nuclear_correlation.out',  status='replace')
  ! open(17,file='out/3d_pes_energy_resolved.out',  status='replace')
  open(18,file='out/3d_pes_asymmetry.out',        status='replace')
  open(19,file='out/3d_pes_momentum_resolved.out',status='replace')
  open(20,file='out/3d_pes_asymmetry_sum.out'    ,status='replace')

  call dfftw_plan_dft_3d(planF,  NR, Nx, Nx, psi, psi, FFTW_FORWARD,  FFTW_MEASURE)
  call dfftw_plan_dft_3d(planB,  NR, Nx, Nx, psi, psi, FFTW_BACKWARD, FFTW_MEASURE)
  call dfftw_plan_dft_1d(planFR, NR, psiR, psiR,       FFTW_FORWARD,  FFTW_MEASURE)
  call dfftw_plan_dft_1d(planBR, NR, psiR, psiR,       FFTW_BACKWARD, FFTW_MEASURE)
  call dfftw_plan_dft_2d(planFx, Nx, Nx, psix, psix,   FFTW_FORWARD,  FFTW_MEASURE)
  call dfftw_plan_dft_2d(planBx, Nx, Nx, psix, psix,   FFTW_BACKWARD, FFTW_MEASURE)

  write(6,"(A)") '   3D Propagation...'

  ! Setting up (nuclear-electronic) kinetic propagator
  do J1 = 1, Nx
  do J2 = 1, Nx
  do  I = 1, NR
      kprop(i,j1,j2) = exp(-im *dt *((PR(i)**2/(2.d0 *mass))+(Px(j1)**2/(2.d0*me))+(Px(j2)**2/(2.d0*me))))          
  end do
  end do
  end do

  ! Setting initial wave function
  do J1 = 1, Nx
  do J2 = 1, Nx
    x1 = x0 + (j1 - 1) * dx
    x2 = x0 + (j2 - 1) * dx

    do I = 1, NR
      R = R0 + (i - 1) * dR
      psi(i,j1,j2) = ewf(i,j1,j2,init_state) * exp(-beta_R *(R-R_init)**2)
    end do

  end do
  end do

  call integ_complex(psi, norm)
  psi = psi / sqrt(norm)        ! Normalizing intial wave function

  ! Preparing for timeloop
  A       = 0d0
  psi_out1e = ( 0d0 , 0d0 )
  psi_out2e = ( 0d0 , 0d0 )

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%% Propagation Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  timeloop: do T = 1, Nt

  time = T *dt

  ! Setting current vector potentials
  if (E0_ir.NE.0 .OR. E0_xuv.NE.0) then

    A_ir  = E0_ir /omega_ir  * exp(-beta_ir  *(time - t0_ir )**2) * sin(omega_ir  *(time - t0_ir ))  
    A_xuv = E0_xuv/omega_xuv * exp(-beta_xuv *(time - t0_xuv)**2) * sin(omega_xuv *(time - t0_xuv) + cep*pi)

    ! Calculating electric field (for storage only)
    E_tot = - (A_ir + A_xuv - A)/dt

    ! Setting vector potential
    A = A_ir + A_xuv

  end if


  !%%%%%%%%%% Split-Step Methode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  !%%%%% 1st potential Step %%%%%
  do J1 = 1, Nx
  do J2 = 1, Nx
  do  I = 1, NR
    psi(i,j1,j2) = psi(i,j1,j2) * exp(-0.5d0 * im * dt * pot(i,j1,j2) )
  end do
  end do
  end do

  !%%%%%%%% Kinetic Step %%%%%%%%

  call dfftw_execute_dft(planF, psi, psi)   ! FT to momentum space

  ! Propagation due to momentum
  psi = psi * kprop

  ! Propagation due to external fields
  if (A.NE.0) then
    do J1 = 1, Nx
    do J2 = 1, Nx
    do  I = 1, NR
      psi(i,j1,j2) = psi(i,j1,j2) * exp(-im * dt * A * ( Px(j1) + Px(j2) - PR(i)/mass ) )
    end do
    end do
    end do
  end if

  call dfftw_execute_dft(planB, psi, psi)   ! FT to real space

  psi = psi / (NR * Nx**2)                  ! re-normalization due to FT

  !%%%%% 2nd potential Step %%%%%
  do J1 = 1, Nx
  do J2 = 1, Nx
  do  I = 1, NR
      psi(i,j1,j2) = psi(i,j1,j2) * exp(-0.5d0 * im * dt * pot(i,j1,j2) )
  end do
  end do
  end do        


  !%%%%%%%%%%%%%%%% Propagation of (existing) Outrun Part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! 1st half-step of (ionic) potential propagation
  do I = 1, NR
    do J1 = 1, Nx
      psi_out1e(i,j1,:) = psi_out1e(i,j1,:) * exp(-0.5d0 * im * dt * pot1e(i,j1))
    end do
    do J1 = 1, Nx
      psi_out1e(i,:,j1) = psi_out1e(i,:,j1) * exp(-0.5d0 * im * dt * pot1e(i,j1))
    end do
    psi_out2e(i,:,:) = psi_out2e(i,:,:) * exp(-0.5d0 * im * dt * ionicpot(i))
  end do

  ! Forward FT of nuclear coordinate from psi_out(R,p_x) to psi_out(p_R,p_x)
  do J1 = 1, Nx
  do J2 = 1, Nx

    do  I = 1, NR
      psiR(i) = psi_out1e(i,j1,j2)
    end do
    
    call dfftw_execute_dft(planFR, psiR, psiR)
    
    do I = 1, NR
      psi_out1e(i,j1,j2) = psiR(i) / sqrt(dble(NR))
    end do


    do  I = 1, NR
      psiR(i) = psi_out2e(i,j1,j2)
    end do
    
    call dfftw_execute_dft(planFR, psiR, psiR)
    
    do I = 1, NR
      psi_out2e(i,j1,j2) = psiR(i) / sqrt(dble(NR))
    end do

  end do
  end do

  ! kinetic & field propagation
  if (A.NE.0) then
    do J1 = 1, Nx
    do J2 = 1, Nx
    do  I = 1, NR
      psi_out1e(i,j1,j2) = psi_out1e(i,j1,j2) * exp(-im * dt * A * (Px(j1) + Px(j2) - PR(i)/mass) )   
      psi_out2e(i,j1,j2) = psi_out2e(i,j1,j2) * exp(-im * dt * A * (Px(j1) + Px(j2) - PR(i)/mass) )    
    end do
    end do
    end do
  end if

  psi_out1e = psi_out1e * kprop
  psi_out2e = psi_out2e * kprop

  ! Backward FT of nuclear coordinate from psi_out(p_R,p_x) to psi_out(R,p_x)
  do J1 = 1, Nx
  do J2 = 1, Nx

    do I = 1, NR
      psiR(i) = psi_out1e(i,j1,j2)
    end do
    
    call dfftw_execute_dft(planBR, psiR, psiR)
    
    do I = 1, NR
      psi_out1e(i,j1,j2) = psiR(i) / sqrt(dble(NR))
    end do


    do I = 1, NR
      psiR(i) = psi_out2e(i,j1,j2)
    end do
    
    call dfftw_execute_dft(planBR, psiR, psiR)
    
    do I = 1, NR
      psi_out2e(i,j1,j2) = psiR(i) / sqrt(dble(NR))
    end do

  end do
  end do

  ! 2nd half-step of (ionic) potential propagation
  do I = 1, NR
    do J1 = 1, Nx
      psi_out1e(i,j1,:) = psi_out1e(i,j1,:) * exp(-0.5d0 * im * dt * pot1e(i,j1))
    end do
    do J1 = 1, Nx
      psi_out1e(i,:,j1) = psi_out1e(i,:,j1) * exp(-0.5d0 * im * dt * pot1e(i,j1))
    end do
    psi_out2e(i,:,:) = psi_out2e(i,:,:) * exp(-0.5d0 * im * dt * ionicpot(i))
  end do


  !%%%%%%%% Apply Cut-Off Function // Collect & add (new) Outrunning Part %%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  ! New single ionization
  do J1 = 1, Nx
  do J2 = 1, Nx
    psi_out1e_new(:,j1,j2) = psi(:,j1,j2) *( (1d0 - cofx(j1)) *cofx(j2) + cofx(j1) *(1d0 - cofx(j2)) )
  end do
  end do

  ! New double ionization
  do J1 = 1, Nx
  do J2 = 1, Nx
    psi_out2e_new(:,j1,j2) = psi(:,j1,j2) *(1d0 - cofx(j1)) *(1d0 - cofx(j2))
  end do
  end do

  ! New remaining inner part 
  do J1 = 1, Nx
  do J2 = 1, Nx
    psi(:,j1,j2) = psi(:,j1,j2) *cofx(j1) *cofx(j2)
  end do
  end do

  ! Forward FT of elec. coordinate from psi_out_new(R,x) to psi_out_new(R,p_x)
  do  I = 1, NR

    do J1 = 1, Nx
    do J2 = 1, Nx
      psix(j1,j2) = psi_out1e_new(i,j1,j2)
    end do
    end do

    call dfftw_execute_dft(planFx, psix, psix)

    do J1 = 1, Nx
    do J2 = 1, Nx
      psi_out1e_new(i,j1,j2) = psix(j1,j2) / dble(Nx)
    end do
    end do


    do J1 = 1, Nx
    do J2 = 1, Nx
      psix(j1,j2) = psi_out2e_new(i,j1,j2)
    end do
    end do

    call dfftw_execute_dft(planFx, psix, psix)

    do J1 = 1, Nx
    do J2 = 1, Nx
      psi_out2e_new(i,j1,j2) = psix(j1,j2) / dble(Nx)
    end do
    end do

  end do

  ! Adding of new outrun part to propagated existing part
  psi_out1e = psi_out1e + psi_out1e_new
  psi_out2e = psi_out2e + psi_out2e_new


  !======================================================================================
  !%%%%%%%%%%%%%%%%%% Analysis of dynamics without ionizing XUV pulse %%%%%%%%%%%%%%%%%%%
  !======================================================================================
  if (mod(T,2).NE.0) go to 10

  call pop_analysis(psi, time, T)                  ! calculates population of elec. & superposition states & nuclear correlation
  call density(psi, densR, densx)                  ! calculates nuclear & electronic density
  call integ_complex(psi, norm)
  call integ_complex_p(psi_out1e, norm1e)
  call integ_complex_p(psi_out2e, norm2e)

  ! Resetting expected values
  erR  = 0d0
  erx  = 0d0
  epR  = 0d0
  epx  = 0d0
  epot = 0d0

  ! R & x expectation value
  do I = 1, NR
    R = R0 + (i-1) *dR
    erR = erR + dble(R *densR(i))
  end do
  erR = erR *dR /norm

  do J1 = 1, Nx
    x1 = x0 + (j1 - 1) *dx
    erx  = erx  + dble(x1 * densx(j1) )
  end do
  erx  = erx  *dx /norm

  ! p_R & p_x expectation value
  call dfftw_execute_dft(planF, psi, psi)
  psi = psi / sqrt(dble(NR * Nx**2))
      
  do J2 = 1, Nx
  do J1 = 1, Nx
  do  I = 1, NR
    epR = epR + dble(PR(i) * abs(psi(i,j1,j2))**2)
    epx = epx + dble(2 * Px(j1) * abs(psi(i,j1,j2))**2)
  end do
  end do
  end do
  epR = epR *dpR /norm
  epx  = epx  *dpx**2 /norm
      
  call dfftw_execute_dft(planB, psi, psi)
  psi = psi / sqrt(dble(NR * Nx**2))

!  ! epot | expectation value of potential energy
!  do J2 = 1, Nx
!  do J1 = 1, Nx
!  do  I = 1, NR
!    epot = epot + pot(i,j1,j2) * abs(psi(i,j1,j2))**2
!  end do
!  end do
!  end do
!  epot = epot *dx**2 *dR /norm
!
!  ! Norm of outrunning part
!  ! Backward FT of elec. coordinate from psi_out(R,p_x) to psi_out(R,x)
!  do i = 1, NR
!
!    do j2 = 1, Nx
!    do j1 = 1, Nx
!      psix(j1,j2) = psi_out(i,j1,j2)
!    end do
!    end do
!
!    call dfftw_execute_dft(planBx, psix, psix)
!
!    do j2 = 1, Nx
!    do j1 = 1, Nx
!      psix_out(i,j1,j2) = psix(j1,j2) / dble(Nx)
!    end do
!    end do
!
!  end do

  !======================================================================================
  !%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  write(3,*) sngl(time *au2fs), sngl(norm) , sngl(norm1e), sngl(norm2e)
  write(4,*) sngl(time *au2fs), sngl(erR *au2A)
  write(5,*) sngl(time *au2fs), sngl(erx *au2A)
  write(7,*) sngl(time *au2fs), sngl(epR) 
  write(8,*) sngl(time *au2fs), sngl(epx)
!  write(9,*) sngl(time *au2fs), sngl(epot *au2eV)

  do I = 1, NR
    R = R0 + (i - 1) *dR
    write(10,*) sngl(time *au2fs), sngl(R *au2a), sngl(densR(i))
  end do
  write(10,*)
  write(10,*)

  do J1 = 1, Nx
    x1 = x0 + (j1 - 1) *dx
    write(11,*) sngl(time *au2fs), sngl(x1 *au2a), sngl(densx(j1))
  end do
  write(11,*)
  write(11,*)


  10 continue

  write(2,*) sngl(time *au2fs), sngl(A), sngl(E_tot)

  end do timeloop

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Photeelectron Spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (E0_xuv.NE.0) then
  
    call integ_complex_p(psi_out1e, norm1e)
    psi_out1e = psi_out1e / sqrt(norm1e)
  
    call spectrum(psi_out1e, pes1e)  ! psi_out1e(R,p_x1,p_x2)

    call integ_complex_p(psi_out2e, norm2e)
    psi_out2e = psi_out2e / sqrt(norm2e)
  
    call spectrum(psi_out2e, pes2e)  ! psi_out2e(R,p_x1,p_x2)
   
    ! do j1 = 1, Nx/2
    ! do j2 = 1, Nx/2
    !   E_grid(1) = px(j1)**2 / 2.d0       ! setting energy value E = p^2/(2m)
    !   E_grid(2) = px(j2)**2 / 2.d0
  
    !   if ( j1.eq.1 .AND. j2.eq.1 ) then
    !     pes_1e = pes(j1,j2)
    !   else if ( j1.eq.1 ) then
    !     pes_1e = pes(j1,j2) + pes(j1,Nx+2-j2)
    !   else if ( j2.eq.1 ) then
    !     pes_1e = pes(j1,j2) + pes(Nx+2-j1,j2)
    !   else
    !     pes_1e = pes(j1,j2) + pes(Nx+2-j1,j2) + pes(j1,Nx+2-j2) + pes(Nx+2-j1,Nx+2-j2)
    !   end if
  
    !   write(17,*) sngl(E_grid *au2eV), sngl(pes_1e)
  
    ! end do
    ! end do
    ! write(17,*)
    ! write(17,*)
 
  !%%%%%%% Writing momentum resolved spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    do j1 = Nx/2+1, Nx
    do j2 = Nx/2+1, Nx
 
      write(19,*) sngl(px(j1)), sngl(px(j2)), sngl(pes1e(j1,j2)), sngl(pes2e(j1,j2))
 
    end do
    do j2 = 1, Nx/2
    
      write(19,*) sngl(px(j1)), sngl(px(j2)), sngl(pes1e(j1,j2)), sngl(pes2e(j1,j2))
    
    end do
    write(19,*)
    write(19,*)
    end do
 
    do j1 = 1, Nx/2
    do j2 = Nx/2+1, Nx
 
      write(19,*) sngl(px(j1)), sngl(px(j2)), sngl(pes1e(j1,j2)), sngl(pes2e(j1,j2))
 
    end do
    do j2 = 1, Nx/2
    
      write(19,*) sngl(px(j1)), sngl(px(j2)), sngl(pes1e(j1,j2)), sngl(pes2e(j1,j2))
    
    end do
    write(19,*)
    write(19,*)
    end do
 
  !%%%%%%% Writing asymmetry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2D asymmetry (sector I minus sector II)
    asy2d_1e = 0; asy2d_2e = 0
    do j1 = 2, Nx/2
    do j2 = 2, Nx/2
      asy2d_1e(j1,j2) = (pes1e(j1,j2) - pes1e(Nx+2-j1,Nx+2-j2)) / (pes1e(j1,j2) + pes1e(Nx+2-j1,Nx+2-j2) + 1d-10)    ! 2D: sector I minus sector III (solely pos. momenta - solely neg. momenta)
      asy2d_2e(j1,j2) = (pes2e(j1,j2) - pes2e(Nx+2-j1,Nx+2-j2)) / (pes2e(j1,j2) + pes2e(Nx+2-j1,Nx+2-j2) + 1d-10)    ! 2D: sector I minus sector III (solely pos. momenta - solely neg. momenta)
    end do
    end do

    do j1 = 1, Nx/2
    do j2 = 1, Nx/2
      write(18,*) sngl(px(j1)), sngl(px(j2)), sngl(asy2d_1e(j1,j2)), sngl(asy2d_2e(j1,j2))
    end do
    write(18,*)
    write(18,*)
    end do

    ! Condensed asymmetry
    sum1_1e = 0; sum2_1e = 0
    sum1_2e = 0; sum2_2e = 0

    do j1 = 2, Nx/2
    do j2 = 2, Nx/2
      ! Summing up sectors I & III
      sum1_1e = sum1_1e + pes1e(j1,j2)
      sum2_1e = sum2_1e + pes1e(Nx+2-j1,Nx+2-j2)
      sum1_2e = sum1_2e + pes2e(j1,j2)
      sum2_2e = sum2_2e + pes2e(Nx+2-j1,Nx+2-j2)
    end do
    end do
    asy_opt_1e = (sum1_1e-sum2_1e)/(sum1_1e+sum2_1e)
    asy_opt_2e = (sum1_2e-sum2_2e)/(sum1_2e+sum2_2e)

    write(20,*) "I/III sum-wise    ", sngl(asy_opt_1e), sngl(asy_opt_2e)

  end if 

  write(6,*) '  3D Propagation...   done.'

  deallocate( densR , densx , pes1e , pes2e , asy2d_1e , asy2d_2e , psi , psiR , psix , kprop )
  deallocate( psi_out1e , psi_out1e_new , psi_out2e , psi_out2e_new , psix_out )


  close( 1,status='keep')
  close( 2,status='keep')
  close( 3,status='keep')
  close( 4,status='keep')
  close( 5,status='keep')
  close( 7,status='keep')
  close( 8,status='keep')
  ! close( 9,status='keep')
  close(10,status='keep')
  close(11,status='keep')
  close(12,status='keep')
  close(13,status='keep')
  close(14,status='keep')
  close(15,status='keep')
  close(16,status='keep')
  ! close(17,status='keep')
  close(18,status='keep')
  close(19,status='keep')
  close(20,status='keep')


  call dfftw_destroy_plan(planF)
  call dfftw_destroy_plan(planB)
  call dfftw_destroy_plan(planFR)
  call dfftw_destroy_plan(planBR)
  call dfftw_destroy_plan(planFx)
  call dfftw_destroy_plan(planBx)

  return
end subroutine prop3D


end module