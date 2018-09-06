! module MODprop1D
!   contains

!   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   !%%%%%% Propagation in 1D: Nuclear coordinate - BO case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   subroutine prop1D

!     use global_arrays
!     use pot_param
!     use MODsubroutines
!     use FFTW3

!     implicit none

!   integer I,J,T
!   integer*8 planF, planB

!   double precision R, time
!   double precision norm(2), norm_ion  ! Norm of wave function & ionic wave function
!   double precision en                 ! energy
!   double precision pes(NE)            ! Photoelectron spectrum
!   double precision E_ir, E_xuv(NE)    ! electric field for pump (IR) & probe (XUV) pulses
!   double precision evV(2), evT(2)     ! expectation values for pot. & kinetic energy
!   double precision evP(2), evR(2)     ! expectation values momentum and nuclear position
!   double precision pdens(NR,2)        ! momentum density

!   complex*16 psi(NR), kprop(NR)   ! auxiallary wave function, kinetic propagator
!   complex*16 psi_ges(NR,2)        ! wave function
!   complex*16 psi_ion(NR,NE)       ! ionic wave function
!   complex*16 tout(2,2)            ! diagonal matrix


!   open(1,file='out/1d_initial_wavefunction.out',  status='replace')
!   open(2,file='out/1d_norm.out',                  status='replace')
!   open(3,file='out/1d_expectation_values.out',    status='replace')
!   open(4,file='out/1d_nuclear_density.out',       status='replace')
!   open(5,file='out/1d_momentum_density.out',      status='replace')
!   open(7,file='out/1d_photoelectron_spectrum.out',status='replace')


!   call dfftw_plan_dft_1d(planF, NR, psi, psi, FFTW_FORWARD, FFTW_MEASURE)
!   call dfftw_plan_dft_1d(planB, NR, psi, psi, FFTW_BACKWARD,FFTW_MEASURE)

!   call ionic_dipole_e

!   psi_ges = (0.d0,0.d0)    
!   psi_ion = (0.d0,0.d0)

!   ! Setting up (nuclear) kinetic propagator & initial wave packet
!   do I = 1, NR
!     R = R0 + (i-1)* dR
!     kprop(I) = exp(-im *dt * PR(I)**2/(4.d0*mass))    ! kinetic propagator
!     psi_ges(i,1) = nwf(i,1) !exp(-beta_R *(R-R_init)**2)          ! initial wave packet
!         psi_ges(i,2) = 0d0 !exp(-beta_R *(R-R_init)**2)
!   end do 

!   call integ_complex_v2(dR, psi_ges, norm)
!   do J = 1, 2
!     psi_ges(:,j) = psi_ges(:,j) / sqrt(2*norm(j))                 ! normalizing wave packet
!   end do

!   do I = 1, NR
!     R = R0 + (i-1)* dR
!     write(1,*) sngl(R *au2A), sngl(abs(psi_ges(i,1))**2), sngl(abs(psi_ges(i,2))**2)
!   end do


!   !%%%%%%%%%% Propagation Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
!   print*
!   write(6,"(A,A,$)") '  1D Propagation...', char(13)

!   timeloop: do T = 1, Nt

!       time = T * dt

!       evV = 0d0
!       evT = 0d0
!       evP = 0d0
!       evR = 0d0

!       ! IR Laser field, gaussian envelope
!       E_ir = E0_ir *exp(-beta_ir * (time - t0_ir)**2)* sin(omega_ir * (time- t0_ir) + cep*pi)

!       !%%%%%%%%%% Split-Step Methode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
!       !%%%%% 1st kinetic Step %%%%%
!       ! Propagation of bound portion of wavefunction
!       do J = 1, 2
!         do I = 1, NR                  ! 1D transformation of R-coordinate
!           psi(i) = psi_ges(i,j)       ! psi: auxillary variable
!         end do

!         call dfftw_execute(planF)
!         psi = psi * kprop             ! 1st kintetic
!         call dfftw_execute(planB)
!         psi = psi / dble(NR)

!         do I = 1, NR
!           psi_ges(i,j) = psi(i)
!         end do
!       end do

!       ! Propagation of ionized portion of wavefunction
!       do J = 1, NE
!         do I = 1, NR                  ! 1D transformation of R-coordinate
!           psi(i) = psi_ion(i,j)       ! psi: auxillary array
!         end do

!         call dfftw_execute(planF)
!         psi = psi * kprop             ! 1st kintetic propagation
!         call dfftw_execute(planB)
!         psi = psi / dble(NR)

!         do I = 1, NR
!           psi_ion(i,j) = psi(i)
!         end do
!       end do

!       !%%%%% Potential Propagation %%%%%
!       ! Field interaction - diagonalizing of interaction matrix
!       do I = 1, NR
!         call pulse(tout, mu(i,:), E_ir)   
!         psi_ges(i,:) = matmul(psi_ges(i,:),tout(:,:))
!       end do

!       ! Potential propagation of nuclear wave function
!       do J = 1, 2
!         do I = 1, NR
!           psi_ges(i,j) = psi_ges(i,j) * exp(-im * dt * adb(i,j))
!         end do
!       end do

!       ! Definition of probe pulse
!       do J = 1, NE
!         en = E0 + (j-1) *dE
!         E_xuv(j) = exp(-beta_xuv * (time - t0_xuv)**2 - im * (omega_xuv - en) * time)
!       end do

!       do J = 1, NE
!         do I = 1, NR
!           ! Estimating ionic wave function via perturbation theory
!           psi_ion(i,j) = psi_ion(i,j) + (1.d0/im *psi_ges(i,1) * mu_ionic_e(i,j) + psi_ges(i,2)) *dt *E_xuv(j)
!           ! Potential propagation of ionic wave function
!           psi_ion(i,j) = psi_ion(i,j) * exp(-im *dt *ionicpot(i))
!         end do
!       end do

!       !%%%%% 2nd kinetic Step %%%%%
!       ! Propagation of bound portion of wavefunction
!       do J = 1, 2
!         do I = 1, NR                  ! 1D transformation of R-coordinate
!           psi(i) = psi_ges(i,j)       ! psi: auxillary variable
!         end do

!         call dfftw_execute(planF)
!         psi = psi * kprop             ! 2nd kintetic propagation
!         psi = psi / sqrt(dble(NR))

!           do I = 1, NR
!             evT(j) = evT(j) + abs(psi(i))**2 *(PR(i)**2/(2.d0*mass))
!             evP(j) = evP(j) + abs(psi(i))**2 *PR(i)
!             pdens(i,j) = abs(psi(i))**2
!           end do

!             evT = evT *dpR
!             evP = evP *dpR

!         call dfftw_execute(planB)
!         psi = psi / sqrt(dble(NR))

!         do I = 1, NR
!           psi_ges(i,j) = psi(i)
!         end do
!       end do

!       ! Propagation of ionized portion of wavefunction
!       do J = 1, NE
!         do I = 1, NR                  ! 1D transformation of R-coordinate
!           psi(i) = psi_ion(i,j)       ! psi: auxillary array
!         end do

!         call dfftw_execute(planF)
!         psi = psi * kprop             ! 1st kintetic
!         call dfftw_execute(planB)
!         psi = psi / dble(NR)

!         do I = 1, NR
!           psi_ion(i,j) = psi(i)
!         end do
!       end do

!       !%%%%%%%%%% End of Propagation | Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
!       if(mod(T,100).EQ.0) then

!       ! Norm of wave function
!       call integ_complex_v2(dR, psi_ges, norm)

!       ! Norm if ionic wave function
!       call integ_complex(dE, psi_ion, norm_ion)

!       ! Calculation of expectation values
!       do I = 1, NR
!         R = R0 + (i-1) *dR

!         evR(1) = abs(psi_ges(i,1))**2 * R
!         evR(2) = abs(psi_ges(i,2))**2 * R

!         evV(1) = evV(1)+ abs(psi_ges(I,1))**2 * adb(i,1)
!         evV(2) = evV(2)+ abs(psi_ges(I,2))**2 * adb(i,2)
!       end do

!         evR = evR *dR
!         evV = evV *dR

!        do J = 1, 2
!          if (norm(j).GE.1d-8) then
!            evV(j) = evV(j) / norm(j)
!            evT(j) = evT(j) / norm(j)
!            evP(j) = evP(j) / norm(j)
!            evR(j) = evR(j) / norm(j)
!          end if
!       end do

!       ! Output of Norm
!       write(2,*) sngl(time *au2fs), sngl(norm), sngl(norm_ion)
!       ! Output of expectation values
!       write(3,*) sngl(time *au2fs), sngl(evT *au2eV), sngl(evV *au2eV), sngl(evP), sngl(evR *au2A)

!       ! Nuclear density as function of time
!       do I = 1, NR
!         R = R0 + (i-1) *dR
!         write(4,*) sngl(time *au2fs), sngl(R *au2A), sngl(abs(psi_ges(i,1)**2)), sngl(abs(psi_ges(i,2)**2))
!       end do
!         write(4,*)

!       ! Momentum density as function of time
!       do I = NR/2+1, NR
!         write(5,*) sngl(time*au2fs), sngl(PR(i)), sngl(pdens(i,1)), sngl(pdens(i,2))
!       end do
!       do I = 1, NR/2
!         write(5,*) sngl(time*au2fs), sngl(PR(i)), sngl(pdens(i,1)), sngl(pdens(i,2))
!       end do
!         write(5,*)

!       end if


!       !%%%%%%%%%% Multiplication with Cutoff-function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       !%%%%%%%%%% in order to avoid reflections at grid end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       do I = 1, NR
!         psi_ges(i,:) = psi_ges(i,:) * cofR(i)  
!             psi_ion(i,:) = psi_ion(i,:) * cofR(i)
!           end do

!   end do timeloop


!   ! Calculating spectral distribution
!   call spectrum(dE, psi_ion, pes)

!   do J = 1, NE
!     en = E0 + (j-1) *dE
!     write(7,*) sngl(en *au2eV), sngl(pes(j))
!   end do

!   write(6,*) '  1D Propagation...   done.'
!   print*

!   call dfftw_destroy_plan(planF)
!   call dfftw_destroy_plan(planB)

!   close(1,status='keep')
!   close(2,status='keep')
!   close(3,status='keep')
!   close(4,status='keep')
!   close(5,status='keep')
!   close(7,status='keep')

!   return
!   end subroutine

! end module