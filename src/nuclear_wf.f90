module MODnuclear_wavefct
  contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Calculation of nuclear wavefunction (vibrational states) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine nuclear_wavefct

    !%%%%% Subroutine calculates solution of nuclear SE (1D) in the elec. groundstate (adb(:,1)) %%

    use pot_param
    use global_arrays
    use MODsubroutines
    use FFTW3

    implicit none
    integer I,N,K,M
    integer*8 planF, planB, istep

    double precision:: psi(NR), psi_stored(NR)
    double precision:: E, E_stored
    double precision:: norm
    double precision:: R
    double precision:: vprop(NR)

    open(1,file='out/vib_states.dat'  ,status='replace')
    open(2,file='out/vib_energies.dat',status='replace')

    istep = 1.d6  ! maximum number of steps in ITP

  call dfftw_plan_r2r_1d(planF, NR, psi, psi, FFTW_R2HC, FFTW_ESTIMATE)
  call dfftw_plan_r2r_1d(planB, NR, psi, psi, FFTW_HC2R, FFTW_ESTIMATE)

  print*, '  Nuclear states'

  ! Defining potential propagator using elec. groundstate
  do I = 1, NR
    vprop(i) = exp(-0.5d0 *dt_vib *adb(i,1))
  end do


    !%%%%%%%%%% Loop over nuclear (vibrational) states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    states: do N = 1, Vstates


      ! initial wavefunction estimated as gaussian distribution
      do I = 1, NR

        R = R0 + (i-1) *dR

        psi(i) = exp(-beta_nucl * (R - Req)**2) + (-1.d0)**(N - 1) * exp(-beta_nucl * (R + Req)**2)

      end do

      ! Normalization of initial wavefunction
      call integ_real_1D(psi, psi, norm)
      psi = psi / sqrt(norm)

      E = 0.d0  ! setting initial eigenvalue to zero

      !%%%%%%%%%% Loop for Imaginary Time Propagation (ITP) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ITP: do K = 1, istep

      ! Storage of previous step
        psi_stored = psi  ! storing of psi for eigenvalue calculation
          E_stored = E  ! storing of E for checking exit criterion

        ! Extracting already calculated states
        if (N.GT.1) then
          do M = 1, (N-1)

            call integ_real_1D(nwf(:,M), psi, norm) ! projecting psi on existing states

            psi = psi - norm * nwf(:,M)           ! extracting projection

          end do
        end if

        !%%%%%%%%%% Propagation via split-step method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! 1/2 potential propagation
        psi = psi * vprop

        ! kinetic propagation in Fourier domain (= momentum space)
        call dfftw_execute(planF)
        psi = psi * exp(-dt_vib * PR**2 / (2.d0 * mass))
        call dfftw_execute(planB)
        psi = psi / dble(NR)    ! normalization due to Fourier Transform

        ! 1/2 potential propagation
        psi = psi * vprop


        !%%%%%%%%%% Checking for convergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Calculating new eigenvalue
        call eigenvalue_real_1D(psi, psi_stored, E, dt_vib)

      ! Normalization of wavefunction for next step
        call integ_real_1D(psi, psi, norm)
        psi = psi / sqrt(norm) 

        ! Comparing new eigenvalue with previous eigenvalue
        if ( abs(E - E_stored) .LE. thresh ) then
        goto 1    ! jump to closing step
        end if

      end do ITP

      print*,'Iteration not converged!'
      print*,'Program stopped!'
      print*
      print*,'E =',           E        / cm2au
      print*,'E_stored =',    E_stored / cm2au
      print*,'thresh =',    thresh   / cm2au
      print*,'step =',        K

      STOP

      !%%%%%%%%%% Storing eigenenergies & eigenfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      1 continue

      ! storing eigenfunctions
      nwf(:,N) = psi  ! storing in array

      ! storing eigenenergies
      write(2,*) N, sngl(E *au2eV)
      write(6,"(A,I1,A,F7.2,A)") '    Vib. State ', N, ': ', sngl(E *au2eV), ' eV'

    end do states


    do I = 1, NR  ! storing in vib_states.dat
      R = R0 + (i - 1) *dR
      write(1,*) sngl(R *au2A), (sngl(nwf(i,n)), N = 1, Vstates)
    end do


    call dfftw_destroy_plan(planF)
    call dfftw_destroy_plan(planB)

    close(1,status='keep')
    close(2,status='keep')

    return
  end subroutine nuclear_wavefct


end module