module MODadiabatic
  contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Calculation of electronic wavefunction, BO curves & dipole moment %%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine adiabatic_surface

    use pot_param
    use global_arrays
    use MODsubroutines
    use FFTW3
    use options
    use omp_lib

    implicit none
    integer I,J1,J2,K,M,N
    integer*8 planF, planB, istep
    double precision psi(Nx,Nx), psi_stored(Nx,Nx)
    double precision:: R, x1, x2
    double precision:: E, E_stored
    double precision:: norm
    double precision:: kprop(Nx,Nx), vprop(Nx,Nx)  
    double precision:: wtime
    

      !%%%%%%%%%% If files already exist, reading instead of recalculating %%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (option_adiabatic) then
        call adiabatic_read
        goto 2
      end if    

    ! open(1,file='out/elec_states.dat'       ,status='replace')
    ! open(2,file='out/elec_dipole_moment.dat',status='replace')
    open(3,file='out/elec_energies.dat'     ,status='replace')

    istep = 1d6 ! maximum number of steps in ITP

    call dfftw_plan_r2r_2d(planF, Nx, Nx, psi, psi, FFTW_R2HC, FFTW_R2HC, FFTW_ESTIMATE)
    call dfftw_plan_r2r_2d(planB, Nx, Nx, psi, psi, FFTW_HC2R, FFTW_HC2R, FFTW_ESTIMATE)

    print*
    print*, '  Surface scan'

    ! Defining kinetic propagator
    do J1 = 1, Nx
    do J2 = 1, Nx
      kprop(j1,j2) = exp(-dt_adia * (Px(j1)**2 + Px(j2)**2) / (2.d0 * me) )
    end do
    end do

    !%%%%%%%%%% Loop over electronic states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    states: do N = 1, Nstates

      write(6,"(A,I1,A)") '    Imaginary Time Propagation, State ', N, '                  '

      ! initial wavefunction estimated as gaussian distribution
      do J1 = 1, Nx
        x1 = x0 + (j1 - 1) *dx

      do J2 = 1, Nx
        x2 = x0 + (j2 - 1) *dx    

        ! Initial guess for 2e wave function
        psi(j1,j2) = + exp(-beta_elec *((x1 - xeq1)**2 + (x2 + xeq2)**2)) &
                  &  - exp(-beta_elec *((x1 + xeq1)**2 + (x2 - xeq2)**2))

      end do
      end do


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%% TEST 2D FT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! open(5,file='ft_test_initial.out',status='replace')

!   do j1 = 1, Nx
!   do j2 = 1, Nx

!         x1 = x0 + (j1 - 1) *dx
!         x2 = x0 + (j2 - 1) *dx 
! !     if (((j1-128)**2+(j2-128)**2).LE.2) then
! !       psi(j1,j2) = 1
! !     else
! !       psi(j1,j2) = 0
! !     end if

!     write(5,*) x1*au2A, x2*au2a, psi(j1,j2)

!   end do
!   end do

! close(5,status='keep')  

! !   call dfftw_execute(planF)

! ! open(7,file='ft_test_transformed.out',status='replace')

! !   do j1 = 1, Nx
! !   do j2 = 1, Nx

! !     write(7,*) j1, j2, psi(j1,j2)

! !   end do
! !   end do

! ! close(7,status='keep') 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%% TEST 2D FT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    ! Normalization of initial wavefunction
      call integ_real_2D(psi, psi, norm)
      psi = psi / sqrt(norm)

      !%%%%%%%%%% Loop over nuclear coordinate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      R_coordinate: do I = 1, NR

        R = R0 + (i - 1) *dR

        ! Defining potential propagator
        do J1 = 1, Nx
        do J2 = 1, Nx
          vprop(j1,j2) = exp(-0.5d0 * dt_adia * pot(i,j1,j2))
        end do
        end do

        E = 0.d0  ! setting initial eigenvalue to zero


        !%%%%%%%%%% Loop for Imaginary Time Propagation (ITP) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ITP: do K = 1, istep

          ! if (mod(K,1000) == 0 ) print*, 'K = ', K

          ! Storage of previous step
          psi_stored = psi      ! storing of psi for eigenvalue calculation
            E_stored = E        ! storing of E for checking exit criterion

          ! Extracting already calculated states
          if (N.GT.1) then
            do M = 1, (N-1)

              call integ_real_2D(ewf(i,:,:,M), psi, norm)   ! projecting psi on existing states

              psi = psi - norm * ewf(i,:,:,M)             ! extracting projection

            end do
          end if

          call integ_real_2D(psi, psi, norm)
          psi = psi / sqrt(norm)

          !%%%%%%%%%% Propagation via split-step method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ! 1/2 potential propagation
          psi = psi * vprop

          ! kinetic propagation in Fourier domain (= momentum space)
          call dfftw_execute(planF)
          psi = psi * kprop
          call dfftw_execute(planB)
          psi = psi / dble(Nx**2)   ! normalization due to Fourier Transform

          ! 1/2 potential propagation
          psi = psi * vprop

          !%%%%%%%%%% Checking for convergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ! Calculating new eigenvalue
          call eigenvalue_real_2D(psi, psi_stored, E, dt_adia)

          ! Enforcing symmetry (or anti-symmetry)
          do j1 = 1, Nx
            do j2 = j1+1, Nx
              psi(j1,j2) = sym * psi(j2,j1)
            end do
          end do

          ! Normalization of wavefunction for next step
          call integ_real_2D(psi, psi, norm)
          psi = psi / sqrt(norm) 

          ! Comparing new eigenvalue with previous eigenvalue
          if ( abs(E - E_stored) .LE. thresh ) then

            ! print 10, 'NR: ', i, ' / ', NR, '      imag. steps: ', k, char(13)
            ! 10 format (2X, A, I3, A, I3, A, I7, A, $)

            goto 1    ! jump to closing step

          end if


      end do ITP


      print*,'Iteration not converged!'
      print*,'Program stopped!'
      print*
      print*,'E =',           E / cm2au
      print*,'E_stored =',   E_stored / cm2au
      print*,'thresh =', thresh / cm2au
      print*,'step =',        K

      STOP


      !%%%%%%%%%% Storing eigenenergies & eigenfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      1 continue

      ! storing eigenenergy (= BO adiabatic potential surfaces)
      adb(i,N) = E

      ! storing eigenfunctions
      ewf(i,:,:,N) = psi    ! Note that this array format is not well suited (ewf(j,i,n) would be better), but is chosen for coherence reasons with rest of code


    end do R_coordinate


  end do states


  !%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  mu = 0.d0 ! setting dipole moment to zero

  do I = 1, NR
    R = R0 + (i - 1) * dR

    ! do J1 = 1, Nx
    ! do J2 = 1, Nx
      
    !   x1 = x0 + (j1 - 1) * dx
    !   x2 = x0 + (j2 - 1) * dx

    !   mu(i,1) = mu(i,1) + ewf(i,j1,j2,1) * 2*(R - x1) * ewf(i,j1,j2,1)    ! permanent dipole moment
    !   mu(i,2) = mu(i,2) + ewf(i,j1,j2,2) * 2*(R - x1) * ewf(i,j1,j2,2)    ! permanent dipole moment
    !   mu(i,3) = mu(i,3) + ewf(i,j1,j2,2) * 2*(R - x1) * ewf(i,j1,j2,1)    ! transition dipole moment

    !   if ( (i.eq.1).OR.(mod(i,32).eq.0) ) then
    !    write(1,*) sngl(R *au2A), sngl(x1 *au2A), sngl(x2 *au2A), ewf(i,j1,j2,:)     ! output elec. states
    !   end if
      
    ! end do
    ! end do

    ! if ( (i.eq.1).OR.(mod(i,32).eq.0) ) then
    !   write(1,*)
    !   write(1,*)
    ! end if

    ! mu(i,:) = mu(i,:) * dx**2

    ! write(2,*) sngl(R *au2A), mu(i,:)                   ! output dipole moment
    write(3,*) sngl(R *au2A), adb(i,:) *au2eV           ! output BO curves

  end do


    call dfftw_destroy_plan(planF)
    call dfftw_destroy_plan(planB)

    ! close(1,status='keep')
    ! close(2,status='keep')
    close(3,status='keep')

    2 continue

    print*, '                                                        ' ! clears terminal line

    return
  end subroutine adiabatic_surface



  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Reading initial wavefunction, BO curves & dipole moment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine adiabatic_read

    use data_au
    use global_arrays

    integer I, J1, J2
    double precision dummy

    open(1,file='out/states_elec.dat'  ,status='old')
    open(2,file='out/dipole_moment.dat',status='old')
    open(3,file='out/BO_curves.dat'    ,status='old')


    print*, 'Reading wavefunctions...'
    print*

    do I = 1, NR

      do J1 = 1, Nx
      do J2 = 1, Nx

        read(1,*) dummy, dummy, ewf(i,j1,j2,:)

      end do
      end do

      read(2,*) dummy, mu(i,:)
      read(3,*) dummy, adb(i,:)

    end do


    adb = adb / au2eV


    close(1,status='keep')
    close(2,status='keep')
    close(3,status='keep')


    return
  end subroutine adiabatic_read


end module