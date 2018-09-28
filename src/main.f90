program metiu_model_system

  use MODinput
  use MODsubroutines
  use MODpotential
  use MODadiabatic
  use MODnuclear_wavefct
  use MODprop3D
  use omp_lib

  double precision wtime  ! walltime
  open(100,file='time.out',status='replace')
  wtime = omp_get_wtime()

!======================================================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!======================================================================================
  call input              ! reading input-file
  call p_grid             ! generating momentum grids
  call make_potential     ! calculating potential & ionic potential
  ! call make_cofs          ! setup of cofR & cofx (cutoff functions)
  call make_new_cofs          ! setup of cofR & cofx (cutoff functions)
  ! open(101,file='time1.out',status='replace')
  ! write(101,*) sngl((omp_get_wtime ( ) - wtime)/3600)       ! Calculation time in hours
  ! close(101,status='keep')

  call adiabatic_surface  ! calculating elec. wavefunction, BO curves & dipole moment
  open(102,file='time_adia.out',status='replace')
  write(102,*) sngl((omp_get_wtime ( ) - wtime)/3600)       ! Calculation time in hours
  close(102,status='keep')

  ! call nuclear_wavefct    ! calculating nuclear wavefunction
  ! open(103,file='time3.out',status='replace')
  ! write(103,*) sngl((omp_get_wtime ( ) - wtime)/3600)       ! Calculation time in hours
  ! close(103,status='keep')
  
  call prop3D             ! conducting propagation with x & R coordinate
!======================================================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!======================================================================================

  wtime = omp_get_wtime ( ) - wtime
  write(100,*) sngl(wtime/3600)       ! Calculation time in hours
  close(100,status='keep')

  print*
  print*, '================================================================================'
  print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'  
  print*, '================================================================================'
  print*

 stop
end program