module MODpotential
  contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Setup of potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine make_potential

    use global_arrays
    use pot_param

    implicit none
    integer:: i, j1, j2
    double precision:: v12, v13, v23     ! repulsion of cores
    double precision:: e1e2              ! repulsion of electrons
    double precision:: v1e1, v2e1, v3e1  ! attraction of electron1 and cores
    double precision:: v1e2, v2e2, v3e2  ! attraction of electron2 and cores
    double precision:: R, x1, x2         ! coordinates of nucleus and electrons
    double precision:: erf               ! error function

    ! open(1,file='out/pot.out',status='replace')
    ! open(2,file='out/pot_ionic.out',status='replace')

    write(6,'(12X,A)'), 'Setting up potential.'

    do i = 1, NR
  
        R = R0 + (i - 1) *dR
  
        v12 = 0.d0                    ! core-core repulsion (left & right core, here: constant)
        v13 = (z1 * z3) / abs(rl - R) ! core-core repulsion (with left core)
        v23 = (z2 * z3) / abs(rr - R) ! core-core repulsion (with right core)
  
  
        do j1 = 1, Nx
        do j2 = 1, Nx

            x1 = x0 + (j1 - 1) *dx
            x2 = x0 + (j2 - 1) *dx
    
            if (abs(rl - x1).LE.1d-5) then   ! core-electron attraction (left core)
              v1e1 = -z1 * (2.d0 / sqrt(pi))
              v1e1 = v1e1 / rf
            else      
              v1e1 = -z1 * (erf(abs(rl - x1) / rf)/ abs(rl - x1))
            end if

            if (abs(rl - x2).LE.1d-5) then   ! core-electron attraction (left core)
              v1e2 = -z1 * (2.d0 / sqrt(pi))
              v1e2 = v1e2 / rf
            else      
              v1e2 = -z1 * (erf(abs(rl - x2) / rf)/ abs(rl - x2))
            end if

    
            if (abs(rr - x1).LE.1d-5) then   ! core-electron attraction (right core)
              v2e1 = -z2 * (2.d0 / sqrt(pi))
              v2e1 = v2e1 / rf
            else
              v2e1 = -z2 * (erf(abs(rr - x1) / rf)/ abs(rr - x1))
            end if
  
            if (abs(rr - x2).LE.1d-5) then   ! core-electron attraction (right core)
              v2e2 = -z2 * (2.d0 / sqrt(pi))
              v2e2 = v2e2 / rf
            else
              v2e2 = -z2 * (erf(abs(rr - x2) / rf)/ abs(rr - x2))
            end if


            if (abs(R - x1).LE.1d-5) then    ! core-electron attraction (moving core)
              v3e1 = -z3 * (2.d0 / sqrt(pi))
              v3e1 = v3e1 / rc
            else
              v3e1 = -z3 * (erf(abs(R - x1) / rc ) / abs(R - x1))
            end if

            if (abs(R - x2).LE.1d-5) then    ! core-electron attraction (moving core)
              v3e2 = -z3 * (2.d0 / sqrt(pi))
              v3e2 = v3e2 / rc
            else
              v3e2 = -z3 * (erf(abs(R - x2) / rc ) / abs(R - x2))
            end if


            if (abs(x1 - x2).LE.1d-5) then    ! electron-electron repulsion
              e1e2 = 2.d0 / sqrt(pi)
              e1e2 = e1e2 / re
            else
              e1e2 = erf(abs(x1 - x2) / re ) / abs(x1 - x2)
            end if

            ! total potential
            pot(i,j1,j2) = v12 + v13 + v23 + v1e1 + v2e1 + v3e1 + v1e2 + v2e2 + v3e2 + e1e2

        end do

            pot1e(i,j2) = v12 + v13 + v23 + v1e1 + v2e1 + v3e1
        end do

        ionicpot(i) = v12 + v13 + v23

    end do
    

    ! do I = 1, NR
    !   R = R0 + (i - 1) *dR

    !   do j1 = 1, Nx
    !   do j2 = 1, Nx

    !     x1 = x0 + (j1 - 1) *dx
    !     x2 = x0 + (j2 - 1) *dx

    !     if ( (i.eq.1).OR.(mod(i,32).eq.0) ) then
    !       write(1,*) sngl(R *au2A), sngl(x1 *au2A), sngl(x2 *au2A), sngl(pot(i,j1,j2)*au2eV)
    !     end if

    !   end do
    !   end do

    !   if ( (i.eq.1).OR.(mod(i,32).eq.0) ) then
    !     write(1,*)
    !     write(1,*)
    !   end if

    !   write(2,*) sngl(R *au2a), sngl(ionicpot(i) *au2eV)

    ! end do


    ! close(1,status='keep')
    ! close(2,status='keep')


    return
  end subroutine

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Error function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function erf(X)
    double precision x, erf

      erf = 1.d0 - erfcc(x)

    return
  end function


  function erfcc(X)
    double precision x, z, t
    double precision erfcc
 
      z = abs(x)

      t = 1.d0 / (1.d0 + 0.5d0 *z)

      erfcc=t*dexp(-Z*Z-1.26551223+t*(1.00002368+t*(.37409196+&
            &    t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+&
            &    t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
  
      if (x.LT.0.d0) erfcc = 2.d0 - erfcc

    return
  end function


end module