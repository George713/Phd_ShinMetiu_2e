module MODsubroutines
  contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Setup of cutoff functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! subroutine make_cofs

  !   use pot_param
  !   use global_arrays

  !   implicit none
  !   integer I
  !   double precision:: R, x

  !   open(1,file='out/cofR.out',status='replace')
  !   open(2,file='out/cofx.out',status='replace')

  !   write(6,'(12X,A)'), 'Defining cutoff-functions...'

  !   !%%%%%% Make R cutoff function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !   do I = 1, NR
  !     R = R0 + (i-1) * dR
  !       if ( R .GE. 0.d0 ) then
  !         cofR(i) = 1.d0 / (1.d0 + exp( (R - (Rend *3/4) ) * 7 ))   ! cofR = .5 @ R = Rend*3/4, '7' specifies steepness
  !       end if
  !   end do
  
  !   do I = 1, NR/2
  !     cofR(i) = cofR(NR-i+1)  ! making nuclear cutoff symmetric
  !   end do

  !   cofR = 1.

  !   do I = 1, NR
  !     R = R0 + (i-1) * dR
  !     write(1,*) sngl(R *au2A), sngl(cofR(i))
  !   end do


  !   !%%%%%% Make x cutoff function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !   do I = 1, Nx
  !     x = x0 + (i-1) * dx

  !     if ( x .GE. 0.d0 ) then
  !       cofx(i) = 1.d0 / (1.d0 + exp( (x - (xend *3/4)) *Rend/xend * 7 ))
  !     end if
  !   end do

  !   do I = 1, Nx/2
  !     cofx(i) = cofx(Nx-i+1)  ! making electronic cutoff symmetric
  !   end do

  !   do I = 1, Nx
  !     x = x0 + (i-1) * dx
  !     write(2,*) sngl(x *au2A), sngl(cofx(i))
  !   end do    


  !   close(1,status='keep')
  !   close(2,status='keep')


  !   return
  ! end subroutine make_cofs



  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Setup of new cutoff masks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



  subroutine make_new_cofs

    use pot_param
    use global_arrays

    implicit none
    integer j1,j2
    double precision:: x, y
    double precision:: boundary_1, boundary_2

    open(1,file='out/cof_masks.out',status='replace')

    write(6,'(12X,A)'), 'Defining cutoff-functions...'

    boundary_1 = 80d0
    boundary_2 = 50d0

    do j1 = 1, Nx
    do j2 = 1, Nx
      x = x0 + (j1-1) * dx
      y = x0 + (j2-1) * dx

      cofx(j1,j2) = (1d0 - cutoff_1d(j1,boundary_1)) * cutoff_1d(j2,boundary_2)

      cofy(j1,j2) = (1d0 - cutoff_1d(j2,boundary_1)) * cutoff_1d(j1,boundary_2)


      if (x*y.GE.0) then
        cofxy_same(j1,j2) = (1d0 - cutoff_1d(j1,boundary_1)) * (1d0 - cutoff_1d(j2,boundary_2)) &
          &               + (1d0 - cutoff_1d(j2,boundary_1)) * (1d0 - cutoff_1d(j1,boundary_2)) &
          &               - (1d0 - cutoff_1d(j1,boundary_1)) * (1d0 - cutoff_1d(j2,boundary_1))
      else
        cofxy_oppo(j1,j2) = (1d0 - cutoff_1d(j1,boundary_1)) * (1d0 - cutoff_1d(j2,boundary_2)) &
          &               + (1d0 - cutoff_1d(j2,boundary_1)) * (1d0 - cutoff_1d(j1,boundary_2)) &
          &               - (1d0 - cutoff_1d(j1,boundary_1)) * (1d0 - cutoff_1d(j2,boundary_1))
      end if


      cofin(j1,j2) = cutoff_1d(j1,boundary_1) * cutoff_1d(j2,boundary_1)

      if (abs(x-y).LE.15d0/au2A) mask(j1,j2) = cofxy_same(j1,j2)

      write(1,*) x*au2A, y*au2A, cofx(j1,j2), cofy(j1,j2), &
        &                      cofxy_same(j1,j2), cofxy_oppo(j1,j2),   & 
        &                      cofin(j1,j2), mask(j1,j2)

    end do
    end do


    return
  end subroutine make_new_cofs

  ! 1D cutoff function. 1/2 @ center
  function cutoff_1d(idx,center)
    use pot_param
    use global_arrays
    integer, intent(in):: idx
    double precision center, x, cutoff_1d

      x = x0 + (idx-1) * dx

      cutoff_1d = 1d0 / (1d0 + exp((abs(x)-center/au2A)/8d0))

    return
  end function
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Setup of momentum grids %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine p_grid

    use global_arrays

    implicit none
    integer I
      
    write(6,'(12X,A)'), 'Defining momenta grids.'

    do I = 1, NR  
      if (i.LE.(NR / 2)) then    
          PR(i) = (i - 1) * dpR    
      else    
          PR(i) = - (NR + 1 - i) * dpR    
      end if    
    end do

    do I = 1, Nx    
      if (i.LE.(Nx / 2)) then    
          Px(i) = (i - 1) * dpx    
      else    
          Px(i) = - (Nx + 1 - i) * dpx    
      end if    
    end do


    return
  end subroutine p_grid


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Norm calculation of real wavefunction / overlap of two wavefunctions %%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine integ_real_1D(psi1, psi2, norm)

    use data_grid

    implicit none
    integer I, N
    double precision, dimension(:), intent(in) :: psi1, psi2
    double precision, intent(out):: norm

      norm = 0.d0

      do I = 1, NR
        norm = norm + psi1(i) * psi2(i)
      end do

      norm = norm * dR

    return
  end subroutine integ_real_1D


  subroutine integ_real_2D(psi1, psi2, norm)

    use data_grid

    implicit none
    integer J1, J2, N
    double precision, dimension(:,:), intent(in) :: psi1, psi2
    double precision, intent(out):: norm

      norm = 0.d0

      do j2 = 1, Nx
      do j1 = 1, Nx
        norm = norm + psi1(j1,j2) * psi2(j1,j2)
      end do
      end do

      norm = norm *dx**2

    return
  end subroutine integ_real_2D


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Norm calculation of complex wavefunction 2d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine integ_complex_2d(psi, norm)

    use data_grid

    implicit none
    integer J1,J2
    complex*16, dimension(:,:), intent(in):: psi
    double precision, intent(out):: norm

      norm = 0.d0

      do J1 = 1, Nx
      do J2 = 1, Nx
        norm = norm + abs(psi(j1,j2))**2
      end do
      end do

      norm = norm *dx**2

    return
  end subroutine integ_complex_2d


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Norm calculation of complex wavefunction 3d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine integ_complex(psi, norm)

    use data_grid

    implicit none
    integer I,J1,J2
    complex*16, dimension(:,:,:), intent(in):: psi
    double precision, intent(out):: norm

      norm = 0.d0

      do J2 = 1, Nx
      do J1 = 1, Nx
      do  I = 1, NR
        norm = norm + abs(psi(i,j1,j2))**2
      end do
      end do
      end do

      norm = norm *dR *dx**2

    return
  end subroutine integ_complex


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Norm calculation of complex wavefunction 3d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine integ_complex_p(psi, norm)

    use data_grid

    implicit none
    integer I,J1,J2
    complex*16, dimension(:,:,:), intent(in):: psi
    double precision, intent(out):: norm

      norm = 0.d0

      do J2 = 1, Nx
      do J1 = 1, Nx
      do  I = 1, NR
        norm = norm + abs(psi(i,j1,j2))**2
      end do
      end do
      end do

      norm = norm *dpR *dpx**2

    return
  end subroutine integ_complex_p


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Norm calculation of complex wavefunction 3d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine integ_complex_xxp(psi, norm)

    use data_grid

    implicit none
    integer I,J1,J2
    complex*16, dimension(:,:,:), intent(in):: psi
    double precision, intent(out):: norm

      norm = 0.d0

      do J2 = 1, Nx
      do J1 = 1, Nx
      do  I = 1, NR
        norm = norm + abs(psi(i,j1,j2))**2
      end do
      end do
      end do

      norm = norm *dR *dx *dpx

    return
  end subroutine integ_complex_xxp


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Calculation of Eigenenergy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine eigenvalue_real_1D(psi1, psi2, E, dt_aux)

    use data_grid

    implicit none
    double precision, dimension(:), intent(in):: psi1, psi2
    double precision, intent(in) :: dt_aux
    double precision, intent(out):: E
    double precision:: e1, e2, norm

      call integ_real_1D(psi1, psi1, norm)
      e1 = norm

      call integ_real_1D(psi2, psi2, norm)
      e2 = norm

      E = (-0.5d0/dt_aux) * log(e1/e2)

    return
  end subroutine eigenvalue_real_1D


  subroutine eigenvalue_real_2D(psi1, psi2, E, dt_aux)

    use data_grid

    implicit none
    double precision, dimension(:,:), intent(in):: psi1, psi2
    double precision, intent(in) :: dt_aux
    double precision, intent(out):: E
    double precision:: e1, e2, norm

      call integ_real_2D(psi1, psi1, norm)
      e1 = norm

      call integ_real_2D(psi2, psi2, norm)
      e2 = norm

      E = (-0.5d0/dt_aux) * log(e1/e2)

    return
  end subroutine eigenvalue_real_2D


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Calculation of spectral distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine spectrum(psi, pes)

    use data_grid

    implicit none
    integer I,J1,J2
    double precision,intent(out),dimension(:,:)  :: pes
    complex*16,      intent(in) ,dimension(:,:,:):: psi

    pes = 0.d0

    do J2 = 1, Nx
    do J1 = 1, Nx
    do I = 1, NR      
        pes(j1,j2) = pes(j1,j2) + abs(psi(i,j1,j2))**2              
    end do
    end do
    end do
 
    pes = pes * dR

    return
  end subroutine spectrum


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Calculation of population distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine pop_analysis(psi, time, T)

    use pot_param
    use global_arrays

    implicit none
    integer I,J1,J2,N

    complex*16,        intent(in):: psi(NR,Nx,Nx)
    double precision,  intent(in):: time
    integer,           intent(in):: T

    double precision R
    double precision pop(Nstates), pop_super(2) ! population of regular & superposition states

    double precision sup_state_plus(NR,Nx,Nx)  ! superposition state: plus
    double precision sup_state_minus(NR,Nx,Nx) ! superposition state: minus

    complex*16 nwf_elec(NR,Nstates)  ! nuclear wave function in electronic states
    complex*16 nwf_sup(NR,2)         ! nuclear wave function in superposition states
    complex*16 c01, c12              ! <c0|c1> | <c1|c2> : nuclear correlation function


    pop       = 0d0
    pop_super = 0d0
    nwf_elec  = ( 0d0 , 0d0 )
    nwf_sup   = ( 0d0 , 0d0 )
    c01       = 0d0
    c12       = 0d0

    ! Defining superposition states
    sup_state_plus  = 1d0/sqrt(2d0) * ( ewf(:,:,:,1) + ewf(:,:,:,2) )
    sup_state_minus = 1d0/sqrt(2d0) * ( ewf(:,:,:,1) - ewf(:,:,:,2) )

    ! Projecting wf onto electronic eigenstates (= nuclear wf in electronic states)
    do N = 1, Nstates
      do J2 = 1, Nx
      do J1 = 1, Nx
      do  I = 1, NR
        nwf_elec(i,n) = nwf_elec(i,n) + ewf(i,j1,j2,n) * psi(i,j1,j2)
      end do
      end do
      end do
    end do
    nwf_elec = nwf_elec *dx**2

    ! Projecting wf onto superposition states
    do J2 = 1, Nx
    do J1 = 1, Nx
    do  I = 1, NR
        nwf_sup(i,1) = nwf_sup(i,1) + sup_state_plus(i,j1,j2)  * psi(i,j1,j2)
        nwf_sup(i,2) = nwf_sup(i,2) + sup_state_minus(i,j1,j2) * psi(i,j1,j2)
    end do
    end do
    end do
    nwf_sup = nwf_sup *dx**2

    ! Calculation of population in electronic states
    do N = 1, Nstates
      do I = 1, NR
        pop(n) = pop(n) + abs(nwf_elec(i,n))**2
      end do
    end do
    pop = pop *dR

    ! Calculation of population in superposition states
    do N = 1, 2
      do I = 1, NR
        pop_super(n) = pop_super(n) + abs(nwf_sup(i,n))**2
      end do
    end do
    pop_super = pop_super *dR

    ! Output of populations
    write(12,'(*(E15.7))') sngl(time *au2fs), (sngl(pop(n)),       N = 1, Nstates)
    ! write(13,'(*(E15.7))') sngl(time *au2fs), (sngl(pop_super(n)), N = 1, 2)

    ! Output of nuclear wavefunctions
    if (mod(T,100).EQ.0) then
      do I = 1, NR
        R = R0 + (i - 1) *dR
        write(14,'(*(E15.7))') sngl(time *au2fs), sngl(R *au2A), (sngl(abs(nwf_elec(i,n))**2), N = 1, Nstates)
        ! write(15,'(*(E15.7))') sngl(time *au2fs), sngl(R *au2A), (sngl(abs(nwf_sup(i,n))**2),  N = 1, 2)
      end do
      write(14,*)
      write(15,*)
    end if

    ! Calculation of nuclear correlation function between first two elec. eigenstates
    do I = 1, NR
      c01 = c01 + nwf_elec(i,1) * conjg(nwf_elec(i,2))
      c12 = c12 + nwf_elec(i,2) * conjg(nwf_elec(i,3))
    end do
    c01 = c01 *dR
    c12 = c12 *dR

    write(16,*) sngl(time *au2fs), sngl(abs(c01)), sngl(abs(c12))

    return
  end subroutine pop_analysis


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Density calculation for nuclear and electronic part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine density(psi, densR, densx)

    use data_grid

    implicit none
    integer I,J1,J2
    complex*16,       intent(in) :: psi(NR,Nx,Nx)
    double precision, intent(out):: densR(NR), densx(Nx)
      densR = 0d0
      densx = 0d0

      do J2 = 1, Nx
      do J1 = 1, Nx
      do  I = 1, NR
        densR(i) = densR(i) + abs(psi(i,j1,j2))**2
      end do
      end do
      end do
      densR = densR *dx**2

      do J2 = 1, Nx
      do J1 = 1, Nx
      do  I = 1, NR
        densx(j2) = densx(j2) + abs(psi(i,j1,j2))**2
      end do
      end do
      end do
      densx = densx *dR *dx

    return
  end subroutine density


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Density calculation electronic part in momentum space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine density_p(psi_p, densx_p)

    use data_grid

    implicit none
    integer I,J1,J2
    complex*16,       intent(in) :: psi_p(NR,Nx,Nx)
    double precision, intent(out):: densx_p(Nx)
      densx_p = 0d0

      do J2 = 1, Nx
      do J1 = 1, Nx
      do  I = 1, NR
        densx_p(j2) = densx_p(j2) + abs(psi_p(i,j1,j2))**2
      end do
      end do
      end do
      densx_p = densx_p *dpR *dpx

    return
  end subroutine density_p


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Calculation of interaction matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine pulse(tout,mu,E)

    use data_grid, only:dt
    use data_au, only:im

    implicit none

    integer:: J, nrot
    double precision, intent(in):: mu(3), E
    double precision:: w, u(2,2), d(2), v(2,2)
    complex*16:: b(2,2), z(2,2)
    complex*16,intent(out):: tout(2,2)

      ! Dieses Programm berechnet die Diagonalmatrix fuer die potentielle
      ! Propagation. Dabei ist die 2x2 Matrix u die Potentialmatrix, wie sie
      ! im Hamiltonoperator auftaucht:
          ! (V_1(R)-mu_g(R)* E(t)   -mu_ge(R)*E(t))
          ! (-mu_ge(R)*E(t)               V_2(R) -mu_u(R)* E(t) ) 
      ! Diese 2x2 Matrix wird mithilfe einer Jacobi-Diagonalisierung in die
      ! Diagonalform umgewandelt. (Alternative: analytische Loesung)
      ! Diagonalisierung der WW-Matrix 
   
      u(1,1) = -mu(1) * E
      u(1,2) = -mu(3) * E
      
      u(2,1) = -mu(3) * E
      u(2,2) = -mu(2) * E

 
      ! Man erhaelt die Matrix der Eigenvektoren, D, auch eine 2x2 Matrix. 
      ! Die Diagonalmatrix wird berechnet, indem man die Matrixmultiplikation
      ! D * U * D^(-1) durchfuehrt. D^(-1) ist die inverse Matrix.
      
      ! Das gilt natuerlich analog, wenn die Matrix D im Exponenten steht -
      ! dazu wird hier eine neue 2x2 Matrix B definiert, die den Propagationsschritt
      ! exp(-i * dt * D) darstellt.
      
      ! Die resultierende 2x2 Diagonalmatrix wird dann auf die Wellenfunktion 
      ! multipliziert - das beendet den Schritt mit der potentiellen Propagation
      
      call jacobi(u,2,d) 
 

      b= (0.d0,0.d0)
       
      do J = 1,2
        b(J,J) = cdexp(-im * dt * d(J))   ! Diagonalmatrix im Exponenten
      end do
       
      
      z = matmul(u,b)
      tout = matmul(z,transpose(u))

    return
  end subroutine


  subroutine jacobi (mat,dim,ewerte)

    implicit       none

    real*8         genau
    parameter      (genau=1.d-15)

    integer        Jmax,mmax
    parameter      (Jmax=15,mmax=18)
    integer        matdim
    parameter      (matdim=2)

    real*8         mat(matdim,matdim)
    integer        dim
    real*8         ewerte(matdim)

    real*8         s(matdim,matdim)
    integer        ca,cb,p,q
    real*8         c1,c2,t1,t2,t3,v1,v2,v3
    real*8         tmp,l,n,t,m1,w,m
    logical        flag

      s= 0.d0

      do 1 ca=1,dim,1
        s(ca,ca)=1.d0
1           continue

         l=0.d0
         do 2 ca=2,dim,1
            do 2 cb=1,dim,1
               tmp=mat(ca,cb)
               l=l+2.d0*tmp*tmp
2              continue

         n=dsqrt(l)
         m=genau*n/dim
         t=n

3        t=t/dim
4           do 6 q=2,dim,1
               do 6 p=1,q-1,1
                  flag=.false.
                  if (dabs(mat(p,q)).gt.t) then
                     flag=.true.
                     v1=mat(p,p)
                     v2=mat(p,q)
                     v3=mat(q,q)
                     m1=(v1-v3)/2.d0
                     if (m1.eq.0.d0) then
                           w=-1.d0
                        else
                           if (m1.gt.0.d0) then
                                 w=-v2/(dsqrt(v2*v2+m1*m1))
                              else
                                 w=v2/(dsqrt(v2*v2+m1*m1))
                              endif
                        endif

                     t1=w/dsqrt(2.d0*(1+dsqrt(1.d0-w/2.d0)))
                     t2=t1*t1
                     c1=dsqrt(1.d0-t2)
                     c2=c1*c1
                     t3=t1*c1

                     do 7 ca=1,dim,1
                        l=mat(ca,p)*c1-mat(ca,q)*t1
                        mat(ca,q)=mat(ca,p)*t1+mat(ca,q)*c1
                        mat(ca,p)=l
                        l=s(ca,p)*c1-s(ca,q)*t1
                        s(ca,q)=s(ca,p)*t1+s(ca,q)*c1
                        s(ca,p)=l
7                       continue
                     do 8 ca=1,dim,1
                        mat(p,ca)=mat(ca,p)
                        mat(q,ca)=mat(ca,q)
8                       continue
                     mat(p,p)=v1*c2+v3*t2-2*v2*t3
                     mat(q,q)=v1*t2+v3*c2+2*v2*t3
                     tmp=(v1-v3)*t3+v2*(c2-t2)
                     mat(p,q)=tmp
                     mat(q,p)=tmp
                     end if
6                 continue
               if (flag) go to 4
            if (m.lt.t) go to 3
ewerte=0.d0
         do 9 ca=1,dim,1
            ewerte(ca)=mat(ca,ca)
9           continue
         do 10 ca=1,dim,1
            do 10 cb=1,dim,1
               mat(ca,cb)=s(ca,cb)
10             continue

    return
  end subroutine jacobi


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Calculation of ionic dipole moment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine ionic_dipole

    ! use data_au
    ! use pot_param
    ! use global_arrays

    ! implicit none

    ! integer i,j,k
    ! double precision:: x,R

    ! open(1,file='out/ionic_dipole.out',status='replace')

    ! mu_ionic = 0d0

    ! do i = 1, NR
    !   R = R0 + (i-1) * dR

    !   do k = Nx/2+1, Nx
    !   do j = 1, Nx
    !     x = x0 + (j-1) * dx
    !     mu_ionic(i,k) = mu_ionic(i,k) + ewf(i,j,1) * x * nwf(i,Vstates+1) * exp(im * px(k) * x)
    !   end do
    !   write(1,*) sngl(R*au2A), sngl(Px(k)), real(real(mu_ionic(i,k))), real(aimag(mu_ionic(i,k)))
    !   end do

    !   do k = 1, Nx/2
    !   do j = 1, Nx
    !     x = x0 + (j-1) * dx
    !     mu_ionic(i,k) = mu_ionic(i,k) + ewf(i,j,1) * x * nwf(i,Vstates+1) * exp(im * px(k) * x)
    !   end do
    !   write(1,*) sngl(R*au2A), sngl(Px(k)), real(real(mu_ionic(i,k))), real(aimag(mu_ionic(i,k)))
    !   end do
    !   write(1,*)
    !   write(1,*)

    ! end do

    ! close(1,status='keep')

  end subroutine ionic_dipole


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%% Calculation of ionic dipole moment (mu_ionic(R,E) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine ionic_dipole_e

    ! use data_au
    ! use pot_param
    ! use global_arrays

    ! implicit none

    ! integer i,j
    ! integer idx
    ! double precision dif, difnew, Ekin(Nx/2), en

    ! do j = 1, Nx/2
    !   Ekin(j) = Px(j)**2 / 2
    ! end do

    ! do i = 1, Nx/2
    !   dif = 10
    !   en = E0 + (i-1) * dE

    !   do j = 1, Nx/2
    !     difnew = abs(Ekin(j) - en)
    !     if ( difnew .LE. dif ) then
    !       dif = difnew
    !       idx = j
    !     end if
    !   end do

    !   mu_ionic_e(:,i) = mu_ionic(:,idx)
    ! end do

  end subroutine ionic_dipole_e

end module