program sswe1d

  ! Simple Shallow water equation 1D
	
  implicit none

  ! Parameters
  integer, parameter :: nx = 201        ! number of spatial points
  integer, parameter :: nt = 400        ! number of time steps
  real(8), parameter :: dx = 1.0d0      ! spatial resolution
  real(8), parameter :: dt = 0.1d0      ! time step
  real(8), parameter :: g = 9.81d0      ! gravitational acceleration
  real(8), parameter :: de = 10.0d0      ! constant water depth
  
  real(8) :: detadta,dudta, detadtn,dudtn

  ! Arrays to store the full time history
  real(8), dimension(nx, nt) :: eta, u

  ! Temporary arrays for predictor step
  real(8), dimension(nx) :: etap, up, detadt, dudt, dd

  ! Loop variables
  integer :: i, n
  real(8) :: x

  ! File output
  integer :: outf

  ! ------------------------------
  ! Initialization at t = 0
  ! ------------------------------
  !eta(:,1) = 0.0d0
  !u(:,1) = 0.0d0

  do i = 1, nx
     x = (i - 1) * dx
     eta(i,1)=0.
	 u(i,1)=0.
	 dd(i)=5-i*0.02
	 !eta(i,1) = exp( -((x - 50.0d0)**2) / (2.0d0 * 25.0d0) )
  end do

  ! ------------------------------
  ! Time stepping loop
  ! ------------------------------
 
!boundary cond. 
 do n = 1, nt - 1
     eta(1,n)=0.5*sin(2*3.14/5*n*dt)
	 !print*, eta(1,n)
enddo 
	 
  do n = 1, nt-1

	 ! Predictor step (forward difference)
     do i = 1, nx
        detadt(i) 	= -dd(i)*(u(i+1,n) - u(i,n))/ dx
		dudt(i) 	= -g*(eta(i+1,n) - eta(i,n))/dx
		up(i) 		= u(i,n) + dudt(i)* dt
		etap(i) 	= eta(i,n) + detadt(i)* dt
		
		!print*, dd(i)
     end do

		etap(nx) 	= etap(nx - 1)
		up(nx)   	= up(nx - 1)
	
	 ! print*, up(i),etap(i)
		
     ! Corrector step (backward difference)
     do i=2,nx-1
		detadtn		= -dd(i)*(up(i)-up(i-1))/dx
		dudtn		= -g*(etap(i)-etap(i-1))/dx 
		
		detadta		= (detadt(i)+detadtn)*0.5
		dudta		= (dudt(i)+dudtn)*0.5
		
		eta(i,n+1)	= eta(i,n)+detadta*dt
		u(i,n+1)	= u(i,n)+dudta*dt
	 enddo
	 
! 	Neumann boundary condition
     u(nx,n+1) = u(nx - 1,n+1)
     eta(nx,n+1) = eta(nx - 1,n+1)


  end do

  ! ------------------------------
  ! Write u(x,t) to file
  ! ------------------------------
  outf=10
  open(unit=outf, file="uoutput.txt", status="replace", action="write")

  do n = 1, nt
     !write(outf,*) 't =', (n-1)*dt
     do i = 1, nx
        write(outf,'(F10.2,2x,F10.3,2x,f8.2,3X,F10.3,3X,F10.3)') (n-1)*dt,(i-1)*dx,dd(i),u(i,n),eta(i,n)
     end do
     !write(outf,*)  ! blank line
  end do

  close(outf)

  print *, 'Finished. Output written to uoutput.txt.'

end program sswe1d
