program sswe1d

  ! Simple Shallow water equation 1D
  ! Fatkhurrozi - Ocean Engineering ITERA - Jul 2025
	
  implicit none

  ! Parameters
  integer, parameter :: nx = 201        
  integer, parameter :: nt = 400        
  real(8), parameter :: dx = 1.0d0      
  real(8), parameter :: dt = 0.1d0      
  real(8), parameter :: g = 9.81d0      
  real(8), parameter :: de = 10.0d0      ! if proper
  
  real(8) :: detadta,dudta, detadtn,dudtn
  real(8), dimension(nx, nt) :: eta, u
  real(8), dimension(nx) :: etap, up, detadt, dudt, dd
  integer :: i, n
  real(8) :: x
  integer :: outf

! -- initialization
  do i = 1, nx
     x = (i - 1) * dx
     eta(i,1)=0.
	u(i,1)=0.
	dd(i)=5-i*0.02		! variying depth
  end do


! -- time integration
 
!boundary cond. 
 do n = 1, nt - 1
     eta(1,n)=0.5*sin(2*3.14/5*n*dt)
	 !print*, eta(1,n)
enddo 
print *, ' -- SSWE1D --'
print *, 'Ocean Engineering ITERA - Jul 2025'

  do n = 1, nt-1
  ! Predictor  
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
		
! Corrector s
     do i=2,nx-1
		detadtn		= -dd(i)*(up(i)-up(i-1))/dx
		dudtn		= -g*(etap(i)-etap(i-1))/dx 
		
		detadta		= (detadt(i)+detadtn)*0.5
		dudta		= (dudt(i)+dudtn)*0.5
		
		eta(i,n+1)	= eta(i,n)+detadta*dt
		u(i,n+1)	= u(i,n)+dudta*dt
	 enddo
	 
! Neumann bc
     u(nx,n+1) = u(nx - 1,n+1)
     eta(nx,n+1) = eta(nx - 1,n+1)

  end do

  print *, '.'
  print *, '.'
  
! -- Write u(x,t) to file
 
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

  print *, ' -- Success! -- Printed to uoutput.txt.'

end program sswe1d
