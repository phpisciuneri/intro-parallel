!
! First order upwind method for the one dimensional convection-diffusion equation
!
program main
  implicit none

  integer, parameter :: PREC = kind( 1.0d0 ) ! double precision 
  integer, parameter :: NX = 5001             ! grid points
  
  real(PREC) :: pi = 4 * atan( 1._PREC ) ! pi
  real(PREC) :: L                        ! domain length
  real(PREC) :: nu = 0._PREC             ! viscosity
  real(PREC) :: U = 1._PREC              ! velocity
  real(PREC) :: dt                       ! time step
  real(PREC) :: x                        ! x location
  real(PREC) :: dx                       ! grid spacing
  real(PREC) :: c = .6_PREC              ! CFL number
  real(PREC) :: d                        ! fourier number
  real(PREC) :: R                        ! cell reynolds number

  real(PREC) :: Tn(NX)   = 0._PREC       ! T at time step n
  real(PREC) :: Tnp1(NX) = 0._PREC       ! T at time step n+1

  real(PREC) :: eps = epsilon( 1._PREC ) ! machine precision

  integer :: i                           ! spatial loop counter
  integer :: n                           ! temporal loop counter
  integer :: max_iter                    ! max iterations
  
  ! +----------------+
  ! | INITIALIZATION |
  ! +----------------+
  L       = 2*pi
  dx      = L / ( NX - 1 )
  dt      = ( c - eps )*dx / U
  d       = nu*dt / ( dx*dx )
  max_iter = floor( L / ( U*dt ) + eps )
  
  ! tophat begins when x=1, and is 25% of the domain length
  do i=1,NX
     x = ( i - 1 )*dx
     if ( x > ( 1 - eps ) .and. x < ( 1 + .25_PREC*L + eps ) ) then
        Tn(i) = 1._PREC
     end if
  end do

  ! +---------------+
  ! | UPWIND SOLVER |
  ! +---------------+
  do n=1,max_iter
     
     ! left-most point
     Tnp1(1) = Tn(1) - c*( Tn(1) - Tn(NX-1) ) &
          + d*( Tn(2) - 2*Tn(1) + Tn(NX-1) )     

     ! interior
     do i=2,NX-1
        Tnp1(i) = Tn(i) - c*( Tn(i) - Tn(i-1) ) &
             + d*( Tn(i+1) - 2*Tn(i) + Tn(i-1) )
     end do

     ! right-most point (periodic)
     Tnp1(NX) = Tnp1(1)

     ! prepare for next iteration
     do i=1,NX
        Tn(i) = Tnp1(i)
     end do

  end do

  call write_profile( Tn )


contains

  ! +---------------+
  ! | WRITE_PROFILE |
  ! +---------------+
  subroutine write_profile( T )
    real(PREC), intent(in) :: T(NX)

    real(PREC) :: x
    integer    :: i
    
    open( unit=101, file='profile.txt' )
    do i=1,NX
       x = ( i - 1 )*dx
       write( 101, * ) x, T(i)
    end do
    close(101)
    
  end subroutine write_profile
  
end program main
