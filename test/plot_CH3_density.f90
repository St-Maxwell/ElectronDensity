module m_tmp
  use iso_fortran_env, only: r8 => real64
  implicit none
  contains
  function linspace(a,b,n_elements)
    ! from https://github.com/kookma/ogpf
    !..............................................................................
    !   returns a linearly spaced vector with n points in [a, b]
    !   if n is omitted, 100 points will be considered
    !..............................................................................
  
    real(r8), intent(in)           :: a
    real(r8), intent(in)           :: b
    integer,  intent(in), optional :: n_elements
    real(r8), allocatable          :: linspace(:)
  
    !   Local vars
    real(r8) :: dx
    integer  :: i
    integer  :: n
    integer  :: ierr
  
    if (present(n_elements)) then
        if (n_elements <=1 ) then
            print*, "linspace procedure: Error: wrong value of n_elements, use an n_elements > 1"
            stop
        end if
        n=n_elements
    else
        n=100
    end if
  
    allocate(linspace(n), stat=ierr)
    if (ierr /= 0) then
        print*, "linspace procedure: Fatal Error, Allocation failed in linspace function"
        stop
    end if
  
    dx=(b-a)/real((n-1),r8)
    linspace=[(i*dx+a, i=0,n-1)]
  
  end function linspace

end module
program main
  use m_basis_func
  use m_routines
  use m_tmp
  use iso_fortran_env, only: r8 => real64
  implicit none
  type(molecule) :: mol
  character(len=80) :: filename
  integer, parameter :: iounit = 12
  real(kind=r8), dimension(:), allocatable :: x, y
  real(kind=r8), parameter :: ang2br = 1.8897260_r8
  real(kind=r8) :: z = 0.0_r8
  real(kind=r8), dimension(:,:), allocatable :: rho
  integer :: iter, jter

  write(*,*) "Input .fch file"
  read(*,*) filename
  open(unit=iounit, file=filename, status='old', action='read')
  call read_basis_func(iounit, mol)
  call read_mo_coeff(iounit, mol)
  close(iounit)

  x = linspace(-2.0_r8, 2.0_r8, 100) * ang2br
  y = linspace(-2.0_r8, 2.0_r8, 100) * ang2br
  allocate( rho(size(x), size(y)) )

  do iter = 1, size(x)
    do jter = 1, size(y)
      rho(iter, jter) = mol%density(x(iter), y(jter), z, 't')
    end do
  end do

  open(unit=iounit, file='CH3_density.txt', status='new', action='write')
  write(iounit,'(100ES15.7)') ((rho(iter,jter),jter=1,100),iter=1,100)
  close(iounit)

end program