module m_basis_func
  use iso_fortran_env, only: r8 => real64
  use m_general, f => factorial
  implicit none
  private
  public :: molecule

  type :: gtf
  ! gtf = N * x^i * y^j * z^k * exp( -zeta * r^2 )
    real(kind=r8) :: x, y, z
    integer :: lx, ly, lz ! i, j, k 
    real(kind=r8) :: zeta
    real(kind=r8) :: Norm ! N: normalization factor

    contains
      procedure :: initialize => initialize_gtf
      procedure :: v => value_gtf
  end type

  type :: basis_function
    type(gtf), dimension(:), allocatable :: gf ! gauss function
    real(kind=r8), dimension(:), allocatable :: cf !contraction factor
  end type

  type :: molecule
    integer :: num_elec ! number of electrons
    integer :: num_alpha ! number of alpha electrons
    integer :: num_beta ! number of beta electrons
    character(len=2) :: wf_type ! wave function type: 'R', 'U', 'RO', 'NO'
    type(basis_function), dimension(:), allocatable :: basis
    real(kind=r8), dimension(:), allocatable :: num_occ
    real(kind=r8), dimension(:,:), allocatable :: amo_coeff
    real(kind=r8), dimension(:,:), allocatable :: bmo_coeff
    contains
    !  procedure :: density
  end type

  contains
  
  !real function density(this, x, y, z, dtype)
  !  class(molecule), intent(in) :: this
  !  real(kind=r8), intent(in) :: x, y, z
  !  character(len=1), intent(in), optional :: dtype ! density type
  !
  !  if ( present(dtype) ) then
  !    if (dtype == 'a') continue ! calculate alpha density
  !    if (dtype == 'b') continue ! calculate beta density
  !  else
  !    continue ! calculate total density
  !  end if
  !
  !end function
  
  subroutine initialize_gtf(this, locate, zeta, lx, ly, lz)
    class(gtf) :: this
    real(kind=r8), intent(in) :: zeta
    real(kind=r8), dimension(3), intent(in) :: locate
    integer, intent(in) :: lx, ly, lz

    this%x = locate(1)
    this%y = locate(2)
    this%z = locate(3)
    this%zeta = zeta
    this%lx = lx
    this%ly = ly
    this%lz = lz
    
    this%Norm = (2.*zeta/pi)**0.75 * sqrt( &
                    (8.*zeta)**(lx+ly+lz) &
                    * f(lx)*f(ly)*f(lz) &
                    / (f(2*lx)*f(2*ly)*f(2*lz)) )

  end subroutine
  
  real(kind=r8) function value_gtf(this, x, y, z)
    class(gtf) :: this
    real(kind=r8) :: x, y, z

    associate(rx => x-this%x, &
              ry => y-this%y, &
              rz => z-this%z)

    value_gtf = this%Norm * rx**this%lx &
                  * ry**this%ly * rz**this%lz &
                  * exp( -this%zeta * (rx**2 + ry**2 + rz**2) )
    
    end associate

  end function

end module

module m_general
  use iso_fortran_env, only: r8 => real64
  implicit none
  real(kind=r8) :: pi = 3.141592653589793_r8
  
  contains

  integer function factorial(n)
    integer, intent(in) :: n
    integer :: iter

    factorial = 1
    do iter = 2, n
      factorial = factorial * iter
    end do

  end function

  subroutine locate_label(iounit, label)
    integer, intent(in) :: iounit
    character(len=*), intent(in) :: label
    character(len=80) :: char_tmp
    integer :: ierr

    rewind(iounit)
    do while (.true.)
      read(iounit, '(A80)', iostat=ierr) char_tmp
      if ( ierr /= 0) stop "label not found"
      if ( index(adjustl(char_tmp), label) == 1 ) exit
    end do
    backspace(iounit)

  end subroutine


end module