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
      procedure :: v => gtf_value
      procedure :: dx => gtf_derivx
      procedure :: dy => gtf_derivy
      procedure :: dz => gtf_derivz
  end type

  type :: basis_function
    type(gtf), dimension(:), allocatable :: gf ! gauss function
    real(kind=r8), dimension(:), allocatable :: cf !contraction factor

    contains
      procedure :: v => basis_value
      procedure :: dx => basis_derivx
      procedure :: dy => basis_derivy
      procedure :: dz => basis_derivz
  end type

  type :: molecule
    integer :: num_elec ! number of electrons
    integer :: num_alpha ! number of alpha electrons
    integer :: num_beta ! number of beta electrons
    character(len=2) :: wf_type ! wave function type: 'R', 'U', 'RO', 'NO'
    type(basis_function), dimension(:), allocatable :: basis
    real(kind=r8), dimension(:), allocatable :: num_occ
    real(kind=r8), dimension(:,:), allocatable :: amo_coeff ! C(imo,ibss)
    real(kind=r8), dimension(:,:), allocatable :: bmo_coeff
    contains
      procedure :: density
      procedure :: den_grad
  end type

  contains
  
  real(kind=r8) function density(this, x, y, z, dtype)
    class(molecule), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z
    character(len=1), intent(in) :: dtype ! density type
    !---------------------------------------------------------------------------
    real(kind=r8) :: rhoa, rhob
    
    if (dtype == 't') then

      if (this%wf_type == 'U ' .or. this%wf_type == 'RO') then
        rhoa = rho(this, x, y, z, 'a')
        rhob = rho(this, x, y, z, 'b')
        density = rhoa + rhob
      else if (this%wf_type == 'R ') then
        rhoa = rho(this, x, y, z, 'a')
        density = 2.0_r8 * rhoa
      else if (this%wf_type == 'NO') then
        density = rho(this, x, y, z, 'n')
      end if

      return

    end if

    if (this%wf_type == 'NO') stop "routine density: alpha/beta density is not available for natural orbital"
    if (dtype == 'a') then
      density = rho(this, x, y, z, 'a')
    else if (dtype == 'b') then
      density = rho(this, x, y, z, 'b')
    end if
  
  end function

  real(kind=r8) function rho(this, x, y, z, dtype)
    class(molecule), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z
    character(len=1), intent(in) :: dtype
    !---------------------------------------------------------------------------
    real(kind=r8), dimension(:), allocatable :: amo, bmo ! alpha/beta occupied orbitals
    integer :: imo
  
    if (dtype == 'a') then

      allocate( amo(this%num_alpha) )
      do imo = 1, this%num_alpha
        amo(imo) = sum( this%amo_coeff(imo,:) * this%basis(:)%v(x,y,z) )
      end do
      rho = sum( amo**2 )

    else if (dtype == 'b') then

      allocate( bmo(this%num_beta) )
      do imo = 1, this%num_beta
        bmo(imo) = sum( this%bmo_coeff(imo,:) * this%basis(:)%v(x,y,z) )
      end do
      rho = sum( bmo**2 )

    else if (dtype == 'n') then

      allocate( amo( size(this%num_occ) ) )
      do imo = 1, size(this%num_occ)
        amo(imo) = sum( this%amo_coeff(imo,:) * this%basis(:)%v(x,y,z) )
      end do
      rho = sum( amo**2 * this%num_occ )

    end if

  end function

  real(kind=r8) function den_grad(this, x, y, z, gtype)
    class(molecule), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z
    character(len=1), intent(in) :: gtype ! component type
    !-------------------------------------------------------
    real(kind=r8) :: gradx, grady, gradz

    if (gtype == 'm') then
      gradx = grad(this, x, y, z, 'x')
      grady = grad(this, x, y, z, 'y')
      gradz = grad(this, x, y, z, 'z')
      den_grad = sqrt(gradx**2 + grady**2 + gradz**2)
    else if (gtype == 'x') then
      den_grad = grad(this, x, y, z, 'x')
    else if (gtype == 'y') then
      den_grad = grad(this, x, y, z, 'y')
    else if (gtype == 'z') then
      den_grad = grad(this, x, y, z, 'z')
    end if

  end function

  real(kind=r8) function grad(this, x, y, z, gtype)
    class(molecule), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z
    character(len=1), intent(in) :: gtype ! component type
    !-------------------------------------------------------
    real(kind=r8), dimension(:), allocatable :: amo, bmo ! alpha/beta occupied orbitals
    real(kind=r8), dimension(:), allocatable :: amox, bmox
    real(kind=r8), dimension(:), allocatable :: amoy, bmoy
    real(kind=r8), dimension(:), allocatable :: amoz, bmoz
    integer :: imo

    if (this%wf_type == 'NO') then
      
      allocate( amo( size(this%num_occ) ) )
      do imo = 1, size(this%num_occ)
        amo(imo) = sum( this%amo_coeff(imo,:) * this%basis(:)%v(x,y,z) )
      end do

      if (gtype == 'x') then
        allocate( amox( size(this%num_occ) ) )
        do imo = 1, size(this%num_occ)
          amox(imo) = sum( this%amo_coeff(imo,:) * this%basis(:)%dx(x,y,z) )
        end do

        grad = sum( 2.0 * this%num_occ * amo * amox )

      else if (gtype == 'y') then
        allocate( amoy( size(this%num_occ) ) )
        do imo = 1, size(this%num_occ)
          amoy(imo) = sum( this%amo_coeff(imo,:) * this%basis(:)%dy(x,y,z) )
        end do

        grad = sum( 2.0 * this%num_occ * amo * amoy )

      else if (gtype == 'z') then
        allocate( amoz( size(this%num_occ) ) )
        do imo = 1, size(this%num_occ)
          amoz(imo) = sum( this%amo_coeff(imo,:) * this%basis(:)%dz(x,y,z) )
        end do

        grad = sum( 2.0 * this%num_occ * amo * amoz )

      end if

    else

      allocate( amo(this%num_alpha) )
      do imo = 1, this%num_alpha
        amo(imo) = sum( this%amo_coeff(imo,:) * this%basis(:)%v(x,y,z) )
      end do

      allocate( bmo(this%num_beta) )
      do imo = 1, this%num_beta
        bmo(imo) = sum( this%bmo_coeff(imo,:) * this%basis(:)%v(x,y,z) )
      end do

      if (gtype == 'x') then
        allocate( amox(this%num_alpha) )
        do imo = 1, this%num_alpha
          amox(imo) = sum( this%amo_coeff(imo,:) * this%basis(:)%dx(x,y,z) )
        end do
  
        allocate( bmox(this%num_beta) )
        do imo = 1, this%num_beta
          bmox(imo) = sum( this%bmo_coeff(imo,:) * this%basis(:)%dx(x,y,z) )
        end do

        grad = sum( 2.0 * amo * amox ) + sum( 2.0 * bmo * bmox )

      else if (gtype == 'y') then
        allocate( amoy(this%num_alpha) )
        do imo = 1, this%num_alpha
          amoy(imo) = sum( this%amo_coeff(imo,:) * this%basis(:)%dy(x,y,z) )
        end do

        allocate( bmoy(this%num_beta) )
        do imo = 1, this%num_beta
          bmoy(imo) = sum( this%bmo_coeff(imo,:) * this%basis(:)%dy(x,y,z) )
        end do

        grad = sum( 2.0 * amo * amoy ) + sum( 2.0 * bmo * bmoy )

      else if (gtype == 'z') then
        allocate( amoz(this%num_alpha) )
        do imo = 1, this%num_alpha
          amoz(imo) = sum( this%amo_coeff(imo,:) * this%basis(:)%dz(x,y,z) )
        end do

        allocate( bmoz(this%num_beta) )
        do imo = 1, this%num_beta
          bmoz(imo) = sum( this%bmo_coeff(imo,:) * this%basis(:)%dz(x,y,z) )
        end do

        grad = sum( 2.0 * amo * amoz ) + sum( 2.0 * bmo * bmoz )

      end if

    end if

  end function

  elemental real(kind=r8) function basis_value(this, x, y, z)
    class(basis_function), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z

    basis_value = sum( this%gf(:)%v(x,y,z) * this%cf(:) )

  end function

  elemental real(kind=r8) function basis_derivx(this, x, y, z)
    class(basis_function), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z

    basis_derivx = sum( this%gf(:)%dx(x,y,z) * this%cf(:) )

  end function

  elemental real(kind=r8) function basis_derivy(this, x, y, z)
    class(basis_function), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z

    basis_derivy = sum( this%gf(:)%dy(x,y,z) * this%cf(:) )

  end function

  elemental real(kind=r8) function basis_derivz(this, x, y, z)
    class(basis_function), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z

    basis_derivz = sum( this%gf(:)%dz(x,y,z) * this%cf(:) )

  end function

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
    
    this%Norm = (2.*zeta/pi)**0.75 * sqrt( (8.*zeta)**(lx+ly+lz) &
                    * f(lx)*f(ly)*f(lz) / (f(2*lx)*f(2*ly)*f(2*lz)) )

  end subroutine
  
  elemental real(kind=r8) function gtf_value(this, x, y, z)
    class(gtf), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z

    associate(rx => x-this%x, &
              ry => y-this%y, &
              rz => z-this%z)

    gtf_value = this%Norm * rx**this%lx &
              * ry**this%ly * rz**this%lz &
              * exp( -this%zeta * (rx**2 + ry**2 + rz**2) )
    
    end associate

  end function

  elemental real(kind=r8) function gtf_derivx(this, x, y, z)
    class(gtf), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z
  
    associate(rx => x-this%x, &
              ry => y-this%y, &
              rz => z-this%z)

    if (this%lx == 0) then
      gtf_derivx = -2.0*this%zeta * rx &
                 * this%Norm * ry**this%ly * rz**this%lz &
                 * exp( -this%zeta * (rx**2 + ry**2 + rz**2) )
    else if (this%lx > 0) then
      gtf_derivx = (this%lx * rx**(this%lx-1) - 2.0*this%zeta * rx**(this%lx+1)) &
                 * this%Norm * ry**this%ly * rz**this%lz &
                 * exp( -this%zeta * (rx**2 + ry**2 + rz**2) )
    end if

    end associate

  end function

  elemental real(kind=r8) function gtf_derivy(this, x, y, z)
    class(gtf), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z
  
    associate(rx => x-this%x, &
              ry => y-this%y, &
              rz => z-this%z)
  
    if (this%ly == 0) then
      gtf_derivy = -2.0*this%zeta * ry &
                 * this%Norm * rx**this%lx * rz**this%lz &
                 * exp( -this%zeta * (rx**2 + ry**2 + rz**2) )
    else if (this%ly > 0) then
      gtf_derivy = (this%ly * ry**(this%ly-1) - 2.0*this%zeta * ry**(this%ly+1)) &
                 * this%Norm * rx**this%lx * rz**this%lz &
                 * exp( -this%zeta * (rx**2 + ry**2 + rz**2) )
    end if
    
    end associate

  end function

  elemental real(kind=r8) function gtf_derivz(this, x, y, z)
    class(gtf), intent(in) :: this
    real(kind=r8), intent(in) :: x, y, z
  
    associate(rx => x-this%x, &
              ry => y-this%y, &
              rz => z-this%z)
  
    if (this%lz == 0) then
      gtf_derivz = -2.0*this%zeta * rz &
                 * this%Norm * rx**this%lx * ry**this%ly &
                 * exp( -this%zeta * (rx**2 + ry**2 + rz**2) )
    else if (this%lz > 0) then
      gtf_derivz = (this%lz * rz**(this%lz-1) - 2.0*this%zeta * rz**(this%lz+1)) &
                 * this%Norm * rx**this%lx * ry**this%ly &
                 * exp( -this%zeta * (rx**2 + ry**2 + rz**2) )
    end if
    
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

    if (n < 0) stop "routine factorial: negative number input"
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
      if ( ierr /= 0) stop "routinte locate_label: label not found"
      if ( index(adjustl(char_tmp), label) == 1 ) exit
    end do
    backspace(iounit)

  end subroutine

end module