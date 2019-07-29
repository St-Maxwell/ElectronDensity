module m_routines
  use m_basis_func
  use iso_fortran_env, only: r8 => real64
  use m_general, only: locate_label
  implicit none
  private
  public :: read_basis_func

  interface
    module subroutine read_basis_func(iounit, mol)
      integer, intent(in) :: iounit
      type(molecule), intent(inout) :: mol
    end subroutine
  end interface

end module

submodule (m_routines) fch_routine
  !-------------------------------------------------------
  ! local variables
  ! but are global for submodule subroutines
  integer :: num_basis
  integer :: num_contra_shl
  integer :: num_primitive
  integer, dimension(:), allocatable :: shl_type
  integer, dimension(:), allocatable :: shl_ctr_odr
  real(kind=r8), dimension(:), allocatable :: prim_exp
  real(kind=r8), dimension(:), allocatable :: ctr_cff
  real(kind=r8), dimension(:), allocatable :: SP_ctr_cff
  real(kind=r8), dimension(:), allocatable :: shl_coord
  !-------------------------------------------------------

  contains
  module subroutine read_basis_func(iounit, mol)
    integer, intent(in) :: iounit
    type(molecule), intent(inout) :: mol
    !-------------------------------------------------------
    integer, dimension(-3:3), parameter :: type2nbss = [7,5,4,1,3,6,10]
    ! reference Sobereva's code: http://sobereva.com/55
    integer, dimension(-3:3,30) :: s2f
    integer, dimension(20), parameter :: type2lx = [0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,0,1,1,0,1]
    integer, dimension(20), parameter :: type2ly = [0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1]
    integer, dimension(20), parameter :: type2lz = [0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,1,0,2,2,1]
    !-------------------------------------------------------
    character(len=80) :: char_tmp
    integer :: lx, ly, lz
    real(kind=r8), dimension(3) :: locate
    integer :: ibss, ishl, ibsshl, ictr, iprim
    !-------------------------------------------------------
    s2f(0,1) = 1
    s2f(-1,1:4) = [1,2,3,4]
    s2f(1,1:3) = [2,3,4]
    s2f(2,1:6) = [5,6,7,8,9,10]
    s2f(3,1:10) = [11,12,13,17,14,15,18,19,16,20]
    !-------------------------------------------------------

    ! read number of basis functions
    call locate_label(iounit, 'Number of basis functions')
    read(iounit, '(50X,I11)') num_basis

    ! read number of contracted shells
    call locate_label(iounit, 'Shell types')
    read(iounit, '(50X,I11)') num_contra_shl

    ! read type of each shell
    allocate( shl_type(num_contra_shl) )
    read(iounit, *) shl_type ! just utilize list-directed i/o

    ! read contraction order of each shell   
    read(iounit, *) ! skip one line
    allocate( shl_ctr_odr(num_contra_shl) )
    read(iounit, *) shl_ctr_odr

    ! read number of contracted shells
    call locate_label(iounit, 'Primitive exponents')
    read(iounit, '(50X,I11)') num_primitive

    ! read exponents of primitives (GTF)
    allocate( prim_exp(num_primitive) )
    read(iounit, *) prim_exp

    ! read contraction coefficients
    read(iounit, *) ! skip one line
    allocate( ctr_cff(num_primitive) )
    read(iounit, *) ctr_cff

    ! read SP(P) contraction coefficients if exist
    read(iounit, '(A80)') char_tmp
    if ( index(adjustl(char_tmp), &
         'P(S=P) Contraction coefficients') == 1) then
      allocate( SP_ctr_cff(num_primitive) )
      read(iounit, *) SP_ctr_cff
    end if

    ! read coordinates of contracted shells
    call locate_label(iounit, 'Coordinates of each shell ')
    read(iounit, *) 
    allocate( shl_coord(3*num_contra_shl) )
    read(iounit, *) shl_coord

    !-------------------------------------------------------

    allocate( mol%basis(num_basis) )
    ibss = 1; iprim = 1

    do ishl = 1, num_contra_shl ! loop contracted shells

      do ibsshl = 1, type2nbss( shl_type(ishl) ) ! loop basis functions in a shell

        associate( basis => mol%basis(ibss), &
                   sco => shl_ctr_odr(ishl) )

        allocate( basis%cf( sco ) )
        allocate( basis%gf( sco ) )
        
        locate = shl_coord(ishl*3-2:ishl*3)
        lx = type2lx( s2f(shl_type(ishl), ibsshl) )
        ly = type2ly( s2f(shl_type(ishl), ibsshl) )
        lz = type2lz( s2f(shl_type(ishl), ibsshl) )

        do ictr = 1, sco ! loop primitive GTFs in a basis function

          if ( shl_type(ishl) == -1 .and. ibsshl /= 1 ) then
            ! sp type shell, and the last three GTFs are p-functions
            basis%cf(ictr) = SP_ctr_cff(ictr+iprim-1) 
          else
            ! sp type shell, the first is s-function
            basis%cf(ictr) = ctr_cff(ictr+iprim-1)
          end if

          call basis%gf(ictr)%initialize(locate, prim_exp(ictr+iprim-1), lx, ly, lz)
 
        end do

        end associate

        ibss = ibss + 1

      end do

      iprim = iprim + shl_ctr_odr(ishl)

    end do

  end subroutine

  !module subroutine init_mol_basis(mol)
  !  type(molecule), intent(inout) :: mol 
  !
  !end subroutine

end submodule