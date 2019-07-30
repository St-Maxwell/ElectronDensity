program main
  use m_basis_func
  use m_routines
  implicit none
  type(molecule) :: mol
  integer, parameter :: iounit = 12

  open(unit=iounit, file='NO_6-31Gd.fch', status='old', action='read')
  call read_basis_func(iounit, mol)
  call read_mo_coeff(iounit, mol)
  close(iounit)

end program