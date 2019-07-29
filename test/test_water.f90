program main
  use m_basis_func
  use m_routines
  implicit none
  type(molecule) :: mol
  integer, parameter :: iounit = 12

  open(unit=iounit, file='water_6-31Gdp.fch', status='old', action='read')
  call read_basis_func(iounit, mol)
  close(iounit)

end program