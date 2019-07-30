program main
  use m_basis_func
  use m_routines
  use iso_fortran_env, only: r8 => real64
  implicit none
  type(molecule) :: mol
  integer, parameter :: iounit = 12
  real(kind=r8) :: x, y, z

  open(unit=iounit, file='water_6-31Gdp.fch', status='old', action='read')
  call read_basis_func(iounit, mol)
  call read_mo_coeff(iounit, mol)
  close(iounit)

  do while(.true.)
		write(*,*) "Input the coordinate in Bohr   e.g. 1.0,1.2,-3.5"
		read(*,*) x,y,z
		write(*,"(' The density is',f16.10,' e/Bohr^3',/)") mol%density(x, y, z, 't')
  end do
  
end program