module settings_mod
use, intrinsic :: iso_fortran_env
implicit none
public
logical, parameter :: LOMPTEST = .false. ! boolean to set whether to test for thread count in openmp
logical, parameter :: LHIPREC = .true. ! boolean to set whether to print energies to file in higher precision
integer, parameter :: LWORK = 132  ! size of work arrays for lapack; you may want to optimize this for a given Hamiltonian
integer, parameter :: N_DIM = 2 ! number of dimensions for kk vector; set this to 3 for 3D crystals
integer, parameter :: NUMK = 2000 ! sets number of kk vectors to NUMK along the GKMG path, and to NUMK**2 in the 2D grid
contains

end module settings_mod
