!------------------------------------------------------------------------------
! CU Boulder, Jila
!------------------------------------------------------------------------------
!
! MODULE:  simulation_parametersf90.F90
!
!> Spencer.Walker@colorado.edu
!> Spencer Walker
!
! DESCRIPTION: 
!>  This module will be shared between all simulation programs and contins all
!   shared parameters that the program will need.
!
! REVISION HISTORY:
! 01 22 2019 - Initial Version
!------------------------------------------------------------------------------

module simulation_parametersf90
#include <slepc/finclude/slepceps.h>
  use slepceps
  implicit none
  integer,   parameter :: dp = kind(0.d0) 
  real(dp),  parameter :: grid_space = 0.065d0  
  real(dp),  parameter :: R_min = 30.d0    
  real(dp),  parameter :: R_max = 1000.d0   
  integer,   parameter :: n_max = 500
  integer,   parameter :: l_max = 50
  real(dp),  parameter :: binary_search_tol = 1d-5  
  real(dp),  parameter :: refinement_tol = 1d-14  
  real(dp),  parameter :: E_max  = 1000.d0    
  real(dp),  parameter :: E_min  = -1.d0      
  integer,   parameter :: Z_nuc  = 1 

  PetscReal, parameter :: time_resolution = 0.05d0
  PetscReal, parameter :: envelope_phase  = 0.d0
  PetscReal, parameter :: num_cycles = 2.d0
  PetscReal, parameter :: Electric_field_strength = 0.053375290941998d0

  PetscReal, parameter :: omega =   0.056953098011833d0
  PetscReal, parameter :: max_time = 2.206441976474781d+02

  character(len = 15), parameter :: hdf5_file_label = 'H' 
  PetscScalar, parameter :: two = 2.d0
  PetscScalar, parameter :: one = 1.d0
  PetscReal, parameter   :: zero = 0.d0
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197169
end module simulation_parametersf90
