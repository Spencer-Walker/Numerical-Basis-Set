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
  integer,   parameter :: dp = kind(1.d0) 

  ! Grid detials
  real(dp),  parameter :: grid_space = 0.01d0  
  real(dp),  parameter :: R_min = 30.d0    
  real(dp),  parameter :: R_max = 100.d0   

  ! Max energy level and angular momentum
  integer,   parameter :: n_max = 150
  integer,   parameter :: l_max = 50
  
  ! Tolerances needed for searching for energies 
  real(dp),  parameter :: binary_search_tol = 1d-4  
  real(dp),  parameter :: refinement_tol = 1d-12
  real(dp),  parameter :: E_upper  = 100.d0    
  real(dp),  parameter :: E_lower  = -1d0  
  
  ! Nuclear charge 
  ! Coulomb
  real(dp),  parameter :: sae_c0 = 1.d0 
  ! Yukawa coef
  real(dp),  parameter :: sae_Zc = 0d0
  ! Yukawa exp
  real(dp),  parameter :: sae_c = 0d0
  ! decay coef
  real(dp),  parameter :: sae_a1 = 0d0
  ! decay exp
  real(dp),  parameter :: sae_b1 = 0d0

  ! additional decay 
  real(dp),  parameter :: sae_a2 = 0d0
  real(dp),  parameter :: sae_b2 = 0d0
  real(dp),  parameter :: sae_a3 = 0d0
  real(dp),  parameter :: sae_b3 = 0d0
  real(dp),  parameter :: sae_a4 = 0d0
  real(dp),  parameter :: sae_b4 = 0d0
  real(dp),  parameter :: sae_a5 = 0d0
  real(dp),  parameter :: sae_b5 = 0d0

  ! Parameters for time propagation
  PetscReal, parameter :: time_resolution = 0.05d0
  PetscReal, parameter :: num_cycles = 2.d0
  logical,   parameter :: custom_envalope_phase = .false.
  ! If above is set to true the below values are used.   
  PetscReal, parameter :: envelope_phase  = 0.d0
  PetscReal, parameter :: time_envelope_phase_set = 0.d0

  ! Field free matrix 
  logical,  parameter :: ecs_present = .true.
  real(dp), parameter :: gobbler = 0.7d0
  logical,  parameter :: masking_function_present = .true.
  integer,  parameter :: n_abs = 75
  integer,  parameter :: l_abs = 51
  logical,  parameter :: energy_absorber_present = .false.    
  real(dp), parameter :: energy_absorber = 100d0
  real(dp), parameter :: V_max = 2.d0

  ! Laster perameters
  PetscReal, parameter :: electric_field_strength =  0.053375290941998d0
  PetscReal, parameter :: omega_electric_field = 0.056953098011833d0

  ! A name for the hdf5 file that the wfns and energy will be saved 2
  character(len = 15), parameter :: hdf5_file_label = 'H'
  ! Pulse envelope function must be either 'sin2' or 'gaussian'
  character(len = 15), parameter :: envelope_function = 'sin2' 

  ! Some useful constants that show up often
  PetscScalar, parameter :: two = 2.d0
  PetscScalar, parameter :: one = 1.d0
  PetscReal, parameter   :: zero = 0.d0
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197169
  end module simulation_parametersf90
