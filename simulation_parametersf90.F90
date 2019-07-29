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
!>  This module will be shared between all simulation programs and
!contins all
!   shared parameters that the program will need.
!
! REVISION HISTORY:
! 01 22 2019 - Initial Version
!------------------------------------------------------------------------------

module simulation_parametersf90
#include <slepc/finclude/slepceps.h>
  use slepceps
  implicit none
  PetscInt,   parameter :: dp = kind(1.d0) 

  ! Initial State
  PetscInt,   parameter :: n_init = 1
  PetscInt,   parameter :: l_init = 0

  ! Blocking Pathways 
  logical,   parameter :: block_pathways = .false.
  logical,   parameter :: block_individual = .false.
  integer,   parameter :: number_to_block = 0

  ! Laser perameters
  PetscReal,  parameter :: electric_field_strength = 0.053375290941998d0 
  PetscReal,  parameter :: omega_electric_field = 0.056953098011833d0

  ! Parameters for time propagation
  PetscReal,  parameter :: time_resolution = 0.05d0
  PetscReal,  parameter :: num_cycles = 10.d0
  PetscBool,  parameter :: custom_envalope_phase = .false.

  ! If above is set to true the below values are used.   
  PetscReal,  parameter :: envelope_phase  = 0.d0
  PetscReal,  parameter :: time_envelope_phase_set = 0.d0

  ! A name for the hdf5 file that the wfns and energy will be saved 2
  character(len = 15), parameter :: hdf5_file_label = 'H'
  ! Pulse envelope function must be either 'sin2' or 'gaussian'
  character(len = 15), parameter :: envelope_function = 'sin2' 

  ! Grid detials
  PetscReal,  parameter :: grid_space = 0.01d0  
  PetscReal,  parameter :: R_max = 150.d0   
  PetscReal,  parameter :: eps_tol = 1d-10
  PetscInt,   parameter :: eps_max_its = 100000

  ! Max energy level and angular momentum
  PetscInt,   parameter :: n_max = 200
  PetscInt,   parameter :: l_max = 120

  ! Nuclear charge 
  ! Coulomb
  PetscReal,  parameter :: sae_c0 = 1.d0 
  ! Yukawa coef
  PetscReal,  parameter :: sae_Zc = 9d0
  ! Yukawa exp
  PetscReal,  parameter :: sae_c =  0.887d0
  ! decay coef
  PetscReal,  parameter :: sae_a1 = -9.9286d0
  ! decay exp
  PetscReal,  parameter :: sae_b1 = 1.3746d0

  ! additional decay 
  PetscReal,  parameter :: sae_a2 = -5.995d0
  PetscReal,  parameter :: sae_b2 = 3.7963d0
  PetscReal,  parameter :: sae_a3 = 0d0
  PetscReal,  parameter :: sae_b3 = 0d0
  PetscReal,  parameter :: sae_a4 = 0d0
  PetscReal,  parameter :: sae_b4 = 0d0
  PetscReal,  parameter :: sae_a5 = 0d0
  PetscReal,  parameter :: sae_b5 = 0d0

  ! Field free matrix 
  PetscBool,  parameter :: eps_ecs_present = .true.
  PetscReal,  parameter :: eps_theta = 0.0d0
  PetscReal,  parameter :: eps_eta   = 0.5d0
  PetscReal,  parameter :: eps_gobbler = 0.9d0
  EPSWhich,   parameter :: eps_which = EPS_SMALLEST_REAL 
  EPSType,    parameter :: eps_type  = EPSKRYLOVSCHUR     
  EPSProblemType, parameter :: eps_problem = EPS_NHEP      
  PetscInt,   parameter :: eps_ncv = PETSC_DEFAULT_INTEGER
  PetscInt,   parameter :: eps_mpd = PETSC_DEFAULT_INTEGER
  PetscBool,  parameter :: eps_two_sided = PETSC_FALSE
  PetscBool,  parameter :: eps_world = PETSC_FALSE
  PetscScalar,parameter :: eps_target = (-0.5d0,0d0)
  MatType,    parameter :: eps_mat_type = MATSBAIJ    
  PetscBool,  parameter :: masking_function_present = .false.
  character(len=20), parameter :: mask_type = 'cosine'
  PetscReal,  parameter :: mask_pow = 1d0
  PetscInt,   parameter :: n_abs = 1000
  PetscInt,   parameter :: l_abs = 600

  ! Some useful constants that show up often
  PetscReal,  parameter :: pi = 3.141592653589793238462643383279502884197169
end module simulation_parametersf90
