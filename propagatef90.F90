!------------------------------------------------------------------------------
! CU Boulder, Jila
!------------------------------------------------------------------------------
!
! PROGRAM:  propagatef90.F90
!
!> Spencer.Walker@colorado.edu
!> Spencer Walker
!
! DESCRIPTION: 
!>  This program propagates the TDSE in time for a hydrogenlike atom in the 
!   presence of a linearly polarized laser in the z direction. 
!
! REVISION HISTORY:
! 22 01 2019 - Initial Version
!------------------------------------------------------------------------------
module time_propagation_module
#include <petsc/finclude/petscts.h>
  use petscts
  implicit none
  Vec :: tmp
  PetscReal :: told 
  Mat :: Z_scale,H0_scale,A
  PetscReal :: omega_vector_potential
  PetscReal, allocatable :: dipoleA(:)
end module time_propagation_module

! --------------------------------------------------------------------------
! Functions
! --------------------------------------------------------------------------
! --------------------------- Evaluate E(t) --------------------------------
function E(t)
  ! This function computes the electric field at some given time 
  use simulation_parametersf90
  use time_propagation_module 
  implicit none 
  PetscReal       :: E0,wa,cep,T0,t,tcep
  PetscScalar     :: E

  E0 = electric_field_strength
  wa  = omega_vector_potential
  T0 = num_cycles*2d0*pi/wa
  if (custom_envalope_phase == .true.) then
    cep = envelope_phase
    tcep = time_envelope_phase_set 
  else 
    tcep = T0/2d0
    cep  = tcep*wa
  end if
  
  if (trim(envelope_function) == 'sin2') then
    E = E0*(-(cos(cep + (t - tcep)*wa)*sin((pi*t)/T0)**2d0) - &
    & (2d0*pi*cos((pi*t)/T0)*sin((pi*t)/T0)*sin(cep + (t - tcep)*wa))/(T0*wa))

  else if ( trim(envelope_function) == 'gaussian') then
    E = -E0*cos(wa*(t-tcep)+cep)*exp(-log(2d0)*((2d0*(t-tcep))/T0)**2d0) -  &
    & (8d0*E0*(t-tcep)*log(2d0)/T0**2d0)*sin(wa*(t-tcep)+cep)*  &
    & exp(-log(2d0)*((2d0*(t-tcep))/T0)**2d0) 
  end if

  return
endfunction E

! --------------------------------------------------------------------------
! Subroutines 
! --------------------------------------------------------------------------
! ----------------------- RHSMatrixSchrodinger -----------------------------
subroutine RHSMatrixSchrodinger(ts,t,psi,J,BB,user,ierr)
  use petscts
  use simulation_parametersf90
  use time_propagation_module
  implicit none
  TS              :: ts
  PetscReal       :: t,real_part
  Mat             :: J,BB
  integer         :: user
  Vec             :: psi
  PetscScalar     :: E,scale
  PetscErrorCode  :: ierr
  PetscInt        :: size1,size2,step

  print*,t
  call TSGetStepNumber(ts,step,ierr)
  CHKERRA(ierr)
  call MatMult(A,psi,tmp,ierr)
  CHKERRA(ierr)
  call VecDotRealPart(psi,tmp,real_part,ierr)
  CHKERRA(ierr)
  dipoleA(step) =  real_part
  call MatCopy(H0_scale,J,SUBSET_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)
  call MatAXPY(J,+(E(t)),Z_scale,SUBSET_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)

  told  = t
  
  return
endsubroutine RHSMatrixSchrodinger

! --------------------------------------------------------------------------
! Main
! --------------------------------------------------------------------------
program Main
use time_propagation_module 
use simulation_parametersf90
#include <petsc/finclude/petscts.h>
  use petscts
  implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  KSP               :: ksp
  PC                :: pc
  SNES              :: snes
  TS                :: ts
  PetscErrorCode    :: ierr
  Mat               :: J
  PetscBool         :: flg
  PetscInt          :: i,i_start,i_end,failures,size1,size2
  PetscInt          :: max_Steps
  Vec               :: psi
  PetscViewer       :: viewer
  PetscMPIInt       :: rank 
  PetscScalar       :: norm, scale, E
  integer           :: nmax,lmax,size,num_grid_points
  PetscReal         :: E0,wa,cep,T0,intensity,wavelength,numcycles,dt
  PetscReal         :: h,r0,maxtime,we,mu,t
  character(len=15) :: label
  external          :: RHSMatrixSchrodinger


! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    print*,'Unable to initialize PETSc'
    stop
  endif

  nmax = n_max
  lmax = l_max
  failures    = -1
  label = hdf5_file_label
  ! If nmax <= lmax+1 we chose the normal basis size but if it is large 
  ! the basis is truncated in l
  if(nmax .le. lmax+1) then
     size = (nmax - 1)*nmax/2 + lmax + 1
  else
     size = (lmax + 1)*(lmax + 2)/2 + (nmax - lmax - 2)*(lmax + 1) + &
          lmax + 1
  endif
  print*,size
  h               = r0/num_grid_points 

  ! These are the default values in atomic units the user is free to change 
  ! them as an extra argument when running the program 
  dt   = time_resolution
  cep  = envelope_phase
  numcycles  = num_cycles
  E0 = electric_field_strength
  we = omega_electric_field

  if ( trim(envelope_function) == 'sin2' ) then
    mu = 4d0*asin(exp(-0.25d0))**2d0
  else if( trim(envelope_function) == 'gaussian') then
    mu = 8d0*log(2d0)/pi**2d0
    print*,'only sin2 is supported'
    stop
  else 
    print*,'only sin2 is are supported'
    stop
  end if 

  omega_vector_potential = 2d0*we/(1d0 + sqrt(1d0 + mu/num_cycles**2d0))
  wa = omega_vector_potential

  T0 = num_cycles*2d0*pi/wa

  maxtime = T0

  print*, 'E0 =',E0
  print*, 'we =',we
  print*, 'wa =',wa
  print*, 'T0 =',T0
  print*, 'dt =',dt 

  scale =  cmplx(0.d0,-1.d0)

! --------------------------------------------------------------------------
! Create Z_scale matrix
! --------------------------------------------------------------------------

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  & trim(label)//'_dipoleMatrix.bin',FILE_MODE_READ,viewer,ierr)
  CHKERRA(ierr)
  call MatCreate(MPI_COMM_WORLD,Z_scale,ierr)
  CHKERRA(ierr)
  call MatSetSizes(Z_scale,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetType(Z_scale,MATMPISBAIJ,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(Z_scale,ierr)
  CHKERRA(ierr)
  call MatSetUp(Z_scale,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(Z_scale,i_start,i_end,ierr)
  CHKERRA(ierr)
  call MatLoad(Z_scale,viewer,ierr)
  CHKERRA(ierr)
  call MatScale(Z_scale,scale,ierr)
  CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr)
  CHKERRA(ierr)

! --------------------------------------------------------------------------
! Create dipole acceleration matrix
! --------------------------------------------------------------------------

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  & trim(label)//'_dipoleAccelerationMatrix.bin',FILE_MODE_READ,viewer,ierr)
  CHKERRA(ierr)
  call MatCreate(MPI_COMM_WORLD,A,ierr)
  CHKERRA(ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetType(A,MATMPISBAIJ,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(A,ierr)
  CHKERRA(ierr)
  call MatSetUp(A,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(A,i_start,i_end,ierr)
  CHKERRA(ierr)
  call MatLoad(A,viewer,ierr)
  CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr)
  CHKERRA(ierr)

! --------------------------------------------------------------------------
! Create the H0_scale matrix
! --------------------------------------------------------------------------

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  & trim(label)//'_fieldFreeMatrix.bin',FILE_MODE_READ,viewer,ierr)
  CHKERRA(ierr)
  call MatCreate(MPI_COMM_WORLD,H0_scale,ierr)
  CHKERRA(ierr)
  call MatSetSizes(H0_scale,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(H0_scale,ierr)
  CHKERRA(ierr)
  call MatSetUp(H0_scale,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(H0_scale,i_start,i_end,ierr)
  CHKERRA(ierr)
  call MatLoad(H0_scale,viewer,ierr)
  CHKERRA(ierr)
  call MatScale(H0_scale,scale,ierr)
  CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr)
  CHKERRA(ierr)

! --------------------------------------------------------------------------
! Create the Jacobian Matrix
! --------------------------------------------------------------------------
  call MatCreate(PETSC_COMM_WORLD,J,ierr)
  CHKERRA(ierr)
  call MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetType(J,MATMPISBAIJ,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(J,ierr)
  CHKERRA(ierr)
  call MatSetUp(J,ierr)
  CHKERRA(ierr)
  call MatCopy(H0_scale,J,DIFFERENT_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)
  told = 0.d0
  call MatAXPY(J,+E(told),Z_scale,DIFFERENT_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)

! --------------------------------------------------------------------------
! Propagate
! --------------------------------------------------------------------------
  ! Create the solution, and initial condition vector
  call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,size,psi,ierr)
  CHKERRA(ierr)
  call VecDuplicate(psi,tmp,ierr)
  CHKERRA(ierr)

  ! We want the electron to start in the 1s state so we set psi0(0) = 1
  call VecSetValue(psi,0,one,INSERT_VALUES,ierr)
  CHKERRA(ierr)
  call VecAssemblyBegin(psi,ierr)
  CHKERRA(ierr)
  call VecAssemblyEnd(psi,ierr)
  CHKERRA(ierr)

  ! Create timestepper context where to compute solutions
  call TSCreate(PETSC_COMM_WORLD,ts,ierr)
  CHKERRA(ierr)
  call TSSetProblemType(ts,TS_LINEAR,ierr)
  CHKERRA(ierr)
  
  ! Tell the timestepper context where to compute solutions
  call TSSetSolution(ts,psi,ierr)
  CHKERRA(ierr)

  ! Tell the timestepper context what type of solver to use. I'll use 
  ! Crank-Nicolson, but this can be changed via a command line option
  call TSSetType(ts,TSCN,ierr)
  CHKERRA(ierr)
  
  ! We will provide the Jacobian matrix only. Petsc will find the RHS. 
  call TSSetRHSFunction(ts,PETSC_NULL_VEC,TSComputeRHSFunctionLinear, &
  & 0, ierr)
  CHKERRA(ierr)

  ! Now we will provide the Jacobian matrix
  call TSSetRHSJacobian(ts,J,J,RHSMatrixSchrodinger,0,ierr)
  CHKERRA(ierr)

  ! I'll set the initial time to be t = 0
  call TSSetTime(ts,zero,ierr)
  CHKERRA(ierr)

  ! Here we set the timestep 
  call TSSetTimeStep(ts,dt,ierr)
  CHKERRA(ierr)

  ! Here we set the maximum time maxtime
  call TSSetMaxTime(ts,maxtime,ierr)
  CHKERRA(ierr)

  ! If maxtime isnt an exact multiple of dt we just interpolate backwards
  ! to get the value at this time
  call TSSetExactFinalTime(ts,TS_EXACTFINALTIME_INTERPOLATE,ierr)
  CHKERRA(ierr)

  ! Here we set the KSP solver in SNES in TS
  call TSGetSNES(ts,snes,ierr)
  CHKERRA(ierr)
  call SNESSetType(snes,SNESKSPONLY,ierr)
  CHKERRA(ierr)
  call SNESGetKSP(snes,ksp,ierr)
  CHKERRA(ierr)
  call KSPSetType(ksp,KSPGMRES,ierr)
  CHKERRA(ierr)
  call KSPGetPC(ksp,pc,ierr)
  CHKERRA(ierr)
  call KSPSetPC(ksp,pc,ierr)
  CHKERRA(ierr)
  call SNESSetKSP(snes,ksp,ierr)
  CHKERRA(ierr)
  call TSSetSNES(ts,snes,ierr)
  CHKERRA(ierr)
  
  call TSGetMaxSteps(ts,max_steps,ierr)
  CHKERRA(ierr)
  allocate(dipoleA(0:int(maxtime/dt)))
  ! Now we finally solve the system 
  call TSSolve(ts,psi,ierr)
  CHKERRA(ierr)

  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(label)//'_psi.output',&
  & viewer,ierr)
  CHKERRA(ierr)
  call VecView(psi,viewer,ierr)
  CHKERRA(ierr)
  call VecAbs(psi,ierr)
  CHKERRA(ierr)
  call VecPointwiseMult(psi,psi,psi,ierr)
  CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr)
  CHKERRA(ierr)

  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(label)//'_rho.output',&
  & viewer,ierr)
  CHKERRA(ierr)
  call VecView(psi,viewer,ierr)
  CHKERRA(ierr)
  call VecSum(psi,norm,ierr)
  CHKERRA(ierr)
  print*,'Norm is',norm


  open(10,file = trim(label)//'_dipoleAcceleration.output')
  write(10,*) dipoleA(:)
  close(10)
  deallocate(dipoleA)
! --------------------------------------------------------------------------
! Clear memory
! --------------------------------------------------------------------------
  call PetscViewerDestroy(viewer,ierr)
  CHKERRA(ierr)
  call MatDestroy(H0_scale,ierr)
  CHKERRA(ierr)
  call MatDestroy(Z_scale,ierr)
  CHKERRA(ierr)
  call MatDestroy(A,ierr)
  CHKERRA(ierr)
  call VecDestroy(psi,ierr)
  CHKERRA(ierr)
  call VecDestroy(tmp,ierr)
  CHKERRA(ierr)
  call PetscFinalize(ierr)
end program Main




