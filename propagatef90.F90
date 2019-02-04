! This program loads the dipole matrix Z, and the field free hamiltonion  
! and then propagates it with Crank-Nicolson method
module dipole
#include <petsc/finclude/petscts.h>
  use petscts
  implicit none
  PetscReal, allocatable :: dipoleA(:)
end module dipole


program main
use dipole 
use simulation_parametersf90
#include <petsc/finclude/petscts.h>
  use petscts
  implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  KSP             :: ksp
  PC              :: pc
  SNES            :: snes
  TS              :: ts
  PetscErrorCode  :: ierr
  Mat             :: user(3),J
  PetscBool       :: flg
  PetscInt        :: i,i_start,i_end,i_Z,i_H0,i_A,failures,size1,size2
  PetscInt        :: max_Steps
  Vec             :: psi,r
  PetscViewer     :: viewer
  PetscMPIInt     :: rank 
  PetscScalar     :: norm
  integer         :: nmax,lmax,size,num_grid_points
  PetscReal       :: E0,w,cep,T0,intensity,wavelength,numcycles,dt
  PetscReal       :: h,r0,max_time
  character(len=15) :: label
  external        :: RHSMatrixSchrodinger
  parameter         (i_Z = 1, i_H0 = 2,i_A = 3)
  
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
  failures        = -1
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
  E0 = Electric_field_strength
  w  = omega
  T0 = period
  max_time = numcycles*T0/2.d0
  print*,'E0 =',E0
  print*,' w =',w
  print*,'T0 =',T0
  print*,'dt =',dt 

  
! --------------------------------------------------------------------------
! Create Z matrix
! --------------------------------------------------------------------------

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  & trim(label)//'_dipoleMatrix.bin',FILE_MODE_READ,viewer,ierr)
  CHKERRA(ierr)
  call MatCreate(MPI_COMM_WORLD,user(i_Z),ierr)
  CHKERRA(ierr)
  call MatSetSizes(user(i_Z),PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(user(i_Z),ierr)
  CHKERRA(ierr)
  call MatSetUp(user(i_Z),ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(user(i_Z),i_start,i_end,ierr)
  CHKERRA(ierr)
  call MatLoad(user(i_Z),viewer,ierr)
  CHKERRA(ierr)

  call MatGetSize(user(i_Z),size1,size2,ierr)
  
! --------------------------------------------------------------------------
! Create dipole acceleration matrix
! --------------------------------------------------------------------------

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  & trim(label)//'_dipoleAccelerationMatrix.bin',FILE_MODE_READ,viewer,ierr)
  CHKERRA(ierr)
  call MatCreate(MPI_COMM_WORLD,user(i_A),ierr)
  CHKERRA(ierr)
  call MatSetSizes(user(i_A),PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(user(i_A),ierr)
  CHKERRA(ierr)
  call MatSetUp(user(i_A),ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(user(i_A),i_start,i_end,ierr)
  CHKERRA(ierr)
  call MatLoad(user(i_A),viewer,ierr)
  CHKERRA(ierr)

  call MatGetSize(user(i_A),size1,size2,ierr)
  


! --------------------------------------------------------------------------
! Create the H0 matrix
! --------------------------------------------------------------------------

  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  & trim(label)//'_fieldFreeMatrix.bin',FILE_MODE_READ,viewer,ierr)
  CHKERRA(ierr)
  call MatCreate(MPI_COMM_WORLD,user(i_H0),ierr)
  CHKERRA(ierr)
  call MatSetSizes(user(i_H0),PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(user(i_H0),ierr)
  CHKERRA(ierr)
  call MatSetUp(user(i_H0),ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(user(i_H0),i_start,i_end,ierr)
  CHKERRA(ierr)
  call MatLoad(user(i_H0),viewer,ierr)
  CHKERRA(ierr)

! --------------------------------------------------------------------------
! Create the Jacobian Matrix
! --------------------------------------------------------------------------
  call MatCreate(PETSC_COMM_WORLD,J,ierr)
  CHKERRA(ierr)
  call MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(J,ierr)
  CHKERRA(ierr)
  call MatSetUp(J,ierr)
  CHKERRA(ierr)

! --------------------------------------------------------------------------
! Propagate
! --------------------------------------------------------------------------
  ! Create the solution, and initial condition vector
  call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,size,psi,ierr)
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
  call TSSetRHSFunction(ts,PETSC_NULL_VEC,TSComputeRHSFunctionLinear,user,&
  & ierr)
  CHKERRA(ierr)

  ! Now we will provide the Jacobian matrix
  call TSSetRHSJacobian(ts,J,J,RHSMatrixSchrodinger,user,ierr)
  CHKERRA(ierr)

  ! I'll set the initial time to be t = 0
  call TSSetTime(ts,zero,ierr)
  CHKERRA(ierr)

  ! Here we set the timestep 
  call TSSetTimeStep(ts,dt,ierr)
  CHKERRA(ierr)

  ! Here we set the maximum time maxtime
  call TSSetMaxTime(ts,max_time,ierr)
  CHKERRA(ierr)

  ! If max_time isnt an exact multiple of dt we just interpolate backwards
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
  allocate(dipoleA(0:int(max_time/dt)))
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
  call MatDestroy(user(i_H0),ierr)
  CHKERRA(ierr)
  call MatDestroy(user(i_Z),ierr)
  CHKERRA(ierr)
  call MatDestroy(user(i_A),ierr)
  CHKERRA(ierr)
  call PetscFinalize(ierr)
end program main


! --------------------------------------------------------------------------
! Functions
! --------------------------------------------------------------------------
! --------------------------- Evaluate E(t) --------------------------------
function E(t)
  ! This function computes the electric field at some given time 
  use simulation_parametersf90
  implicit none 
  PetscReal       :: E0,w,cep,T0,t
  PetscScalar     :: E

  E0 = Electric_field_strength
  w  = omega
  cep = envelope_phase
  T0 = period
  
  E = E0*sin(w*t+cep)*sin(pi*t/T0)**2.d0
  return

endfunction E

! --------------------------------------------------------------------------
! Subroutines 
! --------------------------------------------------------------------------
! ----------------------- RHSMatrixSchrodinger -----------------------------
subroutine RHSMatrixSchrodinger(ts,t,psi,J,BB,user,ierr)
  use petscts
  use simulation_parametersf90
  use dipole
  implicit none
  TS              :: ts
  PetscReal       :: t,T0,real_part
!  PetscReal, allocatable :: dipoleA(:)
  Mat             :: user(3),Z,H0,A,J,BB
  Vec             :: psi,tmp
  PetscScalar     :: E,scale
  PetscErrorCode  :: ierr
  PetscInt        :: i_Z,i_H0,i_A,size1,size2,step
  parameter (i_Z = 1, i_H0 = 2,i_A = 3)
  T0    =  period
  Z     =  user(i_Z)
  H0    =  user(i_H0)
  A     =  user(i_A)
  scale =  cmplx(0.d0,-1.d0)
  print*,t
  call TSGetStepNumber(ts,step,ierr)
  CHKERRA(ierr)
  call VecDuplicate(psi,tmp,ierr)
  CHKERRA(ierr)
  call MatMult(A,psi,tmp,ierr)
  CHKERRA(ierr)
  call VecDotRealPart(psi,tmp,real_part,ierr)
  CHKERRA(ierr)
  dipoleA(step) =  real_part
  call MatCopy(H0,J,DIFFERENT_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)
  call MatAXPY(J,-E(t),Z,DIFFERENT_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)
  call MatScale(J,scale,ierr)
  CHKERRA(ierr)
  
  call VecDestroy(tmp,ierr)
  CHKERRA(ierr)
  return
endsubroutine RHSMatrixSchrodinger


