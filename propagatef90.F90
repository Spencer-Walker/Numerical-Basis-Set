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
  Vec :: tmp, mask_vector, dipoleA
  PetscReal :: told 
  Mat :: Z_scale,H0_scale,A
  PetscReal :: omega_vector_potential
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
  if (custom_envalope_phase .eqv. .true.) then
    cep = envelope_phase
    tcep = time_envelope_phase_set 
  else 
    tcep = T0/2d0
    cep  = tcep*wa
  end if
  
  if (trim(envelope_function) .eq. 'sin2') then
    E = E0*(-(dcos(cep + (t - tcep)*wa)*dsin((pi*t)/T0)**2d0) - &
    & (2d0*pi*dcos((pi*t)/T0)*dsin((pi*t)/T0)*dsin(cep + (t - tcep)*wa))/(T0*wa))

  else if ( trim(envelope_function) .eq. 'gaussian') then
    E = -E0*dcos(wa*(t-tcep)+cep)*dexp(-dlog(2d0)*((2d0*(t-tcep))/T0)**2d0) -  &
    & (8d0*E0*(t-tcep)*dlog(2d0)/T0**2d0)*dsin(wa*(t-tcep)+cep)*  &
    & dexp(-dlog(2d0)*((2d0*(t-tcep))/T0)**2d0) 
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
  PetscReal       :: t
  Mat             :: J,BB
  integer         :: user
  Vec             :: psi
  PetscScalar     :: E,dotProduct
  PetscErrorCode  :: ierr
  PetscInt        :: step

  print*,t
  call TSGetStepNumber(ts,step,ierr)
  CHKERRA(ierr)
  call MatMult(A,psi,tmp,ierr)
  CHKERRA(ierr)
  call VecDot(psi,tmp,dotProduct,ierr)
  CHKERRA(ierr)
  call VecSetValue(dipoleA,step,dotProduct,INSERT_VALUES,ierr)
  CHKERRA(ierr)
  call MatCopy(H0_scale,J,SUBSET_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)
  call MatAXPY(J,+(E(t)),Z_scale,SUBSET_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)

  told  = t
  
  call VecPointwiseMult(tmp, psi, mask_vector, ierr)
  CHKERRA(ierr)

  call TSSetSolution(ts, tmp, ierr)
  CHKERRA(ierr)

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
  PetscInt          :: i_start,i_end,failures
  PetscInt          :: max_steps
  Vec               :: psi
  PetscViewer       :: viewer
  PetscScalar       :: norm, scale, E, val
  integer           :: nmax,lmax,size,num_grid_points,nabs,labs,index
  integer           :: n,l,ninit,linit
  PetscReal         :: E0,wa,cep,T0,numcycles,dt
  PetscReal         :: h,r0,maxtime,we,mu
  character(len=15) :: label
  external          :: RHSMatrixSchrodinger
  real(dp) :: start_time, end_time

! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------

  call CPU_TIME(start_time)

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    print*,'Unable to initialize PETSc'
    stop
  endif

  ninit = n_init
  linit = l_init
  nmax = n_max
  lmax = l_max

  if (masking_function_present) then
    nabs = n_abs
    labs = l_abs
  end if 

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

  if ( trim(envelope_function) .eq. 'sin2' ) then
    mu = 4d0*dasin(dexp(-0.25d0))**2d0
  else if( trim(envelope_function) .eq. 'gaussian') then
    mu = 8d0*dlog(2d0)/pi**2d0
    print*,'only sin2 is supported'
    stop
  else 
    print*,'only sin2 is are supported'
    stop
  end if 

  omega_vector_potential = 2d0*we/(1d0 + dsqrt(1d0 + mu/num_cycles**2d0))
  wa = omega_vector_potential

  T0 = num_cycles*2d0*pi/wa

  maxtime = T0

  print*, 'E0 =',E0
  print*, 'we =',we
  print*, 'wa =',wa
  print*, 'T0 =',T0
  print*, 'dt =',dt 

  scale =  dcmplx(0.d0,-1.d0)

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
  ! We want the electron to start in the 1s state so we set psi0(0) = 1
  index =  -1 + ninit - (linit*(1 + linit - 2*nmax))/2
  call VecSetValue(psi,index,one,INSERT_VALUES,ierr)
  CHKERRA(ierr)
  call VecAssemblyBegin(psi,ierr)
  CHKERRA(ierr)
  call VecAssemblyEnd(psi,ierr)
  CHKERRA(ierr)
  call VecDuplicate(psi,tmp,ierr)
  CHKERRA(ierr)
  call VecDuplicate(psi,mask_vector,ierr)
  CHKERRA(ierr)

  ! Create the masking vector
  if (masking_function_present) then
    do l = 0,lmax
      do n=l+1,nmax
        index =  -1 + n - (l*(1 + l - 2*nmax))/2
        ! We want the electron to start in the 1s state so we set psi0(0) = 1
        if (l>=labs .and. n>=nabs) then 
          val = abs(dcos((dble(l + (labs - lmax))*pi)/(2.d0*dble(labs))))**0.125d0*  &
          & abs(dcos((dble(n + nabs - nmax)*pi)/(2.d0*dble(nabs))))**0.125d0 
        else if (l>=labs) then
          val = abs(dcos(((l + labs - lmax)*pi)/(2.d0*labs)))**0.125d0
        else if (n>=nabs) then 
          val = abs(dcos(((n + nabs - nmax)*pi)/(2.d0*nabs)))**0.125d0 
        else 
          val = 1d0
        end if 
        call VecSetValue(mask_vector,index,val,INSERT_VALUES,ierr)
        CHKERRA(ierr)
      end do
    end do
    
    call VecAssemblyBegin(mask_vector,ierr) 
    call VecAssemblyEnd(mask_vector,ierr)
  end if 

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
  call TSSetType(ts,TSTHETA,ierr)
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
  call TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP,ierr)
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
  call PCSetType(pc,PCJACOBI,ierr)
  CHKERRA(ierr)
  call KSPSetPC(ksp,pc,ierr)
  CHKERRA(ierr)
  call KSPSetPCSide(ksp,PC_SYMMETRIC,ierr)
  CHKERRA(ierr)
  call SNESSetKSP(snes,ksp,ierr)
  CHKERRA(ierr)
  call TSSetSNES(ts,snes,ierr)
  CHKERRA(ierr)

  
  max_steps = ceiling(T0/dt) 

  call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,max_steps,dipoleA,ierr)
  CHKERRA(ierr)
  call VecSet(dipoleA,0d0*one,ierr)
  CHKERRA(ierr)
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

  call VecAssemblyBegin(dipoleA,ierr)
  call VecAssemblyEnd(dipoleA,ierr)

  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(label)//'_dipoleAcceleration.output',&
  & viewer,ierr)
  CHKERRA(ierr)

  call VecView(dipoleA,viewer,ierr)
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



! --------------------------------------------------------------------------
! Clear memory
! --------------------------------------------------------------------------
  print*,'before clear memory'
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
  call VecDestroy(dipoleA,ierr)
  CHKERRA(ierr)
  call PetscFinalize(ierr)

  call CPU_TIME(end_time)
    
  print*, 'time   :', end_time-start_time


end program Main