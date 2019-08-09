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
  PetscScalar, allocatable :: E(:)
  PetscInt  :: masking_function_present
end module time_propagation_module

! --------------------------------------------------------------------------
! Subroutines 
! --------------------------------------------------------------------------
! ----------------------- RHSMatrixSchrodinger -----------------------------
subroutine RHSMatrixSchrodinger(ts,t,psi,J,BB,user,ierr)
  use petscts
  use time_propagation_module
  implicit none
  TS              :: ts
  PetscReal       :: t,dt
  Mat             :: J,BB
  PetscInt        :: user
  Vec             :: psi
  PetscScalar     :: dotProduct
  PetscErrorCode  :: ierr
  PetscInt        :: step
  character(9)   :: tstring
  
  call TSGetStepNumber(ts,step,ierr)
  CHKERRA(ierr)

  dt = t-told

  if (mod(floor(t/dt),100) .eq. 0) then
    write(tstring, "(ES9.2)") t
    call PetscPrintf(MPI_COMM_WORLD, "t = "//tstring//"\n", ierr)
    CHKERRA(ierr)
  end if 

  call MatCopy(H0_scale,J,SUBSET_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)

  call MatMult(A,psi,tmp,ierr)
  CHKERRA(ierr)
  call VecDot(psi,tmp,dotProduct,ierr)
  CHKERRA(ierr)
  call VecSetValue(dipoleA,step,dotProduct,INSERT_VALUES,ierr)
  CHKERRA(ierr)
  call MatCopy(H0_scale,J,SUBSET_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)
  call MatAXPY(J,E(step),Z_scale,SUBSET_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)

  if (masking_function_present .eq. 1) then
    call VecPointwiseMult(tmp, psi, mask_vector, ierr)
    CHKERRA(ierr)
    call TSSetSolution(ts, tmp, ierr)
    CHKERRA(ierr)
  end if 

  told  = t

  return
endsubroutine RHSMatrixSchrodinger

! --------------------------------------------------------------------------
! Main
! --------------------------------------------------------------------------
program Main
use time_propagation_module 
use petscts
use hdf5
#include <petsc/finclude/petscts.h>
  implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  KSP                 :: ksp
  PC                  :: pc
  SNES                :: snes
  TS                  :: ts
  PetscErrorCode      :: ierr
  Mat                 :: J
  PetscInt            :: i_start,i_end, i
  PetscInt            :: max_steps, proc_id, num_proc
  Vec                 :: psi
  PetscViewer         :: viewer
  PetscScalar         :: norm, val
  PetscInt            :: nmax,lmax,size,nabs,labs,index
  PetscInt            :: n,l, h5_err, mask_pow, num_states
  PetscInt            :: operators_local
  integer(HSIZE_T)    :: dims(1)
  PetscReal           :: dt
  PetscReal           :: maxtime
  character(len=15)   :: mask
  character(len = 15) :: label ! File name without .h5 extension
  external            :: RHSMatrixSchrodinger
  PetscReal           :: start_time, end_time
  integer(HID_T)      :: param_file_id
  integer(HID_T)      :: eps_group_id, memtype, eps_dat_id
  integer(HID_T)      :: operators_group_id, operators_dat_id
  integer(HID_T)      :: start_group_id, start_dat_id
  integer(HID_T)      :: laser_group_id, laser_dat_id
  integer(HID_T)      :: tdse_group_id, tdse_dat_id
  character(len = 300) :: tmp_character, operator_directory
  MatType             :: mat_type
  integer(SIZE_T), parameter :: sdim = 300 
  PetscReal,  parameter :: pi = 3.141592653589793238462643383279502884197169
  PetscReal,  allocatable :: EE(:)
  PetscInt,   allocatable  :: init_n(:), init_l(:)
  PetscReal,  allocatable :: init_amp(:), init_phase(:)
! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------
  call CPU_TIME(start_time)

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    call PetscPrintf(MPI_COMM_WORLD, 'Unable to initialize PETSC\n', ierr)
    CHKERRA(ierr)
    stop
  endif

  ! Opens hdf5 to read in the label.h5 file 
  call h5open_f( h5_err)
  if ( h5_err /= 0 ) then
    call PetscPrintf(MPI_COMM_WORLD, 'h5open_f failed\n', ierr)
    CHKERRA(ierr)
    stop
  end if

  call MPI_Comm_rank(MPI_COMM_WORLD,proc_id,ierr)
  CHKERRA(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,num_proc,ierr)
  CHKERRA(ierr)

  dims(1) = 1
  call h5fopen_f( trim("parameters.h5"), H5F_ACC_RDONLY_F, param_file_id, h5_err)
  call h5gopen_f(param_file_id, "TDSE", tdse_group_id, h5_err)

  call h5dopen_f(tdse_group_id, "l_max", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_INTEGER, lmax, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call h5dopen_f(tdse_group_id, "n_max", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_INTEGER, nmax, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call H5Tcopy_f(H5T_FORTRAN_S1, memtype, h5_err)
  call H5Tset_size_f(memtype, sdim, h5_err)
  
  call h5dopen_f(tdse_group_id, "mask_type", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, memtype, mask, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call h5dopen_f(tdse_group_id, "mask_present", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_INTEGER, masking_function_present, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call h5dopen_f(tdse_group_id, "mask_n_abs", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_INTEGER, nabs, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call h5dopen_f(tdse_group_id, "mask_l_abs", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_INTEGER, labs, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call h5dopen_f(tdse_group_id, "mask_pow", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_DOUBLE, mask_pow, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call h5gopen_f(param_file_id, "EPS", eps_group_id, h5_err)

  call h5dopen_f(eps_group_id, "label", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, memtype, label, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5gopen_f(param_file_id, "start_state", start_group_id, h5_err)
  dims(1) = 1
  call h5dopen_f(start_group_id, "num_states", start_dat_id, h5_err)
  call h5dread_f(start_dat_id, H5T_NATIVE_INTEGER, num_states, dims, h5_err)
  call h5dclose_f(start_dat_id, h5_err)

  allocate(init_n(num_states))
  allocate(init_l(num_states))
  allocate(init_amp(num_states))
  allocate(init_phase(num_states))
  
  dims(1) = num_states

  call h5dopen_f(start_group_id, "n_index", start_dat_id, h5_err)
  call h5dread_f(start_dat_id, H5T_NATIVE_INTEGER, init_n, dims, h5_err)
  call h5dclose_f(start_dat_id, h5_err)

  call h5dopen_f(start_group_id, "l_index", start_dat_id, h5_err)
  call h5dread_f(start_dat_id, H5T_NATIVE_INTEGER, init_l, dims, h5_err)
  call h5dclose_f(start_dat_id, h5_err)

  call h5dopen_f(start_group_id, "amplitude", start_dat_id, h5_err)
  call h5dread_f(start_dat_id, H5T_NATIVE_DOUBLE, init_amp, dims, h5_err)
  call h5dclose_f(start_dat_id, h5_err)

  call h5dopen_f(start_group_id, "phase", start_dat_id, h5_err)
  call h5dread_f(start_dat_id, H5T_NATIVE_DOUBLE, init_phase, dims, h5_err)
  call h5dclose_f(start_dat_id, h5_err)

  ! If nmax <= lmax+1 we chose the normal basis size but if it is large 
  ! the basis is truncated in l
  if(nmax .le. lmax+1) then
      size = (nmax - 1)*nmax/2 + lmax + 1
  else
      size = (lmax + 1)*(lmax + 2)/2 + (nmax - lmax - 2)*(lmax + 1) + &
          lmax + 1
  endif

  ! These are the default values in atomic units the user is free to change 
  ! them as an extra argument when running the program 
  dims(1) = 1
  call h5dopen_f(tdse_group_id, "delta_t", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_DOUBLE, dt, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  write(tmp_character,"(ES9.2)") dt
  call PetscPrintf(MPI_COMM_WORLD, 'dt ='//trim(tmp_character)//"\n", ierr)
  CHKERRA(ierr)

  call h5gopen_f(param_file_id, "operators", operators_group_id, h5_err)

  call h5dopen_f(operators_group_id, "mat_type", operators_dat_id, h5_err)
  call h5dread_f(operators_dat_id, memtype, tmp_character,dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  if (trim(tmp_character) .eq. "MATSBAIJ") then
    mat_type = MATSBAIJ
  else if (trim(tmp_character) .eq. "MATAIJ") then
    mat_type = MATAIJ
  else 
    call PetscPrintf(MPI_COMM_WORLD,  "mat_type not supported defaulting to MATAIJ\n", ierr)
    CHKERRA(ierr)
    mat_type = MATAIJ   
  end if 

  call h5dopen_f(operators_group_id, "location", operators_dat_id, h5_err)
  call h5dread_f(operators_dat_id, memtype, operator_directory, dims, h5_err)
  call h5dclose_f(operators_dat_id, h5_err)

  call h5dopen_f(operators_group_id, "local", operators_dat_id, h5_err)
  call h5dread_f(operators_dat_id, H5T_NATIVE_INTEGER, operators_local, dims, h5_err)
  call h5dclose_f(operators_dat_id, h5_err)

! --------------------------------------------------------------------------
! Create Z_scale matrix
! --------------------------------------------------------------------------
  if ( operators_local .eq. 0) then
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
    & trim(operator_directory)//'/'//&
    & trim(label)//'_dipoleMatrix.bin',FILE_MODE_READ,viewer,ierr)
    CHKERRA(ierr)
  else 
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
    & trim(label)//'_dipoleMatrix.bin',FILE_MODE_READ,viewer,ierr)
    CHKERRA(ierr)
  end if 

  call MatCreate(MPI_COMM_WORLD,Z_scale,ierr)
  CHKERRA(ierr)
  call MatSetSizes(Z_scale,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetType(Z_scale,mat_type,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(Z_scale,ierr)
  CHKERRA(ierr)
  call MatSetUp(Z_scale,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(Z_scale,i_start,i_end,ierr)
  CHKERRA(ierr)
  call MatLoad(Z_scale,viewer,ierr)
  CHKERRA(ierr)
  call MatScale(Z_scale,(0.d0,-1.d0),ierr)
  CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr)
  CHKERRA(ierr)

! --------------------------------------------------------------------------
! Create dipole acceleration matrix
! --------------------------------------------------------------------------
  if ( operators_local .eq. 0) then
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
    & trim(operator_directory)//'/'//&
      trim(label)//'_dipoleAccelerationMatrix.bin',FILE_MODE_READ,viewer,ierr)
    CHKERRA(ierr)
  else 
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
      trim(label)//'_dipoleAccelerationMatrix.bin',FILE_MODE_READ,viewer,ierr)
    CHKERRA(ierr)
  end if 
  
  call MatCreate(MPI_COMM_WORLD,A,ierr)
  CHKERRA(ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetType(A,mat_type,ierr)
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
  if ( operators_local .eq. 0) then
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
    & trim(operator_directory)//'/'//&
    & trim(label)//'_fieldFreeMatrix.bin',FILE_MODE_READ,viewer,ierr)
    CHKERRA(ierr)
  else 
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
    & trim(label)//'_fieldFreeMatrix.bin',FILE_MODE_READ,viewer,ierr)
    CHKERRA(ierr)
  end if 

  call MatCreate(MPI_COMM_WORLD,H0_scale,ierr)
  CHKERRA(ierr)
  call MatSetSizes(H0_scale,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetType(H0_scale,mat_type,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(H0_scale,ierr)
  CHKERRA(ierr)
  call MatSetUp(H0_scale,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(H0_scale,i_start,i_end,ierr)
  CHKERRA(ierr)
  call MatLoad(H0_scale,viewer,ierr)
  CHKERRA(ierr)
  call MatScale(H0_scale,(0.d0,-1.d0),ierr)
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
  call MatSetType(J,mat_type,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(J,ierr)
  CHKERRA(ierr)
  call MatSetUp(J,ierr)
  CHKERRA(ierr)
  call MatCopy(H0_scale,J,DIFFERENT_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)
  told = 0.d0
  call MatAXPY(J,(0.0d0,0.0d0),Z_scale,DIFFERENT_NONZERO_PATTERN,ierr)
  CHKERRA(ierr)

! --------------------------------------------------------------------------
! Propagate
! --------------------------------------------------------------------------
  ! Create the solution, and initial condition vector
  call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,size,psi,ierr)
  CHKERRA(ierr)
  ! We want the electron to start in the 1s state so we set psi0(0) = 1
  if (proc_id .eq. 0) then
    do i = 1, num_states 
      index =  -1 + init_n(i) - (init_l(i)*(1 + init_l(i) - 2*nmax))/2
      call VecSetValue(psi,index,init_amp(i)*zexp(dcmplx(0d0,init_phase(i))),ADD_VALUES,ierr)
      CHKERRA(ierr)
    end do 
  end if
  call VecAssemblyBegin(psi,ierr)
  CHKERRA(ierr)
  call VecAssemblyEnd(psi,ierr)
  CHKERRA(ierr)
  call VecDuplicate(psi,tmp,ierr)
  CHKERRA(ierr)
  call VecDuplicate(psi,mask_vector,ierr)
  CHKERRA(ierr)

  ! Create the masking vector
  if (masking_function_present .eq. 1) then
    do l = 0,lmax
      do n=l+1,nmax
        index =  -1 + n - (l*(1 + l - 2*nmax))/2
        if (trim(mask) == 'cosine') then
          ! We want the electron to start in the 1s state so we set psi0(0) = 1
          if (l>=labs .and. n>=nabs) then 
            val = (dabs(dcos(((dble(l - labs))/(dble(lmax-labs)))*pi/2d0))**(1d0/mask_pow)*  &
            & dabs(dcos(((dble(n - nabs)*pi)/(dble(nmax-nabs)))*pi/2d0))**(1d0/mask_pow))**(dt/maxtime)
          else if (l>=labs) then
            val = (dabs(dcos(((dble(l -labs))/(dble(lmax-labs)))*pi/2d0))**(1d0/mask_pow))**(dt/maxtime)
          else if (n>=nabs) then 
            val = (dabs(dcos(((dble(n - nabs ))/(dble(nmax-nabs)))*pi/2d0))**(1d0/mask_pow))**(dt/maxtime)
          else 
            val = 1d0
          end if 
        else if (trim(mask) == 'raised cosine') then
          val = (0.5d0*(1d0 + dcos(pi*dble(n)/dble(nmax))))**(dt/maxtime)
        else 
          val = 1.d0
        end if
        call VecSetValue(mask_vector,index,val,INSERT_VALUES,ierr)
        CHKERRA(ierr)
      end do
    end do
    call VecAssemblyBegin(mask_vector,ierr) 
    CHKERRA(ierr)
    call VecAssemblyEnd(mask_vector,ierr)
    CHKERRA(ierr)
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
  & proc_id, ierr)
  CHKERRA(ierr)
  ! Now we will provide the Jacobian matrix
  call TSSetRHSJacobian(ts,J,J,RHSMatrixSchrodinger,0,ierr)
  CHKERRA(ierr)

  call h5gopen_f(param_file_id, "laser", laser_group_id, h5_err)
  
  call h5dopen_f(laser_group_id, "E_length", laser_dat_id, h5_err)
  call h5dread_f(laser_dat_id, H5T_NATIVE_INTEGER, max_steps, dims, h5_err)
  call h5dclose_f( laser_dat_id, h5_err)
  
  allocate(E(0:max_steps-1))
  allocate(EE(0:max_steps-1))
  
  dims(1) = max_steps
  call h5dopen_f(laser_group_id, "E", laser_dat_id, h5_err)
  call h5dread_f(laser_dat_id, H5T_NATIVE_DOUBLE, EE, dims, h5_err)
  call h5dclose_f( laser_dat_id, h5_err)

  E = dcmplx(EE,0d0)

  ! I'll set the initial time to be t = 0
  call TSSetTime(ts,0d0,ierr)
  CHKERRA(ierr)
  ! Here we set the timestep 
  call TSSetTimeStep(ts,dt,ierr)
  CHKERRA(ierr)

  maxtime = max_steps*dt

  write(tmp_character,"(ES9.2)") maxtime

  call PetscPrintf(MPI_COMM_WORLD, 'maxtime ='//trim(tmp_character)//"\n", ierr)
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

  call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,max_steps,dipoleA,ierr)
  CHKERRA(ierr)
  call VecSet(dipoleA,(0d0,0d0),ierr)
  CHKERRA(ierr)
  ! Now we finally solve the system 
  call TSSolve(ts,psi,ierr)
  CHKERRA(ierr)

  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(label)//'_psi.output',&
  & viewer,ierr)
  CHKERRA(ierr)
  call PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_COMMON,ierr)
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
  CHKERRA(ierr)
  call VecAssemblyEnd(dipoleA,ierr)
  CHKERRA(ierr)
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(label)//'_dipoleAcceleration.output',&
  & viewer,ierr)
  CHKERRA(ierr)
  call PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_COMMON,ierr)
  CHKERRA(ierr)
  call VecView(dipoleA,viewer,ierr)
  CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr)
  CHKERRA(ierr)
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD,trim(label)//'_rho.output',&
  & viewer,ierr)
  CHKERRA(ierr)
  call PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_COMMON,ierr)
  CHKERRA(ierr)
  call VecView(psi,viewer,ierr)
  CHKERRA(ierr)
  call VecSum(psi,norm,ierr)
  CHKERRA(ierr)

  write(tmp_character, "(ES9.2)")  real(real(norm))
  call PetscPrintf(MPI_COMM_WORLD, 'Norm ='//trim(tmp_character)//"\n", ierr)
  CHKERRA(ierr)
  
! --------------------------------------------------------------------------
! Clear memory
! --------------------------------------------------------------------------
  deallocate(E)
  deallocate(EE)
  deallocate(init_n)
  deallocate(init_l)
  deallocate(init_amp)
  deallocate(init_phase)

  call h5gclose_f( eps_group_id, h5_err)
  call h5gclose_f( operators_group_id, h5_err)
  call h5gclose_f( tdse_group_id, h5_err)
  call h5fclose_f( param_file_id, h5_err)
  call h5close_f( h5_err)

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
  call CPU_TIME(end_time)
  write(tmp_character, "(ES9.2)")  end_time-start_time
  call PetscPrintf(MPI_COMM_WORLD, 'time   :'//trim(tmp_character)//"\n", ierr)
  CHKERRA(ierr)
  write(tmp_character,"(ES9.2)") maxtime
  call PetscPrintf(MPI_COMM_WORLD, 'maxtime ='//trim(tmp_character)//"\n", ierr)
  CHKERRA(ierr)
  call PetscFinalize(ierr)
  CHKERRA(ierr)

end program Main