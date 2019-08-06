! This program generates the field free hamiltonion in energy basis
program main
#include <slepc/finclude/slepceps.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
use hdf5
use slepceps
  implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  PetscInt,   parameter :: dp = kind(1.d0)
  PetscErrorCode      :: ierr
  Mat                 :: H0
  PetscInt            :: i_start, i_end, index, l_i_end, l_i_start
  PetscViewer         :: viewer
  PetscMPIInt         :: proc_id, num_proc, comm 
  PetscInt            :: l, n, nmax, lmax, size, left_index, h5_err
  PetscInt            :: tdse_nmax, tdse_lmax
  PetscScalar, allocatable :: val(:)
  PetscInt,    allocatable :: col(:)
  PetscReal, allocatable   :: El(:,:)
  PetscScalar, allocatable :: E(:,:)
  PetscReal  :: start_time, end_time
  integer(HSIZE_T)    :: dims(1), ener_dims(1:2)
  integer(SIZE_T), parameter :: sdim = 50 
  integer(HID_T)      :: file_id, ener_id, h5_kind, param_file_id
  integer(HID_T)      :: eps_group_id, memtype, eps_dat_id
  integer(HID_T)      :: operators_group_id, operators_dat_id
  integer(HID_T)      :: tdse_group_id, tdse_dat_id
  character(len = 15) :: label ! File name without .h5 extension
  character(len = 3)  :: strl  ! file number (l converted to a string)
  character(len = 6)  :: fmt   ! format descriptor
  character(len = 30) :: file_name
  character(len = 12) :: ener_name
  character(len = 50) :: tmp_character
  MatType :: mat_type
! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------
  call CPU_TIME(start_time)

  10 format(A2,I4)

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    print*,'Unable to initialize PETSc'
    stop
  endif

  ! Opens hdf5 to read in the label.h5 file 
  call h5open_f( h5_err)
  if ( h5_err /= 0 ) then
    print*,'h5open_f failed'
    stop
  end if

  comm  = MPI_COMM_WORLD

  call h5fopen_f( trim("parameters.h5"), H5F_ACC_RDONLY_F, param_file_id, h5_err)
  call h5gopen_f(param_file_id, "EPS", eps_group_id, h5_err)

  dims(1) = 1
  CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, h5_err)
  CALL H5Tset_size_f(memtype, sdim, h5_err)

  call h5dopen_f(eps_group_id, "label", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, memtype, label, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "l_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, lmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "n_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, nmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5gopen_f(param_file_id, "TDSE", tdse_group_id, h5_err)

  call h5dopen_f(tdse_group_id, "l_max", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_INTEGER, tdse_lmax, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call h5dopen_f(tdse_group_id, "n_max", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_INTEGER, tdse_nmax, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call h5gopen_f(param_file_id, "operators", operators_group_id, h5_err)
  call h5dopen_f(operators_group_id, "mat_type", operators_dat_id, h5_err)
  call h5dread_f(operators_dat_id, memtype, tmp_character,dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  if (trim(tmp_character) .eq. "MATSBAIJ") then
    mat_type = MATSBAIJ
  else if (trim(tmp_character) .eq. "MATAIJ") then
    mat_type = MATAIJ
  else 
    print*, "mat_type not supported defaulting to MATAIJ"
    mat_type = MATAIJ   
  end if 

  call MPI_Comm_rank(comm,proc_id,ierr)
  CHKERRA(ierr)
  call MPI_Comm_size(comm,num_proc,ierr)
  CHKERRA(ierr)

  allocate(E(nmax,0:tdse_lmax))
  allocate(El(nmax,2))
  
  ! If nmax <= lmax+1 then we use the normal rule but if its larger we 
  ! need to truncate our basis
  if(tdse_nmax .le. tdse_lmax+1) then
    size = (tdse_nmax - 1)*tdse_nmax/2 + tdse_lmax + 1
  else
    size = (tdse_lmax + 1)*(tdse_lmax + 2)/2 + (tdse_nmax - tdse_lmax - 2)*(tdse_lmax + 1) +&
    tdse_lmax + 1
  endif 

  ! Adds the .h5 extension to the input file 
  file_name = trim(label)//'.h5'

  ! Opens the above file
  call h5fopen_f( trim(file_name), H5F_ACC_RDWR_F, file_id, h5_err)

  ! Converts fortrans double kind to an hdf5 double type 
  h5_kind = h5kind_to_type( dp, H5_REAL_KIND)
! --------------------------------------------------------------------------
! Create H0 matrix
! --------------------------------------------------------------------------
  call MatCreate(comm,H0,ierr)
  CHKERRA(ierr)
  call MatSetSizes(H0,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetType(H0,mat_type,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(H0,ierr)
  CHKERRA(ierr)
  call MatSetUp(H0,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(H0,i_start,i_end,ierr)
  CHKERRA(ierr)

  allocate(col(tdse_nmax),val(tdse_nmax))


  l_i_start = -1
  l_i_end = -1
  do l = 0,tdse_lmax
    do n=l+1,tdse_nmax
      index = -1 + n - (l*(1 + l - 2*tdse_nmax))/2
      if (index>=i_start .and. l_i_start == -1) then
        l_i_start = l
      end if
      if (index<=i_end) then
        l_i_end = l
      end if 
    end do
  end do


  ! Itterate over angular momentum to calculate the matrix elements
  do l = l_i_start,l_i_end  
    print 10, 'l0', l
    ! Figure out how many digits are needed for the input file
    if ( l .le. 9 ) then
      fmt = '(I1.1)'
    elseif (l .le. 99) then
      fmt = '(I2.2)'
    else
      fmt = '(I3.3)'
    endif
    ! l converted to a string
    write(strl,fmt) l

    ! Energy dataset name 
    ener_name ='Energy_l'//trim(strl)

    ! Opens the Energy_l#l dataset
    call h5dopen_f(file_id, ener_name, ener_id, h5_err)
    
    ! Sets the dimensions of the Energy data in the dataset
    ener_dims(1) = nmax-l
    ener_dims(2) = 2
    
    ! Fills the lth column of E(:,:) with the energy data from the dataset
    call h5dread_f( ener_id, h5_kind, El(1:nmax-l,:), ener_dims, h5_err) 
    E(1:nmax-l,l) = dcmplx(El(1:nmax-l,1),El(1:nmax-l,2))
    ! Closes the dataset
    call h5dclose_f( ener_id, h5_err)
    ! For each l value we compute the corresponding matrix elements of H0
    
    do n=l+1,tdse_nmax    
      ! Here I convert the n and l value to its corresponding left_index 
      left_index =  -1 + n - (l*(1 + l - 2*tdse_nmax))/2
      ! Convert the real energy to a PetscScalar
      val(1) = E(n-l,l)  
      ! Insert the energy into the field free matrix 
      call MatSetValue(H0,left_index,left_index,val(1),INSERT_VALUES,ierr)
      CHKERRA(ierr)
    enddo
  enddo

  deallocate(col,val)

  ! Build the matrix
  call MatAssemblyBegin(H0,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)
  call MatAssemblyEnd(H0,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)

  ! Write the matrix to a binary file
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  trim(label)//"_fieldFreeMatrix.bin",FILE_MODE_WRITE,viewer,ierr)
  CHKERRA(ierr)
  call MatView(H0,viewer,ierr)
  CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr)
  CHKERRA(ierr)

  ! Clear memory 
  call MatDestroy(H0,ierr)
  CHKERRA(ierr)
  deallocate(E,El)
  call h5fclose_f( file_id, h5_err)
  call h5gclose_f( eps_group_id, h5_err)
  call h5gclose_f( operators_group_id, h5_err)
  call h5gclose_f( tdse_group_id, h5_err)
  call h5fclose_f( param_file_id, h5_err)
  call h5close_f( h5_err)

  call PetscFinalize(ierr)
  call CPU_TIME(end_time)
  
  if(proc_id .eq. 0) then
    print 20, 'time   :', end_time-start_time
    20 format(A8,ES9.2)
  end if 

end program main