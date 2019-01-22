! This program generates the field free hamiltonion in energy basis
program main
use hdf5
#include <slepc/finclude/slepceps.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

  use slepceps
  implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  integer,  parameter  :: dp = kind(0.d0) ! double precision double precision
  PetscErrorCode      :: ierr
  Mat                 :: H0
  PetscInt            :: i_start, i_end
  integer(HID_T)      :: file_id, ener_id, h5_kind
  PetscViewer         :: viewer, hdf5v
  PetscMPIInt         :: proc_id, num_proc, comm 
  PetscInt            :: one, i
  integer             :: l, n, nmax, lmax, size, index, h5_err
  character(len = 15) :: label ! File name without .h5 extension
  character(len = 3)  :: strl  ! file number (l converted to a string)
  character(len = 6)  :: fmt   ! format descriptor
  character(len = 30) :: file_name
  character(len = 12) :: ener_name
  character(len = 25) :: data
  PetscScalar         :: val(1)
  integer(HSIZE_T)    :: ener_dims(1:1)
  PetscReal, allocatable   :: E(:,:)
! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------
  comm             = MPI_COMM_WORLD
  lmax             = 10
  nmax             = 20
  one              = 1
  label =  'H_test'

  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  if ( ierr /= 0 ) then
    print*,'SlepcInitialize failed'
    stop
  endif

  call MPI_Comm_rank(comm,proc_id,ierr)
  CHKERRA(ierr)
  call MPI_Comm_size(comm,num_proc,ierr)
  CHKERRA(ierr)

  allocate(E(nmax,0:lmax))
  
  ! If nmax <= lmax+1 then we use the normal rule but if its larger we 
  ! need to truncate our basis
  if(nmax .le. lmax+1) then
    size = (nmax - 1)*nmax/2 + lmax + 1
  else
    size = (lmax + 1)*(lmax + 2)/2 + (nmax - lmax - 2)*(lmax + 1) +	&
    &	lmax + 1
  endif 

  ! Opens hdf5 to read in the label.h5 file 
  call h5open_f( h5_err)
  if ( h5_err /= 0 ) then
    print*,'h5open_f failed'
    stop
  end if

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
  call MatSetFromOptions(H0,ierr)
  CHKERRA(ierr)
  call MatSetUp(H0,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(H0,i_start,i_end,ierr)
  CHKERRA(ierr)

  ! Itterate over angular momentum to calculate the matrix elements
  do l = 0,lmax    
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
    
    ! Fills the lth column of E(:,:) with the energy data from the dataset
    call h5dread_f( ener_id, h5_kind, E(1:nmax-l,l), ener_dims, h5_err) 

    ! Closes the dataset
    call h5dclose_f( ener_id, h5_err)
    
    ! For each l value we compute the corresponding matrix elements of H0
    do n=l+1,nmax
      
      ! Here I convert the n and l value to its corresponding index 
      if(n .le. lmax+1) then
          index = (n - 1)*n/2 + l
      else
          index = (lmax + 1)*(lmax + 2)/2 + (n - lmax - 2)*(lmax + 1) &
          &	+ l
      endif    
      
      ! Convert the real energy to a PetscScalar
      val = E(n-l,l)
      ! Insert the energy into the field free matrix 
      call MatSetValue(H0,index,index,val,INSERT_VALUES,ierr)
      CHKERRA(ierr)
    enddo
    close(20)
  enddo

  ! Build the matrix
  call MatAssemblyBegin(H0,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)
  call MatAssemblyEnd(H0,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)

  ! Write the matrix to a binary file
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  &	label//"_fieldFreeMatrix.bin",FILE_MODE_WRITE,viewer,ierr)
  CHKERRA(ierr)
  call MatView(H0,viewer,ierr)
  CHKERRA(ierr)

  ! Clear memory 
  call MatDestroy(H0,ierr)
  CHKERRA(ierr)
  call h5fclose_f( file_id, h5_err)
  call h5close_f( h5_err)
  call SlepcFinalize(ierr)
end program main
