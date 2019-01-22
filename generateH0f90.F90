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
  PetscViewer         :: viewer
  PetscMPIInt         :: proc_id, num_proc 
  PetscInt            :: one, size1, size2, i
	integer             :: l, n, nmax, lmax, size, index, h5_err
	character(len = 15) :: lable ! File name without .h5 extension
  character(len = 3)  :: strl  ! file number (l converted to a string)
  character(len = 6)  :: fmt   ! format descriptor
	character(len = 24) :: file_name
	character(len = 12) :: ener_name
  character(len = 25) :: data
	PetscScalar         :: val(1)
	integer(HSIZE_T)    :: ener_dims(1:1)
  PetscReal, allocatable   :: E(:,:)
! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------
  lmax             = 10
  nmax             = 20
	one              = 1
	lable =  'H_test'
  allocate(E(nmax,0:lmax))
	
	! If nmax <= lmax+1 then we use the normal rule but if its larger we 
  ! need to truncate our basis
  if(nmax .le. lmax+1) then
    size = (nmax - 1)*nmax/2 + lmax + 1
  else
		size = (lmax + 1)*(lmax + 2)/2 + (nmax - lmax - 2)*(lmax + 1) +	&
		&	lmax + 1
  endif 
 
  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  if ( ierr /= 0 ) then
    print*,'SlepcInitialize failed'
    stop
  endif

	call MPI_Comm_rank(MPI_COMM_WORLD,proc_id,ierr)
	CHKERRA(ierr)
	call MPI_Comm_size(MPI_COMM_WORLD,num_proc,ierr)
	CHKERRA(ierr)

	call h5open_f( h5_err)
  if ( h5_err /= 0 ) then
		print*,'h5open_f failed'
		stop
  end if

	file_name = trim(lable)//'.h5'

	call h5fopen_f( file_name, H5F_ACC_RDWR_F, file_id, h5_err)

	h5_kind = h5kind_to_type( dp, H5_REAL_KIND)
! --------------------------------------------------------------------------
! Create H0 matrix
! --------------------------------------------------------------------------
	call MatCreate(MPI_COMM_WORLD,H0,ierr)
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

    write(strl,fmt) l

    ener_name ='Energy_l'//trim(strl)

		call h5dopen_f(file_id, ener_name, ener_id, h5_err)
		
		ener_dims(1) = nmax-l
		call h5dread_f( ener_id, h5_kind, E(1:nmax-l,l), ener_dims, h5_err) 
		call h5dclose_f( ener_id, h5_err)
		
    ! For each l value we compute the corresponding matrix elements of H0
		do n=l+1,nmax
			
			! Here I convert the n and l value to its corresponding index 
			if(n .le. lmax+1) then
					index = (n - 1)*n/2 + l
			else
					index = (lmax + 1)*(lmax + 2)/2 + (n - lmax - 2)*(lmax + 1) &
							+ l
			endif    
			
			! Convert the real energy to a PetscScalar
			val = E(n-l,l)
			print*, val
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
	call MatGetSize(H0,size1,size2,ierr)
	CHKERRA(ierr)
  print*, size1,size2
  ! Write the matrix to a binary file
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
      &	"small_fieldFreeMatrix.bin",FILE_MODE_WRITE,viewer,ierr)
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
