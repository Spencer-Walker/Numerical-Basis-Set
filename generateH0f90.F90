! This program generates the field free hamiltonion in energy basis
program main
#include <slepc/finclude/slepceps.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

  use slepceps
  implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  PetscErrorCode      :: ierr
  Mat                 :: H0
  PetscInt            :: i_start, i_end
  PetscViewer         :: viewer
  PetscMPIInt         :: rank 
  PetscInt            :: one, size1, size2, i
  integer             :: l, n, n_max, l_max, size, index, h5_err
  character(len = 3)  :: file! file number
  character(len = 6)  :: fmt ! format descriptor
  character(len = 24) :: file_name
  character(len = 25) :: data
  PetscScalar         :: val(1)
  PetscReal, allocatable   :: E(:,:)
! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------
  l_max            = 20
  n_max            = 20
  one              = 1
  allocate(E(n_max,0:l_max))
  ! If n_max <= l_max+1 then we use the normal rule but if its larger we 
  ! need to truncate our basis
  if(n_max .le. l_max+1) then
     size = (n_max - 1)*n_max/2 + l_max + 1
  else
     size = (l_max + 1)*(l_max + 2)/2 + (n_max - l_max - 2)*(l_max + 1) + &
          l_max + 1
  endif 
 
  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  if ( ierr /= 0 ) then
     print*,'SlepcInitialize failed'
     stop
  endif

  call h5open_f( h5_err)
  if ( h5_err /= 0 ) then
    print*,'h5open_f failed'
  end if


  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);CHKERRA(ierr)

  
! --------------------------------------------------------------------------
! Create H0 matrix
! --------------------------------------------------------------------------
  call MatCreate(MPI_COMM_WORLD,H0,ierr);CHKERRA(ierr)
  call MatSetSizes(H0,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr);&
       CHKERRA(ierr)
  call MatSetFromOptions(H0,ierr);CHKERRA(ierr)
  call MatSetUp(H0,ierr);CHKERRA(ierr)
  call MatGetOwnershipRange(H0,i_start,i_end,ierr);CHKERRA(ierr)

  ! Itterate over angular momentum to calculate the matrix elements
  do l = 0,l_max
     
     ! Figure out how many digits are needed for the input file
      if     (l .le. 9 ) then
        fmt = '(I1.1)'
     elseif (l .le. 99) then
        fmt = '(I2.2)'
     else
        fmt = '(I3.3)'
     endif

     write(file,fmt) l

     file_name ='small_val'//trim(file)//'.bin'
     open(10, file=file_name, form='unformatted',access='sequential')
     do i = 1,n_max-l
        read(10) E(i,l)
        print*, i
     end do
     close(10)

     ! For each l value we compute the corresponding matrix elements of H0
      do n=l+1,n_max
       
        ! Here I convert the n and l value to its corresponding index 
        if(n .le. l_max+1) then
           index = (n - 1)*n/2 + l
        else
           index = (l_max + 1)*(l_max + 2)/2 + (n - l_max - 2)*(l_max + 1) &
                + l
        endif    
        
        ! Convert the real energy to a PetscScalar
        val = E(n-l,l)
        print*, val
        ! Insert the energy into the field free matrix 
        call MatSetValue(H0,index,index,val,INSERT_VALUES,ierr);&
             CHKERRA(ierr)
     enddo
     close(20)
  enddo

  ! Build the matrix
  call MatAssemblyBegin(H0,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyEnd(H0,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatGetSize(H0,size1,size2,ierr);CHKERRA(ierr)
  print*, size1,size2
  ! Write the matrix to a binary file
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
       "small_fieldFreeMatrix.bin",FILE_MODE_WRITE,viewer,ierr);&
       CHKERRA(ierr)
  call MatView(H0,viewer,ierr);CHKERRA(ierr)
  ! Clear memory 
  call MatDestroy(H0,ierr);CHKERRA(ierr)
  call SlepcFinalize(ierr)
end program main
