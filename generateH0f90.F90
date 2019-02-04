! This program generates the field free hamiltonion in energy basis
program main
use hdf5
use simulation_parametersf90
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
  integer(HID_T)      :: file_id, ener_id, h5_kind, psi_id
  PetscViewer         :: viewer, hdf5v
  PetscMPIInt         :: proc_id, num_proc, comm 
  PetscInt            :: i
  integer             :: l, n, nmax, lmax, size, index, h5_err, num_points
  integer             :: ces_point
  character(len = 15) :: label ! File name without .h5 extension
  character(len = 3)  :: strl  ! file number (l converted to a string)
  character(len = 6)  :: fmt   ! format descriptor
  character(len = 30) :: file_name
  character(len = 12) :: ener_name, psi_name
  character(len = 25) :: data
  PetscScalar         :: val(1)
  integer(HSIZE_T)    :: ener_dims(1:1), psi_dims(1:2)
  PetscReal, allocatable   :: u(:,:),E(:,:),El(:)
  PetscScalar, allocatable :: M(:,:)
  real(dp)  :: start_time, end_time, Rmax, h, Rces
! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------
  call CPU_TIME(start_time)

  10 format(A1,I4)
  20 format(A8,ES9.2)

  comm  = MPI_COMM_WORLD
  h     = grid_space
  lmax  = l_max
  nmax  = n_max
  Rmax  = R_max
  label = hdf5_file_label 
  num_points = int(Rmax/h)
  Rces  = Rmax
  ces_point = int(Rces/h)
  print*, 'num_points', num_points
  

  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  if ( ierr /= 0 ) then
    print*,'SlepcInitialize failed'
    stop
  endif

  call MPI_Comm_rank(comm,proc_id,ierr)
  CHKERRA(ierr)
  call MPI_Comm_size(comm,num_proc,ierr)
  CHKERRA(ierr)

  allocate(E(nmax,0:lmax),M(ces_point-1:num_points,ces_point-1:num_points))
  allocate(u(num_points,nmax),El(nmax))

  M(:,:) = 0.d0
 
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

  M(ces_point,ces_point-1) = -(2.d0/(1.d0+exp(cmplx(0.d0,pi/4.d0))) - 1.d0)&
  & *0.5d0/h**2.d0

  M(ces_point,ces_point)   =  (exp(cmplx(0.d0,-pi/4.d0))  - 1.d0  )*1.d0/h**2.d0  
  M(ces_point,ces_point+1) = -(2.d0*exp(cmplx(0.d0,-pi/4.d0))/(1.d0 &
  & +exp(cmplx(0.d0,pi/4.d0))) - 1.d0)*0.5d0/h**2.d0

  M(num_points,num_points-1) = -( cmplx(0.d0,-1.d0) - 1.d0  )*0.5d0/h**2.d0  
  M(num_points,num_points)   =  ( cmplx(0.d0,-1.d0) - 1.d0  )*1.d0/h**2.d0  


  do i= ces_point+1,num_points-1
    M(i,i-1) = -( cmplx(0.d0,-1.d0) - 1.d0  )*0.5d0/h**2.d0  
    M(i,i) =  ( cmplx(0.d0,-1.d0) - 1.d0  )*1.d0/h**2.d0 
    M(i,i+1) = -( cmplx(0.d0,-1.d0) - 1.d0  )*0.5d0/h**2.d0  
  end do



  ! Itterate over angular momentum to calculate the matrix elements
  do l = 0,lmax  
    print 10, 'l', l
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
    call h5dread_f( ener_id, h5_kind, El(1:nmax-l), ener_dims, h5_err) 
    E(1:nmax-l,l) = El
    ! Closes the dataset
    call h5dclose_f( ener_id, h5_err)
    
    psi_name ='Psi_l'//trim(strl)
    
    call h5dopen_f(file_id, psi_name, psi_id, h5_err)

    psi_dims(1) = int(Rmax/h)
    psi_dims(2) = nmax-l
    call h5dread_f( psi_id, h5_kind, u(1:num_points,1:nmax-l),  &
    & psi_dims, h5_err)
    call h5dclose_f( psi_id, h5_err)
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
      val = E(n-l,l) !+ DOT_PRODUCT(u(ces_point-1:,n-l),MATMUL(M(ces_point-1:,ces_point-1:),u(ces_point-1:,n-l)))
      ! Insert the energy into the field free matrix 
      call MatSetValue(H0,index,index,val,INSERT_VALUES,ierr)
      CHKERRA(ierr)
    enddo
    close(20)
  enddo

  print*,'buildig matrix'
  ! Build the matrix
  call MatAssemblyBegin(H0,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)
  call MatAssemblyEnd(H0,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)

  ! Write the matrix to a binary file
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  &	trim(label)//"_fieldFreeMatrix.bin",FILE_MODE_WRITE,viewer,ierr)
  CHKERRA(ierr)
  call MatView(H0,viewer,ierr)
  CHKERRA(ierr)

  do n = 1, 26
    do l = 0,n-1
      print*, n,E(n-l,l)
    end do
  end do 

  ! Clear memory 
  deallocate(u,E,M,El)
  call MatDestroy(H0,ierr)
  CHKERRA(ierr)
  call h5fclose_f( file_id, h5_err)
  call h5close_f( h5_err)
  call SlepcFinalize(ierr)

  call CPU_TIME(end_time)
  
  print 20, 'time   :', end_time-start_time

end program main
