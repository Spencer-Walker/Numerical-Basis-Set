! This program generates the field free hamiltonion in energy basis
program main
#include <slepc/finclude/slepceps.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
use hdf5
use simulation_parametersf90
  use slepceps
  implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  PetscErrorCode      :: ierr
  Vec                 :: yl,yr,tmp
  Mat                 :: H0,W
  PetscInt            :: i_start, i_end, i, n, index, l_i_end, l_i_start
  integer(HID_T)      :: file_id, ener_id, h5_kind, psi_id
  PetscViewer         :: viewer
  PetscMPIInt         :: proc_id, num_proc, comm 
  integer             :: l, n1, n2, nmax, lmax, size, left_index, h5_err, num_points
  integer             :: ecs_point, itter
  character(len = 15) :: label ! File name without .h5 extension
  character(len = 3)  :: strl  ! file number (l converted to a string)
  character(len = 6)  :: fmt   ! format descriptor
  character(len = 30) :: file_name
  character(len = 12) :: ener_name, psi_name
  PetscScalar         :: absorber
  PetscReal           :: sig, width
  PetscScalar, allocatable :: val(:),vl(:),vr(:)
  PetscInt,    allocatable :: col(:),indicies(:)
  integer(HSIZE_T)    :: ener_dims(1:1), psi_dims(1:2)
  PetscReal, allocatable   :: u(:,:),E(:,:),El(:)
  real(dp)  :: start_time, end_time, Rmax, h, Recs
! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------
  call CPU_TIME(start_time)

  10 format(A1,I4)
  20 format(A8,ES9.2)

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    print*,'Unable to initialize PETSc'
    stop
  endif

  comm  = MPI_COMM_WORLD
  h     = grid_space
  lmax  = l_max
  nmax  = n_max
  Rmax  = R_max
  label = hdf5_file_label 
  num_points = int(Rmax/h)
  Recs  = gobbler*Rmax
  ecs_point = int(Recs/h)

  print*, 'num_points', num_points
  

  call MPI_Comm_rank(comm,proc_id,ierr)
  CHKERRA(ierr)
  call MPI_Comm_size(comm,num_proc,ierr)
  CHKERRA(ierr)

  if (ecs_present) then
    absorber = 1.d0
  else 
    absorber = 0.d0
  end if 

  call MatCreate(PETSC_COMM_SELF,W,ierr)
  CHKERRA(ierr)
  call MatSetSizes(W,PETSC_DECIDE,PETSC_DECIDE,num_points,num_points,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(W,ierr)
  CHKERRA(ierr)
  call MatSetUp(W,ierr)
  CHKERRA(ierr)

  call VecCreateSeq(PETSC_COMM_SELF,num_points,yl,ierr)
  CHKERRA(ierr)
  call VecDuplicate(yl,yr,ierr)
  CHKERRA(ierr)
  call VecDuplicate(yl,tmp,ierr)
  CHKERRA(ierr)

  if (ecs_present) then
    allocate(col(1),val(1))

    width = (R_max - Recs)/3.d0
    sig = 1.859d0

    do i = 1, num_points

      col(1) = i - 1

      val(1) = dcmplx(0,-V_max*(dtanh(sig*(dble(i)*h-Recs)/width-sig)+1.d0)/2.d0)
      
      call MatSetValues(W,1,i-1,1,col,val,INSERT_VALUES,ierr)
      CHKERRA(ierr)
    end do
    deallocate(col,val)
  end if

  call MatAssemblyBegin(W,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)
  call MatAssemblyEnd(W,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)

  
  allocate(E(nmax,0:lmax))
  allocate(u(num_points,nmax),El(nmax))
  
  ! If nmax <= lmax+1 then we use the normal rule but if its larger we 
  ! need to truncate our basis
  if(nmax .le. lmax+1) then
    size = (nmax - 1)*nmax/2 + lmax + 1
  else
    size = (lmax + 1)*(lmax + 2)/2 + (nmax - lmax - 2)*(lmax + 1) +&
    lmax + 1
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
  call MatSetType(H0,MATMPISBAIJ,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(H0,ierr)
  CHKERRA(ierr)
  call MatSetUp(H0,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(H0,i_start,i_end,ierr)
  CHKERRA(ierr)

  allocate(col(nmax),val(nmax),indicies(num_points))
  allocate(vl(num_points),vr(num_points))

  do i = 1,num_points
    indicies(i) = i-1
  end do 

  l_i_start = -1
  l_i_end = -1
  do l = 0,lmax
    do n=l+1,nmax
      index = -1 + n - (l*(1 + l - 2*nmax))/2
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
    
    do n1=l+1,nmax    
      vl = u(:,n1-l)
      call  VecSetValues(yl,num_points,indicies,vl,INSERT_VALUES,ierr)
      CHKERRA(ierr)
      call VecAssemblyBegin(yl,ierr)
      CHKERRA(ierr)
      call VecAssemblyEnd(yl,ierr)
      CHKERRA(ierr)

      ! Here I convert the n1 and l value to its corresponding left_index 
      left_index =  -1 + n1 - (l*(1 + l - 2*nmax))/2
      ! Convert the real energy to a PetscScalar
      
      if (ecs_present) then 
        itter = 1
        n2 = n1
        col(itter) =  -1 + n2 - (l*(1 + l - 2*nmax))/2
        vr = u(:,n2-l)
        call VecSetValues(yr,num_points,indicies,vr,INSERT_VALUES,ierr)
        CHKERRA(ierr)
        call VecAssemblyBegin(yr,ierr)
        CHKERRA(ierr)
        call VecAssemblyEnd(yr,ierr)
        CHKERRA(ierr)
        call MatMult(W,yr,tmp,ierr)
        CHKERRA(ierr)
        call VecDot(tmp,yl,val(itter),ierr)
        CHKERRA(ierr)
        val(itter) = val(itter) + E(n1-l,l)
        itter = itter + 1
        do n2 =n1+1,nmax 
          col(itter) =  -1 + n2 - (l*(1 + l - 2*nmax))/2
          vr = u(:,n2-l)
          call VecSetValues(yr,num_points,indicies,vr,INSERT_VALUES,ierr)
          CHKERRA(ierr)
          call VecAssemblyBegin(yr,ierr)
          CHKERRA(ierr)
          call VecAssemblyEnd(yr,ierr)
          CHKERRA(ierr)
          call MatMult(W,yr,tmp,ierr)
          CHKERRA(ierr)
          call VecDot(tmp,yl,val(itter),ierr)
          CHKERRA(ierr)
          itter = itter + 1
        end do 
        call MatSetValues(H0,1,left_index,itter-1,col,val,INSERT_VALUES,ierr)
        CHKERRA(ierr)
      else 
        if (E(n1-l,l) > energy_absorber) then
          val(1) = dcmplx(0,-1d0)*E(n1-l,l)
        else 
          val(1) = E(n1-l,l) 
        end if 
          ! Insert the energy into the field free matrix 
        call MatSetValue(H0,left_index,left_index,val(1),INSERT_VALUES,ierr)
        CHKERRA(ierr)
      end if 
    enddo
  enddo

  deallocate(col,val,indicies,vl,vr)

  print*,'buildig matrix'
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
  call MatDestroy(W,ierr)
  CHKERRA(ierr)
  call VecDestroy(yl,ierr)
  CHKERRA(ierr)
  call VecDestroy(yr,ierr)
  CHKERRA(ierr)
  call VecDestroy(tmp,ierr)
  CHKERRA(ierr)
  deallocate(u,E,El)
  call h5fclose_f( file_id, h5_err)
  call h5close_f( h5_err)
  print*,'blah'
  call PetscFinalize(ierr)
  call CPU_TIME(end_time)
  
  print 20, 'time   :', end_time-start_time

end program main
