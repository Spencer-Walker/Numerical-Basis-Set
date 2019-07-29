! This program generates the dipole acceleration operator 
program main
#include <slepc/finclude/slepceps.h>
  use slepceps
  use simulation_parametersf90
  use hdf5
  implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  PetscErrorCode      :: ierr
  Mat                 :: Z
  PetscReal           :: h,Rmax
  PetscViewer         :: viewer
  PetscMPIInt         :: proc_id, num_proc, comm 
  PetscScalar         :: rint
  PetscInt            :: i_start, i_end, left_index, right_index, index
  PetscInt            :: l_i_start, l_i_end, n, l, i
  PetscInt            :: itter, n1, n2, nmax, l1, l2, lmax, size
  PetscInt            :: num_points, h5_err
  PetscReal           :: start_time, end_time
  PetscReal,      allocatable :: r(:), y(:,:,:)
  PetscScalar,    allocatable :: val_right(:), val_left(:), v1(:), v2(:)
  PetscScalar,    allocatable :: u_right(:,:,:), u_left(:,:,:)
  PetscInt,       allocatable :: col(:)
  PetscReal,      allocatable :: clebsch_gordan(:)
  integer(HID_T)      :: file_id, psi_id, h5_kind
  integer(HSIZE_T)    :: psi_dims(1:3) 
  character(len = 15) :: label ! File name without .h5 extension
  character(len = 3)  :: strl! file number
  character(len = 12) :: psi_name
  character(len = 6)  :: fmt ! format descriptor
  character(len = 30) :: file_name


! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------
  call CPU_TIME(start_time)
  
  10 format(A1,I4)
  20 format(A8,ES9.2)
  
  comm = MPI_COMM_WORLD 
  h = grid_space
  lmax = l_max
  nmax = n_max
  Rmax = R_max
  num_points = int(Rmax/h)
  label = hdf5_file_label 


  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    print*,'SlepcInitialize failed'
    stop
  endif

  call MPI_Comm_rank(comm,proc_id,ierr)
  CHKERRA(ierr)
  call MPI_Comm_size(comm,num_proc,ierr)
  CHKERRA(ierr)

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

  ! If we truncate the basis (set nmax > lmax + 1) then we set the size of
  ! Z differently  
  if(nmax .le. lmax+1) then 
    size = (nmax - 1)*nmax/2 + lmax+ 1
  else
      size = (lmax + 1)*(lmax + 2)/2 + (nmax - lmax - 2)*(lmax + 1) + &
          lmax + 1
  endif

  allocate(clebsch_gordan(lmax))
  allocate(col(nmax))
  allocate(val_right(nmax))
  allocate(val_left(nmax))


! --------------------------------------------------------------------------
! Create Z matrix
! --------------------------------------------------------------------------
  call MatCreate(PETSC_COMM_WORLD,Z,ierr)
  CHKERRA(ierr)
  call MatSetSizes(Z,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetType(Z,eps_mat_type,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(Z,ierr)
  CHKERRA(ierr)
  call MatSetUp(Z,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(Z,i_start,i_end,ierr)
  CHKERRA(ierr)

  l_i_start = -1
  l_i_end = -1
  do l = 0,lmax-1
    do n=l+1,nmax
      index = -1 + n - (l*(1 + l - 2*nmax))/2
      if (index>=i_start .and. l_i_start == -1) then
        l_i_start = l
      end if
      if (index<i_end) then
        l_i_end = l
      end if 
    end do
  end do

  allocate(u_right(num_points,nmax,l_i_start:l_i_end+1))
  allocate(u_left(num_points,nmax,l_i_start:l_i_end+1))
  allocate(v1(num_points))
  allocate(v2(num_points))
  allocate(r(num_points))

  do l = l_i_start,l_i_end+1
    allocate(y(1:num_points,1:nmax-l,2))
    if     (l .le. 9 ) then
      fmt = '(I1.1)'
    elseif (l .le. 99) then
      fmt = '(I2.2)'
    else
      fmt = '(I3.3)'
    endif

    write(strl,fmt) l

    Psi_name = 'Psi_r_l'//trim(strl)

    call h5dopen_f(file_id, psi_name, psi_id, h5_err)

    psi_dims(1) = int(Rmax/h)
    psi_dims(2) = nmax-l
    psi_dims(3) = 2

    call h5dread_f( psi_id, h5_kind, y(1:num_points,1:nmax-l,:),  &
    & psi_dims, h5_err)
    
    u_right(1:num_points,1:nmax-l,l) = dcmplx(y(1:num_points,1:nmax-l,1),y(1:num_points,1:nmax-l,2))
    call h5dclose_f( psi_id, h5_err)
    deallocate(y)
  end do

  do l = l_i_start,l_i_end+1
    allocate(y(1:num_points,1:nmax-l,2))
    if     (l .le. 9 ) then
      fmt = '(I1.1)'
    elseif (l .le. 99) then
      fmt = '(I2.2)'
    else
      fmt = '(I3.3)'
    endif

    write(strl,fmt) l

    Psi_name = 'Psi_l_l'//trim(strl)

    call h5dopen_f(file_id, psi_name, psi_id, h5_err)

    psi_dims(1) = int(Rmax/h)
    psi_dims(2) = nmax-l
    psi_dims(3) = 2

    call h5dread_f( psi_id, h5_kind, y(1:num_points,1:nmax-l,:),  &
    & psi_dims, h5_err)
    
    u_left(1:num_points,1:nmax-l,l) = dcmplx(y(1:num_points,1:nmax-l,1),y(1:num_points,1:nmax-l,2))
    call h5dclose_f( psi_id, h5_err)
    deallocate(y)
  end do

  
! --------------------------------------------------------------------------
! Create r Vector
! --------------------------------------------------------------------------
  r(:) = (/(i*h,i = 1,num_points)/)
  r(:) = r(:)**(-2.d0)
! --------------------------------------------------------------------------
! Load in Clebsch Gordan Coefficients from .bin File
! --------------------------------------------------------------------------

  open(5, file='clebsch_gordan.bin', form='unformatted',access='stream')
  
  do itter = 1,lmax
      read(5, pos=8*itter - 7) clebsch_gordan(itter)
  enddo
  close(5)

  ! Itterate over angular momentum to calculate half the matrix elements
  L_LP_LOOP: do l = l_i_start,l_i_end
    ! For linearly polarized pulses the only non-zero matrix elements are 
    ! for l+/-1
    ! we will compute the (l,l+1) elements only and take care of the rest 
    ! using the fact that Z is hermetian. 
    print 10, 'l', l
    l1 = l
    l2 = l + 1
    ! For each l value we compute the corresponding matrix elements of Z
    do n1=l1+1,nmax
      ! Here I convert the n1 and l1 value to its corresponding index
      left_index = -1 + n1 - (l1*(1 + l1 - 2*nmax))/2
      if (left_index < i_end .and. left_index >= i_start) then 
        ! Here I convert the n2 and l2 value to its corresponding index
        itter = 1
        do n2=l2+1,nmax
          right_index = -1 + n2 - (l2*(1 + l2 - 2*nmax))/2
          ! I create a vector of indicies corresponding to what columns I 
          ! am setting for each row
          col(itter) = right_index
          
          ! Here I compute matrix elements of the r operator as a weighted 
          ! inner product
          v1 = u_left(:,n1-l1,l1)
          v2 = u_right(:,n2-l2,l2)
          rint = sum(v1*r*v2)

          ! Here I compute matrix elements corresponding to the indicies in 
          ! the col vector
          val_right(itter) = 2.d0*((pi/3.d0)**0.5d0)*clebsch_gordan(l1+1)*rint
          if (eps_mat_type .ne. MATSBAIJ) then
            v1 = u_left(:,n2-l2,l2)
            v2 = u_right(:,n1-l1,l1)
            rint = sum(v1*r*v2)

            val_left(itter) = 2.d0*((pi/3.d0)**0.5d0)*clebsch_gordan(l1+1)*rint
          end if 

          itter = itter + 1
        end do

        ! Now that all matrix elements have been computed for the n1 state for
        ! all n2 states I need to 
        ! start again with the first n2 state and compute inner products with 
        ! the nex n1 state
        
        call MatSetValues(Z,1,left_index,itter-1,col,val_right,INSERT_VALUES,ierr)
        CHKERRA(ierr)
        if (eps_mat_type .ne. MATSBAIJ) then
          call MatSetValues(Z,itter-1,col,1,left_index,val_left,INSERT_VALUES,ierr)
          CHKERRA(ierr)
        end if

      else if (left_index >= i_end) then
        exit L_LP_LOOP
      end if 
    end do

    ! Here I destroy the petsc viewers since we are going to move onto the 
    ! next l values
  end do L_LP_LOOP

  ! We finish building Z now that we've finished adding elements
  call MatAssemblyBegin(Z,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)
  call MatAssemblyEnd(Z,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)

! We now save the complete AZ matrix to a binary file on the disk
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  trim(label)//"_dipoleAccelerationMatrix.bin",FILE_MODE_WRITE,viewer,ierr)
  CHKERRA(ierr)

    
  call MatView(Z,viewer,ierr)
  CHKERRA(ierr)

  ! Now that Z is saved to memory we clear out memory and exit
  call MatDestroy(Z,ierr)
  CHKERRA(ierr)
  deallocate(val_right,val_left)
  deallocate(col)

  deallocate(clebsch_gordan)
  deallocate(v1)
  deallocate(v2)
  deallocate(u_right,u_left)
  call h5fclose_f( file_id, h5_err )
  call h5close_f( h5_err)
  call SlepcFinalize(ierr)

  call CPU_TIME(end_time)
  
  print 20, 'time   :', end_time-start_time

end program main
    
    