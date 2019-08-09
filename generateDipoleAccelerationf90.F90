! This program generates the dipole acceleration operator 
program main
#include <slepc/finclude/slepceps.h>
use slepceps
use hdf5
#if defined(__INTEL_COMPILER)
use ifport
#endif
implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  PetscInt,   parameter :: dp = kind(1.d0)
  PetscErrorCode      :: ierr
  Mat                 :: Z
  PetscReal           :: h,Rmax
  PetscViewer         :: viewer
  PetscMPIInt         :: proc_id, num_proc, comm 
  PetscScalar         :: rint
  PetscInt            :: i_start, i_end, left_index, right_index, index
  PetscInt            :: l_i_start, l_i_end, n, l, i
  PetscInt            :: itter, n1, n2, nmax, l1, l2, lmax, size
  PetscInt            :: num_points, h5_err, tdse_nmax, tdse_lmax
  PetscInt            :: status, basis_local
  PetscReal           :: start_time, end_time
  PetscReal,      allocatable :: r(:), y(:,:,:)
  PetscScalar,    allocatable :: val_right(:), val_left(:), v1(:), v2(:)
  PetscScalar,    allocatable :: u_right(:,:,:), u_left(:,:,:)
  PetscInt,       allocatable :: col(:)
  PetscReal,      allocatable :: clebsch_gordan(:)
  logical             :: skip
  integer(HID_T)      :: psi_id, h5_kind
  integer(HSIZE_T)    :: psi_dims(1:3), dims(1)
  integer(SIZE_T), parameter :: sdim = 300 
  integer(HID_T)      :: file_id, param_file_id
  integer(HID_T)      :: eps_group_id, memtype, eps_dat_id
  integer(HID_T)      :: operators_group_id, operators_dat_id
  integer(HID_T)      :: tdse_group_id, tdse_dat_id
  character(len = 15) :: label ! File name without .h5 extension
  character(len = 3)  :: strl! file number
  character(len = 12) :: psi_name
  character(len = 6)  :: fmt ! format descriptor
  character(len = 30) :: file_name
  character(len = 300) :: tmp_character, basis_directory, working_directory
  MatType :: mat_type
  PetscReal,  parameter :: pi = 3.141592653589793238462643383279502884197169
  PetscInt,   allocatable  :: block_n(:), block_l(:)
  integer(HID_T)      :: block_group_id, block_dat_id
  PetscInt            :: num_block
! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------
  call CPU_TIME(start_time)
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    call PetscPrintf(MPI_COMM_WORLD, 'Unable to initialize PETSc\n', ierr)
    CHKERRA(ierr)
    stop
  endif

  ! Opens hdf5 to read in the label.h5 file 
  call h5open_f( h5_err)
  if ( h5_err .ne. 0 ) then
    call PetscPrintf(MPI_COMM_WORLD, 'h5open_f failed\n', ierr)
    CHKERRA(ierr)
    stop
  end if

  comm  = MPI_COMM_WORLD

  call h5fopen_f( trim("parameters.h5"), H5F_ACC_RDONLY_F, param_file_id, h5_err)
  call h5gopen_f(param_file_id, "EPS", eps_group_id, h5_err)

  dims(1) = 1
  call H5Tcopy_f(H5T_FORTRAN_S1, memtype, h5_err)
  call H5Tset_size_f(memtype, sdim, h5_err)

  call h5dopen_f(eps_group_id, "label", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, memtype, label, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "R_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, Rmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)
  
  call h5dopen_f(eps_group_id, "delta_x", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_DOUBLE, h, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "l_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, lmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "n_max", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, nmax, dims, h5_err)
  call h5dclose_f( eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "location", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, memtype, basis_directory, dims, h5_err)
  call h5dclose_f(eps_dat_id, h5_err)

  call h5dopen_f(eps_group_id, "local", eps_dat_id, h5_err)
  call h5dread_f(eps_dat_id, H5T_NATIVE_INTEGER, basis_local, dims, h5_err)
  call h5dclose_f(eps_dat_id, h5_err)

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
    call PetscPrintf(MPI_COMM_WORLD, "mat_type not supported defaulting to MATAIJ\n", ierr)
    CHKERRA(ierr)
    mat_type = MATAIJ   
  end if 

  call h5gopen_f(param_file_id, "block_state", block_group_id, h5_err)
  dims(1) = 1
  call h5dopen_f(block_group_id, "num_block", block_dat_id, h5_err)
  call h5dread_f(block_dat_id, H5T_NATIVE_INTEGER, num_block, dims, h5_err)
  call h5dclose_f(block_dat_id, h5_err)

  allocate(block_n(num_block))
  allocate(block_l(num_block))
  
  dims(1) = num_block

  call h5dopen_f(block_group_id, "n_index", block_dat_id, h5_err)
  call h5dread_f(block_dat_id, H5T_NATIVE_INTEGER, block_n, dims, h5_err)
  call h5dclose_f(block_dat_id, h5_err)

  call h5dopen_f(block_group_id, "l_index", block_dat_id, h5_err)
  call h5dread_f(block_dat_id, H5T_NATIVE_INTEGER, block_l, dims, h5_err)
  call h5dclose_f(block_dat_id, h5_err)

  num_points = nint(Rmax/h)

  call MPI_Comm_rank(comm,proc_id,ierr)
  CHKERRA(ierr)
  call MPI_Comm_size(comm,num_proc,ierr)
  CHKERRA(ierr)
  
  status = getcwd(working_directory)
  status = chdir(trim(basis_directory))
  ! Adds the .h5 extension to the input file 
  file_name = trim(label)//'.h5'

  ! Opens the above file 
  call h5fopen_f( trim(file_name), H5F_ACC_RDWR_F, file_id, h5_err)

  ! Converts fortrans double kind to an hdf5 double type 
  h5_kind = h5kind_to_type( dp, H5_REAL_KIND)

  ! If we truncate the basis (set nmax > lmax + 1) then we set the size of
  ! Z differently  
  if(tdse_nmax .le. tdse_lmax+1) then 
    size = (tdse_nmax - 1)*tdse_nmax/2 + tdse_lmax+ 1
  else
      size = (tdse_lmax + 1)*(tdse_lmax + 2)/2 + (tdse_nmax - tdse_lmax - 2)*(tdse_lmax + 1) + &
        tdse_lmax + 1
  endif

  allocate(clebsch_gordan(tdse_lmax))
  allocate(col(tdse_nmax))
  allocate(val_right(tdse_nmax))
  allocate(val_left(tdse_nmax))


! --------------------------------------------------------------------------
! Create Z matrix
! --------------------------------------------------------------------------
  call MatCreate(PETSC_COMM_WORLD,Z,ierr)
  CHKERRA(ierr)
  call MatSetSizes(Z,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetType(Z,mat_type,ierr)
  CHKERRA(ierr)
  call MatSetFromOptions(Z,ierr)
  CHKERRA(ierr)
  call MatSetUp(Z,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(Z,i_start,i_end,ierr)
  CHKERRA(ierr)

  l_i_start = -1
  l_i_end = -1
  do l = 0,tdse_lmax-1
    do n=l+1,tdse_nmax
      index = -1 + n - (l*(1 + l - 2*tdse_nmax))/2
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
  
  do itter = 1,tdse_lmax
      read(5, pos=8*itter - 7) clebsch_gordan(itter)
  enddo
  close(5)

  ! Itterate over angular momentum to calculate half the matrix elements
  L_LP_LOOP: do l = l_i_start,l_i_end
    ! For linearly polarized pulses the only non-zero matrix elements are 
    ! for l+/-1
    ! we will compute the (l,l+1) elements only and take care of the rest 
    ! using the fact that Z is hermetian. 
    write(tmp_character, "(A2,I4)")  'l1', l
    call PetscPrintf(MPI_COMM_SELF, trim(tmp_character)//"\n", ierr)
    CHKERRA(ierr)
    l1 = l
    l2 = l + 1
    ! For each l value we compute the corresponding matrix elements of Z
    do n1=l1+1,tdse_nmax
      ! Here I convert the n1 and l1 value to its corresponding index
      left_index = -1 + n1 - (l1*(1 + l1 - 2*tdse_nmax))/2
      if (left_index < i_end .and. left_index >= i_start) then 
        ! Here I convert the n2 and l2 value to its corresponding index
        itter = 1
        skip = .false.
        do n2=l2+1,tdse_nmax
          if (num_block .ne. 0) then
            do i = 1,num_block
              if ( (block_n(i) .eq. n1) .and. block_l(i) .eq. l1) then
                skip = .true.
              end if
              if ( (block_n(i) .eq. n2) .and. block_l(i) .eq. l2) then
                skip = .true.
              end if
            end do 
          end if 
          if (.not. skip) then
            right_index = -1 + n2 - (l2*(1 + l2 - 2*tdse_nmax))/2
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
            if (mat_type .ne. MATSBAIJ) then
              v1 = u_left(:,n2-l2,l2)
              v2 = u_right(:,n1-l1,l1)
              rint = sum(v1*r*v2)

              val_left(itter) = 2.d0*((pi/3.d0)**0.5d0)*clebsch_gordan(l1+1)*rint
            end if 

            itter = itter + 1
          end if 
          skip = .false.
        end do

        ! Now that all matrix elements have been computed for the n1 state for
        ! all n2 states I need to 
        ! start again with the first n2 state and compute inner products with 
        ! the nex n1 state
        
        call MatSetValues(Z,1,left_index,itter-1,col,val_right,INSERT_VALUES,ierr)
        CHKERRA(ierr)
        if (mat_type .ne. MATSBAIJ) then
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

  status = chdir(trim(working_directory))
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
  
  deallocate(val_right)
  deallocate(val_left)
  deallocate(col)
  deallocate(clebsch_gordan)
  deallocate(v1)
  deallocate(v2)
  deallocate(u_right)
  deallocate(u_left)
  deallocate(r)
  deallocate(block_l)
  deallocate(block_n)

  call h5fclose_f( file_id, h5_err)
  call h5gclose_f( eps_group_id, h5_err)
  call h5gclose_f( operators_group_id, h5_err)
  call h5gclose_f( tdse_group_id, h5_err)
  call h5gclose_f( block_group_id, h5_err)
  call h5fclose_f( param_file_id, h5_err)
  call h5close_f( h5_err)
  call CPU_TIME(end_time)
  write(tmp_character, "(ES9.2)")  end_time-start_time
  call PetscPrintf(MPI_COMM_WORLD, 'time   :'//trim(tmp_character)//"\n", ierr)
  CHKERRA(ierr)
  call SlepcFinalize(ierr)

end program main