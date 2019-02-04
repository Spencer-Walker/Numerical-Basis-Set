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
    Mat                 :: Z,ZH
    PetscReal           :: h,Rmax
    PetscInt            :: i_start,i_end,left_index,right_index
    PetscViewer         :: view_l1,view_l2
    PetscMPIInt         :: proc_id, num_proc, comm 
    PetscReal           :: rint
    PetscInt            :: i,j
    integer(HID_T)      :: file_id, psi_id, h5_kind
    integer             :: l,itter,n1,n2,nmax,l1,l2,lmax,size
    integer             :: num_points,fd_l1,fd_l2, h5_err
    character(len = 15) :: label ! File name without .h5 extension
    character(len = 3)  :: strl! file number
    character(len = 12) :: psi_name
    character(len = 6)  :: fmt ! format descriptor
    character(len = 30) :: file_name
    PetscReal,      allocatable :: u(:,:,:),v1(:),v2(:),r(:),y(:,:)
    PetscInt,       allocatable :: col(:)
    PetscScalar,    allocatable :: val(:)
    real(dp),       allocatable :: clebsch_gordan(:)
    integer(HSIZE_T) :: psi_dims(1:2) 
    real(dp) :: start_time, end_time
  
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
  
    allocate(u(num_points,nmax,0:lmax))
    allocate(v1(num_points))
    allocate(v2(num_points))
    allocate(r(num_points))
  
    do l = 0,lmax
      allocate(y(1:num_points,1:nmax-l))
      if     (l .le. 9 ) then
        fmt = '(I1.1)'
      elseif (l .le. 99) then
        fmt = '(I2.2)'
      else
        fmt = '(I3.3)'
      endif
  
      write(strl,fmt) l
  
      Psi_name = 'Psi_l'//trim(strl)
  
      call h5dopen_f(file_id, psi_name, psi_id, h5_err)
  
      psi_dims(1) = int(Rmax/h)
      psi_dims(2) = nmax-l
  
      call h5dread_f( psi_id, h5_kind, y(1:num_points,1:nmax-l),  &
      & psi_dims, h5_err)
      
      u(1:num_points,1:nmax-l,l) = y(1:num_points,1:nmax-l)
      call h5dclose_f( psi_id, h5_err)
      deallocate(y)
    end do
  
  
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
    allocate(val(nmax))
   
  
  
   
  ! --------------------------------------------------------------------------
  ! Create r Vector
  ! --------------------------------------------------------------------------
    r(:) = (/(i*h,i = 1,num_points)/)
    r(:) = r(:)**-2.d0
  ! --------------------------------------------------------------------------
  ! Load in Clebsch Gordan Coefficients from .bin File
  ! --------------------------------------------------------------------------
  
    open(5, file='clebsch_gordan.bin', form='unformatted',access='stream')
    
    do itter = 1,lmax
       read(5, pos=8*itter - 7) clebsch_gordan(itter)
    enddo
    close(5)
  
  ! --------------------------------------------------------------------------
  ! Create Z matrix
  ! --------------------------------------------------------------------------
    call MatCreate(PETSC_COMM_WORLD,Z,ierr)
    CHKERRA(ierr)
    call MatSetSizes(Z,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
    CHKERRA(ierr)
    call MatSetFromOptions(Z,ierr)
    CHKERRA(ierr)
    call MatSetUp(Z,ierr)
    CHKERRA(ierr)
    call MatGetOwnershipRange(Z,i_start,i_end,ierr)
    CHKERRA(ierr)
   
    ! Itterate over angular momentum to calculate half the matrix elements
    do l = 0,lmax-1
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
        if(n1 .le. lmax+1) then
          left_index = (n1 - 1)*n1/2 + l1
        else
          left_index = (lmax + 1)*(lmax + 2)/2 + (n1 - lmax - 2)*(lmax &
          & +1) + l1
        endif
        if (left_index <= i_end .and. left_index >= i_start) then 
          ! Here I convert the n2 and l2 value to its corresponding index
          itter = 1
          do n2=l2+1,nmax
            if(n2 .le. lmax+1) then
              right_index = (n2 - 1)*n2/2 + l2
            else
              right_index = (lmax + 1)*(lmax + 2)/2 + (n2 - lmax - 2)* &
              & (lmax + 1) + l2
            endif
            ! I create a vector of indicies corresponding to what columns I 
            ! am setting for each row
            col(itter) = right_index
            
            ! Here I compute matrix elements of the r operator as a weighted 
            ! inner product
            v1 = u(:,n1-l1,l1)
            v2 = u(:,n2-l2,l2)
            rint = DOT_PRODUCT(v1,r*v2)
            ! Here I compute matrix elements corresponding to the indicies in 
            ! the col vector
            val(itter) = dcmplx(2.0*((pi/3.0)**0.5)*clebsch_gordan(l1+1)*rint,0)
            itter = itter + 1
          enddo
          ! Now that all matrix elements have been computed for the n1 state for
          ! all n2 states I need to 
          ! start again with the first n2 state and compute inner products with 
          ! the nex n1 state
          
          call MatSetValues(Z,1,left_index,itter-1,col,val,INSERT_VALUES,ierr)
          CHKERRA(ierr)
        end if 
      enddo
  
      ! Here I destroy the petsc viewers since we are going to move onto the 
      ! next l values
    enddo
  
    ! We finish building Z now that we've finished adding elements
    call MatAssemblyBegin(Z,MAT_FINAL_ASSEMBLY,ierr)
    CHKERRA(ierr)
    call MatAssemblyEnd(Z,MAT_FINAL_ASSEMBLY,ierr)
    CHKERRA(ierr)
  
    ! Since we only added half the elements we can get the rest by adding Z to 
    ! its conjugate transpose 
    call MatDuplicate(Z,MAT_DO_NOT_COPY_VALUES,ZH,ierr)
    CHKERRA(ierr)
    call MatHermitianTranspose(Z,MAT_INITIAL_MATRIX,ZH,ierr)
    CHKERRA(ierr)
    call MatAXPY(Z,one,ZH,DIFFERENT_NONZERO_PATTERN,ierr)
    CHKERRA(ierr)
  
    ! We now save the complete Z matrix to a binary file on the disk
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
      & trim(label)//"_dipoleAccelerationMatrix.bin",FILE_MODE_WRITE,view_l1,ierr);&
      CHKERRA(ierr)
    call MatView(Z,view_l1,ierr)
    CHKERRA(ierr)
    call MatDestroy(ZH,ierr)
    CHKERRA(ierr)
    ! Now that Z is saved to memory we clear out memory and exit
    
    call MatDestroy(Z,ierr)
    CHKERRA(ierr)
    deallocate(val)
    deallocate(col)
  
    deallocate(clebsch_gordan)
    deallocate(v1)
    deallocate(v2)
    deallocate(u)
    call h5fclose_f( file_id, h5_err )
    call h5close_f( h5_err)
    call SlepcFinalize(ierr)
  
    call CPU_TIME(end_time)
    
    print 20, 'time   :', end_time-start_time
  
  end program main