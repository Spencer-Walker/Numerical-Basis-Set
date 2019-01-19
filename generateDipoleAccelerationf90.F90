! This program generates the dipole matrix Z.
program main
#include <slepc/finclude/slepceps.h>
  use slepceps
  implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  PetscErrorCode      :: ierr
  Mat                 :: Z,ZH
  PetscReal           :: h,R0
  PetscInt            :: i_start,i_end,left_index,right_index
  PetscViewer         :: view_l1,view_l2
  PetscReal           :: rint
  PetscScalar         :: one
  PetscInt            :: i,j
  integer             :: l,itter,n1,n2,n_max,l1,l2,l_max,size
  integer             :: num_grid_points,fd_l1,fd_l2
  real(8)             :: pi
  character(len = 3)  :: file! file number
  character(len = 6)  :: fmt ! format descriptor
  character(len = 24) :: file_name
  PetscReal,      allocatable :: u(:,:,:),v1(:),v2(:),r(:)
  PetscInt,       allocatable :: col(:)
  PetscScalar,    allocatable :: val(:)
  real(8),        allocatable :: clebsch_gordan(:)
! --------------------------------------------------------------------------
! Beginning of Program
! --------------------------------------------------------------------------
  h                = 0.065d0
  one              = 1.d0
  pi               = 4.d0*datan(1.d0)
  l_max            = 20
  n_max            = 20
  R0               = 100d0
  num_grid_points  = int(R0/h)
  R0               = num_grid_points*h
  allocate(u(num_grid_points,n_max,0:l_max))
  allocate(v1(num_grid_points))
  allocate(v2(num_grid_points))
  allocate(r(num_grid_points))

  do l = 0,l_max
     if     (l .le. 9 ) then
        fmt = '(I1.1)'
     elseif (l .le. 99) then
        fmt = '(I2.2)'
     else
        fmt = '(I3.3)'
     endif

     write(file,fmt) l
     file_name ='small_vec'//trim(file)//'.bin'
     open(10, file=file_name, form='unformatted',access='sequential')
     do i = 1,n_max-l
        do j = 1,num_grid_points
           read(10) u(j,i,l)
        end do
     end do
     close(10)
  end do


  ! If we truncate the basis (set n_max > l_max + 1) then we set the size of
  ! Z differently  
  if(n_max .le. l_max+1) then 
     size = (n_max - 1)*n_max/2 + l_max+ 1
  else
     size = (l_max + 1)*(l_max + 2)/2 + (n_max - l_max - 2)*(l_max + 1) + &
          l_max + 1
  endif
  
  allocate(clebsch_gordan(l_max))
  allocate(col(n_max))
  allocate(val(n_max))
 

  call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
     print*,'SlepcInitialize failed'
     stop
  endif
 
! --------------------------------------------------------------------------
! Create r Vector
! --------------------------------------------------------------------------
  r(:) = (/(i*h,i = 1,num_grid_points)/)
  r(:) = r(:)**-2.d0
! --------------------------------------------------------------------------
! Load in Clebsch Gordan Coefficients from .bin File
! --------------------------------------------------------------------------

  open(5, file='clebsch_gordan.bin', form='unformatted',access='stream')
  
  do itter = 1,l_max
     read(5, pos=8*itter - 7) clebsch_gordan(itter)
  enddo
  close(5)

! --------------------------------------------------------------------------
! Create Z matrix
! --------------------------------------------------------------------------
  call MatCreate(PETSC_COMM_WORLD,Z,ierr);CHKERRA(ierr)
  call MatSetSizes(Z,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr);CHKERRA(ierr)
  call MatSetFromOptions(Z,ierr);CHKERRA(ierr)
  call MatSetUp(Z,ierr);CHKERRA(ierr)
  call MatGetOwnershipRange(Z,i_start,i_end,ierr);CHKERRA(ierr)
 
  ! Itterate over angular momentum to calculate half the matrix elements
  do l = 0,l_max-1
     ! For linearly polarized pulses the only non-zero matrix elements are 
     ! for l+/-1
     ! we will compute the (l,l+1) elements only and take care of the rest 
     ! using the fact that Z is hermetian. 
     print*,l
     l1 = l
     l2 = l + 1
     ! For each l value we compute the corresponding matrix elements of Z
     do n1=l1+1,n_max
        ! Here I convert the n1 and l1 value to its corresponding index 
        if(n1 .le. l_max+1) then
           left_index = (n1 - 1)*n1/2 + l1
        else
           left_index = (l_max + 1)*(l_max + 2)/2 + (n1 - l_max - 2)*(l_max &
                +1) + l1
        endif
        ! Here I convert the n2 and l2 value to its corresponding index
        itter = 1
        do n2=l2+1,n_max
           if(n2 .le. l_max+1) then
              right_index = (n2 - 1)*n2/2 + l2
           else
              right_index = (l_max + 1)*(l_max + 2)/2 + (n2 - l_max - 2)* &
                   (l_max + 1) + l2
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
        
        call MatSetValues(Z,1,left_index,itter-1,col,val,INSERT_VALUES,ierr);&
             CHKERRA(ierr)
     enddo

     ! Here I destroy the petsc viewers since we are going to move onto the 
     ! next l values
  enddo

  ! We finish building Z now that we've finished adding elements
  call MatAssemblyBegin(Z,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
  call MatAssemblyEnd(Z,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)

  ! Since we only added half the elements we can get the rest by adding Z to 
  ! its conjugate transpose 
  call MatDuplicate(Z,MAT_DO_NOT_COPY_VALUES,ZH,ierr);CHKERRA(ierr)
  call MatHermitianTranspose(Z,MAT_INITIAL_MATRIX,ZH,ierr);CHKERRA(ierr)
  call MatAXPY(Z,one,ZH,DIFFERENT_NONZERO_PATTERN,ierr);CHKERRA(ierr)

  ! We now save the complete Z matrix to a binary file on the disk
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
       "small_dipoleAccelerationMatrix.bin",FILE_MODE_WRITE,view_l1,ierr);&
       CHKERRA(ierr)
  call MatView(Z,view_l1,ierr);CHKERRA(ierr)
  call MatDestroy(ZH,ierr);CHKERRA(ierr)
  ! Now that Z is saved to memory we clear out memory and exit
  
  call MatDestroy(Z,ierr);CHKERRA(ierr)
  deallocate(val)
  deallocate(col)

  deallocate(clebsch_gordan)
  deallocate(v1)
  deallocate(v2)
  deallocate(u)
  call SlepcFinalize(ierr)

end program main
