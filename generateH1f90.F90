! This program generates the dipole acceleration operator 
function Indicies(n,l,m,nmax,mmax) 
  implicit none
  integer, intent(in)  :: n,l,m,nmax,mmax
  integer :: Indicies
  integer :: first, rest
  if (l .le. mmax + 1) then
    first = (2*nmax-1)*l*(l-1)/2+ nmax*l - 2*(l-1)*l*(2*l-1)/6
  else 
    first = (2*nmax-1)*mmax*(mmax+1)/2+ nmax*(mmax+1) - 2*mmax*(mmax+1)*(2*mmax+1)/6
    first = first + (2*mmax+1)*( nmax*(l-1-mmax) - l*(l-1)/2 + mmax*(mmax + 1)/2)
  end if
  rest =  (m + min(l,mmax))*(nmax-l) + n - l - 1;
  Indicies = first + rest
  return
end function Indicies

function Transfer_Coeff(m1,m2)
#include <slepc/finclude/slepceps.h>
  use slepceps
  implicit none
  integer, intent(in) :: m1, m2
  PetscScalar :: Transfer_Coeff
  if (abs(m1) .ne. abs(m2)) then
    Transfer_Coeff = 0d0 
    return
  else if ( m1 .eq. 0 .and. m2 .eq. 0 ) then 
    Transfer_Coeff = 1d0
    return
  else if ( m1 .gt. 0 .and. m2 .gt. 0 ) then
    Transfer_Coeff = (-1d0)**dble(m1)/dsqrt(2.d0) 
    return
  else if ( m1 .gt. 0 .and. m2 .lt. 0 ) then
    Transfer_Coeff = 1.d0/dsqrt(2.d0)
    return
  else if ( m1 .lt. 0 .and. m2 .gt. 0 ) then
    Transfer_Coeff = dcmplx(0d0,-1d0)*(-1d0)**dble(m1)/dsqrt(2d0)
    return
  else if ( m1 .lt. 0 .and. m2 .lt. 0 ) then
    Transfer_Coeff = dcmplx(0d0,1d0)/dsqrt(2d0)
    return 
  end if
end function Transfer_Coeff

function ThreeJ(l1,l2,l3,m1,m2,m3)
  use fgsl
  use iso_c_binding
  implicit none 
  integer, intent(in) :: l1, l2, l3, m1, m2, m3
  real(fgsl_double) :: ThreeJ 
  ThreeJ = fgsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3)
  return 
end function ThreeJ

function Angular(l1,k,l2,m1,q,m2)
#include <slepc/finclude/slepceps.h>
  use fgsl
  use iso_c_binding
  use slepceps
  implicit none
  real(fgsl_double)   :: ThreeJ
  PetscScalar :: Transfer_Coeff, Angular
  integer, intent(in) :: l1,m1,k,q,l2,m2
  integer :: mm, qq
  Angular = 0d0 

  if ( m1 .eq. 0 .and. q .eq. 0 ) then
    Angular = Angular + ThreeJ(l1,k,l2,-m1,0,m2)
  elseif ( m1 .eq. 0 ) then
    do qq = -abs(q),abs(q), 2*abs(q)
      Angular = Angular + Transfer_Coeff(m2,-qq)*Transfer_Coeff(q,qq)*ThreeJ(l1,k,l2,0,qq,-qq)
    enddo
  elseif ( q .eq. 0) then
    do mm = -abs(m1),abs(m1), 2*abs(m1)
      Angular = Angular + dconjg(Transfer_Coeff(m1,mm))*Transfer_Coeff(m2,mm)*(-1d0)**dble(mm)*ThreeJ(l1,k,l2,-mm,0,mm)
    enddo
  else
    do mm = -abs(m1),abs(m1), 2*abs(m1)
      do qq = -abs(q),abs(q), 2*abs(q)
        Angular = Angular + dconjg(Transfer_Coeff(m1,mm))*Transfer_Coeff(m2,mm-qq)*Transfer_Coeff(q,qq)*(-1d0)**dble(mm)*ThreeJ(l1,k,l2,-mm,qq,mm-qq)
      enddo
    enddo
  end if

  Angular =  dsqrt(dble((2*l1+1)*(2*l2+1)))*ThreeJ(l1,k,l2,0,0,0)*Angular
  
  return
end function Angular 

program main
#include <slepc/finclude/slepceps.h>
use slepceps
use hdf5
use fgsl
use iso_c_binding
#if defined(__INTEL_COMPILER)
use ifport
#endif
implicit none
! --------------------------------------------------------------------------
! Declarations
! --------------------------------------------------------------------------
  PetscInt,   parameter :: dp = kind(1.d0)
  PetscErrorCode      :: ierr
  Mat                 :: X, Y, Z
  PetscReal           :: h,Rmax
  PetscViewer         :: viewer
  PetscMPIInt         :: proc_id, num_proc, comm 
  PetscScalar, allocatable  :: rint_left(:,:),rint(:,:)
  PetscInt            :: i_start, i_end, left_index, right_index, index
  PetscInt            :: l_i_start, l_i_end, m_i_start, m_i_end, n, l, m, m1, m2, i
  PetscInt            :: m1_low, m1_high
  PetscInt            :: itter_x,itter_y,itter_z, n1, n2, nmax, l1, l2, lmax, size
  PetscInt            :: num_points, h5_err, tdse_nmax, tdse_lmax, tdse_mmax
  PetscInt            :: status, basis_local
  PetscReal           :: start_time, end_time
  PetscReal,      allocatable :: r(:), wfn(:,:,:)
  PetscScalar,    allocatable :: val_z(:), v1(:), v2(:)
  PetscScalar,    allocatable :: val_x(:)
  PetscScalar,    allocatable :: val_y(:)
  PetscScalar,    allocatable :: u(:,:,:)
  PetscInt,       allocatable :: col_x(:),col_y(:),col_z(:)
  logical             :: skip
  integer(HID_T)      :: psi_id, h5_kind
  integer(HSIZE_T)    :: psi_dims(1:3), dims(1)
  integer(SIZE_T), parameter :: sdim = 300 
  integer(HID_T)      :: file_id, param_file_id
  integer(HID_T)      :: eps_group_id, memtype, eps_dat_id
  integer(HID_T)      :: install_dat_id
  integer(HID_T)      :: operators_group_id, operators_dat_id
  integer(HID_T)      :: tdse_group_id, tdse_dat_id
  real(fgsl_double)   :: threej_m0, ThreeJ
  character(len = 15) :: label ! File name without .h5 extension
  character(len = 3)  :: strl! file number
  character(len = 12) :: psi_name
  character(len = 6)  :: fmt ! format descriptor
  character(len = 30) :: file_name
  character(len = 300):: tmp_character, basis_directory, working_directory
  character(len = 300):: install_directory
  MatType :: mat_type
  PetscReal,  parameter :: pi = 3.141592653589793238462643383279502884197169
  PetscInt,   allocatable  :: block_n(:), block_l(:), block_m(:)
  integer(HID_T)      :: block_group_id, block_dat_id
  PetscInt            :: num_block
  integer :: Indicies, nnzx, nnzy, nnzz
  PetscScalar :: Angular, blah
  complex(16) :: blah2
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

  dims(1) = 1
  call H5Tcopy_f(H5T_FORTRAN_S1, memtype, h5_err)
  call H5Tset_size_f(memtype, sdim, h5_err)

  call h5dopen_f(param_file_id, "install_directory", install_dat_id, h5_err)
  call h5dread_f(install_dat_id, memtype, install_directory, dims, h5_err)
  call h5dclose_f(install_dat_id, h5_err)

  call h5gopen_f(param_file_id, "EPS", eps_group_id, h5_err)

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

  call h5dopen_f(tdse_group_id, "m_max", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_INTEGER, tdse_mmax, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call h5dopen_f(tdse_group_id, "n_max", tdse_dat_id, h5_err)
  call h5dread_f(tdse_dat_id, H5T_NATIVE_INTEGER, tdse_nmax, dims, h5_err)
  call h5dclose_f( tdse_dat_id, h5_err)

  call h5gopen_f(param_file_id, "operators", operators_group_id, h5_err)
  call h5dopen_f(operators_group_id, "mat_type", operators_dat_id, h5_err)
  call h5dread_f(operators_dat_id, memtype, tmp_character,dims, h5_err)
  call h5dclose_f(operators_dat_id, h5_err)
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
  allocate(block_m(num_block))

  if (num_block .ne. 0) then 
    dims(1) = num_block
  
    call h5dopen_f(block_group_id, "n_index", block_dat_id, h5_err)
    call h5dread_f(block_dat_id, H5T_NATIVE_INTEGER, block_n, dims, h5_err)
    call h5dclose_f(block_dat_id, h5_err)

    call h5dopen_f(block_group_id, "l_index", block_dat_id, h5_err)
    call h5dread_f(block_dat_id, H5T_NATIVE_INTEGER, block_l, dims, h5_err)
    call h5dclose_f(block_dat_id, h5_err)
    
    call h5dopen_f(block_group_id, "m_index", block_dat_id, h5_err)
    call h5dread_f(block_dat_id, H5T_NATIVE_INTEGER, block_m, dims, h5_err)
    call h5dclose_f(block_dat_id, h5_err)

  end if 
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
  call h5fopen_f( trim(file_name), H5FD_MPIO_INDEPENDENT_F, file_id, h5_err)

  ! Converts fortrans double kind to an hdf5 double type 
  h5_kind = h5kind_to_type( dp, H5_REAL_KIND)

  ! Total dim of the hilbert space for this problem.
  size = Indicies(tdse_nmax,tdse_lmax,tdse_mmax,tdse_nmax,tdse_mmax) + 1 
  
  allocate(col_x(tdse_nmax))
  allocate(col_y(tdse_nmax))
  allocate(col_z(tdse_nmax))
  allocate(val_z(tdse_nmax))
  allocate(val_x(tdse_nmax))
  allocate(val_y(tdse_nmax))

! --------------------------------------------------------------------------
! Create Z, X, Y matrix
! --------------------------------------------------------------------------
  call MatCreate(PETSC_COMM_WORLD,Z,ierr)
  CHKERRA(ierr)
  call MatSetSizes(Z,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
  CHKERRA(ierr)
  call MatSetType(Z,mat_type,ierr)
  CHKERRA(ierr)
  if (tdse_mmax .eq. 0) then
    call MatMPISBAIJSetPreallocation(Z,1,tdse_nmax+1,PETSC_NULL_INTEGER,tdse_nmax,PETSC_NULL_INTEGER,ierr)
    CHKERRA(ierr)
  else 
    call MatMPISBAIJSetPreallocation(Z,1,3*tdse_nmax+1,PETSC_NULL_INTEGER,3*tdse_nmax,PETSC_NULL_INTEGER,ierr)
    CHKERRA(ierr)
  end if 
  call MatSetFromOptions(Z,ierr)
  CHKERRA(ierr)
  call MatSetUp(Z,ierr)
  CHKERRA(ierr)
  call MatGetOwnershipRange(Z,i_start,i_end,ierr)
  CHKERRA(ierr)

  if (tdse_mmax .ne. 0) then
    call MatCreate(PETSC_COMM_WORLD,X,ierr)
    CHKERRA(ierr)
    call MatSetSizes(X,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
    CHKERRA(ierr)
    call MatSetType(X,mat_type,ierr)
    CHKERRA(ierr)
    call MatMPISBAIJSetPreallocation(X,1,3*tdse_nmax+1,PETSC_NULL_INTEGER,3*tdse_nmax,PETSC_NULL_INTEGER,ierr)
    CHKERRA(ierr)
    call MatSetFromOptions(X,ierr)
    CHKERRA(ierr)
    call MatSetUp(X,ierr)
    CHKERRA(ierr)

    call MatCreate(PETSC_COMM_WORLD,Y,ierr)
    CHKERRA(ierr)
    call MatSetSizes(Y,PETSC_DECIDE,PETSC_DECIDE,size,size,ierr)
    CHKERRA(ierr)
    call MatSetType(Y,mat_type,ierr)
    CHKERRA(ierr)
    call MatMPISBAIJSetPreallocation(Y,1,3*tdse_nmax+1,PETSC_NULL_INTEGER,3*tdse_nmax,PETSC_NULL_INTEGER,ierr)
    CHKERRA(ierr)
    call MatSetFromOptions(Y,ierr)
    CHKERRA(ierr)
    call MatSetUp(Y,ierr)
    CHKERRA(ierr)
  end if 
  
  l_i_start = -1
  l_i_end = -1
  do l = 0,tdse_lmax-1
    do m = -min(tdse_mmax,l),min(tdse_mmax,l)
      do n=l+1,tdse_nmax
        index = Indicies(n,l,m,tdse_nmax,tdse_mmax)
        if (index>=i_start .and. l_i_start == -1) then
          l_i_start = l
          m_i_start = m
        end if
        if (index<=i_end) then
          l_i_end = l
          m_i_end = m
        end if 
      end do
    end do
  end do 

  if (l_i_start .eq. -1) then
    l_i_start = tdse_lmax - 1
    m_i_start = -min(tdse_mmax,tdse_lmax-1)
  end if 

  allocate(u(num_points,nmax,l_i_start:l_i_end+1))
  allocate(v1(num_points))
  allocate(v2(num_points))
  allocate(r(num_points))

  do l = l_i_start,l_i_end+1
    allocate(wfn(1:num_points,1:nmax-l,2))
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

    call h5dread_f( psi_id, h5_kind, wfn(1:num_points,1:nmax-l,:),  &
    & psi_dims, h5_err)
    
    u(1:num_points,1:nmax-l,l) = dcmplx(wfn(1:num_points,1:nmax-l,1),wfn(1:num_points,1:nmax-l,2))
    call h5dclose_f( psi_id, h5_err)
    deallocate(wfn)
  end do
 
! --------------------------------------------------------------------------
! Create r Vector
! --------------------------------------------------------------------------
  r(:) = (/(i*h,i = 1,num_points)/)
! --------------------------------------------------------------------------
  ! Itterate over angular momentum to calculate half the matrix elements
  !open(1, file = 'data1.dat', status = 'new')
  nnzx = 0
  nnzy = 0
  nnzz = 0

  L_LOOP: do l = l_i_start,l_i_end
    l1 = l
    l2 = l+1
    
    threej_m0 = ThreeJ(l1,1,l2,0,0,0)
    
    LOW_CONDITION: if (l1 .eq. l_i_start) then
      m1_low  = m_i_start
    else
      m1_low = -min(tdse_mmax,l1)
    end if LOW_CONDITION

    HIGH_CONDITION: if (l1 .eq. l_i_end) then 
      m1_high = m_i_end 
    else 
      m1_high = min(tdse_mmax,l1)
    end if HIGH_CONDITION
    
    allocate(rint(tdse_nmax-l1,tdse_nmax-l2))
    allocate(rint_left(tdse_nmax-l2,tdse_nmax-l1))

    N1_INT_LOOP: do n1 = l1+1, tdse_nmax
      N2_INT_LOOP: do n2 = l2+1, tdse_nmax
        v1 = u(:,n1-l1,l1)
        v2 = u(:,n2-l2,l2)
        rint(n1-l1,n2-l2) = sum(v1*r*v2)

      end do N2_INT_LOOP
    end do N1_INT_LOOP

    M1_LOOP: do m1 = m1_low,m1_high
      write(tmp_character, "(A2,I4,A4,I4)")  'l1', l1, ', m1', m1
      call PetscPrintf(MPI_COMM_SELF, trim(tmp_character)//"\n", ierr)
      CHKERRA(ierr)
      N1_LOOP: do n1 = l1+1,tdse_nmax
        left_index = Indicies(n1,l1,m1,tdse_nmax,tdse_mmax)
        RANGE_CONDITION: if (left_index < i_end .and. left_index >= i_start) then
          M2_LOOP: do m2 = -max(tdse_mmax,l2),min(tdse_mmax,l2)
            itter_x = 1
            itter_y = 1
            itter_z = 1
            N2_LOOP: do n2 = l2+1,tdse_nmax
              skip = .false.  
              BLOCK_CONDITION: if ( num_block .ne. 0 )  then
                BLOCK_LOOP: do i = 1,num_block
                  ROW_CONDITION: if ( (block_n(i) .eq. n1) .and. block_l(i) .eq. l1 .and. block_m(i) .eq. m1) then
                    skip = .true.
                  end if ROW_CONDITION
                  COLUMN_CONDITION: if ( (block_n(i) .eq. n2) .and. block_l(i) .eq. l2 .and. block_m(i) .eq. m2) then
                    skip = .true.
                  end if COLUMN_CONDITION
                end do BLOCK_LOOP                
              end if BLOCK_CONDITION
              SKIP_CONDITION: if (skip .eqv. .false.) then
                right_index = Indicies(n2,l2,m2,tdse_nmax,tdse_mmax)
                col_x(itter_x) = right_index
                col_y(itter_y) = right_index
                col_z(itter_z) = right_index
                val_z(itter_z) = Angular(l1,1,l2,m1,0,m2)*rint(n1-l1,n2-l2)
                
                if (zabs(val_z(itter_z)) .gt. 1d-16) then 
                  itter_z = itter_z + 1
                  nnzz = nnzz + 1
                endif

                if (tdse_mmax .ne. 0) then 
                  val_x(itter_x) = Angular(l1,1,l2,m1,1,m2)*rint(n1-l1,n2-l2)
                  if (zabs(val_x(itter_x)) .gt. 1d-16 ) then 
                    itter_x = itter_x + 1
                    nnzx = nnzx + 1
                  endif
                  val_y(itter_y) =Angular(l1,1,l2,m1,-1,m2)*rint(n1-l1,n2-l2)
                  if (zabs(val_y(itter_y)) .gt. 1d-16 ) then 
                    itter_y = itter_y + 1
                    nnzy = nnzy + 1
                  endif
                end if 
                !write(1,*) l,m1,m2,dble(real( Angular(l1,1,l2,m1,0,m2))),dble(real( Angular(l1,1,l2,m1,1,m2))),dble(real( Angular(l1,1,l2,m1,-1,m2)))
              end if SKIP_CONDITION
            end do N2_LOOP
            call MatSetValues(Z,1,left_index,itter_z-1,col_z,val_z,INSERT_VALUES,ierr)
            CHKERRA(ierr)

            if (tdse_mmax .ne. 0) then 
              call MatSetValues(X,1,left_index,itter_x-1,col_x,val_x,INSERT_VALUES,ierr)
              CHKERRA(ierr)

              call MatSetValues(Y,1,left_index,itter_y-1,col_y,val_y,INSERT_VALUES,ierr)
              CHKERRA(ierr)
            end if 
          end do M2_LOOP
        else if (left_index >= i_end) then
          deallocate(rint_left,rint)
          exit L_LOOP
        end if RANGE_CONDITION
      end do N1_LOOP
    end do M1_LOOP
    deallocate(rint_left,rint)
  end do L_LOOP

  print*, nnzx,nnzy,nnzz

  status = chdir(trim(working_directory))

  call MatAssemblyBegin(Z,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)
  call MatAssemblyEnd(Z,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)



  if (tdse_mmax .ne. 0) then
    call MatAssemblyBegin(X,MAT_FINAL_ASSEMBLY,ierr)
    CHKERRA(ierr)
    call MatAssemblyEnd(X,MAT_FINAL_ASSEMBLY,ierr)
    CHKERRA(ierr)

    call MatAssemblyBegin(Y,MAT_FINAL_ASSEMBLY,ierr)
    CHKERRA(ierr)
    call MatAssemblyEnd(Y,MAT_FINAL_ASSEMBLY,ierr)
    CHKERRA(ierr)
  end if 


  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
  trim(label)//"_Z.bin",FILE_MODE_WRITE,viewer,ierr)
  CHKERRA(ierr)
  call MatView(Z,viewer,ierr)
  CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr)
  CHKERRA(ierr)
  call MatDestroy(Z,ierr)
  CHKERRA(ierr)

  if (tdse_mmax .ne. 0) then
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
    trim(label)//"_X.bin",FILE_MODE_WRITE,viewer,ierr)
    CHKERRA(ierr)
    call MatView(X,viewer,ierr)
    CHKERRA(ierr)
    call PetscViewerDestroy(viewer,ierr)
    CHKERRA(ierr)
    call MatDestroy(X,ierr)
    CHKERRA(ierr)

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,&
    trim(label)//"_Y.bin",FILE_MODE_WRITE,viewer,ierr)
    CHKERRA(ierr)
    call MatView(Y,viewer,ierr)
    CHKERRA(ierr)
    call PetscViewerDestroy(viewer,ierr)
    CHKERRA(ierr)
    call MatDestroy(Y,ierr)
    CHKERRA(ierr)
  end if 

  
  deallocate(val_z)
  deallocate(val_x)
  deallocate(val_y)
  deallocate(col_x)
  deallocate(col_y)
  deallocate(col_z)
  deallocate(v1)
  deallocate(v2)
  deallocate(u)
  deallocate(r)
  deallocate(block_l)
  deallocate(block_n)
  deallocate(block_m)
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
