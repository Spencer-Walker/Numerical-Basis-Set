program main
#include <slepc/finclude/slepceps.h>
  use slepceps
  use hdf5
  use mpi
  use simulation_parametersf90
    implicit none
    PetscInt  :: nmax, INFO
    PetscReal :: Rmax
    character(len = 15) :: label ! file name w/o .h5
    integer(HID_T), allocatable :: psi_space_right(:), psi_dset_right(:)
    integer(HID_T), allocatable :: ener_space(:), ener_dset(:)
    integer(HID_T), allocatable :: psi_space_left(:), psi_dset_left(:)
    integer(HID_T) :: file_id
    PetscInt :: h5_err, num_proc, proc_id
    integer(HSIZE_T)  :: ener_dims(1:2)
    integer(HSIZE_T)  :: psi_dims(1:3) 
    integer(HID_T)    :: h5_kind, plist_id
    character(len = 24)   :: file_name 
    character(len = 3)    :: strl ! file number
    character(len = 6)    :: fmt  ! format descriptor  
    Mat            :: H
    Vec            :: Vi
    EPS            :: eps
    EPSType        :: tname
    PetscInt       :: num_points, i, j, Istart, Iend, l_start, l_stop, l_stride
    PetscInt       :: nev, Iabs, ncv, mpd, ncon, maxits, lmax, l
    PetscInt       :: row(1), col(3)
    PetscInt, allocatable :: ix(:), IPIV(:)
    PetscErrorCode :: ierr
    PetscScalar    :: value(3), rint 
    PetscScalar, allocatable :: u_right(:,:), u_left(:,:), E_right(:)
    PetscScalar, allocatable ::V(:),S(:,:),UL(:,:),UR(:,:)
    PetscReal      :: dr, gobbler, theta, tol, start_time, end_time, minimum, maximum
    PetscReal, allocatable :: EE_right(:,:), uu_right(:,:,:)
    PetscReal, allocatable :: uu_left(:,:,:)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     Beginning of program
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
    call CPU_TIME(start_time)

    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
    if (ierr .ne. 0) then
      print*,'SlepcInitialize failed'
      stop
    end if
    call MPI_Comm_rank(PETSC_COMM_WORLD,proc_id,ierr);CHKERRA(ierr)
    call MPI_Comm_size(PETSC_COMM_WORLD,num_proc,ierr);CHKERRA(ierr)
    
    minimum = 1d0
    maximum = 0d0
    Rmax  = R_max
    dr = grid_space
    nmax  = n_max
    maxits = eps_max_its
    tol = eps_tol
    lmax = l_max
    label = hdf5_file_label

    if ( eps_ecs_present ) then
      gobbler = eps_gobbler
      theta  = eps_theta
    else 
      gobbler = 1.0d0
      theta = 0.0d0
    end if 
  
    ! Add the .h5 extension to the file label
    file_name = trim(label)//'.h5'
    allocate(psi_space_right(0:lmax),psi_dset_right(0:lmax))
    allocate(ener_space(0:lmax),ener_dset(0:lmax))
    allocate(psi_space_left(0:lmax),psi_dset_left(0:lmax))

   
    ncv = eps_ncv
    mpd = eps_mpd
    num_points = nint(Rmax/dr)
    Rmax = num_points*dr
    Iabs = floor(gobbler*Rmax/dr)
  
    print*,'h5open_f'
    ! Initialize hdf5 
    call h5open_f( h5_err)
  
    print*,'h5pcreate_f'
    ! Creates the hdf5 property list that will be used for the file 
    ! (like an interface or sae_c template)
    call h5pcreate_f( H5P_FILE_ACCESS_F, plist_id, h5_err)
  
    print*,'h5pset_fapl_mpio_f'
    ! Stores MPI IO communicator information to the file access property 
    ! list. 
    call h5pset_fapl_mpio_f( plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, h5_err)
  
    print*,'h5create_f'
    ! Create the hdf5 file with the property list and parallel acess
    call h5fcreate_f( file_name, H5F_ACC_TRUNC_F, file_id, h5_err, &
    & access_prp = plist_id)
  
    ! Converts fortran real(dp) kind to mpi real numbers
    h5_kind = h5kind_to_type( dp, H5_REAL_KIND)
  
    ! Iterate over all l and create a file for each energy and wfn
    do l = 0, lmax
      ! Figures out how many characters are needed to represent l
      ! (if l = 0 1 character if l = 10 2 characters, ... etc. )
      if ( l .le. 9 ) then
        fmt = '(I1.1)'
      else if ( l .le. 99 ) then
        fmt = '(I2.2)'
      else
        fmt = '(I3.3)'
      end if
      
      ! Create string strl for number l
      write(strl,fmt) l
  
      ! For each l there will be an array of nmax-l energies
      ener_dims(1) = nmax-l
      ener_dims(2) = 2
      ! There will be nmax-l wfns, each with dim int(Rmax/dr)
      psi_dims(1) = int(Rmax/dr)
      psi_dims(2) = nmax-l
      psi_dims(3) = 2
      ! Creates an hdf5 interface that knows information about the data
      call h5screate_simple_f( 2, ener_dims, ener_space(l), h5_err)
      
      ! Connects the datasapce to the hdf5 file and gives it a name 
      call h5dcreate_f( file_id, "Energy_l"//trim(strl), h5_kind, ener_space(l), &
      & ener_dset(l), h5_err)
  
      ! Does the same as before for the energy but with the wfns 
      call h5screate_simple_f( 3, psi_dims, psi_space_right(l), h5_err)
      call h5dcreate_f( file_id, "Psi_r_l"//trim(strl), h5_kind, psi_space_right(l), &
      & psi_dset_right(l), h5_err)
  
      ! Does the same as before for the energy but with the wfns 
      call h5screate_simple_f( 3, psi_dims, psi_space_left(l), h5_err)
      call h5dcreate_f( file_id, "Psi_l_l"//trim(strl), h5_kind, psi_space_left(l), &
      & psi_dset_left(l), h5_err)
  
    end do 
 
    if (eps_world .eqv. PETSC_FALSE) then
      l_start = proc_id
      l_stop = lmax
      l_stride = num_proc
    else
      l_start = 1
      l_stop = lmax
      l_stride = 1 
    end if 
    
    do l = l_start,l_stop,l_stride
      print*,l
      nev = nmax - l
      allocate(V(0:num_points-1))
 
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !     Compute the operator matrix that defines the eigensystem, Ax=kx
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (eps_world .eqv. PETSC_TRUE) then
        call MatCreate(PETSC_COMM_WORLD,H,ierr);CHKERRA(ierr)
      else
        call MatCreate(PETSC_COMM_SELF,H,ierr);CHKERRA(ierr)
      end if
      call MatSetSizes(H,PETSC_DECIDE,PETSC_DECIDE,num_points,num_points,ierr);CHKERRA(ierr)
      call MatSetFromOptions(H,ierr);CHKERRA(ierr)
      call MatSetUp(H,ierr);CHKERRA(ierr)
  
      call MatGetOwnershipRange(H,Istart,Iend,ierr);CHKERRA(ierr)

      do i = 0,num_points-1 
        V(i) = -sae_c0/(dble(i+1)*dr) - (sae_Zc*dexp(-sae_c*(dble(i+1)*dr)))/(dble(i+1)*dr) - &
        sae_a1*dexp(-sae_b1*(dble(i+1)*dr)) - sae_a2*dexp(-sae_b2*(dble(i+1)*dr)) - &
        sae_a3*dexp(-sae_b3*(dble(i+1)*dr)) - sae_a4*dexp(-sae_b4*(dble(i+1)*dr)) - &
        sae_a5*dexp(-sae_b5*(dble(i+1)*dr)) + 0.5d0*dble(l)*(dble(l) + 1.d0)/(dble(i+1)*dr)**2.d0
      end do 

      do i = iabs,num_points-1
        V(i) = V(i) - dcmplx(0d0,eps_eta)*dsin(pi*(dble(i+1)*dr-iabs*dr)/(2d0*(Rmax-iabs*dr)))**2d0
      end do

      if (Istart .eq. 0) then
        row(1) = 0
        col(1) = 0
        col(2) = 1
        value(1) =  1.0d0/(dr**2.0d0) + V(Istart)
        value(2) = -0.5d0/(dr**2.0d0)
        call MatSetValues(H,1,row,2,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        Istart = Istart+1
      end if

      if (Iend .eq. num_points) then
        row(1) = num_points-1
        col(1) = num_points-2
        col(2) = num_points-1
        value(1) = -0.5d0/(dr**2.0d0)
        value(2) =  1.0d0/(dr**2.0d0) + V(num_points-1) 
        call MatSetValues(H,1,row,2,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
        Iend = Iend-1
      end if

      value(1) = -0.5d0/(dr**2.0d0)
      value(2) =  1.0d0/(dr**2.0d0) 
      value(3) = -0.5d0/(dr**2.0d0)
      do i=Istart,Iend-1
        row(1) = i
        col(1) = i-1
        col(2) = i
        col(3) = i+1
        value(2) =  1.0d0/(dr**2.0d0) + V(i)
        call MatSetValues(H,1,row,3,col,value,INSERT_VALUES,ierr);CHKERRA(ierr)
      end do

      call MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY,ierr);CHKERRA(ierr)
      call MatCreateVecs(H,Vi,PETSC_NULL_VEC,ierr);CHKERRA(ierr)
  
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !     Create the eigensolver and display info
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
      !     ** Create eigensolver context
      if (eps_world .eqv. PETSC_TRUE) then
        call EPSCreate(PETSC_COMM_WORLD,eps,ierr);CHKERRA(ierr)
      else
        call EPSCreate(PETSC_COMM_SELF,eps,ierr);CHKERRA(ierr)
      end if 

      call EPSSetTwoSided(eps,PETSC_FALSE ,ierr);CHKERRA(ierr)

      !     ** Set operators. In this case, it is a standard eigenvalue problem
      call EPSSetOperators(eps,H,PETSC_NULL_MAT,ierr);CHKERRA(ierr)
      call EPSSetWhichEigenpairs(eps,eps_which,ierr);CHKERRA(ierr)
      call EPSSetType(eps,eps_type,ierr);CHKERRA(ierr)

      if (eps_which .eq. EPS_TARGET_MAGNITUDE ) then
        call EPSSetTarget(eps,eps_target,ierr);CHKERRA(ierr)
      end if 

      call EPSSetTolerances(eps,tol,maxits,ierr);CHKERRA(ierr)
      call EPSSetProblemType(eps,eps_problem,ierr);CHKERRA(ierr)
      call EPSSetDimensions(eps,nev,ncv,mpd,ierr);CHKERRA(ierr)
  
      !     ** Set solver parameters at runtime
      call EPSSetFromOptions(eps,ierr);CHKERRA(ierr)
  
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !     Solve the eigensystem
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
      call EPSSolve(eps,ierr);CHKERRA(ierr)
  
      !     ** Optional: Get some information from the solver and display it
      call EPSGetType(eps,tname,ierr);CHKERRA(ierr)

      call EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
  
      call EPSGetConverged(eps,ncon,ierr);CHKERRA(ierr)
  
      if (ncon < nev) then
        print*, 'ncon < nev'
        print*, 'ncon = ', ncon, ' nev = ',nev
        stop
      end if

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !     Store solution
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      allocate(u_right(num_points,nev),E_right(nev),ix(num_points))
      allocate(uu_right(num_points,nev,2),EE_right(nev,2))
      allocate(u_left(num_points,nev))
      allocate(uu_left(num_points,nev,2))
      allocate(S(nev,nev),UL(nev,nev),UR(nev,nev),IPIV(nev))

      do i = 1, num_points
        ix(i) = i
      end do 
  
      do i = 0, nev-1
        call EPSGetEigenpair(eps,i,E_right(i+1),PETSC_NULL_SCALAR,Vi,PETSC_NULL_VEC,ierr);CHKERRA(ierr)
        call VecGetValues(Vi,num_points,ix,u_right(1:num_points,i+1),ierr);CHKERRA(ierr) 

        if (abs(u_right(num_points,i+1)) .gt. 0.1 ) then
          print*, u_right(num_points,i+1),l
        end if
 

        u_right(num_points,i+1) = (0d0,0d0)
        rint = sum(u_right(:,i+1)*u_right(:,i+1))
        u_right(:,i+1) = u_right(:,i+1)/zsqrt(rint)

        if (eps_two_sided .eqv. PETSC_TRUE) then
          call EPSGetLeftEigenvector(eps,i,Vi,PETSC_NULL_VEC,ierr);CHKERRA(ierr)
          call VecGetValues(Vi,num_points,ix,u_left(1:num_points,i+1),ierr);CHKERRA(ierr) 
          u_left(num_points,i+1) = (0d0,0d0)
          rint  = sum(u_left(:,i+1)*u_left(:,i+1))
          u_left(:,i+1) = dconjg(u_left(:,i+1)/zsqrt(rint))

          rint  = sum(u_left(:,i+1)*u_right(:,i+1))
          u_left(:,i+1) = u_left(:,i+1)/zsqrt(rint)
          u_right(:,i+1) = u_right(:,i+1)/zsqrt(rint)
        else 
          do j = 0,i-1
            u_right(:,i+1) = u_right(:,i+1) - sum(u_right(:,i+1)*u_right(:,j+1))*u_right(:,j+1)
          end do
          rint = sum(u_right(:,i+1)*u_right(:,i+1))
          u_right(:,i+1) = u_right(:,i+1)/zsqrt(rint)
        end if
      end do 

      if(eps_two_sided .eqv. PETSC_TRUE ) then
        S = matmul(transpose(u_left),u_right)

        call zgetrf(nev,nev,S,nev,IPIV,INFO)

        if (INFO .ne. 0) then
          print*, "'LU' factorization failed"
          print*, " exit code ", INFO
        end if 

        do i = 1, nev
          if (IPIV(i) .ne. i) then
            print*,'Pivits are hapening',l
            !stop
          end if 
        end do 

        UR = 0.d0
        UL = 0.d0
        do i = 1,nev
          do j = i,nev
            UR(i,j) = S(i,j)
          end do
        end do 
        
        do j = 1, nev
          do i = j+1,nev
            UL(i,j) = S(i,j)
          end do
        end do 

        do i = 1,nev
          UL(i,i) = 1.d0
        end do 

        call ztrtri('U','N',nev,UR,nev,INFO)
        call ztrtri('L','N',nev,UL,nev,INFO)
        
        if (INFO .ne. 0) then
          print*, 'Inversion failed'
          stop
        end if 
  

        u_left = matmul(u_left,transpose(UL))

        u_right = matmul(u_right,UR)
      else
        u_left = u_right
      end if 

      ! There will be nmax-l eigenstates for each l 
      ener_dims(1) = nmax-l
      ener_dims(2) = 2
      ! Each of the nmax-l eigenstates will have a wfn of length int(Rmax/h)
      psi_dims(1) = int(Rmax/dr)
      psi_dims(2) = nmax-l
      psi_dims(3) = 2
      
      uu_right(:,:,1) = real(real(u_right))
      uu_right(:,:,2) = real(aimag(u_right))
  
      EE_right(:,1) = real(real(E_right))
      EE_right(:,2) = real(aimag(E_right))
  
      ! Writes E_right to Energy_l#l (l#l can be l0, l1, ... etc ) 
      call h5dwrite_f( ener_dset(l), h5_kind, EE_right, ener_dims, h5_err)
      
      ! Writes u_right to Psi_l#l
      call h5dwrite_f( psi_dset_right(l), h5_kind, uu_right, psi_dims, h5_err)


      uu_left(:,:,1) = real(real(u_left))
      uu_left(:,:,2) = real(aimag(u_left))
  
      ! Writes u_left to Psi_l_l#l
      call h5dwrite_f( psi_dset_left(l), h5_kind, uu_left, psi_dims, h5_err)
  
      deallocate(uu_right,EE_right,ix)
      deallocate(u_left,u_right,E_right,uu_left)
      deallocate(S,UL,UR,IPIV)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !     Display solution and clean up
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call EPSDestroy(eps,ierr);CHKERRA(ierr)
      call MatDestroy(H,ierr);CHKERRA(ierr)
      call MatDestroy(H,ierr);CHKERRA(ierr)

      deallocate(V)
    end do 
  
     ! For each l close the hdf5 objects that were created in initialize
    do l = 0, lmax
      ! Closes the energy dataset
      call h5dclose_f( ener_dset(l), h5_err)
      
      ! Closes the energy dataspace
      call h5sclose_f( ener_space(l), h5_err)
  
      ! Closes the psi dataset
      call h5dclose_f( psi_dset_right(l), h5_err)
      
      ! Closes the psi dataspace
      call h5sclose_f( psi_space_right(l), h5_err)
  
      ! Closes the psi dataset
      call h5dclose_f( psi_dset_left(l), h5_err)
      
      ! Closes the psi dataspace
      call h5sclose_f( psi_space_left(l), h5_err)

    end do 
    
    ! Closes the hdf5 file
    call h5fclose_f( file_id, h5_err)
  
    call SlepcFinalize(ierr)

    call CPU_TIME(end_time)

    print 20, 'time   :', end_time-start_time
    20 format(A8,ES9.2)

  end program main 
  
  