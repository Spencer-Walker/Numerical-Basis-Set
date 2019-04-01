!------------------------------------------------------------------------------
! CU Boulder, Jila
!------------------------------------------------------------------------------
!
! PROGRAM:  schrodinger1Df90.F90
!
!> Spencer.Walker@colorado.edu
!> Spencer Walker
!
! DESCRIPTION: 
!>  This program generates the basis for numerov method and writes it to an .h5
!   file.
!
! REVISION HISTORY:
! 22 01 2019 - Initial Version
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Gram_schmidt
!------------------------------------------------------------------------------
subroutine Gram_schmidt( n, l, nmax, length, u)
  ! This subroutine orthogonalizes the set of eigenfunctions produced
  ! from this program.
  implicit none
  integer, parameter :: dp = kind(1.d0)! double precision
  ! Input
  integer,  intent(in) :: n,l,nmax,length
  ! Output
  real(dp), dimension(length,nmax-l), intent(out) :: u
  ! Derived
  real(dp)  :: dot
  integer   :: i

  ! Orthogonalize the new basis vector for the nth state
  do i = l+1,n-1
    dot =  dot_product(u(:,n-l),u(:,i-l))
    u(:,n-l) = u(:,n-l) - dot*u(:,i-l)

  end do

  u(:,n-l) = u(:,n-l)/norm2(u(:,n-l))

end subroutine Gram_schmidt

!------------------------------------------------------------------------------
! Init_mesh
!------------------------------------------------------------------------------
subroutine Init_mesh( num_points, R0, r)
  ! This subroutine initializes a uniformly spaced mesh with num_points
  ! between 0 and R0.
  implicit none
  integer, parameter     :: dp = kind(1.d0)! double precision
  ! Input
  integer,  intent(in)   :: num_points
  real(dp), intent(in)   :: R0
  ! Output
  real(dp), dimension(0:num_points+1), intent(out) :: r
  ! Derived
  real(dp)   :: h
  integer    :: i

  ! Grid Spacing
  h = R0/dble(num_points+1)
  ! Forms mesh
  do i = 0,num_points+1
    r(i) = dble(i)*h

  end do

end subroutine Init_mesh

!------------------------------------------------------------------------------
! Init_pot
!------------------------------------------------------------------------------
subroutine Init_pot( num_points, l, Znuc, r, V)
  ! This subroutine initializes the effective potential for the radial
  ! schrodinger equation for hydrogen.
  use simulation_parametersf90
  implicit none
  ! Input
  integer,  intent(in)  :: num_points, l, Znuc
  real(dp), dimension(0:num_points+1), intent(in)  :: r
  ! Output
  real(dp), dimension(0:num_points+1), intent(out) :: V
  ! Derived
  integer :: i

  ! Sets the potential
  if(l /= 0) then
    do i = 1,num_points+1
      V(i) = -sae_c0/r(i) - (sae_Zc*dexp(-sae_c*r(i)))/r(i) - &
      sae_a1*dexp(-sae_b1*r(i)) - sae_a2*dexp(-sae_b2*r(i)) - &
      sae_a3*dexp(-sae_b3*r(i)) - sae_a4*dexp(-sae_b4*r(i)) - &
      sae_a5*dexp(-sae_b5*r(i)) + 0.5d0*dble(l)*(dble(l) + 1.d0)/r(i)**2.d0 

    end do
  else
    do i = 1, num_points + 1
      V(i) = -sae_c0/r(i) - (sae_Zc*dexp(-sae_c*r(i)))/r(i) - &
      sae_a1*dexp(-sae_b1*r(i)) - sae_a2*dexp(-sae_b2*r(i)) - &
      sae_a3*dexp(-sae_b3*r(i)) - sae_a4*dexp(-sae_b4*r(i)) - &
      sae_a5*dexp(-sae_b5*r(i))
    end do
  end if

  ! Extend the potential
  V(0) = 0

end subroutine Init_pot

!------------------------------------------------------------------------------
! Turning_pt
!------------------------------------------------------------------------------
subroutine Turning_pt( num_points, r, l, E, Rm, im)
  ! Finds the classical turning point 
  ! This point will be used to match the left and right wfns 
  use simulation_parametersf90
  implicit none
  ! Input
  integer,  intent(in)  :: num_points, l
  real(dp), intent(in)  :: E
  real(dp), dimension(0:num_points+1), intent(in) :: r
  ! Output
  real(dp), intent(out) :: Rm
  integer,  intent(out) :: im
  ! Derived 
  real(dp) :: desc,denom

  ! Finds the turning point
  denom = (sae_a1 + sae_a2 + sae_a3 + sae_a4 + sae_a5 + 2d0*E - 2d0*sae_c*sae_Zc)
  
  desc = (sae_c0 + sae_Zc)**2d0 + dble(l*(1 + l))*denom
  
  if (desc >= 0) then
    Rm  = max(-((sae_c0 + sae_Zc - Sqrt(desc))/denom),-((sae_c0 + sae_Zc + Sqrt(desc))/denom))
  else
    Rm = r(num_points+1)/2
  end if

  ! Find index of midpoint
  im = 1
  do while ((r(im)<Rm).and.(im<num_points))
    im = im + 1

  end do

end subroutine Turning_pt

!------------------------------------------------------------------------------
! Get_error
!------------------------------------------------------------------------------
subroutine Get_error( num_points, R0, im, yl, yr, error)
  !This function computes the logarithmic error between the two wfns
  !error = ( d (log u) - d (log v) )/h
  implicit none
  integer, parameter :: dp = kind(1.d0)! double precision
  ! Input
  integer,  intent(in)  :: num_points, im
  real(dp), intent(in)  :: R0
  real(dp), dimension(0:num_points+1),intent(in) :: yl, yr
  ! Output
  real(dp),  intent(out) :: error
  ! Derived
  real(dp) :: h, Dy

  ! Computes the distance between neighboring pts 
  h = R0/dble(num_points+1)

  ! Compute log error
  error =(Dy(yl,im,num_points)/yl(im)-Dy(yr,im,num_points)/yr(im))/h

end subroutine Get_error

!------------------------------------------------------------------------------
! Dy
!------------------------------------------------------------------------------
function Dy( y, im, num_points)
  ! Finds the difference of y between points near im 
  implicit none
  integer, parameter :: dp = kind(1.d0)! double precision
  ! Input
  integer,  intent(in) :: num_points, im
  real(dp), dimension(0:num_points+1), intent(in) :: y
  ! Output
  real(dp) :: Dy

  ! Computes the difference
  Dy = (y(im+1)-y(im-1))/2d0
end function Dy

!------------------------------------------------------------------------------
! Match_and_normalize
!------------------------------------------------------------------------------
subroutine Match_and_normalize( num_points, im, yl, yr, y)
  ! Matches the left and right wfns and normalizes the overall wfn y
  implicit none
  integer, parameter :: dp = kind(1.d0)! double precision
  ! Input
  integer,  intent(in) :: im, num_points
  real(dp), dimension(0:num_points+1), intent(in)  :: yr
  ! Output
  real(dp), dimension(0:num_points+1), intent(out) :: y
  ! In+Out
  real(dp), dimension(0:num_points+1) :: yl
  ! Derived
  real(dp)  :: norm

  ! Match Left and Right wfn
  yl = yr(im)*yl/yl(im)

  ! Form the combined wavefunction
  y(0:im)              = yl(0:im)
  y(im+1:num_points+1) = yr(im+1:num_points+1)

  ! Now we normalize
  norm = norm2(y)
  y = y/norm

end subroutine Match_and_normalize

!------------------------------------------------------------------------------
! Itterate
!------------------------------------------------------------------------------
subroutine Itterate(l, num_points, R0, Rm, im, E, r, V, yl, yr, &
  & num_nodes, error)
  ! This subroutine (1) finds the classical (call) turning pt for given
  ! E. (2) integrates the TISE from 0 to the turning pt (yl) and from R
  ! to the turning pt (yr) and (3) finds the error btwn yl and yr and
  ! the turning pt.
  implicit none
  integer, parameter :: dp = kind(1.d0)! double precision
  ! Input
  integer,  intent(in)  :: l, num_points
  real(dp), intent(in)  :: E, R0
  real(dp), dimension(0:num_points+1), intent(in) :: r, V
  ! Output
  real(dp), intent(out) :: error
  real(dp), dimension(0:num_points+1), intent(out) ::yl, yr
  integer,  intent(out) :: num_nodes
  ! Derived
  integer  :: im, i
  real(dp) :: Rm, h
  real(dp), dimension(0:num_points+1) :: Q
  ! External
  external :: Get_error

  ! Computes the distance between nodes
  h = R0/dble(num_points+1)

  ! Computes the (outer) classical turning point 
  ! This is where the left and right wfns will be matched
  call Turning_pt(num_points,r,l,E,Rm,im)

  ! Apply y(0) = 0 boundary condition
  yl(0) = 0.d0

  ! Normalize the arbituary y(1) value st things don't blow up
  if(abs(yl(im)) < 1d10) then
    yl(1) = 1.d-32
  else
    yl(1) = yl(1)/yl(im)
  end if

  ! Apply y(end) = 0 boundary condition for excited states
  yr(num_points+1) = 0

  ! Normalize things appropriately
  if(abs(yr(im))< 1.d10) then
    yr(num_points)   = 1.d-16
  else
    yr(num_points)   = yr(num_points)/yr(im)
  end if

  ! If we're searching for a bound state use the asymtotic solution
  if(E<0) then
    yr(num_points+1)= yr(num_points)*exp(-h*sqrt(-2.d0*E*h))
  end if

  ! Define a RHS operator for the problem
  Q = (/( 2.d0*( E - V(i) ), i = 0,num_points + 1  )/)
  
  ! Count the nodes as we integrate in space
  num_nodes = 0

  10 continue
  
  ! Here I use a cusp correction for points near r = 0
  ! Look up "Numerov method for singular potentials"
  ! A few papers should come up.
  yl(2) = (-24.d0 + h**4.d0*Q(1)*Q(3) + h**2.d0*( &
  & 13.d0*Q(1) -3.d0*Q(3) ))*yl(1)/(-12.d0 + h**2.d0*(&
  & 2.d0*Q(2) - 3.d0*Q(3) ) + h**4.d0*Q(2)*Q(3))

  yl(3) = (-36.d0 -11.d0*h**4.d0*Q(1)*Q(2) + 9.d0*h**2.d0*( &
  & 3.d0*Q(1) + 2.d0*Q(2) ))*yl(1)/(-12.d0 + h**2.d0*(&
  & 2.d0*Q(2) - 3.d0*Q(3) ) + h**4.d0*Q(2)*Q(3))

  ! Integrate the rest of the values with the standard Numerov method
  do i = 3,im
    ! Finds the next value of yl
    yl(i+1) = ( 2d0*( 1.d0 - 5d0*Q(i)*h**2d0/12d0 )*yl(i)-&
    & ( 1.d0 + Q(i-1)*h**2d0/12d0)*yl(i-1))/&
    & (1.d0+Q(i+1)*h**2d0/12d0)

    ! Counts the number of nodes
    if((yl(i+1)>0.and.(yl(i)<0.and.(yl(i-1)<0.and.(yl(i-2)<0.and.&
    & (yl(i-3)<0))))).or.&
    & (yl(i+1)<0.and.(yl(i)>0.and.(yl(i-1)>0.and.(yl(i-2)>0.and.&
    & (yl(i-3)>0)))))) then
      if(i>3.and.i<im) then
        num_nodes = num_nodes+1
      end if
    end if

  end do

  ! Generate values for the right wfn from i = num_pts,im
  do i = num_points,im,-1
    ! Finds the next value of yr
    yr(i-1) = ( 2d0*( 1.d0 - 5d0*Q(i)*h**2d0/12d0 )*yr(i)-&
    & (1.d0+Q(i+1)*h**2d0/12d0)*yr(i+1))/&
    & (1.d0+Q(i-1)*h**2d0/12d0)
    ! Counts the number of nodes
    if((yr(i-1)>0.and.(yr(i)<0.and.(yr(i+1)<0.and.(yr(i+2)<0.and.&
    & (yr(i+3)<0))))).or.&
    & (yr(i-1)<0.and.(yr(i)>0.and.(yr(i+1)>0.and.(yr(i+2)>0.and.&
    & (yr(i+3)>0)))))) then
      if(i<num_points-2 .and. i>im) then
        num_nodes = num_nodes+1
      end if
    end if
    if (isnan(yl(im+1))) then
      yl(1) = yl(1)*1d-16
      goto 10
    end if 
    if (isnan(yr(im-1))) then
      yr(num_nodes) = yr(num_nodes)*1d-16
      goto 10
    end if
   
  end do

  ! Gets the log error between the left and right wfns 
  call Get_error(num_points,R0,im,yl,yr,error)

end subroutine Itterate

!------------------------------------------------------------------------------
! Get_bounds
!------------------------------------------------------------------------------
subroutine Get_bounds( n, l, num_points, R0, Rm, im, r, V, yl, yr, dE, &
  & Emin, Emax, min, max)
  ! This subroutine finds an upper and lower bound for possible energies
  ! which yield the correct number of nodes n-l-1 for the wfns
  use simulation_parametersf90
  implicit none
  ! Input
  integer,   intent(in)  :: n, l, num_points, im
  real(dp),  intent(in)  :: R0, Rm, dE
  real(dp),  dimension(0:num_points+1) ::r, V, yl, yr
  ! Output
  real(dp),  intent(out) :: Emin, Emax, min, max
  ! Derived
  integer  :: nodes, num_nodes
  real(dp) :: Eold
  ! Code
  nodes = n-l-1
  ! Finds an appropriate bound for Emin
  call Itterate( l, num_points, R0, Rm, im, Emin, r, V, yl, yr, num_nodes,&
  & min)

  ! Increase Emin till an appropriate number of nodes are formed
  do while(num_nodes<=nodes)
    Emin = Emin+sqrt(dE)
    call Itterate( l, num_points, R0, Rm, im, Emin, r, V, yl, yr, &
    & num_nodes, min)

  end do

  ! Decrease Emin till an appropriate number of nodes are formed
  do while(num_nodes>nodes)
    Emin = Emin-sqrt(dE)
    call Itterate( l, num_points, R0, Rm, im, Emin, r, V, yl, yr, &
    & num_nodes, min)

  end do

  do while(num_nodes>nodes)
    Emin = Emin-dE
    call Itterate( l, num_points, R0, Rm, im, Emin, r, V, yl, yr, &
    & num_nodes, min)

  end do

  ! Decrease Emin until the bottom of the appropiate doMain is found
  if (Emin <E_lower-dE) then
    Emin = E_lower -dE

  end if

  Eold = Emin
  do while(num_nodes == nodes .and. Emin>-0.51d0)
    Eold = Emin
    Emin = Emin-dE

    call Itterate( l, num_points, R0, Rm, im, Emin, r, V, yl, yr, &
    & num_nodes, min)

  end do

  ! Verify that the condition is met
  Emin = Eold
  call Itterate( l, num_points, R0, Rm, im, Emin, r, V, yl, yr, num_nodes,&
  & min)

  ! Finds the appropriate bound for Emax
  Emax = Emin+dE
  call Itterate( l, num_points, R0, Rm, im, Emax, r, V, yl, yr, num_nodes,&
  & max)

  ! Increase Emax until the node constraint is satsified
  do while(num_nodes<nodes+1)
    Emax = Emax+sqrt(dE)
    call Itterate( l, num_points, R0, Rm, im, Emax, r, V, yl, yr, &
    & num_nodes, max)

  end do

  ! Decrease Emax until the node constraint is satsified
  do while(num_nodes>nodes)
    Emax = Emax-dE
    call Itterate( l, num_points, R0, Rm, im, Emax, r, V, yl, yr, &
    & num_nodes, max)

  end do

  ! Increase Emax to the top of the appropriate doMain
  Eold = Emax
  do while(num_nodes == nodes .and. Emax < 15)
    Eold = Emax
    Emax = Emax+dE
    call Itterate( l, num_points, R0, Rm, im, Emax, r, V, yl, yr, &
    & num_nodes, max)

  end do

  ! Set Emax to the top most energy that still satisfies the nodes
  Emax = Eold

  ! Generate left and right wfns for Emax = Eold 
  call Itterate( l, num_points, R0, Rm, im, Emax, r, V, yl, yr, num_nodes,&
  & max)

end subroutine Get_bounds

!------------------------------------------------------------------------------
! Get_correction
!------------------------------------------------------------------------------
subroutine Get_correction( yl, yr, num_points, im, R0, E1)
  ! We compute energy corrections by taking variations of the
  ! logarithm of the left and right wfns and matching them at
  ! the turning point.(See The Calculation of Atomic Structures by Hartree)
  implicit none
  integer, parameter :: dp = kind(1.d0)! double precision
  ! Input
  real(dp)  :: R0
  real(dp),   dimension(0:num_points+1), intent(in) :: yl, yr
  integer,    intent(in) :: num_points, im
  ! Output
  real(dp),   intent(out) :: E1
  ! Derived
  real(dp)                :: Il, Ir, dvl, dvr, Dy, h
  real(dp),   dimension(0:num_points+1) :: vl, vr

  ! Get gridspacing
  h = R0/dble(num_points+1)

  ! Here I normalize the data s.t. things don't blow up
  vr = yr/yr(im)
  vl = yl/yl(im)

  ! Integrate the square of the left wfn
  Il = sum(vl(1:im)**2)*h
  ! Integrate the square of the right wfn
  Ir = sum(vr(im:num_points)**2)*h

  ! Normalized derivatives
  dvl = Dy(vl,im,num_points)
  dvr = Dy(vr,im,num_points)

  ! Compute the energy update
  ! (given in Hartree's atomic structure text )
  E1 = -((dvr - dvl)/(Ir + Il))/2.d0
  if (isnan(E1)) then
    E1 = 0.d0
  end if 
end subroutine Get_correction

!------------------------------------------------------------------------------
! Refine
!------------------------------------------------------------------------------
subroutine Refine( n, l, num_points, R0, r, V, E, tol, y, error)
  ! Once the energy is in some neighborhood dE of the "exact" energy
  ! we are allowed to treat the left and right wfns as variations of
  ! eachother, and are able to iteratively compute better energy estimates
  implicit none
  integer, parameter :: dp = kind(1.d0)! double precision
  ! Input
  integer,  intent(in)  :: n, l, num_points
  real(dp), intent(in)  :: R0, tol
  real(dp), dimension(0:num_points+1), intent(in) :: r, V
  ! Outout
  real(dp),   intent(out) :: E, error
  real(dp),   dimension(0:num_points+1) :: y
  ! Derived
  real(dp) :: Rm,E1
  real(dp),   dimension(0:num_points+1) :: yl, yr
  integer  :: im, nodes, num_nodes
  ! External
  external  :: Get_correction, Match_and_normalize

  ! Compute the number of nodes the wfn must have
  nodes = n-l-1

  ! Initialize the energy correction to start at zero
  E1 = 0

  ! This loop computes itterative updates to the wfn using the method 
  ! outlined in Hartree's text on atomic structure 
  REFINEMENT_LOOP: do
    ! Update the old energy
    E = E+E1
    if(isnan(E)) then
      print*,'REFINEMENT ERROR: E1 ISNAN'
      print*,'im',im,'num_points',num_points
      stop
    end if 
    ! Compute the left and right wfns using E
    call Itterate( l, num_points, R0, Rm, im, E, r, V, yl, yr, num_nodes,&
    & error)

    ! Get the energy correction from the left and right wfns
    call Get_correction(yl,yr,num_points,im,R0,E1)

    ! If the energy update is under some tolerance or the produced wfn
    ! does not have the correct number of nodes exit the loop
    if (abs(E1) < tol .or. nodes /= num_nodes) exit REFINEMENT_LOOP

  end do REFINEMENT_LOOP

  ! Take the final left and right wfns, match them and normalize the wfn
  call Match_and_normalize(num_points,im,yl,yr,y)

end subroutine Refine

!------------------------------------------------------------------------------
! Search
!------------------------------------------------------------------------------
subroutine Search( n, l, Znuc, num_points, R0, Emin, Emax, dE, tol, y, E)
  ! We search for the radial wavefunctions and energy by first using
  ! bisection to get a rough estimate and then use variations of the
  ! left and right wfns to compute additional refinements
  implicit none
  integer, parameter :: dp = kind(1.d0)! double precision
  ! Input
  integer,    intent(in)  :: n, l, num_points, Znuc
  real(dp),   intent(in)  :: dE, R0, tol
  ! Outout
  real(dp),   intent(out) :: E
  real(dp),   intent(out) :: y(0:num_points+1)
  ! Derived
  real(dp) :: Rm, error, start, finish, range, left, right, Emin, Emax
  real(dp) :: min, max
  real(dp),   dimension(0:num_points+1) :: yl, yr, r, V
  integer  :: im, nodes, num_nodes
  ! External
  external :: Init_mesh, Init_pot, Itterate, Match_and_normalize
  
  ! Compute the number of nodes that the n,l wfn must have
  nodes = n-l-1

  ! Initialize the uniformally spaced mesh
  call Init_mesh(num_points,R0,r)

  ! Initialize the effective radial potential
  call Init_pot(num_points,l, Znuc, r,V)

  ! Find appropriate bounds for the energy for the binary Search
  call Get_bounds(n,l,num_points,R0,Rm,im,r,V,yl,yr,dE,Emin,Emax,min,max)

  ! We will use a binary search to find an energy estimate here we chose 
  ! the bounds + range + midpoint of the binary search method 
  start  = Emin
  finish = Emax
  range  = finish - start
  E      = (start + finish)/2

  ! We will keep bisecting some energy interval until the range used 
  ! for the search is under some tolerance
  do while( range > dE )
    ! Look at an energy higher than E 
    call Itterate( l, num_points, R0, Rm, im, E+dE, r, V, yl, yr, &
    & num_nodes, right)

    ! Look at an energy less than E
    call Itterate( l, num_points, R0, Rm, im, E-dE, r, V, yl, yr, &
    & num_nodes, left)

    ! If once Energy E+/-dE produces wfns that match better split the
    ! domain in half 
    if(abs(right)>abs(left)) then
      finish  = E
    else
      start = E
    end if

    ! Compute the new range
    range = finish - start

    ! Compute the new midpoint
    E = (start + finish)/2

  end do

  ! Refine the prior rough estimate for the wfn 
  call Refine(n,l,num_points,R0,r,V,E,tol,y,error)

end subroutine Search

!------------------------------------------------------------------------------
! Initialize
!------------------------------------------------------------------------------
subroutine Initialize( h, nmax, lmax, Rmax, label, file_id, ener_space, &
  & ener_dset, psi_space, psi_dset, h5_err, num_proc, proc_id,  mpi_err)
  ! Starts mpi + hdf5 and also makes the hdf5 file and structure that will 
  ! be used to output all of the data form this program 
  use hdf5
  use mpi
  implicit none 
  integer,  parameter  :: dp = kind(1.d0)! double precision double precision
  ! Input
  integer,  intent(in) :: nmax, lmax
  real(dp), intent(in) :: h, Rmax
  character(len = 15), intent(in) :: label ! file name w/o .h5
  ! Output
  integer(HID_T), dimension(0:lmax), intent(out) :: psi_space, psi_dset
  integer(HID_T), dimension(0:lmax), intent(out) :: ener_space, ener_dset
  integer(HID_T), intent(out) :: file_id
  integer, intent(out) :: mpi_err, h5_err, num_proc, proc_id
  ! Derived
  integer(HSIZE_T)  :: ener_dims(1:1)
  integer(HSIZE_T)  :: psi_dims(1:2) 
  integer(HID_T)    :: h5_kind, plist_id
  character(len = 24)   :: file_name 
  character(len = 3)    :: strl ! file number
  character(len = 6)    :: fmt  ! format descriptor
  integer :: comm, info, l 

  ! Add the .h5 extension to the file label
  file_name = trim(label)//'.h5'

  ! Abbreviate mpi comm+info
  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  ! Initialize mpi
  call MPI_Init( mpi_err)
  
  ! Find how many workers are available 
  call MPI_Comm_size( comm, num_proc, mpi_err)

  ! Figure out what the worker id is 
  call MPI_Comm_rank( comm, proc_id, mpi_err)

  ! Initialize hdf5 
  call h5open_f( h5_err)

  ! Creates the hdf5 property list that will be used for the file 
  ! (like an interface or sae_c template)
  call h5pcreate_f( H5P_FILE_ACCESS_F, plist_id, h5_err)

  ! Stores MPI IO communicator information to the file access property 
  ! list. 
  call h5pset_fapl_mpio_f( plist_id, comm, info, h5_err)

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

    ! There will be nmax-l wfns, each with dim int(Rmax/h)
    psi_dims(1) = int(Rmax/h)
    psi_dims(2) = nmax-l

    ! Creates an hdf5 interface that knows information about the data
    call h5screate_simple_f( 1, ener_dims, ener_space(l), h5_err)
    
    ! Connects the datasapce to the hdf5 file and gives it a name 
    call h5dcreate_f( file_id, "Energy_l"//trim(strl), h5_kind, ener_space(l), &
    & ener_dset(l), h5_err)

    ! Does the same as before for the energy but with the wfns 
    call h5screate_simple_f( 2, psi_dims, psi_space(l), h5_err)
    call h5dcreate_f( file_id, "Psi_l"//trim(strl), h5_kind, psi_space(l), &
    & psi_dset(l), h5_err)

  end do 

end subroutine Initialize

!------------------------------------------------------------------------------
! Finalize
!------------------------------------------------------------------------------
subroutine Finalize( lmax, file_id, ener_space, ener_dset, &
  & psi_space, psi_dset, h5_err, mpi_err)
  ! Closes all hdf5 and mpi objects 
  use mpi
  use hdf5
  implicit none
  ! Input
  integer, intent(in) :: lmax
  integer(HID_T), intent(in) :: file_id
  integer(HID_T), dimension(0:lmax), intent(in) :: ener_space, ener_dset
  integer(HID_T), dimension(0:lmax), intent(in) :: psi_space, psi_dset
  ! Output
  integer, intent(out) :: h5_err, mpi_err
  ! Derived 
  integer :: l

  ! For each l close the hdf5 objects that were created in initialize
  do l = 0, lmax
    ! Closes the energy dataset
    call h5dclose_f( ener_dset(l), h5_err)
    
    ! Closes the energy dataspace
    call h5sclose_f( ener_space(l), h5_err)

    ! Closes the psi dataset
    call h5dclose_f( psi_dset(l), h5_err)
    
    ! Closes the psi dataspace
    call h5sclose_f( psi_space(l), h5_err)
    
  end do 
  
  ! Closes the hdf5 file
  call h5fclose_f( file_id, h5_err)

  ! Closes all mpi 
  call MPI_Finalize( mpi_err)

end subroutine Finalize

!------------------------------------------------------------------------------
! Save_spectrum
!------------------------------------------------------------------------------
subroutine Save_spectrum(l, nmax, lmax, Rmax, h, u, E, ener_dset, psi_dset,&
  & h5_err)
  ! Saves all of the data for each angular momentum block 
  use mpi
  use hdf5
  implicit none 
  integer,  parameter   :: dp = kind(1.d0)! double precision
  ! Input
  integer,  intent(in)  :: l, nmax, lmax
  real(dp), intent(in)  :: Rmax, h
  real(dp), dimension(1:int(Rmax/h),1:nmax-l), intent(in)  :: u
  real(dp), dimension(1:nmax-l), intent(in) :: E
  integer(HID_T), dimension(0:lmax), intent(in) :: psi_dset, ener_dset
  ! Output
  integer, intent(out) :: h5_err
  ! Derived 
  integer(HID_T) :: h5_kind
  integer(HSIZE_T)    :: ener_dims(1:1)
  integer(HSIZE_T)    :: psi_dims(1:2) 

  ! Lets the user know that l is being written to lable.h5
  10 format(A8,I4)
  print 10, 'writing ', l

  ! Converts fortran dp kind to mpi real type
  h5_kind = h5kind_to_type( dp, H5_REAL_KIND)

  ! There will be nmax-l eigenstates for each l 
  ener_dims(1) = nmax-l

  ! Each of the nmax-l eigenstates will have a wfn of length int(Rmax/h)
  psi_dims(1) = int(Rmax/h)
  psi_dims(2) = nmax-l
  
  ! Writes E to Energy_l#l (l#l can be l0, l1, ... etc ) 
  call h5dwrite_f( ener_dset(l), h5_kind, E, ener_dims, h5_err)
  
  ! Writes u to Psi_l#l
  call h5dwrite_f( psi_dset(l), h5_kind, u, psi_dims, h5_err)

end subroutine Save_spectrum

!------------------------------------------------------------------------------
! Process_vec
!------------------------------------------------------------------------------
subroutine Process_vec( n, l, nmax, num_points, Rmax, h, y, u)
  ! As each new energy En is computed the wfn yn is orthogonalized and 
  ! added to the matrix u 
  implicit none 
  integer,  parameter     :: dp = kind(1.d0)! double precision
  ! Input
  integer,  intent(in)    :: n, nmax, num_points, l
  real(dp), intent(in)    :: h, Rmax
  real(dp), dimension(0:num_points+1), intent(in) :: y
  ! Input+Output
  real(dp), dimension(int(Rmax/h),nmax-l), intent(inout) :: u
  ! Derived 
  integer                 :: i
  
  ! Store the final eigenfunction
  u(1:num_points + 1,n - l) = y(1:num_points + 1)

  ! Make the leftover space 0
  if ( num_points+1 < int(Rmax/h) ) then
    u(num_points + 2:int(Rmax/h),n - l) = &
    & (/( 0,i = num_points + 2, int(Rmax/h) )/)
  end if

  ! Orthogonalize the basis every time a new vector is added.
  call Gram_schmidt( n, l, nmax, int(Rmax/h), u)

end subroutine Process_vec

!------------------------------------------------------------------------------
! Process_vec
!------------------------------------------------------------------------------
subroutine Get_basis( h, Rmin, nmax, lmax, Rmax, dE, tol, Emax, Emin, Znuc, &
  & proc_id, num_proc, ener_dset, psi_dset, h5_err)
  use hdf5
  use mpi
  implicit none 
  integer, parameter   :: dp = kind(1.d0)! double precision
  ! Input
  integer,  intent(in) :: nmax, lmax, Znuc, proc_id, num_proc
  real(dp), intent(in) :: h, Rmin, Rmax, dE, tol
  integer(HID_T), dimension(0:lmax) :: ener_dset, psi_dset
  ! Output
  integer :: h5_err
  ! Input + Output
  real(dp), intent(inout) :: Emax, Emin
  ! Derived
  integer               :: n, num_points, l, stride, mpi_err
  real(dp)              :: R0, Eold
  real(dp), allocatable :: u(:,:), E(:), y(:)
  external       :: search

  ! Creates formats for the n,l output table
  10 format(A47)
  20 format(A3,I4,A3,I4,A3,ES11.4,A7,ES11.4,A1)

  ! Top of n,l table
  if (proc_id == 0) then
    print 10,'-----------------------------------------------'      
  end if 
  call MPI_Barrier( MPI_COMM_WORLD, mpi_err)

  if( num_proc < lmax+1) then
    stride = ceiling( real(lmax+1)/num_proc)
  else 
    stride = 1
  end if 

  do l = proc_id*stride, min( (proc_id+1)*stride-1,lmax)
    ! Alocate space for the wavefunction
    allocate( u(int(Rmax/h), nmax-l))
    ! Allocate space for the energies
    allocate( E(nmax-l))

    ! Iterate through allowed quantum numbers
    do n = l+1, nmax
      ! The num_proc of the grid increases for every energy level
      R0 = Rmin*n**2.d0
      ! Truncate the doMain at some R0
      if ( R0 > Rmax) then
        R0 = Rmax
      end if

      ! Calculate the number of points from h and R0
      num_points = int(R0/h)-1;

      ! Allocate memory for a wavefunction in the new grid space
      allocate( y(0:num_points+1))

      ! Search for the correct n,l eigenfunction
      call Search( n, l, Znuc, num_points, R0, Emin, Emax, dE, tol, y, &
      & E(n-l))
      
      ! Print information about each n,l computed 
      print 20, '|n|', n, '|l|', l, '|E|', E(n-l), '|ΔE/E|', &
      & (E(n-l)+ 0.5d0/dble(n)**2.d0)/E(n-l), '|'
      
      ! Use old calculation for new energy estimates
      Eold = Emin
      Emin = Emax + dE
      Emax = Emax + ( Emax-Eold)

      ! Orthogonalize y and add it to u 
      call Process_vec( n, l, nmax, num_points, Rmax, h, y, u)

      ! Remove y from memory
      deallocate( y)

    end do
  
    ! Save the l block to our label.h5 file 
    call Save_spectrum( l, nmax, lmax, Rmax, h, u, E, ener_dset, &
    & psi_dset, h5_err)
    
    ! Remove u and E from memory
    deallocate( u)
    deallocate( E)

  end do 
  call MPI_Barrier( MPI_COMM_WORLD, mpi_err)

  ! Bottom of n,l table 
  if (proc_id == 0) then
    print 10,'-----------------------------------------------'     
  end if 
end subroutine Get_basis

!------------------------------------------------------------------------------
! Main
!------------------------------------------------------------------------------
program Main
  use mpi
  use simulation_parametersf90
  use hdf5
  implicit none
  ! Declarations
  integer             :: lmax, nmax, Znuc
  real(dp)            :: Rmin, dE, Emin, Emax, tol, h, Rmax
  integer :: num_proc, mpi_err, proc_id
  character(len = 15) :: label
  integer(HID_T) :: file_id
  integer(HID_T), allocatable ::  ener_space(:), ener_dset(:)
  integer(HID_T), allocatable ::  psi_space(:), psi_dset(:)
  integer        :: h5_err
  external       :: search
  real(dp)       :: start_time, end_time 


  ! Start CPU time 
  call CPU_TIME(start_time)

  ! Parameters are set in simulation_parametersf90
  h     = grid_space
  Rmin  = R_min
  nmax  = n_max
  Rmax  = R_max
  dE    = binary_search_tol
  tol   = refinement_tol
  lmax  = l_max
  Emax  = E_upper
  Emin  = E_lower
  label = hdf5_file_label

  ! If lmax is greater than nmax the program stops since these parameters 
  ! do not make any sense 
  if (lmax > nmax) then
    ! Prints an error message before stopping the program 
    print*,'ERROR: nmax must be greater than lmax'
    stop 
  end if

  ! Fomats for the terminal output 
  10 format(A47)
  20 format(A8,I4)
  30 format(A8,ES9.2)

  ! Allocates memory for bothe energy and wavefunction space + dset 
  allocate( ener_space(0:lmax), ener_dset(0:lmax))
  allocate( psi_space(0:lmax), psi_dset(0:lmax))

  ! Initializes mpi and hdf5 and sets up the file structure 
  call Initialize( h, nmax, lmax, Rmax, label, file_id, &
  & ener_space, ener_dset, psi_space, psi_dset, h5_err, num_proc, &
  & proc_id,  mpi_err)

  ! Only processessor 0 will print this so duplicates are not created
  if ( proc_id == 0 ) then
    print 10,'-----------------------------------------------'
    print 20, 'num_proc', num_proc
    print 10,'-----------------------------------------------'
    if (num_proc > lmax + 1) then
      print*, "Warning: Each l is assigned a processor having more &
      &than lmax+1 workers will not speed up code!"
    end if
  end if

  ! I stop all workers until all hdf5 objects are created this probably 
  ! isnt needed but is safer regardless, and garentees that the print 
  ! statement above is printed first 
  call MPI_Barrier( MPI_COMM_WORLD, mpi_err)

  ! Prints that proc_id #proc_id has been initialized properly
  print 20, 'start   ', proc_id

  ! Generates the basis + energy for this system and writes it to 
  ! label.h5
  call Get_basis( h, Rmin, nmax, lmax, Rmax, dE, tol, Emax, Emin, Znuc, &
  & proc_id, num_proc, ener_dset, psi_dset, h5_err)

  ! Prints when proc_id #proc_id has finished the above step 
  ! this lets the user know whether or not any processor has gotten 
  ! stuck in the prior step. If it has gotten stuck 9/10 times it 
  ! happens during the get bounds step.
  print 20, 'end     ', proc_id

  ! Closes all of the mpi and hdf5 structures 
  call Finalize( lmax, file_id, ener_space, ener_dset, psi_space, &
  & psi_dset, h5_err, mpi_err)

  ! Deallocates all remaining varaibles in the program 
  deallocate( ener_space, ener_dset, psi_space, psi_dset)

  call CPU_TIME(end_time)
  
  print 30, 'time   :', end_time-start_time

end program Main