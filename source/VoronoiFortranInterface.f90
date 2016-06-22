!==============================================================================
! Module VoronoiFortranInterface
! D. A. Hubber - 27/08/2013
! Fortran interface for Voronoi tessellation generating functions.
! Currently used as a front-end for qhull and Voro++ libraries.
!==============================================================================
MODULE VoronoiFortranInterface
  use iso_c_binding
  use common_mod  !, only : Connection, Vertex
  implicit none



  ! Precision of float variables
  !----------------------------------------------------------------------------
  !integer,parameter :: ILP = selected_int_kind(r=15) ! Integer long precision
  integer,parameter :: DP = selected_real_kind(p=15) ! Double precision
  integer,parameter :: QP = selected_real_kind(p=33) ! Quadruple precision
  integer,parameter :: SP = selected_real_kind(p=6)  ! Single precision
#if defined(QUADRUPLE_PRECISION)
  !integer,parameter :: PR = QP
#elif defined(DOUBLE_PRECISION)
  !integer,parameter :: PR = DP
#else
  !integer,parameter :: PR = SP
#endif


  ! Dimensionality of simulation
  !----------------------------------------------------------------------------
#if NDIM==1
  !integer,parameter :: ndim = 1
#elif NDIM==2
  !integer,parameter :: ndim = 2
#elif NDIM==3
  !integer,parameter :: ndim = 3
#endif
  integer,parameter :: ndim = 3


  ! Other constants and parameters
  !----------------------------------------------------------------------------
  !integer :: Npackets
  !real(kind=PR),parameter :: opacity = 2.0_PR
  !real(kind=PR),parameter :: pi = 3.14159265358979_PR
  !real(kind=PR),parameter :: Lstar = 1.0_PR


  !integer,parameter :: Nconnectmax = 32         ! Max. no. of connections
                                                ! per vertex
  logical :: allocated_voronoi                  ! Is Voronoi memory allocated?
  integer :: Nvertex                            ! No. of points
  integer :: Nvertexmax                         ! Max no. of points
  integer :: Nline                              ! No. of lines/connections
  integer :: Nlinemax                           ! Max. no. of lines/connections
  real(c_double) :: c_boxsize                   ! ..

  !type Connection                               ! Connection data type
  !   integer :: i1                              ! id of point 1
  !   integer :: i2                              ! id of point 2
  !   real(kind=4) :: length                     ! Length of connection
  !   real(kind=4) :: invlength                  ! 1 / length
  !   real(kind=4) :: runit(1:3)                 ! Unit vector (from i1 to i2)
  !end type Connection
  type(Connection),allocatable :: lines(:)      ! Array of lines/connections

  !type Vertex                                   ! Vertex data type
  !   real(kind=4) :: r(1:3)                     ! Position of vertex
  !   real(kind=4) :: mass                       ! Mass contained in cell
  !   real(kind=4) :: volume                     ! Voronoi cell volume
  !   real(kind=4) :: density                    ! Density of gas in cell
  !   real(kind=4) :: intensity                  ! Mean intensity
  !   integer :: Nconnect                        ! No. of connections
  !   integer :: iconnect(1:Nconnectmax)         ! ids of connecting particles
  !   type(Connection) :: connect(1:Nconnectmax) ! Connection data
  !end type Vertex
  !type(Vertex),allocatable :: vertices(:)       ! Array of verticies


CONTAINS



  !============================================================================
  ! InitialiseVoronoiVariables
  ! Initialise all variables related to generating the Voronoi tessellation.
  !============================================================================
  SUBROUTINE InitialiseVoronoiVariables
    implicit none

    Nline = 0
    Nlinemax = 0
    Nvertex = 0
    Nvertexmax = 0
    allocated_voronoi = .false.

    return
  END SUBROUTINE InitialiseVoronoiVariables



  !============================================================================
  ! AllocateVoronoiMemory
  ! Allocate memory for all arrays required to store data on the 
  ! Voronoi tessellation.
  !============================================================================
  SUBROUTINE AllocateVoronoiMemory(N)
    implicit none

    integer,intent(in) :: N             ! No. of vertices

    if (N > Nvertexmax) then
       !if (allocated_voronoi) deallocate(vertices)
       Nvertexmax = N
       Nlinemax = 20*N
       !allocate(vertices(1:Nvertexmax))
       allocated_voronoi = .true.
    end if

    return
  END SUBROUTINE AllocateVoronoiMemory



  !============================================================================
  ! GenerateVoronoiTessellation
  ! Generate the Delaunay triangulation from a set of N points.
  !============================================================================
  SUBROUTINE GenerateVoronoiTessellation(algorithm,N,vertices,&
       &xmin,xmax,ymin,ymax,zmin,zmax,boxsizeV)
    implicit none

    character(len=*),intent(in) :: algorithm  ! User-chosen algorithm
    integer,intent(in) :: N                   ! No. of points for triangulation
    real(kind=4),intent(in) :: xmin           ! ..
    real(kind=4),intent(in) :: xmax           ! ..
    real(kind=4),intent(in) :: ymin           ! ..
    real(kind=4),intent(in) :: ymax           ! ..
    real(kind=4),intent(in) :: zmin           ! ..
    real(kind=4),intent(in) :: zmax           ! ..
    real(kind=DP),intent(inout) :: boxsizeV   ! ..

    !integer :: Nvertexmax          ! Max. no. of Voronoi vertices
    !real(kind=4),intent(in) :: r(1:ndim,1:N)  ! Positions of points
    !integer  :: Nvertex            ! No. of vertices
    type(Vertex),intent(inout) :: vertices(1:N)  ! ..


    integer :: i,j                            ! Aux. counter
    integer :: p1,p2                          ! Aux. vertex ids
    integer :: Nconnecttot                    ! ..
    logical :: doneflag                       ! Connection already been done?
    real(kind=DP) :: dr(1:ndim)               ! Relative position vector
    real(kind=DP) :: drmag                    ! Distance
    real(kind=DP) :: drsqd                    ! Distance squared
    real(kind=DP) :: invdrmag                 ! 1 / distance
    real(kind=DP) :: volaux

    integer(c_int) :: c_N                     ! No. of points (c-type)
    integer(c_int) :: c_ndim                  ! No. of dimensiones (c-type)
    integer(c_int) :: c_Nline                 ! No. of lines (c-type)
    integer(c_int) :: c_Nlinemax              ! Max. no. of lines (c-type)
    integer(c_int) :: Ntempline               ! Temp. var for no. of lines
    real(c_double) :: c_rmin(1:ndim)          ! ..
    real(c_double) :: c_rmax(1:ndim)          ! ..
    real(c_double) :: c_boxsize               ! ..
    integer(c_int),allocatable :: i1(:)       ! Temp. particle pair array
    integer(c_int),allocatable :: i2(:)       ! Temp. particle pair array
    real(c_double),allocatable :: c_vol(:)    ! Voronoi volumes (c-type)
    real(c_double),allocatable :: c_r(:,:)    ! Positions of points (c-type)

    !interface
    !   subroutine qhull_generate_voronoi_tessellation(ndimaux,Naux,&
    !        &Nlinemaxaux,Nlineaux,i1,i2,raux,volaux) bind(c)
    !     use iso_c_binding
    !     integer(c_int),value :: ndimaux
    !     integer(c_int),value :: Naux
    !     integer(c_int),value :: Nlinemaxaux
    !     integer(c_int) :: Nlineaux
    !     integer(c_int) :: i1(1:Nlinemaxaux)
    !     integer(c_int) :: i2(1:Nlinemaxaux)
    !     real(c_double) :: raux(1:ndimaux,1:Naux)
    !     real(c_double) :: volaux(1:Naux)
    !   end subroutine qhull_generate_voronoi_tessellation
    !end interface

    interface
       subroutine voropp_generate_voronoi_tessellation(ndimaux,Naux,&
            &Nlinemaxaux,Nlineaux,i1,i2,raux,volaux,rminaux,rmaxaux) bind(c)
         use iso_c_binding
         integer(c_int),value :: ndimaux
         integer(c_int),value :: Naux
         integer(c_int),value :: Nlinemaxaux
         integer(c_int) :: Nlineaux
         integer(c_int) :: i1(1:Nlinemaxaux)
         integer(c_int) :: i2(1:Nlinemaxaux)
         real(c_double) :: raux(1:ndimaux,1:Naux)
         real(c_double) :: volaux(1:Naux)
         real(c_double) :: rminaux(1:ndimaux)
         real(c_double) :: rmaxaux(1:ndimaux)
       end subroutine voropp_generate_voronoi_tessellation

       subroutine voropp_delete_container() bind(c)
         use iso_c_binding
       end subroutine voropp_delete_container
    end interface


    if (taskid == 0)     write(6,*) "[GenerateVoronoiTessellation]"
    
    Nline = 0
    Nconnecttot = 0
    Nvertex = N

    ! Allocate memory to store all particle/vertex/line information
    call AllocateVoronoiMemory(N)
    allocate(c_r(1:ndim,1:N))
    allocate(c_vol(1:N))
    allocate(i1(1:Nlinemax))
    allocate(i2(1:Nlinemax))

    ! Convert all quantities to 'c-types' for calling c routines
    c_ndim = ndim
    c_N = N
    c_Nline = Nline
    c_Nlinemax = Nlinemax
    c_rmin(1) = xmin
    c_rmin(2) = ymin
    c_rmin(3) = zmin
    c_rmax(1) = xmax
    c_rmax(2) = ymax
    c_rmax(3) = zmax

    ! Re-size box for precision purposes
    c_boxsize = maxval(c_rmax(1:ndim) - c_rmin(1:ndim))
    c_rmin(1:ndim) = c_rmin(1:ndim) / c_boxsize
    c_rmax(1:ndim) = c_rmax(1:ndim) / c_boxsize
    boxsizeV = c_boxsize

    do i=1,N
!       c_r(1:ndim,i) = r(1:ndim,i)
       c_r(1:ndim,i) = vertices(i)%r(1:ndim)/c_boxsize
    end do


    ! Select Delaunay triangulation algorithm
    !--------------------------------------------------------------------------
    !if (algorithm == "qhull") then
    !   write(6,*) "Using qhull algorithm to generate Delaunay triangulation"
    !   
    !   call qhull_generate_voronoi_tessellation(c_ndim,c_N,&
    !        &c_Nlinemax,c_Nline,i1,i2,c_r,c_vol)
    !
    !--------------------------------------------------------------------------
    if (algorithm == "voro++") then
       if (taskid == 0) write(6,*) "Using Voro++ algorithm to generate Delaunay triangulation"
       
       call voropp_generate_voronoi_tessellation(c_ndim,c_N,&
            &c_Nlinemax,c_Nline,i1,i2,c_r,c_vol,c_rmin,c_rmax)

    end if
    !--------------------------------------------------------------------------
    

    Ntempline = c_Nline
    Nline = 0
    
    ! Write line information to screen
    if (taskid == 0) write(6,*) "Found ",Ntempline," lines"
    write(6,*) "Average number of temp lines per vertex : ",Ntempline/Nvertex


    ! Initialise data structures before processing data from tessellation
    do i=1,Nvertexmax
       vertices(i)%Nconnect = 0
    end do
    do i=1,Nvertex
       !volaux = real(c_vol(i),DP)*(c_boxsize**ndim)
       !vertices(i)%volume = volaux  !c_vol(i)*(c_boxsize**ndim)
       vertices(i)%volume = (c_vol(i)/1.d45)*(c_boxsize**ndim)
       !if (taskid == 0)  write(6,*) "VOLUME : ",i,c_vol(i),vertices(i)%volume
       if (vertices(i)%volume <= 0.) vertices(i)%volume = 1.d-15
    end do
    

    ! Now process lines to generate required data structures for ray tracing
    !--------------------------------------------------------------------------
    do i=1,Ntempline
       doneflag = .false.
       p1 = i1(i) + 1
       p2 = i2(i) + 1

       ! Check first that line pair has not already been included
       do j=1,vertices(p1)%Nconnect
          if (vertices(p1)%iconnect(j) == p2) doneflag = .true.
       end do
       do j=1,vertices(p2)%Nconnect
          if (vertices(p2)%iconnect(j) == p1) doneflag = .true.
       end do
       if (doneflag) cycle

       ! Check here that we've not exceeded max. no. of connections
       if (vertices(p1)%Nconnect == Nconnectmax .or. &
            & vertices(p2)%Nconnect == Nconnectmax) then

          ! Calculate average no. of connections
          write(6,*) "Average no. of connections : ",Nconnecttot/Nvertex
          write(6,*) "Checking connections for : ",p1,"   : ",&
               &vertices(p1)%Nconnect,vertices(p1)%r(1:ndim)
          do j=1,vertices(p1)%Nconnect
             write(6,*) "connection ",j,vertices(p1)%iconnect(j),&
                  &vertices(vertices(p1)%iconnect(j))%r(1:ndim),&
                  &vertices(p1)%connect(j)%length
          end do

          write(6,*) "Checking connections for : ",p2,"   : ",&
               &vertices(p2)%Nconnect,vertices(p2)%r(1:ndim)
          do j=1,vertices(p2)%Nconnect
             write(6,*) "connection ",j,vertices(p2)%iconnect(j),&
                  &vertices(vertices(p2)%iconnect(j))%r(1:ndim),&
                  &vertices(p2)%connect(j)%length
          end do

          if (taskid == 0) write(6,*) "Exceeded maximum number of connections for : ",&
               &p1,p2,vertices(p1)%Nconnect,vertices(p1)%Nconnect,Nconnectmax
          stop
       end if

       ! Compute distance and unit vector
       dr(1:ndim) = vertices(p2)%r(1:ndim) - vertices(p1)%r(1:ndim)
       drsqd = dot_product(dr(1:ndim),dr(1:ndim))
       drmag = sqrt(drsqd)
       invdrmag = 1.0/drmag
       dr(1:ndim) = dr(1:ndim)*invdrmag

       ! For new connections, calculate and store all important variables
       vertices(p1)%Nconnect = vertices(p1)%Nconnect + 1
       vertices(p1)%iconnect(vertices(p1)%Nconnect) = p2
       vertices(p1)%connect(vertices(p1)%Nconnect)%i1 = p1
       vertices(p1)%connect(vertices(p1)%Nconnect)%i2 = p2
       vertices(p1)%connect(vertices(p1)%Nconnect)%length = drmag
       vertices(p1)%connect(vertices(p1)%Nconnect)%invlength = invdrmag
       vertices(p1)%connect(vertices(p1)%Nconnect)%runit(1:ndim) = dr(1:ndim)

       vertices(p2)%Nconnect = vertices(p2)%Nconnect + 1
       vertices(p2)%iconnect(vertices(p2)%Nconnect) = p1
       vertices(p2)%connect(vertices(p2)%Nconnect)%i1 = p2
       vertices(p2)%connect(vertices(p2)%Nconnect)%i2 = p1
       vertices(p2)%connect(vertices(p2)%Nconnect)%length = drmag
       vertices(p2)%connect(vertices(p2)%Nconnect)%invlength = invdrmag
       vertices(p2)%connect(vertices(p2)%Nconnect)%runit(1:ndim) = -dr(1:ndim)

       Nconnecttot = Nconnecttot + 2

    end do
    !--------------------------------------------------------------------------


    write(6,*) "Average no. of connections : ",Nconnecttot/Nvertex


    deallocate(i2)
    deallocate(i1)
    deallocate(c_vol)
    deallocate(c_r)

    !call voropp_delete_container
    

    return
  END SUBROUTINE GenerateVoronoiTessellation



  !============================================================================
  ! FindVoronoiCell
  ! Find i.d. of Voronoi cell containing point (x,y,z).  If point is outside 
  ! tessellation, return -1.
  !============================================================================
  SUBROUTINE FindVoronoiCell(x,y,z,iVoronoi,boxsizeV)
    implicit none

    real, intent(in) :: x,y,z                 ! ..
    integer,intent(inout) :: iVoronoi         ! ..
    real(kind=DP),intent(in) :: boxsizeV      ! ..

    real(c_double) :: c_x,c_y,c_z,c_boxsize   ! ..
    integer(c_int) :: c_iVoronoi              ! ..

    interface
       function voropp_find_voronoi_cell(xaux,yaux,zaux) bind(c)
         use iso_c_binding
         integer(c_int) :: voropp_find_voronoi_cell
         !integer(c_int),value :: iaux
         real(c_double) :: xaux,yaux,zaux
       end function voropp_find_voronoi_cell
    end interface

    c_boxsize = boxsizeV
    c_x = real(x,DP)/c_boxsize
    c_y = real(y,DP)/c_boxsize
    c_z = real(z,DP)/c_boxsize
    
    if (taskid == 0) then
       write(6,*) "Original pos : ",x,y,z
       write(6,*) "boxsize      : ",c_boxsize
       write(6,*) "New pos      : ",c_x,c_y,c_z
    end if

    !if (algorithm == "voro++") then
    c_iVoronoi = voropp_find_voronoi_cell(c_x,c_y,c_z)
    iVoronoi = c_iVoronoi + 1
    !end if

    if (taskid == 0) then
       write(6,*) "New pos2     : ",c_x,c_y,c_z
       write(6,*) "iVoronoi     : ",iVoronoi,c_iVoronoi
    end if

  END SUBROUTINE FindVoronoiCell



  !============================================================================
  ! FindVoronoiCell
  ! Find i.d. of Voronoi cell containing point (x,y,z).  If point is outside 
  ! tessellation, return -1.
  !============================================================================
  SUBROUTINE DeleteContainer
    implicit none

    interface
       subroutine voropp_delete_container() bind(c)
         use iso_c_binding
       end subroutine voropp_delete_container
    end interface

    write(6,*) "[DeleteContainer]"

    call voropp_delete_container

    return
  END SUBROUTINE DeleteContainer



END MODULE VoronoiFortranInterface
