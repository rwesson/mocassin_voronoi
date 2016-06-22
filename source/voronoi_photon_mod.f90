!------------------------------------------------------------------------------
! voronoi_photon_mod
! Version of 'photon_mod' module for use in walking Voronoi tessellation to 
! propagate photon packets around computational domain.
! David Hubber & Barbara Ercolano - 30/10/2013
!------------------------------------------------------------------------------
module voronoi_photon_mod

  use common_mod
  use constants_mod
  use continuum_mod
  use voronoi_grid_mod
  use interpolation_mod
  use pathIntegration_mod
  use vector_mod
  
  ! common variables
  integer,parameter      :: safeLim = 10000          ! safety limit for no. of loops ??
  type(vector),parameter :: origin=vector(0.,0.,0.)  ! origin of cartesian grid axis

  integer                :: iPhot                    ! counter
  integer                :: totalEscapedV            ! total no. of escaped photon packets??
  real                   :: Qphot = 0.0              ! ??

contains
  
  
  !----------------------------------------------------------------------------
  subroutine energyPacketDriverV(iStar, n, grid, plot, gpLoc, cellLoc)
    implicit none
    
    integer,intent(in) :: n                        ! no. of energy packets 
    integer,intent(in) :: iStar                    ! central star index
    type(grid_type),intent(inout) :: grid(:)       ! main voronoi grid
    type(plot_type),intent(inout),optional :: plot ! only used in the mocassinPlot version
    integer,intent(inout),optional :: gpLoc        ! local grid (only used for extra diffuse sources) (Needed??) 
    integer,intent(inout),optional :: cellLoc      ! local cell (only used for extra diffuse sources)


    ! local variables
    integer          :: dt(8)             ! date and time values
    integer          :: freqP             ! pointer to frequency
    integer          :: gPIn              ! ..
    integer          :: i,iG,ii,j,k       ! counters        
    integer          :: ian               ! angle counter
    integer          :: iCell             ! cell counter
    integer          :: ierr              ! allocation error status
    integer          :: ifreq             ! freq counter
    integer          :: igp               ! always 1 for voronoi grid
    integer          :: igrid             ! ..
    integer          :: inV               ! voronoi cell id
    integer          :: iview             ! viewing angle counter
    integer          :: iV                ! voronoi cell id
    integer          :: msec              ! millisecs of the sec
    integer          :: plotNum           ! counter
    integer          :: reRun             ! flag to re-run photon packet
    integer          :: seedSize          ! pseudo random number generator seed 
    integer          :: trapped           ! no. of trapped photons
    integer, pointer :: seed(:)           ! seed array
    real             :: JDifTot           ! tot JDif
    real             :: JsteTot           ! tot Jste
    real             :: radius            ! radius
    character(len=7) :: chTypeD           ! character type for driver
    character(len=7) :: chTypeIn          ! character type  
    type(vector)     :: absPosition       ! position of packet absorption
    type(vector)     :: positionIn        ! position of packet absorption
    type(vector)     :: posDiff           ! initial position vector for diff ext
    type(vector)     :: posVector         ! initial position vector for dust emi


    ! Diffuse source
    if (iStar == 0) then
       deltaE(0) = grid(gpLoc)%LdiffuseLoc(grid(gpLoc)%activeV&
            &(cellLoc))/NphotonsDiffuseLoc
    end if
    
    ! Obtain timing information (used to generate 'true' random number)
    call date_and_time(values=dt)
    msec = dt(8)
    
    ! Allocate memory for containing random number data
    call random_seed(seedSize) 
    allocate(seed(1:seedSize), stat= ierr)
    if (ierr /= 0) then
       print*, "energyPacketDriverV: can't allocate array memory: seed"
       stop
    end if

    ! Initial random number seed using seed + clock time info
    seed = 0
    call random_seed(get = seed)
    seed = seed + msec + taskid

    ! to force the same random sequence de-comment the following 2 lines
    ! seed(1) = -2147482802 
    ! seed(2) = -2147482966
    call random_seed(put = seed)
    if (associated(seed)) deallocate(seed)
    
    ! Initialise arrays and counters
    Qphot = 0.0
    trapped = 0


    ! Main photon packet loop
    !--------------------------------------------------------------------------
    do iPhot = 1, n

       ! If photon is emitted directly from a star
       !-----------------------------------------------------------------------
       if (iStar >= 1) then
          
          chTypeD = "stellar"
          igp = 1
          inV = starIndecesV(iStar)
          chTypeIn = chTypeD
          positionIn = starPosition(iStar)  !!DAH
          

       ! If photon is a diffuse (scattered or emitted) photon
       !-----------------------------------------------------------------------
       else if (iStar == 0) then

          chTypeD = "diffExt"              
          igp = 1
          inV = cellLoc
          posDiff%x = grid(gpLoc)%voronoi(cellLoc)%r(1)
          posDiff%y = grid(gpLoc)%voronoi(cellLoc)%r(2)
          posDiff%z = grid(gpLoc)%voronoi(cellLoc)%r(3)
          chTypeIn = chTypeD
          positionIn = posDiff
          gPIn = gPLoc


       ! Stop if something has gone wrong
       !-----------------------------------------------------------------------
       else
          
          print*, '! energyPacketDriverV: insanity in iStar value'
          stop
          
       end if
       !-----------------------------------------------------------------------

       reRun = 0
       do i = 1, recursionLimit
          call energyPacketRunV(chTypeIn, positionIn, inV, reRun)
          if (rerun == 0) exit
       end do
       if (i >= recursionLimit) trapped = trapped + 1
       
    end do
    !--------------------------------------------------------------------------

    if (iStar > 0) then
       print*, 'Star: ', iStar
       print*, 'Qphot = ', Qphot
    end if


    ! ..
    !--------------------------------------------------------------------------
    if (lgDust .and. convPercent>=resLinesTransfer .and.&
         & .not. lgResLinesFirst&
         & .and. (.not. nIterateMC==1) .and. .not.lgVoronoi) then
       
       print*, "! energyPacketDriverV: starting resonance line packets transfer"
       
       
       iCell = 0
       igrid = 1


       ! Loop over all Voronoi cells 
       !-----------------------------------------------------------------------
       do iV = 1,grid(igrid)%nCellsV
          iCell = iCell+1
                   
          !--------------------------------------------------------------------
          if (mod(iCell-(taskid+1),numtasks) == 0 .and. &
               & grid(igrid)%activeV(iV) > 0) then
                         
             !-----------------------------------------------------------------
             do iPhot = 1, grid(igrid)%resLinePackets(grid(igrid)%activeV(iV))
                
                chTypeD = "diffuse"
                inV = iV                                
                posVector%x = grid(igrid)%voronoi(iV)%r(1)
                posVector%y = grid(igrid)%voronoi(iV)%r(2)
                posVector%z = grid(igrid)%voronoi(iV)%r(3)                            
                chTypeIn = chTypeD
                positionIn = posVector
                gPIn = igrid
                reRun = 0
                
                do i = 1, recursionLimit      
                   call energyPacketRunV(chTypeIn, positionIn, inV, reRun)
                   if (rerun == 0) exit                      
                end do
                if (i >= recursionLimit) trapped = trapped + 1
                
             end do
             !-----------------------------------------------------------------
             
          end if
          !--------------------------------------------------------------------

       end do
       !-----------------------------------------------------------------------
       
       print*, "! energyPacketDriverV: ending resonance line packets transfer"
       
    end if
    !--------------------------------------------------------------------------    

    
    if (iStar>0) print*, 'Qphot = ', Qphot
    print*, 'Packets trapped = ', trapped, taskid
    
    ! evaluate Jste and Jdif
    ! NOTE : deltaE is in units of [E36 erg/s] however we also need a factor of
    ! 1.e-45 from the calculations of the volume of the cell hence these
    ! two factors cancel each other out giving units of [E-9erg/s] so we need to
    ! multiply by 1.E-9
    ! NOTE : Jste and Jdif calculated this way are in units of
    ! [erg sec^-1 cm^-2] -> no Hz^-1 as they are summed over separate bins (see
    ! Lucy A&A (1999)                                                                                                                                   
    
    if (iStar > 0.0) then
       print*, 'Lstar', Lstar(iStar)           
    end if
    
    

  contains



    !--------------------------------------------------------------------------
    subroutine energyPacketRunV(chType, position, vP, rR)
      implicit none
      
      character(len=7),intent(inout) :: chType    ! stellar or diffuse?
      type(vector),intent(inout) :: position      ! the position of the photon
      integer,intent(inout) :: vP                 ! voronoi cell index
      integer,intent(inout) :: rR                 ! rerun?

      integer                 :: difSourceV       ! ..
      integer                 :: err              ! allocation error status
      integer                 :: gP               ! ..
      integer                 :: i, j             ! counters  
      integer                 :: idirP, idirT     ! direction cosines
      integer                 :: igpr             ! grid pointer 1= mother 2=sub
      type(photon_packet)     :: enPacket         ! the energu packet


      ! Hard-wire no re-run??
      rR = 0
      gP = 1


      ! create a new photon packet
      !-------------------------------------------------------------------------   
      select case (chType)
      
   
      ! if the energy packet is stellar
      !-------------------------------------------------------------------------   
      case ("stellar")

         ! check for errors in the sources position
         if( position /= starPosition(iStar) ) then
            print*, "! energyPacketRun: stellar energy packet must&
                 & start at the stellar position"
            stop
         end if
         
         ! create the packet        
         enPacket = newPhotonPacketV(chType,position, vP, noCellLocV)
         
      ! if the photon is from an extra source of diffuse radiation
      !-------------------------------------------------------------------------   
      case ("diffExt")
         
         difSourceV = vP         
         enPacket = newPhotonPacketV(chType, position, vP, difSourceV)
         
         
      ! if the photon is diffuse
      !-------------------------------------------------------------------------   
      case ("diffuse")
         
         ! create the packet
         enPacket = newPhotonPacketV(chType=chType, position=position, &
              &iVP=vP, difSourceV=noCellLocV)
         
      !-------------------------------------------------------------------------   
      case ("dustEmi")

         ! create the packet
         enPacket = newPhotonPacketV(chType=chType, position=position, &
              &iVP=vP, difSourceV = noCellLocV)
         
      end select
      !-------------------------------------------------------------------------   

      !-------------------------------------------------------------------------   
      if (.not.lgDust .and. enPacket%nu < ionEdge(1) .and. .not.enPacket%lgLine) then

         call photonPacketEscapes(enPacket)
         return

      end if

      
      ! if the new packet is capable of re-ionizing or we have dust
      !-------------------------------------------------------------------------   
      if (.not.enPacket%lgLine) then

         ! compute the next segment of trajectory
         call pathSegmentV(enPacket)
         
         return
         
      else ! if the packet is a line packet 
         ! add to respective line packet bin
         ! this was just used for debugging

      end if
      
    end subroutine energyPacketRunV



    ! this function initializes a photon packet
    !--------------------------------------------------------------------------
    function initPhotonPacketV(nuP, position, direction, lgLine, lgStellar, iVP, lgHg)
      implicit none

      integer, intent(in) :: nuP               ! the frequency of the photon packet
      integer, intent(in) :: iVP               ! ..
      logical, intent(in) :: lgLine, lgStellar ! line, stellar packet
      logical, intent(in) :: lgHG              ! Heyney-Greenstein
      type(vector), intent(in) :: position     ! position at which photon packet is created
      type(vector), intent(in) :: direction    ! direction of photon packet??

      integer             :: iboundary         ! ..
      integer             :: igpi              ! grid pointer 1=mother, 2=sub
      integer             :: i, irepeat        ! counter
      integer             :: gP                ! grid number (always 1 on Voronoi)
      real                :: random            ! random number     
      double precision    :: Smin              ! ..
      type(vector)        :: pn, pp            ! ..
      type(photon_packet) :: initPhotonPacketV ! the photon packet
      
      
      initPhotonPacketV%direction = direction
      
      if (.not. (direction%x >= 0. .or. direction%x < 0)) then
         print*, '! initPhotonPacketV: [0] insane direction%x [direction%x]'
         print*, direction%x
         stop
      end if
      
      
      gP = 1
      igpi = 1
      
      initPhotonPacketV%position  = position
      initPhotonPacketV%iG        = gP
      initPhotonPacketV%nuP       = nuP       
      initPhotonPacketV%lgStellar = lgStellar

      ! check if photon packet is line or continuum photon
      if ( lgLine ) then
         initPhotonPacketV%nu       = 0.
         initPhotonPacketV%lgLine   = .true.
      else
         initPhotonPacketV%nu       = nuArray(nuP)
         initPhotonPacketV%lgLine   = .false.          
      end if
      
      initPhotonPacketV%iVoronoi = iVP
      

      !DAVID : Probably needs substantial re-writing for Voronoi
      !------------------------------------------------------------------------
      if (.not.lgHG .or. lgIsotropic .or. initPhotonPacketV%lgStellar) then

         ! cater for plane parallel ionization case
         !---------------------------------------------------------------------
         if (initPhotonPacketV%lgStellar .and. lgPlaneIonization) then
      !      
      !      ! get position
      !      
      !      ! x-direction
      !      call notrandom2(random)
      !      random = 1. - random
      !      initPhotonPacketV%position%x = &
      !           & -(grid(gP)%xAxis(2)-grid(gP)%xAxis(1))/2. + random*( &
      !           & (grid(gP)%xAxis(2)-grid(gP)%xAxis(1))/2.+&
      !           & (grid(gP)%xAxis(grid(gP)%nx)-grid(gP)%xAxis(grid(gP)%nx-1))/2.+&
      !           & grid(gP)%xAxis(grid(gP)%nx))
      !      if (initPhotonPacketV%position%x<grid(gP)%xAxis(1)) &
      !           & initPhotonPacketV%position%x=grid(gP)%xAxis(1)
      !      if (initPhotonPacketV%position%x>grid(gP)%xAxis(grid(gP)%nx))&
      !           & initPhotonPacketV%position%x=grid(gP)%xAxis(grid(gP)%nx)
      !      
      !      call locate(grid(gP)%xAxis, initPhotonPacketV%position%x,&
      !           & initPhotonPacketV%xP(igpi))
      !      if (initPhotonPacketV%xP(igpi) < grid(gP)%nx) then
      !         if (initPhotonPacketV%position%x >= &
      !              &(grid(gP)%xAxis(initPhotonPacketV%xP(igpi))+&
      !              & grid(gP)%xAxis(initPhotonPacketV%xP(igpi)+1))/2.) &
      !              & initPhotonPacketV%xP(igpi) = initPhotonPacketV%xP(igpi)+1
      !      end if
      !      
      !      ! y-direction
      !      initPhotonPacketV%position%y = 0.
      !      initPhotonPacketV%yP(igpi) = 1
      !      
      !      ! z-direction
      !      call notrandom2(random)
      !      random = 1. - random
      !      initPhotonPacketV%position%z = -&
      !           & (grid(gP)%zAxis(2)-grid(gP)%zAxis(1))/2. &
      !           & + random*( &
      !           & (grid(gP)%zAxis(2)-grid(gP)%zAxis(1))/2.+&
      !           & (grid(gP)%zAxis(grid(gP)%nz)-&
      !           & grid(gP)%zAxis(grid(gP)%nz-1))/2.+&
      !           & grid(gP)%zAxis(grid(gP)%nz))
      !      if (initPhotonPacketV%position%z<grid(gP)%zAxis(1)) & 
      !           & initPhotonPacketV%position%z=grid(gP)%zAxis(1)
      !      if (initPhotonPacketV%position%z>&
      !           & grid(gP)%zAxis(grid(gP)%nz)) & 
      !           & initPhotonPacketV%position%z=&
      !           & grid(gP)%zAxis(grid(gP)%nz)
      !      
      !      call locate(grid(gP)%zAxis, &
      !           & initPhotonPacketV%position%z, initPhotonPacketV%zP(igpi))
      !      if (initPhotonPacketV%zP(igpi) < grid(gP)%nz) then 
      !         if (initPhotonPacketV%position%z >= &
      !              & (grid(gP)%xAxis(initPhotonPacketV%zP(igpi))+&
      !              & grid(gP)%zAxis(initPhotonPacketV%zP(igpi)+1))&
      !              & /2.) initPhotonPacketV%zP(igpi) =& 
      !              & initPhotonPacketV%zP(igpi)+1
      !      end if
      !      
      !      if (initPhotonPacketV%xP(igpi)<1) &
      !           & initPhotonPacketV%xP(igpi)=1             
      !      if (initPhotonPacketV%zP(igpi)<1) &
      !           & initPhotonPacketV%zP(igpi)=1
      !      
      !      ! direction is parallel to y-axis direction
      !      initPhotonPacketV%direction%x = 0.
      !      initPhotonPacketV%direction%y = 1.
      !      initPhotonPacketV%direction%z = 0.
      !      
      !      if (initPhotonPacketV%xP(igpi) >  &
      !           & grid(gP)%xAxis(grid(gP)%nx) .or. &
      !           & initPhotonPacketV%zP(igpi) >  &
      !           & grid(gP)%zAxis(grid(gP)%nz)) then
      !         print*, "! initPhotonPacketV: insanity in &
      !              & planeIonisation init"
      !         print*, igpi, initPhotonPacketV%xP(igpi),  &
      !              & grid(gP)%xAxis(grid(gP)%nx), &
      !              & initPhotonPacketV%zP(igpi), &
      !              & grid(gP)%zAxis(grid(gP)%nz),  random, &
      !              & initPhotonPacketV%position%z
      !         
      !         stop
      !      end if
      !      
      !      planeIonDistribution(initPhotonPacketV%xP(igpi),&
      !           & initPhotonPacketV%zP(igpi)) = &
      !           & planeIonDistribution(initPhotonPacketV%xP(igpi),&
      !           & initPhotonPacketV%zP(igpi)) + 1
      !      
         !---------------------------------------------------------------------
         else
            
            do irepeat = 1, 1000000
               ! get a random direction
               initPhotonPacketV%direction = randomUnitVector(iphot)
               
               if (.not. (initPhotonPacketV%direction%x >= 0. .or. initPhotonPacketV%direction%x < 0)) then
                  print*, '! initPhotonPacketV: [1] insane initPhotonPacketV%direction%x [initPhotonPacketV%direction%x]'
                  print*, initPhotonPacketV%direction%x
                  stop
               end if
               
               if (initPhotonPacketV%direction%x/=0. .and. & 
                    & initPhotonPacketV%direction%y/=0. .and. & 
                    & initPhotonPacketV%direction%z/=0.) exit
            end do
            
         end if
         !---------------------------------------------------------------------

         
         if ((lgSymmetricXYZ) .and. initPhotonPacketV%lgStellar &
              & .and. .not.lgMultistars) then
!            print*, '! initPhotonPacketV: SymmetricXYZ is not yet implemented in Voronoi'
!            stop

            if (initPhotonPacketV%direction%x<0.) &
                 & initPhotonPacketV%direction%x = -initPhotonPacketV%direction%x
            if (initPhotonPacketV%direction%y<0.) &
                 & initPhotonPacketV%direction%y = -initPhotonPacketV%direction%y
            if (initPhotonPacketV%direction%z<0.) &
                 & initPhotonPacketV%direction%z = -initPhotonPacketV%direction%z
         end if
         
      end if
      !------------------------------------------------------------------------

      
      initPhotonPacketV%origin(1) = gP
      initPhotonPacketV%origin(2) = grid(gP)%activeV(initPhotonPacketV%iVoronoi)

      ! Calculate path to next boundary (mirror or open)
      call calculateBoundaryPath(initPhotonPacketV%direction,&
           &initPhotonPacketV%position,Smin,iboundary)

      initPhotonPacketV%STot = 0.0
      initPhotonPacketV%SEscape = Smin

      if (lgSymmetricXYZ) then
         initPhotonPacketV%iboundary = iboundary
      else
         initPhotonPacketV%iboundary = 0
      end if

    end function initPhotonPacketV

    

    ! this function initializes a photon packet
    !--------------------------------------------------------------------------
    subroutine calculateBoundaryPath(direction,position,Smin,iboundary)
      implicit none

      type(vector),intent(in)      :: direction   ! ..
      type(vector),intent(in)      :: position    ! ..
      double precision,intent(out) :: Smin        ! ..
      integer,intent(out)          :: iboundary   ! ..

      double precision             :: Saux        ! Minimum plane intersection distance
      type(vector)                 :: n0,pn,pp    ! ..


      ! Find out the distance to the first boundary along the photon path
      ! and record for future checking if the photon intercepts the boundary
      !------------------------------------------------------------------------
      Smin = 9.9e30 ! BE was 9.9e40
      pn%x = boxXn;  pn%y = boxYn;  pn%z = boxZn
      pp%x = boxXp;  pp%y = boxYp;  pp%z = boxZp

      ! Check if photon packet strikes left or right x-boundary
      if (direction%x > 0.0) then
         n0%x = 1.0;  n0%y = 0.0;   n0%z = 0.0
         Saux = ((pp - position).dot.n0)/&
              &(direction.dot.n0)
         if (Saux < Smin) then
            Smin = Saux
            iboundary = 0
         end if
      else
         n0%x = -1.0;  n0%y = 0.0;   n0%z = 0.0
         Saux = ((pn - position).dot.n0)/&
              &(direction.dot.n0)
         if (Saux < Smin) then
            Smin = Saux
            iboundary = 1
         end if
      end if

      if (direction%y > 0.0) then
         n0%x = 0.0;  n0%y = 1.0;   n0%z = 0.0
         Saux = ((pp - position).dot.n0)/&
              &(direction.dot.n0)
         if (Saux < Smin) then
            Smin = Saux
            iboundary = 0
         end if
      else
         n0%x = 0.0;  n0%y = -1.0;   n0%z = 0.0
         Saux = ((pn - position).dot.n0)/&
              &(direction.dot.n0)
         if (Saux < Smin) then
            Smin = Saux
            iboundary = 2
         end if
      end if

      if (direction%z > 0.0) then
         n0%x = 0.0;  n0%y = 0.0;   n0%z = 1.0
         Saux = ((pp - position).dot.n0)/&
              &(direction.dot.n0)
         if (Saux < Smin) then
            Smin = Saux
            iboundary = 0
         end if
      else
         n0%x = 0.0;  n0%y = 0.0;   n0%z = -1.0
         Saux = ((pn - position).dot.n0)/&
              &(direction.dot.n0)
         if (Saux < Smin) then
            Smin = Saux
            iboundary = 3
         end if
      end if

      if (Smin < 0.0) then
         print*, '! initPhotonPacketV: Negative path length to exit surface'
         stop
      end if

      return
    end subroutine calculateBoundaryPath

    

    ! this subroutine determines the frequency of a newly created photon packet
    ! according to the given probability density
    !---------------------------------------------------------------------------
    subroutine getNuV(probDen, nuP)
      
      real, dimension(:), intent(in) :: probDen    ! probability density function   
      integer, intent(out)           :: nuP        ! frequency index of the new
      
      real                           :: random     ! random number
      
      ! get a random number
      call random_number(random)
      random = 1.0 - random

      
      ! see what frequency random corresponds to 
      call locate(probDen, random, nuP)
      if (nuP <= 0) nuP = 1
      
      if (nuP<nbins) then
         if (random>=(probDen(nuP+1)+probDen(nuP))/2.) nuP=nuP+1
      end if
      
    end subroutine getNuV


    
    ! this subroutine determines the frequency of a newly created photon packet
    ! according to the given probability density
    ! does not use bisection to locate nu on array
    !---------------------------------------------------------------------------
    subroutine getNu2V(probDen, nuP)
      
      real, dimension(:), intent(in) :: probDen    ! probability density function
      integer, intent(out)           :: nuP        ! frequency index of the new
      
      integer                        :: isearch,i  ! ..
      real                           :: random     ! random number
      
            
      ! get a random number
      call random_number(random)
      
      ! ensure random number does not have 'special values' (e.g. 0 or 1)
      do i = 1, 10000
         if (random==0 .or. random==1. .or. random == 0.9999999) then
            call random_number(random)
         else
            exit
         end if
      end do
      if (i>=10000) then
         print*, '! getNu2: problem with random number generator', random, i
         stop
      end if
      
      ! see what frequency random corresponds to 
      nuP=1
      do isearch = 1, nbins
         if (random >= probDen(isearch)) then
            nuP = isearch
         else 
            exit                  
         end if
      end do
      
      if (nuP < nbins - 1) then
         nuP = nuP + 1
      end if
      
      if (nuP>=nbins) then
         print*, 'random: ', random
         print*, 'probDen: ', probDen
      end if
      
    end subroutine getNu2V
    
    

    ! this function creates a new photon packet
    !---------------------------------------------------------------------------
    function newPhotonPacketV(chType, position, iVP, difSourceV)

      character(len=7), intent(in) :: chType         ! stellar or diffuse?
      integer, intent(in)          :: iVP            ! voronoi cell index of photon
      integer, intent(in)          :: difSourceV     ! Grid/Voronoi indicies
      type(vector), intent(in), optional :: position ! the position of the photon

      logical             :: lgLine_loc=.false.   ! line photon?
      integer             :: gP                   ! ..
      integer             :: igpn                 ! grid pointe 1=motehr, 2=sub
      integer             :: nuP                  ! frequency index of the photon packet
      integer             :: orV                  ! ..
      real                :: random               ! random number
      type(vector)        :: positionLoc          ! position of the photon packet
      type(photon_packet) :: newPhotonPacketV     ! photon packet to be created
      

      gP = 1
      igpn = 1

      
      !-------------------------------------------------------------------------
      select case (chType)
      
   
      ! if the photon is stellar
      !-------------------------------------------------------------------------
      case ("stellar")
         
         ! check for errors in the sources position
         if( position /= starPosition(iStar) ) then
            print*, "! newPhotonPacketV: stellar photon packet must&
                 & start at the stellar position"
            stop
         end if
         
         
         !gP = starIndeces(iStar,4)
         igpn = 1
         
         ! determine the frequency of the newly created photon packet
         call getNu2V(inSpectrumProbDen(iStar,1:nbins), nuP)
         
         if (nuP >= nbins) then
            print*, "! newPhotonPacketV: insanity occured in stellar photon &
                 &nuP assignment (nuP,iVP,activeP)", nuP, iVP, &
                 & grid(starIndeces(iStar,4))%activeV(iVP)
            print*, "inSpectrumProbDen: ",iStar,inSpectrumProbDen(iStar,:), nuP
            stop
         end if
         
         if (nuP < 1) then
            print*, "! newPhotonPacketV: insanity occured in stellar photon &
                 &                               nuP assignment"
            stop
         end if
         
         ! initialize the new photon packet
         orV = starIndecesV(iStar)

         if (grid(igpn)%activeV(orV) < 0) then
            print*, "! newPhotonPacketV: new packet cannot be emitted from re-mapped cell -1-" 
            print*, "chType, nuP, starPosition(iStar), .false., .true., orV"
            print*, chType, nuP, starPosition(iStar), .false., .true., orV
            stop
         end if
         

         newPhotonPacketV = initPhotonPacketV(nuP, starPosition(iStar), &
              & nullUnitVector, .false., .true., orV, .false.)
         
         
         if (newPhotonPacketV%nu>1.) then
            Qphot = Qphot + deltaE(iStar)/(2.1799153e-11*newPhotonPacketV%nu)
         end if
      
   
      ! if the photon is from an extra diffuse source
      !-------------------------------------------------------------------------
      case ("diffExt")
         
         call getNu2V(inSpectrumProbDen(0,1:nbins), nuP)
         
         if (nuP>=nbins) then
            print*, "! newPhotonPacketV: insanity occured in extra diffuse photon &
                 & nuP assignment (nuP,gp,activeP)", nuP, gp
            print*, "difSpectrumProbDen: ", inSpectrumProbDen(0,:)
            stop
         end if
         
         if (nuP < 1) then
            print*, "! newPhotonPacketV: insanity occured in extra diffuse photon &
                 & nuP assignment (nuP,gp,activeP)", nuP, gp,grid(gP)%activeV(iVP)
            print*, "difSpectrumProbDen: ", inSpectrumProbDen(0,:)
            stop
         end if
         
         positionLoc%x = grid(gP)%voronoi(difSourceV)%r(1)
         positionLoc%y = grid(gP)%voronoi(difSourceV)%r(2)
         positionLoc%z = grid(gP)%voronoi(difSourceV)%r(3)
         
         ! initialize the new photon packet
         orV = difSourceV
         
         if (grid(gP)%activeV(orV) < 0) then
            print*, "! newPhotonPacketV: new packet cannot be emitted from re-mapped cell -3-"
            print*, "chType, nuP, starPosition(iStar), .false., .true., orV, gp"
            print*, chType, nuP, starPosition(iStar), .false., .true., orV
            stop
         end if
         
         newPhotonPacketV = initPhotonPacketV(nuP, positionLoc, &
              & nullUnitVector, .false., .false., orV, .false.)
         
         
      ! if the photon is diffuse
      !-------------------------------------------------------------------------
      case ("diffuse")
         
         ! check that gas is present in the grid
         if (.not.lgGas) then
            print*, "! newPhotonPacketV: diffuse packet cannot be created in a no gas grid"
            stop
         end if
         
         ! check that the position is not inside the inner region
         if (grid(gP)%activeV(iVP)<= 0) then                   
            print*, "! newPhotonPacketV: position of the new diffuse &
                 & photon packet is inside the inner region", iVP
            stop
         end if
         
         ! decide whether continuum or line photon
         call random_number(random)
         random = 1.0 - random
         

         !----------------------------------------------------------------------
         if (random <= grid(gP)%totalLines(grid(gP)%activeV(iVP))) then 

            nuP = 0
            
            ! initialize the new photon packet
            if (grid(gP)%activeV(iVP) < 0.) then
               print*, "! newPhotonPacketV: new packet cannot be emitted from re-mapped cell -4-"
               print*, "chType, nuP, starPosition(iStar), .false., .true., iVP"
               print*, chType, nuP, starPosition(iStar), .false., .true., iVP
               stop
            end if
            
            newPhotonPacketV = initPhotonPacketV(nuP, position, &
                 & nullUnitVector, .true., .false., iVP, .false.)

         ! continuum photon
         !----------------------------------------------------------------------
         else 
            
            call getNu2V(grid(gP)%recPDF(grid(gP)%activeV(iVP),1:nbins), nuP)
            
            if (nuP>=nbins) then
               print*, "! newPhotonPacketV: insanity occured in diffuse photon &
                    & nuP assignment (nuP,iVP,activeP)", nuP, iVP,&
                    &  grid(gP)%activeV(iVP)
               print*, "recPDF: ", grid(gP)%recPDF(grid(gP)%activeV(iVP),1:nbins)
               stop
            end if
            
            if (nuP < 1) then
               print*, "! newPhotonPacketV: insanity occured in diffuse photon &
                    & nuP assignment", nuP, iVP, & 
                    & grid(gP)%activeV(iVP)
               print*, "recPDF: ", grid(gP)%recPDF(grid(gP)%activeV(iVP),:)
               stop
            end if
            
            ! initialize the new photon packet
            if (grid(gP)%activeV(iVP) < 0) then
               print*, "! newPhotonPacketV: new packet cannot be emitted from re-mapped cell -5-"
               print*, "chType, nuP, starPosition(iStar), .false., .true., iVP"
               stop
            end if
            
            
            newPhotonPacketV = initPhotonPacketV(nuP, position, nullUnitVector,&
                 & .false., .false., iVP, .false.)
            
         end if
         !----------------------------------------------------------------------
         

      !-------------------------------------------------------------------------
      case ("dustEmi")
         
         ! check dust is present
         if (.not.lgDust) then
            print*, "! newPhotonPacketV: dust emitted packet cannot be created in a &
                 &no dust grid."
            stop
         end if
         
         if (lgGas) then
            print*, "! newPhotonPacketV: dustEmi-type packet should be created in a &
                 & grid containing gas."
            stop
         end if
         
         ! check that the position is not inside the inner region
         if (grid(gP)%activeV(iVP)<= 0) then
            print*, "! newPhotonPacketV: position of the new dust emitted &
                 &photon packet is inside the inner region", iVP
            stop
         end if
         
         call getNu2V(grid(gP)%dustPDF(grid(gP)%activeV(iVP),1:nbins), nuP)
         
         if (nuP >= nbins) then
            print*, "! newPhotonPacketV: insanity occured in dust emitted photon &
                 &nuP assignment (iphot, nuP, iVP,activeP)", iphot, &
                 & nuP, iVP, &
                 & grid(gP)%activeV(iVP)
            print*, "dustPDF: ", grid(gP)%dustPDF(grid(gP)%activeV(iVP),1:nbins)
            print*, "grain T: ", grid(gP)%Tdust(:,0,grid(gP)%activeV(iVP))
            stop
         end if
         
         if (nuP < 1) then
            print*, "! newPhotonPacketV: insanity occured in dust emitted photon &
                 &nuP assignment", nuP,iVP,grid(gP)%activeV(iVP)
            print*, "dustPDF: ", grid(gP)%dustPDF(grid(gP)%activeV(iVP),1:nbins)
            print*, "grain T: ", grid(gP)%Tdust(:, 0, grid(gP)%activeV(iVP))
            stop
         end if
         
         ! initialize the new photon packet
         if (grid(gP)%activeV(iVP) < 0.) then
            print*, "! newPhotonPacketV: new packet cannot be emitted from re-mapped cell -6-"
            print*, "chType, nuP, starPosition(iStar), .false., .true., iVP"
            print*, chType, nuP, starPosition(iStar), .false., .true., iVP
            stop
         end if
         
         newPhotonPacketV = initPhotonPacketV(nuP, position, nullUnitVector,&
              & .false., .false., iVP, .false.)
         
         
         
      ! if the photon packet type is wrong or missing
      !-------------------------------------------------------------------------
      case default
         
         print*, "! newPhotonPacketV: wrong photon packet type - 'stellar', 'diffuse' and &
              & dust emitted types allowed-"
         stop
      end select
      !-------------------------------------------------------------------------

      
    end function newPhotonPacketV
    


    !--------------------------------------------------------------------------
    subroutine pathSegmentV(enPacket)
      implicit none

      type(photon_packet),intent(inout) :: enPacket   ! energy packet

      type(vector)         :: rVec2,vHat2             ! position vector
      type(Vertex),pointer :: v                       ! pointer to current voronoi cell
      character(len=7)     :: packetType              ! photon packet type
      logical              :: lgScattered             ! is the packet scattering with dust?
      logical              :: lgReturn                ! flag to return from routine
      integer              :: gP                      ! grid index
      integer              :: i,j                     ! ..
      integer              :: iboundary               ! ..
      integer              :: igpp,ihg                ! ..
      integer              :: iierr                   ! ..
      integer              :: inext                   ! id of next Voronoi cell
      integer              :: iprev                   ! id of previous Voronoi cell
      integer              :: iVoronoi                ! id of Voronoi cell
      integer              :: jcell                   ! id of neighbour Voronoi cell
      integer              :: jj                      ! aux. neighbour cell counter
      integer              :: jjmin                   ! ..
      integer              :: nS                      ! ..
      integer              :: safeLimit               ! ..
      double precision     :: absTau                  ! optical depth
      double precision     :: dlLoc                   ! local displacement
      double precision     :: dr(1:3)                 ! relative position vector
      double precision     :: dS                      ! distance to next voronoi cell face
      double precision     :: dSaux                   ! ..
      double precision     :: dV                      ! volume of Voronoi cell
      double precision     :: invr12comp              ! 1 / r12comp
      double precision     :: probSca                 ! ..
      double precision     :: random                  ! random number
      double precision     :: rVec(1:3)               ! position vector
      double precision     :: r1sqd                   ! distance (squared) of point 1
      double precision     :: r2sqd                   ! distance (squared) of point 2
      double precision     :: r12comp                 ! component of relative position vector
      double precision     :: Smin                    ! ..
      double precision     :: tauCell                 ! optical depth through cell
      double precision     :: totTau                  ! tot optical depth travelled by photon
      double precision     :: vHat(1:3)               ! direction vector

      ! define rVec and vHat (position and direction of photon packet)
      rVec(1) = enPacket%position%x
      rVec(2) = enPacket%position%y
      rVec(3) = enPacket%position%z
      vHat(1) = enPacket%direction%x
      vHat(2) = enPacket%direction%y
      vHat(3) = enPacket%direction%z
      iVoronoi = enPacket%iVoronoi
      gp = 1  !enPacket%iG
      v => grid(gp)%voronoi(iVoronoi)  
      iprev = iVoronoi


      ! sanity check for valid values and compatible options
      if (enPacket%iG /= 1) then
         print *, "! pathSegment: insane grid index : ",enPacket%iG
         stop
      end if
      if (lg1D) then
         print *, "1D simulations not yet implemented for Voronoi cell walk"
         stop
      end if
!      if (lgSymmetricXYZ) then
!         print *, "3D-symmetric simulations not yet implemented for Voronoi cell walk"
!         stop
!      end if
      if (lgPlaneIonization) then
         print *, "Plane ionization not yet implemented for Voronoi cell walk"
         stop
      end if


      ! initialise optical depth
      absTau = 0.0
      
      ! calculate the total optical depth to be travelled by the photon packet

      call random_number(random)
!      random = random2

      totTau = -log(1.0 - random)

      ! speed up photons that may be trapped
      !if (lgPlaneIonization) then
      !   safeLimit = 5000
      !else
      !   safeLimit = 50000
      !end if
      safeLimit = 50000


      ! Loop through photon propagation and scatterings until either photon 
      ! escapes computational domain, or reaches maximum limit (safeLimit) 
      ! for optically thick media.
      !------------------------------------------------------------------------
      do i=1,safeLimit

            
         ! verify that ray is inside a valid Voronoi cell
         ! (need to find calls to relevant Voronoi library routines)
         if (iVoronoi < 1 .and. iVoronoi > grid(gp)%nCellsV) then
            print*, "! pathSegmentV: insanity [gp,iVoronoi]", gp,iVoronoi
            stop
         end if
         
         ! (DAVID : Missing code here for symmetric walls; implement in the future)
         ! (Update : Not needed; walls done in different way)
         
         ! initialise variables before finding voronoi path length
         dr(1:3) = v%r(1:3) - rVec(1:3)
         r1sqd = dot_product(dr,dr)
         dS = 9.9e30
         inext = -1
         jjmin = -1
         
         
         ! loop over all voronoi connections to find minimum face intersection distance
         !---------------------------------------------------------------------
         do jj=1,v%Nconnect
            jcell = v%iconnect(jj)
            if (jcell == iprev) cycle
            r12comp = dot_product(v%connect(jj)%runit,vHat)
            dr = grid(gP)%voronoi(jcell)%r - rVec
            r2sqd = dot_product(dr,dr)
            dSaux = r2sqd - r1sqd
            
            ! Check special case
            if (dSaux < 0.0) then
               invr12comp = 1.0/r12comp
               dS = 0.0
               inext = jcell
               jjmin = jj
               exit
            end if
            
            if (r12comp*dSaux > 0.0 .and. dSaux < dS*v%connect(jj)%length*r12comp) then
               invr12comp = 1.0/r12comp
               dS = dSaux*v%connect(jj)%invlength*invr12comp
               inext = jcell
               jjmin = jj
            end if

         end do
         !---------------------------------------------------------------------
         
         ! Set path-length fractionally beyond cell boundary to ensure boundary 
         ! is crossed by photon packet
         dS = 0.5*dS
         dS = 1.00001*dS
         dV = v%volume

         
         ! Check if we've reached the edge of the domain
         if (enPacket%STot + dS >= enPacket%SEscape) then
            dS = enPacket%SEscape - enPacket%STot
            inext = -1
         end if
         
         
         ! calculate the optical depth to the next cell wall 
         tauCell = dS*grid(gP)%opacity(grid(gP)%activeV(iVoronoi), enPacket%nuP)
         
         ! check if photon packet interacets within this cell
         !---------------------------------------------------------------------
         if ((absTau + tauCell > totTau) .and. grid(gP)%activeV(iVoronoi) > 0) then

            ! calculate where within this cell the packet is absorbed
            dlLoc = (totTau - absTau)/&
                 &grid(gP)%opacity(grid(gP)%activeV(iVoronoi),enPacket%nuP)
            
            ! update packet's position
            rVec = rVec + dlLoc*vHat
            
            ! add contribution of the packet to the radiation field (BOTH THE SAME!?)

            grid(gP)%Jste(grid(gP)%activeV(iVoronoi),enPacket%nuP) = &
                 & grid(gP)%Jste(grid(gP)%activeV(iVoronoi),enPacket%nuP) + dlLoc*deltaE(iStar) / dV

               
            ! check if the position within the cell is still withing the Outer radius
            ! (Not sure how this applies to Voronoi yet; Ask Barbara)
            !------------------------------------------------------------------
            !if (dlLoc > dS .and. inext == -1) then

            if ( sqrt(dot_product(rvec,rvec)) >= R_out .and. R_out > 0.0) then
               call photonPacketEscapes(enPacket)
               return
            end if
            !------------------------------------------------------------------
            
            
            ! check if the packet is absorbed or scattered
            !------------------------------------------------------------------
            if (lgDust) then
               
               probSca = grid(gP)%scaOpac(grid(gP)%activeV(iVoronoi),enPacket%nuP)/&
                    & (grid(gP)%opacity(grid(gP)%activeV(iVoronoi),enPacket%nuP))
               
               call random_number(random)
!               random = random2
               random = 1.0 - random
               
               if (random > probSca) then
                  lgScattered = .false.         
               else if (random <= probSca) then
                  lgScattered = .true.         
               else
                  print*, '! pathSegment: insanity occured at the scattering/absorption&
                       & decision stage.'
                  print*,'scaOpac,opacity,Ndust'
                  print*,&
                       &grid(gP)%scaOpac(grid(gP)%activeV(iVoronoi),enPacket%nuP),&
                       &grid(gP)%opacity(grid(gP)%activeV(iVoronoi),enPacket%nuP),&
                       &grid(gP)%Ndust(grid(gP)%activeV(iVoronoi))
                  print*,'random number=',random,'probSca=',probSca
                  stop
               end if
               
               
               !---------------------------------------------------------------
               if (.not. lgScattered) then
                  
                  absInt = absInt + 1.
                  
                  ! packet is absorbed by the dust
                  if (.not.lgGas) then
                     packetType = "dustEmi"
                     exit                         
                     
                  ! packet is absorbed by the dust+gas
                  else
                     packetType = "diffuse"
                     exit    
                     
                  end if
                  
                  
               !---------------------------------------------------------------
               else
                  
                  scaInt = scaInt + 1.                           
                  
                  if (lgMultiDustChemistry) then
                     do nS = 1, nSpeciesPart(grid(gP)%dustAbunIndex(grid(gP)%activeV(iVoronoi)))
                        if (grainabun(grid(gP)%dustAbunIndex(grid(gP)%activeV(iVoronoi)),nS)>0. &
                             &.and. grid(gP)%Tdust(nS, 0, & 
                             & grid(gP)%activeV(iVoronoi))<TdustSublime(dustComPoint(&
                             &grid(gP)%dustAbunIndex(grid(gP)%activeV(iVoronoi)))-1+nS)) exit
                     end do
                     if (nS>nSpeciesPart(grid(gP)%dustAbunIndex(grid(gP)%activeV(iVoronoi)))) then
                        print*, "! pathSegment: packet scatters with dust at position where all &
                             &grains have sublimed -1-."
                        print*, iVoronoi, grid(gP)%activeV(iVoronoi), tauCell, absTau, totTau
                        stop
                     end if
                  else
                     do nS = 1, nSpeciesPart(1)
                        if (grainabun(1,nS)>0. &
                             &.and. grid(gP)%Tdust(nS, 0, & 
                             & grid(gP)%activeV(iVoronoi))<TdustSublime(dustComPoint(&
                             &1)-1+nS)) exit
                     end do
                     if (nS>nSpeciesPart(1)) then
                        print*, "! pathSegment: packet scatters with dust at position where all &
                             &grains have sublimed -2-."
                        print*, iVoronoi, grid(gP)%activeV(iVoronoi), tauCell, absTau, totTau
                        stop
                     end if
                  end if
                  
                  ! packet is scattered by the grain                         
                  ! calculate new direction
                  enPacket%iVoronoi = iVoronoi
                  rVec2%x = rVec(1)
                  rVec2%y = rVec(2)
                  rVec2%z = rVec(3)
             
                  
                  if (grid(gP)%activeV(iVoronoi) < 0) then
                     print*, "! pathSegment: new packet cannot be emitted from re-mapped cell -1-"
                     print*, "nuP, starPosition(iStar), .false., .true., iVoronoi"
                     print*, nuP, starPosition(iStar), .false., .true.,  iVoronoi
                     stop
                  end if
                  
                  enPacket = initPhotonPacketV(enPacket%nuP, rVec2, &
                       & enPacket%direction, .false., .false., iVoronoi, .true.)
                  
                  if (.not.lgIsotropic .and. .not.enPacket%lgStellar) then
                     do ihg = 1,10
                        call hgV(enPacket,iierr)
                        if(iierr==0) exit
                     end do
                     if (ihg >=10) then
                        print*, '! pathSegment [warning]: hg called ten times [enPacket]', enPacket
                     end if
                  end if
                  
                  vHat(1) = enPacket%direction%x
                  vHat(2) = enPacket%direction%y
                  vHat(3) = enPacket%direction%z
                  
                  
                  if (.not. (enPacket%direction%x >= 0. .or. enPacket%direction%x < 0)) then
                     print*, '! pathSegment: [1] insane direction%x [direction%x]'
                     print*, enPacket%direction%x
                     stop
                  end if
                  
                  
                  ! initialize optical depth for next photon packet
                  absTau = 0.0
                  call random_number(random)
!                  random = random2

                  totTau = -log(1.0 - random)
                  
               end if
               !---------------------------------------------------------------
               
               
            else
               
               absInt = absInt + 1.
               
               if (.not.lgGas) then
                  print*, "! pathSegment: Insanity occured - no gas present when no dust interaction"
                  stop
               end if
               
               ! packet interacts with gas
               packetType = "diffuse"
               exit
               
            end if
            !------------------------------------------------------------------
            
            
         ! the packet is not absorbed within this cell;
         ! add contribution of the packet to the radiation field
         !---------------------------------------------------------------------
         else
            
            grid(gP)%Jste(grid(gP)%activeV(iVoronoi),enPacket%nuP) = &
                 & grid(gP)%Jste(grid(gP)%activeV(iVoronoi),enPacket%nuP) + &
                 & dS*deltaE(iStar) / dV


            ! update absTau
            absTau = absTau + tauCell

            ! update packet's position
            rVec = rVec + dS*vHat
            enPacket%STot = enPacket%STot + dS

            ! update cell index
            iprev = iVoronoi
            iVoronoi = inext
            if (iVoronoi /= -1) then
               v => grid(gp)%voronoi(iVoronoi)
            end if


            ! be 6/6/06
            !------------------------------------------------------------------
            if(.not.lgPlaneIonization.and..not.lgSymmetricXYZ) then
               lgReturn=.false.

               ! DAVID : Some extra code to add here in future (ask Barbara)
               ! (Update : Also not needed I think)

               if (inext == -1) then
                  lgReturn = .true.
                  call photonPacketEscapes(enPacket)
                  return
               end if
               
            end if
            !------------------------------------------------------------------


            ! DAVID : Code here for lgPlaneIonisation, for future
            if (lgPlaneIonization) then
               lgReturn = .false.

               print*, ' pathSegmentV: plane ionisation not yet implemented for Voronoi'
               stop

               ! DAVID : Some extra code to add here in future (ask Barbara)

               if (lgReturn) then
                  call photonPacketEscapes(enPacket)
                  return
               end if
               
            end if
            !------------------------------------------------------------------

            
            ! check if path is still within the simulation region
            radius = sqrt(dot_product(rVec,rVec))

! BE 25.4.14
!            if (.not. lgPlaneIonization) then
!
!               if (iVoronoi == -1) then
!                  call photonPacketEscapes(enPacket)
!                  return
!               end if
!
!            end if


            ! DAVID : stuff here also
            if (lgSymmetricXYZ) then
               !print*, '! pathSegmentV: SymmetricXYZ not yet implemented in Voronoi'
               !stop

               if (inext == -1) then
                  if (enPacket%iboundary == 1) then
                     vHat(1) = -vHat(1)
                     iVoronoi = iprev
                     inext = iprev
                  else if (enPacket%iboundary == 2) then
                     vHat(2) = -vHat(2)
                     iVoronoi = iprev
                     inext = iprev
                  else if (enPacket%iboundary == 3) then
                     vHat(3) = -vHat(3)
                     iVoronoi = iprev
                     inext = iprev
                  end if

                  ! Make sure photon packet remains in bounding box
                  if (rVec(1) < boxXn) rVec(1) = rVec(1) + 2.0*(boxXn - rVec(1))
                  if (rVec(2) < boxYn) rVec(2) = rVec(2) + 2.0*(boxYn - rVec(2))
                  if (rVec(2) < boxZn) rVec(3) = rVec(3) + 2.0*(boxZn - rVec(3))

                  ! If packet has been reflected, recalculate next boundary values
                  if (enPacket%iboundary /= 0) then
                     rVec2%x = rVec(1)
                     rVec2%y = rVec(2)
                     rVec2%z = rVec(3)
                     vHat2%x = vHat(1)
                     vHat2%y = vHat(2)
                     vHat2%z = vHat(3)
                     call calculateBoundaryPath(vHat2,rVec2,Smin,iboundary)

                     enPacket%STot = 0.0
                     enPacket%SEscape = Smin
                     
                     if (lgSymmetricXYZ) then
                        enPacket%iboundary = iboundary
                     else
                        enPacket%iboundary = 0
                     end if
                  end if

               end if

            end if
            
         end if
         !---------------------------------------------------------------------
         
! 3.4.
         if (inext == -1) then
            call photonPacketEscapes(enPacket)
            return
         end if
     

    
      end do
      !------------------------------------------------------------------------


      if (i >= safeLimit) then
         if (.not.lgPlaneIonization) then
            print*, '! pathSegment: [warning] packet trajectory has exceeded&
                 &  maximum number of events', safeLimit, iVoronoi, &
                 & grid(gP)%activeV(iVoronoi), rvec, vhat, enPacket, iphot
            !stop
         end if
         return
         
      end if
      
      ! the energy packet has been absorbed - reemit a new packet from this position
      enPacket%iVoronoi = iVoronoi
      chTypeIn          = packetType
      positionIn%x      = rVec(1)
      positionIn%y      = rVec(2)
      positionIn%z      = rVec(3)
      inV               = enPacket%iVoronoi
      gPIn              = 1
      reRun             = 1

      return

    end subroutine pathSegmentV



    !--------------------------------------------------------------------------
    subroutine photonPacketEscapes(enPacket)
      implicit none

      type(photon_packet),intent(in) :: enPacket  ! escaping energy packet

      integer :: idirT,idirP                      ! anglular bin variables (theta and phi)

      ! the packet escapes without further interaction
      ! the packet escapes without further interaction
      if (lgSymmetricXYZ) then
         idirT = int(acos(abs(enPacket%direction%z))/dTheta)+1
      else
         idirT = int(acos(enPacket%direction%z)/dTheta)+1
      end if
      
      !                      idirT = int(acos(enPacket%direction%z)/dTheta)+1
      if (idirT>totangleBinsTheta) then
         idirT=totangleBinsTheta
      end if
      if (idirT<1 .or. idirT>totAngleBinsTheta) then
         print*, '! pathSegment: error in theta direction cosine assignment',&
              &  idirT, enPacket, dTheta, totAngleBinsTheta
         stop
      end if
      
      
      if (abs(enPacket%direction%x)<1.e-35) then
         idirP = 0
      else
         if (lgSymmetricXYZ) then
            idirP = int(atan(abs(enPacket%direction%y)/abs(enPacket%direction%x))/dPhi)
         else
            idirP = int(atan(enPacket%direction%y/enPacket%direction%x)/dPhi)
         end if
      end if
      if (idirP<0) idirP=totAngleBinsPhi+idirP
      
      idirP=idirP+1
      
      if (idirP>totangleBinsPhi) then
         idirP=totangleBinsPhi
      end if
      if (idirP<1 .or. idirP>totAngleBinsPhi) then
         print*, '! photonPacketEscapes: error in phi direction cosine assignment -3',&
              &  idirP, enPacket, dPhi, totAngleBinsPhi
         stop
      end if
      
      if (nAngleBins>0) then
         if (viewPointPtheta(idirT) > 0 .and. viewPointPhi(viewPointPtheta(idirT)) < 0) then
            
            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,viewPointPtheta(idirT)) = &
                 &grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                 & enPacket%nuP,viewPointPtheta(idirT)) + deltaE(iStar)                     
            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), & 
                 & enPacket%nuP,0) = &
                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                 & enPacket%nuP,0) +  deltaE(iStar)
            
         elseif (viewPointPtheta(idirT) == viewPointPphi(idirP).or. &
              & (viewPointTheta(viewPointPphi(idirP))==viewPointTheta(viewPointPtheta(idirT))) .or. &
              & (viewPointPhi(viewPointPtheta(idirT))==viewPointPhi(viewPointPphi(idirP))) ) then
            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), enPacket%nuP,& 
                 & viewPointPtheta(idirT)) =grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2),& 
                 & enPacket%nuP,viewPointPtheta(idirT)) +  deltaE(iStar)
            
            if (viewPointPtheta(idirT)/=0) grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                 &enPacket%nuP,0) = &
                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                 & enPacket%nuP,0) +  deltaE(iStar)
            
         else
            grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                 & enPacket%nuP,0) = &
                 & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
                 & enPacket%nuP,0) +  deltaE(iStar)
            
         end if
         
      else
         
         grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
              enPacket%nuP,0) = &
              & grid(enPacket%origin(1))%escapedPackets(enPacket%origin(2), &
              & enPacket%nuP,0) +  deltaE(iStar)
      end if

      return
    end subroutine photonPacketEscapes



    !--------------------------------------------------------------------------
    subroutine hgV(inpacket,ierr)
      
      implicit none
      
      type(photon_packet), intent(inout) :: inpacket
      
      real :: nxp,nyp,nzp
      real :: sint,cost,sinp,cosp,phi
      real :: costp,sintp,phip
      real :: bmu,b,ri1,ri3,cosi3,sini3,cosb2,sinbt,sini2,bott,cosdph
      real :: cosi2,sin2i3,sin2i2,cos2i3,cos2i2,sin2,cos2,sin2cos1
      real :: cos2sin1,cosi1,sini1,sin2i1,cos2i1
      real :: random, random0, vin(3), vout(3), s, denom
      real :: hgg, g2, znorm
      
      integer :: ierr
      
      ierr = 0
      
      ! BEKS 2010
      vin(1) = inpacket%direction%x
      vin(2) = inpacket%direction%y
      vin(3) = inpacket%direction%z
      
      hgg = gSca(inpacket%nuP)
      
      ! henyey-greenstein http://robertoreif.com/Documents/Chapter3.pdf
      call random_number(random0) ! this is the theta dependence
      s=2.*random0-1.
      if (hgg.ge.0.0001) then
         cost=0.5/hgg*(1.+hgg**2-((1.-hgg**2)/(1.+hgg*s))**2)
      else
         cost=s !+1.5*hgg*(1.-s**2)-2*hgg**2*s*(1.-s**2)
      endif
      if (cost .ge. 1.0) then 
         cost=1.0
         sint=0.
      elseif (cost .lt. -1.0) then 
         cost=-1.0
         sint=0.
      else
         sint=sqrt(1.-cost**2)
      endif
      
      call random_number(random) ! this is the phi dependence 
      phi = twoPi*random 
      cosp=cos(phi)
      sinp=sin(phi)
      denom=sqrt(1.-vin(3)**2)
      
      if (denom.gt.0.001) then
         vout(1)=sint/denom*(vin(1)*vin(3)*cosp-vin(2)*sinp) + vin(1)*cost
         vout(2)=sint/denom*(vin(2)*vin(3)*cosp+vin(1)*sinp) + vin(2)*cost
         vout(3)=-sint*cosp*denom+vin(3)*cost
      else
         vout(1)=sint*cosp
         vout(2)=sint*sinp
         if (vin(3).ge.0.) then
            vout(3)=cost
         else
            vout(3)=-cost
         endif
      endif
      
      
      
      if ( ( (abs(vout(1)) <= 1.) .and. (abs(vout(2)) <= 1.) .and. (abs(vout(3)) <= 1.) )&
           & .and. (vout(1) >= 0. .or. vout(1) < 0.) .and. &
           & (vout(2) >= 0. .or. vout(2) < 0.) .and. & 
           & (vout(3) >= 0. .or. vout(3) < 0.)) then
         
         inpacket%direction%x =  vout(1)
         inpacket%direction%y =  vout(2)
         inpacket%direction%z =  vout(3)
         
         
         ierr = 0
      elseif ((abs(vout(1))>=1.).or.(abs(vout(2))>=1.).or.(abs(vout(3))>=1.)&
           & .and. (vout(1) >= 0. .or. vout(1) < 0.) .and. &
           & (vout(2) >= 0. .or. vout(2) < 0.) .and. & 
           & (vout(3) >= 0. .or. vout(3) < 0.)) then
         znorm=sqrt(vout(1)**2 + vout(2)**2 + vout(3)**2)
         vout=vout/znorm
      else
         ierr = 1
         print*,'FAIL',random0,random,s
         print*,' ',vin
         print*,' ',cost,sint
         print*,' ',cosp,sinp,denom
         print*,' ',vout
         
         vout(1) = vin(1)
         vout(2) = vin(2)
         vout(3) = vin(3)
         
         !      print*, vout
         !      print*, vin
         !      print*, random0, random
         !      print*, hgg
         !      print*, inpacket%nuP
         !      print*, sint, cost, sinp, cosp
         
      end if
      
      return
      
                        
    end subroutine hgV

  end subroutine energyPacketDriverV
  
  
end module voronoi_photon_mod
