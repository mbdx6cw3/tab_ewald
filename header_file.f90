      MODULE header_file
!     ***********************************************************************
!     This module provides the global variables for the elec_pot.f90 program.
!     Includes a namelist for the input parameters. 
!     ***********************************************************************
      implicit none
      integer, parameter :: dp = SELECTED_REAL_KIND(15,99)

!     main program
!     ---------------------------------------------------------------------
      integer       :: n,natm,i,j,k,p
      real(kind=dp) :: vol, u, rc, rcsq
      real(kind=dp) :: qtot, qO, qSi, qH, qN, qC, qC_1, qH_N, qH_C, qO_s !charges
      real(kind=dp), allocatable, dimension (:) :: cell
      character(LEN=30) :: name
      integer, allocatable, dimension (:) :: aindex
      real(kind=dp), allocatable, dimension (:) :: rxi, ryi, rzi  !adsorbent atom positions
      real(kind=dp), allocatable, dimension (:) :: rx, ry, rz !vector distances
      real(kind=dp) :: vxi, vyi, vzi  !adsorbent atom velocities
      real(kind=dp) :: fxi, fyi, fzi  !adsorbent atom forces
      real(kind=dp) :: rxg, ryg, rzg, grid_res  !grid positions and resolution
      real(kind=dp), allocatable, dimension (:) :: qi, q
      character(LEN=2), allocatable, dimension (:) :: atype  !adsorbent atoms types
      integer :: CK, BK   !config key, boundary key
      real(kind=dp) :: ax, ay, az, lx, lx2   !x,y,z components of a cell vector
      real(kind=dp) :: bx, by, bz, ly, ly2   !x,y,z components of b cell vector
      real(kind=dp) :: cx, cy, cz, lz, lz2   !x,y,z components of c cell vector
      real(kind=dp) :: zint    ! interval between grid points

!     global ewald summation variables
!     ------------------------------------------------------------------------
      integer, parameter :: maxk = 100000
      integer :: elec_switch, kmax1, kmax2, kmax3
      real(kind=dp) :: ewald_neut,ewald_self,recip_sum,real_sum
      real(kind=dp) :: kappa, kcut, ksqcut
      real(kind=dp), allocatable, dimension (:) :: ugrid !tabulated grid

!     parameters
!     ------------------------------------------------------------------------
      real(kind=dp), parameter :: qE = 1.602E-19
      real(kind=dp), parameter :: R = 8.314              ! J K**-1 mol**-1
      real(kind=dp), parameter :: kb = 1.381E-23         ! J K**-1
      real(kind=dp), parameter :: Na = 6.022E23          ! mol**-1
      real(kind=dp), parameter :: pi = 3.1415927 
      real(kind=dp), parameter :: twopi = 6.2831854
      real(kind=dp), parameter :: root_pi = 1.772453851  
      real(kind=dp), parameter :: e0 = 8.854E-12         ! J-1 C2 m-1
      real(kind=dp), parameter :: r4pie0 = 1/ (4.0D0 * e0 * pi) !J m C-2 
      real(kind=dp), parameter :: qEsq = qE * qE         ! C2

      namelist /input_deck/qO,qSi,qH,qN,qH_N,qC,qC_1,qH_C,qO_s,rc,&
                         & grid_res,kappa,kmax1,kmax2,kmax3
!     ************************************************************************
      END MODULE header_file
