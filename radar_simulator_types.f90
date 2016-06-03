  module radar_simulator_types

! Collection of common variables and types
! Part of QuickBeam v1.1a by John Haynes
! http://cloudsat.atmos.colostate.edu/radarsim

  integer, parameter ::       &
  maxhclass = 50             ,& ! max number of hydrometeor classes
  nd = 85                       ! number of discrete particles  

  real*8, parameter ::        &
  dmin = 0.1                 ,& ! min size of discrete particle
  dmax = 10000.                 ! max size of discrete particle
   
! ---- hydrometeor class type -----  
  
  type class_param
    real*8, dimension(maxhclass) :: p1,p2,p3,dmin,dmax,apm,bpm,rho
    integer, dimension(maxhclass) :: dtype,col,cp,phase
    logical, dimension(maxhclass) :: scaled
    real*8, dimension(maxhclass,nd) :: fc, rho_eff
    integer, dimension(maxhclass,nd) :: ifc
    integer, dimension(maxhclass) :: idd
  end type class_param

! ----- mie table variables -----
  
  integer, parameter :: &
  mt_nfreq = 5              , &
  mt_ntt = 39               , & ! num temperatures in table
  mt_nf = 14                , & ! number of ice fractions in table  
  mt_nd = 85                    ! num discrete mode-p drop sizes in table

  type mie
    real*8 :: freq(mt_nfreq), tt(mt_ntt), f(mt_nf), D(mt_nd)
    real*8, dimension(mt_nd,mt_ntt,mt_nf,mt_nfreq) :: qext, qbsca
    integer :: phase(mt_ntt)
    integer*4 :: iff
  end type mie

  real*8, dimension(:), allocatable :: &
    mt_ttl, &                   ! liquid temperatures (C)
    mt_tti                      ! ice temperatures (C)

  integer*4 :: &
    cnt_liq, &                  ! liquid temperature count
    cnt_ice                     ! ice temperature count

  end module radar_simulator_types
