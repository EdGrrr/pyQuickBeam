  subroutine optical_sphere(freq,D,nsizes,tt,ice,rho_e,qext,qbsca)
  use math_lib
  use optics_lib, only : mieint, m_wat, m_ice
  implicit none
  
! Purpose:
!   Returns optical properties of spheres using Mie calculations
!   Part of QuickBeam v1.1a by John Haynes
!   http://cloudsat.atmos.colostate.edu/radarsim
!
! Inputs:
!   [freq]      radar frequency (GHz)
!   [D]         discrete drop sizes (um)
!   [nsizes]    number of discrete drop sizes
!   [tt]        hydrometeor temperature (C)
!   [ice]       indicates volume consists of ice
!   [rho_e]     medium effective density (kg m^-3) (-1 = pure)
!
! Outputs:
!   [qext]      extinction efficiency
!   [qbsca]     backscatter efficiency
!
! Created:
!   11/28/05  John Haynes (haynes@atmos.colostate.edu)  
  
! ----- INPUTS -----  
  integer, intent(in) :: ice
  integer*4, intent(in) :: nsizes  
  real*8, intent(in) :: freq,D(nsizes),tt,rho_e(nsizes)
  
! ----- OUTPUTS -----
  real*8, intent(out) :: qext(nsizes),qbsca(nsizes)
  
! ----- INTERNAL -----
  integer :: &
  correct_for_rho               ! correct for density flag
  real*8, dimension(nsizes) :: &
  D0, &                         ! D in (m)
  sizep, &                      ! size parameter
  rho_ice, &                    ! bulk density ice (kg m^-3)
  f, &                          ! ice fraction
  beta                          ! Bruggeman factor
  real*8 :: &
  wl                            ! wavelength (m)
  complex*16 :: &
  m                             ! complex index of refraction of bulk form
  complex*16, dimension(nsizes) :: &
  m0, &                         ! complex index of refraction
  e0                            ! complex dielectric factor
  
  integer*4 :: i
  real*8 :: pi
  real*8 :: n_r, n_i, dqv(1), dqsc, dg, dph(1)
  integer*4 :: err
  complex*16 :: Xs1(1), Xs2(1)
  
  pi = acos(-1.0)
  rho_ice(:) = 917

! // conversions
  D0 = D*1E-6                   ! m
  wl = 2.99792458/(freq*10)     ! m
  
! // get the refractive index of the bulk hydrometeors
  if (ice == 0) then
    call m_wat(freq,tt,n_r,n_i)
  else
    call m_ice(freq,tt,n_r,n_i)
  endif
  m = cmplx(n_r,-n_i)
  m0(:) = m
    
  correct_for_rho = 0
  if ((ice == 1) .and. (minval(rho_e) >= 0)) correct_for_rho = 1
    
! // correct refractive index for ice density if needed
! // using the Bruggeman formula
  if (correct_for_rho == 1) then
    f = rho_e/rho_ice
    e0 = m0**2.
    beta = 3.*f*e0-e0+2.-3.*f
    m0 = (0.25*(beta+sqrt(beta**2.+8.*e0)))**0.5
  endif

! // Mie calculations
  sizep = (pi*D0)/wl
  dqv(1) = 0.
  do i=1,nsizes
    call mieint(sizep(i), m0(i), 1, dqv, qext(i), dqsc, qbsca(i), &
      dg, xs1, xs2, dph, err)
  end do

  end subroutine optical_sphere
