  subroutine zeff(freq,D,N,nsizes,k2,xr,qext,qbsca,z_eff,z_ray,kr)
  use math_lib
  implicit none
  
! Purpose:
!   Simulates radar return of a volume given DSD of spheres
!   Part of QuickBeam v1.1a by John Haynes
!   http://cloudsat.atmos.colostate.edu/radarsim
!
! Inputs:
!   [freq]      radar frequency (GHz)
!   [D]         discrete drop sizes (um)
!   [N]         discrete concentrations (cm^-3 um^-1)
!   [nsizes]    number of discrete drop sizes
!   [k2]        |K|^2, -1=use frequency dependent default 
!   [xr]        perform Rayleigh calculations?
!   [qext]      extintion efficiency
!   [qbsca]     backscatter efficiency
!
! Outputs:
!   [z_eff]     unattenuated effective reflectivity factor (mm^6/m^3)
!   [z_ray]     reflectivity factor, Rayleigh only (mm^6/m^3)
!   [kr]        attenuation coefficient (db km^-1)
!
! Created:
!   11/28/05  John Haynes (haynes@atmos.colostate.edu)

! ----- INPUTS -----  
  integer, intent(in) :: xr
  integer*4, intent(in) :: nsizes
  real*8, intent(in) :: freq,D(nsizes),N(nsizes),qext(nsizes),qbsca(nsizes)
  real*8, intent(inout) :: k2
  
! ----- OUTPUTS -----
  real*8, intent(out) :: z_eff,z_ray,kr
    
! ----- INTERNAL -----
  real*8, dimension(nsizes) :: &
  D0, &                         ! D in (m)
  N0                            ! N in m^-3 m^-1
  real*8 :: &
  wl, &                         ! wavelength (m)
  cr                            ! kr(dB/km) = cr * kr(1/km)
  
  real*8 :: pi
  real*8 :: eta_sum, eta_mie, const, z0_eff, z0_ray, k_sum

  pi = acos(-1.0)

! // conversions
  D0 = d*1E-6                   ! m
  N0 = n*1E12                   ! 1/(m^3 m)
  wl = 2.99792458/(freq*10)     ! m
  
! // dielectric constant |k^2| defaults
  if (k2 < 0) then
    k2 = 0.933
    if (abs(94.-freq) < 3.) k2=0.75
    if (abs(35.-freq) < 3.) k2=0.88
    if (abs(13.8-freq) < 3.) k2=0.925
  endif  
 
! // eta_mie = 0.25*sum[qbsca*pi*D^2*N(D)*deltaD]
!                   <--------- eta_sum --------->
! // z0_eff = (wl^4/!pi^5)*(1./k2)*eta_mie
  eta_sum = 0.
  if (size(D0) == 1) then
    eta_sum = qbsca(1)*(n(1)*1E6)*D0(1)**2
  else
    call avint(qbsca*N0*D0**2,D0,nsizes,D0(1),D0(size(D0,1)),eta_sum)
  endif
 
  eta_mie = eta_sum*0.25*pi
  const = (wl**4/pi**5)*(1./k2)
  z0_eff = const*eta_mie

! // kr = 0.25*cr*sum[qext*pi*D^2*N(D)*deltaD]
!                 <---------- k_sum --------->  
  k_sum = 0.
  if (size(D0) == 1) then
    k_sum = qext(1)*(n(1)*1E6)*D0(1)**2
  else
    call avint(qext*N0*D0**2,D0,nsizes,D0(1),D0(size(D0,1)),k_sum)
  endif
  cr = 10./log(10.)
  kr = k_sum*0.25*pi*(1000.*cr)
  
! // z_ray = sum[D^6*N(D)*deltaD]
  if (xr == 1) then
    z0_ray = 0.
    if (size(D0) == 1) then
      z0_ray = (n(1)*1E6)*D0(1)**6
    else
      call avint(N0*D0**6,D0,nsizes,D0(1),D0(size(D0)),z0_ray)
    endif
  endif
  
! // convert to mm^6/m^3
  z_eff = z0_eff*1E18 !  10.*alog10(z0_eff*1E18)
  z_ray = z0_ray*1E18 !  10.*alog10(z0_ray*1E18)
  
  end subroutine zeff
