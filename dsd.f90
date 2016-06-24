  subroutine dsd(Q,D,nsizes,dtype,rho_a,tc, &
             dmin,dmax,p1,p2,p3,fc,scaled,apm,bpm,N)
  use array_lib
  use math_lib 
  implicit none

! Purpose:
!   Create a discrete drop size distribution
!   Part of QuickBeam v1.1a by John Haynes
!   http://cloudsat.atmos.colostate.edu/radarsim
!
! Inputs:
!   [Q]        hydrometeor mixing ratio (g/kg)
!   [D]        discrete drop sizes (um)
!   [nsizes]   number of elements of [D]
!   [dtype]    distribution type
!   [rho_a]    ambient air density (kg m^-3)
!   [tc]       temperature (C)
!   [dmin]     minimum size cutoff (um)
!   [dmax]     maximum size cutoff (um)
!   [rho_c]    alternate constant density (kg m^-3)
!      Now removed
!   [p1],[p2],[p3]  distribution parameters
!
! Input/Output:
!   [fc]       scaling factor for the distribution
!   [scaled]   has this hydrometeor type been scaled?
!   [apm]      a parameter for mass (kg m^[-bpm])
!   [bmp]      b params for mass
!
! Outputs:
!   [N]        discrete concentrations (cm^-3 um^-1)
!              or, for monodisperse, a constant (1/cm^3)
!
! Requires:
!   function infind
!
! Created:
!   11/28/05  John Haynes (haynes@atmos.colostate.edu)
! Modified:
!   01/31/06  Port from IDL to Fortran 90
!   07/07/06  Rewritten for variable DSD's
!   10/02/06  Rewritten using scaling factors (Roger Marchand and JMH)
!   11/08/06  Corrected problem with gamma scaling factor (1.03a)
 
! ----- INPUTS -----  
  
  integer*4, intent(in) :: nsizes
  integer, intent(in) :: dtype
  real*8, intent(in) :: Q,D(nsizes),rho_a,tc,dmin,dmax, &
    p1,p2,p3,apm,bpm
    
! ----- INPUT/OUTPUT -----

  logical, intent(inout) :: scaled  
  real*8, intent(inout) :: fc(nsizes)
    
! ----- OUTPUTS -----

  real*8, intent(out) :: N(nsizes)
  
! ----- INTERNAL -----
  
  real*8 :: &
  N0,D0,vu,np,dm,ld, &                  ! gamma, exponential variables
  dmin_mm,dmax_mm,ahp,bhp, &            ! power law variables
  rg,log_sigma_g, &                     ! lognormal variables
  rho_e                                 ! particle density (kg m^-3)
  
  real*8 :: tmp1, tmp2
  real*8 :: pi

  ! Parameters for Field 2005 distribution
  REAL*8 :: mob, mo2, mo3, loga_, a_, b_, cse, mrat, mo0, slam1, slam2
  REAL, DIMENSION(10), PARAMETER:: &
       sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
       0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
  REAL, DIMENSION(10), PARAMETER:: &
       sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
                     0.060366,  0.000079,  0.000594, 0.0,      -0.003577/) 
  REAL*8, PARAMETER :: mu_s = 0.6357
  REAL*8, PARAMETER :: kap0 = 490.6
  REAL*8, PARAMETER :: kap1 = 17.46
  REAL*8, PARAMETER :: lam0 = 20.78
  REAL*8, PARAMETER :: lam1 = 3.29

  ! Parameters for Thompson graupel distribution
  REAL*8, PARAMETER :: mu_g = 0.0
  REAL, PARAMETER :: gonv_max = 3.e6
  REAL, PARAMETER :: gonv_min = 1.e4
  REAL*8 :: cgg1, cgg2, cgg3, cge1, cge2, cge3, &
       org2, lamr, mvd_r, N0_min, &
       ygra1, zans1, N0_exp, lam_exp, lamg, ilamg, &
       N0_g, xslw1, crg3, mu_r, apm_r, bpm_r
     
  pi = acos(-1.0)
  
! // if density is constant, store equivalent values for apm and bpm
  !Now done at readin
  !if ((rho_c > 0) .and. (apm < 0)) then
  !  apm = (pi/6)*rho_c
  !  bpm = 3.
  !endif
  
  select case(dtype)
  
! ---------------------------------------------------------!
! // modified gamma                                        !
! ---------------------------------------------------------!
! :: N0 = total number concentration (m^-3)
! :: np = fixed number concentration (kg^-1)
! :: D0 = characteristic diameter (um)
! :: dm = mean diameter (um)
! :: vu = distribution width parameter

  case(1)  
    if (vals_equal(p1,-1)) then

!     // D0, vu are given    
      dm = p2
      vu = p3
      D0 = gamma(vu)/gamma(vu+1)*dm
      
      if (scaled .eqv. .false.) then
      
        fc = ( &
             ((D*1E-6)**(vu-1)*exp(-1*D/D0)) / &
             (apm*((D0*1E-6)**(vu+bpm))*gamma(vu+bpm)) &
             ) * 1E-12
        scaled = .true.

      endif        
           
      N = fc*rho_a*(Q*1E-3)
    
   elseif (vals_equal(p3, -1)) then
      !     // vu isn't given
      print *, 'Error: Must specify a value for vu'
      stop
      
   else
       !     // N0, vu are given
       if (vals_equal(p2, -2)) then
          !units of p1 is in per m3
          np = p1/rho_a
       else
          np = p1
       endif
       
       N0 = np*rho_a

       if (vals_equal(p3,-2)) then
          !If p3 is -2, use Thompson size parameter - check units! 
          vu = MIN(15., (1000.E6/N0 + 2.)) -1
       elseif (vals_equal(p3, -3)) then
          ! Morrison shape parameter
          vu = (1/(0.0005714*(N0/1.E6)+0.2714))**2 -1
          vu = MIN(MAX(vu,2.),10.)
       else
          vu = p3
       endif
       !write (*,*) vu
      
      !write (*,*) N0
      tmp1 = (Q*1E-3)**(1./bpm)
      
      if (scaled .eqv. .false.) then

        fc = (D*1E-6 / (gamma(vu)/(apm*np*gamma(vu+bpm)))** &
             (1./bpm))**vu
             
        scaled = .true.

      endif

      N = ( &
          (rho_a*np*fc*(D*1E-6)**(-1.))/(gamma(vu)*tmp1**vu) * &
          exp(-1.*fc**(1./vu)/tmp1) &
          ) * 1E-12
             
    endif
    
! ---------------------------------------------------------!
! // exponential                                           !
! ---------------------------------------------------------!
! :: N0 = intercept parameter (m^-4)
! :: ld = slope parameter (um)

  case(2)
    if (vals_differ(p1,-1)) then

!     // N0 has been specified, determine ld
      N0 = p1
      tmp1 = 1./(1.+bpm)
      
      if (scaled .eqv. .false.) then
      
        fc = ((apm*gamma(1.+bpm)*N0)**tmp1)*(D*1E-6)
        scaled = .true.
        
      endif
      
      N = ( &
          N0*exp(-1.*fc*(1./(rho_a*Q*1E-3))**tmp1) &
          ) * 1E-12

    elseif (vals_differ(p2,-1)) then

!     // ld has been specified, determine N0
      ld = p2

      if (scaled .eqv. .false.) then

        fc = (ld*1E6)**(1.+bpm)/(apm*gamma(1+bpm))* &
             exp(-1.*(ld*1E6)*(D*1E-6))*1E-12
        scaled = .true.

      endif

      N = fc*rho_a*(Q*1E-3)

    else

!     // ld will be determined from temperature, then N0 follows
      ld = 1220*10.**(-0.0245*tc)*1E-6
      N0 = ((ld*1E6)**(1+bpm)*Q*1E-3*rho_a)/(apm*gamma(1+bpm))
    
      N = ( &
          N0*exp(-1*ld*D) &
          ) * 1E-12
    
    endif
  
! ---------------------------------------------------------!
! // power law                                             !
! ---------------------------------------------------------!
! :: ahp = Ar parameter (m^-4 mm^-bhp)
! :: bhp = br parameter
! :: dmin_mm = lower bound (mm)
! :: dmax_mm = upper bound (mm)

  case(3)

!   :: br parameter
    if (vals_equal(p1,-2)) then
!     :: if p1=-2, bhp is parameterized according to Ryan (2000),
!     :: applicatable to cirrus clouds
      if (tc < -30) then
        bhp = -1.75+0.09*((tc+273)-243.16)
      elseif ((tc >= -30) .and. (tc < -9)) then
        bhp = -3.25-0.06*((tc+273)-265.66)
      else
        bhp = -2.15
      endif
    elseif (vals_equal(p1,-3)) then
!     :: if p1=-3, bhp is parameterized according to Ryan (2000),
!     :: applicable to frontal clouds
      if (tc < -35) then
        bhp = -1.75+0.09*((tc+273)-243.16)
      elseif ((tc >= -35) .and. (tc < -17.5)) then
        bhp = -2.65+0.09*((tc+273)-255.66)
      elseif ((tc >= -17.5) .and. (tc < -9)) then
        bhp = -3.25-0.06*((tc+273)-265.66)
      else
        bhp = -2.15
      endif    
    else
!     :: otherwise the specified value is used
      bhp = p1
    endif

!   :: Ar parameter
    dmin_mm = dmin*1E-3
    dmax_mm = dmax*1E-3
!   :: commented lines are original method with constant density
!      rc = 500.                ! (kg/m^3)
!      tmp1 = 6*rho_a*(bhp+4)
!      tmp2 = pi*rc*(dmax_mm**(bhp+4))*(1-(dmin_mm/dmax_mm)**(bhp+4))
!      ahp = (Q*1E-3)*1E12*tmp1/tmp2
!   :: new method is more consistent with the rest of the distributions
!   :: and allows density to vary with particle size
    tmp1 = rho_a*(Q*1E-3)*(bhp+bpm+1)
    tmp2 = apm*(dmax_mm**bhp*dmax**(bpm+1)-dmin_mm**bhp*dmin**(bpm+1))
    ahp = tmp1/tmp2
    
    N = ( &
        ahp*(D*1E-3)**bhp &
        ) * 1E-12    

! ---------------------------------------------------------!
! // monodisperse                                          !
! ---------------------------------------------------------!
! :: D0 = particle diameter (um)

  case(4)
  
    if (scaled .eqv. .false.) then
    
      D0 = p1
      rho_e = (6/pi)*apm*(D0*1E-6)**(bpm-3)
      fc(1) = (6./(pi*D0**3*rho_e))*1E12
      scaled = .true.
      
    endif
    
    N(1) = fc(1)*rho_a*(Q*1E-3)
    
! ---------------------------------------------------------!
! // lognormal                                             !
! ---------------------------------------------------------!
! :: N0 = total number concentration (m^-3)
! :: np = fixed number concentration (kg^-1)
! :: rg = mean radius (um)
! :: log_sigma_g = ln(geometric standard deviation)

  case(5)
    if (vals_equal(p1,-1)) then

!     // rg, log_sigma_g are given
      rg = p2
      log_sigma_g = p3
      tmp2 = (bpm*log_sigma_g)**2.

      if (scaled .eqv. .false.) then
            
        fc = 0.5 * ( &
             (1./((2.*rg*1E-6)**(bpm)*apm*(2.*pi)**(0.5) * &
             log_sigma_g*D*0.5*1E-6)) * &
             exp(-0.5*((log(0.5*D/rg)/log_sigma_g)**2.+tmp2)) &
             ) * 1E-12
        scaled = .true.
             
      endif
                
      N = fc*rho_a*(Q*1E-3)
      
    elseif (vals_equal(p2,-1)) then

!     // N0, log_sigma_g are given    
      Np = p1
      log_sigma_g = p3
      N0 = np*rho_a
      tmp1 = (rho_a*(Q*1E-3))/(2.**bpm*apm*N0)
      tmp2 = exp(0.5*bpm**2.*(log_sigma_g))**2.      
      rg = ((tmp1/tmp2)**(1/bpm))*1E6
      
      N = 0.5*( &
        N0 / ((2.*pi)**(0.5)*log_sigma_g*D*0.5*1E-6) * &
        exp((-0.5*(log(0.5*D/rg)/log_sigma_g)**2.)) &
        ) * 1E-12      
      
    else

!     // vu isn't given
      print *, 'Error: Must specify a value for sigma_g'
      stop
    
    endif

  case(6)  
    ! Dual gamma distribution from Field et al. 2005
    ! Adapted from Thompson microphysics code in WRF
     cse = bpm + 1
     
     mob = (Q * 1e-3 * rho_a)/apm

     if (vals_equal(bpm, 2)) then
        mo2 = mob
     else
        loga_ = sa(1) + sa(2)*tc + sa(3)*bpm &
             + sa(4)*tc*bpm + sa(5)*tc*tc &
             + sa(6)*bpm*bpm + sa(7)*tc*tc*bpm &
             + sa(8)*tc*bpm*bpm + sa(9)*tc*tc*tc &
             + sa(10)*bpm*bpm*bpm
        a_ = 10.0**loga_
        b_ = sb(1) + sb(2)*tc + sb(3)*bpm &
             + sb(4)*tc*bpm + sb(5)*tc*tc &
             + sb(6)*bpm*bpm + sb(7)*tc*tc*bpm &
             + sb(8)*tc*bpm*bpm + sb(9)*tc*tc*tc &
             + sb(10)*bpm*bpm*bpm
        mo2 = (mob/a_)**(1./b_)
     endif
     
     loga_ = sa(1) + sa(2)*tc + sa(3)*cse &
          + sa(4)*tc*cse + sa(5)*tc*tc &
          + sa(6)*cse*cse + sa(7)*tc*tc*cse &
          + sa(8)*tc*cse*cse + sa(9)*tc*tc*tc &
          + sa(10)*cse*cse*cse
     a_ = 10.0**loga_
     b_ = sb(1)+ sb(2)*tc + sb(3)*cse + sb(4)*tc*cse &
          + sb(5)*tc*tc + sb(6)*cse*cse &
          + sb(7)*tc*tc*cse + sb(8)*tc*cse*cse &
          + sb(9)*tc*tc*tc + sb(10)*cse*cse*cse
     mo3 = a_ * mo2**b_   
     
     mrat = mo2*(mo2/mo3)*(mo2/mo3)*(mo2/mo3)
     mo0 = (mo2/mo3)**mu_s
     slam1 = (mo2/mo3)*lam0
     slam2 = (mo2/mo3)*lam1

     N = (mrat*kap0*exp(-1*slam1*D*1e-6) &
          +kap1*mrat*mo0* (D*1e-6)**mu_s * exp(-1*slam2*D*1e-6)) * 1e-12


  case(7)
     !Thompson graupel size distribution
     ! p1 - Qrain (g/kg)
     ! p2 - Nrain (kg-1)

     !rain values fixed - pretty unliekly you would want different ones
     ! anyway, as rain is just water....?
     mu_r = 0
     apm_r = 524
     bpm_r = 3
     
     cge1 = bpm+1
     cge2 = mu_g + 1
     cge3 = bpm + mu_g + 1
     cgg1 = gamma(cge1)
     cgg2 = gamma(cge2)
     cgg3 = gamma(cge3)
     crg3 = gamma(bpm_r) 
     org2 = 1./gamma(mu_r+1.)
     
     !Calc mvd_rain from p1 (Qrain), p2 (QNRAIN)
     lamr = (apm_r * crg3 * org2 * (MAX(1e-12, p2*rho_a)/(p1*1e-3*rho_a)))**(1/bpm_r)
     mvd_r = (3.0 + mu_r + 0.672) / lamr

     !write (*,*) "mvd_r ",mvd_r
     
     ! Thompson graupel distribution
     N0_min = gonv_max
     if (tc.lt. -2.5 .and. p1.gt.1e-12 .and. mvd_r.gt.100.E-6) then
        xslw1 = 4.01 + log10(mvd_r)
     else
        xslw1 = 0.01
     endif
     ygra1 = 4.31 + log10(max(5.E-5, Q*1e-3*rho_a))

     zans1 = 3.1 + (100./(300.*xslw1*ygra1/(10./xslw1+1.+0.25*ygra1)+30.+10.*ygra1))
     N0_exp = 10.**(zans1)
     N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
     N0_min = MIN(N0_exp, N0_min)
     N0_exp = N0_min
     !write (*,*) "N0_exp ",N0_exp

     lam_exp = (N0_exp*apm*cgg1/(Q*1e-3*rho_a))**(1./cge1)
     lamg = lam_exp * (cgg3*(1./cgg2)*(1./cgg1))**(1./bpm)
     N0_g = N0_exp/(cgg2*lam_exp) * lamg**cge2

     !write (*,*) "lam_exp ",lam_exp
     !write (*,*) "lamg ",lamg
     !write (*,*) "N0_g ",N0_g
     
     N = N0_g * exp(-1*lamg*(D*1e-6)) * 1e-12
     
  end select
  
  end subroutine dsd
