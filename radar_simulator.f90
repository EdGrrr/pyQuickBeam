  subroutine radar_simulator(freq,k2,do_ray,use_gas_abs,use_mie_table,mt, &
    nhclass,hp,nprof,ngate,nsizes,D,hgt_matrix,hm_matrix,mf_matrix, &
    p_matrix,t_matrix,rh_matrix,Ze_non,Ze_ray,h_atten_to_vol, &
    g_atten_to_vol,dBZe)
  
  use array_lib
  use math_lib
  use optics_lib
  use radar_simulator_types
  implicit none
  
! Purpose:
!   Simulates a vertical profile of radar reflectivity
!   Part of QuickBeam v1.1a by John Haynes
!   http://cloudsat.atmos.colostate.edu/radarsim
!
! Inputs:
!   [freq]            radar frequency (GHz), can be anything unless
!                     use_mie_table=1, in which case one of 94,35,13.8,9.6,3
!   [k2]              |K|^2, the dielectric constant, set to -1 to use the
!                     frequency dependent default
!   [do_ray]          1=do Rayleigh calcs, 0=not
!   [use_gas_abs]     1=do gaseous abs calcs, 0=not,
!                     2=use same as first profile (undocumented)
!   [use_mie_table]   1=use Mie tables, 0=not
!   [mt]              Mie look up table
!   [nhclass]         number of hydrometeor types
!   [hp]              structure that defines hydrometeor types
!   [nprof]           number of hydrometeor profiles
!   [ngate]           number of vertical layers
!   [nsizes]          number of discrete particles in [D]
!   [D]               array of discrete particles (um)
!
!   (The following 5 arrays must be in order from closest to the radar
!    to farthest...)
!   [hgt_matrix]      height of hydrometeors (km)
!   [hm_matrix]       table of hydrometeor mixing rations (g/kg)
!   [mf_matrix]       table of melted fractions
!   [p_matrix]        pressure profile (hPa)
!   [t_matrix]        temperature profile (C)
!   [rh_matrix]       relative humidity profile (%)
!
! Outputs:
!   [Ze_non]          radar reflectivity without attenuation (dBZ)
!   [Ze_ray]          Rayleigh reflectivity (dBZ)
!   [h_atten_to_vol]  attenuation by hydromets, radar to vol (dB)
!   [g_atten_to_vol]  gaseous atteunation, radar to vol (dB)
!   [dBZe]            effective radar reflectivity factor (dBZ)
!
! Created:
!   11/28/2005  John Haynes (haynes@atmos.colostate.edu)
! Modified:
!   09/2006  placed into subroutine form, scaling factors (Roger Marchand,JMH)

! ----- INPUTS -----  
  type(mie), intent(in) :: mt
  type(class_param), intent(inout) :: hp
  real*8, intent(in) :: freq,k2
  integer, intent(in) ::  do_ray,use_gas_abs,use_mie_table, &
    nhclass,nprof,ngate,nsizes
  real*8, dimension(nsizes), intent(in) :: D
  real*8, dimension(nprof,ngate), intent(in) :: hgt_matrix, p_matrix, &
    t_matrix,rh_matrix
  real*8, dimension(nhclass,nprof,ngate), intent(in) :: hm_matrix, &
    mf_matrix
    
! ----- OUTPUTS -----
  real*8, dimension(nprof,ngate), intent(out) :: Ze_non,Ze_ray, &
    h_atten_to_vol,g_atten_to_vol,dBZe
    
! ----- INTERNAL -----
  integer :: &
  ns                            ! number of discrete drop sizes

  integer*4, dimension(ngate) :: &
  hydro                         ! 1=hydrometeor in vol, 0=none
  real*8 :: &
  rho_a, &                      ! air density (kg m^-3)
  gases                         ! function: 2-way gas atten (dB/km)

  real*8, dimension(:), allocatable :: &
  Di, &                         ! discrete drop sizes (um)
  Ni, &                         ! discrete concentrations (cm^-3 um^-1)
  qext, &                       ! discrete extinction efficiencies
  qbsca                         ! discrete backscatter efficiencies
  
  real*8, dimension(ngate) :: &
  z_vol, &                      ! effective reflectivity factor (mm^6/m^3)
  z_ray, &                      ! reflectivity factor, Rayleigh only (mm^6/m^3)
  kr_vol, &                     ! attenuation coefficient hydro (dB/km)
  g_vol, &                      ! attenuation coefficient gases (dB/km)
  a_to_vol, &                   ! integrated atten due to hydometeors, r>v (dB)
  g_to_vol                      ! integrated atten due to gases, r>v (dB)

  real*8, dimension(nhclass) :: &
  hm                            ! hydrometeor mixing ratio (g/kg)
  
  real*8 :: kr, ze, zr, pi
  integer*4 :: tp, i, j, k, pr, itt

  pi = acos(-1.0)
  
! // loop over each profile (nprof)
  do pr=1,nprof

  ! ----- calculations for each volume ----- 
    z_vol(:) = 0
    z_ray(:) = 0
    kr_vol(:) = 0
    hydro(:) = 0

  ! // loop over each range gate (ngate)
    do k=1,ngate

      rho_a = (p_matrix(pr,k)*100.)/(287*(t_matrix(pr,k)+273.15))
    
    ! :: if mixing ratios are in (g/m^3), convert to (g/kg)
      hm = hm_matrix(:,pr,k)
      where(hp%cp == 1) hm = hm/rho_a
    
    ! :: determine if hydrometeor(s) present in volume
      hydro(k) = 0
      j = 1
      do while ((hydro(k) == 0) .and. (j <= nhclass))
        if ((hm(j) > 1E-12) .and. (hp%dtype(j) > 0)) hydro(k) = 1
        j = j + 1
      enddo

      if (hydro(k) == 1) then
      ! :: if there is hydrometeor in the volume            

      ! // loop over hydrometeor type
        do tp=1,nhclass

          if ((hm(tp) <= 1E-12) .or. (hp%dtype(tp) < 0)) cycle
          
        ! // setup arrays containing D,N
          select case(hp%dtype(tp))
          case(4)
            ns = 1
            allocate(Di(ns),Ni(ns),qext(ns),qbsca(ns))
            Di = hp%p1(tp)
            Ni = 0.
          case default
            ns = nsizes            
            allocate(Di(ns),Ni(ns),qext(ns),qbsca(ns)) 
            Di = D
            Ni = 0.
          end select

        ! // calculate particle density for distribution
          if ((hp%rho_eff(tp,1) < 0) .and. (hp%phase(tp) == 1)) then

            if (vals_equal(hp%rho(tp),-1)) then     
          ! :: density from size and mass coefficients
              hp%rho_eff(tp,1:ns) = (6/pi)*hp%apm(tp)*(Di*1E-6)** &
                (hp%bpm(tp)-3)
            else
          ! :: density is constant
              hp%rho_eff(tp,1:ns) = hp%rho(tp)
            endif

          ! :: limit density to reasonable values
            where(hp%rho_eff(tp,:) < 5) hp%rho_eff(tp,:) = 5
            where(hp%rho_eff(tp,:) > 917) hp%rho_eff(tp,:) = 917

          ! :: Note that density is used to calulate a modified index of
          ! :: refraction, but the particle mass is unchanged for
          ! :: reflectivity and attenuation calculations. This
          ! :: approximation results in some (small) error but is
          ! :: chosen over the time-consuming alternative of
          ! :: integrating the PSD over finite limits.

          endif
          
        ! // get electromagnetic properties of volume, setup N
          if (mf_matrix(tp,pr,k) < -8.) then

            if (use_mie_table == 1) then
          
            ! :: electromagentic properties from Mie table
              if ((hp%dtype(tp) == 4) .and. (hp%idd(tp) < 0)) then
                hp%idd(tp) = infind(mt%D,Di(1))
              endif
            
              if (hp%phase(tp) == 0) then
                itt = infind(mt_ttl,t_matrix(pr,k))
                select case(hp%dtype(tp))
                case(4)
                  qext(1) = mt%qext(hp%idd(tp),itt,1,mt%iff)
                  qbsca(1) = mt%qbsca(hp%idd(tp),itt,1,mt%iff)
                case default
                  qext = mt%qext(:,itt,1,mt%iff)
                  qbsca = mt%qbsca(:,itt,1,mt%iff)
                end select
              else
                itt = infind(mt_tti,t_matrix(pr,k))
                select case(hp%dtype(tp))
                case(4)
                  if (hp%ifc(tp,1) < 0) then
                    hp%ifc(tp,1) = infind(mt%f,hp%rho_eff(tp,1)/917.)
                  endif               
                  qext(1) = &
                    mt%qext(hp%idd(tp),itt+cnt_liq,hp%ifc(tp,1),mt%iff)
                  qbsca(1) = &
                    mt%qbsca(hp%idd(tp),itt+cnt_liq,hp%ifc(tp,1),mt%iff)              
                case default
                  do i=1,ns
                    if (hp%ifc(tp,i) < 0) then
                      hp%ifc(tp,i) = infind(mt%f,hp%rho_eff(tp,i)/917.)
                    endif             
                    qext(i) = mt%qext(i,itt+cnt_liq,hp%ifc(tp,i),mt%iff)
                    qbsca(i) = mt%qbsca(i,itt+cnt_liq,hp%ifc(tp,i),mt%iff)
                  enddo
                end select
              endif

            else
            
            ! :: electromagentic properties on-the-fly
              call optical_sphere(freq,Di,ns,t_matrix(pr,k),hp%phase(tp), &
                hp%rho_eff(tp,1:ns),qext,qbsca)
                
            endif

          ! :: create DSD of hydrometeors
            call dsd(hm(tp),Di,ns,hp%dtype(tp),rho_a, &
              t_matrix(pr,k),hp%dmin(tp),hp%dmax(tp),hp%rho(tp),hp%p1(tp), &
              hp%p2(tp),hp%p3(tp),hp%fc(tp,1:ns),hp%scaled(tp),hp%apm(tp), &
              hp%bpm(tp),Ni)
          
          else
          
          ! :: electromagnetic properties from melting layer model
            call optical_melt()

          ! :: create DSD of hydrometeors from melting layer model
            call dsd_melt()

          endif

        ! // calculate effective reflectivity factor of volume    
          call zeff(freq,Di,Ni,ns,k2,do_ray,qext,qbsca,ze,zr,kr)          
          kr_vol(k) = kr_vol(k) + kr
          z_vol(k) = z_vol(k) + ze
          z_ray(k) = z_ray(k) + zr
          
          deallocate(Di,Ni,qext,qbsca)
        
        enddo   ! end loop of tp (hydrometeor type)

      else
      ! :: volume is hydrometeor-free

        kr_vol(k) = 0
        z_vol(k) = -999
        z_ray(k) = -999
        
      endif

    ! // attenuation due to hydrometeors between radar and volume
      a_to_vol(k) = 2*path_integral(kr_vol,hgt_matrix(pr,:),1,k-1)
      
    ! // attenuation due to gaseous absorption between radar and volume
      if ( &
         (use_gas_abs == 1) .or. &
         ((use_gas_abs == 2) .and. (pr == 1)) &
         ) &
      then
        g_vol(k) = gases(p_matrix(pr,k),t_matrix(pr,k)+273.15, &
          rh_matrix(pr,k),freq)
        g_to_vol(k) = path_integral(g_vol,hgt_matrix(pr,:),1,k-1)
      elseif (use_gas_abs == 0) then
        g_to_vol(k) = 0
      endif      
      
    ! // store results in matrix for return to calling program
      h_atten_to_vol(pr,k)=a_to_vol(k)
      g_atten_to_vol(pr,k)=g_to_vol(k)
      if ((do_ray == 1) .and. (z_ray(k) > 0)) then
        Ze_ray(pr,k) = 10*log10(z_ray(k))
      else
        Ze_ray(pr,k) = -999
      endif
      if (z_vol(k) > 0) then
        dBZe(pr,k) = 10*log10(z_vol(k))-a_to_vol(k)-g_to_vol(k)
        Ze_non(pr,k) = 10*log10(z_vol(k))
      else
        dBZe(pr,k) = -999
        Ze_non(pr,k) = -999
      endif
      
    enddo       ! end loop of k (range gate)
  enddo         ! end loop over pr (profile)  

  end subroutine radar_simulator
  
