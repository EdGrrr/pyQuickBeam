! QuickBeam v1.1a
! Created 2005/11/28 John Haynes
! 2006/08/23  placed into subroutine form (Roger Marchand)
! v1.1a released 2008/09/10
!
! haynes@atmos.colostate.edu
! http://cloudsat.atmos.colostate.edu/radarsim

  program driver
  use array_lib
  use atmos_lib
  use format_input
  use radar_simulator_types
  implicit none

! ----- variables -----  

  integer :: & 
  sfc_radar, &                  ! 1=surface radar, 0=spaceborne
  do_ray, &                     ! 1=do Rayleigh calcs, 0=not
  use_gas_abs, &                ! 1=do gaseous abs calcs, 0=not
  use_mie_table, &              ! 1=use Mie tables, 0=not
  nhclass, &                    ! num of hydrometeor classes
  nlev, &                       ! num of vertical layers in sounding profile
  nsizes, &                     ! num of discrete drop sizes
  sonde_format, &               ! 1=separate, 2=combined, 3=MLS, 4=MLW, 5=TRO
  output_format, &              ! 1=full output, 2=simple output
  output_disp, &                ! 1=output to display, 0=file only
  melt_lay, &                   ! 1=activate melting layer model, 0=not
  nprof, &                      ! num of hydrometeor profiles
  ngate                         ! num of vertical layers in hydrometeor prof

  real*8 :: &
  freq, &                       ! radar frequency (GHz)
  k2                            ! |K|^2, -1=use frequency dependent default
  
  real*8, dimension(:), allocatable :: &
  D                             ! array of discrete drop sizes (um)
  
  real*8, dimension(:,:), allocatable :: &
  env_hgt_matrix, &             ! height (km)
  env_p_matrix, &               ! pressure (hPa)
  env_t_matrix, &               ! temperature (C)
  env_rh_matrix, &              ! relative humidity (%)
  hgt_matrix, &                 ! as above, gridded to radar res (km)
  t_matrix, &                   ! as above, gridded to radar res (hPa)
  p_matrix, &                   ! as above, gridded to radar res (C)
  rh_matrix, &                  ! as above, gridded to radar res (%)
  Ze_non, &                     ! radar reflectivity withOUT attenuation (dBZ)
  Ze_ray, &                     ! Rayleigh reflectivity (dBZ)
  h_atten_to_vol, &             ! attenuation by hydromets, radar to vol (dB)
  g_atten_to_vol, &             ! gaseous atteunation, radar to vol (dB)
  dBZe                          ! effective radar reflectivity factor (dBZ)
  
  real*8, dimension(:,:,:), allocatable :: &
  hm_matrix, &                  ! hydrometeor mixing ratio (g kg^-1)
  mf_matrix                     ! melted fraction
  
  character*200 :: &
  sonde_file, &                 ! sounding file
  hydromet_file, &              ! hydrometeor input file
  melt_file, &                  ! melted fraction file
  results_file, &               ! results file
  hclass_file, &                ! hydrometeor classes file
  mie_table_name                ! Mie table name  
                 
  type(mie) ::          mt
  type(class_param) ::  hp
 
  real*8, dimension(:), allocatable :: prof_in
  real*8  :: o1, o2, o3, o4, o5, delt, deltp
  integer :: i1,i2,i3
  integer, parameter :: one = 1
  logical :: fexist, hgt_reversed
  integer :: pr,i,k
 
! ----- main program settings ------
  
! // does a settings file exists -- if so load in values
  inquire(file='settings.dat',exist=fexist)
  if (fexist) then      
    
    open(9,file='settings.dat',action='read')
    read(9,*) i, freq
    read(9,*) i, sfc_radar
    read(9,*) i, use_mie_table
    read(9,*) i, use_gas_abs
    read(9,*) i, sonde_format
    read(9,*) i, do_ray
    read(9,*) i, melt_lay
    read(9,*) i, output_format
    read(9,*) i, output_disp
    read(9,*) i, k2
    read(9,*) i, sonde_file
    read(9,*) i, hydromet_file
    read(9,*) i, hclass_file
    read(9,*) i, results_file
    read(9,*) i, mie_table_name
    close(9)
    
  ! // the melting layer model is not yet operational
    melt_file = 'meltf.dat'
    melt_lay = 0

  else
  
    print *, 'Error: settings.dat does not exist'
    stop
    
  endif
  
! ----- Mie tables -----

  if (use_mie_table == 1) call load_mie_table(mie_table_name,freq,mt)

! ----- hydrometeor class parameters -----
  
  call load_hydrometeor_classes(hclass_file,hp,nhclass)
  allocate(prof_in(nhclass))
  
! ----- create array of discrete drop sizes -----
  
  if (use_mie_table == 1) then
 
  ! :: D specified by table
    nsizes = mt_nd
    allocate(D(nsizes))
    D = mt%D
    
  else
  
  ! :: D created on a log-linear scale
    nsizes = nd
    delt = (log(dmax)-log(dmin))/(nsizes-1)
    deltp = exp(delt)
    allocate(D(nsizes))
    D(1) = dmin
    do i=2,nsizes
      D(i) = D(i-1)*deltp
    enddo   
    
  endif  

! ----- read sounding data if separate from hydrometeor profile -----

  if (sonde_format == 1) then
        
  ! :: read environmental sounding
    open(10,file=sonde_file,action='read')
    read(10,*) nlev

    allocate(env_hgt_matrix(one,nlev))
    allocate(env_p_matrix(one,nlev))
    allocate(env_t_matrix(one,nlev))
    allocate(env_rh_matrix(one,nlev))
                        
    do i=1,nlev
      read(10,*) env_hgt_matrix(one,i), &
        env_p_matrix(one,i), &
        env_t_matrix(one,i), &
        env_rh_matrix(one,i)
    enddo
        
    close(10)

  elseif (sonde_format >= 3) then
        
  ! :: get a standard profile
    nlev = 33
    allocate(env_hgt_matrix(one,nlev))
    allocate(env_p_matrix(one,nlev))
    allocate(env_t_matrix(one,nlev))
    allocate(env_rh_matrix(one,nlev))
        
    call mcclatchey(sonde_format-2,env_hgt_matrix,env_p_matrix, &
      env_t_matrix,env_rh_matrix)
      
  endif  
  
! ----- get number of profiles -----
  
  open(11,file=hydromet_file,action='read')
  if (melt_lay == 1) open(14,file=melt_file,action='read')

! // read in number of profiles
  read(11,*) nprof  
  
! ----- loop over profiles -----

  do pr=1,nprof
  
    read(11,*) ngate
    if (melt_lay == 1) read(14,*) i1
    allocate(hgt_matrix(one,ngate))
    allocate(t_matrix(one,ngate),p_matrix(one,ngate),rh_matrix(one,ngate))
    allocate(hm_matrix(nhclass,one,ngate),mf_matrix(nhclass,one,ngate))

  ! ----- read hydrometeor data (and sounding data if on same grid) -----

    if (sonde_format /= 2) then

    ! // read from file containing only hydrometeor data

      do i=1,ngate
        read(11,*) hgt_matrix(one,i), prof_in
        hm_matrix(:,one,i) = prof_in
      enddo
      
    ! // sounding to radar height conversion

      call irreg_to_grid(hgt_matrix,t_matrix,p_matrix,rh_matrix, &
        env_hgt_matrix,env_t_matrix,env_p_matrix,env_rh_matrix)
      
    else
        
    ! // read from file with both hydrometeor and sounding data
      
      do i=1,ngate
        read(11,*) hgt_matrix(one,i), p_matrix(one,i), &
          t_matrix(one,i), rh_matrix(one,i), prof_in
        hm_matrix(:,one,i) =  prof_in
      enddo
      
    endif
    
    hgt_matrix = hgt_matrix/1000.
    if (maxval(t_matrix) > 75) t_matrix = t_matrix-273.15

  ! ----- read melted fraction data (if used) -----

    if (melt_lay == 1) then
      read(14,*) i1
      do i=1,ngate
        read(14,*) o1, prof_in
          mf_matrix(:,one,i) = prof_in
      enddo
      do i=1,nhclass
        if (hp%col(i) == 0) mf_matrix(i,one,:) = -9.
      enddo
    else
      mf_matrix(:,one,:) = -9.
    endif
    
  ! ----- order data for call to simulator -----

    call order_data(hgt_matrix,hm_matrix,mf_matrix,p_matrix, &
      t_matrix,rh_matrix,sfc_radar,hgt_reversed)

  ! ----- call radar simulator -----
        
    allocate(Ze_non(one,ngate)) 
    allocate(Ze_ray(one,ngate)) 
    allocate(h_atten_to_vol(one,ngate))                 
    allocate(g_atten_to_vol(one,ngate))                 
    allocate(dBZe(one,ngate))
    
  ! :: In this case, the radar_simulator routine is called nprof times,
  ! :: once for each profile. This may be advantageous if an extremely
  ! :: large number of profiles are to be processed since it is not
  ! :: necessary to store all the profiles and their output results
  ! :: in memory simultaneously. In other circumstances it may be
  ! :: advantageous to pass multiple profiles into radar_simulator with
  ! :: a single call to the routine.

    call radar_simulator(freq,k2,do_ray,use_gas_abs,use_mie_table,mt, &
      nhclass,hp,one,ngate,nsizes,D,hgt_matrix,hm_matrix,mf_matrix, &
      p_matrix,t_matrix,rh_matrix,Ze_non,Ze_ray,h_atten_to_vol, &
      g_atten_to_vol,dBZe)
      
  ! ----- BEGIN output section -----
  
    if (pr == 1) then
      open(12,file=results_file,action='write')
      if (output_disp == 1) then
        if (output_format == 1) then
          write(6,101) 'km', 'z_eff', 'z_ray', 'a_to_vol', &
            'g_to_vol', 'z_cor'
        else
          write(6,103) 'km', 'z_cor'
        endif
      endif
      write(12,*) nprof
    endif

    if (hgt_reversed .eqv. .false.) then
      i1=1
      i2=ngate
      i3=1
    else
      i1=ngate
      i2=1
      i3=-1
    endif          
    
    do k=i1,i2,i3
      
      select case(output_format)
      
    ! // full output      
      case(1)
        o1 = Ze_non(one,k)
        o2 = Ze_ray(one,k)
        o3 = h_atten_to_vol(one,k)
        o4 = g_atten_to_vol(one,k)
        o5 = dBZe(one,k)
        if (output_disp == 1) &
        write(6,102) hgt_matrix(one,k), o1, o2, o3, o4, o5
        write(12,102) hgt_matrix(one,k), o1, o2, o3, o4, o5
        
    ! // abridged output
      case(2)
        o5 = dBZe(one,k)
        if (output_disp == 1) &
        write(6,104) hgt_matrix(one,k), o5
        write(12,104) hgt_matrix(one,k), o5
      
      end select
     
    enddo
     
    deallocate(hgt_matrix,hm_matrix,mf_matrix,p_matrix,t_matrix, &
      rh_matrix,Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe)
     
  enddo
  
  if (use_mie_table == 1) deallocate(mt_ttl,mt_tti)

  close(12)       ! close output file
  close(11)       ! close input file
  if (melt_lay == 1) close(14)       ! close melt fraction file

! ----- file formats -----
  101 format (A5,5A10)
  102 format (F5.2,5F10.2)
  103 format (A5,A10)
  104 format (F5.2,F10.2)

  end program driver
