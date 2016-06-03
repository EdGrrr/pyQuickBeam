  subroutine load_mie_table(mie_table_name,freq,mt)
  use array_lib
  use radar_simulator_types
  implicit none
  
! Purpose:
!   Loads the Mie table data
!   Part of QuickBeam v1.1a by John Haynes
!   http://cloudsat.atmos.colostate.edu/radarsim
!
! Inputs:  
!   [mie_table_name]   Mie table input file
!
! Outputs:
!   [mt]            structure of Mie table data
!
! Created from QuickBeam v1.02 08/24/2006 by Roger Marchand  

! ----- INPUT -----
  character*200, intent(in) :: mie_table_name
  real*8, intent(in) :: freq

! ----- OUTPUT -----
  type(mie), intent(out) :: mt

! ----- INTERNAL -----  
  integer :: i,tmp
  real*8 :: dst

! // read table  
  open(50,file=mie_table_name,form='unformatted',action='read')
  read(50) tmp,tmp,tmp,tmp ! mt_nfreq, mt_ntt, mt_nf, mt_nd
  read(50) mt%freq, mt%tt, mt%f, mt%phase, mt%D, mt%qext, mt%qbsca
  close(50)
   
! // create arrays of liquid/ice temperature
  cnt_liq = 0
  cnt_ice = 0
  do i=1,mt_ntt
    if (mt%phase(i) == 0) cnt_liq = cnt_liq + 1
    if (mt%phase(i) == 1) cnt_ice = cnt_ice + 1
  enddo
  allocate(mt_ttl(cnt_liq),mt_tti(cnt_ice))
  do i=1,cnt_liq
    mt_ttl(i) = mt%tt(i)
  enddo
  do i=1,cnt_ice
    mt_tti(i) = mt%tt(cnt_liq+i)
  enddo
  
! // determine frequency index
  mt%iff = infind(mt%freq,freq,sort=1,dist=dst)
  if (abs(dst) > 2) then
    print *, 'Error: Radar frequency not found in Mie tables, use ' // &
              'on-the-fly calculations'
    stop
   endif

  end subroutine load_mie_table
