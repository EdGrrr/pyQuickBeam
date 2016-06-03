  subroutine load_hydrometeor_classes(hclass_file,hp,nhclass)
  use radar_simulator_types
  implicit none
  
! Purpose:
!   Loads the hydrometeor classes to be used in calculations
!   Part of QuickBeam v1.1a by John Haynes
!   http://cloudsat.atmos.colostate.edu/radarsim
!
! Inputs:  
!   [hclass_file]   hydrometeor class input file
!
! Outputs:
!   [nhclass]       number of hydrometeor types
!   [hp]            structure that define hydrometeor types
!
! Modified:
!   08/23/2006  placed into subroutine form (Roger Marchand)
   
! ----- INPUT -----
  character*200, intent(in) :: hclass_file
  
! ----- OUTPUTS -----  
  integer, intent(out) :: nhclass
  type(class_param), intent(out) :: hp
  
! ----- INTERNAL -----  
  integer :: i
  real*8  :: val
  logical :: fexist 
  character :: comment*20
        
  inquire(file=hclass_file,exist=fexist) 
  if (fexist) then 
   
    open(9,file=hclass_file,action='read')
    read(9,*) comment
    nhclass = 0
    hp%rho(:) = -1
    hp%cp(:) = 0
    do i = 1,maxhclass
      read(9,*,end=500) val,hp%dtype(i),hp%col(i),hp%phase(i), &
        hp%cp(i),hp%dmin(i),hp%dmax(i),hp%apm(i),hp%bpm(i),hp%rho(i), &
        hp%p1(i),hp%p2(i),hp%p3(i)
      nhclass = nhclass + 1
    enddo
    500 close(9)
    
    hp%rho_eff(:,:) = -999.
    hp%ifc(:,:) = -9
    hp%idd(:) = -9
    
!   // setup scaling arrays
    hp%fc(:,:) = -999.
    hp%scaled(:) = .false.
         
  else
  
    print *, 'Error: hclass_file does not exist'
    stop
  
  endif
  
  end subroutine load_hydrometeor_classes
