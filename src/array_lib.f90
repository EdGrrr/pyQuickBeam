! ARRAY_LIB: Array procedures for F90
! Compiled/Modified:
!   07/01/06  John Haynes (haynes@atmos.colostate.edu)
!
! infind (function)
! lin_interpolate (function)
  
  module array_lib
  implicit none

  contains

! ----------------------------------------------------------------------------
! function INFIND
! ----------------------------------------------------------------------------
  function infind(list,val,sort,dist)
  use m_mrgrnk, only : mrgrnk
  implicit none
!
! Purpose:
!   Finds the index of an array that is closest to a value, plus the
!   difference between the value found and the value specified
!
! Inputs:
!   [list]   an array of sequential values
!   [val]    a value to locate
! Optional input:
!   [sort]   set to 1 if [list] is in unknown/non-sequential order
!
! Returns:
!   index of [list] that is closest to [val]
!
! Optional output:
!   [dist]   set to variable containing [list([result])] - [val]
!
! Requires:
!   mrgrnk library
!
! Created:
!   10/16/03  John Haynes (haynes@atmos.colostate.edu)
! Modified:
!   01/31/06  IDL to Fortran 90
 
! ----- INPUTS -----
  real*8, dimension(:), intent(in) :: list
  real*8, intent(in) :: val  
  integer, intent(in), optional :: sort
  
! ----- OUTPUTS -----
  integer*4 :: infind
  real*8, intent(out), optional :: dist

! ----- INTERNAL -----
  real*8, dimension(size(list)) :: lists
  integer*4 :: nlist, result, tmp(1), sort_list
  integer*4, dimension(size(list)) :: mask, idx

  if (present(sort)) then
    sort_list = sort
  else
    sort_list = 0
  endif  

  nlist = size(list)
  if (sort_list == 1) then
    call mrgrnk(list,idx)
    lists = list(idx)
  else
    lists = list
  endif

  if (val >= lists(nlist)) then
    result = nlist
  else if (val <= lists(1)) then
    result = 1
  else
    mask(:) = 0
    where (lists < val) mask = 1
      tmp = minloc(mask,1)
      if (abs(lists(tmp(1)-1)-val) < abs(lists(tmp(1))-val)) then
        result = tmp(1) - 1
      else
        result = tmp(1)
      endif
  endif
  if (present(dist)) dist = lists(result)-val
  if (sort_list == 1) then
    infind = idx(result)
  else
    infind = result
  endif

  end function infind

! ----------------------------------------------------------------------------
! function LINEAR_INTERPOLATE
! ----------------------------------------------------------------------------  
  function linear_interpolate(X,Y,Z,IND)
!
! Purpose:
!   Linearly interpolate a value
!
! Inputs:
!   [X]   tabulated values of the independent variable,
!         in ascending order
!   [Y]   tabulated values of the dependent variable
!   [Z]   point at which the function is to be estimated
!
! Returns:
!   Linearly interpolated value corresponding to [Z]
!
! Optional arguments:
!   [IND]   - use this optional argument if you have multiple [Y] vectors
!             to do interpolation over, for identical [X] and [Z]
!           - call linear_interpolate the first time with [IND] set to
!             a negative integer; [IND] will then be returned with
!             bounding index information that can speed up subsequent
!             interpolations using the same [X] and [Z]
!           - call linear_interpolate again with the updated value
!             of [IND]
!
! Notes:
!   - If [Z] is outside the range of [X], then the corresponding
!     lower or upper value of [Y] is returned
!
! Coded:
!   Based on binary search interpolation code by Denis O'Brien
!   08/23/07  to F90, John Haynes (haynes@atmos.colostate.edu)
!   09/13/07  IND keyword  (JMH)

! ----- INPUTS -----
  real*8, intent(in), dimension(:) ::  X,Y
  real*8, intent(in) ::  Z
  integer, intent(inout), dimension(2), optional :: IND

! ----- OUTPUTS -----
  real*8 :: linear_interpolate

! ----- INTERNAL -----
  real*8 :: XMin,XMax
  integer :: IMin,IMax,I,N
  logical :: store_IND
   
      
  store_IND = .false.

  N = size(X)
  if (X(1) < X(N)) then
    XMin=X(1)
    XMax=X(N)
    IMin=  1
    IMax=  N
  else
    XMin=X(N)
    XMax=X(1)
    IMin=  N
    IMax=  1
  endif

  if      (Z <= XMin) then
    linear_interpolate=Y(IMin)
  elseif (Z >= XMax) then
    Linear_Interpolate=Y(IMax)
  else

    if ( present(ind) ) then
      if (maxval(ind) > 0) then
        IMin = IND(1)
        IMax = IND(2)
        Linear_Interpolate=Y(IMin)+(Y(IMax)-Y(IMin))*&
                                   (Z      -X(IMin))/&
                                   (X(IMax)-X(IMin))
        return
      else
        store_ind = .true.
      endif
    endif

  ! // find indices encompassing the point Z
    do while (abs(IMax-IMin) > 1)
      I=(IMin+IMax)/2
      if (Z >= X(I)) then
        IMin=I
      else
        IMax=I
      end if
    enddo

  ! // estimate the function by linear interpolation
    linear_interpolate=Y(IMin)+(Y(IMax)-Y(IMin))*&
                               (Z      -X(IMin))/&
                               (X(IMax)-X(IMin))
                                     
    if (store_ind) then
      IND(1) = IMin
      IND(2) = IMax
    endif
                                    
  endif
  return

  end function linear_interpolate

  end module array_lib
