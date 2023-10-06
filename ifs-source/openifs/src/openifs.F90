!
!    OpenIFS specific code
!
!	Purpose
!	-------
!	OpenIFS is a subset of the main IFS code. This file holds
!	all the OpenIFS variables and code. In order to keep
!	the interface between IFS and OpenIFS as clean as possible,
!	all OpenIFS specific code should go in here. No OpenIFS 
!	variables should be used directly by the model; access should
!	>always< be provided by functions/subroutines.
!
!	Currently this is not a module to avoid inserting 'use openifs'
!	in all the code needed, though it would be cleaner to make it one.

!	Interface
!	---------
!	The interface between the model and OpenIFS should be only via
!	the functions defined in here. 
!	All variables/functions/subroutines should be prefixed with 'oifs_'.

!	Author
!	------
!		Glenn Carver. *ECMWF* 2012-2021

!	Modifications
!	-------------
!		Original version: feb 2012. GC

!	-------------------------------------------------------------------

subroutine oifs_dummy( sname, fname )

  implicit none

  character(len=*), intent(in) :: sname ! subroutine call this dummy replaced
  character(len=*), intent(in) :: fname ! file in which the call was replaced

  ! OIFS_DUMMY_ACTION environment variable
  ! If set to 'quiet'                 then do/print nothing.
  !           'verbose'               print when subroutine called
  !           'abort'                 print message then call abor1
  ! If unset, default is 'abort'
  character(len=20), save :: action_value = ''
  logical, save           :: value_set = .false.

  if ( .not. value_set ) then
     call get_environment_variable( 'OIFS_DUMMY_ACTION', action_value, trim_name=.true. )
     if ( len_trim(action_value) == 0 ) action_value = 'abort'
     value_set = .true.
  endif

  if ( trim(action_value) == 'quiet' )  return

  if ( trim(action_value) == 'verbose' .or. trim(action_value) == 'abort' ) then
     write(0,*) 'oifs_dummy called as : ',sname, ' from file : ', fname
     call flush(0)
  endif
  if ( trim(action_value) == 'abort' ) then
     ! kills the run
     call abor1('Error! OpenIFS dummy routine '//trim(sname)//' was called. Contact support.')
  endif
  return
end subroutine 


subroutine oifs_dummympi( sname, fname )

  implicit none

  character(len=*), intent(in) :: sname ! subroutine call this dummy replaced
  character(len=*), intent(in) :: fname ! file in which the call was replaced

  ! OIFS_DUMMYMPI_ACTION environment variable
  ! If set to 'quiet'                 then do/print nothing.
  !           'verbose'               print when subroutine called
  !           'abort'                 print message then call abor1
  ! If unset, default is 'quiet'
  character(len=20), save :: action_value = ''
  logical, save           :: value_set = .false.

  if ( .not. value_set ) then
     call get_environment_variable( 'OIFS_DUMMYMPI_ACTION', action_value, trim_name=.true. )
     if ( len_trim(action_value) == 0 ) action_value = 'quiet'
     value_set = .true.
  endif

  if ( trim(action_value) == 'quiet' )  return

  if ( trim(action_value) == 'verbose' .or. trim(action_value) == 'abort' ) then
     write(0,*) 'oifs_dummympi called as : ',sname, ' from file : ', fname
     call flush(0)
  endif
  if ( trim(action_value) == 'abort' ) then
     ! kills the run
     call abor1('Error! OpenIFS dummy mpi subroutine was called. Contact support.')
  endif
  return
end subroutine 



subroutine oifs_print_copyright(kout)
  implicit none
  integer, intent(in) :: kout
  character(len=*), parameter :: oifs_cycle = 'CY43R3'
  character(len=*), parameter :: oifs_version = '2'
  
  write(kout,'(//70(''='')/5(10x,a/)70(''='')//)' ) &
        'ECMWF OpenIFS model', &
        'OpenIFS is licensed software (c) 2011-2021', &
        'Version: '//oifs_get_vsn(), &
        'Website: '//oifs_get_weburl(), &
        'Support: '//oifs_get_supportemail()

  return

contains

character(len=64) function oifs_get_vsn()
  implicit none
  oifs_get_vsn = 'OIFS'//oifs_cycle//'v'//oifs_version
end function

character(len=64) function oifs_get_supportemail()
  implicit none
  oifs_get_supportemail = 'openifs-support@ecmwf.int'
end function

character(len=64) function oifs_get_weburl()
  implicit none
  oifs_get_weburl = 'https://software.ecmwf.int/oifs/'
end function

end subroutine oifs_print_copyright
