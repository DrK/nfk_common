!     Last change: UOSF 25/08/2011 20:43:00
Module KIND_PARAMETERS
!<cat>
!Norman Kirkby 13/2/2004
!Common Version created 7/1/06

! all the compiler dependent stuff including my_system
! This is a base class that almost evert other module and routine should use
!
!I removed the command line stuff, and used the new standard intrinsics
!in COMMAND_LINE_SWITCHES 29/4/11 Royal Wedding Day
!
!at v8 I introduced the error_message_queue
!which the user call pull from to determine if there have been errors
!
!</cat>
use ERROR_MESSAGE_QUEUE
implicit none

PUBLIC KIND_PARAMETERS_initialise,          &
       KIND_PARAMETERS_system,              &
       My_system,                           & !old version available while upgrading
       KIND_PARAMETERS_get_kind_name,       &
       KIND_PARAMETERS_test,                &
       KIND_PARAMETERS_test_get_reply_size, &
       KIND_PARAMETERS_get_test_reply_dim,  & !for the test system
       KIND_PARAMETERS_get_test_reply_len

INTEGER, parameter :: i_kind=4 !IVF
INTEGER, parameter :: r_kind=8 !IVF
INTEGER, parameter :: l_kind=1 !IVF

!gFINTEGER, parameter :: i_kind=4
!gFINTEGER, parameter :: r_kind=8
!gFINTEGER, parameter :: l_kind=1

!lf95INTEGER, parameter :: i_kind=4
!lf95INTEGER, parameter :: r_kind=8
!lf95INTEGER, parameter :: l_kind=1

!ftn95INTEGER, parameter :: i_kind=2
!ftn95INTEGER, parameter :: r_kind=2
!ftn95INTEGER, parameter :: l_kind=1

INTEGER, parameter :: kindr1=4  !IVF
INTEGER, parameter :: kindr2=8  !IVF
INTEGER, parameter :: kindr3=16 !IVF

!gFINTEGER, parameter :: kindr1=4  !gF
!gFINTEGER, parameter :: kindr2=8  !gF
!gFINTEGER, parameter :: kindr3=16 !gF

!lf95INTEGER, parameter :: kindr1=4
!lf95INTEGER, parameter :: kindr2=8
!lf95INTEGER, parameter :: kindr3=16

!ftn95INTEGER, parameter :: kindr1=1
!ftn95INTEGER, parameter :: kindr2=2
!ftn95INTEGER, parameter :: kindr3=3

INTEGER, parameter :: kindc1=4  !IVF
INTEGER, parameter :: kindc2=8  !IVF
INTEGER, parameter :: kindc3=16 !IVF

!gFINTEGER, parameter :: kindc1=4
!gFINTEGER, parameter :: kindc2=8
!gFINTEGER, parameter :: kindc3=16

!lf95INTEGER, parameter :: kindc1=4
!lf95INTEGER, parameter :: kindc2=8
!lf95INTEGER, parameter :: kindc3=16

!FTN95INTEGER, parameter :: kindc1=1
!FTN95INTEGER, parameter :: kindc2=2
!FTN95INTEGER, parameter :: kindc3=3

INTEGER, parameter :: kindi1=1 !IVF
INTEGER, parameter :: kindi2=2 !IVF
INTEGER, parameter :: kindi3=4 !IVF
INTEGER, parameter :: kindi4=8 !IVF

!gFINTEGER, parameter :: kindi1=1
!gFINTEGER, parameter :: kindi2=2
!gFINTEGER, parameter :: kindi3=4
!gFINTEGER, parameter :: kindi4=8

!lf95INTEGER, parameter :: kindi1=1
!lf95INTEGER, parameter :: kindi2=2
!lf95INTEGER, parameter :: kindi3=4
!lf95INTEGER, parameter :: kindi4=8

!ftn95INTEGER, parameter :: kindi1=1
!ftn95INTEGER, parameter :: kindi2=2
!ftn95INTEGER, parameter :: kindi3=3
!ftn95INTEGER, parameter :: kindi4=4

INTEGER, parameter :: kindl1=1 !IVF
INTEGER, parameter :: kindl2=2 !IVF
INTEGER, parameter :: kindl3=4 !IVF
INTEGER, parameter :: kindl4=8 !IVF

!gFINTEGER, parameter :: kindl1=1
!gFINTEGER, parameter :: kindl2=2
!gFINTEGER, parameter :: kindl3=4
!gFINTEGER, parameter :: kindl4=8

!lf95INTEGER, parameter :: kindl1=1
!lf95INTEGER, parameter :: kindl2=2
!lf95INTEGER, parameter :: kindl3=4
!lf95INTEGER, parameter :: kindl4=8

!only three logical kinds for silverfrost - a pain in the arse later

!ftn95INTEGER, parameter :: kindl1=1
!ftn95INTEGER, parameter :: kindl2=2
!ftn95INTEGER, parameter :: kindl3=3

CHARACTER (LEN=5), parameter :: kind_name_IVF   = 'IVF'    ! Intel Visual Fortran
CHARACTER (LEN=5), parameter :: kind_name_LF95  = 'LF95'   ! Lahey
CHARACTER (LEN=5), parameter :: kind_name_FTN95 = 'FTN95'  ! Silverfrost
CHARACTER (LEN=5), parameter :: kind_name_g95   = 'g95'    ! g95
CHARACTER (LEN=5), parameter :: kind_name_gF    = 'gF'     ! gFortran
CHARACTER (LEN=5), parameter :: kind_name_ABS   = 'ABS'    ! Absoft

!exactly one of the next statements must be active
CHARACTER (len=5)            :: kind_name = kind_name_IVF   !IVF
!LF95CHARACTER (len=5)           :: kind_name = kind_name_LF95
!FTN95CHARACTER (len=5)          :: kind_name = kind_name_FTN95
!G95CHARACTER (len=5)            :: kind_name = kind_name_G95
!GFCHARACTER (len=5)             :: kind_name = kind_name_GF
!ABSCHARACTER (len=5)            :: kind_name = kind_name_ABS

LOGICAL, PRIVATE             :: initialised = .false.
INTEGER, private                            :: my_iostat
character (len=256), private                :: my_errmsg, my_iomsg

INTEGER, private, parameter  :: reply_dim=27    !for reply in test
INTEGER, private, parameter  :: reply_len=100

!****************************************************************************************
contains
!****************************************************************************************
pure function KIND_PARAMETERS_get_test_reply_dim() RESULT(ii)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'KIND_PARAMETERS_get_test_reply_dim'
INTEGER :: ii
ii = reply_dim
return
END function KIND_PARAMETERS_get_test_reply_dim
!****************************************************************************************
pure function KIND_PARAMETERS_get_test_reply_len() RESULT(ii)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'KIND_PARAMETERS_get_test_reply_len'
INTEGER :: ii
ii = reply_len
return
END function KIND_PARAMETERS_get_test_reply_len
!****************************************************************************************
subroutine My_System(command) !backwards compatible
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'My_System'

CHARACTER (LEN=*), INTENT(IN) :: command
! to handle non-standard operating system calls

call KIND_PARAMETERS_system(command)
!lf has no way to determine the return code

return
end subroutine My_System
!****************************************************************************************
subroutine KIND_PARAMETERS_system(command)
!<cat>
!
!</cat>
use IFPORT                              !IVF
implicit none
character (:), parameter:: place = 'KIND_PARAMETERS_system'
CHARACTER (LEN=*), INTENT(IN) :: command
! to handle non-standard operating system calls
Logical                       :: result !IVF

!LF95call system(trim(command))   
!FTN95call system(trim(command))
!lf has no way to determine the return code
result = systemqq(command)              !IVF
!result is true if successful
!The standard intrinsic will go in here as it is adopted by the major compilers
return
end subroutine KIND_PARAMETERS_System
!****************************************************************************************
subroutine KIND_PARAMETERS_initialise(message)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'KIND_PARAMETERS_initialise'
CHARACTER (len=*), OPTIONAL, INTENT(OUT) :: message

initialised = .true.
if(present(message)) message = 'KIND_PARAMETERS_initialise initialised for '//kind_name !IVF
!LF95if(present(message))  message = 'KIND_PARAMETERS_initialise initialised for '//kind_name
!GFif(present(message))    message = 'KIND_PARAMETERS_initialise initialised for '//kind_name
!FTN95if(present(message)) message = 'KIND_PARAMETERS_initialise initialised for Silverfrost'

return
end subroutine kind_parameters_initialise
!****************************************************************************************
subroutine KIND_PARAMETERS_get_kind_name(ch)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'KIND_PARAMETERS_get_kind_name'
CHARACTER (len=*), intent(out) :: ch
ch = kind_name
return
end subroutine KIND_PARAMETERS_get_kind_name
!****************************************************************************************
subroutine KIND_PARAMETERS_test(reply)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'KIND_PARAMETERS_test'
CHARACTER (LEN=reply_len), DIMENSION(reply_dim), INTENT(OUT) :: reply !to transfer answers to the calling routine

WRITE(reply(1),fmt='(A)'   , iostat=my_iostat, iomsg=my_iomsg) 'Kind Parameters supplies the following constants'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(2),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'r_kind = ', r_kind
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(3),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'i_kind = ', i_kind
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(4),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'l_kind = ', l_kind
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(5),fmt='(A)'   , iostat=my_iostat, iomsg=my_iomsg)    'kind_name = '//TRIM(kind_name)&
               //' with the following auxilliary types'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(6),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg)  'kindr1 = ', kindr1
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(7),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg)  'kindr2 = ', kindr2
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(8),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg)  'kindr3 = ', kindr3
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(9),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg)  'kindc1 = ', kindc1
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(10),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'kindc2 = ', kindc2
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(11),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'kindc3 = ', kindc3
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(12),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'kindi1 = ', kindi1
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(13),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'kindi2 = ', kindi2
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(14),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'kindi3 = ', kindi3
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(15),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'kindi4 = ', kindi4
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(16),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'kindl1 = ', kindl1
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(17),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'kindl2 = ', kindl2
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(18),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg) 'Kind_Parameters has two public functions                     '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(19),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg) 'KIND_PARAMETERS_get_test_reply_dim                           '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(20),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg) 'KIND_PARAMETERS_get_test_reply_len                           '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(21),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg) 'Kind_Parameters has the following public subroutines         '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(22),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg) 'KIND_PARAMETERS_initialise(message)                          '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(23),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg) 'KIND_PARAMETERS_system(command)                              '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(24),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg) 'My_system(command)  -- legacy code use KIND_PARAMETERS_system'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(25),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg) 'KIND_PARAMETERS_get_kind_name(name)                          '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(26),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg) 'KIND_PARAMETERS_test(reply)                                  '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(27),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg) 'KIND_PARAMETERS_test_get_reply_size(dim,len)                 '
if(my_iostat /= 0) call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write')

return
end subroutine KIND_PARAMETERS_test
!******************************************************************************
subroutine KIND_PARAMETERS_test_get_reply_size(mydim, mylen)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'KIND_PARAMETERS_test_get_reply_size'
INTEGER, INTENT(OUT) :: mydim
INTEGER, INTENT(OUT) :: mylen

mydim = reply_dim
mylen = reply_len

return
end subroutine KIND_PARAMETERS_test_get_reply_size
!****************************************************************************************
end module KIND_PARAMETERS
