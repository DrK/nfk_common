!     Last change: UOSF 04/09/2011 10:50:18
module SPECIAL_CONSTANTS

USE KIND_PARAMETERS
USE VERBOSITY

implicit none

public  SPECIAL_CONSTANTS_initialise_verbose,  &
        SPECIAL_CONSTANTS_get_verbose,         &
        SPECIAL_CONSTANTS_initialise,          &
        IsSmall, IsExpBig, IsBetween,          &
        SPECIAL_CONSTANTS_test,                &
        SPECIAL_CONSTANTS_test_get_reply_size, &
        SPECIAL_CONSTANTS_colour_to_number8,   &
        SPECIAL_CONSTANTS_colour_to_number24,  &
        SPECIAL_CONSTANTS_next_colour24,       &
        SPECIAL_CONSTANTS_html_constant_found

REAL (KIND=R_KIND), PARAMETER :: ZERO     = 0.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: ONE      = 1.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: TWO      = 2.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: THREE    = 3.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: FOUR     = 4.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: FIVE     = 5.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: SIX      = 6.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: SEVEN    = 7.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: EIGHT    = 8.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: NINE     = 9.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: TEN      = 10.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: ELEVEN   = 11.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: TWELVE   = 12.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: THIRTEEN = 13.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: TWENTY  = 20.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: THIRTY  = 30.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: FORTY   = 40.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: FIFTY   = 50.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: SIXTY   = 60.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: SEVENTY = 70.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: EIGHTY  = 80.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: NINETY  = 90.0_r_kind

REAL (KIND=R_KIND), PARAMETER :: HALF     = 0.5_r_kind

REAL (KIND=R_KIND), PARAMETER :: HUNDRED  = 100.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: THOUSAND = 1000.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: MILLION  = 1000000.0_r_kind
REAL (KIND=R_KIND), PARAMETER :: INFINITY = HUGE(zero)

REAL (KIND=R_KIND), DIMENSION(2), PARAMETER    :: POSITIVE = (/ZERO,INFINITY/)
REAL (KIND=R_KIND), DIMENSION(2), PARAMETER    :: NEGATIVE = (/-INFINITY, ZERO/)

INTEGER (KIND=I_KIND), DIMENSION(2), parameter :: POSITIVE_INTEGER_I_KIND = (/INT(0,i_kind), HUGE(INT(0,i_kind))/)
INTEGER (KIND=KINDI1), DIMENSION(2), parameter :: POSITIVE_INTEGER_KINDI1 = (/INT(0,kindi1), HUGE(INT(0,kindi1))/)
INTEGER (KIND=KINDI2), DIMENSION(2), parameter :: POSITIVE_INTEGER_KINDI2 = (/INT(0,kindi2), HUGE(INT(0,kindi2))/)
INTEGER (KIND=KINDI3), DIMENSION(2), parameter :: POSITIVE_INTEGER_KINDI3 = (/INT(0,kindi3), HUGE(INT(0,kindi3))/)
!INTEGER (KIND=KINDI4), DIMENSION(2), parameter :: POSITIVE_INTEGER_KINDI4 = (/INT(0,kindi4), HUGE(INT(0,kindi4))/)!LF95
INTEGER              , DIMENSION(2), parameter :: POSITIVE_INTEGER_DEFAULT = (/0, HUGE(0)/)

INTEGER (KIND=I_KIND), DIMENSION(2), parameter :: NEGATIVE_INTEGER_I_KIND = (/-HUGE(INT(0,i_kind)), INT(0,i_kind)/)
INTEGER (KIND=KINDI1), DIMENSION(2), parameter :: NEGATIVE_INTEGER_KINDI1 = (/-HUGE(INT(0,kindi1)), INT(0,kindi1)/)
INTEGER (KIND=KINDI2), DIMENSION(2), parameter :: NEGATIVE_INTEGER_KINDI2 = (/-HUGE(INT(0,kindi2)), INT(0,kindi2)/)
INTEGER (KIND=KINDI3), DIMENSION(2), parameter :: NEGATIVE_INTEGER_KINDI3 = (/-HUGE(INT(0,kindi3)), INT(0,kindi3)/)
INTEGER (KIND=KINDI4), DIMENSION(2), parameter :: NEGATIVE_INTEGER_KINDI4 = (/-HUGE(INT(0,kindi4)), INT(0,kindi4)/)
INTEGER              , DIMENSION(2), parameter :: NEGATIVE_INTEGER_DEFAULT = (/-HUGE(0), 0/)

REAL (KIND=R_KIND), protected                 :: PI       !accessible outside but cannot change
REAL (KIND=R_KIND), protected                 :: two_pi
REAL (KIND=R_KIND), protected                 :: log2
REAL (KIND=R_KIND), protected                 :: log_infinity
REAL (KIND=R_KIND), protected                 :: one_by_sqrt_2Pi
REAL (KIND=r_kind), protected                 :: sqrt2
REAL (KIND=r_kind), protected                 :: pi_by_2
REAL (KIND=r_kind), protected                 :: root_pi_by_2

!The next two may be a mistake
LOGICAL, PARAMETER                            :: T = .TRUE., F = .false.
CHARACTER (LEN=1), PARAMETER                  :: Y = 'Y',    N = 'N'

CHARACTER (LEN=1), PARAMETER                  :: equals  = '='
CHARACTER (LEN=2), PARAMETER                  :: CRLF    = achar(13)//achar(11)  !but they do not
CHARACTER (LEN=3), PARAMETER                  :: degC    = ' '//char(167)//'C'
CHARACTER (len=1), PARAMETER                  :: squared = char(253)
                                 !work in conventional output
CHARACTER (LEN=3), PARAMETER                  :: OFF     = 'off'
CHARACTER (LEN=2), PARAMETER                  :: ON      = 'on'
CHARACTER (LEN=3), DIMENSION(2), PARAMETER    :: ON_OFF  = (/on//' ',off/)
CHARACTER (LEN=3), DIMENSION(2), PARAMETER    :: YES_NO  = (/'yes', 'no '/)
CHARACTER (LEN=1), DIMENSION(2), PARAMETER    :: Y_N     = (/Y,N/)
CHARACTER (LEN=1), PARAMETER                  :: path_separator = '\' ! for Window and DOS but not UNIX??

!colours for html output

CHARACTER (len=7), parameter  :: red         = '#ff0000'
CHARACTER (len=7), parameter  :: black       = '#000000'
CHARACTER (len=7), parameter  :: white       = '#ffffff'
CHARACTER (len=7), parameter  :: yellow      = '#ffff00'
CHARACTER (len=7), parameter  :: blue        = '#3333ff'
CHARACTER (len=7), parameter  :: grey        = '#999999'
CHARACTER (len=7), parameter  :: green       = '#33ff33'
CHARACTER (len=7), parameter  :: light_blue  = '#33ffff'
CHARACTER (len=7), parameter  :: dark_blue   = '#000066'
CHARACTER (len=7), parameter  :: dark_green  = '#009900'
CHARACTER (len=7), parameter  :: light_green = '#99ff99'
CHARACTER (len=7), parameter  :: pink        = '#ff99ff'
CHARACTER (len=7), parameter  :: violet      = '#993399'
CHARACTER (len=7), parameter  :: orange      = '#ffcc00'
CHARACTER (len=7), parameter  :: dark_brown  = '#993300'
CHARACTER (len=7), parameter  :: light_brown = '#cc9933'
CHARACTER (len=7), parameter  :: brown       = '#cc6600'
CHARACTER (len=7), parameter  :: turquoise   = '#66cccc'
CHARACTER (len=7), parameter  :: lemon       = '#ffffcc'
CHARACTER (len=7), parameter  :: salmon      = '#ffcccc'
CHARACTER (len=7), parameter  :: charcoal    = '#333333'
CHARACTER (len=7), parameter  :: purple      = '#cc33cc'
CHARACTER (len=7), parameter  :: maroon      = '#660033'
CHARACTER (len=7), parameter  :: mauve       = '#663399'
CHARACTER (len=7), parameter  :: olive       = '#666600'

!special characters for html output

CHARACTER (len=6), parameter  :: html_copy   = '&#169;' !copyright = char(169)(in Arial)
CHARACTER (len=6), parameter  :: html_reg    = '&#174;' !registered = char(174)( " )
CHARACTER (len=6), parameter  :: html_deg    = '&#176;' !as in deg C = char(248)
CHARACTER (len=6), parameter  :: html_pmi    = '&#177;' !plus minus = char(242)
CHARACTER (len=6), parameter  :: html_sqd    = '&#178;' !squared = char(178)
CHARACTER (len=6), parameter  :: html_cbd    = '&#179;' !cubed = char(179)
CHARACTER (len=6), parameter  :: html_sup1   = '&#185;' !superscript 1 char(185)
CHARACTER (len=6), parameter  :: html_qtr    = '&#188;' !quarter = char(172)
CHARACTER (len=6), parameter  :: html_hlf    = '&#189;' !half = char(171)
CHARACTER (len=6), parameter  :: html_3qtr   = '&#190;' !3/4 = char(243)
CHARACTER (len=7), parameter  :: html_pd     = '&#8706;' !partial d
CHARACTER (len=7), parameter  :: html_sqrt   = '&#8730;' !square root
CHARACTER (len=7), parameter  :: html_inf    = '&#8734;' !infinity
CHARACTER (len=7), parameter  :: html_int    = '&#8747;' !integral sign
CHARACTER (len=7), parameter  :: html_aeq    = '&#8773;' !approx =
CHARACTER (len=7), parameter  :: html_neq    = '&#8800;' !not =
CHARACTER (len=7), parameter  :: html_leq    = '&#8804;' !less than or =
CHARACTER (len=7), parameter  :: html_geq    = '&#8805;' !gt than or =
CHARACTER (len=6), parameter  :: html_sum    = '&#931;'  !summation
CHARACTER (len=7), parameter  :: html_sep    = '&#8260;' !fraction separator

!other useful characters in dos
! much greater than = char(175)
! much less than    = char(174)
! superscript minus = char(238)

!Greek alphabet for html

CHARACTER (len=7), parameter  :: html_alpha   = '&alpha;'
CHARACTER (len=6), parameter  :: html_beta    = '&beta;'
CHARACTER (len=7), parameter  :: html_gamma   = '&gamma;'
CHARACTER (len=7), parameter  :: html_delta   = '&delta;'
CHARACTER (len=9), parameter  :: html_epsilon = '&epsilon;'
CHARACTER (len=6), parameter  :: html_zeta    = '&zeta;'
CHARACTER (len=5), parameter  :: html_eta     = '&eta;'
CHARACTER (len=7), parameter  :: html_theta   = '&theta;'
CHARACTER (len=6), parameter  :: html_iota    = '&iota;'
CHARACTER (len=7), parameter  :: html_kappa   = '&kappa;'
CHARACTER (len=8), parameter  :: html_lambda  = '&lambda;'
CHARACTER (len=4), parameter  :: html_mu      = '&mu;'
CHARACTER (len=4), parameter  :: html_nu      = '&nu;'
CHARACTER (len=4), parameter  :: html_xi      = '&xi;'
CHARACTER (len=9), parameter  :: html_omicron = '&omicron;'
CHARACTER (len=4), parameter  :: html_pi      = '&pi;'
CHARACTER (len=5), parameter  :: html_rho     = '&rho;'
CHARACTER (len=7), parameter  :: html_sigma   = '&sigma;'
CHARACTER (len=5), parameter  :: html_tau     = '&tau;'
CHARACTER (len=9), parameter  :: html_upsilon = '&upsilon;'
CHARACTER (len=5), parameter  :: html_phi     = '&phi;'
CHARACTER (len=5), parameter  :: html_chi     = '&chi;'
CHARACTER (len=5), parameter  :: html_psi     = '&psi;'
CHARACTER (len=7), parameter  :: html_omega   = '&omega;'

CHARACTER (len=7), parameter  :: html_u_alpha   = '&Alpha;'
CHARACTER (len=6), parameter  :: html_u_beta    = '&Beta;'
CHARACTER (len=7), parameter  :: html_u_gamma   = '&Gamma;'
CHARACTER (len=7), parameter  :: html_u_delta   = '&Delta;'
CHARACTER (len=9), parameter  :: html_u_epsilon = '&Epsilon;'
CHARACTER (len=6), parameter  :: html_u_zeta    = '&Zeta;'
CHARACTER (len=5), parameter  :: html_u_eta     = '&Eta;'
CHARACTER (len=7), parameter  :: html_u_theta   = '&Theta;'
CHARACTER (len=6), parameter  :: html_u_iota    = '&Iota;'
CHARACTER (len=7), parameter  :: html_u_kappa   = '&Kappa;'
CHARACTER (len=8), parameter  :: html_u_lambda  = '&Lambda;'
CHARACTER (len=4), parameter  :: html_u_mu      = '&Mu;'
CHARACTER (len=4), parameter  :: html_u_nu      = '&Nu;'
CHARACTER (len=4), parameter  :: html_u_xi      = '&Xi;'
CHARACTER (len=9), parameter  :: html_u_omicron = '&Omicron;'
CHARACTER (len=4), parameter  :: html_u_pi      = '&Pi;'
CHARACTER (len=5), parameter  :: html_u_rho     = '&Rho;'
CHARACTER (len=7), parameter  :: html_u_sigma   = '&Sigma;'
CHARACTER (len=5), parameter  :: html_u_tau     = '&Tau;'
CHARACTER (len=9), parameter  :: html_u_upsilon = '&Upsilon;'
CHARACTER (len=5), parameter  :: html_u_phi     = '&Phi;'
CHARACTER (len=5), parameter  :: html_u_chi     = '&Chi;'
CHARACTER (len=5), parameter  :: html_u_psi     = '&Psi;'
CHARACTER (len=7), parameter  :: html_u_omega   = '&Omega;'

! The official colours for the web

CHARACTER (len=7), parameter  :: web_aqua    = 'aqua'
CHARACTER (len=7), parameter  :: web_black   = 'black'
CHARACTER (len=7), parameter  :: web_white   = 'white'
CHARACTER (len=7), parameter  :: web_silver  = 'silver'
CHARACTER (len=7), parameter  :: web_yellow  = 'yellow'
CHARACTER (len=7), parameter  :: web_teal    = 'teal'
CHARACTER (len=7), parameter  :: web_fushia  = 'fushia'
CHARACTER (len=7), parameter  :: web_red     = 'red'
CHARACTER (len=7), parameter  :: web_blue    = 'blue'
CHARACTER (len=7), parameter  :: web_green   = 'green'
CHARACTER (len=7), parameter  :: web_gray    = 'gray'
CHARACTER (len=7), parameter  :: web_lime    = 'lime'
CHARACTER (len=7), parameter  :: web_maroon  = 'maroon'
CHARACTER (len=7), parameter  :: web_navy    = 'navy'
CHARACTER (len=7), parameter  :: web_olive   = 'olive'
CHARACTER (len=7), parameter  :: web_purple  = 'purple'

! Some useful fonts

CHARACTER (len=20), parameter :: Arial              = '"Arial"'
CHARACTER (len=20), parameter :: Courier_New        = '"Courier New"'
CHARACTER (len=20), parameter :: Times_New_Roman    = '"Times New Roman"'
CHARACTER (len=20), parameter :: Symbol             = '"Symbol"'
CHARACTER (len=20), parameter :: Helvetica_95_Black = '"Helvetica 95 Black"'
CHARACTER (len=20), parameter :: Helvetica_45_Light = '"Helvetica 45 Light"'

!
! Finally the local privates
!

logical, private, dimension(verbose_levels) :: verbose = verbose_defaults ! The local verbose

interface IsSmall
   module PROCEDURE IsSmall4, IsSmall8, IsSmall16
end interface

interface IsExpBig  ! are we going to overflow if we exponentiate
   module PROCEDURE IsExpBig4, IsExpBig8, IsExpBig16
end interface

interface IsBetween
   module PROCEDURE IsBetweenr4, IsBetweenr8, IsBetweenr16
   module PROCEDURE IsBetweeni1, IsBetweeni2, IsBetweeni4, IsBetweeni8
end interface

INTEGER, private, parameter :: reply_dim=179   !for reply in test
INTEGER, private, parameter :: reply_len=100

INTEGER, private                :: my_iostat
character (len=256), private    :: my_errmsg, my_iomsg

!******************************************************************************
contains
!******************************************************************************
function check_html(ch, html, found, pos, letter_in, letter_out, lench) result(gotcha)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'check_html'
logical                      :: gotcha
character(len=*), intent(in) :: ch, html
integer, intent(out)         :: pos, lench
character(len=1), intent(in) :: letter_in
character(len=1), intent(out):: Letter_out
logical, intent(inout)       :: found

gotcha = .false.
found  = .false.
pos    = index(ch,html)
if(pos /= 0) then
   gotcha = .true.
   found  = .true.
   letter_out = Letter_in
   lench      = len_trim(html)
endif

return
end function check_html
!******************************************************************************
function SPECIAL_CONSTANTS_html_constant_found(ch, myposition, mylen, myletter) &
                                               result(found)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'SPECIAL_CONSTANTS_html_constant_found'
character(len=*), intent(in) :: ch
character(len=1), optional   :: myletter
integer, optional            :: myposition, mylen
integer                      :: pos, lench
character(len=1)             :: letter
logical                      :: found

found  = .False.
letter = ' '

if(check_html(ch, html_alpha,    found, pos, 'a', letter, lench)) goto 100
if(check_html(ch, html_u_alpha,  found, pos, 'A', letter, lench)) goto 100
if(check_html(ch, html_beta,     found, pos, 'b', letter, lench)) goto 100
if(check_html(ch, html_u_beta,   found, pos, 'B', letter, lench)) goto 100
if(check_html(ch, html_gamma,    found, pos, 'g', letter, lench)) goto 100
if(check_html(ch, html_u_gamma,  found, pos, 'G', letter, lench)) goto 100
if(check_html(ch, html_delta,    found, pos, 'd', letter, lench)) goto 100
if(check_html(ch, html_u_delta,  found, pos, 'D', letter, lench)) goto 100
if(check_html(ch, html_epsilon,  found, pos, 'e', letter, lench)) goto 100
if(check_html(ch, html_u_epsilon,found, pos, 'E', letter, lench)) goto 100
if(check_html(ch, html_zeta,     found, pos, 'z', letter, lench)) goto 100
if(check_html(ch, html_u_zeta,   found, pos, 'Z', letter, lench)) goto 100
if(check_html(ch, html_eta,      found, pos, 'h', letter, lench)) goto 100
if(check_html(ch, html_u_eta,    found, pos, 'H', letter, lench)) goto 100
if(check_html(ch, html_theta,    found, pos, 'q', letter, lench)) goto 100
if(check_html(ch, html_u_theta,  found, pos, 'Q', letter, lench)) goto 100
if(check_html(ch, html_iota,     found, pos, 'i', letter, lench)) goto 100
if(check_html(ch, html_u_iota,   found, pos, 'I', letter, lench)) goto 100
if(check_html(ch, html_kappa,    found, pos, 'k', letter, lench)) goto 100
if(check_html(ch, html_u_kappa,  found, pos, 'K', letter, lench)) goto 100
if(check_html(ch, html_lambda,   found, pos, 'l', letter, lench)) goto 100
if(check_html(ch, html_u_lambda, found, pos, 'L', letter, lench)) goto 100
if(check_html(ch, html_mu,       found, pos, 'm', letter, lench)) goto 100
if(check_html(ch, html_u_mu,     found, pos, 'M', letter, lench)) goto 100
if(check_html(ch, html_nu,       found, pos, 'n', letter, lench)) goto 100
if(check_html(ch, html_u_nu,     found, pos, 'N', letter, lench)) goto 100
if(check_html(ch, html_xi,       found, pos, 'x', letter, lench)) goto 100
if(check_html(ch, html_u_xi,     found, pos, 'X', letter, lench)) goto 100
if(check_html(ch, html_omicron,  found, pos, 'o', letter, lench)) goto 100
if(check_html(ch, html_u_omicron,found, pos, 'O', letter, lench)) goto 100
if(check_html(ch, html_pi,       found, pos, 'p', letter, lench)) goto 100
if(check_html(ch, html_u_pi,     found, pos, 'P', letter, lench)) goto 100
if(check_html(ch, html_rho,      found, pos, 'r', letter, lench)) goto 100
if(check_html(ch, html_u_rho,    found, pos, 'R', letter, lench)) goto 100
if(check_html(ch, html_sigma,    found, pos, 's', letter, lench)) goto 100
if(check_html(ch, html_u_sigma,  found, pos, 'S', letter, lench)) goto 100
if(check_html(ch, html_tau,      found, pos, 't', letter, lench)) goto 100
if(check_html(ch, html_u_tau,    found, pos, 'T', letter, lench)) goto 100
if(check_html(ch, html_upsilon,  found, pos, 'u', letter, lench)) goto 100
if(check_html(ch, html_u_upsilon,found, pos, 'U', letter, lench)) goto 100
if(check_html(ch, html_phi,      found, pos, 'f', letter, lench)) goto 100
if(check_html(ch, html_u_phi,    found, pos, 'F', letter, lench)) goto 100
if(check_html(ch, html_chi,      found, pos, 'c', letter, lench)) goto 100
if(check_html(ch, html_u_chi,    found, pos, 'C', letter, lench)) goto 100
if(check_html(ch, html_psi,      found, pos, 'y', letter, lench)) goto 100
if(check_html(ch, html_u_psi,    found, pos, 'Y', letter, lench)) goto 100
if(check_html(ch, html_omega,    found, pos, 'w', letter, lench)) goto 100
if(check_html(ch, html_u_omega,  found, pos, 'W', letter, lench)) goto 100

!and some special characters

if(check_html(ch, html_sqd,  found, pos, char(178), letter, lench)) goto 100 !use FFHelvetica
if(check_html(ch, html_cbd,  found, pos, char(179), letter, lench)) goto 100 !use FFHelvetica
if(check_html(ch, html_copy, found, pos, char(169), letter, lench)) goto 100
if(check_html(ch, html_reg,  found, pos, char(174), letter, lench)) goto 100
if(check_html(ch, html_deg,  found, pos, char(176), letter, lench)) goto 100
if(check_html(ch, html_pmi,  found, pos, char(177), letter, lench)) goto 100
if(check_html(ch, html_sup1, found, pos, char(185), letter, lench)) goto 100
if(check_html(ch, html_qtr,  found, pos, char(188), letter, lench)) goto 100
if(check_html(ch, html_hlf,  found, pos, char(189), letter, lench)) goto 100
if(check_html(ch, html_3qtr, found, pos, char(190), letter, lench)) goto 100
if(check_html(ch, html_pd,   found, pos, char(182), letter, lench)) goto 100
if(check_html(ch, html_sqrt, found, pos, char(214), letter, lench)) goto 100
if(check_html(ch, html_inf,  found, pos, char(165), letter, lench)) goto 100
if(check_html(ch, html_int,  found, pos, char(242), letter, lench)) goto 100
if(check_html(ch, html_aeq,  found, pos, char(187), letter, lench)) goto 100
if(check_html(ch, html_neq,  found, pos, char(185), letter, lench)) goto 100
if(check_html(ch, html_leq,  found, pos, char(163), letter, lench)) goto 100
if(check_html(ch, html_geq,  found, pos, char(179), letter, lench)) goto 100 
if(check_html(ch, html_sum,  found, pos, char(229), letter, lench)) goto 100 
if(check_html(ch, html_sep,  found, pos, char(164), letter, lench)) goto 100

100 continue      !goto justified in this rare case?

if(present(myletter)  ) myletter   = letter
if(present(myposition)) myposition = pos
if(present(mylen)     ) mylen      = lench

return
end function SPECIAL_CONSTANTS_html_constant_found
!******************************************************************************
function SPECIAL_CONSTANTS_font_for_html_constant(ch) result(font)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'SPECIAL_CONSTANTS_font_for_html_constant'
character (len=*):: ch
integer          :: font
integer          :: pos, lench
character(len=1) :: letter
logical          :: found
integer          :: my_arial  = 204
integer          :: my_symbol = 207

!the default font is FFSymbol

!  FFSoftware  (1) Current software font
!  FFDriver    (2) Current TrueType font
!  FFUser      (100) User defined GDI font
!  FFCourier   (101) Courier
!  FFHelvetica (102) Helvetica/Arial
!  FFTimes     (103) Times Roman
!  FFSymbol    (104) Symbol

! Check with font manager, but...
! Use 204 for FFuser Arial
! Use 207 for FFUser Symbol


font = my_symbol !this will work for all the greek characters

if(    check_html(ch, html_sqd,  found, pos, char(178), letter, lench)) then
   font = my_arial
elseif(check_html(ch, html_cbd,  found, pos, char(179), letter, lench)) then
   font = my_arial
elseif(check_html(ch, html_copy, found, pos, char(169), letter, lench)) then
   font = my_arial
elseif(check_html(ch, html_reg,  found, pos, char(174), letter, lench)) then
   font = my_arial
elseif(check_html(ch, html_deg,  found, pos, char(176), letter, lench)) then
   font = my_arial
elseif(check_html(ch, html_pmi,  found, pos, char(177), letter, lench)) then
   font = my_arial
elseif(check_html(ch, html_sup1, found, pos, char(185), letter, lench)) then
   font = my_arial
elseif(check_html(ch, html_qtr,  found, pos, char(188), letter, lench)) then
   font = my_arial
elseif(check_html(ch, html_hlf,  found, pos, char(189), letter, lench)) then
   font = my_arial
elseif(check_html(ch, html_3qtr, found, pos, char(190), letter, lench)) then
   font = my_arial
elseif(check_html(ch, html_pd,   found, pos, char(182), letter, lench)) then
   font = my_symbol
elseif(check_html(ch, html_sqrt, found, pos, char(214), letter, lench)) then
   font = my_symbol
elseif(check_html(ch, html_inf,  found, pos, char(165), letter, lench)) then
   font = my_symbol
elseif(check_html(ch, html_int,  found, pos, char(242), letter, lench)) then
   font = my_symbol
elseif(check_html(ch, html_aeq,  found, pos, char(187), letter, lench)) then
   font = my_symbol
elseif(check_html(ch, html_neq,  found, pos, char(185), letter, lench)) then
   font = my_symbol
elseif(check_html(ch, html_leq,  found, pos, char(163), letter, lench)) then
   font = my_symbol
elseif(check_html(ch, html_geq,  found, pos, char(179), letter, lench)) then
   font = my_symbol
elseif(check_html(ch, html_sum,  found, pos, char(229), letter, lench)) then
   font = my_symbol
elseif(check_html(ch, html_sep,  found, pos, char(164), letter, lench)) then
   font = my_symbol
endif

return
end function SPECIAL_CONSTANTS_font_for_html_constant
!******************************************************************************
LOGICAL function IsBetweenr4(x,x1,x2)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsBetweenr4'
real(kind=kindr1) :: x,x1,x2
IsBetweenr4 = (x.ge.x1 .and. x.le.x2) .OR. (x.ge.x2 .AND. x.le.x1)
return
end function IsBetweenr4
!******************************************************************************
LOGICAL function IsBetweenr8(x,x1,x2)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsBetweenr8'
real(kind=kindr2) :: x,x1,x2
IsBetweenr8 = (x.ge.x1 .and. x.le.x2) .OR. (x.ge.x2 .AND. x.le.x1)
return
end function IsBetweenr8
!******************************************************************************
LOGICAL function IsBetweenr16(x,x1,x2)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsBetweenr16'
real(kind=kindr3) :: x,x1,x2
IsBetweenr16 = (x.ge.x1 .and. x.le.x2) .OR. (x.ge.x2 .AND. x.le.x1)
return
end function IsBetweenr16
!******************************************************************************
LOGICAL function IsBetweeni1(x,x1,x2)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsBetweeni1'
integer(kind=kindi1) :: x,x1,x2
IsBetweeni1 = (x.ge.x1 .and. x.le.x2) .OR. (x.ge.x2 .AND. x.le.x1)
return
end function IsBetweeni1
!******************************************************************************
LOGICAL function IsBetweeni2(x,x1,x2)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsBetweeni2'
integer(kind=kindi2) :: x,x1,x2
IsBetweeni2 = (x.ge.x1 .and. x.le.x2) .OR. (x.ge.x2 .AND. x.le.x1)
return
end function IsBetweeni2
!******************************************************************************
LOGICAL function IsBetweeni4(x,x1,x2)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsBetweeni4'
integer(kind=kindi3) :: x,x1,x2
IsBetweeni4 = (x.ge.x1 .and. x.le.x2) .OR. (x.ge.x2 .AND. x.le.x1)
return
end function IsBetweeni4
!******************************************************************************
LOGICAL function IsBetweeni8(x,x1,x2)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsBetweeni8'
integer(kind=kindi4) :: x,x1,x2
IsBetweeni8 = (x.ge.x1 .and. x.le.x2) .OR. (x.ge.x2 .AND. x.le.x1)
return
end function IsBetweeni8
!******************************************************************************
LOGICAL function IsExpBig4(x)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsExpBig4'
REAL (kind=kindr1) :: x,y

y = log(huge(x))

IsExpBig4 = ( x .ge. y )

return
end function IsExpBig4
!******************************************************************************
LOGICAL function IsExpBig8(x)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsExpBig8'
REAL (kind=kindr2) :: x,y

y = log(huge(x))

IsExpBig8 = ( x .ge. y )

return
end function IsExpBig8
!******************************************************************************
LOGICAL function IsExpBig16(x)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsExpBig16'
REAL (kind=kindr3) :: x,y

y = log(huge(x))

IsExpBig16 = ( x .ge. y )

return
end function IsExpBig16
!******************************************************************************
LOGICAL function IsSmall4(x)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsSmall4'
REAL (kind=kindr1) :: x,y

y = ten*tiny(x)

!IsSmall4 = (abs(x) .lt. y)

IsSmall4 = ( x .lt. y ) .AND. (x .GT. -y)

return
end function IsSmall4
!******************************************************************************
LOGICAL function IsSmall8(x)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsSmall8'
REAL (kind=kindr2) :: x,y

y = ten*tiny(x)

IsSmall8 = ( x .lt. y ) .AND. (x .GT. -y)

return
end function IsSmall8
!******************************************************************************
LOGICAL function IsSmall16(x)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'IsSmall16'
REAL (kind=kindr3) :: x,y

y = ten*tiny(x)

IsSmall16 = ( x .lt. y ) .AND. (x .GT. -y)

return
end function IsSmall16
!******************************************************************************
function SPECIAL_CONSTANTS_colour_to_number8(colour_name) RESULT (colour)
!<cat>
!
!</cat>
!Upcase must be called before entry
!These are named Winteracter colours
implicit none
character (:), parameter:: place = 'SPECIAL_CONSTANTS_colour_to_number8'
INTEGER           :: colour
CHARACTER (LEN=*) :: colour_name

select case (colour_name)
    case ('BLACK');   colour = 223
    case ('RED');     colour =  31
    case ('YELLOW');  colour =  63
    case ('GREEN');   colour =  95
    case ('CYAN');    colour =  127
    case ('BLUE');    colour =  159
    case ('MAGENTA'); colour =  191
    case ('WHITE');   colour =  0
    case default
       colour = 223
end select

return
END function SPECIAL_CONSTANTS_colour_to_number8
!******************************************************************************
function SPECIAL_CONSTANTS_colour_to_number24(colour_name) RESULT (colour)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'SPECIAL_CONSTANTS_colour_to_number24'

!Winteracter RGB
!WRGB = red + green*256 + blue*256*256
!where red, green and blue are 0 - 255
!The result is a 24 bit colour
!Upcase must be called before entry

INTEGER           :: colour
CHARACTER (LEN=*) :: colour_name

select case (colour_name)
    case ('BLACK');   colour =  0
    case ('RED');     colour =  255
    case ('YELLOW');  colour =  65535
    case ('GREEN');   colour =  65280
    case ('CYAN');    colour =  16776960
    case ('BLUE');    colour =  16711680
    case ('MAGENTA'); colour =  16711935
    case ('WHITE');   colour =  16777215
    case default
       colour = 0
end select

return
END function SPECIAL_CONSTANTS_colour_to_number24
!******************************************************************************
recursive function SPECIAL_CONSTANTS_next_colour24(colour_name, not_colour) RESULT (colour_next)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'SPECIAL_CONSTANTS_next_colour24'

!Winteracter RGB
!WRGB = red + green*256 + blue*256*256
!where red, green and blue are 0 - 255
!The result is a colour name for use above
!Upcase must be called before entry

CHARACTER (LEN=*), INTENT(IN)           :: colour_name
CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: not_colour   !e.g. do not return the background colour
CHARACTER (LEN=7)                       :: colour_next

select case (colour_name)
    case ('BLACK');   colour_next =  'RED    '
    case ('RED');     colour_next =  'YELLOW '
    case ('YELLOW');  colour_next =  'GREEN  '
    case ('GREEN');   colour_next =  'CYAN   '
    case ('CYAN');    colour_next =  'BLUE   '
    case ('BLUE');    colour_next =  'MAGENTA'
    case ('MAGENTA'); colour_next =  'WHITE  '  !bit dangerous on a white background!
    case ('WHITE');   colour_next =  'BLACK  '
    case default
       colour_next = 'BLACK'
end select

IF(PRESENT(not_colour)) then
   IF(colour_next == not_colour) colour_next = SPECIAL_CONSTANTS_next_colour24(colour_next)
endif

return
END function SPECIAL_CONSTANTS_next_colour24
!******************************************************************************
subroutine SPECIAL_CONSTANTS_test(reply)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'SPECIAL_CONSTANTS_test'
CHARACTER (LEN=reply_len), DIMENSION(reply_dim), INTENT(OUT) :: reply
CHARACTER (len=1)                                            :: letter
INTEGER                                                      :: pos, lenc

WRITE(reply(1),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)        'Special_Constants supplies the following constant'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(2),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg)  'ZERO     = ', ZERO
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(3),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg)  'ONE      = ', ONE
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(4),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg)  'TWO      = ', TWO
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(5),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg)  'THREE    = ', THREE
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(6),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg)  'FOUR     = ', FOUR
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(7),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg)  'FIVE     = ', FIVE
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(8),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg)  'SIX      = ', SIX
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(9),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg)  'SEVEN    = ', SEVEN
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(10),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'EIGHT    = ', EIGHT
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(11),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'NINE     = ', NINE
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(12),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'TEN      = ', TEN
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(13),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'HALF     = ', HALF
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(14),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'SIXTY    = ', SIXTY
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(15),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'HUNDRED  = ', HUNDRED
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(16),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'INFINITY = ', INFINITY
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(17),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'PI              = ', PI
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(18),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'two_pi          = ', two_pi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(19),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'log2            = ', log2
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(20),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'log_infinity    = ', log_infinity
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(21),FMT='(A,F11.4)', iostat=my_iostat, iomsg=my_iomsg) 'one_by_sqrt_2Pi = ', one_by_sqrt_2Pi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(22),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)       'equals  = '//equals
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(23),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)       'CRLF    = '//CRLF
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(24),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)       'degC    = '//degC
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(25),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)       'squared = '//squared
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(26),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)       'OFF     = '//OFF
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(27),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)       'ON      = '//ON
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(28),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)       'ON_OFF  = '//ON_OFF(1)//' '//ON_OFF(2)
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(29),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)       'YES_NO  = '//YES_NO(1)//' '//YES_NO(2)
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(30),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)       'Y_N     = '//Y_N(1)//' '//Y_N(2)
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(31),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'red         = '// red
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(32),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'black       = '// black
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(33),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'white       = '// white
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(34),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'yellow      = '// yellow
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(35),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'blue        = '// blue
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(36),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'grey        = '// grey
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(37),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'green       = '// green
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(38),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'light_blue  = '// light_blue
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(39),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'dark_blue   = '// dark_blue
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(40),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'dark_green  = '// dark_green
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(41),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'light_green = '// light_green
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(42),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'pink        = '// pink
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(43),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'violet      = '// violet
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(44),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'orange      = '// orange
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(45),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'dark_brown  = '// dark_brown
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(46),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'light_brown = '// light_brown
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(47),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'brown       = '// brown
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(48),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'turquoise   = '// turquoise
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(49),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'lemon       = '// lemon
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(50),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'salmon      = '// salmon
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(51),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'charcoal    = '// charcoal
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(52),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'purple      = '// purple
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(53),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'maroon      = '// maroon
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(54),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'mauve       = '// mauve
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(55),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'olive       = '// olive
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(56),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'Special characters for HTML'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(57),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_copy   = '//html_copy
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(58),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_reg    = '//html_reg
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(59),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_deg    = '//html_deg
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(60),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_pmi    = '//html_pmi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(61),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_sqd    = '//html_sqd
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(62),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_cbd    = '//html_cbd
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(63),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_qtr    = '//html_qtr
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(64),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_hlf    = '//html_hlf
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(65),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_3qtr   = '//html_3qtr
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(66),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_pd     = '//html_pd
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(67),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_sqrt   = '//html_sqrt
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(68),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_inf    = '//html_inf
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(69),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_int    = '//html_int
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(70),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_aeq    = '//html_aeq
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(71),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_neq    = '//html_neq
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(72),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_leq    = '//html_leq
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(73),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_geq    = '//html_geq
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(74),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_sum    = '//html_sum
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(75),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_sep    = '//html_sep
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(76),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'HTML Greek Alphabet'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(77),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_alpha   = '//html_alpha
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(78),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_beta    = '//html_beta
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(79),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_gamma   = '//html_gamma
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(80),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_delta   = '//html_delta
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(81),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_epsilon = '//html_epsilon
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(82),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_zeta    = '//html_zeta
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(83),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_eta     = '//html_eta
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(84),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_theta   = '//html_theta
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(85),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_iota    = '//html_iota
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(86),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_kappa   = '//html_kappa
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(87),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_lambda  = '//html_lambda
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(88),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_mu      = '//html_mu
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(89),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_nu      = '//html_nu
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(90),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_xi      = '//html_xi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(91),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_omicron = '//html_omicron
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(92),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_pi      = '//html_pi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(93),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_rho     = '//html_rho
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(94),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_sigma   = '//html_sigma
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(95),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_tau     = '//html_tau
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(96),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_upsilon = '//html_upsilon
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(97),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_phi     = '//html_phi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(98),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_chi     = '//html_chi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(99),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_psi     = '//html_psi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(100),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_omega   = '//html_omega
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(101),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_alpha   = '//html_u_alpha
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(102),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_beta    = '//html_u_beta
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(103),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_gamma   = '//html_u_gamma
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(104),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_delta   = '//html_u_delta
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(105),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_epsilon = '//html_u_epsilon
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(106),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_zeta    = '//html_u_zeta
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(107),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_eta     = '//html_u_eta
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(108),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_theta   = '//html_u_theta
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(109),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_iota    = '//html_u_iota
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(110),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_kappa   = '//html_u_kappa
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(111),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_lambda  = '//html_u_lambda
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(112),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_mu      = '//html_u_mu
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(113),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_nu      = '//html_u_nu
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(114),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_xi      = '//html_u_xi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(115),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_omicron = '//html_u_omicron
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(116),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_pi      = '//html_u_pi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(117),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_rho     = '//html_u_rho
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(118),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_sigma   = '//html_u_sigma
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(119),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_tau     = '//html_u_tau
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(120),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_upsilon = '//html_u_upsilon
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(121),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_phi     = '//html_u_phi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(122),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_chi     = '//html_u_chi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(123),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_psi     = '//html_u_psi
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(124),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'html_u_omega   = '//html_u_omega
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(125),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'Official Web Colours'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(126),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_aqua    = '//web_aqua
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(127),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_black   = '//web_black
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(128),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_white   = '//web_white
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(129),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_silver  = '//web_silver
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(130),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_yellow  = '//web_yellow
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(131),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_teal    = '//web_teal
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(132),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_fushia  = '//web_fushia
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(133),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_red     = '//web_red
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(134),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_blue    = '//web_blue
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(135),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_green   = '//web_green
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(136),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_gray    = '//web_gray
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(137),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_lime    = '//web_lime
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(138),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_maroon  = '//web_maroon
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(139),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_navy    = '//web_navy
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(140),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_olive   = '//web_olive
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(141),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'web_purple  = '//web_purple
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(142),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'Some useful fonts'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(143),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'Arial              = '//Arial
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(144),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'Courier_New        = '//Courier_New
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(145),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'Times_New_Roman    = '//Times_New_Roman
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(146),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'Symbol             = '//Symbol
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(147),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'Helvetica_95_Black = '//Helvetica_95_Black
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(148),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg) 'Helvetica_45_Light = '//Helvetica_45_Light
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(149),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'The Public Functions are:'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(150),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'IsSmall(x)               '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(151),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'IsExpBig(x)              '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(152),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'IsBetween(x,x1,x2)       '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(153),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'The Public Subroutines are:'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(154),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'SPECIAL_CONSTANTS_initialise_verbose'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(155),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'SPECIAL_CONSTANTS_initialise        '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
WRITE(reply(156),FMT='(A)', iostat=my_iostat, iomsg=my_iomsg)    'SPECIAL_CONSTANT_test               '
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif

!now test the functions

write(reply(157),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg)    'Some function tests...'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
write(reply(158),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg)    'SPECIAL_CONSTANTS_html_constant_found'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
write(reply(159),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg)    'input = " sadfkjhskdjfh*^%&^$%$££ "'
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif
write(reply(160),fmt='(A,L3)', iostat=my_iostat, iomsg=my_iomsg) 'output = ', &
                            SPECIAL_CONSTANTS_html_constant_found(' sadfkjhskdjfh*^%&^$%$££ ')
if(my_iostat /= 0) then
   call ERROR_MESSAGE_QUEUE_push(my_iomsg, place, my_iostat, 'write'); return
endif

write(reply(161),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg)    'input = "test &beta" '
write(reply(162),fmt='(A,L3)', iostat=my_iostat, iomsg=my_iomsg) 'output = ', &
                            SPECIAL_CONSTANTS_html_constant_found('test &beta', pos, lenc, letter)
write(reply(163),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg)    'letter = b = '//letter
write(reply(164),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'pos = ', pos
write(reply(165),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg)    'input = "&gamma; test" '
write(reply(166),fmt='(A,L3)', iostat=my_iostat, iomsg=my_iomsg) 'output = ', &
                            SPECIAL_CONSTANTS_html_constant_found('&gamma; test', pos, lenc, letter)
write(reply(167),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg)    'letter = g = '//letter
write(reply(168),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'pos = ', pos
write(reply(169),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'lenc = ', lenc

write(reply(170),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg)    'input = "&#178; test" '
write(reply(171),fmt='(A,L3)', iostat=my_iostat, iomsg=my_iomsg) 'output = ', &
                            SPECIAL_CONSTANTS_html_constant_found('&#178;; test', pos, lenc, letter)
write(reply(172),fmt='(A)', iostat=my_iostat, iomsg=my_iomsg)    'letter = '//char(253)//' = '//letter
write(reply(173),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'pos = ', pos
write(reply(174),fmt='(A,I0)', iostat=my_iostat, iomsg=my_iomsg) 'lenc = ', lenc
reply(175) = ' '
reply(176) = 'Calling SPECIAL_CONSTANTS_font_for_html_constant'
write(reply(177),*, iostat=my_iostat, iomsg=my_iomsg) 'SPECIAL_CONSTANTS_font_for_html_constant(html_sqd) = ', &
                     SPECIAL_CONSTANTS_font_for_html_constant(html_sqd)
write(reply(178),*, iostat=my_iostat, iomsg=my_iomsg) 'SPECIAL_CONSTANTS_font_for_html_constant(html_alpha) = ', &
                     SPECIAL_CONSTANTS_font_for_html_constant(html_alpha)
reply(179) = ' '

return
END subroutine SPECIAL_CONSTANTS_test
!******************************************************************************
subroutine SPECIAL_CONSTANTS_test_get_reply_size(mydim, mylen)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'SPECIAL_CONSTANTS_test_get_reply_size'
INTEGER, INTENT(OUT) :: mydim
INTEGER, INTENT(OUT) :: mylen

mydim = reply_dim
mylen = reply_len

return
end subroutine SPECIAL_CONSTANTS_test_get_reply_size
!******************************************************************************
subroutine SPECIAL_CONSTANTS_initialise_verbose(verb)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'SPECIAL_CONSTANTS_initialise_verbose'
logical, intent(IN), dimension(verbose_levels) :: verb

verbose = verb

return
end subroutine SPECIAL_CONSTANTS_initialise_verbose
!******************************************************************************
subroutine SPECIAL_CONSTANTS_get_verbose(verb)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'SPECIAL_CONSTANTS_get_verbose'
logical, intent(out), dimension(verbose_levels) :: verb

verb = verbose

return
end subroutine SPECIAL_CONSTANTS_get_verbose
!******************************************************************************
subroutine SPECIAL_CONSTANTS_initialise(message)
!<cat>
!
!</cat>
implicit none
character (:), parameter:: place = 'SPECIAL_CONSTANTS_initialise'
CHARACTER (len=*), OPTIONAL :: message

log2            = LOG(two)
pi              = asin(one)*two
two_pi          = two*pi
pi_by_2         = pi/two
log_infinity    = log(infinity)-one   !just changed to -1 28/2/06 NFK
one_by_sqrt_2Pi = one/sqrt(two_Pi)
sqrt2           = SQRT(two)
root_pi_by_2    = SQRT(pi)/two

if(verbose(v_top).and.present(message)) message='SPECIAL_CONSTANTS_initialise initialised'

return
end subroutine SPECIAL_CONSTANTS_initialise
!******************************************************************************

end module special_constants
