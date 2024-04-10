
module mod_kndconsts

    implicit none

    integer, parameter :: kind_short = 2
    integer, parameter :: kind_long  = 4
#ifdef INTEGER_SHORT
    integer, parameter :: kind_int = kind_short
#else
    integer, parameter :: kind_int = kind_long
#endif

    integer, parameter :: kind_single = 4
    integer, parameter :: kind_double = 8
#ifdef REAL_SINGLE
    integer, parameter :: kind_real = kind_single
#else
    integer, parameter :: kind_real = kind_double
#endif

#ifdef FORTRAN_2003_STANDARD
    integer, parameter :: kind_char = 4
#endif

    integer, parameter :: len_char_name = 32
    integer, parameter :: len_char_msg  = 128
    integer, parameter :: len_char_file = 128

end module mod_kndconsts
