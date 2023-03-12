#ifndef __DISABLEDEEP_UPC

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  int   upc_check_digit_calculate ( char *s )
/******************************************************************************/
/*
  Purpose:
    UPC_CHECK_DIGIT_CALCULATE returns the check digit of a UPC.
  Discussion:
    UPC stands for Universal Product Code.
    A full UPC is a string of 12 digits, in groups of size 1, 5, 5, and 1,
    of the form P-LLLLL-RRRRR-C, where:
      P is the one-digit product type code.
      L is the five-digit manufacturer code.
      R is the five_digit product code
      C is the check digit.
  Example:
    0-72890-00011-8
    0-12345-67890-5
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 May 2001
  Author:
    John Burkardt
  Reference:
    David Savir, George Laurer,
    The Characteristics and Decodability of the Universal Product Code,
    IBM Systems Journal,
    Volume 14, Number 1, pages 16-34, 1975.
  Parameters:
    Input, char *S, a string containing at least 11 digits.
    Dashes and other characters will be ignored.  A 12th digit may be
    included, but it will be ignored.
    Output, int UPC_CHECK_DIGIT_CALCULATE, the check digit.
*/
{
    int d;
    int *dvec = s_to_digits ( s, 11 );
    d = 3 * ( dvec[0] + dvec[2] + dvec[4] + dvec[6] + dvec[8] + dvec[10] )+ ( dvec[1] + dvec[3] + dvec[5] + dvec[7] + dvec[9] );
    d = ( d % 10 );
    d = ( ( 10 - d ) % 10 );
    free ( dvec );
    return d;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  bool   upc_is_valid ( char *s )
/******************************************************************************/
/*
  Purpose:
    UPC_IS_VALID reports whether a UPC is valid.
  Discussion:
    UPC stands for Universal Product Code.
    A full UPC is a string of 12 digits, in groups of size 1, 5, 5, and 1,
    of the form P-LLLLL-RRRRR-C, where:
      P is the one-digit product type code.
      L is the five-digit manufacturer code.
      R is the five_digit product code
      C is the check digit.
  Example:
    0-72890-00011-8
    0-12345-67890-5
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 September 2015
  Author:
    John Burkardt
  Reference:
    David Savir, George Laurer,
    The Characteristics and Decodability of the Universal Product Code,
    IBM Systems Journal,
    Volume 14, Number 1, pages 16-34, 1975.
  Parameters:
    Input, char *S, a string containing 12 digits.
    Dashes and other characters will be ignored.
    Output, int UPC_IS_VALID, is TRUE if the string
    is a valid UPC.
*/
{
    int d1;
    int d2;
    int *dvec;
    int value;

    dvec = s_to_digits ( s, 12 );

    d1 = upc_check_digit_calculate ( s );
    d2 = dvec[11];
    free ( dvec );
    return d1==d2;
}
#endif
