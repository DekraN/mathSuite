#ifndef __DISABLEDEEP_ISBN

#include "../dutils.h"



/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ch_is_digit ( void * data)
/******************************************************************************/
/*
  Purpose:
    CH_IS_DIGIT returns TRUE if a character is a decimal digit.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 January 2014
  Author:
    John Burkardt
  Parameters:
    Input, char C, the character to be analyzed.
    Output, int CH_IS_DIGIT, is TRUE if C is a digit.
*/
{
	static bool result = 2;
	
	register char c = *(char *) data;
	
	result = '0' <= c && c <= '9';
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ch_to_digit ( void * data)
/******************************************************************************/
/*
  Purpose:
    CH_TO_DIGIT returns the integer value of a base 10 digit.
  Example:
     CH  DIGIT
    ---  -----
    '0'    0
    '1'    1
    ...  ...
    '9'    9
    ' '    0
    'X'   -1
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    13 June 2003
  Author:
    John Burkardt
  Parameters:
    Input, char CH, the decimal digit, '0' through '9'.
    Output, int CH_TO_DIGIT, the corresponding integer value.  If the
    character was 'illegal', then DIGIT is -1.
*/
{
	static int result = INT_MAX;
	
	register char ch = *(char *) data;
	
	result = '0' <= ch && ch <= '9' ? ch-'0' : -1;
    return &result;
}


/******************************************************************************/
__MATHSUITE __JBURKARDT  void   * _s_to_digits ( void * data)
/******************************************************************************/
/*
  Purpose:
    S_TO_DIGITS extracts N digits from a string.
  Discussion:
    The string may include spaces, letters, and dashes, but only the
    digits 0 through 9 will be extracted.
  Example:
    S  => 34E94-70.6
    N  => 5
    D <= (/ 3, 4, 9, 4, 7 /)
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    09 September 2015
  Author:
    John Burkardt
  Parameters:
    Input, char *S, the string.
    Input, int N, the number of digits to extract.
    Output, int S_TO_DIGITS[N], the extracted digits.
*/
{
	const pcdt * const s_data = data;
	char * s = s_data->a0;
	const register dim_typ n = s_data->a1;
	
    char c;
    dim_typ d;
    dim_typ d_pos;
    int *dvec;
    dim_typ s_len;
    dim_typ s_pos;

    dvec = ( int * ) malloc ( n * sizeof ( int ) );

    s_len = strlen ( s );

    s_pos = 0;
    d_pos = 0;

    while ( d_pos < n )
    {

        if ( s_len <= s_pos )
            return NULL;

        c = s[s_pos];
        ++ s_pos;

        if ( ch_is_digit ( c ) )
        {
            d = ch_to_digit ( c );
            dvec[d_pos] = d;
            ++ d_pos;
        }
    }

    return dvec;
}

#endif
