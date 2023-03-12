#ifndef __DISABLEDEEP_ASA183

#include "../dutils.h"

/*******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8_random ( void * data)
/*******************************************************************************/
/*
  Purpose:
    r8_RANDOM returns a pseudorandom number between 0 and 1.
  Discussion:
    This function returns a pseudo-random number rectangularly distributed
    between 0 and 1.   The cycle length is 6.95E+12. (See page 123
    of Applied Statistics (1984) volume 33), not as claimed in the
    original article.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 July 2008
  Author:
    Original FORTRAN77 version by Brian Wichman, David Hill.
    C++ version by John Burkardt.
  Reference:
    Brian Wichman, David Hill,
    Algorithm AS 183: An Efficient and Portable Pseudo-Random
    Number Generator,
    Applied Statistics,
    Volume 31, Number 2, 1982, pages 188-190.
  Parameters:
    Input/output, int *S1, *S2, *S3, three values used as the
    seed for the sequence.  These values should be positive
    integers between 1 and 30,000.
    Output, ityp r8_RANDOM, the next value in the sequence.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ ** a_data = data;
	dim_typ *s1 = a_data[0];
	dim_typ *s2 = a_data[1];
	dim_typ *s3 = a_data[2];
	
	*s1 = ( ( 171 * *s1 ) % 30269 );
	*s2 = ( ( 172 * *s2 ) % 30307 );
	*s3 = ( ( 170 * *s3 ) % 30323 );
	result = fmod ( ( ityp ) ( *s1 ) / 30269.0  + ( ityp ) ( *s2 ) / 30307.0  + ( ityp ) ( *s3 ) / 30323.0, 1.0 );
	return &result;
}

/*******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _r8_uni ( void * data)
/*******************************************************************************/
/*
  Purpose:
    r8_UNI returns a pseudorandom number between 0 and 1.
  Discussion:
    This function generates uniformly distributed pseudorandom numbers
    between 0 and 1, using the 32-bit generator from figure 3 of
    the article by L'Ecuyer.
    The cycle length is claimed to be 2.30584E+18.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 July 2008
  Author:
    Original PASCAL version by Pierre L'Ecuyer.
    C++ version by John Burkardt.
  Reference:
    Pierre LEcuyer,
    Efficient and Portable Combined Random Number Generators,
    Communications of the ACM,
    Volume 31, Number 6, June 1988, pages 742-751.
  Parameters:
    Input/output, int *S1, *S2, two values used as the
    seed for the sequence.  On first call, the user should initialize
    S1 to a value between 1 and 2147483562;  S2 should be initialized
    to a value between 1 and 2147483398.
    Output, ityp r8_UNI, the next value in the sequence.
*/
{
	static ityp result = MAX_VAL;
	
	dim_typ ** a_data = data;
	dim_typ * s1 = a_data[0]; 
	dim_typ * s2 = a_data[1];
	
	dim_typ k;
	unsigned z;

	k = *s1 / 53668;
	*s1 = 40014 * ( *s1 - k * 53668 ) - k * 12211;
	if ( *s1 < 0 )
		*s1 += 2147483563;
	k = *s2 / 52774;
	*s2 = 40692 * ( *s2 - k * 52774 ) - k * 3791;
	if ( *s2 < 0 )
		*s2 += 2147483399;
	z = *s1 - *s2;
	if ( z < 1 )
		z = 2147483562;
		
	result = ( ityp ) ( z ) / 2147483563.00;
    return &result;
}

#endif
