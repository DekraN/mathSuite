#ifndef __DISABLEDEEP_WTIME

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _wtime ( void * data)
/******************************************************************************/
/*
  Purpose:

    WTIME reports the elapsed wallclock time.
  Discussion:
    The reliability of this function depends in part on the value of
    CLOCKS_PER_SECOND.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 April 2009
  Author:
    John Burkardt
  Parameters:
    Output, double WTIME, the a reading of the wall clock timer,
    in seconds.
*/
{
	static ityp result = MAX_VAL;
	
	result = ( ityp ) clock ( ) / ( ityp ) CLOCKS_PER_SEC;
    return &result;
}

#endif
