#ifndef __DISABLEDEEP_NORMAL

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  double   complex c8_normal_01 ( int *seed )
/******************************************************************************/
/*
  Purpose:
    C8_NORMAL_01 returns a unit pseudonormal C8.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    03 March 2015
  Author:
    John Burkardt
  Parameters:
    Input/output, int *SEED, a seed for the random
    number generator.
    Output, double complex C8_NORMAL_01, a unit pseudonormal value.
*/
{
    ityp v1;
    ityp v2;
    v1 = r8_uniform_01 ( seed );
    v2 = r8_uniform_01 ( seed );
    return sqrt ( - 2.00 * log ( v1 ) ) * cos ( M_2TPI * v2 ) + sqrt ( - 2.00 * log ( v1 ) ) * sin ( M_2TPI * v2 ) * I;
}
#endif
