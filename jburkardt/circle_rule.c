#ifndef __DISABLEDEEP_CIRCLERULE

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _circle_rule ( void * data)
/******************************************************************************/
/*
  Purpose:
    CIRCLE_RULE computes a quadrature rule for the unit circle.
  Discussion:
    The unit circle is the region:
      x * x + y * y = 1.
    The integral I(f) is then approximated by
      Q(f) = 2 * M_PI * sum ( 1 <= i <= NT ) W(i) * F ( cos(T(i)), sin(T(i)) ).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    06 April 2014
  Author:
    John Burkardt
  Parameters:
    Input, int NT, the number of angles to use.
    Output, double W[NT], the weights for the rule.
    Output, double T[NT], the angles for the rule.
*/
{
	const dt2pit * const s_data = data;
	const register dim_typ nt = s_data->a0;
	ityp * w = s_data->a1;
	ityp * t = s_data->a2;
	
    for (dim_typ it = 0; it < nt; ++it )
    {
        w[it] = 1.00 / ( ityp ) ( nt );
        t[it] = M_2TPI * ( ityp ) ( it ) / ( ityp ) ( nt );
    }
    return NULL;
}

#endif
