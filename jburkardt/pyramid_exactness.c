#ifndef __DISABLEDEEP_PYRAMIDEXACTNESS

#include "../dutils.h"

 /******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _pyra_unit_monomial ( void * data)
/******************************************************************************/
/*
  Purpose:
    PYRA_UNIT_MONOMIAL: monomial integral in a unit pyramid.
  Discussion:
    This function returns the value of the integral of X^ALPHA Y^BETA Z^GAMMA
    over the unit pyramid.
    The integration region is:
    - ( 1 - Z ) <= X <= 1 - Z
    - ( 1 - Z ) <= Y <= 1 - Z
              0 <= Z <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 March 2008
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Input, int EXPON[3], the exponents.
    Output, double PYRA_UNIT_MONOMIAL, the integral of the monomial
    over the pyramid.
*/
{
	static ityp result = MAX_VAL;
	
	int * expon = data;
	
    dim_typ i;
    int i_hi;
    ityp value = 0.00;

    if ( ( expon[0] % 2 ) == 0 && ( expon[1] % 2 ) == 0 )
    {
        i_hi = 2 + expon[0] + expon[1];

        for ( i = 0; i <= i_hi; ++i)
            value += r8_mop ( i ) * r8_choose ( i_hi, i )/ ( ityp ) ( i + expon[2] + 1 );

        value *= 4.00 * (1.00 /(( ityp ) ( expon[0] + 1 ) * ( ityp ) ( expon[1] + 1 )));
    }

	result = value; 
    return &result;
}

#endif
