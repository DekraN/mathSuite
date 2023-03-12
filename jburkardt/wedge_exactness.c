#ifndef __DISABLEDEEP_WEDGEEXACTNESS

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _wedge01_integral ( void * data)
/******************************************************************************/
/*
  Purpose:
    WEDGE01_INTEGRAL returns the integral of a monomial in the unit wedge in 3D.
  Discussion:
    This routine returns the integral of
      product ( 1 <= I <= 3 ) X(I)^E(I)
    over the unit wedge.
    The integration region is:
      0 <= X
      0 <= Y
      X + Y <= 1
      -1 <= Z <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 August 2014
  Author:
    John Burkardt
  Reference:
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971,
    ISBN: 0130438936,
    LC: QA311.S85.
  Parameters:
    Input, int E[3], the exponents.
    Output, double WEDGE01_INTEGRAL, the integral of the monomial.
*/
{
	static ityp result = MAX_VAL;
	
	int * e = data;
	
    ityp value = 1.0;
    dim_typ i, k = e[0];

    for ( i = 1; i <= e[1]; ++i)
    {
        ++ k;
        value *= ( ityp ) ( i ) / ( ityp ) ( k );
    }

    ++ k;
    value /= ( ityp ) ( k );

    ++ k;
    value /= ( ityp ) ( k );
    /*
    Now account for integration in Z.
    */
    
    result = e[2] == -1 ? MAX_VAL : ( e[2] % 2 ) == 1 ? 0.00 : value*(2.00 / ( ityp ) ( e[2] + 1 ));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _wedge01_volume ( void * data)
/******************************************************************************/
/*
  Purpose:
    WEDGE01_VOLUME returns the volume of the unit wedge in 3D.
  Discussion:
    The unit wedge is:
      0 <= X
      0 <= Y
      X + Y <= 1
      -1 <= Z <= 1.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 August 2014
  Author:
    John Burkardt
  Parameters:
    Output, double WEDGE01_VOLUME, the volume of the unit wedge.
*/
{
	static ityp result = 1.00;
    return &result;
}

#endif
