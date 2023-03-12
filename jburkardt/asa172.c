#ifndef __DISABLEDEEP_ASA172

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _revers( void * data)
/******************************************************************************/
/*
  Purpose:
    REVERS reorders the subscript vector, if required.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 July 2008
  Author:
    M O'Flaherty, G MacKenzie
  Reference:
    M O'Flaherty, G MacKenzie,
    Algorithm AS 172:
    Direct Simulation of Nested Fortran DO-LOOPS,
    Applied Statistics,
    Volume 31, Number 1, 1982, pages 71-74.
  Parameters:
    Input/output, int IVEC[KDIM], the subscript vector.
    Input, int KDIM, the dimension of the subscript vector.
*/
{
	const dtpit * const s_data = data;
	const register dim_typ kdim = s_data->a0;
	ityp * ivec = s_data->a1;
	
	dim_typ itemp;
	for (dim_typ i = 0; i < (kdim>>1); ++i )
	{
		itemp          = ivec[i];
		ivec[i]        = ivec[kdim-1-i];
		ivec[kdim-1-i] = itemp;
	}
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _simdo ( void * data)
/******************************************************************************/
/*
  Purpose:
    SIMDO generates multi-indices, simulating nested DO-loops.
  Discussion:
    The loops are assumed to be nested to a depth of K.
    The R-th loop is assumed to have upper limit N(R) and increment Inc(R).
    The total number of executions of the innermost loop is
      N = product ( 1 <= R <= K ) N(R).
    Let these executions be indexed by the single integer J, which
    we call the index subscript.
    Each value of J corresponds to a particular set of loop indices,
    which we call the subscript vector I(J).
    This routine can start with J and find I(J), or determine
    J from I(J).
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    27 July 2008
  Author:
    M O'Flaherty, G MacKenzie
  Reference:
    M O'Flaherty, G MacKenzie,
    Algorithm AS 172:
    Direct Simulation of Nested Fortran DO-LOOPS,
    Applied Statistics,
    Volume 31, Number 1, 1982, pages 71-74.
  Parameters:
    Input, int QIND.
    TRUE to convert an index subscript J to the subscript vector I(J).
    FALSE to convert the subscript vector I(J) to the index subscript J.
    Input, int QFOR,
    TRUE if conversion is required in standard Fortran subscripting order,
    FALSE otherwise.
    Input, int IPROD[KDIM], contains the partial products.
    If QFOR is FALSE, then
      IPROD(S) = product ( 1 <= R <= S ) N(R).
    If QFOR is TRUE, then
      IPROD(S) = product ( 1 <= R <= S ) N(KDIM+1-R).
    Input, int KDIM, the nesting depth of the loops.
    Input/output, int *JSUB.
    If QIND is TRUE, then JSUB is an input quantity, an index subscript J
    to be converted into the subscript vector I(J).
    If QIND is FALSE, then JSUB is an output quantity, the index subscript J
    corresponding to the subscript vector I(J).
    Input/output, int IVEC[KDIM].
    if QIND is TRUE, then IVEC is an output quantity, the subscript vector I(J)
    corresponding to the index subscript J.
    If QIND is FALSE, then IVEC is an input quantity, a subscript vector I(J)
    for which the corresponding index subscript J is to be computed.
    Output, int SIMDO, error flag.
    0, no error was detected.
    1, if QIND is TRUE, and the input value of JSUB exceeds IPROD(KDIM).
    2, if QIND is FALSE, and IVEC contains an illegal component.
*/
{
	static sel_typ result = UCHAR_MAX;
	
	const _2bdtpitpdtpit * const s_data = data;
	const bool qind = s_data->a0;
	const bool qfor = s_data->a1;
	const register dim_typ kdim = s_data->a2;
	ityp * iprod = s_data->a3;
	dim_typ * jsub = s_data->a4;
	ityp * ivec = s_data->a5;
	
	dim_typ i;
	int ik;
	int itempv;

	/*
	Index subscript to subscript vector conversion.
	*/
	if ( qind )
	{
		if ( iprod[kdim-1] < *jsub )
		{
			result = SIMDO_JSUMEXCEEDEDIPROD;
			return &result;
		}

		itempv = *jsub - 1;

		for ( i = 0; i < kdim - 1; i++ )
		{
			ik = kdim - 2 - i;
			ivec[i] = itempv / iprod[ik];
			itempv -= iprod[ik] * ivec[i];
			++ ivec[i];
		}

		ivec[kdim-1] = itempv + 1;
		if ( qfor )
			revers (kdim, ivec);
	}
	/*
	Subscript vector to index subscript conversion.
	*/
	else
	{
		if ( !qfor )
			revers(kdim, ivec);

	if ( iprod[0] < ivec[0] )
	{
		result = SIMDO_IVECILLEGALCOMPO;
		return &result;
	}

		for ( i = 1; i < kdim; ++i )
			if ( iprod[i] / iprod[i-1] < ivec[i] )
			{
				result = SIMDO_IVECILLEGALCOMPO;
				return &result;
			}

		*jsub = ivec[0];
		for(i = 1; i < kdim; ++i)
			*jsub += ( ivec[i] - 1 ) * iprod[i-1];
		/*
		As a courtesy to the caller, UNREVERSE the IVEC vector
		if you reversed it.
		*/
		if ( !qfor )
			revers (kdim, ivec);
	}

	result = SIMDO_SUCCESS;
	return &result;
}

#endif
