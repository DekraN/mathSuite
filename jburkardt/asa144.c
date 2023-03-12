#ifndef __DISABLEDEEP_ASA144

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rcont ( void * data) 
/******************************************************************************/
/*
  Purpose:
    RCONT generates a random two-way table with given marginal totals.
  Discussion:
    Each time the program is called, another table will be randomly
    generated.
    Note that it should be the case that the sum of the row totals
    is equal to the sum of the column totals.  However, this program
    does not check for that condition.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 November 2010
  Author:
    Original FORTRAN77 version by James Boyett.
    C version by John Burkardt.
  Reference:
    James Boyett,
    Algorithm AS 144:
    Random R x C Tables with Given Row and Column Totals,
    Applied Statistics,
    Volume 28, Number 3, pages 329-332, 1979.
  Parameters:
    Input, int NROW, the number of rows in the observed matrix.
    Input, int NCOL, the number of columns in the observed matrix.
    Input, int NROWT[NROW], the row totals of the observed matrix.
    Input, int NCOLT[NCOL], the column totals of the observed matrix.
    Input/output, int NSUBT[NCOL], used by RCONT for partial column sums.
    Must not be changed by the calling program.
    Output, int MATRIX[NROW*NCOL], the randomly generated matrix.
    Input/output, int *KEY, should be set to FALSE by the user before
    the initial call.  RCONT will reset it to TRUE, and it should be left
    at that value for subsequent calls in which the same values of NROW,
    NCOL, NROWT and NCOLT are being used.
    Output, int *IFAULT, fault indicator.
    0, no error occurred.
    1, NROW <= 0.
    2, NCOL <= 1.
    3, some entry of NROWT is less than 0.
    4, some entry of NCOLT is less than 0.
*/
{
	static bool result = 2;
	
	const _4pdtpitpb * const s_data = data;
	const dim_typ * dim = s_data->a0;
	dim_typ * nrowt = s_data->a1;
	dim_typ * ncolt = s_data->a2;
	dim_typ * nsubt = s_data->a3;
	ityp * matrix = s_data->a4;
	bool * key = s_data->a5;
	
	dim_typ i;
	dim_typ ii;
	dim_typ j;
	dim_typ k;
	dim_typ limit;
	dim_typ noct;
	dim_typ ntemp;
	static dim_typ ntotal;
	static dim_typ *nvect = NULL;
	static unsigned seed = 0;

	if ( !(*key) )
	{
		/*
		Set KEY for subsequent calls.
		*/
		*key = true;
		seed = 123456789;
		/*
		Check for faults and prepare for future calls.
		*/
		if ( dim[ROWS] <= 0 || dim[COLUMNS] <= 1 || ncolt[0] <= 0)
		{
			result = RCONT_INVALIDPARAMS;
			return &result;
		}

		for (i = 0; i < dim[ROWS]; ++i)
			if ( nrowt[i] <= 0 )
			{
				result = RCONT_INVALIDPARAMS;
				return &result;
			}

		nsubt[0] = ncolt[0];

		for ( j = 1; j < dim[COLUMNS]; ++j )
		{
			if(ncolt[j] <= 0 )
			{
				result = RCONT_INVALIDPARAMS;
				return &result;
			}
			nsubt[j] = nsubt[j-1] + ncolt[j];
		}

		ntotal = nsubt[dim[COLUMNS]-1];
		nvect = ( dim_typ * ) malloc ( ntotal * sizeof ( dim_typ ) );
		/*
		Initialize vector to be permuted.
		*/
		#pragma omp parallel for
		for (i = 0; i < ntotal; ++i )
			nvect[i] = i + 1;
	}
	/*
	Initialize vector to be permuted.
	*/
	dim_typ nnvect[ntotal];

	#pragma omp parallel for
	for ( i = 0; i < ntotal; ++i )
		nnvect[i] = nvect[i];
	/*
	Permute vector.
	*/
	ntemp = ntotal;

	for (i = 0; i < ntotal; ++i )
	{
		noct = (dim_typ) ( r8_uniform_01 ( &seed ) * ( ityp ) ( ntemp ) + 1.00 );
		nvect[i] = nnvect[noct-1];
		nnvect[noct-1] = nnvect[ntemp-1];
		-- ntemp;
	}
	/*
	Construct random matrix.
	*/
	#pragma omp parallel for
	for (j = 0; j < dim[COLUMNS]; ++j)
		#pragma omp parallel for
		for (i = 0; i < dim[ROWS]; ++i )
			matrix[i*dim[ROWS]+j] = 0.00;

	ii = 0;

	for ( i = 0; i < dim[ROWS]; ++i )
	{
		limit = nrowt[i];

		for ( k = 0; k < limit; ++k )
		{
			for ( j = 0; j < dim[COLUMNS]; ++j )
			{
				if ( nvect[ii] <= nsubt[j] )
				{
					++ ii;
					matrix[i+j*dim[ROWS]] = matrix[i+j*dim[ROWS]] + 1;
					break;
				}
			}
		}
	}

	result = RCONT_SUCCESS;
	return &result;
}

#endif
