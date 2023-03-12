#ifndef __DISABLEDEEP_ASA314

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _invmod ( void * data)
/******************************************************************************/
/*
  Purpose:
    INVMOD inverts a matrix using modulo arithmetic.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2013
  Author:
    Original FORTRAN77 version by Roger Payne.
    C version by John Burkardt.
  Reference:
    Roger Payne,
    Inversion of matrices with contents subject to modulo arithmetic,
    Applied Statistics,
    Volume 46, Number 2, 1997, pages 295-298.
  Parameters:
    Input/output, int MAT[NROW*NROW].
    On input, the matrix to be inverted.
    On output, the product of the input matrix and IMAT.
    Output, int IMAT[NROW*NROW], the inverse matrix.
    If IFAULT = -1 on output, then IMAT is only a left inverse.
    Input, int RMOD[NROW], the modulus for values in each row.
    Input, int CMOD[NROW], the modulus for values
    in each column.
    Input, int NROW, the order of the matrix.
    Output, int *IFAULT, an error flag.
    0, no error was detected.
    -1, only a left inverse could be formed.
    1, the matrix contains elements that are negative, or too large.
    2, the matrix contains nonzero elements in mixed modulus positions.
    3, the matrix cannot be inverted.
*/
{
	static bool result = 10;
	
	const dt2pit2pdt * const s_data = data;
	const register dim_typ nrow = s_data->a0;
	ityp * mat = s_data->a1;
	ityp * imat =  s_data->a2;
	dim_typ * rmod = s_data->a3;
	dim_typ * cmod = s_data->a4;
	
	int all_zero;
	int i;
	int ir;
	int j;
	int k;
	int kir;
	int kjr;
	int n;
	/*
	Check that elements in 'mixed-moduli' positions are all zero.
	*/
	n = 0;
	for ( i = 1; i <= nrow; ++i )
	{
		for ( j = 1; j <= nrow; ++j )
		{
			++ n;
			if ( ( rmod[i-1] != cmod[j-1] ) && ( 0 < mat[(n-1)*nrow] ) )
			{
				result = INVMOD_NONZEROELMSINMIXEDMODULUSPOSITIONS;
				return &result;	
			}
			if ( ( rmod[i-1] < mat[n-1] ) ||( mat[n-1] < 0 ) )
			{
				result = INVMOD_INVALIDELEMENTS;
				return &result;
			}	
		}
	}

	n = 0;
	for ( i = 1; i <= nrow; ++i )
	{
		for ( j = 1; j <= nrow; ++j )
		{
			++ n;
			imat[(n-1)*nrow] = 0;
		}
	}
	/*
	Sort rows and columns into ascending order of moduli
	*/

	ityp rsort[nrow];
	ityp csort[nrow];


	msort ( nrow, mat, imat, rmod, cmod, rsort, csort );
	/*
	Complete initialization of inverse matrix
	*/
	#pragma omp parallel for
	for ( n = 1; n <= nrow * nrow; n += nrow + 1 )
		imat[(n-1)*nrow] = 1;
	/*
	Invert the matrix.
	*/
	for ( ir = 1; ir <= nrow; ++ir )
	{
		kir = ( ir - 1 ) * nrow;

		if ( !mat[(kir+ir-1)*nrow])
		{
			/*
			Find a row JR below IR such that K(JR,IR)>0
			*/
			all_zero = 1;

			for ( kjr = kir + nrow + ir; kjr <= nrow * nrow; kjr += nrow )
			{
				if ( 0 < mat[(kjr-1)*nrow] )
				{
					all_zero = 0;
					break;
				}
			}
			/*
			Column IR contains all zeros in rows IR or below:
			look for a row above with zeros to left of column IR
			and K(JR,IR)>0
			*/
			if ( all_zero )
			{
				for ( kjr = ir; kjr <= kir; kjr += nrow )
				{
					if ( 0 < mat[(kjr-1)*nrow] )
					{
						for ( i = kjr - ir + 1; i < kjr; ++i )
						{
							if ( 0 < mat[(i-1)*nrow] )
							{
								result = INVMOD_NOTINVERTIBLEMATRIX;
								return &result;
							}	
						}
						all_zero = 0;
						break;
					}
				}
			}
			/*
			Column IR contains all zeros
			*/
			if ( all_zero )
				continue;
			/*
			Switch row JR with row IR
			*/
			kjr -= ir;

			for ( i = 1; i <= nrow; ++i )
			{
				k = mat[(kir+i-1)*nrow];
				mat[(kir+i-1)*nrow] = mat[(kjr+i-1)*nrow];
				mat[(kjr+i-1)*nrow] = k;

				k = imat[(kir+i-1)*nrow];
				imat[(kir+i-1)*nrow] = imat[(kjr+i-1)*nrow];
				imat[(kjr+i-1)*nrow] = k;
			}
		}
		/*
		Find a multiplier N such that N*MAT(IR,IR)=1 mod(P{IR})
		*/
		k = mat[(kir+ir-1)*nrow];
		for ( n = 1; n < rmod[ir-1]; n++ )
			if ( ( n * k ) % rmod[ir-1] == 1 )
				break;
		/*
		Multiply row IR by N.
		*/
		if ( 1 < n )
		{
			for ( i = kir + 1; i <= ir * nrow; i++ )
			{
			mat[i-1] = mat[i-1] * n;
			imat[i-1] = imat[i-1] * n;
			}
		}
		/*
		Subtract MAT(JR,IR) * row IR from each row JR
		*/
		for ( kjr = 0; kjr < nrow * nrow; kjr = kjr + nrow )
		{
			n = rmod[ir-1] - mat[(kjr+ir-1)*nrow];
			if ( ( kjr != kir ) && ( n != 0 ) )
				for ( i = 1; i <= nrow; ++i )
				{
					mat[(kjr+i-1)*nrow]  = (int64_t)(  mat[(kjr+i-1)*nrow] + n *  mat[(kir+i-1)*nrow] ) % cmod[i-1];
					imat[(kjr+i-1)*nrow] = (int64_t)( imat[(kjr+i-1)*nrow] + n * imat[(kir+i-1)*nrow] ) % cmod[i-1];
				}
		}

	}
	/*
	Check inversion was possible - that result has
	non-zero elements only on diagonal.
	*/
	/*
	If we encounter a zero diagonal element, then only a left inverse
	will be formed.
	*/
	#pragma omp parallel for
	for ( n = 1; n <= nrow * nrow; n += nrow + 1 )
		mat[n-1] *= -1;

	for ( n = 1; n <= nrow * nrow; ++n )
		if ( 0 < mat[n-1] )
		{
			result = INVMOD_NOTINVERTIBLEMATRIX;
			return &result;
		}
	#pragma omp parallel for
	for ( n = 1; n <= nrow * nrow; n += nrow + 1 )
		mat[n-1] *= -1;
	/*
	Unsort the rows and columns back into their original order.
	*/
	musort ( nrow, mat, imat, rmod, cmod, rsort, csort );

	result = INVMOD_SUCCESS;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _msort ( void * data)
/******************************************************************************/
/*
  Purpose:
    MSORT sorts matrix rows and columns in ascending order of moduli.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2013
  Author:
    Original FORTRAN77 version by Roger Payne.
    C version by John Burkardt.
  Reference:
    Roger Payne,
    Inversion of matrices with contents subject to modulo arithmetic,
    Applied Statistics,
    Volume 46, Number 2, 1997, pages 295-298.
  Parameters:
    Input/output, int MAT[NROW*NROW].
    On output, the matrix has been sorted.
    Ignoreput, int IMAT[NROW*NROW].
    This quantity is ignored.
    Input/output, int RMOD[NROW], the modulus for values in
    each row.  On output, these have been rearranged according to the sorting.
    Input/output, int CMOD[NROW], the modulus for values in
    each column.  On output, these have been rearranged according to the
    sorting.
    Output, int RSORT[NROW], the sorted row indices.
    Output, int CSORT[NROW], the sorted column indices.
    Input, int NROW, the order of the matrix.
*/
{
	const dt2pit2pdt2pit * const s_data = data;
	const register dim_typ nrow = s_data->a0;
	ityp * mat = s_data->a1;
	ityp * imat = s_data->a2;
	dim_typ * rmod = s_data->a3;
	dim_typ * cmod = s_data->a4;
	ityp * rsort = s_data->a5;
	ityp * csort = s_data->a6;
	
	dim_typ i;
	dim_typ irc;
	dim_typ j;
	dim_typ jrc;
	dim_typ kirc;
	dim_typ kjrc;
	int p;
	/*
	Initialize row and column addresses.
	*/
	#pragma omp parallel for
	for ( i = 1; i <= nrow; ++i )
		rsort[i-1] = csort[i-1] = i;
	/*
	Sort the rows.
	*/
	for ( irc = 1; irc <= nrow; ++irc )
	{
		/*
		Find the next row.
		*/
		jrc = irc;
		p = rmod[irc-1];

		for ( i = irc + 1; i <= nrow; ++i)
		{
			if ( rmod[i-1] < p )
			{
				p = rmod[i-1];
				jrc = i;
			}
		}

		if ( irc != jrc )
		{
			i = rmod[irc-1];
			rmod[irc-1] = rmod[jrc-1];
			rmod[jrc-1] = i;

			i = rsort[irc-1];
			rsort[irc-1] = rsort[jrc-1];
			rsort[jrc-1] = i;
			/*
			Switch the rows.
			*/
			kirc = ( irc - 1 ) * nrow;
			kjrc = ( jrc - 1 ) * nrow;

			for ( j = 1; j <= nrow; ++j )
			{
				i = mat[(kirc+j-1)*nrow];
				mat[(kirc+j-1)*nrow] = mat[(kjrc+j-1)*nrow];
				mat[(kjrc+j-1)*nrow] = i;
			}
		}
	}
	/*
	Sort the columns.
	*/
	for ( irc = 1; irc <= nrow; ++irc )
	{
		/*
		Find the next column.
		*/
		jrc = irc;
		p = cmod[irc-1];

		for ( i = irc + 1; i <= nrow; ++i )
		{
			if ( cmod[i-1] < p )
			{
				p = cmod[i-1];
				jrc = i;
			}
		}

		if ( irc != jrc )
		{
			i = cmod[irc-1];
			cmod[irc-1] = cmod[jrc-1];
			cmod[jrc-1] = i;

			i = csort[irc-1];
			csort[irc-1] = csort[jrc-1];
			csort[jrc-1] = i;
			/*
			Switch the columns.
			*/
			for ( j = 0; j < nrow * nrow; j += nrow )
			{
				i = mat[(irc+j-1)*nrow];
				mat[(irc+j-1)*nrow] = mat[(jrc+j-1)*nrow];
				mat[(jrc+j-1)*nrow] = i;
			}
		}
	}

	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _musort ( void * data)
/******************************************************************************/
/*
  Purpose:
    MUSORT unsorts the inverse matrix rows and columns into the original order.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    08 December 2013
  Author:
    Original FORTRAN77 version by Roger Payne.
    C version by John Burkardt.
  Reference:
    Roger Payne,
    Inversion of matrices with contents subject to modulo arithmetic,
    Applied Statistics,
    Volume 46, Number 2, 1997, pages 295-298.
  Parameters:
    Input/output, int MAT[NROW*NROW].
    On output, the matrix has been "unsorted".
    Input/output, int IMAT[NROW*NROW].
    On output, the matrix has been "unsorted".
    Input/output, int RMOD[NROW], the modulus for values in
    each row.  On output, these have been restored to their original ordering.
    Input/output, int CMOD[NROW], the modulus for values in
    each column.  On output, these have been restored to their original
    ordering.
    Input/output, int RSORT[NROW], the sorted row indices.
    Input/output, int CSORT[NROW], the sorted column indices.
    Input, int NROW, the order of the matrix.
*/
{
	const dt2pit2pdt2pit * const s_data = data;
	const register dim_typ nrow = s_data->a0;
	ityp * mat = s_data->a1;
	ityp * imat = s_data->a2;
	dim_typ * rmod = s_data->a3;
	dim_typ * cmod = s_data->a4;
	ityp * rsort = s_data->a5;
	ityp * csort = s_data->a6;
	
	dim_typ i;
	dim_typ irc;
	dim_typ j;
	dim_typ jrc;
	dim_typ kirc;
	dim_typ kjrc;
	/*
	Sort rows of inverse (= columns of original).
	*/
	for ( irc = 1; irc <= nrow; ++irc )
	{
		/*
		Find next row.
		*/
		if ( csort[irc-1] != irc )
		{
			for ( jrc = irc + 1; jrc <= nrow; ++jrc )
				if ( csort[jrc-1] == irc )
					break;

			i = cmod[irc-1];
			cmod[irc-1] = cmod[jrc-1];
			cmod[jrc-1] = i;

			i = csort[irc-1];
			csort[irc-1] = csort[jrc-1];
			csort[jrc-1] = i;
			/*
			Switch rows.
			*/
			kirc = ( irc - 1 ) * nrow;
			kjrc = ( jrc - 1 ) * nrow;

			for ( j = 1; j <= nrow; ++j )
			{
				i = imat[(kirc+j-1)*nrow];
				imat[(kirc+j-1)*nrow] = imat[(kjrc+j-1)*nrow];
				imat[(kjrc+j-1)*nrow] = i;
			}
		}
	}
	/*
	Sort the columns of the inverse (= rows of original).
	*/
	for ( irc = 1; irc <= nrow; ++irc )
	{
		/*
		Find the next column.
		*/
		if ( rsort[irc-1] != irc )
		{
			for ( jrc = irc + 1; jrc <= nrow; ++jrc )
				if ( rsort[jrc-1] == irc )
					break;

			i = rmod[irc-1];
			rmod[irc-1] = rmod[jrc-1];
			rmod[jrc-1] = i;

			i = rsort[irc-1];
			rsort[irc-1] = rsort[jrc-1];
			rsort[jrc-1] = i;
			/*
			Switch the columns of IMAT.
			*/
			for ( j = 0; j < nrow * nrow; j += nrow )
			{
				i = imat[(irc+j-1)*nrow];
				imat[(irc+j-1)*nrow] = imat[(jrc+j-1)*nrow];
				imat[(jrc+j-1)*nrow] = i;
			}
			/*
			Switch the diagonal elements of MAT (others are zero).
			*/
			kirc = ( irc - 1 ) * nrow + irc;
			kjrc = ( jrc - 1 ) * nrow + jrc;

			i = mat[(kirc-1)*nrow];
			mat[(kirc-1)*nrow] = mat[(kjrc-1)*nrow];
			mat[(kjrc-1)*nrow] = i;
		}
	}
	return NULL;
}

#endif
