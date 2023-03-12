#ifndef __DISABLEDEEP_ASA006

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cholesky (void *data)
/******************************************************************************/
/*
  Purpose:
    CHOLESKY computes the Cholesky factorization of a PDS matrix.
  Discussion:
    For a positive definite symmetric matrix A, the Cholesky factor U
    is an upper triangular matrix such that A = U' * U.
    This routine was originally named "CHOL", but that conflicted with
    a built in MATLAB routine name.
    The missing initialization "II = 0" has been added to the code.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    12 February 2008
  Author:
    Original FORTRAN77 version by Michael Healy.
    Modifications by AJ Miller.
    C version by John Burkardt.
  Reference:
    Michael Healy,
    Algorithm AS 6:
    Triangular decomposition of a symmetric matrix,
    Applied Statistics,
    Volume 17, Number 2, 1968, pages 195-197.
  Parameters:
    Input, double A((N*(N+1))/2), a positive definite matrix
    stored by rows in lower triangular form as a one dimensional array,
    in the sequence
    A(1,1),
    A(2,1), A(2,2),
    A(3,1), A(3,2), A(3,3), and so on.
    Input, int N, the order of A.
    Input, int NN, the dimension of the array used to store A,
    which should be at least (N*(N+1))/2.
    Output, double U((N*(N+1))/2), an upper triangular matrix,
    stored by columns, which is the Cholesky factor of A.  The program is
    written in such a way that A and U can share storage.
    Output, int NULLTY, the rank deficiency of A.  If NULLTY is zero,
    the matrix is judged to have full rank.
    Output, int IFAULT, an error indicator.
    0, no error was detected;
    1, if N < 1;
    2, if A is not positive semi-definite.
    3, if NN < (N*(N+1))/2.
  Local Parameters:
    Local, double ETA, should be set equal to the smallest positive
    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
    accuracy being used.
*/
{
	static sel_typ result = UCHAR_MAX;
	
	const pit2ipitpdt * const s_data = data;
	ityp * a = s_data->a0;
	int n = s_data->a1;
	int nn = s_data->a2;
	ityp * u = s_data->a3;
	dim_typ * nullty = s_data->a4;
	
	dim_typ i;
	dim_typ icol;
	dim_typ ii;
	dim_typ irow;
	dim_typ j;
	dim_typ k;
	dim_typ kk;
	dim_typ l;
	dim_typ m;
	ityp w;
	ityp x;

	*nullty = 0;

	if(n < 1)
	{
		result = CHOLESKY_INVALIDMATRIXDIM;
		return &result;
	}

	if ( nn < (n*(n+1))/2)
	{
		result = CHOLESKY_INVALIDSTORAGEDIM;
		return &result;
	}

	j = 1;
	k = ii = 0;
	/*
	Factorize column by column, ICOL = column number.
	*/
	for ( icol = 1; icol <= n; ++icol)
	{
		ii += icol;
		x = 1.0E-09 * 1.0E-09 * a[ii-1];
		l = 0;
		kk = 0;
		/*
		IROW = row number within column ICOL.
		*/
		for ( irow = 1; irow <= icol; ++irow )
		{
			kk = kk + irow;
			k = k + 1;
			w = a[k-1];
			m = j;

			for ( i = 1; i < irow; i++ )
			{
				++ l;
				w -= u[l-1] * u[m-1];
				++ m;
			}

			++ l;

			if (irow == icol)
				break;

			if (u[l-1])
				u[k-1] = w / u[l-1];
			else
			{
				u[k-1] = 0.00;

				if (fabs ( x*a[k-1]) < w*w)
				{
					result = CHOLESKY_INVALIDMATRIX; 
					return &result;
				}
			}
		}
		/*
		End of row, estimate relative accuracy of diagonal element.
		*/
		if ( fabs ( w ) <= fabs ( 1.0E-09 * a[k-1] ) )
		{
			u[k-1] = 0.00;
			++ (*nullty);
		}
		else
		{
			if ( w < 0.00 )
			{
				result = CHOLESKY_INVALIDMATRIX; 
				return &result;
			}
			u[k-1] = sqrt ( w );
		}
		j += icol;
	}

	result = CHOLESKY_SUCCESS; 
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _subchl(void *data)
/******************************************************************************/
/*
  Purpose:
    SUBCHL computes the Cholesky factorization of a (subset of a) PDS matrix.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    11 February 2008
  Author:
    Original FORTRAN77 version by Michael Healy, PR Freeman.
    C version by  John Burkardt.
  Reference:
    PR Freeman,
    Remark AS R44:
    A Remark on AS 6 and AS7: Triangular decomposition of a symmetric matrix
    and Inversion of a positive semi-definite symmetric matrix,
    Applied Statistics,
    Volume 31, Number 3, 1982, pages 336-339.
    Michael Healy,
    Algorithm AS 6:
    Triangular decomposition of a symmetric matrix,
    Applied Statistics,
    Volume 17, Number 2, 1968, pages 195-197.
  Parameters:
    Input, double A((M*(M+1))/2), a positive definite matrix
    stored by rows in lower triangular form as a one dimensional array,
    in the sequence
    A(1,1),
    A(2,1), A(2,2),
    A(3,1), A(3,2), A(3,3), and so on.
    In the simplest case, M, the order of A, is equal to N.
    Input, int B(N), indicates the order in which the
    rows and columns of A are to be used.  In the simplest case,
    B = (1,2,3...,N).
    Input, int N, the order of the matrix, that is,
    the matrix formed by using B to select N rows and columns of A.
    Output, double U((N*(N+1))/2), an upper triangular matrix,
    stored by columns, which is the Cholesky factor of A.  The program is
    written in such a way that A and U can share storage.
    Output, int *NULLTY, the rank deficiency of A.
    If NULLTY is zero, the matrix is judged to have full rank.
    Output, int *IFAULT, an error indicator.
    0, no error was detected;
    1, if N < 1;
    2, if A is not positive semi-definite.
    Input, int NDIM, the dimension of A and U, which might
    be presumed to be (N*(N+1))/2.
    Output, double *DET, the determinant of the matrix.
*/
{
	static sel_typ result = UCHAR_MAX;
	
	const pitdtpdtpitpdtdtpit * const s_data = data;
	ityp * a = s_data->a0;
	const register dim_typ n = s_data->a1;
	const dim_typ * b = s_data->a2;
	ityp * u = s_data->a3; 
	dim_typ * nullty = s_data->a4;
	const register dim_typ ndim = s_data->a5;
	ityp * det = s_data->a6;
	
	dim_typ i;
	dim_typ icol;
	dim_typ ii;
	dim_typ ij;
	dim_typ irow;
	dim_typ j;
	dim_typ jj;
	dim_typ k;
	dim_typ kk;
	dim_typ l;
	dim_typ m;
	ityp w;
	ityp x;

	*nullty = k = 0;
	*det = 1.00;

	if (n <= 0)
	{
		result = CHOLESKY_INVALIDMATRIXDIM;
		return &result;
	}

	j = 1;

	for ( icol = 1; icol <= n; ++icol )
	{
		ij = ( b[icol-1] * ( b[icol-1] - 1 ) ) / 2;
		ii = ij + b[icol-1];
		x = 1.0E-09 * 1.0E-09 * a[ii-1];
		l = 0;

		for ( irow = 1; irow <= icol; irow++ )
		{
			kk = ( b[irow-1] * ( b[irow-1] + 1 ) ) / 2;
			++ k;
			jj = ij + b[irow-1];
			w = a[jj-1];
			m = j;

			for ( i = 1; i <= irow - 1; ++i )
			{
				++ l;
				w -= u[l-1] * u[m-1];
				++ m;
			}

			++ l;

			if(irow == icol)
				break;

			if (u[l-1])
				u[k-1] = w / u[l-1];
			else
			{
				if ( fabs ( x * a[kk-1] ) < w * w )
				{
					result = CHOLESKY_INVALIDMATRIX;
					return &result;
				}
				u[k-1] = 0.00;
			}
		}

		if ( fabs ( 1.0E-09 * a[kk-1] ) <= fabs ( w ) )
		{
			if ( w < 0.00 )
			{
				result = CHOLESKY_INVALIDMATRIXDIM;
				return &result;
			}

			u[k-1] = sqrt ( w );
		}
		else
		{
			u[k-1] = 0.00;
			++ (*nullty);
		}
		j += icol;
		*det *= u[k-1] * u[k-1];
	}

	result = CHOLESKY_SUCCESS;
	return &result;
}

#endif
