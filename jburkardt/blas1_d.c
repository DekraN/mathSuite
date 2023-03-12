#ifndef __DISABLEDEEP_BLAS1D

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _daxpy ( void * data)
/******************************************************************************/
/*
  Purpose:
    DAXPY computes constant times a vector plus a vector.
  Discussion:
    This routine uses unrolled loops for increments equal to one.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 May 2005
  Author:
    FORTRAN77 original version by Jack Dongarra.
    C++ version by John Burkardt.
  Reference:
    Dongarra, Moler, Bunch, Stewart,
    LINPACK User's Guide,
    SIAM, 1979.
    Lawson, Hanson, Kincaid, Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
  Parameters:
    Input, int N, the number of elements in DX and DY.
    Input, double DA, the multiplier of DX.
    Input, double DX[*], the first vector.
    Input, int INCX, the increment between successive entries of DX.
    Input/output, double DY[*], the second vector.
    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
    Input, int INCY, the increment between successive entries of DY.
*/
{
	const dtitpitdtpitdt * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp da = s_data->a1;
	ityp * dx = s_data->a2;
	const register dim_typ incx = s_data->a3;
	ityp * dy = s_data->a4;
	const register dim_typ incy = s_data->a5;
	
    dim_typ i;
    dim_typ ix;
    dim_typ iy;
    dim_typ m;

    if ( n <= 0 || !da )
        return NULL;
    /*
    Code for unequal increments or equal increments
    not equal to 1.
    */
    if ( incx != 1 || incy != 1 )
    {
        ix = 0<=incx ? 0 : ( - n + 1 ) * incx;
        iy = 0<=incy ? 0 : ( - n + 1 ) * incy;

        for ( i = 0; i < n; ++i )
        {
            dy[iy] = dy[iy] + da * dx[ix];
            ix += incx;
            iy += incy;
        }
    }
    /*
    Code for both increments equal to 1.
    */
    else
    {
        m = n % 4;

        for ( i = 0; i < m; i++ )
            dy[i] += da * dx[i];


        for ( i = m; i < n; i += 4 )
        {
            dy[i  ] += da * dx[i  ];
            dy[i+1] += da * dx[i+1];
            dy[i+2] += da * dx[i+2];
            dy[i+3] += da * dx[i+3];
        }
    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ddot ( void * data)
/******************************************************************************/
/*
  Purpose:
    DDOT forms the dot product of two vectors.
  Discussion:
    This routine uses unrolled loops for increments equal to one.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    02 May 2005
  Author:
    FORTRAN77 original version by Jack Dongarra.
    C version by John Burkardt.
  Reference:
    Dongarra, Moler, Bunch, Stewart,
    LINPACK User's Guide,
    SIAM, 1979.
    Lawson, Hanson, Kincaid, Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input, double DX[*], the first vector.
    Input, int INCX, the increment between successive entries in DX.
    Input, double DY[*], the second vector.
    Input, int INCY, the increment between successive entries in DY.
    Output, double DDOT, the sum of the product of the corresponding
    entries of DX and DY.
*/
{
	static ityp result = MAX_VAL;
	
	const _3dt2pit * const s_data = data;
	
	const register dim_typ n = s_data->a0;
	const register dim_typ incx = s_data->a1;
	const register dim_typ incy = s_data->a2;
	ityp * dx = s_data->a3;
	ityp * dy = s_data->a4;
	
	
    ityp dtemp = 0.00;
    dim_typ i;
    dim_typ ix;
    dim_typ iy;
    dim_typ m;


    if ( n <= 0 )
    {
    	result = dtemp;
        return &result;
    }
    /*
    Code for unequal increments or equal increments
    not equal to 1.
    */
    if ( incx != 1 || incy != 1 )
    {
        ix = 0<=incx ? 0 : ( - n + 1 ) * incx;
        iy = 0<=incy ? 0 : ( - n + 1 ) * incy;

        for ( i = 0; i < n; ++i )
        {
            dtemp += dx[ix] * dy[iy];
            ix += incx;
            iy += incy;
        }
    }
    /*
    Code for both increments equal to 1.
    */
    else
    {
        m = n % 5;

        for ( i = 0; i < m; ++i )
            dtemp += dx[i] * dy[i];

        for ( i = m; i < n; i += 5)
        dtemp += dx[i  ] * dy[i  ]+ dx[i+1] * dy[i+1]+ dx[i+2] * dy[i+2]+ dx[i+3] * dy[i+3]+ dx[i+4] * dy[i+4];

    }
    
    result = dtemp;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dnrm2 ( void * data)
/******************************************************************************/
/*
  Purpose:
    DNRM2 returns the euclidean norm of a vector.
  Discussion:
     DNRM2 ( X ) = sqrt ( X' * X )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2007
  Author:
    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
    David Kincaid, Fred Krogh.
    C version by John Burkardt.
  Reference:
    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.
    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, double X[*], the vector whose norm is to be computed.
    Input, int INCX, the increment between successive entries of X.
    Output, double DNRM2, the Euclidean norm of X.
*/
{
	static ityp result = MAX_VAL;
	
	const dtpiti * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	int incx = s_data->a2;
	
    ityp absxi;
    dim_typ i;
    int ix;
    ityp norm;
    ityp scale;
    ityp ssq;
    ityp value;

    if ( n < 1 || incx < 1 )
        norm = 0.00;
    else if ( n == 1 )
        norm = abs ( x[0] );
    else
    {
        scale = 0.00;
        ssq = 1.00;
        ix = 0;

        for ( i = 0; i < n; ++i )
        {
            if ( x[ix] != 0.00 )
            {
                absxi = abs ( x[ix] );
                if ( scale < absxi )
                {
                    ssq = 1.00 + ssq * ( scale / absxi ) * ( scale / absxi );
                    scale = absxi;
                }
                else
                ssq += ( absxi / scale ) * ( absxi / scale );
            }
            ix += incx;
        }

        norm  = scale * sqrt ( ssq );
    }

	result = norm;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _drot ( void * data)
/******************************************************************************/
/*
  Purpose:
    DROT applies a plane rotation.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2007
  Author:
    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
    David Kincaid, Fred Krogh.
    C version by John Burkardt.
  Reference:
    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.
    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input/output, double X[*], one of the vectors to be rotated.
    Input, int INCX, the increment between successive entries of X.
    Input/output, double Y[*], one of the vectors to be rotated.
    Input, int INCY, the increment between successive elements of Y.
    Input, double C, S, parameters (presumably the cosine and
    sine of some angle) that define a plane rotation.
*/
{
	const dtpitipiti2it * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	int incx = s_data->a2;
	ityp * y = s_data->a3;
	int incy = s_data->a4;
	ityp c = s_data->a5;
	ityp s = s_data->a6;
	
    dim_typ i;
    int ix;
    int iy;
    ityp stemp;

    if ( n == 0 );
    else if ( incx == 1 && incy == 1 )
    {
        for ( i = 0; i < n; ++i )
        {
            stemp = c * x[i] + s * y[i];
            y[i]  = c * y[i] - s * x[i];
            x[i]  = stemp;
        }
    }
    else
    {
        ix = 0 + (0>incx)*( - n + 1 ) * incx;
        iy = 0 + (0>incy)*( - n + 1 ) * incy;

        for ( i = 0; i < n; ++i )
        {
            stemp = c * x[ix] + s * y[iy];
            y[iy] = c * y[iy] - s * x[ix];
            x[ix] = stemp;
            ix += incx;
            iy += incy;
        }

    }

    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _drotg ( void * data)
/******************************************************************************/
/*
  Purpose:
    DROTG constructs a Givens plane rotation.
  Discussion:
    Given values A and B, this routine computes
    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
          = sign ( B ) if abs ( A ) <= abs ( B );
    R     = SIGMA * ( A * A + B * B );
    C = A / R if R is not 0
      = 1     if R is 0;
    S = B / R if R is not 0,
        0     if R is 0.
    The computed numbers then satisfy the equation
 (  C  S ) ( A ) = ( R )
 ( -S  C ) ( B ) = ( 0 )
    The routine also computes
    Z = S     if abs ( A ) > abs ( B ),
      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
      = 1     if C is 0.
    The single value Z encodes C and S, and hence the rotation:
    If Z = 1, set C = 0 and S = 1;
    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2007
  Author:
    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
    David Kincaid, Fred Krogh.
    C version by John Burkardt.
  Reference:
    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.
    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
  Parameters:
    Input/output, double *SA, *SB,  On input, SA and SB are the values
    A and B.  On output, SA is overwritten with R, and SB is
    overwritten with Z.
    Output, double *C, *S, the cosine and sine of the Givens rotation.
*/
{
	ityp ** const a_data = data;
	ityp * sa = a_data[0];
	ityp * sb = a_data[1];
	ityp * c = a_data[2];
	ityp * s = a_data[3];
	
    ityp r;
    ityp roe;
    ityp scale;
    ityp z;

    roe = abs ( *sb ) < abs ( *sa ) ? *sa : *sb;
    scale = abs ( *sa ) + abs ( *sb );

    if ( scale == 0.00 )
    {
        *c = 1.0;
        *s = r = 0.00;
    }
    else
    {
        r = scale * sqrt ( ( *sa / scale ) * ( *sa / scale )+ ( *sb / scale ) * ( *sb / scale ) );
        r *= r8_sign ( roe );
        *c = *sa / r;
        *s = *sb / r;
    }

    z = 0.00 < abs ( *c ) && abs ( *c ) <= *s ? 1.00/ *c : *s;
    *sa = r;
    *sb = z;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dscal ( void * data)
/******************************************************************************/
/*
  Purpose:
    DSCAL scales a vector by a constant.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2007
  Author:
    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
    David Kincaid, Fred Krogh.
    C version by John Burkardt.
  Reference:
    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.
    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
  Parameters:
    Input, int N, the number of entries in the vector.
    Input, double SA, the multiplier.
    Input/output, double X[*], the vector to be scaled.
    Input, int INCX, the increment between successive entries of X.
*/
{
	const dtitpiti * const s_data = data;
	const register dim_typ n = s_data->a0;
	const register ityp sa = s_data->a1;
	ityp * x = s_data->a2;
	int incx = s_data->a3;
	
    dim_typ i;
    int ix;
    dim_typ m;

    if ( incx == 1 )
    {
        m = n % 5;

        for ( i = 0; i < m; ++i )
            x[i] *= sa;

        for ( i = m; i < n; i += 5 )
        {
            x[i]   = sa * x[i];
            x[i+1] = sa * x[i+1];
            x[i+2] = sa * x[i+2];
            x[i+3] = sa * x[i+3];
            x[i+4] = sa * x[i+4];
        }
    }
    else
    {
        ix = 0 + (0>incx)*( - n + 1 ) * incx;

        for ( i = 0; i < n; ++i )
        {
            x[ix] = sa * x[ix];
            ix += incx;
        }
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _dswap ( void * data)
/******************************************************************************/
/*
  Purpose:
    DSWAP interchanges two vectors.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 March 2007
  Author:
    Original FORTRAN77 version by Charles Lawson, Richard Hanson,
    David Kincaid, Fred Krogh.
    C version by John Burkardt.
  Reference:
    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.
    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.
  Parameters:
    Input, int N, the number of entries in the vectors.
    Input/output, double X[*], one of the vectors to swap.
    Input, int INCX, the increment between successive entries of X.
    Input/output, double Y[*], one of the vectors to swap.
    Input, int INCY, the increment between successive elements of Y.
*/
{
	const dtpitipiti * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	int incx = s_data->a2;
	ityp * y = s_data->a3;
	int incy = s_data->a4;
	
    dim_typ i;
    int ix;
    int iy;
    dim_typ m;
    ityp temp;

    if ( n == 0);
    else if ( incx == 1 && incy == 1 )
    {
        m = n % 3;

        for ( i = 0; i < m; ++i )
        {
            temp = x[i];
            x[i] = y[i];
            y[i] = temp;
        }

        for ( i = m; i < n; i += 3 )
        {
            temp = x[i];
            x[i] = y[i];
            y[i] = temp;

            temp = x[i+1];
            x[i+1] = y[i+1];
            y[i+1] = temp;

            temp = x[i+2];
            x[i+2] = y[i+2];
            y[i+2] = temp;
        }
    }
    else
    {
        ix = 0 + (0>incx)*( - n + 1 ) * incx;
        iy = 0 + (0>incy)*( - n + 1 ) * incy;

        for ( i = 0; i < n; ++i )
        {
            temp = x[ix];
            x[ix] = y[iy];
            y[iy] = temp;
            ix += incx;
            iy += incy;
        }

    }

    return NULL;
}

#endif
