#ifndef __DISABLEDEEP_ZERORC

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _r8mat_fs ( void * data)
/******************************************************************************/
/*
  Purpose:
    R8MAT_FS factors and solves a system with one right hand side.
  Discussion:
    This routine uses partial pivoting, but no pivot vector is required.
    This routine differs from R8MAT_FSS in two ways:
    * only one right hand side is allowed;
    * the input matrix A is not modified.
    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    21 January 2013
  Author:
    John Burkardt
  Parameters:
    Input, int LDA, the leading dimension of the matrix.
    Input, int N, the order of the matrix.
    N must be positive.
    Input, double A[N*N], the coefficient matrix of the linear system.
    The matrix is stored in an LDAxN array.
    Input/output, double X[N], on input, the right hand side of the
    linear system.  On output, the solution of the linear system.
*/
{
	const _2dt2pit * const s_data = data; 
	const register dim_typ lda = s_data->a0;
	const register dim_typ n = s_data->a1;
	ityp * a = s_data->a2;
	ityp * x = s_data->a3;
	
	ityp *a2;
	dim_typ i;
	int ipiv;
	dim_typ j;
	dim_typ jcol;
	ityp piv;
	ityp t;
	
	a2 = ( ityp * ) malloc ( n * n * sizeof ( ityp ) );
	for ( j = 0; j < n; ++j )
		for ( i = 0; i < n; ++i )
			a2[i+j*n] = a[i+j*lda];
	
	for ( jcol = 1; jcol <= n; ++jcol )
	{
		/*
		Find the maximum element in column I.
		*/
		piv = abs ( a2[jcol-1+(jcol-1)*n] );
		ipiv = jcol;
		for ( i = jcol+1; i <= n; ++i )
		{
			if ( piv < abs ( a2[i-1+(jcol-1)*n] ) )
			{
				piv = abs ( a2[i-1+(jcol-1)*n] );
				ipiv = i;
			}
		}
	
		if ( piv == 0.00 )
			return NULL;
		/*
		Switch rows JCOL and IPIV, and X.
		*/
		if ( jcol != ipiv )
		{
			for ( j = 1; j <= n; ++j )
			{
				t                  = a2[jcol-1+(j-1)*n];
				a2[jcol-1+(j-1)*n] = a2[ipiv-1+(j-1)*n];
				a2[ipiv-1+(j-1)*n] = t;
			}
			t         = x[jcol-1];
			x[jcol-1] = x[ipiv-1];
			x[ipiv-1] = t;
		}
		/*
		Scale the pivot row.
		*/
		t = a2[jcol-1+(jcol-1)*n];
		a2[jcol-1+(jcol-1)*n] = 1.00;
		for ( j = jcol+1; j <= n; ++j )
			a2[jcol-1+(j-1)*n] /= t;
		x[jcol-1] /= t;
		/*
		Use the pivot row to eliminate lower entries in that column.
		*/
		for ( i = jcol+1; i <= n; ++i )
		{
			if ( a2[i-1+(jcol-1)*n] )
			{
				t = - a2[i-1+(jcol-1)*n];
				a2[i-1+(jcol-1)*n] = 0.00;
				for ( j = jcol+1; j <= n; ++j)
					a2[i-1+(j-1)*n] += t * a2[jcol-1+(j-1)*n];
				x[i-1] += t * x[jcol-1];
			}
		}
	}
	/*
	Back solve.
	*/
	for ( jcol = n; 2 <= jcol; --jcol )
		for ( i = 1; i < jcol; ++i )
		x[i-1] -= a2[i-1+(jcol-1)*n] * x[jcol-1];
	
	free ( a2 );
	
	return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _root_rc ( void * data)
/******************************************************************************/
/*
  Purpose:
    ROOT_RC solves a single nonlinear equation using reverse communication.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    22 January 2013
  Author:
    Original FORTRAN77 version by Gaston Gonnet.
    C version by John Burkardt.
  Reference:
    Gaston Gonnet,
    On the Structure of Zero Finders,
    BIT Numerical Mathematics,
    Volume 17, Number 2, June 1977, pages 170-183.
  Parameters:
    Input, double X, an estimate for the root.  On the first
    call, this must be a value chosen by the user.  Thereafter, it may
    be a value chosen by the user, or the value of ROOT returned on the
    previous call to the function.
    Input, double FX, the value of the function at X.
    Output, double *FERR, the smallest value of F encountered.
    Output, double *XERR, the width of the change-in-sign interval,
    if one was encountered.
    Input/output, double Q[9], storage needed by the function.
    Before the first call, the user must set Q(1) to 0.
    Output, double ROOT_RC, an improved estimate for the root.
*/
{
	static ityp result = MAX_VAL;
	
	const _2it3pit * const s_data = data;
	const register ityp x = s_data->a0;
	const register ityp fx = s_data->a1;
	ityp * ferr = s_data->a2;
	ityp * xerr = s_data->a3;
	ityp * q = s_data->a4;
	
    ityp d;
    ityp decr;
    dim_typ i;
    ityp p;
    ityp r;
    ityp u;
    ityp v;
    ityp w;
    ityp xnew;
    ityp z;
    /*
    If we found an exact zero, there is nothing more to do.
    */
    if ( fx == 0.00 )
    {
        *ferr = *xerr = 0.00;
        result = x;
        return &result;
    }

    *ferr = abs ( fx );
    /*
    If this is the first time, initialize, estimate the first root, and exit.
    */
    if ( q[0] )
    {
        q[0] = fx;
        q[1] = x;
        #pragma omp parallel for num_threads(2)
        for ( i = 2; i < 9; ++i )
            q[i] = 0.00;
        xnew = x + fx;
        *xerr = r8_huge;
        
        result = xnew;
        return &result;
    }
    /*
    This is not the first call.
    */
    ++ q[8];
    /*
    Check for too many iterations.
    */
    /*
    Check for a repeated X value.
    */
    if ( 80.00 < q[8] || ( 2.00 <= q[8] && x == q[3] ) || x == q[1] )
    {
    	result = MAX_VAL;
        return &result;
    }
    /*
    Push X -> A -> B -> C
    */
    for ( i = 5; 2 <= i; --i )
        q[i] = q[i-2];
    q[0] = fx;
    q[1] = x;
    /*
    If we have a change-in-sign interval, store the opposite value.
    */
    if ( r8_sign ( q[0] ) != r8_sign ( q[2] ) )
    {
        q[6] = q[2];
        q[7] = q[3];
    }
    /*
    Calculate XERR.
    */
    *xerr = q[6] != 0.00 ? abs ( q[7] - q[1] ) : r8_huge;
    /*
    If more than 30 iterations, and we have change-in-sign interval, bisect.
    */
    if ( 30.00 < q[8] && q[6] != 0.00 )
    {
    	result = q[1] + ( q[7] - q[1] ) / 2.00;
        return &result;
    }

    v = ( q[2] - q[0] ) / ( q[3] - q[1] );
    /*
    If 3 or more points, try Muller.
    */
    if ( q[4] )
    {
        u = ( q[4] - q[2] ) / ( q[5] - q[3] );
        w = q[3] - q[1];
        z = ( q[5] - q[1] ) / w;
        r = ( z + 1.00 ) * v - u;

        if ( r  )
        {
            p = 2.00 * z * q[0] / r;
            d = 2.00 * p / ( w * r ) * ( v - u );
            if ( -1.00 <= d )
            {
                xnew = q[1] - p / ( 1.00 + sqrt ( 1.00 + d ) );
                if ( q[6] == 0.00 || ( q[1] < xnew && xnew < q[7] ) || ( q[7] < xnew && xnew < q[1] ) )
                {
                	result = xnew;
                    return &result;
                }
            }
        }
    }
    /*
    Try the secant step.
    */
    if ( q[0] != q[2] || q[6] == 0.00 )
    {
        if ( q[0] == q[2] )
        {
        	result = MAX_VAL;
            return &result;
        }
        decr = q[0] / v;
        if ( abs ( decr ) * 4.6E+18 < abs ( q[1] ) )
            decr = 1.74E-18 * abs ( q[1] ) * r8_sign ( decr );
        xnew = q[1] - decr;
        if ( q[6] == 0.00 || ( q[1] < xnew && xnew < q[7] ) || ( q[7] < xnew && xnew < q[1] ) )
        {
        	result = xnew; 
            return &result;
        }
    }
    /*
    Apply bisection.
    */
    xnew = q[1] + ( q[7] - q[1] ) / 2.00;
	
	result = xnew;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _roots_rc ( void * data)
/******************************************************************************/
/*
  Purpose:
    ROOTS_RC solves a system of nonlinear equations using reverse communication.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    26 January 2013
  Author:
    Original FORTRAN77 version by Gaston Gonnet.
    C version by John Burkardt.
  Reference:
    Gaston Gonnet,
    On the Structure of Zero Finders,
    BIT Numerical Mathematics,
    Volume 17, Number 2, June 1977, pages 170-183.
  Parameters:
    Input, int N, the number of equations.
    Input, double X[N].  Before the first call, the user should
    set X to an initial guess or estimate for the root.  Thereafter, the input
    value of X should be the output value of XNEW from the previous call.
    Input, double FX[N], the value of the function at XNEW.
    Output, double *FERR, the function error, that is, the sum of
    the absolute values of the most recently computed function vector.
    Output, double XNEW[N], a new point at which a function
    value is requested.
    Workspace, double Q[(2*N+2)*(N+2)].  Before the first call
    for a given problem, the user must set Q(2*N+1,1) to 0.0.
*/
{
	const dt5pit * const s_data = data;
	const register dim_typ n = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	ityp * ferr = s_data->a3;
	ityp * xnew = s_data->a4;
	ityp * q = s_data->a5;
	
    ityp damp;
    dim_typ i, j;
    dim_typ jsma;
    dim_typ jsus;
    dim_typ lda;
    ityp sump;
    ityp t;

    lda = (n<<1) + 2;

    *ferr = 0.00;
    for ( i = 0; i < n; ++i )
        *ferr += abs ( fx[i] );
    /*
    Initialization if Q(2*N+1,1) = 0.0.
    */
    if ( q[(n<<1)+1+0*lda] == 0.00 )
    {
        for ( i = 1; i <= n; ++i)
        {
            for ( j = 1; j <= n + 1; ++j )
            {
                q[i-1+(j-1)*lda] = 0.00;
                q[i+(j-1)*lda] = 0.00;
            }
            q[i-1+(i-1)*lda] = 100.00;
            q[i+n-1+(i-1)*lda] = 1.00;
        }

        for ( j = 1; j <= n; ++j )
        {
            q[(n<<1)+(j-1)*lda] = r8_huge;
            q[(n<<1)+1+(j-1)*lda] = ( ityp ) ( n );
        }

        for ( i = 1; i <= n; ++i )
        {
            q[i+n-1+n*lda] = x[i-1];
            q[i-1+n*lda] = fx[i-1];
        }

        q[(n<<1)+n*lda] = *ferr;
        q[(n<<1)+1+n*lda] = 0.00;
        damp = 0.99;
    }
    else
    {
        jsus = 1;
        for ( i = 2; i <= n + 1; ++i )
        {
            if ( ( ityp ) ( (n<<1) ) <= q[(n<<1)+1+(i-1)*lda] )
                q[(n<<1)+(i-1)*lda] = r8_huge;

            if ( q[(n<<1)+1+(jsus-1)*lda] < ( n + 3 ) / 2 || ( n + 3 ) / 2 <= q[(n<<1)+1+(i-1)*lda] &&q[(n<<1)+(jsus-1)*lda] < q[(n<<1)+(i-1)*lda])
                jsus = i;
        }

        for ( i = 1; i <= n; ++i )
        {
            q[i+n-1+(jsus-1)*lda] = x[i-1];
            q[i-1+(jsus-1)*lda] = fx[i-1];
        }

        q[(n<<1)+(jsus-1)*lda] = *ferr;
        q[(n<<1)+1+(jsus-1)*lda] = 00;
        jsma = 1;
        damp = 0.00;

        for ( j = 1; j <= n + 1; ++j )
        {
            if ( r8_huge / 10.00 < q[(n<<1)+(j-1)*lda] )
                damp = 0.99;
            if ( q[(n<<1)+(j-1)*lda] < q[(n<<1)+(jsma-1)*lda] )
                jsma = j;
        }

        if ( jsma != n + 1 )
        {
            for ( i = 1; i <= (n<<1) + 2; ++i )
            {
                t = q[i-1+(jsma-1)*lda];
                q[i-1+(jsma-1)*lda] = q[i-1+n*lda];
                q[i-1+n*lda] = t;
            }
        }
    }

    for ( i = 1; i <= n; ++i )
        q[i-1+(n+1)*lda] = q[i-1+n*lda];
    /*
    Call the linear equation solver, which should not destroy the matrix
    in Q(1:N,1:N), and should overwrite the solution into Q(1:N,N+2).
    */
    r8mat_fs ( lda, n, q, q+(n+1)*lda );

    sump = 0.00;
    for ( i = 1; i <= n; ++i )
        sump += q[i-1+(n+1)*lda];

    if ( abs ( 1.0 - sump ) <= 1.0E-10 )
        return NULL;

    for ( i = 1; i <= n; ++i )
    {
        xnew[i-1] = q[i+n-1+n*lda];
        for ( j = 1; j <= n; ++j )
            xnew[i-1] -= q[i+n-1+(j-1)*lda] * q[j-1+(n+1)*lda];
        /*
        If system not complete, damp the solution.
        */
        xnew[i-1] /= ( 1.00 - sump ) * ( 1.00 - damp ) + q[i+n-1+n*lda] * damp;
    }

    for ( j = 1; j <= n + 1; ++j )
        ++ q[(n<<1)+1+(j-1)*lda];

    return NULL;
}

#endif
