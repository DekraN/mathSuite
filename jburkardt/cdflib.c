#ifndef __DISABLEDEEP_CDFLIB

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _algdiv ( void * data)
/******************************************************************************/
//  Purpose:
//    ALGDIV computes ln ( Gamma ( B ) / Gamma ( A + B ) ) when 8 <= B.
//  Discussion:
//    In this algorithm, DEL(X) is the function defined by
//      ln ( Gamma(X) ) = ( X - 0.5 ) * ln ( X ) - X + 0.5 * ln ( 2 * M_PI )
//                      + DEL(X).
//  Parameters:
//    Input, double *A, *B, define the arguments.
//    Output, double ALGDIV, the value of ln(Gamma(B)/Gamma(A+B)).
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * b = a_data[1];
	
    static ityp algdiv;
    static ityp c;
    static ityp c0 =  0.833333333333333e-01;
    static ityp c1 = -0.277777777760991e-02;
    static ityp c2 =  0.793650666825390e-03;
    static ityp c3 = -0.595202931351870e-03;
    static ityp c4 =  0.837308034031215e-03;
    static ityp c5 = -0.165322962780713e-02;
    static ityp d;
    static ityp h;
    static ityp s11;
    static ityp s3;
    static ityp s5;
    static ityp s7;
    static ityp s9;
    static ityp t;
    static ityp T1;
    static ityp u;
    static ityp v;
    static ityp w;
    static ityp x;
    static ityp x2;

    if ( *b <= *a )
    {
        h = *b / *a;
        c = 1.0e0 / ( 1.0e0 + h );
        x = h / ( 1.0e0 + h );
        d = *a + ( *b - 0.5e0 );
    }
    else
    {
        h = *a / *b;
        c = h / ( 1.0e0 + h );
        x = 1.0e0 / ( 1.0e0 + h );
        d = *b + ( *a - 0.5e0 );
    }
    //  SET SN = (1 - X**N)/(1 - X)
    x2 = x * x;
    s3 = 1.0e0 + ( x + x2 );
    s5 = 1.0e0 + ( x + x2 * s3 );
    s7 = 1.0e0 + ( x + x2 * s5 );
    s9 = 1.0e0 + ( x + x2 * s7 );
    s11 = 1.0e0 + ( x + x2 * s9 );
    //  SET W = DEL(B) - DEL(A + B)
    t = pow ( 1.0e0 / *b, 2.0 );

    w = ((((( c5 * s11  * t + c4 * s9 ) * t + c3 * s7 ) * t + c2 * s5 ) * t + c1 * s3 ) * t + c0)*(c/(*b));
    //  Combine the results.
    T1 = *a / *b;
    u = d * alnrel ( &T1 );
    v = *a * ( log ( *b ) - 1.0e0 );

	result = w - (v<u?(v+u):(u+v));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _alnrel ( void * data)
/******************************************************************************/
//  Purpose:
//    ALNREL evaluates the function ln ( 1 + A ).
//  Modified:
//    17 November 2006
//  Reference:
//    Armido DiDinato, Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//  Parameters:
//    Input, double *A, the argument.
//    Output, double ALNREL, the value of ln ( 1 + A ).
{
	static ityp result = MAX_VAL;
	
	ityp * a = (ityp *) data; 
	
    ityp alnrel;
    static ityp p1 = -0.129418923021993e+01;
    static ityp p2 =  0.405303492862024e+00;
    static ityp p3 = -0.178874546012214e-01;
    static ityp q1 = -0.162752256355323e+01;
    static ityp q2 =  0.747811014037616e+00;
    static ityp q3 = -0.845104217945565e-01;
    ityp t;
    ityp t2;
    ityp w;
    ityp x;

    if ( fabs ( *a ) <= 0.375e0 )
    {
        t = *a / ( *a + 2.0e0 );
        t2 = t * t;
        w = (((p3*t2+p2)*t2+p1)*t2+1.0e0)
        / (((q3*t2+q2)*t2+q1)*t2+1.0e0);
        alnrel = 2.0e0 * t * w;
    }
    else
    {
        x = 1.0e0 + *a;
        alnrel = log ( x );
    }
    
    result = alnrel;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _apser ( void * data)
/******************************************************************************/
//  Purpose:
//    APSER computes the incomplete beta ratio I(SUB(1-X))(B,A).
//  Discussion:
//    APSER is used only for cases where
//      A <= MIN ( EPS, EPS * B ),
//      B * X <= 1, and
//      X <= 0.5.
//  Parameters:
//    Input, double *A, *B, *X, the parameters of teh
//    incomplete beta ratio.
//    Input, double *EPS, a tolerance.
//    Output, double APSER, the computed value of the
//    incomplete beta ratio.
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * b = a_data[1];
	ityp * x = a_data[2];
	ityp * eps = a_data[3];
	
    static ityp g = 0.577215664901533e0;
    static ityp apser,aj,bx,c,j,s,t,tol;

    bx = *b**x;
    t = *x-bx;
    if(*b**eps > 2.e-2) goto S10;
    c = log(*x)+psi(b)+g+t;
    goto S20;
    S10:
        c = log(bx)+g+t;
    S20:
        tol = 5.0e0**eps*fabs(c);
        j = 1.0e0;
        s = 0.0e0;
    S30:
        j += 1.0e0;
        t *= (*x-bx/j);
        aj = t/j;
        s += aj;
        if(fabs(aj) > tol) goto S30;

	result = -(*a*(c+s));
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _bcorr ( void * data)
/******************************************************************************/
//  Purpose:
//    BCORR evaluates DEL(A0) + DEL(B0) - DEL(A0 + B0).
//  Discussion:
//    The function DEL(A) is a remainder term that is used in the expression:
//      ln ( Gamma ( A ) ) = ( A - 0.5 ) * ln ( A )
//        - A + 0.5 * ln ( 2 * M_PI ) + DEL ( A ),
//    or, in other words, DEL ( A ) is defined as:
//      DEL ( A ) = ln ( Gamma ( A ) ) - ( A - 0.5 ) * ln ( A )
//        + A + 0.5 * ln ( 2 * M_PI ).
//  Parameters:
//    Input, double *A0, *B0, the arguments.
//    It is assumed that 8 <= A0 and 8 <= B0.
//    Output, double *BCORR, the value of the function.
//
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a0 = a_data[0];
	ityp * b0 = a_data[1];
	
    static ityp c0 =  0.833333333333333e-01;
    static ityp c1 = -0.277777777760991e-02;
    static ityp c2 =  0.793650666825390e-03;
    static ityp c3 = -0.595202931351870e-03;
    static ityp c4 =  0.837308034031215e-03;
    static ityp c5 = -0.165322962780713e-02;
    static ityp bcorr,a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2;

    a = fifdmin1 ( *a0, *b0 );
    b = fifdmax1 ( *a0, *b0 );
    h = a / b;
    c = h / ( 1.0e0 + h );
    x = 1.0e0 / ( 1.0e0 + h );
    x2 = x * x;
    //
    //  SET SN = (1 - X**N)/(1 - X)
    //
    s3 = 1.0e0 + ( x + x2 );
    s5 = 1.0e0 + ( x + x2 * s3 );
    s7 = 1.0e0 + ( x + x2 * s5 );
    s9 = 1.0e0 + ( x + x2 * s7 );
    s11 = 1.0e0 + ( x + x2 * s9 );
    //
    //  SET W = DEL(B) - DEL(A + B)
    //
    t = pow ( 1.0e0 / b, 2.0 );

    w = ((((( c5 * s11  * t + c4 * s9 ) * t + c3 * s7 ) * t + c2 * s5 ) * t + c1 * s3 ) * t + c0)*(c/b);
    //
    //  COMPUTE  DEL(A) + W
    //
    t = pow ( 1.0e0 / a, 2.0 );

	result = ((((( c5 * t + c4 ) * t + c3 ) * t + c2 ) * t + c1 ) * t + c0 ) / a + w;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _beta ( void * data)
/******************************************************************************/
//
//  Purpose:
//    BETA evaluates the beta function.
//  Modified:
//    03 December 1999
//  Author:
//    John Burkardt
//  Parameters:
//    Input, double A, B, the arguments of the beta function.
//    Output, double BETA, the value of the beta function.
//
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	result = ( exp ( beta_log ( &a_data[0], &a_data[1] ) ) );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _beta_asym ( void * data)
/******************************************************************************/
//  Purpose:
//    BETA_ASYM computes an asymptotic expansion for IX(A,B), for large A and B.
//  Parameters:
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that both A and B
//    are greater than or equal to 15.
//    Input, double *LAMBDA, the value of ( A + B ) * Y - B.
//    It is assumed that 0 <= LAMBDA.
//    Input, double *EPS, the tolerance.
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * b = a_data[1];
	ityp * lambda = a_data[2];
	ityp * eps = a_data[3];
	
    static ityp e0 = 1.12837916709551e0;
    static ityp e1 = .353553390593274e0;
    static dim_typ num = 20;
    //  NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
    //            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
    //            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
    //     E0 = 2/SQRT(M_PI)
    //     E1 = 2**(-3/2)
    static dim_typ K3 = 1;
    static ityp value;
    static ityp bsum,dsum,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,t,t0,t1,u,w,w0,z,z0,z2,zn,znm1;
    static dim_typ i,im1,imj,j,m,mm1,mmj,n,np1;
    static ityp a0[21],b0[21],c[21],d[21],T1,T2;

    value = 0.0e0;
    if(*a >= *b) goto S10;
    h = *a/ *b;
    r0 = 1.0e0/(1.0e0+h);
    r1 = (*b-*a)/ *b;
    w0 = 1.0e0/sqrt(*a*(1.0e0+h));
    goto S20;
    S10:
        h = *b/ *a;
        r0 = 1.0e0/(1.0e0+h);
        r1 = (*b-*a)/ *a;
        w0 = 1.0e0/sqrt(*b*(1.0e0+h));
    S20:
        T1 = -(*lambda/ *a);
        T2 = *lambda/ *b;
        f = *a*rlog1(&T1)+*b*rlog1(&T2);
        t = exp(-f);
        if(t == 0.0e0)
		{
			result = value;
			return &result;
		}
        z0 = sqrt(f);
        z = 0.5e0*(z0/e1);
        z2 = f+f;
        a0[0] = 2.0e0/3.0e0*r1;
        c[0] = -(0.5e0*a0[0]);
        d[0] = -c[0];
        j0 = 0.5e0/e0 * error_fc ( &K3, &z0 );
        j1 = e1;
        sum = j0+d[0]*w0*j1;
        s = 1.0e0;
        h2 = h*h;
        hn = 1.0e0;
        w = w0;
        znm1 = z;
        zn = z2;
        for ( n = 2; n <= num; n += 2 )
        {
            hn = h2*hn;
            a0[n-1] = 2.0e0*r0*(1.0e0+h*hn)/((ityp)n+2.0e0);
            np1 = n+1;
            s += hn;
            a0[np1-1] = 2.0e0*r1*s/((ityp)n+3.0e0);
            for ( i = n; i <= np1; ++i )
            {
                r = -(0.5e0*((ityp)i+1.0e0));
                b0[0] = r*a0[0];
                for ( m = 2; m <= i; m++ )
                {
                    bsum = 0.0e0;
                    mm1 = m-1;
                    for ( j = 1; j <= mm1; ++j)
                    {
                        mmj = m-j;
                        bsum += (((ityp)j*r-(ityp)mmj)*a0[j-1]*b0[mmj-1]);
                    }
                    b0[m-1] = r*a0[m-1]+bsum/(ityp)m;
                }
                c[i-1] = b0[i-1]/((double)i+1.0e0);
                dsum = 0.0e0;
                im1 = i-1;
                for ( j = 1; j <= im1; ++j )
                {
                    imj = i-j;
                    dsum += (d[imj-1]*c[j-1]);
                }
                d[i-1] = -(dsum+c[i-1]);
            }
            j0 = e1*znm1+((double)n-1.0e0)*j0;
            j1 = e1*zn+(double)n*j1;
            znm1 = z2*znm1;
            zn = z2*zn;
            w = w0*w;
            t0 = d[n-1]*w*j0;
            w = w0*w;
            t1 = d[np1-1]*w*j1;
            sum += (t0+t1);
            if(fabs(t0)+fabs(t1) <= *eps*sum) goto S80;
        }
        S80:
            u = exp(-bcorr(a,b));
            value = e0*t*u*sum;
            
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _beta_frac ( void * data)
/******************************************************************************/
//  Purpose:
//    BETA_FRAC evaluates a continued fraction expansion for IX(A,B).
//  Parameters:
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that both A and
//    B are greater than 1.
//    Input, double *X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//    Input, double *LAMBDA, the value of ( A + B ) * Y - B.
//    Input, double *EPS, a tolerance.
//    Output, double BETA_FRAC, the value of the continued
//    fraction approximation for IX(A,B).
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * b = a_data[1];
	ityp * x = a_data[2];
	ityp * y = a_data[3];
	ityp * lambda = a_data[4];
	ityp * eps = a_data[5];
	
    static ityp bfrac,alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,t,w,yp1;

    bfrac = beta_rcomp ( a, b, x, y );

    if ( bfrac == 0.0e0 )
    {
    	result = bfrac;
        return &result;
    }

    c = 1.0e0+*lambda;
    c0 = *b/ *a;
    c1 = 1.0e0+1.0e0/ *a;
    yp1 = *y+1.0e0;
    n = 0.0e0;
    p = 1.0e0;
    s = *a+1.0e0;
    an = 0.0e0;
    bn = anp1 = 1.0e0;
    bnp1 = c/c1;
    r = c1/c;
    //  CONTINUED FRACTION CALCULATION
    S10:
        n += 1.0e0;
        t = n/ *a;
        w = n*(*b-n)**x;
        e = *a/s;
        alpha = p*(p+c0)*e*e*(w**x);
        e = (1.0e0+t)/(c1+t+t);
        beta = n+w/s+e*(c+n*yp1);
        p = 1.0e0+t;
        s += 2.0e0;
        //  UPDATE AN, BN, ANP1, AND BNP1
        t = alpha*an+beta*anp1;
        an = anp1;
        anp1 = t;
        t = alpha*bn+beta*bnp1;
        bn = bnp1;
        bnp1 = t;
        r0 = r;
        r = anp1/bnp1;

    if ( fabs(r-r0) <= (*eps) * r )
        goto S20;
    //  RESCALE AN, BN, ANP1, AND BNP1

        an /= bnp1;
        bn /= bnp1;
        anp1 = r;
        bnp1 = 1.0e0;
        goto S10;
    //  TERMINATION
    S20:
        bfrac *=  r;
    
    result = bfrac;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _beta_grat ( void * data)
/******************************************************************************/
//  Purpose:
//    BETA_GRAT evaluates an asymptotic expansion for IX(A,B).
//  Parameters:
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.  It is assumed that 15 <= A
//    and B <= 1, and that B is less than A.
//    Input, double *X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//    Input/output, double *W, a quantity to which the
//    result of the computation is to be added on output.
//    Input, double *EPS, a tolerance.
//    Output, int *IERR, an error flag, which is 0 if no error
//    was detected.
{
	const pdt6pit * const s_data = data;
	
	dim_typ * ierr = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp * w = s_data->a5;
	ityp * eps = s_data->a6;
	
	
    static ityp bm1,bp2n,cn,coef,dj,j,l,lnx,n2,nu,p,q,r,s,sum,t,t2,u,v,z;
    static dim_typ i,n,nm1;
    static ityp c[30],d[30],T1;

    bm1 = *b-0.5e0-0.5e0;
    nu = *a+0.5e0*bm1;
    if(*y > 0.375e0) goto S10;
    T1 = -*y;
    lnx = alnrel(&T1);
    goto S20;
    S10:
        lnx = log(*x);
    S20:
        z = -(nu*lnx);
        if(*b*z == 0.0e0) goto S70;
        //  COMPUTATION OF THE EXPANSION
        //  SET R = EXP(-Z)*Z**B/GAMMA(B)
        r = *b*(1.0e0+gam1(b))*exp(*b*log(z));
        r *= (exp(*a*lnx)*exp(0.5e0*bm1*lnx));
        u = algdiv(b,a)+*b*log(nu);
        u = r*exp(-u);
        if(u == 0.0e0) goto S70;
        gamma_rat1 ( b, &z, &r, &p, &q, eps );
        v = 0.25e0*pow(1.0e0/nu,2.0);
        t2 = 0.25e0*lnx*lnx;
        l = *w/u;
        j = q/r;
        sum = j;
        t = cn = 1.0e0;
        n2 = 0.0e0;
        for ( n = 1; n <= 30; ++n)
        {
            bp2n = *b+n2;
            j = (bp2n*(bp2n+1.0e0)*j+(z+bp2n+1.0e0)*t)*v;
            n2 = n2 + 2.0e0;
            t *= t2;
            cn /= (n2*(n2+1.0e0));
            c[n-1] = cn;
            s = 0.0e0;
            if(n == 1) goto S40;
            nm1 = n-1;
            coef = *b-(ityp)n;
            for ( i = 1; i <= nm1; ++i )
            {
                s += (coef*c[i-1]*d[n-i-1]);
                coef += *b;
            }
            S40:
                d[n-1] = bm1*cn+s/(double)n;
                dj = d[n-1]*j;
                sum = sum + dj;
                if(sum <= 0.0e0) goto S70;
                if(fabs(dj) <= *eps*(sum+l)) goto S60;
        }
        S60:
            //  ADD THE RESULTS TO W
            *ierr = 0;
            *w += (u*sum);
            return NULL;
        S70:
        //  THE EXPANSION CANNOT BE COMPUTED
        *ierr = 1;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _beta_inc ( void * data)
/******************************************************************************/
//  Purpose:
//    BETA_INC evaluates the incomplete beta function IX(A,B).
//  Author:
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//  Parameters:
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.
//    Input, double *X, *Y.  X is the argument of the
//    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
//    Output, double *W, *W1, the values of IX(A,B) and
//    1-IX(A,B).
//    Output, int *IERR, the error flag.
//    0, no error was detected.
//    1, A or B is negative;
//    2, A = B = 0;
//    3, X < 0 or 1 < X;
//    4, Y < 0 or 1 < Y;
//    5, X + Y /= 1;
//    6, X = A = 0;
//    7, Y = B = 0.
{
	const pdt6pit * const s_data = data;
	
	dim_typ * ierr = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp * w = s_data->a5;
	ityp * w1 = s_data->a6;
	
	
    static dim_typ K1 = 1;
    static ityp a0,b0,eps,lambda,t,x0,y0,z;
    static dim_typ ierr1,ind,n;
    static ityp T2,T3,T4,T5;
    //  EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
    //  NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
    eps = dpmpar ( &K1 );
    *w = *w1 = 0.0e0;
    if(*a < 0.0e0 || *b < 0.0e0) goto S270;
    if(*a == 0.0e0 && *b == 0.0e0) goto S280;
    if(*x < 0.0e0 || *x > 1.0e0) goto S290;
    if(*y < 0.0e0 || *y > 1.0e0) goto S300;
    z = *x+*y-0.5e0-0.5e0;
    if(fabs(z) > 3.0e0*eps) goto S310;
    *ierr = 0;
    if(*x == 0.0e0) goto S210;
    if(*y == 0.0e0) goto S230;
    if(*a == 0.0e0) goto S240;
    if(*b == 0.0e0) goto S220;
    eps = fifdmax1(eps,1.e-15);
    if(fifdmax1(*a,*b) < 1.e-3*eps) goto S260;
    ind = 0;
    a0 = *a;
    b0 = *b;
    x0 = *x;
    y0 = *y;
    if(fifdmin1(a0,b0) > 1.0e0) goto S40;
    //  PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
    if(*x <= 0.5e0) goto S10;
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;
    S10:
        if(b0 < fifdmin1(eps,eps*a0)) goto S90;
        if(a0 < fifdmin1(eps,eps*b0) && b0*x0 <= 1.0e0) goto S100;
        if(fifdmax1(a0,b0) > 1.0e0) goto S20;
        if(a0 >= fifdmin1(0.2e0,b0)) goto S110;
        if(pow(x0,a0) <= 0.9e0) goto S110;
        if(x0 >= 0.3e0) goto S120;
        n = 20;
        goto S140;
    S20:
        if(b0 <= 1.0e0) goto S110;
        if(x0 >= 0.3e0) goto S120;
        if(x0 >= 0.1e0) goto S30;
        if(pow(x0*b0,a0) <= 0.7e0) goto S110;
    S30:
        if(b0 > 15.0e0) goto S150;
        n = 20;
        goto S140;
    S40:
        //  PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
        if(*a > *b) goto S50;
        lambda = *a-(*a+*b)**x;
        goto S60;
    S50:
        lambda = (*a+*b)**y-*b;
    S60:
        if(lambda >= 0.0e0) goto S70;
        ind = 1;
        a0 = *b;
        b0 = *a;
        x0 = *y;
        y0 = *x;
        lambda = fabs(lambda);
        S70:
        if(b0 < 40.0e0 && b0*x0 <= 0.7e0) goto S110;
        if(b0 < 40.0e0) goto S160;
        if(a0 > b0) goto S80;
        if(a0 <= 100.0e0) goto S130;
        if(lambda > 0.03e0*a0) goto S130;
        goto S200;
    S80:
        if(b0 <= 100.0e0) goto S130;
        if(lambda > 0.03e0*b0) goto S130;
        goto S200;
    S90:
        //  EVALUATION OF THE APPROPRIATE ALGORITHM
        *w = fpser(&a0,&b0,&x0,&eps);
        *w1 = 0.5e0+(0.5e0-*w);
        goto S250;
    S100:
        *w1 = apser(&a0,&b0,&x0,&eps);
        *w = 0.5e0+(0.5e0-*w1);
        goto S250;
    S110:
        *w = beta_pser(&a0,&b0,&x0,&eps);
        *w1 = 0.5e0+(0.5e0-*w);
        goto S250;
    S120:
        *w1 = beta_pser(&b0,&a0,&y0,&eps);
        *w = 0.5e0+(0.5e0-*w1);
        goto S250;
    S130:
        T2 = 15.0e0*eps;
        *w = beta_frac ( &a0,&b0,&x0,&y0,&lambda,&T2 );
        *w1 = 0.5e0+(0.5e0-*w);
        goto S250;
    S140:
        *w1 = beta_up ( &b0, &a0, &y0, &x0, &n, &eps );
        b0 += (ityp)n;
        S150:
        T3 = 15.0e0*eps;
        beta_grat (&b0,&a0,&y0,&x0,w1,&T3,&ierr1);
        *w = 0.5e0+(0.5e0-*w1);
        goto S250;
        S160:
        n = ( dim_typ ) b0;
        b0 -= (ityp)n;
        if(b0 != 0.0e0) goto S170;
        -- n;
        b0 = 1.0e0;
    S170:
        *w = beta_up ( &b0, &a0, &y0, &x0, &n, &eps );
        if(x0 > 0.7e0) goto S180;
        *w += beta_pser(&a0,&b0,&x0,&eps);
        *w1 = 0.5e0+(0.5e0-*w);
        goto S250;
    S180:
        if(a0 > 15.0e0) goto S190;
        n = 20;
        *w += beta_up ( &a0, &b0, &x0, &y0, &n, &eps );
        a0 += (ityp)n;
    S190:
        T4 = 15.0e0*eps;
        beta_grat ( &a0, &b0, &x0, &y0, w, &T4, &ierr1 );
        *w1 = 0.5e0+(0.5e0-*w);
        goto S250;
    S200:
        T5 = 100.0e0*eps;
        *w = beta_asym ( &a0, &b0, &lambda, &T5 );
        *w1 = 0.5e0+(0.5e0-*w);
        goto S250;
    S210:
        //  TERMINATION OF THE PROCEDURE
        if(*a == 0.0e0) goto S320;
        S220:
        *w = 0.0e0;
        *w1 = 1.0e0;
        return NULL;
    S230:
        if(*b == 0.0e0) goto S330;
        S240:
        *w = 1.0e0;
        *w1 = 0.0e0;
        return NULL;
    S250:
        if(!ind) return NULL;
        t = *w;
        *w = *w1;
        *w1 = t;
        return NULL;
    S260:
        //  PROCEDURE FOR A AND B .LT. 1.E-3*EPS
        *w = *b/(*a+*b);
        *w1 = *a/(*a+*b);
        return NULL;
    S270:
        //  ERROR RETURN
        *ierr = 1;
        return NULL;
    S280:
        *ierr = 2;
        return NULL;
    S290:
        *ierr = 3;
        return NULL;
    S300:
        *ierr = 4;
        return NULL;
    S310:
        *ierr = 5;
        return NULL;
    S320:
        *ierr = 6;
        return NULL;
    S330:
        *ierr = 7;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _beta_log ( void * data)
/******************************************************************************/
//  Purpose:
//    BETA_LOG evaluates the logarithm of the beta function.
//  Reference:
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//  Parameters:
//    Input, double *A0, *B0, the parameters of the function.
//    A0 and B0 should be nonnegative.
//    Output, double *BETA_LOG, the value of the logarithm
//    of the Beta function.
//
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a0 = a_data[0];
	ityp * b0 = a_data[1];
	
    static ityp e = .918938533204673e0;
    static ityp value,a,b,c,h,u,v,w,z;
    static dim_typ i,n;
    static ityp T1;

    a = fifdmin1(*a0,*b0);
    b = fifdmax1(*a0,*b0);
    if(a >= 8.0e0) goto S100;
    if(a >= 1.0e0) goto S20;
    //  PROCEDURE WHEN A .LT. 1
    if(b >= 8.0e0) goto S10;
    T1 = a+b;
    value = gamma_log ( &a )+( gamma_log ( &b )- gamma_log ( &T1 ));
    result = value;
    return &result;
    S10:
        value = gamma_log ( &a )+algdiv(&a,&b);
        result = value;
    	return &result;
    S20:
        //  PROCEDURE WHEN 1 .LE. A .LT. 8
        if(a > 2.0e0) goto S40;
        if(b > 2.0e0) goto S30;
        value = gamma_log ( &a )+ gamma_log ( &b )-gsumln(&a,&b);
        result = value;
    	return &result;
    S30:
        w = 0.0e0;
        if(b < 8.0e0) goto S60;
        value = gamma_log ( &a )+algdiv(&a,&b);
        result = value;
    	return &result;
    S40:
        //  REDUCTION OF A WHEN B .LE. 1000
        if(b > 1000.0e0) goto S80;
        n = ( int ) ( a - 1.0e0 );
        w = 1.0e0;
        for ( i = 1; i <= n; ++i )
        {
            a -= 1.0e0;
            h = a/b;
            w *= (h/(1.0e0+h));
        }
        w = log(w);
        if(b < 8.0e0) goto S60;
        value = w+ gamma_log ( &a )+algdiv(&a,&b);
        result = value;
    	return &result;
    S60:
        //  REDUCTION OF B WHEN B .LT. 8
        n = ( dim_typ ) ( b - 1.0e0 );
        z = 1.0e0;
        for ( i = 1; i <= n; ++i )
        {
            b -= 1.0e0;
            z *= (b/(a+b));
        }
        value = w+log(z)+( gamma_log ( &a )+( gamma_log ( &b )-gsumln(&a,&b)));
        result = value;
    	return &result;
    S80:
        //  REDUCTION OF A WHEN B .GT. 1000
        n = ( dim_typ ) ( a - 1.0e0 );
        w = 1.0e0;
        for ( i = 1; i <= n; ++i )
        {
            a -= 1.0e0;
            w *= (a/(1.0e0+a/b));
        }
        value = log(w)-(ityp)n*log(b)+( gamma_log ( &a )+algdiv(&a,&b));
        result = value;
    	return &result;
    S100:
        //  PROCEDURE WHEN A .GE. 8
        w = bcorr(&a,&b);
        h = a/b;
        c = h/(1.0e0+h);
        u = -((a-0.5e0)*log(c));
        v = b*alnrel(&h);
        if(u <= v) goto S110;
        value = -(0.5e0*log(b))+e+w-v-u;
        result = value;
    	return &result;
    S110:
        value = -(0.5e0*log(b))+e+w-u-v;
        
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _beta_pser ( void * data)
/******************************************************************************/
//
//  Purpose:
//    BETA_PSER uses a power series expansion to evaluate IX(A,B)(X).
//  Discussion:
//    BETA_PSER is used when B <= 1 or B*X <= 0.7.
//  Parameters:
//    Input, double *A, *B, the parameters.
//    Input, double *X, the point where the function
//    is to be evaluated.
//    Input, double *EPS, the tolerance.
//    Output, double BETA_PSER, the approximate value of IX(A,B)(X).
//
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * b = a_data[1];
	ityp * x = a_data[2];
	ityp * eps = a_data[3];
	
    static ityp bpser,a0,apb,b0,c,n,sum,t,tol,u,w,z;
    static dim_typ i,m;

    bpser = 0.0e0;
    if(*x == 0.0e0)
	{
		result = bpser;
		return &result;
	}
    //  COMPUTE THE FACTOR X**A/(A*BETA(A,B))
    a0 = fifdmin1(*a,*b);
    if(a0 < 1.0e0) goto S10;
    z = *a*log(*x)-beta_log(a,b);
    bpser = exp(z)/ *a;
    goto S100;
    S10:
        b0 = fifdmax1(*a,*b);
        if(b0 >= 8.0e0) goto S90;
        if(b0 > 1.0e0) goto S40;
        //  PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1
        bpser = pow(*x,*a);
        if(bpser == 0.0e0)
		{
			result = bpser;
			return &result;
		}
        apb = *a+*b;
        if(apb > 1.0e0) goto S20;
        z = 1.0e0+gam1(&apb);
        goto S30;
    S20:
        u = *a+*b-1.e0;
        z = (1.0e0+gam1(&u))/apb;
    S30:
        c = (1.0e0+gam1(a))*(1.0e0+gam1(b))/z;
        bpser *= (c*(*b/apb));
        goto S100;
    S40:
        //  PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8
        u = gamma_ln1 ( &a0 );
        m = ( dim_typ ) ( b0 - 1.0e0 );
        if(m < 1) goto S60;
        c = 1.0e0;
        for ( i = 1; i <= m; ++i )
        {
            b0 -= 1.0e0;
            c *= (b0/(a0+b0));
        }
        u = log(c)+u;
        S60:
        z = *a*log(*x)-u;
        b0 -= 1.0e0;
        apb = a0+b0;
        if(apb > 1.0e0) goto S70;
        t = 1.0e0+gam1(&apb);
        goto S80;
    S70:
        u = a0+b0-1.e0;
        t = (1.0e0+gam1(&u))/apb;
    S80:
        bpser = exp(z)*(a0/ *a)*(1.0e0+gam1(&b0))/t;
        goto S100;
    S90:
        //  PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8
        u = gamma_ln1 ( &a0 ) + algdiv ( &a0, &b0 );
        z = *a*log(*x)-u;
        bpser = a0/ *a*exp(z);
    S100:
        if(bpser == 0.0e0 || *a <= 0.1e0**eps)
        {
			result = bpser;
			return &result;
		}
        //  COMPUTE THE SERIES
        sum = n = 0.0e0;
        c = 1.0e0;
        tol = *eps/ *a;
    S110:
        n += 1.0e0;
        c *= ((0.5e0+(0.5e0-*b/n))**x);
        w = c/(*a+n);
        sum += w;
        if(fabs(w) > tol) goto S110;
        bpser *= (1.0e0+*a*sum);
        
    result = bpser;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _beta_rcomp ( void * data)
/******************************************************************************/
//
//  Purpose:
//    BETA_RCOMP evaluates X**A * Y**B / Beta(A,B).
//  Parameters:
//    Input, double *A, *B, the parameters of the Beta function.
//    A and B should be nonnegative.
//    Input, double *X, *Y, define the numerator of the fraction.
//    Output, double BETA_RCOMP, the value of X**A * Y**B / Beta(A,B).
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * b = a_data[1];
	ityp * x = a_data[2];
	ityp * y = a_data[3];
	
    static ityp Const = .398942280401433e0;
    static ityp brcomp,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
    static dim_typ i,n;
    //  CONST = 1/SQRT(2*M_PI)
    static ityp T1,T2;

    brcomp = 0.0e0;
    if(*x == 0.0e0 || *y == 0.0e0)
	{
		result = brcomp;
		return &result;
	}
    a0 = fifdmin1(*a,*b);
    if(a0 >= 8.0e0) goto S130;
    if(*x > 0.375e0) goto S10;
    lnx = log(*x);
    T1 = -*x;
    lny = alnrel(&T1);
    goto S30;
    S10:
        if(*y > 0.375e0) goto S20;
        T2 = -*y;
        lnx = alnrel(&T2);
        lny = log(*y);
        goto S30;
    S20:
        lnx = log(*x);
        lny = log(*y);
    S30:
        z = *a*lnx+*b*lny;
        if(a0 < 1.0e0) goto S40;
        z -= beta_log(a,b);
        brcomp = exp(z);
        result = brcomp;
		return &result;
    S40:
        //  PROCEDURE FOR A .LT. 1 OR B .LT. 1
        b0 = fifdmax1(*a,*b);
        if(b0 >= 8.0e0) goto S120;
        if(b0 > 1.0e0) goto S70;
        //  ALGORITHM FOR B0 .LE. 1
        brcomp = exp(z);
        if(brcomp == 0.0e0)
        {
        	result = brcomp;
			return &result;
        }
        apb = *a+*b;
        if(apb > 1.0e0) goto S50;
        z = 1.0e0+gam1(&apb);
        goto S60;
    S50:
        u = *a+*b-1.e0;
        z = (1.0e0+gam1(&u))/apb;
    S60:
        c = (1.0e0+gam1(a))*(1.0e0+gam1(b))/z;
        brcomp = brcomp*(a0*c)/(1.0e0+a0/b0);
        result = brcomp;
		return &result;
    S70:
        //  ALGORITHM FOR 1 .LT. B0 .LT. 8
        u = gamma_ln1 ( &a0 );
        n = ( dim_typ ) ( b0 - 1.0e0 );
        if(n < 1) goto S90;
        c = 1.0e0;
        for ( i = 1; i <= n; ++i )
        {
            b0 -= 1.0e0;
            c *= (b0/(a0+b0));
        }
        u = log(c)+u;
        S90:
        z -= u;
        b0 -= 1.0e0;
        apb = a0+b0;
        if(apb > 1.0e0) goto S100;
        t = 1.0e0+gam1(&apb);
        goto S110;
    S100:
        u = a0+b0-1.e0;
        t = (1.0e0+gam1(&u))/apb;
    S110:
        brcomp = a0*exp(z)*(1.0e0+gam1(&b0))/t;
        result = brcomp;
		return &result;
    S120:
        //  ALGORITHM FOR B0 .GE. 8
        u = gamma_ln1 ( &a0 ) + algdiv ( &a0, &b0 );
        brcomp = a0*exp(z-u);
        result = brcomp;
		return &result;
    S130:
        //  PROCEDURE FOR A .GE. 8 AND B .GE. 8
        if(*a > *b) goto S140;
        h = *a/ *b;
        x0 = h/(1.0e0+h);
        y0 = 1.0e0/(1.0e0+h);
        lambda = *a-(*a+*b)**x;
        goto S150;
    S140:
        h = *b/ *a;
        x0 = 1.0e0/(1.0e0+h);
        y0 = h/(1.0e0+h);
        lambda = (*a+*b)**y-*b;
    S150:
        e = -(lambda/ *a);
        if(fabs(e) > 0.6e0) goto S160;
        u = rlog1(&e);
        goto S170;
    S160:
        u = e-log(*x/x0);
    S170:
        e = lambda/ *b;
        if(fabs(e) > 0.6e0) goto S180;
        v = rlog1(&e);
        goto S190;
    S180:
        v = e-log(*y/y0);
    S190:
        z = exp(-(*a*u+*b*v));
        brcomp = Const*sqrt(*b*x0)*z*exp(-bcorr(a,b));
        
    result = brcomp;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _beta_rcomp1 ( void * data)
/******************************************************************************/
//  Purpose:
//    BETA_RCOMP1 evaluates exp(MU) * X**A * Y**B / Beta(A,B).
//  Parameters:
//    Input, int MU, ?
//    Input, double A, B, the parameters of the Beta function.
//    A and B should be nonnegative.
//    Input, double X, Y, ?
//    Output, double BETA_RCOMP1, the value of
//    exp(MU) * X**A * Y**B / Beta(A,B).
{
	static ityp result = MAX_VAL;
	
	const pdt4pit * const s_data = data;
	dim_typ * mu = s_data->a0;
	ityp * a = s_data->a1;
	ityp * b = s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	
    static ityp Const = .398942280401433e0;
    static ityp brcmp1,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
    static dim_typ i,n;
    //     CONST = 1/SQRT(2*M_PI)
    static ityp T1,T2,T3,T4;

    a0 = fifdmin1(*a,*b);
    if(a0 >= 8.0e0) goto S130;
    if(*x > 0.375e0) goto S10;
    lnx = log(*x);
    T1 = -*x;
    lny = alnrel(&T1);
    goto S30;
    S10:
        if(*y > 0.375e0) goto S20;
        T2 = -*y;
        lnx = alnrel(&T2);
        lny = log(*y);
        goto S30;
    S20:
        lnx = log(*x);
        lny = log(*y);
    S30:
        z = *a*lnx+*b*lny;
        if(a0 < 1.0e0) goto S40;
        z -= beta_log(a,b);
        brcmp1 = esum(mu,&z);
        result = brcmp1;
        return &result;
    S40:
        //   PROCEDURE FOR A .LT. 1 OR B .LT. 1
        b0 = fifdmax1(*a,*b);
        if(b0 >= 8.0e0) goto S120;
        if(b0 > 1.0e0) goto S70;
        //  ALGORITHM FOR B0 .LE. 1
        brcmp1 = esum(mu,&z);
        if(brcmp1 == 0.0e0)
        {
        	result = brcmp1;
        	return &result;
        }
        apb = *a+*b;
        if(apb > 1.0e0) goto S50;
        z = 1.0e0+gam1(&apb);
        goto S60;
    S50:
        u = *a+*b-1.e0;
        z = (1.0e0+gam1(&u))/apb;
    S60:
        c = (1.0e0+gam1(a))*(1.0e0+gam1(b))/z;
        brcmp1 = brcmp1*(a0*c)/(1.0e0+a0/b0);
        result = brcmp1;
        return &result;
    S70:
        //  ALGORITHM FOR 1 .LT. B0 .LT. 8
        u = gamma_ln1 ( &a0 );
        n = ( int ) ( b0 - 1.0e0 );
        if(n < 1) goto S90;
        c = 1.0e0;
        for ( i = 1; i <= n; ++i)
        {
            b0 -= 1.0e0;
            c *= (b0/(a0+b0));
        }
        u = log(c)+u;
    S90:
        z -= u;
        b0 -= 1.0e0;
        apb = a0+b0;
        if(apb > 1.0e0) goto S100;
        t = 1.0e0+gam1(&apb);
        goto S110;
    S100:
        u = a0+b0-1.e0;
        t = (1.0e0+gam1(&u))/apb;
    S110:
        brcmp1 = a0*esum(mu,&z)*(1.0e0+gam1(&b0))/t;
        result = brcmp1;
        return &result;
    S120:
        //  ALGORITHM FOR B0 .GE. 8
        u = gamma_ln1 ( &a0 ) + algdiv ( &a0, &b0 );
        T3 = z-u;
        brcmp1 = a0*esum(mu,&T3);
        result = brcmp1;
        return &result;
    S130:
        //    PROCEDURE FOR A .GE. 8 AND B .GE. 8
        if(*a > *b) goto S140;
        h = *a/ *b;
        x0 = h/(1.0e0+h);
        y0 = 1.0e0/(1.0e0+h);
        lambda = *a-(*a+*b)**x;
        goto S150;
    S140:
        h = *b/ *a;
        x0 = 1.0e0/(1.0e0+h);
        y0 = h/(1.0e0+h);
        lambda = (*a+*b)**y-*b;
    S150:
        e = -(lambda/ *a);
        if(fabs(e) > 0.6e0) goto S160;
        u = rlog1(&e);
        goto S170;
    S160:
        u = e-log(*x/x0);
    S170:
        e = lambda/ *b;
        if(fabs(e) > 0.6e0) goto S180;
        v = rlog1(&e);
        goto S190;
    S180:
        v = e-log(*y/y0);
    S190:
        T4 = -(*a*u+*b*v);
        z = esum(mu,&T4);
        brcmp1 = Const*sqrt(*b*x0)*z*exp(-bcorr(a,b));
        
    result = brcmp1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _beta_up ( void * data)
/******************************************************************************/
//  Purpose:
//    BETA_UP evaluates IX(A,B) - IX(A+N,B) where N is a positive integer.
//  Parameters:
//    Input, double *A, *B, the parameters of the function.
//    A and B should be nonnegative.
//    Input, double *X, *Y, ?
//    Input, int *N, ?
//    Input, double *EPS, the tolerance.
//    Output, double BETA_UP, the value of IX(A,B) - IX(A+N,B).
{
	static ityp result = MAX_VAL;
	
	const _4pitpdtpit * const s_data = data;
	ityp * a = s_data->a0;
	ityp * b = s_data->a1;
	ityp * x = s_data->a2;
	ityp * y = s_data->a3;
	dim_typ * n = s_data->a4;
	ityp * eps = s_data->a5;
	
    static dim_typ K1 = 1;
    static dim_typ K2 = 0;
    static ityp bup,ap1,apb,d,l,r,t,w;
    static dim_typ i,k,kp1,mu,nm1;
    //  OBTAIN THE SCALING FACTOR EXP(-MU) AND
    //  EXP(MU)*(X**A*Y**B/BETA(A,B))/A
    apb = *a+*b;
    ap1 = *a+1.0e0;
    mu = 0;
    d = 1.0e0;
    if(*n == 1 || *a < 1.0e0) goto S10;
    if(apb < 1.1e0*ap1) goto S10;
    mu = ( dim_typ ) fabs ( exparg(&K1) );
    k = ( dim_typ ) exparg ( &K2 );
    if(k < mu) mu = k;
    t = mu;
    d = exp(-t);
    S10:
        bup = beta_rcomp1 ( &mu, a, b, x, y ) / *a;
        if(*n == 1 || bup == 0.0e0)
		{
			result = bup;
			return &result;
		}
        nm1 = *n-1;
        w = d;
        //  LET K BE THE INDEX OF THE MAXIMUM TERM
        k = 0;
        if(*b <= 1.0e0) goto S50;
        if(*y > 1.e-4) goto S20;
        k = nm1;
        goto S30;
    S20:
        r = (*b-1.0e0)**x/ *y-*a;
        if(r < 1.0e0) goto S50;
        t = ( double ) nm1;
        k = nm1;
        if ( r < t ) k = ( dim_typ ) r;
    S30:
        //          ADD THE INCREASING TERMS OF THE SERIES
        for ( i = 1; i <= k; ++i )
        {
            l = i-1;
            d = (apb+l)/(ap1+l)**x*d;
            w += d;
        }
        if(k == nm1) goto S70;
    S50:
        //          ADD THE REMAINING TERMS OF THE SERIES
        kp1 = k+1;
        for ( i = kp1; i <= nm1; ++i )
        {
            l = i-1;
            d = (apb+l)/(ap1+l)**x*d;
            w += d;
            if(d <= *eps*w) goto S70;
        }
    S70:
        //
        //  TERMINATE THE PROCEDURE
        //
        bup *= w;
        
    result = bup;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cdfbet ( void * data)
/******************************************************************************/
//  Purpose:
//    CDFBET evaluates the CDF of the Beta Distribution.
//  Discussion:
//    This routine calculates any one parameter of the beta distribution
//    given the others.
//    The value P of the cumulative distribution function is calculated
//    directly by code associated with the reference.
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//    The beta density is proportional to t^(A-1) * (1-t)^(B-1).
//  Modified:
//    09 June 2004
//  Reference:
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//  Parameters:
//    Input, int *WHICH, indicates which of the next four argument
//    values is to be calculated from the others.
//    1: Calculate P and Q from X, Y, A and B;
//    2: Calculate X and Y from P, Q, A and B;
//    3: Calculate A from P, Q, X, Y and B;
//    4: Calculate B from P, Q, X, Y and A.
//    Input/output, double *P, the integral from 0 to X of the
//    chi-square distribution.  Input range: [0, 1].
//    Input/output, double *Q, equals 1-P.  Input range: [0, 1].
//    Input/output, double *X, the upper limit of integration
//    of the beta density.  If it is an input value, it should lie in
//    the range [0,1].  If it is an output value, it will be searched for
//    in the range [0,1].
//    Input/output, double *Y, equal to 1-X.  If it is an input
//    value, it should lie in the range [0,1].  If it is an output value,
//    it will be searched for in the range [0,1].
//    Input/output, double *A, the first parameter of the beta
//    density.  If it is an input value, it should lie in the range
// (0, +infinity).  If it is an output value, it will be searched
//    for in the range [1D-300,1D300].
//    Input/output, double *B, the second parameter of the beta
//    density.  If it is an input value, it should lie in the range
// (0, +infinity).  If it is an output value, it will be searched
//    for in the range [1D-300,1D300].
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +4, if X + Y /= 1.
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
{
	const pdt6pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * x = s_data->a3;
	ityp * y = s_data->a4;
	ityp * a = s_data->a5;
	ityp * b = s_data->a6;
	short * status = s_data->a7;
	ityp * bound = s_data->a8; 
	
    # define tol (1.0e-8)
    # define atol (1.0e-50)
    # undef zero
    # define zero (1.0e-300)
    # define inf 1.0e300
    # define one 1.0e0

    static dim_typ K1 = 1;
    static ityp K2 = 0.0e0;
    static ityp K3 = 1.0e0;
    static ityp K8 = 0.5e0;
    static ityp K9 = 5.0e0;
    static ityp fx,xhi,xlo,cum,ccum,xy,pq;
    static dim_typ qporq;
    static unsigned long qleft, qhi; 
    static ityp T4,T5,T6,T7,T10,T11,T12,T13,T14,T15;

    *status = 0;
    *bound = 0.00;
    //     Check arguments
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 4.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
        if(!(*p < 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = 1.0e0;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60:
    if(*which == 1) goto S110;
    //     Q
    if(!(*q < 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q < 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
    S80:
        *bound = 1.0e0;
    S90:
        *status = -3;
        return NULL;
    S110:
    S100:
        if(*which == 2) goto S150;
        //     X
        if(!(*x < 0.0e0 || *x > 1.0e0)) goto S140;
        if(!(*x < 0.0e0)) goto S120;
        *bound = 0.0e0;
        goto S130;
    S120:
        *bound = 1.0e0;
    S130:
        *status = -4;
        return NULL;
    S150:
    S140:
        if(*which == 2) goto S190;
        //     Y
        if(!(*y < 0.0e0 || *y > 1.0e0)) goto S180;
        if(!(*y < 0.0e0)) goto S160;
        *bound = 0.0e0;
        goto S170;
    S160:
        *bound = 1.0e0;
    S170:
        *status = -5;
        return NULL;
    S190:
    S180:
        if(*which == 3) goto S210;
        //     A
        if(!(*a <= 0.0e0)) goto S200;
        *bound = 0.0e0;
        *status = -6;
        return NULL;
    S210:
    S200:
        if(*which == 4) goto S230;
        //     B
        if(!(*b <= 0.0e0)) goto S220;
        *bound = 0.0e0;
        *status = -7;
        return NULL;
    S230:
    S220:
        if(*which == 1) goto S270;
        //     P + Q
        pq = *p+*q;
        if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S260;
        if(!(pq < 0.0e0)) goto S240;
        *bound = 0.0e0;
        goto S250;
    S240:
        *bound = 1.0e0;
    S250:
        *status = 3;
        return NULL;
    S270:
    S260:
        if(*which == 2) goto S310;
        //     X + Y
        xy = *x+*y;
        if(!(fabs(xy-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S300;
        if(!(xy < 0.0e0)) goto S280;
        *bound = 0.0e0;
        goto S290;
    S280:
        *bound = 1.0e0;
    S290:
        *status = 4;
        return NULL;
    S310:
    S300:
        if(!(*which == 1)) qporq = *p <= *q;
        //     Select the minimum of P or Q
        //     Calculate ANSWERS
        if(1 == *which)
        {
            //     Calculating P and Q
            cumbet(x,y,a,b,p,q);
            *status = 0;
        }
        else if(2 == *which)
        {
            //     Calculating X and Y
            T4 = atol;
            T5 = tol;
            dstzr(&K2,&K3,&T4,&T5);
            if(!qporq) goto S340;
            *status = 0;
            dzror(status,x,&fx,&xlo,&xhi,&qleft,&qhi);
            *y = one-*x;
            S320:
            if(!(*status == 1)) goto S330;
            cumbet(x,y,a,b,&cum,&ccum);
            fx = cum-*p;
            dzror(status,x,&fx,&xlo,&xhi,&qleft,&qhi);
            *y = one-*x;
            goto S320;
            S330:
                goto S370;
            S340:
                *status = 0;
                dzror(status,y,&fx,&xlo,&xhi,&qleft,&qhi);
                *x = one-*y;
            S350:
                if(!(*status == 1)) goto S360;
                cumbet(x,y,a,b,&cum,&ccum);
                fx = ccum-*q;
                dzror(status,y,&fx,&xlo,&xhi,&qleft,&qhi);
                *x = one-*y;
                goto S350;
            S370:
            S360:
                if(!(*status == -1)) goto S400;
                if(!qleft) goto S380;
                *status = 1;
                *bound = 0.0e0;
                goto S390;
            S380:
                *status = 2;
                *bound = 1.0e0;
            S400:
            S390:
            ;
        }
        else if(3 == *which)
        {
        //     Computing A
        *a = 5.0e0;
        T6 = zero;
        T7 = inf;
        T10 = atol;
        T11 = tol;
        dstinv(&T6,&T7,&K8,&K8,&K9,&T10,&T11);
        *status = 0;
        dinvr(status,a,&fx,&qleft,&qhi);
    S410:
        if(!(*status == 1)) goto S440;
        cumbet(x,y,a,b,&cum,&ccum);
        if(!qporq) goto S420;
        fx = cum-*p;
        goto S430;
    S420:
        fx = ccum-*q;
    S430:
        dinvr(status,a,&fx,&qleft,&qhi);
        goto S410;
    S440:
        if(!(*status == -1)) goto S470;
        if(!qleft) goto S450;
        *status = 1;
        *bound = zero;
        goto S460;
    S450:
        *status = 2;
        *bound = inf;
    S470:
    S460:
    ;
    }
    else if(4 == *which)
    {

        //     Computing B
        *b = 5.0e0;
        T12 = zero;
        T13 = inf;
        T14 = atol;
        T15 = tol;
        dstinv(&T12,&T13,&K8,&K8,&K9,&T14,&T15);
        *status = 0;
        dinvr(status,b,&fx,&qleft,&qhi);
    S480:
        if(!(*status == 1)) goto S510;
        cumbet(x,y,a,b,&cum,&ccum);
        if(!qporq) goto S490;
        fx = cum-*p;
        goto S500;
    S490:
        fx = ccum-*q;
    S500:
        dinvr(status,b,&fx,&qleft,&qhi);
        goto S480;
    S510:
        if(!(*status == -1)) goto S540;
        if(!qleft) goto S520;
        *status = 1;
        *bound = zero;
        goto S530;
    S520:
        *status = 2;
        *bound = inf;
    S530:
    ;
    }
    S540:
        return NULL;
    # undef tol
    # undef atol
    # undef zero
    # undef inf
    # undef one
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cdfbin ( void * data)
/******************************************************************************/
//  Purpose:
//    CDFBIN evaluates the CDF of the Binomial distribution.
//  Discussion:
//    This routine calculates any one parameter of the binomial distribution
//    given the others.
//    The value P of the cumulative distribution function is calculated
//    directly.
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//    P is the probablility of S or fewer successes in XN binomial trials,
//    each trial having an individual probability of success of PR.
//  Modified:
//    09 June 2004
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.24.
//  Parameters:
//    Input, int *WHICH, indicates which of argument values is to
//    be calculated from the others.
//    1: Calculate P and Q from S, XN, PR and OMPR;
//    2: Calculate S from P, Q, XN, PR and OMPR;
//    3: Calculate XN from P, Q, S, PR and OMPR;
//    4: Calculate PR and OMPR from P, Q, S and XN.
//    Input/output, double *P, the cumulation, from 0 to S,
//    of the binomial distribution.  If P is an input value, it should
//    lie in the range [0,1].
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//    Input/output, double *S, the number of successes observed.
//    Whether this is an input or output value, it should lie in the
//    range [0,XN].
//    Input/output, double *XN, the number of binomial trials.
//    If this is an input value it should lie in the range: (0, +infinity).
//    If it is an output value it will be searched for in the
//    range [1.0D-300, 1.0D+300].
//    Input/output, double *PR, the probability of success in each
//    binomial trial.  Whether this is an input or output value, it should
//    lie in the range: [0,1].
//    Input/output, double *OMPR, equal to 1-PR.  Whether this is an
//    input or output value, it should lie in the range [0,1].  Also, it should
//    be the case that PR + OMPR = 1.
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +4, if PR + OMPR /= 1.
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
{
	const pdt6pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * s = s_data->a3;
	ityp * xn = s_data->a4;
	ityp * pr = s_data->a5;
	ityp * ompr = s_data->a6;
	short * status = s_data->a7;
	ityp * bound = s_data->a8; 
	
    # define atol (1.0e-50)
    # define tol (1.0e-8)
    # define zero (1.0e-300)
    # define inf 1.0e300
    # define one 1.0e0

    static dim_typ K1 = 1;
    static ityp K2 = 0.0e0;
    static ityp K3 = 0.5e0;
    static ityp K4 = 5.0e0;
    static ityp K11 = 1.0e0;
    static ityp fx,xhi,xlo,cum,ccum,pq,prompr;
    static dim_typ qporq;
    static unsigned long qleft, qhi;
    static ityp T5,T6,T7,T8,T9,T10,T12,T13;

    *status = 0;
    *bound = 0.00;
    //     Check arguments
    if(!(*which < 1 && *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 4.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
        if(!(*p < 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = 1.0e0;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60:
        if(*which == 1) goto S110;
        //     Q
        if(!(*q < 0.0e0 || *q > 1.0e0)) goto S100;
        if(!(*q < 0.0e0)) goto S80;
        *bound = 0.0e0;
        goto S90;
    S80:
        *bound = 1.0e0;
    S90:
        *status = -3;
        return NULL;
    S110:
    S100:
        if(*which == 3) goto S130;
        //     XN
        if(!(*xn <= 0.0e0)) goto S120;
        *bound = 0.0e0;
        *status = -5;
        return NULL;
    S130:
    S120:
        if(*which == 2) goto S170;
        //     S

        if(!(*s < 0.0e0 || *which != 3 && *s > *xn)) goto S160;
        if(!(*s < 0.0e0)) goto S140;
        *bound = 0.0e0;
        goto S150;
    S140:
        *bound = *xn;
    S150:
        *status = -4;
        return NULL;
    S170:
    S160:
        if(*which == 4) goto S210;
        //     PR
        if(!(*pr < 0.0e0 || *pr > 1.0e0)) goto S200;
        if(!(*pr < 0.0e0)) goto S180;
        *bound = 0.0e0;
        goto S190;
    S180:
        *bound = 1.0e0;
    S190:
        *status = -6;
        return NULL;
        S210:
    S200:
        if(*which == 4) goto S250;
        //     OMPR
        if(!(*ompr < 0.0e0 || *ompr > 1.0e0)) goto S240;
        if(!(*ompr < 0.0e0)) goto S220;
        *bound = 0.0e0;
        goto S230;
    S220:
        *bound = 1.0e0;
    S230:
        *status = -7;
        return NULL;
    S250:
    S240:
        if(*which == 1) goto S290;
        //     P + Q
        pq = *p+*q;
        if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S280;
        if(!(pq < 0.0e0)) goto S260;
        *bound = 0.0e0;
        goto S270;
    S260:
        *bound = 1.0e0;
    S270:
        *status = 3;
        return NULL;
    S290:
    S280:
        if(*which == 4) goto S330;
        //     PR + OMPR
        prompr = *pr+*ompr;
        if(!(fabs(prompr-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S320;
        if(!(prompr < 0.0e0)) goto S300;
        *bound = 0.0e0;
        goto S310;
    S300:
        *bound = 1.0e0;
    S310:
        *status = 4;
        return NULL;
    S330:
    S320:
        if(!(*which == 1)) qporq = *p <= *q;
        //     Select the minimum of P or Q
        //     Calculate ANSWERS
        if(1 == *which)
        {
            //     Calculating P
            cumbin(s,xn,pr,ompr,p,q);
            *status = 0;
        }
        else if(2 == *which)
        {
            //     Calculating S
            *s = 5.0e0;
            T5 = atol;
            T6 = tol;
            dstinv(&K2,xn,&K3,&K3,&K4,&T5,&T6);
            *status = 0;
            dinvr(status,s,&fx,&qleft,&qhi);
        S340:
            if(!(*status == 1)) goto S370;
            cumbin(s,xn,pr,ompr,&cum,&ccum);
            if(!qporq) goto S350;
            fx = cum-*p;
            goto S360;
        S350:
            fx = ccum-*q;
        S360:
            dinvr(status,s,&fx,&qleft,&qhi);
            goto S340;
        S370:
            if(!(*status == -1)) goto S400;
            if(!qleft) goto S380;
            *status = 1;
            *bound = 0.0e0;
            goto S390;
        S380:
            *status = 2;
            *bound = *xn;
        S400:
        S390:
        ;
        }
        else if(3 == *which)
        {
            //     Calculating XN
            *xn = 5.0e0;
            T7 = zero;
            T8 = inf;
            T9 = atol;
            T10 = tol;
            dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
            *status = 0;
            dinvr(status,xn,&fx,&qleft,&qhi);
            S410:
                if(!(*status == 1)) goto S440;
                cumbin(s,xn,pr,ompr,&cum,&ccum);
                if(!qporq) goto S420;
                fx = cum-*p;
                goto S430;
            S420:
                fx = ccum-*q;
                S430:
                dinvr(status,xn,&fx,&qleft,&qhi);
                goto S410;
            S440:
                if(!(*status == -1)) goto S470;
                if(!qleft) goto S450;
                *status = 1;
                *bound = zero;
                goto S460;
            S450:
                *status = 2;
                *bound = inf;
            S470:
            S460:
            ;
        }
        else if(4 == *which)
        {
            //     Calculating PR and OMPR
            T12 = atol;
            T13 = tol;
            dstzr(&K2,&K11,&T12,&T13);
            if(!qporq) goto S500;
            *status = 0;
            dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
            *ompr = one-*pr;
            S480:
                if(!(*status == 1)) goto S490;
                cumbin(s,xn,pr,ompr,&cum,&ccum);
                fx = cum-*p;
                dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
                *ompr = one-*pr;
                goto S480;
            S490:
            goto S530;
            S500:
                *status = 0;
                dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
                *pr = one-*ompr;
            S510:
                if(!(*status == 1)) goto S520;
                cumbin(s,xn,pr,ompr,&cum,&ccum);
                fx = ccum-*q;
                dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
                *pr = one-*ompr;
                goto S510;
            S530:
            S520:
                if(!(*status == -1)) goto S560;
                if(!qleft) goto S540;
                *status = 1;
                *bound = 0.0e0;
                goto S550;
            S540:
                *status = 2;
                *bound = 1.0e0;
            S550:
            ;
        }
        S560:
            return NULL;
    # undef atol
    # undef tol
    # undef zero
    # undef inf
    # undef one
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cdfchi ( void *data)
/******************************************************************************/
//
//  Purpose:
//    CDFCHI evaluates the CDF of the chi square distribution.
//  Discussion:
//    This routine calculates any one parameter of the chi square distribution
//    given the others.
//    The value P of the cumulative distribution function is calculated
//    directly.
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//    The CDF of the chi square distribution can be evaluated
//    within Mathematica by commands such as:
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF [ ChiSquareDistribution [ DF ], X ]
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.4.19.
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//  Parameters:
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X and DF;
//    2: Calculate X from P, Q and DF;
//    3: Calculate DF from P, Q and X.
//    Input/output, double *P, the integral from 0 to X of
//    the chi-square distribution.  If this is an input value, it should
//    lie in the range [0,1].
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//    Input/output, double *X, the upper limit of integration
//    of the chi-square distribution.  If this is an input
//    value, it should lie in the range: [0, +infinity).  If it is an output
//    value, it will be searched for in the range: [0,1.0D+300].
//    Input/output, double *DF, the degrees of freedom of the
//    chi-square distribution.  If this is an input value, it should lie
//    in the range: (0, +infinity).  If it is an output value, it will be
//    searched for in the range: [ 1.0D-300, 1.0D+300].
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +10, an error was returned from CUMGAM.
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
{
	const pdt4pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * x = s_data->a3;
	ityp * df = s_data->a4;
	short * status = s_data->a5;
	ityp * bound = s_data->a6; 
	
    # define tol (1.0e-8)
    # define atol (1.0e-50)
    # define zero (1.0e-300)
    # define inf 1.0e300

    static dim_typ K1 = 1;
    static ityp K2 = 0.0e0;
    static ityp K4 = 0.5e0;
    static ityp K5 = 5.0e0;
    static ityp fx,cum,ccum,pq,porq;
    static dim_typ qporq;
    static unsigned long qleft, qhi;
    static ityp T3,T6,T7,T8,T9,T10,T11;

    *status = 0;
    *bound = 0.00;
    //     Check arguments
    if(!(*which < 1 || *which > 3)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 3.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
        if(!(*p < 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = 1.0e0;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60: 
        if(*which == 1) goto S110;
        //     Q
        if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
        if(!(*q <= 0.0e0)) goto S80;
        *bound = 0.0e0;
        goto S90;
    S80:
        *bound = 1.0e0;
        S90:
        *status = -3;
        return NULL;
    S110:
    S100:
        if(*which == 2) goto S130;
        //     X
        if(!(*x < 0.0e0)) goto S120;
        *bound = 0.0e0;
        *status = -4;
        return NULL;
    S130:
    S120:
        if(*which == 3) goto S150;
        //     DF
        if(!(*df <= 0.0e0)) goto S140;
        *bound = 0.0e0;
        *status = -5;
        return NULL;
    S150:
    S140:
        if(*which == 1) goto S190;
        //     P + Q
        pq = *p+*q;
        if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S180;
        if(!(pq < 0.0e0)) goto S160;
        *bound = 0.0e0;
        goto S170;
    S160:
        *bound = 1.0e0;
    S170:
        *status = 3;
        return NULL;
    S190:
    S180:
    if(*which == 1) goto S220;
    //     Select the minimum of P or Q
    qporq = *p <= *q;
    if(!qporq) goto S200;
    porq = *p;
    goto S210;
    S200:
        porq = *q;
    S220:
    S210:
        //     Calculate ANSWERS
        if(1 == *which)
        {

            //     Calculating P and Q
            *status = 0;
            cumchi(x,df,p,q);
            if(porq > 1.5e0)
            {
                *status = 10;
                return NULL;
            }
        }
        else if(2 == *which)
        {
            //     Calculating X
            *x = 5.0e0;
            T3 = inf;
            T6 = atol;
            T7 = tol;
            dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
            *status = 0;
            dinvr(status,x,&fx,&qleft,&qhi);
            S230:
                if(!(*status == 1)) goto S270;
                cumchi(x,df,&cum,&ccum);
                if(!qporq) goto S240;
                fx = cum-*p;
                goto S250;
            S240:
                fx = ccum-*q;
            S250:
                if(!(fx+porq > 1.5e0)) goto S260;
                *status = 10;
                return NULL;
            S260:
                dinvr(status,x,&fx,&qleft,&qhi);
                goto S230;
            S270:
                if(!(*status == -1)) goto S300;
                if(!qleft) goto S280;
                *status = 1;
                *bound = 0.0e0;
                goto S290;
            S280:
                *status = 2;
                *bound = inf;
            S300:
            S290:
                ;
        }
        else if(3 == *which)
        {
            //  Calculating DF
            *df = 5.0e0;
            T8 = zero;
            T9 = inf;
            T10 = atol;
            T11 = tol;
            dstinv(&T8,&T9,&K4,&K4,&K5,&T10,&T11);
            *status = 0;
            dinvr(status,df,&fx,&qleft,&qhi);
            S310:
                if(!(*status == 1)) goto S350;
                cumchi(x,df,&cum,&ccum);
                if(!qporq) goto S320;
                fx = cum-*p;
                goto S330;
            S320:
                fx = ccum-*q;
            S330:
                if(!(fx+porq > 1.5e0)) goto S340;
                *status = 10;
                return NULL;
            S340:
                dinvr(status,df,&fx,&qleft,&qhi);
                goto S310;
            S350:
                if(!(*status == -1)) goto S380;
                if(!qleft) goto S360;
                *status = 1;
                *bound = zero;
                goto S370;
            S360:
                *status = 2;
                *bound = inf;
            S370:
                ;
        }
        S380:
            return NULL;
    # undef tol
    # undef atol
    # undef zero
    # undef inf
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cdfchn ( void * data)
/******************************************************************************/
//  Purpose:
//    CDFCHN evaluates the CDF of the Noncentral Chi-Square.
//  Discussion:
//    This routine calculates any one parameter of the noncentral chi-square
//    distribution given values for the others.
//    The value P of the cumulative distribution function is calculated
//    directly.
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//    The computation time required for this routine is proportional
//    to the noncentrality parameter (PNONC).  Very large values of
//    this parameter can consume immense computer resources.  This is
//    why the search range is bounded by 10,000.
//    The CDF of the noncentral chi square distribution can be evaluated
//    within Mathematica by commands such as:
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.25.
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//  Parameters: 
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X, DF and PNONC;
//    2: Calculate X from P, DF and PNONC;
//    3: Calculate DF from P, X and PNONC;
//    4: Calculate PNONC from P, X and DF.
//    Input/output, double *P, the integral from 0 to X of
//    the noncentral chi-square distribution.  If this is an input
//    value, it should lie in the range: [0, 1.0-1.0D-16).
//    Input/output, double *Q, is generally not used by this
//    subroutine and is only included for similarity with other routines.
//    However, if P is to be computed, then a value will also be computed
//    for Q.
//    Input, double *X, the upper limit of integration of the
//    noncentral chi-square distribution.  If this is an input value, it
//    should lie in the range: [0, +infinity).  If it is an output value,
//    it will be sought in the range: [0,1.0D+300].
//    Input/output, double *DF, the number of degrees of freedom
//    of the noncentral chi-square distribution.  If this is an input value,
//    it should lie in the range: (0, +infinity).  If it is an output value,
//    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
//    Input/output, double *PNONC, the noncentrality parameter of
//    the noncentral chi-square distribution.  If this is an input value, it
//    should lie in the range: [0, +infinity).  If it is an output value,
//    it will be searched for in the range: [0,1.0D+4]
//    Output, int *STATUS, reports on the calculation.
//    0, if calculation completed correctly;
//    -I, if input parameter number I is out of range;
//    1, if the answer appears to be lower than the lowest search bound;
//    2, if the answer appears to be higher than the greatest search bound.
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
//
{
	const pdt5pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * x = s_data->a3;
	ityp * df = s_data->a4;
	ityp * pnonc = s_data->a5;
	short * status = s_data->a6;
	ityp * bound = s_data->a7; 
	
    # define tent4 1.0e4
    # define tol (1.0e-8)
    # define atol (1.0e-50)
    # define zero (1.0e-300)
    # define one (1.0e0-1.0e-16)
    # define inf 1.0e300

    static ityp K1 = 0.0e0;
    static ityp K3 = 0.5e0;
    static ityp K4 = 5.0e0;
    static ityp fx,cum,ccum;
    static unsigned long qhi,qleft;
    static ityp T2,T5,T6,T7,T8,T9,T10,T11,T12,T13;

    *status = 0;
    *bound = 0.00;
    //     Check arguments
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 4.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p < 0.0e0 || *p > one)) goto S60;
        if(!(*p < 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = one;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60:
        if(*which == 2) goto S90;
        //     X
        if(!(*x < 0.0e0)) goto S80;
        *bound = 0.0e0;
        *status = -4;
        return NULL;
    S90:
    S80:
        if(*which == 3) goto S110;
        //     DF
        if(!(*df <= 0.0e0)) goto S100;
        *bound = 0.0e0;
        *status = -5;
        return NULL;
    S110:
    S100:
        if(*which == 4) goto S130;
        //     PNONC
        if(!(*pnonc < 0.0e0)) goto S120;
        *bound = 0.0e0;
        *status = -6;
        return NULL;
    S130:
    S120:
        //     Calculate ANSWERS
        if(1 == *which)
        {
            //     Calculating P and Q
            cumchn(x,df,pnonc,p,q);
            *status = 0;
        }
        else if(2 == *which)
        {
            //     Calculating X
            *x = 5.0e0;
            T2 = inf;
            T5 = atol;
            T6 = tol;
            dstinv(&K1,&T2,&K3,&K3,&K4,&T5,&T6);
            *status = 0;
            dinvr(status,x,&fx,&qleft,&qhi);
            S140:
                if(!(*status == 1)) goto S150;
                cumchn(x,df,pnonc,&cum,&ccum);
                fx = cum-*p;
                dinvr(status,x,&fx,&qleft,&qhi);
                goto S140;
            S150:
                if(!(*status == -1)) goto S180;
                if(!qleft) goto S160;
                *status = 1;
                *bound = 0.0e0;
                goto S170;
            S160:
                *status = 2;
                *bound = inf;
            S180:
            S170:
                ;
        }
        else if(3 == *which)
        {
            //     Calculating DF
            *df = 5.0e0;
            T7 = zero;
            T8 = inf;
            T9 = atol;
            T10 = tol;
            dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
            *status = 0;
            dinvr(status,df,&fx,&qleft,&qhi);
            S190:
                if(!(*status == 1)) goto S200;
                cumchn(x,df,pnonc,&cum,&ccum);
                fx = cum-*p;
                dinvr(status,df,&fx,&qleft,&qhi);
                goto S190;
            S200:
                if(!(*status == -1)) goto S230;
                if(!qleft) goto S210;
                *status = 1;
                *bound = zero;
                goto S220;
            S210:
                *status = 2;
                *bound = inf;
            S230:
            S220:
                ;
        }
        else if(4 == *which)
        {
            //     Calculating PNONC
            *pnonc = 5.0e0;
            T11 = tent4;
            T12 = atol;
            T13 = tol;
            dstinv(&K1,&T11,&K3,&K3,&K4,&T12,&T13);
            *status = 0;
            dinvr(status,pnonc,&fx,&qleft,&qhi);
            S240:
                if(!(*status == 1)) goto S250;
                cumchn(x,df,pnonc,&cum,&ccum);
                fx = cum-*p;
                dinvr(status,pnonc,&fx,&qleft,&qhi);
                goto S240;
            S250:
                if(!(*status == -1)) goto S280;
                if(!qleft) goto S260;
                *status = 1;
                *bound = zero;
                goto S270;
            S260:
                *status = 2;
                *bound = tent4;
            S270:
                ;
        }
        S280:
            return NULL;
    # undef tent4
    # undef tol
    # undef atol
    # undef zero
    # undef one
    # undef inf
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cdff ( void * data)
/******************************************************************************/
//  Purpose:
//    CDFF evaluates the CDF of the F distribution.
//  Discussion:
//    This routine calculates any one parameter of the F distribution
//    given the others.
//    The value P of the cumulative distribution function is calculated
//    directly.
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//    The value of the cumulative F distribution is not necessarily
//    monotone in either degree of freedom.  There thus may be two
//    values that provide a given CDF value.  This routine assumes
//    monotonicity and will find an arbitrary one of the two values.
//  Modified:
//    14 April 2007
//  Reference:
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.6.2.
//  Parameters:
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from F, DFN and DFD;
//    2: Calculate F from P, Q, DFN and DFD;
//    3: Calculate DFN from P, Q, F and DFD;
//    4: Calculate DFD from P, Q, F and DFN.
//    Input/output, double *P, the integral from 0 to F of
//    the F-density.  If it is an input value, it should lie in the
//    range [0,1].
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//    Input/output, double *F, the upper limit of integration
//    of the F-density.  If this is an input value, it should lie in the
//    range [0, +infinity).  If it is an output value, it will be searched
//    for in the range [0,1.0D+300].
//    Input/output, double *DFN, the number of degrees of
//    freedom of the numerator sum of squares.  If this is an input value,
//    it should lie in the range: (0, +infinity).  If it is an output value,
//    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
//    Input/output, double *DFD, the number of degrees of freedom
//    of the denominator sum of squares.  If this is an input value, it should
//    lie in the range: (0, +infinity).  If it is an output value, it will
//    be searched for in the  range: [ 1.0D-300, 1.0D+300].
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
{
	const pdt5pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * f = s_data->a3;
	ityp * dfn = s_data->a4;
	ityp * dfd = s_data->a5;
	short * status = s_data->a6;
	ityp * bound = s_data->a7; 
	
    # define tol (1.0e-8)
    # define atol (1.0e-50)
    # define zero (1.0e-300)
    # define inf 1.0e300

    static dim_typ K1 = 1;
    static ityp K2 = 0.0e0;
    static ityp K4 = 0.5e0;
    static ityp K5 = 5.0e0;
    static ityp pq,fx,cum,ccum;
    static dim_typ qporq;
    static unsigned long qleft, qhi;
    static ityp T3,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15;

    *status = 0;
    *bound = 0.00;
    //  Check arguments
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 4.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
        if(!(*p < 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = 1.0e0;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60:
        if(*which == 1) goto S110;
        //     Q
        if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
        if(!(*q <= 0.0e0)) goto S80;
        *bound = 0.0e0;
        goto S90;
    S80:
        *bound = 1.0e0;
    S90:
        *status = -3;
        return NULL;
    S110:
    S100:
        if(*which == 2) goto S130;
        //     F
        if(!(*f < 0.0e0)) goto S120;
        *bound = 0.0e0;
        *status = -4;
        return NULL;
    S130:
    S120:
        if(*which == 3) goto S150;
        //     DFN
        if(!(*dfn <= 0.0e0)) goto S140;
        *bound = 0.0e0;
        *status = -5;
        return NULL;
    S150:
    S140:
        if(*which == 4) goto S170;
        //     DFD
        if(!(*dfd <= 0.0e0)) goto S160;
        *bound = 0.0e0;
        *status = -6;
        return NULL;
    S170:
    S160:
        if(*which == 1) goto S210;
        //     P + Q
        pq = *p+*q;
        if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0 * dpmpar ( &K1 ) ) ) goto S200;
        if(!(pq < 0.0e0)) goto S180;
        *bound = 0.0e0;
        goto S190;
    S180:
        *bound = 1.0e0;
    S190:
        *status = 3;
        return NULL;
    S210:
    S200:
        if(!(*which == 1)) qporq = *p <= *q;
        //     Select the minimum of P or Q
        //     Calculate ANSWERS
        if(1 == *which)
        {
            //     Calculating P
            cumf(f,dfn,dfd,p,q);
            *status = 0;
        }
        else if(2 == *which)
        {
            //     Calculating F
            *f = 5.0e0;
            T3 = inf;
            T6 = atol;
            T7 = tol;
            dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
            *status = 0;
            dinvr(status,f,&fx,&qleft,&qhi);
            S220:
                if(!(*status == 1)) goto S250;
                cumf(f,dfn,dfd,&cum,&ccum);
                if(!qporq) goto S230;
                fx = cum-*p;
                goto S240;
            S230:
                fx = ccum-*q;
            S240:
                dinvr(status,f,&fx,&qleft,&qhi);
                goto S220;
                S250:
                if(!(*status == -1)) goto S280;
                if(!qleft) goto S260;
                *status = 1;
                *bound = 0.0e0;
                goto S270;
            S260:
                *status = 2;
                *bound = inf;
            S280:
            S270:
                ;
        }
        //  Calculate DFN.
        //  Note that, in the original calculation, the lower bound for DFN was 0.
        //  Using DFN = 0 causes an error in CUMF when it calls BETA_INC.
        //  The lower bound was set to the more reasonable value of 1.
        //  JVB, 14 April 2007.
        else if ( 3 == *which )
        {

            T8 = 1.0;
            T9 = inf;
            T10 = atol;
            T11 = tol;
            dstinv ( &T8, &T9, &K4, &K4, &K5, &T10, &T11 );

            *status = 0;
            *dfn = 5.0;
            fx = 0.0;

            dinvr ( status, dfn, &fx, &qleft, &qhi );

            while ( *status == 1 )
            {
                cumf ( f, dfn, dfd, &cum, &ccum );
                fx = *p <= *q ? cum - *p : ccum - *q;
                dinvr ( status, dfn, &fx, &qleft, &qhi );
            }

            if ( *status == -1 )
            {
            if ( qleft )
            {
                *status = 1;
                *bound = 1.00;
            }
            else
            {
                *status = 2;
                *bound = inf;
            }
        }
    }
    //  Calculate DFD.
    //  Note that, in the original calculation, the lower bound for DFD was 0.
    //  Using DFD = 0 causes an error in CUMF when it calls BETA_INC.
    //  The lower bound was set to the more reasonable value of 1.
    //  JVB, 14 April 2007.
    else if ( 4 == *which )
    {

        T12 = 1.0;
        T13 = inf;
        T14 = atol;
        T15 = tol;
        dstinv ( &T12, &T13, &K4, &K4, &K5, &T14, &T15 );

        *status = 0;
        *dfd = 5.0;
        fx = 0.0;
        dinvr ( status, dfd, &fx, &qleft, &qhi );

        while ( *status == 1 )
        {
            cumf ( f, dfn, dfd, &cum, &ccum );
            fx = *p <= *q ? cum - *p : ccum - *q;
            dinvr ( status, dfd, &fx, &qleft, &qhi );
        }

        if ( *status == -1 )
        {
            if ( qleft )
            {
                *status = 1;
                *bound = 1.00;
            }
            else
            {
                *status = 2;
                *bound = inf;
            }
        }
    }

    return NULL;
    # undef tol
    # undef atol
    # undef zero
    # undef inf
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cdffnc ( void * data)
/******************************************************************************/
//  Purpose:
//    CDFFNC evaluates the CDF of the Noncentral F distribution.
//  Discussion:
//    This routine originally used 1.0E+300 as the upper bound for the
//    interval in which many of the missing parameters are to be sought.
//    Since the underlying rootfinder routine needs to evaluate the
//    function at this point, it is no surprise that the program was
//    experiencing overflows.  A less extravagant upper bound
//    is being tried for now!
//    This routine calculates any one parameter of the Noncentral F distribution
//    given the others.
//    The value P of the cumulative distribution function is calculated
//    directly.
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//    The computation time required for this routine is proportional
//    to the noncentrality parameter PNONC.  Very large values of
//    this parameter can consume immense computer resources.  This is
//    why the search range is bounded by 10,000.
//    The value of the cumulative noncentral F distribution is not
//    necessarily monotone in either degree of freedom.  There thus
//    may be two values that provide a given CDF value.  This routine
//    assumes monotonicity and will find an arbitrary one of the two
//    values.
//    The CDF of the noncentral F distribution can be evaluated
//    within Mathematica by commands such as:
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF [ NoncentralFRatioDistribution [ DFN, DFD, PNONC ], X ]
//  Modified:
//    15 June 2004
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.6.20.
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//  Parameters:
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from F, DFN, DFD and PNONC;
//    2: Calculate F from P, Q, DFN, DFD and PNONC;
//    3: Calculate DFN from P, Q, F, DFD and PNONC;
//    4: Calculate DFD from P, Q, F, DFN and PNONC;
//    5: Calculate PNONC from P, Q, F, DFN and DFD.
//    Input/output, double *P, the integral from 0 to F of
//    the noncentral F-density.  If P is an input value it should
//    lie in the range [0,1) (Not including 1!).
//    Dummy, double *Q, is not used by this subroutine,
//    and is only included for similarity with the other routines.
//    Its input value is not checked.  If P is to be computed, the
//    Q is set to 1 - P.
//    Input/output, double *F, the upper limit of integration
//    of the noncentral F-density.  If this is an input value, it should
//    lie in the range: [0, +infinity).  If it is an output value, it
//    will be searched for in the range: [0,1.0D+30].
//    Input/output, double *DFN, the number of degrees of freedom
//    of the numerator sum of squares.  If this is an input value, it should
//    lie in the range: (0, +infinity).  If it is an output value, it will
//    be searched for in the range: [ 1.0, 1.0D+30].
//    Input/output, double *DFD, the number of degrees of freedom
//    of the denominator sum of squares.  If this is an input value, it should
//    be in range: (0, +infinity).  If it is an output value, it will be
//    searched for in the range [1.0, 1.0D+30].
//    Input/output, double *PNONC, the noncentrality parameter
//    If this is an input value, it should be nonnegative.
//    If it is an output value, it will be searched for in the range: [0,1.0D+4].
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
{
	const pdt6pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * f = s_data->a3;
	ityp * dfn = s_data->a4;
	ityp * dfd = s_data->a5;
	ityp * phonc = s_data->a6;
	short * status = s_data->a7;
	ityp * bound = s_data->a8; 
	
    # define tent4 1.0e4
    # define tol (1.0e-8)
    # define atol (1.0e-50)
    # define zero (1.0e-300)
    # define one (1.0e0-1.0e-16)
    # define inf 1.0e300

    static ityp K1 = 0.0e0;
    static ityp K3 = 0.5e0;
    static ityp K4 = 5.0e0;
    static ityp fx,cum,ccum;
    static unsigned long qhi,qleft;
    static ityp T2,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17;

    *status = 0;
    *bound = 0.00;
    //     Check arguments
    if(!(*which < 1 || *which > 5)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 5.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p < 0.0e0 || *p > one)) goto S60;
        if(!(*p < 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = one;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60:
        if(*which == 2) goto S90;
        //     F
        if(!(*f < 0.0e0)) goto S80;
        *bound = 0.0e0;
        *status = -4;
        return NULL;
    S90:
    S80:
        if(*which == 3) goto S110;
        //     DFN
        if(!(*dfn <= 0.0e0)) goto S100;
        *bound = 0.0e0;
        *status = -5;
        return NULL;
    S110:
    S100:
        if(*which == 4) goto S130;
        //     DFD
        if(!(*dfd <= 0.0e0)) goto S120;
        *bound = 0.0e0;
        *status = -6;
        return NULL;
    S130:
    S120:
        if(*which == 5) goto S150;
        //     PHONC
        if(!(*phonc < 0.0e0)) goto S140;
        *bound = 0.0e0;
        *status = -7;
        return NULL;
    S150:
    S140:
        //     Calculate ANSWERS
        if(1 == *which)
        {
            //     Calculating P
            cumfnc(f,dfn,dfd,phonc,p,q);
            *status = 0;
        }
        else if(2 == *which)
        {
            //     Calculating F
            *f = 5.0e0;
            T2 = inf;
            T5 = atol;
            T6 = tol;
            dstinv(&K1,&T2,&K3,&K3,&K4,&T5,&T6);
            *status = 0;
            dinvr(status,f,&fx,&qleft,&qhi);
            S160:
                if(!(*status == 1)) goto S170;
                cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
                fx = cum-*p;
                dinvr(status,f,&fx,&qleft,&qhi);
                goto S160;
            S170:
                if(!(*status == -1)) goto S200;
                if(!qleft) goto S180;
                *status = 1;
                *bound = 0.0e0;
                goto S190;
            S180:
                *status = 2;
                *bound = inf;
            S200:
            S190:
                ;
        }
        else if(3 == *which)
        {
            //     Calculating DFN
            *dfn = 5.0e0;
            T7 = zero;
            T8 = inf;
            T9 = atol;
            T10 = tol;
            dstinv(&T7,&T8,&K3,&K3,&K4,&T9,&T10);
            *status = 0;
            dinvr(status,dfn,&fx,&qleft,&qhi);
            S210:
                if(!(*status == 1)) goto S220;
                cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
                fx = cum-*p;
                dinvr(status,dfn,&fx,&qleft,&qhi);
                goto S210;
            S220:
                if(!(*status == -1)) goto S250;
                if(!qleft) goto S230;
                *status = 1;
                *bound = zero;
                goto S240;
            S230:
                *status = 2;
                *bound = inf;
            S250:
            S240:
                ;
        }
        else if(4 == *which)
        {
            //     Calculating DFD

            *dfd = 5.0e0;
            T11 = zero;
            T12 = inf;
            T13 = atol;
            T14 = tol;
            dstinv(&T11,&T12,&K3,&K3,&K4,&T13,&T14);
            *status = 0;
            dinvr(status,dfd,&fx,&qleft,&qhi);
            S260:
                if(!(*status == 1)) goto S270;
                cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
                fx = cum-*p;
                dinvr(status,dfd,&fx,&qleft,&qhi);
                goto S260;
            S270:
                if(!(*status == -1)) goto S300;
                if(!qleft) goto S280;
                *status = 1;
                *bound = zero;
                goto S290;
            S280:
                *status = 2;
                *bound = inf;
            S300:
            S290:
                ;
        }
        else if(5 == *which)
        {
            //     Calculating PHONC
            *phonc = 5.0e0;
            T15 = tent4;
            T16 = atol;
            T17 = tol;
            dstinv(&K1,&T15,&K3,&K3,&K4,&T16,&T17);
            *status = 0;
            dinvr(status,phonc,&fx,&qleft,&qhi);
            S310:
                if(!(*status == 1)) goto S320;
                cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
                fx = cum-*p;
                dinvr(status,phonc,&fx,&qleft,&qhi);
                goto S310;
            S320:
                if(!(*status == -1)) goto S350;
                if(!qleft) goto S330;
                *status = 1;
                *bound = 0.0e0;
                goto S340;
            S330:
                *status = 2;
                *bound = tent4;
            S340:
                ;
        }
    S350:
        return NULL;
    # undef tent4
    # undef tol
    # undef atol
    # undef zero
    # undef one
    # undef inf
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _cdfgam ( void * data)
/******************************************************************************/
//  Purpose:
//    CDFGAM evaluates the CDF of the Gamma Distribution.
//  Discussion:
//    This routine calculates any one parameter of the Gamma distribution
//    given the others.
//    The cumulative distribution function P is calculated directly.
//    Computation of the other parameters involves a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//    The gamma density is proportional to T**(SHAPE - 1) * EXP(- SCALE * T)
//  Reference:
//    Armido DiDinato and Alfred Morris,
//    Computation of the incomplete gamma function ratios and their inverse,
//    ACM Transactions on Mathematical Software,
//    Volume 12, 1986, pages 377-393.
//  Parameters:
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X, SHAPE and SCALE;
//    2: Calculate X from P, Q, SHAPE and SCALE;
//    3: Calculate SHAPE from P, Q, X and SCALE;
//    4: Calculate SCALE from P, Q, X and SHAPE.
//    Input/output, double *P, the integral from 0 to X of the
//    Gamma density.  If this is an input value, it should lie in the
//    range: [0,1].
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//    Input/output, double *X, the upper limit of integration of
//    the Gamma density.  If this is an input value, it should lie in the
//    range: [0, +infinity).  If it is an output value, it will lie in
//    the range: [0,1E300].
//    Input/output, double *SHAPE, the shape parameter of the
//    Gamma density.  If this is an input value, it should lie in the range:
// (0, +infinity).  If it is an output value, it will be searched for
//    in the range: [1.0D-300,1.0D+300].
//    Input/output, double *SCALE, the scale parameter of the
//    Gamma density.  If this is an input value, it should lie in the range
// (0, +infinity).  If it is an output value, it will be searched for
//    in the range: (1.0D-300,1.0D+300].
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +10, if the Gamma or inverse Gamma routine cannot compute the answer.
//    This usually happens only for X and SHAPE very large (more than 1.0D+10.
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
{
	const pdt5pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * x = s_data->a3;
	ityp * shape = s_data->a4;
	ityp * scale = s_data->a5;
	short * status = s_data->a6;
	ityp * bound = s_data->a7; 
	
    # define tol (1.0e-8)
    # define atol (1.0e-50)
    # define zero (1.0e-300)
    # define inf 1.0e300

    static dim_typ K1 = 1;
    static ityp K5 = 0.5e0;
    static ityp K6 = 5.0e0;
    static ityp xx,fx,xscale,cum,ccum,pq,porq;
    static dim_typ ierr;
    static dim_typ qporq;
    static unsigned long qhi, qleft;
    static ityp T2,T3,T4,T7,T8,T9;

    *status = 0;
    *bound = 0.00;
    //     Check arguments
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 4.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
        if(!(*p < 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = 1.0e0;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60:
        if(*which == 1) goto S110;
        //     Q
        if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
        if(!(*q <= 0.0e0)) goto S80;
        *bound = 0.0e0;
        goto S90;
    S80:
        *bound = 1.0e0;
    S90:
        *status = -3;
        return NULL;
    S110:
    S100:
        if(*which == 2) goto S130;
        //     X
        if(!(*x < 0.0e0)) goto S120;
        *bound = 0.0e0;
        *status = -4;
        return NULL;
    S130:
    S120:
        if(*which == 3) goto S150;
        //     SHAPE
        if(!(*shape <= 0.0e0)) goto S140;
        *bound = 0.0e0;
        *status = -5;
        return NULL;
    S150:
    S140:
        if(*which == 4) goto S170;
        //     SCALE
        if(!(*scale <= 0.0e0)) goto S160;
        *bound = 0.0e0;
        *status = -6;
        return NULL;
    S170:
    S160:
        if(*which == 1) goto S210;
        //     P + Q
        pq = *p+*q;
        if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S200;
        if(!(pq < 0.0e0)) goto S180;
        *bound = 0.0e0;
        goto S190;
    S180:
        *bound = 1.0e0;
    S190:
        *status = 3;
        return NULL;
    S210:
    S200:
        if(*which == 1) goto S240;
        //     Select the minimum of P or Q
        qporq = *p <= *q;
        if(!qporq) goto S220;
        porq = *p;
        goto S230;
    S220:
        porq = *q;
    S240:
    S230:
        //     Calculate ANSWERS
        if(1 == *which)
        {
            //     Calculating P
            *status = 0;
            xscale = *x**scale;
            cumgam(&xscale,shape,p,q);
            if(porq > 1.5e0) *status = 10;
        }
        else if(2 == *which)
        {
            //     Computing X
            T2 = -1.0e0;
            gamma_inc_inv ( shape, &xx, &T2, p, q, &ierr );
            if(ierr < 0.0e0)
            {
                *status = 10;
                return NULL;
            }
            else
            {
                *x = xx/ *scale;
                *status = 0;
            }
        }
        else if(3 == *which)
        {
            //     Computing SHAPE
            *shape = 5.0e0;
            xscale = *x**scale;
            T3 = zero;
            T4 = inf;
            T7 = atol;
            T8 = tol;
            dstinv(&T3,&T4,&K5,&K5,&K6,&T7,&T8);
            *status = 0;
            dinvr(status,shape,&fx,&qleft,&qhi);
            S250:
                if(!(*status == 1)) goto S290;
                cumgam(&xscale,shape,&cum,&ccum);
                if(!qporq) goto S260;
                fx = cum-*p;
                goto S270;
            S260:
                fx = ccum-*q;
            S270:
                if(!(qporq && cum > 1.5e0 || !qporq && ccum > 1.5e0)) goto S280;
                *status = 10;
                return NULL;
            S280:
                dinvr(status,shape,&fx,&qleft,&qhi);
                goto S250;
            S290:
                if(!(*status == -1)) goto S320;
                if(!qleft) goto S300;
                *status = 1;
                *bound = zero;
                goto S310;
            S300:
                *status = 2;
                *bound = inf;
            S320:
            S310:
                ;
        }
        else if(4 == *which)
        {
            //     Computing SCALE
            T9 = -1.0e0;
            gamma_inc_inv ( shape, &xx, &T9, p, q, &ierr );
            if(ierr < 0.0e0)
            {
                *status = 10;
                return NULL;
            }
            else
            {
                *scale = xx/ *x;
                *status = 0;
            }
        }
            return NULL;
    # undef tol
    # undef atol
    # undef zero
    # undef inf
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _cdfnbn ( void * data )
/******************************************************************************/
//  Purpose:
//    CDFNBN evaluates the CDF of the Negative Binomial distribution
//  Discussion:
//    This routine calculates any one parameter of the negative binomial
//    distribution given values for the others.
//    The cumulative negative binomial distribution returns the
//    probability that there will be F or fewer failures before the
//    S-th success in binomial trials each of which has probability of
//    success PR.
//    The individual term of the negative binomial is the probability of
//    F failures before S successes and is
//    Choose( F, S+F-1 ) * PR^(S) * (1-PR)^F
//    Computation of other parameters involve a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.26.
//  Parameters:
//    Input, int WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from F, S, PR and OMPR;
//    2: Calculate F from P, Q, S, PR and OMPR;
//    3: Calculate S from P, Q, F, PR and OMPR;
//    4: Calculate PR and OMPR from P, Q, F and S.
//    Input/output, double P, the cumulation from 0 to F of
//    the negative binomial distribution.  If P is an input value, it
//    should lie in the range [0,1].
//    Input/output, double Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//    Input/output, double F, the upper limit of cumulation of
//    the binomial distribution.  There are F or fewer failures before
//    the S-th success.  If this is an input value, it may lie in the
//    range [0,+infinity), and if it is an output value, it will be searched
//    for in the range [0,1.0D+300].
//    Input/output, double S, the number of successes.
//    If this is an input value, it should lie in the range: [0, +infinity).
//    If it is an output value, it will be searched for in the range:
//    [0, 1.0D+300].
//    Input/output, double PR, the probability of success in each
//    binomial trial.  Whether an input or output value, it should lie in the
//    range [0,1].
//    Input/output, double OMPR, the value of (1-PR).  Whether an
//    input or output value, it should lie in the range [0,1].
//    Output, int STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1;
//    +4, if PR + OMPR /= 1.
//    Output, double BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
{
	const pdt6pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * s = s_data->a3;
	ityp * xn = s_data->a4;
	ityp * pr = s_data->a5;
	ityp * ompr = s_data->a6;
	short * status = s_data->a7;
	ityp * bound = s_data->a8; 
	
    # define tol (1.0e-8)
    # define atol (1.0e-50)
    # define inf 1.0e300
    # define one 1.0e0

    static dim_typ K1 = 1;
    static ityp K2 = 0.0e0;
    static ityp K4 = 0.5e0;
    static ityp K5 = 5.0e0;
    static ityp K11 = 1.0e0;
    static ityp fx,xhi,xlo,pq,prompr,cum,ccum;
    static dim_typ qporq;
    static unsigned long qhi, qleft;
    static ityp T3,T6,T7,T8,T9,T10,T12,T13;

    *status = 0;
    *bound = 0.00;
    //     Check arguments
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 4.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
        if(!(*p < 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = 1.0e0;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60:
        if(*which == 1) goto S110;
        //     Q
        if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
        if(!(*q <= 0.0e0)) goto S80;
        *bound = 0.0e0;
        goto S90;
    S80:
        *bound = 1.0e0;
    S90:
        *status = -3;
        return NULL;
    S110:
    S100:
        if(*which == 2) goto S130;
        //     S
        if(!(*s < 0.0e0)) goto S120;
        *bound = 0.0e0;
        *status = -4;
        return NULL;
    S130:
    S120:
        if(*which == 3) goto S150;
        //     XN
        if(!(*xn < 0.0e0)) goto S140;
        *bound = 0.0e0;
        *status = -5;
        return NULL;
    S150:
    S140:
        if(*which == 4) goto S190;
        //     PR
        if(!(*pr < 0.0e0 || *pr > 1.0e0)) goto S180;
        if(!(*pr < 0.0e0)) goto S160;
        *bound = 0.0e0;
        goto S170;
    S160:
        *bound = 1.0e0;
    S170:
        *status = -6;
        return NULL;
    S190:
    S180:
        if(*which == 4) goto S230;
        //     OMPR
        if(!(*ompr < 0.0e0 || *ompr > 1.0e0)) goto S220;
        if(!(*ompr < 0.0e0)) goto S200;
        *bound = 0.0e0;
        goto S210;
    S200:
        *bound = 1.0e0;
    S210:
        *status = -7;
        return NULL;
    S230:
    S220:
        if(*which == 1) goto S270;
        //     P + Q
        pq = *p+*q;
        if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S260;
        if(!(pq < 0.0e0)) goto S240;
        *bound = 0.0e0;
        goto S250;
    S240:
        *bound = 1.0e0;
    S250:
        *status = 3;
        return NULL;
    S270:
    S260:
        if(*which == 4) goto S310;
        //     PR + OMPR
        prompr = *pr+*ompr;
        if(!(fabs(prompr-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S300;
        if(!(prompr < 0.0e0)) goto S280;
        *bound = 0.0e0;
        goto S290;
    S280:
        *bound = 1.0e0;
    S290:
        *status = 4;
        return NULL;
    S310:
    S300:
        if(!(*which == 1)) qporq = *p <= *q;
        //     Select the minimum of P or Q
        //     Calculate ANSWERS
        if(1 == *which)
        {
            //     Calculating P
            cumnbn(s,xn,pr,ompr,p,q);
            *status = 0;
        }
        else if(2 == *which)
        {
            //     Calculating S
            *s = 5.0e0;
            T3 = inf;
            T6 = atol;
            T7 = tol;
            dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
            *status = 0;
            dinvr(status,s,&fx,&qleft,&qhi);
            S320:
                if(!(*status == 1)) goto S350;
                cumnbn(s,xn,pr,ompr,&cum,&ccum);
                if(!qporq) goto S330;
                fx = cum-*p;
                goto S340;
            S330:
                fx = ccum-*q;
            S340:
                dinvr(status,s,&fx,&qleft,&qhi);
                goto S320;
            S350:
                if(!(*status == -1)) goto S380;
                if(!qleft) goto S360;
                *status = 1;
                *bound = 0.0e0;
                goto S370;
            S360:
                *status = 2;
                *bound = inf;
            S380:
            S370:
                ;
        }
        else if(3 == *which)
        {
            //     Calculating XN
            *xn = 5.0e0;
            T8 = inf;
            T9 = atol;
            T10 = tol;
            dstinv(&K2,&T8,&K4,&K4,&K5,&T9,&T10);
            *status = 0;
            dinvr(status,xn,&fx,&qleft,&qhi);
            S390:
                if(!(*status == 1)) goto S420;
                cumnbn(s,xn,pr,ompr,&cum,&ccum);
                if(!qporq) goto S400;
                fx = cum-*p;
                goto S410;
            S400:
                fx = ccum-*q;
            S410:
                dinvr(status,xn,&fx,&qleft,&qhi);
                goto S390;
            S420:
                if(!(*status == -1)) goto S450;
                if(!qleft) goto S430;
                *status = 1;
                *bound = 0.0e0;
                goto S440;
            S430:
                *status = 2;
                *bound = inf;
            S450:
            S440:
            ;
        }
        else if(4 == *which)
        {
            //     Calculating PR and OMPR
            T12 = atol;
            T13 = tol;
            dstzr(&K2,&K11,&T12,&T13);
            if(!qporq) goto S480;
            *status = 0;
            dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
            *ompr = one-*pr;
            S460:
                if(!(*status == 1)) goto S470;
                cumnbn(s,xn,pr,ompr,&cum,&ccum);
                fx = cum-*p;
                dzror(status,pr,&fx,&xlo,&xhi,&qleft,&qhi);
                *ompr = one-*pr;
                goto S460;
            S470:
                goto S510;
            S480:
                *status = 0;
                dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
                *pr = one-*ompr;
            S490:
                if(!(*status == 1)) goto S500;
                cumnbn(s,xn,pr,ompr,&cum,&ccum);
                fx = ccum-*q;
                dzror(status,ompr,&fx,&xlo,&xhi,&qleft,&qhi);
                *pr = one-*ompr;
                goto S490;
            S510:
            S500:
                if(!(*status == -1)) goto S540;
                if(!qleft) goto S520;
                *status = 1;
                *bound = 0.0e0;
                goto S530;
            S520:
                *status = 2;
                *bound = 1.0e0;
            S530:
                ;
        }
        S540:
            return NULL;
    # undef tol
    # undef atol
    # undef inf
    # undef one
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cdfnor ( void * data)
/******************************************************************************/
//  Purpose:
//    CDFNOR evaluates the CDF of the Normal distribution.
//  Discussion:
//    A slightly modified version of ANORM from SPECFUN
//    is used to calculate the cumulative standard normal distribution.
//    The rational functions from pages 90-95 of Kennedy and Gentle
//    are used as starting values to Newton's Iterations which
//    compute the inverse standard normal.  Therefore no searches are
//    necessary for any parameter.
//    For X < -15, the asymptotic expansion for the normal is used  as
//    the starting value in finding the inverse standard normal.
//    The normal density is proportional to
//    exp( - 0.5D+00 * (( X - MEAN)/SD)**2)
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.2.12.
//    William Cody,
//    Algorithm 715: SPECFUN - A Portable FORTRAN Package of
//      Special Function Routines and Test Drivers,
//    ACM Transactions on Mathematical Software,
//    Volume 19, pages 22-32, 1993.
//    Kennedy and Gentle,
//    Statistical Computing,
//    Marcel Dekker, NY, 1980,
//    QA276.4  K46
//  Parameters:
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from X, MEAN and SD;
//    2: Calculate X from P, Q, MEAN and SD;
//    3: Calculate MEAN from P, Q, X and SD;
//    4: Calculate SD from P, Q, X and MEAN.
//    Input/output, double *P, the integral from -infinity to X
//    of the Normal density.  If this is an input or output value, it will
//    lie in the range [0,1].
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//    Input/output, double *X, the upper limit of integration of
//    the Normal density.
//    Input/output, double *MEAN, the mean of the Normal density.
//    Input/output, double *SD, the standard deviation of the
//    Normal density.  If this is an input value, it should lie in the
//    range (0,+infinity).
//    Output, int *STATUS, the status of the calculation.
//    0, if calculation completed correctly;
//    -I, if input parameter number I is out of range;
//    1, if answer appears to be lower than lowest search bound;
//    2, if answer appears to be higher than greatest search bound;
//    3, if P + Q /= 1.
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
{
	const pdt5pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * x = s_data->a3;
	ityp * mean = s_data->a4;
	ityp * sd = s_data->a5;
	short * status = s_data->a6;
	ityp * bound = s_data->a7; 
	
    static dim_typ K1 = 1;
    static ityp z,pq;

    *status = 0;
    *bound = 0.00;
    //     Check arguments
    *status = 0;
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 4.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
        if(!(*p <= 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = 1.0e0;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60:
        if(*which == 1) goto S110;
        //     Q
        if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
        if(!(*q <= 0.0e0)) goto S80;
        *bound = 0.0e0;
        goto S90;
    S80:
        *bound = 1.0e0;
    S90:
        *status = -3;
        return NULL;
    S110:
    S100:
        if(*which == 1) goto S150;
        //     P + Q
        pq = *p+*q;
        if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S140;
        if(!(pq < 0.0e0)) goto S120;
        *bound = 0.0e0;
        goto S130;
    S120:
        *bound = 1.0e0;
    S130:
        *status = 3;
        return NULL;
    S150:
    S140:
        if(*which == 4) goto S170;
        //     SD
        if(!(*sd <= 0.0e0)) goto S160;
        *bound = 0.0e0;
        *status = -6;
        return NULL;
    S170:
    S160:
        //     Calculate ANSWERS
        if(1 == *which)
        {
            //     Computing P
            z = (*x-*mean)/ *sd;
            cumnor(&z,p,q);
        }
        else if(2 == *which)
        {
            //     Computing X
            z = dinvnr(p,q);
            *x = *sd*z+*mean;
        }
        else if(3 == *which)
        {
            //     Computing the MEAN
            z = dinvnr(p,q);
            *mean = *x-*sd*z;
        }
        else if(4 == *which)
        {
            //     Computing SD
            z = dinvnr(p,q);
            *sd = (*x-*mean)/z;
        }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cdfpoi ( void * data)
/******************************************************************************/
//  Purpose:
//    CDFPOI evaluates the CDF of the Poisson distribution.
//  Discussion:
//    This routine calculates any one parameter of the Poisson distribution
//    given the others.
//    The value P of the cumulative distribution function is calculated
//    directly.
//    Computation of other parameters involve a seach for a value that
//    produces the desired value of P.  The search relies on the
//    monotonicity of P with respect to the other parameters.
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.4.21.
//  Parameters:
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1: Calculate P and Q from S and XLAM;
//    2: Calculate A from P, Q and XLAM;
//    3: Calculate XLAM from P, Q and S.
//    Input/output, double *P, the cumulation from 0 to S of the
//    Poisson density.  Whether this is an input or output value, it will
//    lie in the range [0,1].
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//    Input/output, double *S, the upper limit of cumulation of
//    the Poisson CDF.  If this is an input value, it should lie in
//    the range: [0, +infinity).  If it is an output value, it will be
//    searched for in the range: [0,1.0D+300].
//    Input/output, double *XLAM, the mean of the Poisson
//    distribution.  If this is an input value, it should lie in the range
//    [0, +infinity).  If it is an output value, it will be searched for
//    in the range: [0,1E300].
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
{
	const pdt4pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * s = s_data->a3;
	ityp * xlam = s_data->a4;
	short * status = s_data->a5;
	ityp * bound = s_data->a6; 
	
    # define tol (1.0e-8)
    # define atol (1.0e-50)
    # define inf 1.0e300

    static dim_typ K1 = 1;
    static ityp K2 = 0.0e0;
    static ityp K4 = 0.5e0;
    static ityp K5 = 5.0e0;
    static ityp fx,cum,ccum,pq;
    static dim_typ qporq;
    static unsigned long qhi, qleft;
    static ityp T3,T6,T7,T8,T9,T10;

    *status = 0;
    *bound = 0.00;
    //     Check arguments
    if(!(*which < 1 || *which > 3)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 3.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
        if(!(*p < 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = 1.0e0;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60:
        if(*which == 1) goto S110;
        //     Q
        if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
        if(!(*q <= 0.0e0)) goto S80;
        *bound = 0.0e0;
        goto S90;
    S80:
        *bound = 1.0e0;
    S90:
        *status = -3;
        return NULL;
    S110:
    S100:
        if(*which == 2) goto S130;
        //     S
        if(!(*s < 0.0e0)) goto S120;
        *bound = 0.0e0;
        *status = -4;
        return NULL;
    S130:
    S120:
        if(*which == 3) goto S150;
        //     XLAM
        if(!(*xlam < 0.0e0)) goto S140;
        *bound = 0.0e0;
        *status = -5;
        return NULL;
    S150:
    S140:
        if(*which == 1) goto S190;
        //     P + Q
        pq = *p+*q;
        if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S180;
        if(!(pq < 0.0e0)) goto S160;
        *bound = 0.0e0;
        goto S170;
    S160:
        *bound = 1.0e0;
    S170:
        *status = 3;
        return NULL;
    S190:
    S180:
        if(!(*which == 1)) qporq = *p <= *q;
        //     Select the minimum of P or Q
        //     Calculate ANSWERS
        if(1 == *which)
        {
            //     Calculating P
            cumpoi(s,xlam,p,q);
            *status = 0;
        }
        else if(2 == *which)
        {
            //     Calculating S
            *s = 5.0e0;
            T3 = inf;
            T6 = atol;
            T7 = tol;
            dstinv(&K2,&T3,&K4,&K4,&K5,&T6,&T7);
            *status = 0;
            dinvr(status,s,&fx,&qleft,&qhi);
            S200:
                if(!(*status == 1)) goto S230;
                cumpoi(s,xlam,&cum,&ccum);
                if(!qporq) goto S210;
                fx = cum-*p;
                goto S220;
            S210:
                fx = ccum-*q;
            S220:
                dinvr(status,s,&fx,&qleft,&qhi);
                goto S200;
            S230:
                if(!(*status == -1)) goto S260;
                if(!qleft) goto S240;
                *status = 1;
                *bound = 0.0e0;
                goto S250;
            S240:
                *status = 2;
                *bound = inf;
            S260:
            S250:
                ;
        }
        else if(3 == *which)
        {
            //     Calculating XLAM
            *xlam = 5.0e0;
            T8 = inf;
            T9 = atol;
            T10 = tol;
            dstinv(&K2,&T8,&K4,&K4,&K5,&T9,&T10);
            *status = 0;
            dinvr(status,xlam,&fx,&qleft,&qhi);
            S270:
                if(!(*status == 1)) goto S300;
                cumpoi(s,xlam,&cum,&ccum);
                if(!qporq) goto S280;
                fx = cum-*p;
                goto S290;
            S280:
                fx = ccum-*q;
            S290:
                dinvr(status,xlam,&fx,&qleft,&qhi);
                goto S270;
            S300:
                if(!(*status == -1)) goto S330;
                if(!qleft) goto S310;
                *status = 1;
                *bound = 0.0e0;
                goto S320;
            S310:
                *status = 2;
                *bound = inf;
            S320:
                ;
        }
        S330:
            return NULL;
    # undef tol
    # undef atol
    # undef inf
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cdft ( void * data)
/******************************************************************************/
//  Purpose:
//    CDFT evaluates the CDF of the T distribution.
//  Discussion:
//    This routine calculates any one parameter of the T distribution
//    given the others.
//    The value P of the cumulative distribution function is calculated
//    directly.
//    Computation of other parameters involve a seach for a value that
//    produces the desired value of P.   The search relies on the
//    monotonicity of P with respect to the other parameters.
//    The original version of this routine allowed the search interval
//    to extend from -1.0E+300 to +1.0E+300, which is fine until you
//    try to evaluate a function at such a point!
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.27.
//  Parameters:
//    Input, int *WHICH, indicates which argument is to be calculated
//    from the others.
//    1 : Calculate P and Q from T and DF;
//    2 : Calculate T from P, Q and DF;
//    3 : Calculate DF from P, Q and T.
//    Input/output, double *P, the integral from -infinity to T of
//    the T-density.  Whether an input or output value, this will lie in the
//    range [0,1].
//    Input/output, double *Q, equal to 1-P.  If Q is an input
//    value, it should lie in the range [0,1].  If Q is an output value,
//    it will lie in the range [0,1].
//    Input/output, double *T, the upper limit of integration of
//    the T-density.  If this is an input value, it may have any value.
//    It it is an output value, it will be searched for in the range
//    [ -1.0D+30, 1.0D+30 ].
//    Input/output, double *DF, the number of degrees of freedom
//    of the T distribution.  If this is an input value, it should lie
//    in the range: (0 , +infinity).  If it is an output value, it will be
//    searched for in the range: [1, 1.0D+10].
//    Output, int *STATUS, reports the status of the computation.
//     0, if the calculation completed correctly;
//    -I, if the input parameter number I is out of range;
//    +1, if the answer appears to be lower than lowest search bound;
//    +2, if the answer appears to be higher than greatest search bound;
//    +3, if P + Q /= 1.
//    Output, double *BOUND, is only defined if STATUS is nonzero.
//    If STATUS is negative, then this is the value exceeded by parameter I.
//    if STATUS is 1 or 2, this is the search bound that was exceeded.
{
	const pdt4pitpspit * const s_data = data;
	dim_typ * which = s_data->a0;
	ityp * p = s_data->a1;
	ityp * q =  s_data->a2;
	ityp * t = s_data->a3;
	ityp * df = s_data->a4;
	short * status = s_data->a5;
	ityp * bound = s_data->a6;
	
    # define tol (1.0e-8)
    # define atol (1.0e-50)
    # define zero (1.0e-300)
    # define inf 1.0e30
    # define maxdf 1.0e10

    static dim_typ K1 = 1;
    static ityp K4 = 0.5e0;
    static ityp K5 = 5.0e0;
    static ityp fx,cum,ccum,pq;
    static dim_typ qporq;
    static unsigned long qhi, qleft;
    static ityp T2,T3,T6,T7,T8,T9,T10,T11;

    *status = 0;
    *bound = 0.00;
    //     Check arguments
    if(!(*which < 1 || *which > 3)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
        *bound = 3.0e0;
    S20:
        *status = -1;
        return NULL;
    S30:
        if(*which == 1) goto S70;
        //     P
        if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
        if(!(*p <= 0.0e0)) goto S40;
        *bound = 0.0e0;
        goto S50;
    S40:
        *bound = 1.0e0;
    S50:
        *status = -2;
        return NULL;
    S70:
    S60:
        if(*which == 1) goto S110;
        //     Q
        if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
        if(!(*q <= 0.0e0)) goto S80;
        *bound = 0.0e0;
        goto S90;
    S80:
        *bound = 1.0e0;
    S90:
        *status = -3;
        return NULL;
    S110:
    S100:
        if(*which == 3) goto S130;
        //     DF
        if(!(*df <= 0.0e0)) goto S120;
        *bound = 0.0e0;
        *status = -5;
        return NULL;
    S130:
    S120:
        if(*which == 1) goto S170;
        //     P + Q
        pq = *p+*q;
        if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*dpmpar(&K1))) goto S160;
        if(!(pq < 0.0e0)) goto S140;
        *bound = 0.0e0;
        goto S150;
    S140:
        *bound = 1.0e0;
    S150:
        *status = 3;
        return NULL;
    S170:
    S160:
        if(!(*which == 1)) qporq = *p <= *q;
        //     Select the minimum of P or Q
        //     Calculate ANSWERS
        if(1 == *which)
        {
            //     Computing P and Q
            cumt(t,df,p,q);
            *status = 0;
        }
        else if(2 == *which)
        {
            //     Computing T
            //     .. Get initial approximation for T
            *t = dt1(p,q,df);
            T2 = -inf;
            T3 = inf;
            T6 = atol;
            T7 = tol;
            dstinv(&T2,&T3,&K4,&K4,&K5,&T6,&T7);
            *status = 0;
            dinvr(status,t,&fx,&qleft,&qhi);
            S180:
                if(!(*status == 1)) goto S210;
                cumt(t,df,&cum,&ccum);
                if(!qporq) goto S190;
                fx = cum-*p;
                goto S200;
            S190:
                fx = ccum-*q;
            S200:
                dinvr(status,t,&fx,&qleft,&qhi);
                goto S180;
            S210:
                if(!(*status == -1)) goto S240;
                if(!qleft) goto S220;
                *status = 1;
                *bound = -inf;
                goto S230;
            S220:
                *status = 2;
                *bound = inf;
            S240:
            S230:
                ;
            }
            else if(3 == *which)
            {
                //     Computing DF
                *df = 5.0e0;
                T8 = zero;
                T9 = maxdf;
                T10 = atol;
                T11 = tol;
                dstinv(&T8,&T9,&K4,&K4,&K5,&T10,&T11);
                *status = 0;
                dinvr(status,df,&fx,&qleft,&qhi);
                S250:
                    if(!(*status == 1)) goto S280;
                    cumt(t,df,&cum,&ccum);
                    if(!qporq) goto S260;
                    fx = cum-*p;
                    goto S270;
                S260:
                    fx = ccum-*q;
                S270:
                    dinvr(status,df,&fx,&qleft,&qhi);
                    goto S250;
                S280:
                    if(!(*status == -1)) goto S310;
                    if(!qleft) goto S290;
                    *status = 1;
                    *bound = zero;
                    goto S300;
                S290:
                    *status = 2;
                    *bound = maxdf;
                S300:
                    ;
            }
    S310:
        return NULL;
    # undef tol
    # undef atol
    # undef zero
    # undef inf
    # undef maxdf
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cumbet ( void * data)
/******************************************************************************/
//  Purpose:
//    CUMBET evaluates the cumulative incomplete beta distribution.
//  Discussion:
//    This routine calculates the CDF to X of the incomplete beta distribution
//    with parameters A and B.  This is the integral from 0 to x
//    of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
//  Modified:
//    14 March 2006
//  Reference:
//    A R Didonato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios.
//    ACM Transactions on Mathematical Software,
//    Volume 18, Number 3, September 1992, pages 360-373.
//  Parameters:
//    Input, double *X, the upper limit of integration.
//    Input, double *Y, the value of 1-X.
//    Input, double *A, *B, the parameters of the distribution.
//    Output, double *CUM, *CCUM, the values of the cumulative
//    density function and complementary cumulative density function.
{
	ityp ** const a_data = data;
	ityp * x = a_data[0];
	ityp * y = a_data[1];
	ityp * a = a_data[2];
	ityp * b = a_data[3];
	ityp * cum = a_data[4];
	ityp * ccum = a_data[5];
		
    static dim_typ ierr;
    if ( *x <= 0.00 )
    {
        *cum = 0.00;
        *ccum = 1.0;
    }
    else if ( *y <= 0.00 )
    {
        *cum = 1.00;
        *ccum = 0.00;
    }
    else
        beta_inc ( a, b, x, y, cum, ccum, &ierr );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _cumbin ( void * data)
/******************************************************************************/
//  Purpose:
//    CUMBIN evaluates the cumulative binomial distribution.
//  Discussion:
//    This routine returns the probability of 0 to S successes in XN binomial
//    trials, each of which has a probability of success, PR.
//  Modified:
//    14 March 2006
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.24.
//  Parameters:
//    Input, double *S, the upper limit of summation.
//    Input, double *XN, the number of trials.
//    Input, double *PR, the probability of success in one trial.
//    Input, double *OMPR, equals ( 1 - PR ).
//    Output, double *CUM, the cumulative binomial distribution.
//    Output, double *CCUM, the complement of the cumulative
//    binomial distribution.
{
	ityp ** const a_data = data;
	ityp * s = a_data[0];
	ityp * xn = a_data[1];
	ityp * pr = a_data[2];
	ityp * ompr = a_data[3];
	ityp * cum = a_data[4];
	ityp * ccum = a_data[5];
	
    static double T1,T2;

    if ( *s < *xn )
    {
        T1 = *s + 1.00;
        T2 = *xn - *s;
        cumbet ( pr, ompr, &T1, &T2, ccum, cum );
    }
    else
    {
        *cum = 1.00;
        *ccum = 0.00;
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cumchi ( void * data)
/******************************************************************************/
//  Purpose:
//    CUMCHI evaluates the cumulative chi-square distribution.
//  Parameters:
//    Input, double *X, the upper limit of integration.
//    Input, double *DF, the degrees of freedom of the
//    chi-square distribution.
//    Output, double *CUM, the cumulative chi-square distribution.
//    Output, double *CCUM, the complement of the cumulative
//    chi-square distribution.
{
	ityp ** const a_data = data;
	ityp * x = a_data[0];
	ityp * df = a_data[1];
	ityp * cum = a_data[2];
	ityp * ccum = a_data[3];
	
    static ityp a;
    static ityp xx;
    a = *df * 0.50;
    xx = *x * 0.50;
    cumgam ( &xx, &a, cum, ccum );
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cumchn ( void * data)
/******************************************************************************/
//  Purpose:
//    CUMCHN evaluates the cumulative noncentral chi-square distribution.
//  Discussion:
//    Calculates the cumulative noncentral chi-square
//    distribution, i.e., the probability that a random variable
//    which follows the noncentral chi-square distribution, with
//    noncentrality parameter PNONC and continuous degrees of
//    freedom DF, is less than or equal to X.
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.4.25.
//  Parameters:
//    Input, double *X, the upper limit of integration.
//    Input, double *DF, the number of degrees of freedom.
//    Input, double *PNONC, the noncentrality parameter of
//    the noncentral chi-square distribution.
//    Output, double *CUM, *CCUM, the CDF and complementary
//    CDF of the noncentral chi-square distribution.
//  Local Parameters:
//    Local, double EPS, the convergence criterion.  The sum
//    stops when a term is less than EPS*SUM.
//    Local, int NTIRED, the maximum number of terms to be evaluated
//    in each sum.
//    Local, bool QCONV, is TRUE if convergence was achieved, that is,
//    the program did not stop on NTIRED criterion.
{
	ityp ** const a_data = data;
	ityp * x = a_data[0];
	ityp * df = a_data[1];
	ityp * pnonc = a_data[2];
	ityp * cum = a_data[3];
	ityp * ccum = a_data[4];
	
    # define dg(i) (*df+2.0e0*(double)(i))
    # define qsmall(xx) (int)(sum < 1.0e-20 || (xx) < eps*sum)
    # define qtired(i) (int)((i) > ntired)

    static ityp eps = 1.0e-5;
    static dim_typ ntired = 1000;
    static ityp adj,centaj,centwt,chid2,dfd2,lcntaj,lcntwt,lfact,pcent,pterm,sum,sumadj,term,wt,xnonc;
    static dim_typ i,icent,iterb,iterf;
    static ityp T1,T2,T3;

    if(!(*x <= 0.0e0)) goto S10;
    *cum = 0.0e0;
    *ccum = 1.0e0;
    return NULL;
    S10:
        if(!(*pnonc <= 1.0e-10)) goto S20;
        //     When non-centrality parameter is (essentially) zero,
        //     use cumulative chi-square distribution
        cumchi(x,df,cum,ccum);
        return NULL;
    S20:
        xnonc = *pnonc/2.0e0;
        //     The following code calculates the weight, chi-square, and
        //     adjustment term for the central term in the infinite series.
        //     The central term is the one in which the poisson weight is
        //     greatest.  The adjustment term is the amount that must
        //     be subtracted from the chi-square to move up two degrees
        //     of freedom.
        icent = fifidint(xnonc);
        if(icent == 0) icent = 1;
        chid2 = *x/2.0e0;
        //     Calculate central weight term
        T1 = (ityp)(icent+1);
        lfact = gamma_log ( &T1 );
        lcntwt = -xnonc+(ityp)icent*log(xnonc)-lfact;
        centwt = exp(lcntwt);
        //     Calculate central chi-square
        T2 = dg(icent);
        cumchi(x,&T2,&pcent,ccum);
        //     Calculate central adjustment term
        dfd2 = dg(icent)/2.0e0;
        T3 = 1.0e0+dfd2;
        lfact = gamma_log ( &T3 );
        lcntaj = dfd2*log(chid2)-chid2-lfact;
        centaj = exp(lcntaj);
        sum = centwt*pcent;
        //     Sum backwards from the central term towards zero.
        //     Quit whenever either
        //  (1) the zero term is reached, or
        //  (2) the term gets small relative to the sum, or
        //  (3) More than NTIRED terms are totaled.
        iterb = 0;
        sumadj = 0.0e0;
        adj = centaj;
        wt = centwt;
        i = icent;
        goto S40;
    S30:
        if(qtired(iterb) || qsmall(term) || i == 0) goto S50;
    S40:
        dfd2 = dg(i)/2.0e0;
        //     Adjust chi-square for two fewer degrees of freedom.
        //     The adjusted value ends up in PTERM.
        adj *= dfd2/chid2;
        sumadj += adj;
        pterm = pcent+sumadj;
        //     Adjust poisson weight for J decreased by one
        wt *= ((ityp)i/xnonc);
        term = wt*pterm;
        sum += term;
        i -= 1;
        ++ iterb;
        goto S30;
    S50:
        iterf = 0;
        //     Now sum forward from the central term towards infinity.
        //     Quit when either
        //  (1) the term gets small relative to the sum, or
        //  (2) More than NTIRED terms are totaled.
        sumadj = adj = centaj;
        wt = centwt;
        i = icent;
        goto S70;
    S60:
        if(qtired(iterf) || qsmall(term)) goto S80;
    S70:
        //     Update weights for next higher J
        wt *= (xnonc/(ityp)(i+1));
        //     Calculate PTERM and add term to sum
        pterm = pcent-sumadj;
        term = wt*pterm;
        sum += term;
        //  Update adjustment term for DF for next iteration
        ++ i;
        dfd2 = dg(i)/2.0e0;
        adj = adj*chid2/dfd2;
        sumadj = sum + adj;
        ++ iterf;
        goto S60;
    S80:
        *cum = sum;
        *ccum = 0.5e0+(0.5e0-*cum);
    return NULL;
    # undef dg
    # undef qsmall
    # undef qtired
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cumf ( void * data)
/******************************************************************************/
//  Purpose:
//    CUMF evaluates the cumulative F distribution.
//  Discussion:
//    CUMF computes the integral from 0 to F of the F density with DFN
//    numerator and DFD denominator degrees of freedom.
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.28.
//  Parameters:
//    Input, double *F, the upper limit of integration.
//    Input, double *DFN, *DFD, the number of degrees of
//    freedom for the numerator and denominator.
//    Output, double *CUM, *CCUM, the value of the F CDF and
//    the complementary F CDF.
{
	ityp ** const a_data = data;
	ityp * f = a_data[0];
	ityp * dfn = a_data[1];
	ityp * dfd = a_data[2];
	ityp * cum = a_data[3];
	ityp * ccum = a_data[4];
	
    # define half 0.5e0
    # define done 1.0e0

    static ityp dsum,prod,xx,yy;
    static dim_typ ierr;
    static ityp T1,T2;

    if(!(*f <= 0.0e0)) goto S10;
    *cum = 0.0e0;
    *ccum = 1.0e0;
    return NULL;
    S10:
        prod = *dfn**f;
        //     XX is such that the incomplete beta with parameters
        //     DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
        //     YY is 1 - XX
        //     Calculate the smaller of XX and YY accurately
        dsum = *dfd+prod;
        xx = *dfd/dsum;

        if ( xx > half )
        {
            yy = prod/dsum;
            xx = done-yy;
        }
        else
            yy = done - xx;

        T1 = *dfd*half;
        T2 = *dfn*half;
        beta_inc ( &T1, &T2, &xx, &yy, ccum, cum, &ierr );
        return NULL;
    # undef half
    # undef done
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cumfnc ( void * data)
/******************************************************************************/
//  Purpose:
//    CUMFNC evaluates the cumulative noncentral F distribution.
//  Discussion:
//    This routine computes the noncentral F distribution with DFN and DFD
//    degrees of freedom and noncentrality parameter PNONC.
//    The series is calculated backward and forward from J = LAMBDA/2
// (this is the term with the largest Poisson weight) until
//    the convergence criterion is met.
//    The sum continues until a succeeding term is less than EPS
//    times the sum (or the sum is less than 1.0e-20).  EPS is
//    set to 1.0e-4 in a data statement which can be changed.
//    The original version of this routine allowed the input values
//    of DFN and DFD to be negative (nonsensical) or zero (which
//    caused numerical overflow.)  I have forced both these values
//    to be at least 1.
//  Modified:
//    15 June 2004
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.16, 26.6.17, 26.6.18, 26.6.20.
//  Parameters:
//    Input, double *F, the upper limit of integration.
//    Input, double *DFN, *DFD, the number of degrees of freedom
//    in the numerator and denominator.  Both DFN and DFD must be positive,
//    and normally would be integers.  This routine requires that they
//    be no less than 1.
//    Input, double *PNONC, the noncentrality parameter.
//    Output, double *CUM, *CCUM, the noncentral F CDF and
//    complementary CDF.
{
	ityp ** const a_data = data;
	ityp * f = a_data[0];
	ityp * dfn = a_data[1];
	ityp * dfd = a_data[2];
	ityp * pnonc = a_data[3];
	ityp * cum = a_data[4];
	ityp * ccum = a_data[5];
	
    # define qsmall(x) (int)(sum < 1.0e-20 || (x) < eps*sum)
    # define half 0.5e0
    # define done 1.0e0

    static ityp eps = 1.0e-4;
    static ityp dsum,dummy,prod,xx,yy,adn,aup,b,betdn,betup,centwt,dnterm,sum,upterm,xmult,xnonc;
    static dim_typ i,icent,ierr;
    static ityp T1,T2,T3,T4,T5,T6;

    if(!(*f <= 0.0e0)) goto S10;
    *cum = 0.0e0;
    *ccum = 1.0e0;
    return NULL;
    S10:
        if(!(*pnonc < 1.0e-10)) goto S20;
        //  Handle case in which the non-centrality parameter is
        // (essentially) zero.
        cumf(f,dfn,dfd,cum,ccum);
        return NULL;
    S20:
        xnonc = *pnonc/2.0e0;
        //  Calculate the central term of the poisson weighting factor.
        icent = ( dim_typ ) xnonc;
        if(icent == 0) icent = 1;
        //  Compute central weight term
        T1 = (ityp)(icent+1);
        centwt = exp(-xnonc+(ityp)icent*log(xnonc)- gamma_log ( &T1 ) );
        //  Compute central incomplete beta term
        //  Assure that minimum of arg to beta and 1 - arg is computed
        //  accurately.
        prod = *dfn**f;
        dsum = *dfd+prod;
        yy = *dfd/dsum;
        if(yy > half)
        {
            xx = prod/dsum;
            yy = done-xx;
        }
        else
            xx = done-yy;
        T2 = *dfn*half+(ityp)icent;
        T3 = *dfd*half;
        beta_inc ( &T2, &T3, &xx, &yy, &betdn, &dummy, &ierr );
        adn = *dfn/2.0e0+(ityp)icent;
        aup = adn;
        b = *dfd/2.0e0;
        betup = betdn;
        sum = centwt*betdn;
        //  Now sum terms backward from icent until convergence or all done
        xmult = centwt;
        i = icent;
        T4 = adn+b;
        T5 = adn+1.0e0;
        dnterm = exp( gamma_log ( &T4 ) - gamma_log ( &T5 ) - gamma_log ( &b ) + adn * log ( xx ) + b * log(yy));
    S30:
        if(qsmall(xmult*betdn) || i <= 0) goto S40;
        xmult *= ((ityp)i/xnonc);
        i -= 1;
        adn -= 1.00;
        dnterm = (adn+1.00)/((adn+b)*xx)*dnterm;
        betdn += dnterm;
        sum += (xmult*betdn);
        goto S30;
    S40:
        i = icent+1;
        //  Now sum forwards until convergence
        xmult = centwt;
        if(aup-1.0+b == 0) upterm = exp(-gamma_log ( &aup ) - gamma_log ( &b ) + (aup-1.00)*log(xx)+b*log(yy));
        else
        {
            T6 = aup-1.00+b;
            upterm = exp( gamma_log ( &T6 ) - gamma_log ( &aup )- gamma_log ( &b ) + (aup-1.00)*log(xx)+b*log(yy));
        }
        goto S60;
    S50:
        if(qsmall(xmult*betup)) goto S70;
    S60:
        xmult *= (xnonc/(ityp)i);
        i += 1;
        aup += 1.00;
        upterm = (aup+b-2.0e0)*xx/(aup-1.00)*upterm;
        betup -= upterm;
        sum += (xmult*betup);
        goto S50;
    S70:
        *cum = sum;
        *ccum = 0.5e0+(0.5e0-*cum);
    return NULL;
    # undef qsmall
    # undef half
    # undef done
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cumgam ( void * data)
/******************************************************************************/
//  Purpose:
//    CUMGAM evaluates the cumulative incomplete gamma distribution.
//  Discussion:
//    This routine computes the cumulative distribution function of the
//   incomplete gamma distribution, i.e., the integral from 0 to X of
//   (1/GAM(A))*EXP(-T)*T**(A-1) DT
//    where GAM(A) is the complete gamma function of A, i.e.,
//      GAM(A) = integral from 0 to infinity of EXP(-T)*T**(A-1) DT
//  Parameters:
//    Input, double *X, the upper limit of integration.
//    Input, double *A, the shape parameter of the incomplete
//    Gamma distribution.
//    Output, double *CUM, *CCUM, the incomplete Gamma CDF and
//    complementary CDF.
{
	ityp ** const a_data = data;
	ityp * x = a_data[0];
	ityp * a = a_data[1];
	ityp * cum = a_data[2];
	ityp * ccum = a_data[3];
	
    static dim_typ K1 = 0;

    if(!(*x <= 0.0e0)) goto S10;
    *cum = 0.0e0;
    *ccum = 1.0e0;
    return NULL;
    S10:
        gamma_inc ( a, x, cum, ccum, &K1 );
        //     Call gratio routine
    return NULL;
}
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cumnbn ( void * data)
/******************************************************************************/
//  Purpose:
//    CUMNBN evaluates the cumulative negative binomial distribution.
//  Discussion:
//    This routine returns the probability that there will be F or
//    fewer failures before there are S successes, with each binomial
//    trial having a probability of success PR.
//    Prob(# failures = F | S successes, PR)  =
//                  ( S + F - 1 )
//                  (            ) * PR^S * (1-PR)^F
//                  (      F     )
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.5.26.
//  Parameters:
//    Input, double *F, the number of failures.
//    Input, double *S, the number of successes.
//    Input, double *PR, *OMPR, the probability of success on
//    each binomial trial, and the value of (1-PR).
//    Output, double *CUM, *CCUM, the negative binomial CDF,
//    and the complementary CDF.
{
	ityp ** const a_data = data;
	ityp * s = a_data[0];
	ityp * xn = a_data[1];
	ityp * pr = a_data[2];
	ityp * ompr = a_data[3];
	ityp * cum = a_data[4];
	ityp * ccum = a_data[5];
	
    static double T1;

    T1 = *s+1.e0;
    cumbet(pr,ompr,xn,&T1,cum,ccum);
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cumnor ( void * data)
/******************************************************************************/
//  Purpose:
//    CUMNOR computes the cumulative normal distribution.
//  Discussion:
//    This function evaluates the normal distribution function:
//                              / x
//                     1       |       -t*t/2
//          P(x) = ----------- |      e       dt
//                 sqrt(2 M_PI)  |
//                             /-oo
//    This transportable program uses rational functions that
//    theoretically approximate the normal distribution function to
//    at least 18 significant decimal digits.  The accuracy achieved
//    depends on the arithmetic system, the compiler, the intrinsic
//    functions, and proper selection of the machine-dependent
//    constants.
//  Author:
//    William Cody
//    Mathematics and Computer Science Division
//    Argonne National Laboratory
//    Argonne, IL 60439
//  Reference:
//    William Cody,
//    Rational Chebyshev approximations for the error function,
//    Mathematics of Computation,
//    1969, pages 631-637.
//    William Cody,
//    Algorithm 715:
//    SPECFUN - A Portable FORTRAN Package of Special Function Routines
//      and Test Drivers,
//    ACM Transactions on Mathematical Software,
//    Volume 19, 1993, pages 22-32.
//  Parameters:
//    Input, double *ARG, the upper limit of integration.
//    Output, double *CUM, *CCUM, the Normal density CDF and
//    complementary CDF.
//  Local Parameters:
//    Local, double EPS, the argument below which anorm(x)
//    may be represented by 0.5D+00 and above which  x*x  will not underflow.
//    A conservative value is the largest machine number X
//    such that   1.0D+00 + X = 1.0D+00   to machine precision.
{
	ityp ** const a_data = data;
	ityp * arg = a_data[0];
	ityp * result = a_data[1];
	ityp * ccum = a_data[2];
	
    static ityp a[5] =
    {
        2.2352520354606839287e00,
        1.6102823106855587881e02,
        1.0676894854603709582e03,
        1.8154981253343561249e04,
        6.5682337918207449113e-2
    };
    static ityp b[4] =
    {
        4.7202581904688241870e01,
        9.7609855173777669322e02,
        1.0260932208618978205e04,
        4.5507789335026729956e04
    };
    static ityp c[9] =
    {
        3.9894151208813466764e-1,
        8.8831497943883759412e00,
        9.3506656132177855979e01,
        5.9727027639480026226e02,
        2.4945375852903726711e03,
        6.8481904505362823326e03,
        1.1602651437647350124e04,
        9.8427148383839780218e03,
        1.0765576773720192317e-8
    };
    static ityp d[8] =
    {
        2.2266688044328115691e01,
        2.3538790178262499861e02,
        1.5193775994075548050e03,
        6.4855582982667607550e03,
        1.8615571640885098091e04,
        3.4900952721145977266e04,
        3.8912003286093271411e04,
        1.9685429676859990727e04
    };
    static ityp half = 0.5e0;
    static ityp p[6] =
    {
        2.1589853405795699e-1,
        1.274011611602473639e-1,
        2.2235277870649807e-2,
        1.421619193227893466e-3,
        2.9112874951168792e-5,
        2.307344176494017303e-2
    };
    static ityp one = 1.0e0;
    static ityp q[5] =
    {
        1.28426009614491121e00,
        4.68238212480865118e-1,
        6.59881378689285515e-2,
        3.78239633202758244e-3,
        7.29751555083966205e-5
    };
    static ityp sixten = 1.60e0;
    static ityp sqrpi = 3.9894228040143267794e-1;
    static ityp thrsh = 0.66291e0;
    static ityp root32 = 5.656854248e0;
    static ityp zero = 0.0e0;
    static dim_typ K1 = 1;
    static dim_typ K2 = 2;
    static dim_typ i;
    static ityp del,eps,temp,x,xden,xnum,y,xsq,MIN;
    //  Machine dependent constants
    eps = dpmpar(&K1)*0.5e0;
    MIN = dpmpar(&K2);
    x = *arg;
    y = fabs(x);
    if(y <= thrsh)
    {
        //  Evaluate  anorm  for  |X| <= 0.66291
        xsq = zero;
        if(y > eps) xsq = x*x;
        xnum = a[4]*xsq;
        xden = xsq;
        #pragma omp parallel for num_threads(3)
        for ( i = 0; i < 3; ++i)
        {
            xnum = (xnum+a[i])*xsq;
            xden = (xden+b[i])*xsq;
        }
        *result = x*(xnum+a[3])/(xden+b[3]);
        temp = *result;
        *result = half+temp;
        *ccum = half-temp;
    }
    //  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
    else if(y <= root32)
    {
        xnum = c[8]*y;
        xden = y;
        #pragma omp parallel for num_threads(7)
        for ( i = 0; i < 7; ++i )
        {
            xnum = (xnum+c[i])*y;
            xden = (xden+d[i])*y;
        }
        *result = (xnum+c[7])/(xden+d[7]);
        xsq = fifdint(y*sixten)/sixten;
        del = (y-xsq)*(y+xsq);
        *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
        *ccum = one-*result;
        if(x > zero)
        {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
    //  Evaluate  anorm  for |X| > sqrt(32)
    else
    {
        *result = zero;
        xsq = one/(x*x);
        xnum = p[5]*xsq;
        xden = xsq;
        #pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i )
        {
            xnum = (xnum+p[i])*xsq;
            xden = (xden+q[i])*xsq;
        }
        *result = xsq*(xnum+p[4])/(xden+q[4]);
        *result = (sqrpi-*result)/y;
        xsq = fifdint(x*sixten)/sixten;
        del = (x-xsq)*(x+xsq);
        *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
        *ccum = one-*result;
        if(x > zero)
        {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
    	}

    }
    if(*result < MIN) *result = 0.0e0;
    //  Fix up for negative argument, erf, etc.
    if(*ccum < MIN) *ccum = 0.0e0;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _cumpoi (void * data)
/******************************************************************************/
//  Purpose:
//    CUMPOI evaluates the cumulative Poisson distribution.
//  Discussion:
//    CUMPOI returns the probability of S or fewer events in a Poisson
//    distribution with mean XLAM.
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    Formula 26.4.21.
//  Parameters:
//    Input, double *S, the upper limit of cumulation of the
//    Poisson density function.
//    Input, double *XLAM, the mean of the Poisson distribution.
//    Output, double *CUM, *CCUM, the Poisson density CDF and
//    complementary CDF.
{
	ityp ** const a_data = data;
	ityp * s = a_data[0];
	ityp * xlam = a_data[1];
	ityp * cum = a_data[2];
	ityp * ccum = a_data[3];
	
    static ityp chi,df;
    df = 2.0e0*(*s+1.0e0);
    chi = 2.0e0**xlam;
    cumchi(&chi,&df,ccum,cum);
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _cumt ( void * data)
/******************************************************************************/
//  Purpose:
//    CUMT evaluates the cumulative T distribution.
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    Formula 26.5.27.
//  Parameters:
//    Input, double *T, the upper limit of integration.
//    Input, double *DF, the number of degrees of freedom of
//    the T distribution.
//    Output, double *CUM, *CCUM, the T distribution CDF and
//    complementary CDF.
{
	ityp ** const a_data = data;
	ityp * t = a_data[0];
	ityp * df = a_data[1];
	ityp * cum = a_data[2];
	ityp * ccum = a_data[3];
	
    static ityp a;
    static ityp dfptt;
    static ityp K2 = 0.5e0;
    static ityp oma;
    static ityp T1;
    static ityp tt;
    static ityp xx;
    static ityp yy;

    tt = (*t) * (*t);
    dfptt = ( *df ) + tt;
    xx = *df / dfptt;
    yy = tt / dfptt;
    T1 = 0.5e0 * ( *df );
    cumbet ( &xx, &yy, &T1, &K2, &a, &oma );

    if ( *t <= 0.0e0 )
    {
        *cum = 0.5e0 * a;
        *ccum = oma + ( *cum );
    }
    else
    {
        *ccum = 0.5e0 * a;
        *cum = oma + ( *ccum );
    }
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dbetrm ( void * data)
/******************************************************************************/
//  Purpose:
//    DBETRM computes the Sterling remainder for the complete beta function.
//  Discussion:
//    Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
//    where Lgamma is the log of the (complete) gamma function
//    Let ZZ be approximation obtained if each log gamma is approximated
//    by Sterling's formula, i.e.,
//    Sterling(Z) = LOG( SQRT( 2*M_PI ) ) + ( Z-0.5D+00 ) * LOG( Z ) - Z
//    The Sterling remainder is Log(Beta(A,B)) - ZZ.
//  Parameters:
//    Input, double *A, *B, the parameters of the Beta function.
//    Output, double DBETRM, the Sterling remainder.
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * b = a_data[1];
	
    static ityp dbetrm,T1,T2,T3;
    //     Try to sum from smallest to largest
    T1 = *a+*b;
    dbetrm = -dstrem(&T1);
    T2 = fifdmax1(*a,*b);
    dbetrm += dstrem(&T2);
    T3 = fifdmin1(*a,*b);
    dbetrm += dstrem(&T3);
    result = dbetrm;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dexpm1 ( void * data)
/******************************************************************************/
//  Purpose:
//    DEXPM1 evaluates the function EXP(X) - 1.
//  Reference:
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//  Parameters:
//    Input, double *X, the value at which exp(X)-1 is desired.
//    Output, double DEXPM1, the value of exp(X)-1.
{
	static ityp result = MAX_VAL;
	
	ityp * x = (ityp * ) data;
	
    static ityp p1 = .914041914819518e-09;
    static ityp p2 = .238082361044469e-01;
    static ityp q1 = -.499999999085958e+00;
    static ityp q2 = .107141568980644e+00;
    static ityp q3 = -.119041179760821e-01;
    static ityp q4 = .595130811860248e-03;
    static ityp dexpm1;
    ityp w;

    if ( fabs(*x) <= 0.15e0 )
        dexpm1 =   *x * ( ( (p2   * *x+ p1 ) * *x+ 1.0e0 )/((((q4   * *x+ q3 ) * *x+ q2 ) * *x+ q1 ) * *x+ 1.0e0 ) );
    else if ( *x <= 0.0e0 )
    {
        w = exp(*x);
        dexpm1 = w-0.5e0-0.5e0;
    }
    else
    {
        w = exp(*x);
        dexpm1 = w*(0.5e0+(0.5e0-1.0e0/w));
    }
    
    result = dexpm1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _dinvnr ( void * data)
/******************************************************************************/
//  Purpose:
//    DINVNR computes the inverse of the normal distribution.
//  Discussion:
//    Returns X such that CUMNOR(X)  =   P,  i.e., the  integral from -
//    infinity to X of (1/SQRT(2*M_PI)) EXP(-U*U/2) dU is P
//    The rational function on page 95 of Kennedy and Gentle is used as a start
//    value for the Newton method of finding roots.
//  Reference:
//    Kennedy and Gentle,
//    Statistical Computing,
//    Marcel Dekker, NY, 1980,
//    QA276.4  K46
//  Parameters:
//    Input, double *P, *Q, the probability, and the complementary
//    probability.
//    Output, double DINVNR, the argument X for which the
//    Normal CDF has the value P.
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * p = a_data[0];
	ityp * q = a_data[1];
	
    # define maxit 100
    # define eps (1.0e-13)
    # define r2pi 0.3989422804014326e0
    # define nhalf (-0.5e0)
    # define dennor(x) (r2pi*exp(nhalf*(x)*(x)))

    static ityp dinvnr,strtx,xcur,cum,ccum,pp,dx;
    static dim_typ i;
    static dim_typ qporq;

    //     FIND MINIMUM OF P AND Q
    qporq = *p <= *q;
    if(!qporq) goto S10;
    pp = *p;
    goto S20;
    S10:
    pp = *q;
    S20:
    //     INITIALIZATION STEP
    strtx = stvaln(&pp);
    xcur = strtx;
    //     NEWTON INTERATIONS

    for ( i = 1; i <= maxit; ++i )
    {
        cumnor(&xcur,&cum,&ccum);
        dx = (cum-pp)/dennor(xcur);
        xcur -= dx;
        if(fabs(dx/xcur) < eps) goto S40;
    }
    dinvnr = strtx;
    //     IF WE GET HERE, NEWTON HAS FAILED
    if(!qporq) dinvnr = -dinvnr;
    return &dinvnr;
    S40:
        //     IF WE GET HERE, NEWTON HAS SUCCEDED
        dinvnr = xcur;
        if(!qporq) dinvnr *= -1;
        {
        	result = dinvnr;
        	return &result;
    	}
    # undef maxit
    # undef eps
    # undef r2pi
    # undef nhalf
    # undef dennor
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _dinvr ( void * data)
/******************************************************************************/
//  Purpose:
//    DINVR bounds the zero of the function and invokes DZROR.
//  Discussion:
//    This routine seeks to find bounds on a root of the function and
//    invokes ZROR to perform the zero finding.  STINVR must have been
//    called before this routine in order to set its parameters.
//  Reference:
//    J C P Bus and T J Dekker,
//    Two Efficient Algorithms with Guaranteed Convergence for
//      Finding a Zero of a Function,
//    ACM Transactions on Mathematical Software,
//    Volume 1, Number 4, pages 330-345, 1975.
//  Parameters:
//    Input/output, integer STATUS.  At the beginning of a zero finding
//    problem, STATUS should be set to 0 and INVR invoked.  The value
//    of parameters other than X will be ignored on this call.
//    If INVR needs the function to be evaluated, it will set STATUS to 1
//    and return.  The value of the function should be set in FX and INVR
//    again called without changing any of its other parameters.
//    If INVR finishes without error, it returns with STATUS 0, and X an
//    approximate root of F(X).
//    If INVR cannot bound the function, it returns a negative STATUS and
//    sets QLEFT and QHI.
//    Output, double precision X, the value at which F(X) is to be evaluated.
//    Input, double precision FX, the value of F(X) calculated by the user
//    on the previous call, when INVR returned with STATUS = 1.
//    Output, logical QLEFT, is defined only if QMFINV returns FALSE.  In that
//    case, QLEFT is TRUE if the stepping search terminated unsucessfully
//    at SMALL, and FALSE if the search terminated unsucessfully at BIG.
//    Output, logical QHI, is defined only if QMFINV returns FALSE.  In that
//    case, it is TRUE if Y < F(X) at the termination of the search and FALSE
//    if F(X) < Y.
{
	const ps2pit2pul * const s_data = data;
	short * status = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	unsigned long * qleft = s_data->a3;
	unsigned long * qhi = s_data->a4;
	
    E0000(0,status,x,fx,qleft,qhi,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dlanor ( void * data)
/******************************************************************************/
//  Purpose:
//    DLANOR evaluates the logarithm of the asymptotic Normal CDF.
//  Discussion:
//    This routine computes the logarithm of the cumulative normal distribution
//    from abs ( x ) to infinity for  5 <= abs ( X ).
//    The relative error at X = 5 is about 0.5D-5.
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions
//    1966, Formula 26.2.12.
//  Parameters:
//    Input, double *X, the value at which the Normal CDF is to be
//    evaluated.  It is assumed that 5 <= abs ( X ).
//    Output, double DLANOR, the logarithm of the asymptotic
//    Normal CDF.
{
	static ityp result = MAX_VAL;
	
	ityp * x = data;
	
    # define dlsqpi 0.91893853320467274177e0

    static ityp coef[12] =
    {
        -1.0e0,
        3.0e0,
        -15.0e0,
        105.0e0,
        -945.0e0,
        10395.0e0,
        -135135.0e0,
        2027025.0e0,
        -34459425.0e0,
        654729075.0e0,
        13749310575.0e0,
        316234143225.0e0
    };
    static dim_typ K1 = 12;
    static ityp dlanor,approx,correc,xx,xx2,T2;

    xx = fabs(*x);
    if ( xx < 5.0e0 )
    {
    	result = DLANOR_INVALIDRETURNVALUE;
        return &result;
    }
    approx = -dlsqpi-0.5e0*xx*xx-log(xx);
    xx2 = xx*xx;
    T2 = 1.0e0/xx2;
    correc = eval_pol ( coef, &K1, &T2 ) / xx2;
    correc = alnrel ( &correc );
    dlanor = approx+correc;
    result = dlanor;
    return &result;
    # undef dlsqpi
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dpmpar ( void * data)
/******************************************************************************/
//  Purpose:
//    DPMPAR provides machine constants for double precision arithmetic.
//  Discussion:
//     DPMPAR PROVIDES THE double PRECISION MACHINE CONSTANTS FOR
//     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
//     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
//     double PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
//     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
//        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
//        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
//        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
//     WRITTEN BY
//        ALFRED H. MORRIS, JR.
//        NAVAL SURFACE WARFARE CENTER
//        DAHLGREN VIRGINIA
//     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
//     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
//     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
{
	static ityp result = MAX_VAL;
	
	dim_typ * i = data; 
	
    static int K1 = 4;
    static int K2 = 8;
    static int K3 = 9;
    static int K4 = 10;
    static ityp value,b,binv,bm1,one,w,z;
    static dim_typ emax,emin,ibeta,m;

    if(*i > 1) goto S10;
    b = ipmpar(&K1);
    m = ipmpar(&K2);
    value = pow(b,(ityp)(1-m));
    result = value;
    return &result;
    S10:
    if(*i > 2) goto S20;
    b = ipmpar(&K1);
    emin = ipmpar(&K3);
    one = 1.00;
    binv = one/b;
    w = pow(b,(ityp)(emin+2));
    value = w*binv*binv*binv;
    result = value;
    return &result;
    S20:
    ibeta = ipmpar(&K1);
    m = ipmpar(&K2);
    emax = ipmpar(&K4);
    b = ibeta;
    bm1 = ibeta-1;
    one = 1.00;
    z = pow(b,(ityp)(m-1));
    w = ((z-one)*b+bm1)/(b*z);
    z = pow(b,(ityp)(emax-2));
    value = w*z*b*b;
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dstinv (void * data)
/******************************************************************************/
//  Purpose:
//    DSTINV seeks a value X such that F(X) = Y.
//  Discussion:
//      Double Precision - SeT INverse finder - Reverse Communication
//                              Function
//     Concise Description - Given a monotone function F finds X
//     such that F(X) = Y.  Uses Reverse communication -- see invr.
//     This routine sets quantities needed by INVR.
//          More Precise Description of INVR -
//     F must be a monotone function, the results of QMFINV are
//     otherwise undefined.  QINCR must be .TRUE. if F is non-
//     decreasing and .FALSE. if F is non-increasing.
//     QMFINV will return .TRUE. if and only if F(SMALL) and
//     F(BIG) bracket Y, i. e.,
//          QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or
//          QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL)
//     if QMFINV returns .TRUE., then the X returned satisfies
//     the following condition.  let
//               TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
//     then if QINCR is .TRUE.,
//          F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X))
//     and if QINCR is .FALSE.
//          F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
//                              Arguments
//     SMALL --> The left endpoint of the interval to be
//          searched for a solution.
//                    SMALL is DOUBLE PRECISION
//     BIG --> The right endpoint of the interval to be
//          searched for a solution.
//                    BIG is DOUBLE PRECISION
//     ABSSTP, RELSTP --> The initial step size in the search
//          is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
//                    ABSSTP is DOUBLE PRECISION
//                    RELSTP is DOUBLE PRECISION
//     STPMUL --> When a step doesn't bound the zero, the step
//                size is multiplied by STPMUL and another step
//                taken.  A popular value is 2.0
//                    DOUBLE PRECISION STPMUL
//     ABSTOL, RELTOL --> Two numbers that determine the accuracy
//          of the solution.  See function for a precise definition.
//                    ABSTOL is DOUBLE PRECISION
//                    RELTOL is DOUBLE PRECISION
//                              Method
//     Compares F(X) with Y for the input value of X then uses QINCR
//     to determine whether to step left or right to bound the
//     desired x.  the initial step size is
//          MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
//     Iteratively steps right or left until it bounds X.
//     At each step which doesn't bound X, the step size is doubled.
//     The routine is careful never to step beyond SMALL or BIG.  If
//     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
//     after setting QLEFT and QHI.
//     If X is successfully bounded then Algorithm R of the paper
//     'Two Efficient Algorithms with Guaranteed Convergence for
//     Finding a Zero of a Function' by J. C. P. Bus and
//     T. J. Dekker in ACM Transactions on Mathematical
//     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
//     to find the zero of the function F(X)-Y. This is routine
//     QRZERO.
{
	ityp ** const a_data = data;
	ityp * zsmall = a_data[0];
	ityp * zbig = a_data[1];
	ityp * zabsst = a_data[2];
	ityp * zrelst = a_data[3];
	ityp * zstpmu = a_data[4];
	ityp * zabsto = a_data[5];
	ityp * zrelto = a_data[6];
	
    E0000(1,NULL,NULL,NULL,NULL,NULL,zabsst,zabsto,zbig,zrelst,zrelto,zsmall,zstpmu);
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dstrem ( void * data)
/******************************************************************************/
//  Purpose:
//    DSTREM computes the Sterling remainder ln ( Gamma ( Z ) ) - Sterling ( Z ).
//  Discussion:
//    This routine returns
//      ln ( Gamma ( Z ) ) - Sterling ( Z )
//    where Sterling(Z) is Sterling's approximation to ln ( Gamma ( Z ) ).
//    Sterling(Z) = ln ( sqrt ( 2 * M_PI ) ) + ( Z - 0.5 ) * ln ( Z ) - Z
//    If 6 <= Z, the routine uses 9 terms of a series in Bernoulli numbers,
//    with values calculated using Maple.
//    Otherwise, the difference is computed explicitly.
//  Modified:
//    14 June 2004
//  Parameters:
//    Input, double *Z, the value at which the Sterling
//    remainder is to be calculated.  Z must be positive.
//    Output, double DSTREM, the Sterling remainder.
{
	static ityp result = MAX_VAL;
	
	ityp * z = data; 
	
    # define hln2pi 0.91893853320467274178e0
    # define ncoef 10

    static ityp coef[ncoef] =
    {
        0.0e0,0.0833333333333333333333333333333e0,
        -0.00277777777777777777777777777778e0,
        0.000793650793650793650793650793651e0,
        -0.000595238095238095238095238095238e0,
        0.000841750841750841750841750841751e0,
        -0.00191752691752691752691752691753e0,
        0.00641025641025641025641025641026e0,
        -0.0295506535947712418300653594771e0,
        0.179644372368830573164938490016e0
    };
    static dim_typ K1 = 10;
    static ityp dstrem,sterl,T2;
    //    For information, here are the next 11 coefficients of the
    //    remainder term in Sterling's formula
    //            -1.39243221690590111642743221691
    //            13.4028640441683919944789510007
    //            -156.848284626002017306365132452
    //            2193.10333333333333333333333333
    //            -36108.7712537249893571732652192
    //            691472.268851313067108395250776
    //            -0.152382215394074161922833649589D8
    //            0.382900751391414141414141414141D9
    //            -0.108822660357843910890151491655D11
    //            0.347320283765002252252252252252D12
    //            -0.123696021422692744542517103493D14
    if(*z <= 0.0e0)
    {
    	result = DSTREM_INVALIDRETURNVALUE;
        return &result;
    }
    if(!(*z > 6.0e0)) goto S10;
    T2 = 1.0e0/pow(*z,2.0);
    dstrem = eval_pol ( coef, &K1, &T2 )**z;
    goto S20;
    S10:
        sterl = hln2pi+(*z-0.5e0)*log(*z)-*z;
        dstrem = gamma_log ( z ) - sterl;
    S20:
    	result = dstrem;
        return &result;
    # undef hln2pi
    # undef ncoef
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _dstzr ( void * data)
/******************************************************************************/
//  Purpose:
//    DSTXR sets quantities needed by the zero finder.
//  Discussion:
//     Double precision SeT ZeRo finder - Reverse communication version
//                              Function
//     Sets quantities needed by ZROR.  The function of ZROR
//     and the quantities set is given here.
//     Concise Description - Given a function F
//     find XLO such that F(XLO) = 0.
//          More Precise Description -
//     Input condition. F is a double function of a single
//     double argument and XLO and XHI are such that
//          F(XLO)*F(XHI)  .LE.  0.0
//     If the input condition is met, QRZERO returns .TRUE.
//     and output values of XLO and XHI satisfy the following
//          F(XLO)*F(XHI)  .LE. 0.
//          ABS(F(XLO)  .LE. ABS(F(XHI)
//          ABS(XLO-XHI)  .LE. TOL(X)
//     where
//          TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
//     If this algorithm does not find XLO and XHI satisfying
//     these conditions then QRZERO returns .FALSE.  This
//     implies that the input condition was not met.
//                              Arguments
//     XLO --> The left endpoint of the interval to be
//           searched for a solution.
//                    XLO is DOUBLE PRECISION
//     XHI --> The right endpoint of the interval to be
//           for a solution.
//                    XHI is DOUBLE PRECISION
//     ABSTOL, RELTOL --> Two numbers that determine the accuracy
//                      of the solution.  See function for a
//                      precise definition.
//                    ABSTOL is DOUBLE PRECISION
//                    RELTOL is DOUBLE PRECISION
//                              Method
//     Algorithm R of the paper 'Two Efficient Algorithms with
//     Guaranteed Convergence for Finding a Zero of a Function'
//     by J. C. P. Bus and T. J. Dekker in ACM Transactions on
//     Mathematical Software, Volume 1, no. 4 page 330
//  (Dec. '75) is employed to find the zero of F(X)-Y.
{
	ityp ** const a_data = data;
	ityp * zxlo = a_data[0];
	ityp * zxhi = a_data[1];
	ityp * zabstl = a_data[2];
	ityp * zreltl = a_data[3];
	
    E0001(1,NULL,NULL,NULL,NULL,NULL,NULL,NULL,zabstl,zreltl,zxhi,zxlo);
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _dt1 ( void * data)
/******************************************************************************/
//  Purpose:
//    DT1 computes an approximate inverse of the cumulative T distribution.
//  Discussion:
//    Returns the inverse of the T distribution function, i.e.,
//    the integral from 0 to INVT of the T density is P. This is an
//    initial approximation.
//  Parameters:
//    Input, double *P, *Q, the value whose inverse from the
//    T distribution CDF is desired, and the value (1-P).
//    Input, double *DF, the number of degrees of freedom of the
//    T distribution.
//    Output, double DT1, the approximate value of X for which
//    the T density CDF with DF degrees of freedom has value P.
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * p = a_data[0];
	ityp * q = a_data[1];
	ityp * df = a_data[2];
	
    static ityp coef[4][5] =
    {
        1.0e0,
        1.0e0,
        0.0e0,
        0.0e0,
        0.0e0,
        3.0e0,
        16.0e0,
        5.0e0,
        0.0e0,
        0.0e0,
        -15.0e0,
        17.0e0,
        19.0e0,
        3.0e0,
        0.0e0,
        -945.0e0,
        -1920.0e0,
        1482.0e0,
        776.0e0,
        79.0e0
    };
    static ityp denom[4] =
    {
        4.0e0,
        96.0e0,
        384.0e0,
        92160.0e0
    };
    static dim_typ ideg[4] =
    {
        2,
        3,
        4,
        5
    };
    static ityp dt1,denpow,sum,term,x,xp,xx;
    static dim_typ i;

    x = fabs(dinvnr(p,q));
    xx = x*x;
    sum = x;
    denpow = 1.0e0;
    #pragma omp parallel for num_threads(4)
    for ( i = 0; i < 4; ++i )
    {
        term = eval_pol ( &coef[i][0], &ideg[i], &xx ) * x;
        denpow *= *df;
        sum += (term/(denpow*denom[i]));
    }
    if(!(*p >= 0.5e0)) goto S20;
    xp = sum;
    goto S30;
    S20:
        xp = -sum;
    S30:
        dt1 = xp;
        
    result = dt1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _dzror ( void * data)
/******************************************************************************/
//  Purpose:
//    DZROR seeks the zero of a function using reverse communication.
//  Discussion:
//     Performs the zero finding.  STZROR must have been called before
//     this routine in order to set its parameters.
//                              Arguments
//     STATUS <--> At the beginning of a zero finding problem, STATUS
//                 should be set to 0 and ZROR invoked. (The value
//                 of other parameters will be ignored on this call.)
//                 When ZROR needs the function evaluated, it will set
//                 STATUS to 1 and return.  The value of the function
//                 should be set in FX and ZROR again called without
//                 changing any of its other parameters.
//                 When ZROR has finished without error, it will return
//                 with STATUS 0.  In that case (XLO,XHI) bound the answe
//                 If ZROR finds an error (which implies that F(XLO)-Y an
//                 F(XHI)-Y have the same sign, it returns STATUS -1.  In
//                 this case, XLO and XHI are undefined.
//                         INTEGER STATUS
//     X <-- The value of X at which F(X) is to be evaluated.
//                         DOUBLE PRECISION X
//     FX --> The value of F(X) calculated when ZROR returns with
//            STATUS = 1.
//                         DOUBLE PRECISION FX
//     XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
//             inverval in X containing the solution below.
//                         DOUBLE PRECISION XLO
//     XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
//             inverval in X containing the solution above.
//                         DOUBLE PRECISION XHI
//     QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
//                at XLO.  If it is .FALSE. the search terminated
//                unsucessfully at XHI.
//                    QLEFT is LOGICAL
//     QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
//              search and .FALSE. if F(X) .LT. Y at the
//              termination of the search.
//                    QHI is LOGICAL
{
	const ps4pit2pul * const s_data = data;
	short * status = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	ityp * xlo = s_data->a3;
	ityp * xhi = s_data->a4;
	unsigned long * qleft = s_data->a5;
	unsigned long * qhi = s_data->a6;
	
    E0001(0,status,x,fx,xlo,xhi,qleft,qhi,NULL,NULL,NULL,NULL);
    return NULL;	
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void * _E0000 ( void * data)
/******************************************************************************/
//  Purpose:
//    E0000 is a reverse-communication zero bounder.
{

	const ips2pit2pul7pit * const s_data = data;
		
	int IENTRY = s_data->a0;
	short *status = s_data->a1;
	ityp *x = s_data->a2;
	ityp *fx = s_data->a3;
  	unsigned long *qleft = s_data->a4;
	unsigned long *qhi = s_data->a5;
	ityp *zabsst = s_data->a6;
  	ityp *zabsto = s_data->a7;
	ityp *zbig = s_data->a8;
	ityp *zrelst = s_data->a9;
  	ityp *zrelto = s_data->a10;
	ityp *zsmall = s_data->a11;
	ityp *zstpmu = s_data->a12;
	
    # define qxmon(zx,zy,zz) (int)((zx) <= (zy) && (zy) <= (zz))

    static ityp absstp;
    static ityp abstol;
    static ityp big,fbig,fsmall,relstp,reltol,small,step,stpmul,xhi,xlb,xlo,xsave,xub,yy;
    static int i99999;
    static unsigned long qbdd,qcond,qdum1,qdum2,qincr,qlim,qok,qup;
    switch(IENTRY)
    {
        case 0: goto DINVR;
        case 1: goto DSTINV;
    }
    DINVR:
    if(*status > 0) goto S310;
    qcond = !qxmon(small,*x,big);
    if(qcond)
        return NULL;
    xsave = *x;
    //     See that SMALL and BIG bound the zero and set QINCR
    *x = small;
    //     GET-FUNCTION-VALUE
    i99999 = 1;
    goto S300;
    S10:
        fsmall = *fx;
        *x = big;
        //     GET-FUNCTION-VALUE
        i99999 = 2;
        goto S300;
    S20:
        fbig = *fx;
        qincr = fbig > fsmall;
        if(!qincr) goto S50;
        if(fsmall <= 0.0e0) goto S30;
        *status = -1;
        *qleft = *qhi = 1;
        return NULL;
    S30:
        if(fbig >= 0.0e0) goto S40;
        *status = -1;
        *qleft = *qhi = 0;
        return NULL;
    S40:
        goto S80;
    S50:
        if(fsmall >= 0.0e0) goto S60;
        *status = -1;
        *qleft = 1;
        *qhi = 0;
        return NULL;
    S60:
        if(fbig <= 0.0e0) goto S70;
        *status = -1;
        *qleft = 0;
        *qhi = 1;
        return NULL;
    S80:
    S70:
        *x = xsave;
        step = fifdmax1(absstp,relstp*fabs(*x));
        //      YY = F(X) - Y
        //     GET-FUNCTION-VALUE
        i99999 = 3;
        goto S300;
    S90:
        yy = *fx;
        if(!(yy == 0.0e0)) goto S100;
        *status = 0;
        qok = 1;
        return NULL;
    S100:
        qup = qincr && yy < 0.0e0 || !qincr && yy > 0.0e0;
        //     HANDLE CASE IN WHICH WE MUST STEP HIGHER
        if(!qup) goto S170;
        xlb = xsave;
        xub = fifdmin1(xlb+step,big);
        goto S120;
    S110:
        if(qcond) goto S150;
    S120:
        //      YY = F(XUB) - Y
        *x = xub;
        //     GET-FUNCTION-VALUE
        i99999 = 4;
        goto S300;
    S130:
        yy = *fx;
        qbdd = qincr && yy >= 0.0e0 || !qincr && yy <= 0.0e0;
        qlim = xub >= big;
        qcond = qbdd || qlim;
        if(qcond) goto S140;
        step = stpmul*step;
        xlb = xub;
        xub = fifdmin1(xlb+step,big);
        S140:
        goto S110;
    S150:
        if(!(qlim && !qbdd)) goto S160;
        *status = -1;
        *qleft = 0;
        *qhi = !qincr;
        *x = big;
        return NULL;
    S160:
        goto S240;
    S170:
        //     HANDLE CASE IN WHICH WE MUST STEP LOWER
        xub = xsave;
        xlb = fifdmax1(xub-step,small);
        goto S190;
    S180:
        if(qcond) goto S220;
    S190:
        //      YY = F(XLB) - Y
        *x = xlb;
        //     GET-FUNCTION-VALUE
        i99999 = 5;
        goto S300;
    S200:
        yy = *fx;
        qbdd = qincr && yy <= 0.0e0 || !qincr && yy >= 0.0e0;
        qlim = xlb <= small;
        qcond = qbdd || qlim;
        if(qcond) goto S210;
        step = stpmul*step;
        xub = xlb;
        xlb = fifdmax1(xub-step,small);
    S210:
        goto S180;
    S220:
        if(!(qlim && !qbdd)) goto S230;
        *status = -1;
        *qleft = 1;
        *qhi = qincr;
        *x = small;
        return NULL;
    S240:
    S230:
        dstzr(&xlb,&xub,&abstol,&reltol);
    //  IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
        *status = 0;
        goto S260;
    S250:
        if(!(*status == 1)) goto S290;
    S260:
        dzror ( status, x, fx, &xlo, &xhi, &qdum1, &qdum2 );
        if(!(*status == 1)) goto S280;
        //     GET-FUNCTION-VALUE
        i99999 = 6;
        goto S300;
    S280:
    S270:
        goto S250;
        S290:
        *x = xlo;
        *status = 0;
        return NULL;
    DSTINV:
        small = *zsmall;
        big = *zbig;
        absstp = *zabsst;
        relstp = *zrelst;
        stpmul = *zstpmu;
        abstol = *zabsto;
        reltol = *zrelto;
        return NULL;
    S300:
        //     TO GET-FUNCTION-VALUE
        *status = 1;
        return NULL;
    S310:
        switch((int)i99999)
        {
            case 1: goto S10;
            case 2: goto S20;
            case 3: goto S90;
            case 4: goto S130;
            case 5: goto S200;
            case 6: goto S270;
            default: break;
        }
    # undef qxmon
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _E0001 ( void * data)
/******************************************************************************/
//  Purpose:
//    E00001 is a reverse-communication zero finder.
{
	const ips4pit2pul4pit * const s_data = data;
		
	int IENTRY = s_data->a0;
	short *status = s_data->a1;
	ityp *x = s_data->a2;
	ityp *fx = s_data->a3;
	ityp *xlo = s_data->a4;
	ityp *xhi = s_data->a5;
  	unsigned long *qleft = s_data->a6;
	unsigned long *qhi = s_data->a7;
	ityp *zabstl = s_data->a8;
  	ityp *zreltl = s_data->a9;
	ityp *zxhi = s_data->a10;
	ityp *zxlo = s_data->a11;
	
    # define ftol(zx) (0.5e0*fifdmax1(abstol,reltol*fabs((zx))))

    static ityp a,abstol,b,c,d,fa,fb,fc,fd,fda;
    static ityp fdb,m,mb,p,q,reltol,tol,w,xxhi,xxlo;
    static int ext,i99999;
    static unsigned long first,qrzero;
    switch(IENTRY)
    {
        case 0: goto DZROR;
        case 1: goto DSTZR;
    }
    DZROR:
        if(*status > 0) goto S280;
        *xlo = xxlo;
        *xhi = xxhi;
        b = *x = *xlo;
        //     GET-FUNCTION-VALUE
        i99999 = 1;
        goto S270;
    S10:
        fb = *fx;
        *xlo = *xhi;
        a = *x = *xlo;
        //     GET-FUNCTION-VALUE
        i99999 = 2;
        goto S270;
    S20:
        //     Check that F(ZXLO) < 0 < F(ZXHI)  or
        //                F(ZXLO) > 0 > F(ZXHI)
        if(!(fb < 0.0e0)) goto S40;
        if(!(*fx < 0.0e0)) goto S30;
        *status = -1;
        *qleft = *fx < fb;
        *qhi = 0;
        return NULL;
    S40:
    S30:
        if(!(fb > 0.0e0)) goto S60;
        if(!(*fx > 0.0e0)) goto S50;
        *status = -1;
        *qleft = *fx > fb;
        *qhi = 1;
        return NULL;
    S60:
    S50:
        fa = *fx;
        first = 1;
        S70:
        c = a;
        fc = fa;
        ext = 0;
    S80:
        if(!(fabs(fc) < fabs(fb))) goto S100;
        if(!(c != a)) goto S90;
        d = a;
        fd = fa;
    S90:
        a = b;
        fa = fb;
        *xlo = c;
        b = *xlo;
        fb = fc;
        c = a;
        fc = fa;
    S100:
        tol = ftol(*xlo);
        m = (c+b)*.5e0;
        mb = m-b;
        if(!(fabs(mb) > tol)) goto S240;
        if(!(ext > 3)) goto S110;
        w = mb;
        goto S190;
    S110:
        tol = fifdsign(tol,mb);
        p = (b-a)*fb;
        if(!first) goto S120;
        q = fa-fb;
        first = 0;
        goto S130;
    S120:
        fdb = (fd-fb)/(d-b);
        fda = (fd-fa)/(d-a);
        p = fda*p;
        q = fdb*fa-fda*fb;
    S130:
        if(!(p < 0.0e0)) goto S140;
        p = -p;
        q = -q;
    S140:
        if(ext == 3) p *= 2.0e0;
        if(!(p*1.0e0 == 0.0e0 || p <= q*tol)) goto S150;
        w = tol;
        goto S180;
    S150:
        if(!(p < mb*q)) goto S160;
        w = p/q;
        goto S170;
    S160:
        w = mb;
    S190:
    S180:
    S170:
        d = a;
        fd = fa;
        a = b;
        fa = fb;
        b += w;
        *xlo = b;
        *x = *xlo;
        //     GET-FUNCTION-VALUE
        i99999 = 3;
        goto S270;
    S200:
        fb = *fx;
        if(!(fc*fb >= 0.0e0)) goto S210;
        goto S70;
    S210:
        if(!(w == mb)) goto S220;
        ext = 0;
        goto S230;
    S220:
        ext += 1;
    S230:
        goto S80;
    S240:
        *xhi = c;
        qrzero = fc >= 0.0e0 && fb <= 0.0e0 || fc < 0.0e0 && fb >= 0.0e0;
        if(!qrzero) goto S250;
        *status = 0;
        goto S260;
    S250:
        *status = -1;
    S260:
        return NULL;
    DSTZR:
        xxlo = *zxlo;
        xxhi = *zxhi;
        abstol = *zabstl;
        reltol = *zreltl;
        return NULL;
    S270:
        //     TO GET-FUNCTION-VALUE
        *status = 1;
        return NULL;
    S280:
        switch((int)i99999)
        {
            case 1: goto S10;
            case 2: goto S20;
            case 3: goto S200;
            default: break;
        }
    # undef ftol
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _error_f ( void * data)
/******************************************************************************/
//  Purpose:
//    ERROR_F evaluates the error function ERF.
//  Parameters:
//    Input, double *X, the argument.
//    Output, double ERROR_F, the value of the error function at X.
{
	static ityp result = MAX_VAL;
	
	ityp * x = data;
	
    static ityp c = .564189583547756e0;
    static ityp a[5] =
    {
        .771058495001320e-04,
        -.133733772997339e-02,
        .323076579225834e-01,
        .479137145607681e-01,
        .128379167095513e+00
    };
    static ityp b[3] =
    {
        .301048631703895e-02,
        .538971687740286e-01,
        .375795757275549e+00
    };
    static ityp p[8] =
    {
        -1.36864857382717e-07,
        5.64195517478974e-01,
        7.21175825088309e+00,
        4.31622272220567e+01,
        1.52989285046940e+02,
        3.39320816734344e+02,
        4.51918953711873e+02,
        3.00459261020162e+02
    };
    static ityp q[8] =
    {
        1.00000000000000e+00,
        1.27827273196294e+01,
        7.70001529352295e+01,
        2.77585444743988e+02,
        6.38980264465631e+02,
        9.31354094850610e+02,
        7.90950925327898e+02,
        3.00459260956983e+02
    };
    static ityp r[5] =
    {
        2.10144126479064e+00,
        2.62370141675169e+01,
        2.13688200555087e+01,
        4.65807828718470e+00,
        2.82094791773523e-01
    };
    static ityp s[4] =
    {
        9.41537750555460e+01,
        1.87114811799590e+02,
        9.90191814623914e+01,
        1.80124575948747e+01
    };
    static ityp erf1,ax,bot,t,top,x2;

    ax = fabs(*x);
    if(ax > 0.5e0) goto S10;
    t = *x**x;
    top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0e0;
    bot = ((b[0]*t+b[1])*t+b[2])*t+1.0e0;
    erf1 = *x*(top/bot);
    result = erf1;
    return &result;
    S10:
        if(ax > 4.0e0) goto S20;
        top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[7];
        bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[7];
        erf1 = 0.5e0+(0.5e0-exp(-(*x**x))*top/bot);
        if(*x < 0.0e0) erf1 *= -1;
        result = erf1;
   	 return &result;
    S20:
        if(ax >= 5.8e0) goto S30;
        x2 = *x**x;
        t = 1.0e0/x2;
        top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
        bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0e0;
        erf1 = (c-top/(x2*bot))/ax;
        erf1 = 0.5e0+(0.5e0-exp(-x2)*erf1);
        if(*x < 0.0e0) erf1 = -erf1;
        result = erf1;
   		return &result;
    S30:
        erf1 = fifdsign(1.0e0,*x);
        
    result = erf1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _error_fc ( void * data) 
/******************************************************************************/
//  Purpose:
//    ERROR_FC evaluates the complementary error function ERFC.
//  Modified:
//    09 December 1999
//  Parameters:
//    Input, int *IND, chooses the scaling.
//    If IND is nonzero, then the value returned has been multiplied by
//    EXP(X*X).
//    Input, double *X, the argument of the function.
//    Output, double ERROR_FC, the value of the complementary
//    error function.
{
	static ityp result = MAX_VAL;
	
	const pdtpit * const s_data = data;
	dim_typ * ind = s_data->a0;
	ityp * x = s_data->a1; 
	
    static ityp c = .564189583547756e0;
    static ityp a[5] =
    {
        .771058495001320e-04,
        -.133733772997339e-02,
        .323076579225834e-01,
        .479137145607681e-01,
        .128379167095513e+00
    };
    static ityp b[3] =
    {
        .301048631703895e-02,
        .538971687740286e-01,
        .375795757275549e+00
    };
    static ityp p[8] =
    {
        -1.36864857382717e-07,
        5.64195517478974e-01,
        7.21175825088309e+00,
        4.31622272220567e+01,
        1.52989285046940e+02,
        3.39320816734344e+02,
        4.51918953711873e+02,
        3.00459261020162e+02
    };
    static ityp q[8] =
    {
        1.00000000000000e+00,
        1.27827273196294e+01,
        7.70001529352295e+01,
        2.77585444743988e+02,
        6.38980264465631e+02,
        9.31354094850610e+02,
        7.90950925327898e+02,
        3.00459260956983e+02
    };
    static ityp r[5] =
    {
        2.10144126479064e+00,
        2.62370141675169e+01,
        2.13688200555087e+01,
        4.65807828718470e+00,
        2.82094791773523e-01
    };
    static ityp s[4] =
    {
        9.41537750555460e+01,
        1.87114811799590e+02,
        9.90191814623914e+01,
        1.80124575948747e+01
    };
    static dim_typ K1 = 1;
    static ityp erfc1,ax,bot,e,t,top,w;

    //                     ABS(X) .LE. 0.5
    ax = fabs(*x);
    if(ax > 0.5e0) goto S10;
    t = *x**x;
    top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0e0;
    bot = ((b[0]*t+b[1])*t+b[2])*t+1.0e0;
    erfc1 = 0.5e0+(0.5e0-*x*(top/bot));
    if(*ind != 0) erfc1 = exp(t)*erfc1;
    result = erfc1;
    return &result;
    S10:
        //                  0.5 .LT. ABS(X) .LE. 4
        if(ax > 4.0e0) goto S20;
        top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[7];
        bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[7];
        erfc1 = top/bot;
        goto S40;
    S20:
        //                      ABS(X) .GT. 4
        if(*x <= -5.6e0) goto S60;
        if(*ind != 0) goto S30;
        if(*x > 100.0e0) goto S70;
        if(*x**x > -exparg(&K1)) goto S70;
    S30:
        t = pow(1.0e0/ *x,2.0);
        top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
        bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0e0;
        erfc1 = (c-t*top/bot)/ax;
    S40:
        //                      FINAL ASSEMBLY
        if(*ind == 0) goto S50;
        if(*x < 0.0e0) erfc1 = 2.0e0*exp(*x**x)-erfc1;
        result = erfc1;
    	return &result;
        S50:
        w = *x**x;
        t = w;
        e = w-t;
        erfc1 = (0.5e0+(0.5e0-e))*exp(-t)*erfc1;
        if(*x < 0.0e0) erfc1 = 2.0e0-erfc1;
        result = erfc1;
   	 	return &result;
    S60:
        //             LIMIT VALUE FOR LARGE NEGATIVE X
        erfc1 = 2.0e0;
        if(*ind != 0) erfc1 = 2.0e0*exp(*x**x);
        result = erfc1;
    	return &result;
    S70:
        //             LIMIT VALUE FOR LARGE POSITIVE X
        //                       WHEN IND = 0
        erfc1 = 0.0e0;
    
    result = erfc1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _esum ( void * data)
/******************************************************************************/
//  Purpose:
//    ESUM evaluates exp ( MU + X ).
//  Parameters:
//    Input, int *MU, part of the argument.
//    Input, double *X, part of the argument.
//    Output, double ESUM, the value of exp ( MU + X ).
{
	static ityp result = MAX_VAL;
	
	const pdtpit * const s_data = data;
	dim_typ * mu = s_data->a0;
	ityp * x = s_data->a1; 
	
    static ityp esum,w;

    if(*x > 0.0e0) goto S10;
    if(*mu < 0) goto S20;
    w = (ityp)*mu+*x;
    if(w > 0.0e0) goto S20;
    esum = exp(w);
    ityp value = esum;
    result = value;
    return &result;
    S10:
        if(*mu > 0) goto S20;
        w = (ityp)*mu+*x;
        if(w < 0.0e0) goto S20;
        esum = exp(w);
        value = esum;
        result = value;
   		return &result;
    S20:
        w = *mu;
        esum = exp(w)*exp(*x);
        
    value = esum;
    result = value;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _eval_pol ( void * data)
/******************************************************************************/
//  Purpose:
//    EVAL_POL evaluates a polynomial at X.
//  Discussion:
//    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
//  Modified:
//    15 December 1999
//  Parameters:
//    Input, double precision A(0:N), coefficients of the polynomial.
//    Input, int *N, length of A.
//    Input, double *X, the point at which the polynomial
//    is to be evaluated.
//    Output, double EVAL_POL, the value of the polynomial at X.
{
	static ityp result = MAX_VAL;
	
	const pdt2pit * const s_data = data;
	
	dim_typ * n = s_data->a0;
	ityp * a = s_data->a1;
	ityp * x = s_data->a2;
	
    static ityp devlpl,term;
    static dim_typ i;
    term = a[*n-1];
    for ( i = *n-1-1; i >= 0; --i )
        term = a[i]+term**x;
    devlpl = term;
    result = devlpl;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _exparg ( void * data)
/******************************************************************************/
//  Purpose:
//    EXPARG returns the largest or smallest legal argument for EXP.
//  Discussion:
//    Only an approximate limit for the argument of EXP is desired.
//  Modified:
//    09 December 1999
//  Parameters:
//    Input, int *L, indicates which limit is desired.
//    If L = 0, then the largest positive argument for EXP is desired.
//    Otherwise, the largest negative argument for EXP for which the
//    result is nonzero is desired.
//    Output, double EXPARG, the desired value.
{
	static ityp result = MAX_VAL;
	
	dim_typ * l = data;
	
    static int K1 = 4;
    static int K2 = 9;
    static int K3 = 10;
    static ityp exparg,lnb;
    static dim_typ b,m;

    b = ipmpar(&K1);
    if(b != 2) goto S10;
    lnb = .69314718055995e0;
    goto S40;
    S10:
        if(b != 8) goto S20;
        lnb = 2.0794415416798e0;
        goto S40;
    S20:
        if(b != 16) goto S30;
        lnb = 2.7725887222398e0;
        goto S40;
    S30:
        lnb = log((ityp)b);
    S40:
        if(*l == 0) goto S50;
        m = ipmpar(&K2)-1;
        exparg = 0.99999e0*((ityp)m*lnb);
        result = exparg;
        return &result;
    S50:
        m = ipmpar(&K3);
        exparg = 0.99999e0*((ityp)m*lnb);
        
    result = exparg;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *    _fifdint ( void * data)
/******************************************************************************/
//  Purpose:
//    FIFDINT truncates a double number to an integer.
//  Parameters:
// a     -     number to be truncated
{
	static ityp result = MAX_VAL;
	
	result = *(dim_typ *) data;
  	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fifdmax1 ( void * data)
/******************************************************************************/
//
//  Purpose:
//    FIFDMAX1 returns the maximum of two numbers a and b
//  Parameters:
//  a     -      first number
//  b     -      second number
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	ityp a = a_data[0];
	ityp b = a_data[1];
	
	result = a<b?b:a;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fifdmin1 ( void * data)
/******************************************************************************/
//
//  Purpose:
//    FIFDMIN1 returns the minimum of two numbers.
//  Parameters:
//  a     -     first number
//  b     -     second number
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	ityp a = a_data[0];
	ityp b = a_data[1];
	
    result = a<b?b:a;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _fifdsign (void * data)
/******************************************************************************/
//
//  Purpose:
//    FIFDSIGN transfers the sign of the variable "sign" to the variable "mag"
//  Parameters:
//  mag     -     magnitude
//  sign    -     sign to be transfered
{
	static ityp result = MAX_VAL;
	
	ityp * const a_data = data;
	ityp mag = a_data[0];
	ityp sign = a_data[1];
	
    if (mag < 0 || sign < 0) mag *= -1;
    
    result = mag;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fifidint ( void * data)
/******************************************************************************/
//
//  Purpose:
//    FIFIDINT truncates a double number to a long integer
//  Parameters:
//  a - number to be truncated
{
	static long result = LONG_MAX;
	
	ityp a = *(ityp *) data; 
	result = a<1.00 ? (long)0 : (long)a;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT   inline void *   _fifmod ( void * data)
/******************************************************************************/
//
//  Purpose:
//    FIFMOD returns the modulo of a and b
//  Parameters:
//  a - numerator
//  b - denominator
{
	static long result = LONG_MAX;
	
	const long * const a_data = data;
	long a = a_data[0];
	long b = a_data[1];
	
	result = ( a % b );
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _fpser ( void * data)
/******************************************************************************/
//
//  Purpose:
//    FPSER evaluates IX(A,B)(X) for very small B.
//  Discussion:
//    This routine is appropriate for use when
//      B < MIN ( EPS, EPS * A )
//    and
//      X <= 0.5.
//  Parameters:
//    Input, double *A, *B, parameters of the function.
//    Input, double *X, the point at which the function is to
//    be evaluated.
//    Input, double *EPS, a tolerance.
//    Output, double FPSER, the value of IX(A,B)(X).
{
	static ityp result = MAX_VAL;
	
  ityp ** const a_data = data;
  ityp * a = a_data[0];
  ityp * b = a_data[1];
  ityp * x = a_data[2];
  ityp * eps = a_data[3];
  
  	
  static dim_typ K1 = 1;
  static ityp fpser,an,c,s,t,tol;

    fpser = 1.0e0;
    if(*a <= 1.e-3**eps) goto S10;
    fpser = 0.0e0;
    t = *a*log(*x);
    if(t < exparg(&K1))
	{
		result = fpser;
		return &result;
	}
    fpser = exp(t);
    S10:
        //
        //                NOTE THAT 1/B(A,B) = B
        //
        fpser = *b/ *a*fpser;
        tol = *eps/ *a;
        an = *a+1.0e0;
        t = *x;
        s = t/an;
    S20:
        an += 1.0e0;
        t = *x*t;
        c = t/an;
        s += c;
        if(fabs(c) > tol) goto S20;
        fpser *= (1.0e0+*a*s);
        
    result = fpser;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gam1 ( void * data)
/******************************************************************************/
//  Purpose:
//    GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5D+00 <= A <= 1.5
//  Parameters:
//    Input, double *A, forms the argument of the Gamma function.
//    Output, double GAM1, the value of 1 / GAMMA ( A + 1 ) - 1.
{
	static ityp result = MAX_VAL;
	
	ityp * a = data;
	
    static ityp s1 = .273076135303957e+00;
    static ityp s2 = .559398236957378e-01;
    static ityp p[7] =
    {
        .577215664901533e+00,
        -.409078193005776e+00,
        -.230975380857675e+00,
        .597275330452234e-01,
        .766968181649490e-02,
        -.514889771323592e-02,
        .589597428611429e-03
    };
    static ityp q[5] =
    {
        .100000000000000e+01,
        .427569613095214e+00,
        .158451672430138e+00,
        .261132021441447e-01,
        .423244297896961e-02
    };
    static ityp r[9] =
    {
        -.422784335098468e+00,
        -.771330383816272e+00,
        -.244757765222226e+00,
        .118378989872749e+00,
        .930357293360349e-03,
        -.118290993445146e-01,
        .223047661158249e-02,
        .266505979058923e-03,-.132674909766242e-03
    };
    static ityp gam1,bot,d,t,top,w,T1;
    t = *a;
    d = *a-0.5e0;
    if(d > 0.0e0) t = d-0.5e0;
    T1 = t;
    if(T1 < 0) goto S40;
    else if(T1 == 0) goto S10;
    else  goto S20;
    S10:
        gam1 = 0.0e0;
        result = gam1;
        return &result;
    S20:
        top = (((((p[6]*t+p[5])*t+p[4])*t+p[3])*t+p[2])*t+p[1])*t+p[0];
        bot = (((q[4]*t+q[3])*t+q[2])*t+q[1])*t+1.0e0;
        w = top/bot;
        if(d > 0.0e0) goto S30;
        gam1 = *a*w;
        result = gam1;
        return &result;
    S30:
        gam1 = t/ *a*(w-0.5e0-0.5e0);
        result = gam1;
        return &result;
    S40:
        top = (((((((r[8]*t+r[7])*t+r[6])*t+r[5])*t+r[4])*t+r[3])*t+r[2])*t+r[1])*t+
        r[0];
        bot = (s2*t+s1)*t+1.0e0;
        w = top/bot;
        if(d > 0.0e0) goto S50;
        gam1 = *a*(w+0.5e0+0.5e0);
        result = gam1;
        return &result;
    S50:
        gam1 = t*w/ *a;
        
    result = gam1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void  * _gamma_inc_inv ( void * data)
/******************************************************************************/
//  Purpose:
//    GAMMA_INC_INV computes the inverse incomplete gamma ratio function.
//  Discussion:
//    The routine is given positive A, and nonnegative P and Q where P + Q = 1.
//    The value X is computed with the property that P(A,X) = P and Q(A,X) = Q.
//    Schroder iteration is employed.  The routine attempts to compute X
//    to 10 significant digits if this is possible for the particular computer
//    arithmetic being used.
//  Author:
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//  Parameters:
//    Input, double *A, the parameter in the incomplete gamma
//    ratio.  A must be positive.
//    Output, double *X, the computed point for which the
//    incomplete gamma functions have the values P and Q.
//    Input, double *X0, an optional initial approximation
//    for the solution X.  If the user does not want to supply an
//    initial approximation, then X0 should be set to 0, or a negative
//    value.
//    Input, double *P, *Q, the values of the incomplete gamma
//    functions, for which the corresponding argument is desired.
//    Output, int *IERR, error flag.
//    0, the solution was obtained. Iteration was not used.
//    0 < K, The solution was obtained. IERR iterations were performed.
//    -2, A <= 0
//    -3, No solution was obtained. The ratio Q/A is too large.
//    -4, P + Q /= 1
//    -6, 20 iterations were performed. The most recent value obtained
//        for X is given.  This cannot occur if X0 <= 0.
//    -7, Iteration failed. No value is given for X.
//        This may occur when X is approximately 0.
//    -8, A value for X has been obtained, but the routine is not certain
//        of its accuracy.  Iteration cannot be performed in this
//        case. If X0 <= 0, this can occur only when P or Q is
//        approximately 0. If X0 is positive then this can occur when A is
//        exceedingly close to X and A is extremely large (say A .GE. 1.E20).
{
	const _4pitpdtpit * const s_data = data;
	
	ityp *a = s_data->a0;
	ityp *x = s_data->a1;
	ityp *x0 = s_data->a2;
	ityp *p = s_data->a3;
	dim_typ *ierr = s_data->a4;
	ityp *q = s_data->a5;
	
    static ityp a0 = 3.31125922108741e0;
    static ityp a1 = 11.6616720288968e0;
    static ityp a2 = 4.28342155967104e0;
    static ityp a3 = .213623493715853e0;
    static ityp b1 = 6.61053765625462e0;
    static ityp b2 = 6.40691597760039e0;
    static ityp b3 = 1.27364489782223e0;
    static ityp b4 = .036117081018842e0;
    static ityp c = .577215664901533e0;
    static ityp ln10 = 2.302585e0;
    static ityp tol = 1.e-5;
    static ityp amin[2] =
    {
    500.0e0,100.0e0
    };
    static ityp bmin[2] =
    {
        1.e-28,
        1.e-13
    };
    static ityp dmin[2] =
    {
        1.e-06,
        1.e-04
    };
    static ityp emin[2] =
    {
        2.e-03,
        6.e-03
    };
    static ityp eps0[2] =
    {
        1.e-10,
        1.e-08
    };
    static dim_typ K1 = 1;
    static dim_typ K2 = 2;
    static dim_typ K3 = 3;
    static dim_typ K8 = 0;
    static ityp am1,amax,ap1,ap2,ap3,apn,b,c1,c2,c3,c4,c5,d,e,e2,eps,g,h,pn,qg,qn,r,rta,s,s2,sum,t,u,w,xmax,xmin,xn,y,z;
    static dim_typ iop;
    static ityp T4,T5,T6,T7,T9;

    //  E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS.
    //            E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0.
    //            XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE
    //            LARGEST POSITIVE NUMBER.
    e = dpmpar(&K1);
    xmin = dpmpar(&K2);
    xmax = dpmpar(&K3);
    *x = 0.0e0;
    if(*a <= 0.0e0) goto S300;
    t = *p+*q-1.e0;
    if(fabs(t) > e) goto S320;
    *ierr = 0;
    if(*p == 0.0e0) return NULL;
    if(*q == 0.0e0) goto S270;
    if(*a == 1.0e0) goto S280;
    e2 = 2.0e0*e;
    amax = 0.4e-10/(e*e);
    iop = 1;
    if(e > 1.e-10) iop = 2;
    eps = eps0[iop-1];
    xn = *x0;
    if(*x0 > 0.0e0) goto S160;
    //        SELECTION OF THE INITIAL APPROXIMATION XN OF X
    //                       WHEN A .LT. 1
    if(*a > 1.0e0) goto S80;
    T4 = *a+1.0e0;
    g = gamma_x(&T4);
    qg = *q*g;
    if(qg == 0.0e0) goto S360;
    b = qg/ *a;
    if(qg > 0.6e0**a) goto S40;
    if(*a >= 0.30e0 || b < 0.35e0) goto S10;
    t = exp(-(b+c));
    u = t*exp(t);
    xn = t*exp(u);
    goto S160;
    S10:
        if(b >= 0.45e0) goto S40;
        if(b == 0.0e0) goto S360;
        y = -log(b);
        s = 0.5e0+(0.5e0-*a);
        z = log(y);
        t = y-s*z;
        if(b < 0.15e0) goto S20;
        xn = y-s*log(t)-log(1.0e0+s/(t+1.0e0));
        goto S220;
    S20:
        if(b <= 0.01e0) goto S30;
        u = ((t+2.0e0*(3.0e0-*a))*t+(2.0e0-*a)*(3.0e0-*a))/((t+(5.0e0-*a))*t+2.0e0);
        xn = y-s*log(t)-log(u);
        goto S220;
    S30:
        c1 = -(s*z);
        c2 = -(s*(1.0e0+c1));
        c3 = s*((0.5e0*c1+(2.0e0-*a))*c1+(2.5e0-1.5e0**a));
        c4 = -(s*(((c1/3.0e0+(2.5e0-1.5e0**a))*c1+((*a-6.0e0)**a+7.0e0))*c1+((11.0e0**a-46.0)**a+47.0e0)/6.0e0));
        c5 = -(s*((((-(c1/4.0e0)+(11.0e0**a-17.0e0)/6.0e0)*c1+((-(3.0e0**a)+13.0e0)*
        *a-13.0e0))*c1+0.5e0*(((2.0e0**a-25.0e0)**a+72.0e0)**a-61.0e0))*c1+(((25.0e0**a-195.0e0)**a+477.0e0)**a-379.0e0)/12.0e0));
        xn = (((c5/y+c4)/y+c3)/y+c2)/y+c1+y;
        if(*a > 1.0e0) goto S220;
        if(b > bmin[iop-1]) goto S220;
        *x = xn;
        return NULL;
    S40:
        if(b**q > 1.e-8) goto S50;
        xn = exp(-(*q/ *a+c));
        goto S70;
    S50:
        if(*p <= 0.9e0) goto S60;
        T5 = -*q;
        xn = exp((alnrel(&T5)+ gamma_ln1 ( a ) ) / *a );
        goto S70;
    S60:
        xn = exp(log(*p*g)/ *a);
    S70:
        if(xn == 0.0e0) goto S310;
        t = 0.5e0+(0.5e0-xn/(*a+1.0e0));
        xn /= t;
        goto S160;
    S80:
        //        SELECTION OF THE INITIAL APPROXIMATION XN OF X
        //                       WHEN A .GT. 1
        if(*q <= 0.5e0) goto S90;
        w = log(*p);
        goto S100;
    S90:
        w = log(*q);
    S100:
        t = sqrt(-(2.0e0*w));
        s = t-(((a3*t+a2)*t+a1)*t+a0)/((((b4*t+b3)*t+b2)*t+b1)*t+1.0e0);
        if(*q > 0.5e0) s = -s;
        rta = sqrt(*a);
        s2 = s*s;
        xn = *a+s*rta+(s2-1.0e0)/3.0e0+s*(s2-7.0e0)/(36.0e0*rta)-((3.0e0*s2+7.0e0)*
        s2-16.0e0)/(810.0e0**a)+s*((9.0e0*s2+256.0e0)*s2-433.0e0)/(38880.0e0**a*
        rta);
        xn = fifdmax1(xn,0.0e0);
        if(*a < amin[iop-1]) goto S110;
        *x = xn;
        d = 0.5e0+(0.5e0-*x/ *a);
        if(fabs(d) <= dmin[iop-1]) return NULL;
    S110:
        if(*p <= 0.5e0) goto S130;
        if(xn < 3.0e0**a) goto S220;
        y = -(w+ gamma_log ( a ) );
        d = fifdmax1(2.0e0,*a*(*a-1.0e0));
        if(y < ln10*d) goto S120;
        s = 1.0e0-*a;
        z = log(y);
        goto S30;
    S120:
        t = *a-1.0e0;
        T6 = -(t/(xn+1.0e0));
        xn = y+t*log(xn)-alnrel(&T6);
        T7 = -(t/(xn+1.0e0));
        xn = y+t*log(xn)-alnrel(&T7);
        goto S220;
    S130:
        ap1 = *a+1.0e0;
        if(xn > 0.70e0*ap1) goto S170;
        w += gamma_log ( &ap1 );
        if(xn > 0.15e0*ap1) goto S140;
        ap2 = *a+2.0e0;
        ap3 = *a+3.0e0;
        *x = exp((w+*x)/ *a);
        *x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2)))/ *a);
        *x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2)))/ *a);
        *x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2*(1.0e0+*x/ap3))))/ *a);
        xn = *x;
        if(xn > 1.e-2*ap1) goto S140;
        if(xn <= emin[iop-1]*ap1) return NULL;
        goto S170;
    S140:
        apn = ap1;
        t = xn/apn;
        sum = 1.0e0+t;
    S150:
        apn += 1.0e0;
        t *= (xn/apn);
        sum += t;
        if(t > 1.e-4) goto S150;
        t = w-log(sum);
        xn = exp((xn+t)/ *a);
        xn *= (1.0e0-(*a*log(xn)-xn-t)/(*a-xn));
        goto S170;
    S160:
        //                 SCHRODER ITERATION USING P
        if(*p > 0.5e0) goto S220;
        S170:
        if(*p <= 1.e10*xmin) goto S350;
        am1 = *a-0.5e0-0.5e0;
    S180:
        if(*a <= amax) goto S190;
        d = 0.5e0+(0.5e0-xn/ *a);
        if(fabs(d) <= e2) goto S350;
    S190:
        if(*ierr >= 20) goto S330;
        *ierr += 1;
        gamma_inc ( a, &xn, &pn, &qn, &K8 );
        if(pn == 0.0e0 || qn == 0.0e0) goto S350;
        r = rcomp(a,&xn);
        if(r == 0.0e0) goto S350;
        t = (pn-*p)/r;
        w = 0.5e0*(am1-xn);
        if(fabs(t) <= 0.1e0 && fabs(w*t) <= 0.1e0) goto S200;
        *x = xn*(1.0e0-t);
        if(*x <= 0.0e0) goto S340;
        d = fabs(t);
        goto S210;
    S200:
        h = t*(1.0e0+w*t);
        *x = xn*(1.0e0-h);
        if(*x <= 0.0e0) goto S340;
        if(fabs(w) >= 1.0e0 && fabs(w)*t*t <= eps) return NULL;
        d = fabs(h);
    S210:
        xn = *x;
        if(d > tol) goto S180;
        if(d <= eps) return NULL;
        if(fabs(*p-pn) <= tol**p) return NULL;
        goto S180;
    S220:
        //                 SCHRODER ITERATION USING Q
        if(*q <= 1.e10*xmin) goto S350;
        am1 = *a-0.5e0-0.5e0;
    S230:
        if(*a <= amax) goto S240;
        d = 0.5e0+(0.5e0-xn/ *a);
        if(fabs(d) <= e2) goto S350;
    S240:
        if(*ierr >= 20) goto S330;
        *ierr += 1;
        gamma_inc ( a, &xn, &pn, &qn, &K8 );
        if(pn == 0.0e0 || qn == 0.0e0) goto S350;
        r = rcomp(a,&xn);
        if(r == 0.0e0) goto S350;
        t = (*q-qn)/r;
        w = 0.5e0*(am1-xn);
        if(fabs(t) <= 0.1e0 && fabs(w*t) <= 0.1e0) goto S250;
        *x = xn*(1.0e0-t);
        if(*x <= 0.0e0) goto S340;
        d = fabs(t);
        goto S260;
    S250:
        h = t*(1.0e0+w*t);
        *x = xn*(1.0e0-h);
        if(*x <= 0.0e0) goto S340;
        if(fabs(w) >= 1.0e0 && fabs(w)*t*t <= eps) return NULL;
        d = fabs(h);
    S260:
        xn = *x;
        if(d > tol) goto S230;
        if(d <= eps) return NULL;
        if(fabs(*q-qn) <= tol**q) return NULL;
        goto S230;
    S270:
        //                       SPECIAL CASES
        *x = xmax;
        return NULL;
    S280:
        if(*q < 0.9e0) goto S290;
        T9 = -*p;
        *x = -alnrel(&T9);
        return NULL;
    S290:
        *x = -log(*q);
        return NULL;
    S300:
        //                       ERROR RETURN
        *ierr = -2;
        return NULL;
    S310:
        *ierr = -3;
        return NULL;
    S320:
        *ierr = -4;
        return NULL;
    S330:
        *ierr = -6;
        return NULL;
    S340:
        *ierr = -7;
        return NULL;
    S350:
        *x = xn;
        *ierr = -8;
        return NULL;
    S360:
        *x = xmax;
        *ierr = -8;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gamma_ln1 ( void * data)
/******************************************************************************/
//  Purpose:
//    GAMMA_LN1 evaluates ln ( Gamma ( 1 + A ) ), for -0.2 <= A <= 1.25.
//  Parameters:
//    Input, double *A, defines the argument of the function.
//    Output, double GAMMA_LN1, the value of ln ( Gamma ( 1 + A ) ).
{
	static ityp result = MAX_VAL;
	
	ityp * a = data;
	
    static ityp p0 = .577215664901533e+00;
    static ityp p1 = .844203922187225e+00;
    static ityp p2 = -.168860593646662e+00;
    static ityp p3 = -.780427615533591e+00;
    static ityp p4 = -.402055799310489e+00;
    static ityp p5 = -.673562214325671e-01;
    static ityp p6 = -.271935708322958e-02;
    static ityp q1 = .288743195473681e+01;
    static ityp q2 = .312755088914843e+01;
    static ityp q3 = .156875193295039e+01;
    static ityp q4 = .361951990101499e+00;
    static ityp q5 = .325038868253937e-01;
    static ityp q6 = .667465618796164e-03;
    static ityp r0 = .422784335098467e+00;
    static ityp r1 = .848044614534529e+00;
    static ityp r2 = .565221050691933e+00;
    static ityp r3 = .156513060486551e+00;
    static ityp r8 = .170502484022650e-01;
    static ityp r5 = .497958207639485e-03;
    static ityp s1 = .124313399877507e+01;
    static ityp s2 = .548042109832463e+00;
    static ityp s3 = .101552187439830e+00;
    static ityp s4 = .713309612391000e-02;
    static ityp s5 = .116165475989616e-03;
    static ityp gamln1,w,x;

    if(*a >= 0.6e0) goto S10;
    w = ((((((p6**a+p5)**a+p4)**a+p3)**a+p2)**a+p1)**a+p0)/((((((q6**a+q5)**a+
    q4)**a+q3)**a+q2)**a+q1)**a+1.0e0);
    gamln1 = -(*a*w);
    result = gamln1;
    return &result;
    S10:
        x = *a-0.5e0-0.5e0;
        w = (((((r5*x+r8)*x+r3)*x+r2)*x+r1)*x+r0)/(((((s5*x+s4)*x+s3)*x+s2)*x+s1)*x
        +1.0e0);
        gamln1 = x*w;
        
    result = gamln1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gamma_log ( void * data)
/******************************************************************************/
//  Purpose:
//    GAMMA_LOG evaluates ln ( Gamma ( A ) ) for positive A.
//  Author:
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//  Reference:
//    Armido DiDinato and Alfred Morris,
//    Algorithm 708:
//    Significant Digit Computation of the Incomplete Beta Function Ratios,
//    ACM Transactions on Mathematical Software,
//    Volume 18, 1993, pages 360-373.
//  Parameters:
//    Input, double *A, the argument of the function.
//    A should be positive.
//    Output, double GAMMA_LOG, the value of ln ( Gamma ( A ) ).
{
	static ityp result = MAX_VAL;
	
	ityp * a = data;
	
    static ityp c0 = .833333333333333e-01;
    static ityp c1 = -.277777777760991e-02;
    static ityp c2 = .793650666825390e-03;
    static ityp c3 = -.595202931351870e-03;
    static ityp c4 = .837308034031215e-03;
    static ityp c5 = -.165322962780713e-02;
    static ityp d = .418938533204673e0;
    static ityp gamln,t,w;
    static dim_typ i,n;
    static ityp T1;

    if(*a > 0.8e0) goto S10;
    gamln = gamma_ln1 ( a ) - log ( *a );
    result = gamln;
    return &result;
    S10:
        if(*a > 2.25e0) goto S20;
        t = *a-0.5e0-0.5e0;
        gamln = gamma_ln1 ( &t );
        result = gamln;
    	return &result;
    S20:
        if(*a >= 10.0e0) goto S40;
        n = ( int ) ( *a - 1.25e0 );
        t = *a;
        w = 1.0e0;
        for ( i = 1; i <= n; ++i )
        {
            t -= 1.0e0;
            w = t*w;
        }
        T1 = t-1.0e0;
        gamln = gamma_ln1 ( &T1 ) + log ( w );
        result = gamln;
    	return &result;
    S40:
        t = pow(1.0e0/ *a,2.0);
        w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/ *a;
        gamln = d+w+(*a-0.5e0)*(log(*a)-1.0e0);
        
    result = gamln;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gamma_rat1 ( void * data)
/******************************************************************************/
//  Purpose:
//    GAMMA_RAT1 evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
//  Parameters:
//    Input, double *A, *X, the parameters of the functions.
//    It is assumed that A <= 1.
//    Input, double *R, the value exp(-X) * X**A / Gamma(A).
//    Output, double *P, *Q, the values of P(A,X) and Q(A,X).
//    Input, double *EPS, the tolerance.
{
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * x = a_data[1];
	ityp * r = a_data[2];
	ityp * p = a_data[3];
	ityp * q = a_data[4];
	ityp * eps = a_data[5];
	
    static dim_typ K2 = 0;
    static ityp a2n,a2nm1,am0,an,an0,b2n,b2nm1,c,cma,g,h,j,l,sum,t,tol,w,z,T1,T3;

    if(*a**x == 0.0e0) goto S120;
    if(*a == 0.5e0) goto S100;
    if(*x < 1.1e0) goto S10;
    goto S60;
    S10:
        //             TAYLOR SERIES FOR P(A,X)/X**A
        an = 3.0e0;
        c = *x;
        sum = *x/(*a+3.0e0);
        tol = 0.1e0**eps/(*a+1.0e0);
    S20:
        an += 1.0e0;
        c = -(c*(*x/an));
        t = c/(*a+an);
        sum += t;
        if(fabs(t) > tol) goto S20;
        j = *a**x*((sum/6.0e0-0.5e0/(*a+2.0e0))**x+1.0e0/(*a+1.0e0));
        z = *a*log(*x);
        h = gam1(a);
        g = 1.0e0+h;
        if(*x < 0.25e0) goto S30;
        if(*a < *x/2.59e0) goto S50;
        goto S40;
    S30:
        if(z > -.13394e0) goto S50;
    S40:
        w = exp(z);
        *p = w*g*(0.5e0+(0.5e0-j));
        *q = 0.5e0+(0.5e0-*p);
        return NULL;
    S50:
        l = rexp(&z);
        w = 0.5e0+(0.5e0+l);
        *q = (w*j-l)*g-h;
        if(*q < 0.0e0) goto S90;
        *p = 0.5e0+(0.5e0-*q);
        return NULL;
    S60:
        //              CONTINUED FRACTION EXPANSION
        a2nm1 = a2n = 1.0e0;
        b2nm1 = *x;
        b2n = *x+(1.0e0-*a);
        c = 1.0e0;
    S70:
        a2nm1 = *x*a2n+c*a2nm1;
        b2nm1 = *x*b2n+c*b2nm1;
        am0 = a2nm1/b2nm1;
        c += 1.0e0;
        cma = c-*a;
        a2n = a2nm1+cma*a2n;
        b2n = b2nm1+cma*b2n;
        an0 = a2n/b2n;
        if(fabs(an0-am0) >= *eps*an0) goto S70;
        *q = *r*an0;
        *p = 0.5e0+(0.5e0-*q);
        return NULL;
    S80:
        //                SPECIAL CASES
        *p = 0.0e0;
        *q = 1.0e0;
        return NULL;
        S90:
        *p = 1.0e0;
        *q = 0.0e0;
        return NULL;
    S100:
        if(*x >= 0.25e0) goto S110;
        T1 = sqrt(*x);
        *p = error_f ( &T1 );
        *q = 0.5e0+(0.5e0-*p);
        return NULL;
    S110:
        T3 = sqrt(*x);
        *q = error_fc ( &K2, &T3 );
        *p = 0.5e0+(0.5e0-*q);
        return NULL;
    S120:
        if(*x <= *a) goto S80;
        goto S90;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gamma_x ( void * data)
/******************************************************************************/
//  Purpose:
//    GAMMA_X evaluates the gamma function.
//  Discussion:
//    This routine was renamed from "GAMMA" to avoid a conflict with the
//    C/C++ math library routine.
//  Author:
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//  Parameters:
//    Input, double *A, the argument of the Gamma function.
//    Output, double GAMMA_X, the value of the Gamma function.
{
	static ityp result = MAX_VAL;
	
	ityp * a = data;
	
    static ityp d = .41893853320467274178e0;
    static ityp r1 = .820756370353826e-03;
    static ityp r2 = -.595156336428591e-03;
    static ityp r3 = .793650663183693e-03;
    static ityp r8 = -.277777777770481e-02;
    static ityp r5 = .833333333333333e-01;
    static ityp p[7] =
    {
        .539637273585445e-03,
        .261939260042690e-02,
        .204493667594920e-01,
        .730981088720487e-01,
        .279648642639792e+00,
        .553413866010467e+00,
        1.0e0
    };
    static ityp q[7] =
    {
        -.832979206704073e-03,
        .470059485860584e-02,
        .225211131035340e-01,
        -.170458969313360e+00,
        -.567902761974940e-01,
        .113062953091122e+01,
        1.0e0
    };
    static dim_typ K2 = 3;
    static dim_typ K3 = 0;
    static ityp Xgamm,bot,g,lnx,s,t,top,w,x,z;
    static dim_typ i,j,m,n,T1;

    Xgamm = 0.0e0;
    x = *a;
    if(fabs(*a) >= 15.0e0) goto S110;
    //            EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15
    t = 1.0e0;
    m = fifidint(*a)-1;
    //     LET T BE THE PRODUCT OF A-J WHEN A .GE. 2
    T1 = m;
    if(T1 < 0) goto S40;
    else if(T1 == 0) goto S30;
    else  goto S10;
    S10:
        for ( j = 1; j <= m; ++j )
        {
            x -= 1.0e0;
            t = x*t;
        }
    S30:
        x -= 1.0e0;
        goto S80;
    S40:
        //     LET T BE THE PRODUCT OF A+J WHEN A .LT. 1
        t = *a;
        if(*a > 0.0e0) goto S70;
        m = -m-1;
        if(m == 0) goto S60;
        for ( j = 1; j <= m; ++j )
        {
            x += 1.0e0;
            t = x*t;
        }
    S60:
        x += (0.5e0+0.5e0);
        t = x*t;
        if(t == 0.0e0)
		{
			result = Xgamm;
			return &result;
		}
    S70:
        //     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
        //     CODE MAY BE OMITTED IF DESIRED.
        if(fabs(t) >= 1.e-30) goto S80;
        if(fabs(t)*dpmpar(&K2) <= 1.0001e0)
        {
        	result = Xgamm;
			return &result;
        }
        Xgamm = 1.0e0/t;
        result = Xgamm;
		return &result;
    S80:
        //     COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1
        top = p[0];
        bot = q[0];
        #pragma omp parallel for num_threads(6)
        for ( i = 1; i < 7; ++i )
        {
            top = p[i]+x*top;
            bot = q[i]+x*bot;
        }
        Xgamm = top/bot;
        //     TERMINATION
        if(*a < 1.0e0) goto S100;
        Xgamm *= t;
        result = Xgamm;
		return &result;
    S100:
        Xgamm /= t;
        result = Xgamm;
		return &result;
    S110:
        //  EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15
        if(fabs(*a) >= 1.e3)
        {
        	result = Xgamm;
			return &result;
        }
        if(*a > 0.0e0) goto S120;
        x = -*a;
        n = ( dim_typ ) x;
        t = x-(ityp)n;
        if(t > 0.9e0) t = 1.0e0-t;
        s = sin(M_PI*t)/M_PI;
        if(fifmod(n,2) == 0) s = -s;
        if(s == 0.0e0)
        {
        	result = Xgamm;
			return &result;
        }
    S120:
        //     COMPUTE THE MODIFIED ASYMPTOTIC SUM
        t = 1.0e0/(x*x);
        g = ((((r1*t+r2)*t+r3)*t+r8)*t+r5)/x;
        //     ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X)
        //     BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.
        lnx = log(x);
        //  FINAL ASSEMBLY
        z = x;
        g = d+g+(z-0.5e0)*(lnx-1.e0);
        w = g;
        t = g-w;
        if(w > 0.99999e0*exparg(&K3))
        {
        	result = Xgamm;
			return &result;
        }
        Xgamm = exp(w)*(1.0e0+t);
        if(*a < 0.0e0) Xgamm = 1.0e0/(Xgamm*s)/x;
        
    result = Xgamm;
	return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _gsumln ( void * data)
/******************************************************************************/
//  Purpose:
//    GSUMLN evaluates the function ln(Gamma(A + B)).
//  Discussion:
//    GSUMLN is used for 1 <= A <= 2 and 1 <= B <= 2
//  Parameters:
//    Input, double *A, *B, values whose sum is the argument of
//    the Gamma function.
//    Output, double GSUMLN, the value of ln(Gamma(A+B)).
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * b = a_data[1];
	
    static ityp gsumln,x,T1,T2;

    x = *a+*b-2.e0;
    if(x > 0.25e0) goto S10;
    T1 = 1.0e0+x;
    gsumln = gamma_ln1 ( &T1 );
    result = gsumln;
    return &result;
    S10:
        if(x > 1.25e0) goto S20;
        gsumln = gamma_ln1 ( &x ) + alnrel ( &x );
        result = gsumln;
    	return &result;
    S20:
        T2 = x-1.0e0;
        gsumln = gamma_ln1 ( &T2 ) + log ( x * ( 1.0e0 + x ) );
        
    result = gsumln;
    return &result;
}

/// MACHINE DEPENDENT CODE!!! MDC WARNING
/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _ipmpar ( void * data)
/******************************************************************************/
//  Purpose:
//    IPMPAR returns integer machine constants.
//  Discussion:
//    Input arguments 1 through 3 are queries about integer arithmetic.
//    We assume integers are represented in the N-digit, base-A form
//      sign * ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
//    where 0 <= X(0:N-1) < A.
//    Then:
//      IPMPAR(1) = A, the base of integer arithmetic;
//      IPMPAR(2) = N, the number of base A digits;
//      IPMPAR(3) = A**N - 1, the largest magnitude.
//    It is assumed that the single and double precision floating
//    point arithmetics have the same base, say B, and that the
//    nonzero numbers are represented in the form
//      sign * (B**E) * (X(1)/B + ... + X(M)/B**M)
//    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
//    EMIN <= E <= EMAX.
//    Input argument 4 is a query about the base of real arithmetic:
//      IPMPAR(4) = B, the base of single and double precision arithmetic.
//    Input arguments 5 through 7 are queries about single precision
//    floating point arithmetic:
//     IPMPAR(5) = M, the number of base B digits for single precision.
//     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
//     IPMPAR(7) = EMAX, the largest exponent E for single precision.
//    Input arguments 8 through 10 are queries about double precision
//    floating point arithmetic:
//     IPMPAR(8) = M, the number of base B digits for double precision.
//     IPMPAR(9) = EMIN, the smallest exponent E for double precision.
//     IPMPAR(10) = EMAX, the largest exponent E for double precision.
//  Reference:
//    Phyllis Fox, Andrew Hall, and Norman Schryer,
//    Algorithm 528,
//    Framework for a Portable FORTRAN Subroutine Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, 1978, pages 176-188.
//  Parameters:
//    Input, int *I, the index of the desired constant.
//    Output, int IPMPAR, the value of the desired constant.
{
	static int result = INT_MAX;
	
	int * i = data;
	
    static int imach[11];
    static int ipmpar;
    //     MACHINE CONSTANTS FOR AMDAHL MACHINES.
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 16;
    //   imach[5] = 6;
    //   imach[6] = -64;
    //   imach[7] = 63;
    //   imach[8] = 14;
    //   imach[9] = -64;
    //   imach[10] = 63;
    //     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
    //       PC 7300, AND AT&T 6300.
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 2;
    //   imach[5] = 24;
    //   imach[6] = -125;
    //   imach[7] = 128;
    //   imach[8] = 53;
    //   imach[9] = -1021;
    //   imach[10] = 1024;
    //     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
    //   imach[1] = 2;
    //   imach[2] = 33;
    //   imach[3] = 8589934591;
    //   imach[4] = 2;
    //   imach[5] = 24;
    //   imach[6] = -256;
    //   imach[7] = 255;
    //   imach[8] = 60;
    //   imach[9] = -256;
    //   imach[10] = 255;
    //     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
    //   imach[1] = 2;
    //   imach[2] = 39;
    //   imach[3] = 549755813887;
    //   imach[4] = 8;
    //   imach[5] = 13;
    //   imach[6] = -50;
    //   imach[7] = 76;
    //   imach[8] = 26;
    //   imach[9] = -50;
    //   imach[10] = 76;
    //     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
    //   imach[1] = 2;
    //   imach[2] = 39;
    //   imach[3] = 549755813887;
    //   imach[4] = 8;
    //   imach[5] = 13;
    //   imach[6] = -50;
    //   imach[7] = 76;
    //   imach[8] = 26;
    //   imach[9] = -32754;
    //   imach[10] = 32780;
    //     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
    //       60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
    //       ARITHMETIC (NOS OPERATING SYSTEM).
    //   imach[1] = 2;
    //   imach[2] = 48;
    //   imach[3] = 281474976710655;
    //   imach[4] = 2;
    //   imach[5] = 48;
    //   imach[6] = -974;
    //   imach[7] = 1070;
    //   imach[8] = 95;
    //   imach[9] = -926;
    //   imach[10] = 1070;
    //     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
    //       ARITHMETIC (NOS/VE OPERATING SYSTEM).
    //   imach[1] = 2;
    //   imach[2] = 63;
    //   imach[3] = 9223372036854775807;
    //   imach[4] = 2;
    //   imach[5] = 48;
    //   imach[6] = -4096;
    //   imach[7] = 4095;
    //   imach[8] = 96;
    //   imach[9] = -4096;
    //   imach[10] = 4095;
    //     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
    //   imach[1] = 2;
    //   imach[2] = 63;
    //   imach[3] = 9223372036854775807;
    //   imach[4] = 2;
    //   imach[5] = 47;
    //   imach[6] = -8189;
    //   imach[7] = 8190;
    //   imach[8] = 94;
    //   imach[9] = -8099;
    //   imach[10] = 8190;
    //     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
    //   imach[1] = 2;
    //   imach[2] = 15;
    //   imach[3] = 32767;
    //   imach[4] = 16;
    //   imach[5] = 6;
    //   imach[6] = -64;
    //   imach[7] = 63;
    //   imach[8] = 14;
    //   imach[9] = -64;
    //   imach[10] = 63;
    //     MACHINE CONSTANTS FOR THE HARRIS 220.
    //   imach[1] = 2;
    //   imach[2] = 23;
    //   imach[3] = 8388607;
    //   imach[4] = 2;
    //   imach[5] = 23;
    //   imach[6] = -127;
    //   imach[7] = 127;
    //   imach[8] = 38;
    //   imach[9] = -127;
    //   imach[10] = 127;
    //     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000
    //       AND DPS 8/70 SERIES.
    //
    //   imach[1] = 2;
    //   imach[2] = 35;
    //   imach[3] = 34359738367;
    //   imach[4] = 2;
    //   imach[5] = 27;
    //   imach[6] = -127;
    //   imach[7] = 127;
    //   imach[8] = 63;
    //   imach[9] = -127;
    //   imach[10] = 127;
    //     MACHINE CONSTANTS FOR THE HP 2100
    //       3 WORD DOUBLE PRECISION OPTION WITH FTN4
    //   imach[1] = 2;
    //   imach[2] = 15;
    //   imach[3] = 32767;
    //   imach[4] = 2;
    //   imach[5] = 23;
    //   imach[6] = -128;
    //   imach[7] = 127;
    //   imach[8] = 39;
    //   imach[9] = -128;
    //   imach[10] = 127;
    //     MACHINE CONSTANTS FOR THE HP 2100
    //       4 WORD DOUBLE PRECISION OPTION WITH FTN4
    //   imach[1] = 2;
    //   imach[2] = 15;
    //   imach[3] = 32767;
    //   imach[4] = 2;
    //   imach[5] = 23;
    //   imach[6] = -128;
    //   imach[7] = 127;
    //   imach[8] = 55;
    //   imach[9] = -128;
    //   imach[10] = 127;
    //     MACHINE CONSTANTS FOR THE HP 9000.
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 2;
    //   imach[5] = 24;
    //   imach[6] = -126;
    //   imach[7] = 128;
    //   imach[8] = 53;
    //   imach[9] = -1021;
    //   imach[10] = 1024;
    //     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
    //       THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
    //       5/7/9 AND THE SEL SYSTEMS 85/86.
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 16;
    //   imach[5] = 6;
    //   imach[6] = -64;
    //   imach[7] = 63;
    //   imach[8] = 14;
    //   imach[9] = -64;
    //   imach[10] = 63;
    //     MACHINE CONSTANTS FOR THE IBM PC.
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 2;
    //   imach[5] = 24;
    //   imach[6] = -125;
    //   imach[7] = 128;
    //   imach[8] = 53;
    //   imach[9] = -1021;
    //   imach[10] = 1024;
    //     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
    //       MACFORTRAN II.
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 2;
    //   imach[5] = 24;
    //   imach[6] = -125;
    //   imach[7] = 128;
    //   imach[8] = 53;
    //   imach[9] = -1021;
    //   imach[10] = 1024;
    //     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN.
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 2;
    //   imach[5] = 24;
    //   imach[6] = -127;
    //   imach[7] = 127;
    //   imach[8] = 56;
    //   imach[9] = -127;
    //   imach[10] = 127;
    //     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
    //   imach[1] = 2;
    //   imach[2] = 35;
    //   imach[3] = 34359738367;
    //   imach[4] = 2;
    //   imach[5] = 27;
    //   imach[6] = -128;
    //   imach[7] = 127;
    //   imach[8] = 54;
    //   imach[9] = -101;
    //   imach[10] = 127;
    //     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
    //   imach[1] = 2;
    //   imach[2] = 35;
    //   imach[3] = 34359738367;
    //   imach[4] = 2;
    //   imach[5] = 27;
    //   imach[6] = -128;
    //   imach[7] = 127;
    //   imach[8] = 62;
    //   imach[9] = -128;
    //   imach[10] = 127;
    //     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
    //       32-BIT INTEGER ARITHMETIC.
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 2;
    //   imach[5] = 24;
    //   imach[6] = -127;
    //   imach[7] = 127;
    //   imach[8] = 56;
    //   imach[9] = -127;
    //   imach[10] = 127;
    //     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 2;
    //   imach[5] = 24;
    //   imach[6] = -125;
    //   imach[7] = 128;
    //   imach[8] = 53;
    //   imach[9] = -1021;
    //   imach[10] = 1024;
    //     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
    //       SERIES (MIPS R3000 PROCESSOR).
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 2;
    //   imach[5] = 24;
    //   imach[6] = -125;
    //   imach[7] = 128;
    //   imach[8] = 53;
    //   imach[9] = -1021;
    //   imach[10] = 1024;
    //     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
    //       3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
    //       PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).

    imach[1] = 2;
    imach[2] = 31;
    imach[3] = 2147483647;
    imach[4] = 2;
    imach[5] = 24;
    imach[6] = -125;
    imach[7] = 128;
    imach[8] = 53;
    imach[9] = -1021;
    imach[10] = 1024;

    //     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
    //   imach[1] = 2;
    //   imach[2] = 35;
    //   imach[3] = 34359738367;
    //   imach[4] = 2;
    //   imach[5] = 27;
    //   imach[6] = -128;
    //   imach[7] = 127;
    //   imach[8] = 60;
    //   imach[9] = -1024;
    //   imach[10] = 1023;

    //     MACHINE CONSTANTS FOR THE VAX 11/780.
    //   imach[1] = 2;
    //   imach[2] = 31;
    //   imach[3] = 2147483647;
    //   imach[4] = 2;
    //   imach[5] = 24;
    //   imach[6] = -127;
    //   imach[7] = 127;
    //   imach[8] = 56;
    //   imach[9] = -127;
    //   imach[10] = 127;
    //
    ipmpar = imach[*i];
    result = ipmpar;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _psi ( void * data)
/******************************************************************************/
//  Purpose:
//    PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
//  Discussion:
//    The main computation involves evaluation of rational Chebyshev
//    approximations.  PSI was written at Argonne National Laboratory
//    for FUNPACK, and subsequently modified by A. H. Morris of NSWC.
//  Reference:
//    William Cody, Strecok and Thacher,
//    Chebyshev Approximations for the Psi Function,
//    Mathematics of Computation,
//    Volume 27, 1973, pages 123-127.
//  Parameters:
//    Input, double *XX, the argument of the psi function.
//    Output, double PSI, the value of the psi function.  PSI
//    is assigned the value 0 when the psi function is undefined.
{
	static ityp result = MAX_VAL;
	
	ityp * xx = data; 
	
    static ityp dx0 = 1.461632144968362341262659542325721325e0;
    static ityp piov4 = .785398163397448e0;
    static ityp p1[7] =
    {
        .895385022981970e-02,
        .477762828042627e+01,
        .142441585084029e+03,
        .118645200713425e+04,
        .363351846806499e+04,
        .413810161269013e+04,
        .130560269827897e+04
    };
    static ityp p2[4] =
    {
        -.212940445131011e+01,
        -.701677227766759e+01,
        -.448616543918019e+01,
        -.648157123766197e+00
    };
    static ityp q1[6] =
    {
        .448452573429826e+02,
        .520752771467162e+03,
        .221000799247830e+04,
        .364127349079381e+04,
        .190831076596300e+04,
        .691091682714533e-05
    };
    static ityp q2[4] =
    {
        .322703493791143e+02,
        .892920700481861e+02,
        .546117738103215e+02,
        .777788548522962e+01
    };
    static int K1 = 3;
    static dim_typ K2 = 1;
    static ityp psi,aug,den,sgn,upper,w,x,xmax1,xmx0,xsmall,z;
    static dim_typ i,m,n,nq;
    //     MACHINE DEPENDENT CONSTANTS ...
    //        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
    //                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
    //                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
    //                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
    //                 PSI MAY BE REPRESENTED AS ALOG(X).
    //        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH M_PI*COTAN(M_PI*X)
    //                 MAY BE REPRESENTED BY 1/X.
    xmax1 = ipmpar(&K1);
    xmax1 = fifdmin1(xmax1,1.0e0/dpmpar(&K2));
    xsmall = 1.e-9;
    x = *xx;
    aug = 0.0e0;
    if(x >= 0.5e0) goto S50;
    //     X .LT. 0.5,  USE REFLECTION FORMULA
    //     PSI(1-X) = PSI(X) + M_PI * COTAN(M_PI*X)
    if(fabs(x) > xsmall) goto S10;
    if(x == 0.0e0) goto S100;
    //     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
    //     FOR  M_PI*COTAN(M_PI*X)
    aug = -(1.0e0/x);
    goto S40;
    S10:
        //     REDUCTION OF ARGUMENT FOR COTAN
        w = -x;
        sgn = piov4;
        if(w > 0.0e0) goto S20;
        w *= -1;
        sgn *= -1;
    S20:
        //     MAKE AN ERROR EXIT IF X .LE. -XMAX1
        if(w >= xmax1) goto S100;
        nq = fifidint(w);
        w -= (ityp)nq;
        nq = fifidint(w*4.0e0);
        w = 4.0e0*(w-(ityp)nq*.25e0);
        //     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
        //     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
        //     QUADRANT AND DETERMINE SIGN
        n = nq/2;
        if(n+n != nq) w = 1.0e0-w;
        z = piov4*w;
        m = n/2;
        if(m+m != n) sgn = -sgn;
        //     DETERMINE FINAL VALUE FOR  -M_PI*COTAN(M_PI*X)
        n = (nq+1)/2;
        m = n/2;
        m += m;
        if(m != n) goto S30;
        //     CHECK FOR SINGULARITY
        if(z == 0.0e0) goto S100;
        //     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
        //     SIN/COS AS A SUBSTITUTE FOR TAN
        aug = sgn*(cos(z)/sin(z)*4.0e0);
        goto S40;
    S30:
        aug = sgn*(sin(z)/cos(z)*4.0e0);
    S40:
        x = 1.0e0-x;
    S50:
        if(x > 3.0e0) goto S70;
        //     0.5 .LE. X .LE. 3.0
        den = x;
        upper = p1[0]*x;
        #pragma omp parallel for num_threads(6)
        for ( i = 1; i <= 5; ++i )
        {
            den = (den+q1[i-1])*x;
            upper = (upper+p1[i+1-1])*x;
        }
        den = (upper+p1[6])/(den+q1[5]);
        xmx0 = x-dx0;
        psi = den*xmx0+aug;
        result = psi;
        return &result;
    S70:
        //     IF X .GE. XMAX1, PSI = LN(X)
        if(x >= xmax1) goto S90;
        //     3.0 .LT. X .LT. XMAX1
        w = 1.0e0/(x*x);
        den = w;
        upper = p2[0]*w;
        #pragma omp parallel for num_threads(3)
        for ( i = 1; i <= 3; ++i )
        {
            den = (den+q2[i-1])*w;
            upper = (upper+p2[i+1-1])*w;
        }
        aug = upper/(den+q2[3])-0.5e0/x+aug;
    S90:
        psi = aug+log(x);
        result = psi;
        return &result;
    S100:
        //     ERROR RETURN
        psi = 0.0e0;
        
    result = psi;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rcomp ( void * data)
/******************************************************************************/
//  Purpose:
//    RCOMP evaluates exp(-X) * X**A / Gamma(A).
//  Parameters:
//    Input, double *A, *X, arguments of the quantity to be computed.
//    Output, double RCOMP, the value of exp(-X) * X**A / Gamma(A).
//  Local parameters:
//    RT2PIN = 1/SQRT(2*M_PI)
{
	static ityp result = MAX_VAL;
	
	ityp ** const a_data = data;
	ityp * a = a_data[0];
	ityp * x = a_data[1];
	
    static ityp rt2pin = .398942280401433e0;
    static ityp rcomp,t,t1,u;
    rcomp = 0.0e0;
    if(*a >= 20.0e0) goto S20;
    t = *a*log(*x)-*x;
    if(*a >= 1.0e0) goto S10;
    rcomp = *a*exp(t)*(1.0e0+gam1(a));
    result = rcomp;
    return &result;
    S10:
        rcomp = exp(t)/ gamma_x(a);
        result = rcomp;
    return &result;
    S20:
        u = *x/ *a;
        result = rcomp;
        if(u == 0.0e0) return &result;
        t = pow(1.0e0/ *a,2.0);
        t1 = (((0.75e0*t-1.0e0)*t+3.5e0)*t-105.0e0)/(*a*1260.0e0);
        t1 -= (*a*rlog(&u));
        rcomp = rt2pin*sqrt(*a)*exp(t1);
        
    result = rcomp;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _rexp ( void * data)
/******************************************************************************/
//  Purpose:
//    REXP evaluates the function EXP(X) - 1.
//  Modified:
//    09 December 1999
//  Parameters:
//    Input, double *X, the argument of the function.
//    Output, double REXP, the value of EXP(X)-1.
{
	static ityp result = MAX_VAL;
	
	ityp * x = data;
	
    static ityp p1 = .914041914819518e-09;
    static ityp p2 = .238082361044469e-01;
    static ityp q1 = -.499999999085958e+00;
    static ityp q2 = .107141568980644e+00;
    static ityp q3 = -.119041179760821e-01;
    static ityp q4 = .595130811860248e-03;
    static ityp rexp,w;

    if(fabs(*x) > 0.15e0) goto S10;
    rexp = *x*(((p2**x+p1)**x+1.0e0)/((((q4**x+q3)**x+q2)**x+q1)**x+1.0e0));
    result = rexp;
    return &result;
    S10:
        w = exp(*x);
        if(*x > 0.0e0) goto S20;
        rexp = w-0.5e0-0.5e0;
        result = rexp;
    return &result;
    S20:
        rexp = w*(0.5e0+(0.5e0-1.0e0/w));
        
    result = rexp;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _rlog ( void * data)
/******************************************************************************/
//  Purpose:
//    RLOG computes  X - 1 - LN(X).
//  Modified:
//    09 December 1999
//  Parameters:
//    Input, double *X, the argument of the function.
//    Output, double RLOG, the value of the function.
{
	static ityp result = MAX_VAL;
	
	ityp * x = data;
	
    static ityp a = .566749439387324e-01;
    static ityp b = .456512608815524e-01;
    static ityp p0 = .333333333333333e+00;
    static ityp p1 = -.224696413112536e+00;
    static ityp p2 = .620886815375787e-02;
    static ityp q1 = -.127408923933623e+01;
    static ityp q2 = .354508718369557e+00;
    static ityp rlog,r,t,u,w,w1;

    if(*x < 0.61e0 || *x > 1.57e0) goto S40;
    if(*x < 0.82e0) goto S10;
    if(*x > 1.18e0) goto S20;
    //              ARGUMENT REDUCTION
    u = *x-0.5e0-0.5e0;
    w1 = 0.0e0;
    goto S30;
    S10:
        u = *x-0.7e0;
        u /= 0.7e0;
        w1 = a-u*0.3e0;
    goto S30;
    S20:
        u = 0.75e0**x-1.e0;
        w1 = b+u/3.0e0;
    S30:
        //               SERIES EXPANSION
        r = u/(u+2.0e0);
        t = r*r;
        w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0e0);
        rlog = 2.0e0*t*(1.0e0/(1.0e0-r)-r*w)+w1;
        result = rlog;
        return &result;
    S40:
        r = *x-0.5e0-0.5e0;
        rlog = r-log(*x);
        
    result = rlog;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _rlog1 ( void * data)
/******************************************************************************/
//  Purpose:
//    RLOG1 evaluates the function X - ln ( 1 + X ).
//  Parameters:
//    Input, double *X, the argument.
//    Output, double RLOG1, the value of X - ln ( 1 + X ).
{
	static ityp result = MAX_VAL;
	
	ityp * x = data;
	
    static ityp a = .566749439387324e-01;
    static ityp b = .456512608815524e-01;
    static ityp p0 = .333333333333333e+00;
    static ityp p1 = -.224696413112536e+00;
    static ityp p2 = .620886815375787e-02;
    static ityp q1 = -.127408923933623e+01;
    static ityp q2 = .354508718369557e+00;
    static ityp rlog1,h,r,t,w,w1;

    if(*x < -0.39e0 || *x > 0.57e0) goto S40;
    if(*x < -0.18e0) goto S10;
    if(*x > 0.18e0) goto S20;
    //
    //              ARGUMENT REDUCTION
    //
    h = *x;
    w1 = 0.0e0;
    goto S30;
        S10:
        h = *x+0.3e0;
        h /= 0.7e0;
        w1 = a-h*0.3e0;
    goto S30;
    S20:
        h = 0.75e0**x-0.25e0;
        w1 = b+h/3.0e0;
    S30:
        //
        //               SERIES EXPANSION
        //
        r = h/(h+2.0e0);
        t = r*r;
        w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0e0);
        rlog1 = 2.0e0*t*(1.0e0/(1.0e0-r)-r*w)+w1;
        result = rlog1;
        return &result;
    S40:
        w = *x+0.5e0+0.5e0;
        rlog1 = *x-log(w);
        
    result = rlog1;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _stvaln ( void * data)
/******************************************************************************/
//  Purpose:
//    STVALN provides starting values for the inverse of the normal distribution.
//  Discussion:
//    The routine returns X such that
//      P = CUMNOR(X),
//    that is,
//      P = Integral from -infinity to X of (1/SQRT(2*M_PI)) EXP(-U*U/2) dU.
//  Reference:
//    Kennedy and Gentle,
//    Statistical Computing,
//    Marcel Dekker, NY, 1980, page 95,
//    QA276.4  K46
//  Parameters:
//    Input, double *P, the probability whose normal deviate
//    is sought.
//    Output, double STVALN, the normal deviate whose probability
//    is P.
{
	static ityp result = MAX_VAL;
	
	ityp * p = data;
	
    static ityp xden[5] =
    {
        0.993484626060e-1,
        0.588581570495e0,
        0.531103462366e0,
        0.103537752850e0,
        0.38560700634e-2
    };
    static ityp xnum[5] =
    {
        -0.322232431088e0,
        -1.000000000000e0,
        -0.342242088547e0,
        -0.204231210245e-1,
        -0.453642210148e-4
    };
    static dim_typ K1 = 5;
    static ityp stvaln,sign,y,z;

    if(!(*p <= 0.5e0)) goto S10;
    sign = -1.0e0;
    z = *p;
    goto S20;
    S10:
        sign = 1.0e0;
        z = 1.0e0-*p;
    S20:
        y = sqrt(-(2.0e0*log(z)));
        stvaln = y+ eval_pol ( xnum, &K1, &y ) / eval_pol ( xden, &K1, &y );
        stvaln = sign*stvaln;
        
    result = stvaln;
    return &result;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _binomial_cdf_values ( void * data)
/******************************************************************************/
//  Purpose:
//    BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
//  Discussion:
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
//  Modified:
//    31 May 2004
//  Author:
//    John Burkardt
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, int *A, double *B, the parameters of the function.
//    Output, int *X, the argument of the function.
//    Output, double *FX, the value of the function.
{
	const _3pdt2pit * const s_data = data;
	
	dim_typ * n_data = s_data->a0;
	dim_typ * a = s_data->a1;
	dim_typ *x = s_data->a2;
	ityp * b = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 17

    dim_typ a_vec[N_MAX] =
    {
        2,  2,  2,  2,
        2,  4,  4,  4,
        4, 10, 10, 10,
        10, 10, 10, 10,
        10
    };
    ityp b_vec[N_MAX] =
    {
        0.05E+00, 0.05E+00, 0.05E+00, 0.50E+00,
        0.50E+00, 0.25E+00, 0.25E+00, 0.25E+00,
        0.25E+00, 0.05E+00, 0.10E+00, 0.15E+00,
        0.20E+00, 0.25E+00, 0.30E+00, 0.40E+00,
        0.50E+00
    };
    ityp fx_vec[N_MAX] =
    {
        0.9025E+00, 0.9975E+00, 1.0000E+00, 0.2500E+00,
        0.7500E+00, 0.3164E+00, 0.7383E+00, 0.9492E+00,
        0.9961E+00, 0.9999E+00, 0.9984E+00, 0.9901E+00,
        0.9672E+00, 0.9219E+00, 0.8497E+00, 0.6331E+00,
        0.3770E+00
    };
    dim_typ x_vec[N_MAX] =
    {
        0, 1, 2, 0,
        1, 0, 1, 2,
        3, 4, 4, 4,
        4, 4, 4, 4,
        4
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *a = *x = 0;
        *b = *fx = 0.0E+00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }
    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _chi_noncentral_cdf_values ( void * data)
/******************************************************************************/
//  Purpose:
//    CHI_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
//  Discussion:
//    The CDF of the noncentral chi square distribution can be evaluated
//    within Mathematica by commands such as:
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF [ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
//  Modified:
//    12 June 2004
//  Author:
//    John Burkardt
//  Reference:
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, double *X, the argument of the function.
//    Output, double *LAMBDA, the noncentrality parameter.
//    Output, int *DF, the number of degrees of freedom.
//    Output, double *CDF, the noncentral chi CDF.
{
	const pdt2pitpdtpit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * lambda = s_data->a2;
	dim_typ * df = s_data->a3;
	ityp * cdf = s_data->a4;
	
    # define N_MAX 27

    ityp cdf_vec[N_MAX] =
    {
        0.839944E+00, 0.695906E+00, 0.535088E+00,
        0.764784E+00, 0.620644E+00, 0.469167E+00,
        0.307088E+00, 0.220382E+00, 0.150025E+00,
        0.307116E-02, 0.176398E-02, 0.981679E-03,
        0.165175E-01, 0.202342E-03, 0.498448E-06,
        0.151325E-01, 0.209041E-02, 0.246502E-03,
        0.263684E-01, 0.185798E-01, 0.130574E-01,
        0.583804E-01, 0.424978E-01, 0.308214E-01,
        0.105788E+00, 0.794084E-01, 0.593201E-01
    };
    dim_typ df_vec[N_MAX] =
    {
        1,   2,   3,
        1,   2,   3,
        1,   2,   3,
        1,   2,   3,
        60,  80, 100,
        1,   2,   3,
        10,  10,  10,
        10,  10,  10,
        10,  10,  10
    };
    ityp lambda_vec[N_MAX] =
    {
        0.5E+00,  0.5E+00,  0.5E+00,
        1.0E+00,  1.0E+00,  1.0E+00,
        5.0E+00,  5.0E+00,  5.0E+00,
        20.0E+00, 20.0E+00, 20.0E+00,
        30.0E+00, 30.0E+00, 30.0E+00,
        5.0E+00,  5.0E+00,  5.0E+00,
        2.0E+00,  3.0E+00,  4.0E+00,
        2.0E+00,  3.0E+00,  4.0E+00,
        2.0E+00,  3.0E+00,  4.0E+00
    };
    ityp x_vec[N_MAX] =
    {
        3.000E+00,  3.000E+00,  3.000E+00,
        3.000E+00,  3.000E+00,  3.000E+00,
        3.000E+00,  3.000E+00,  3.000E+00,
        3.000E+00,  3.000E+00,  3.000E+00,
        60.000E+00, 60.000E+00, 60.000E+00,
        0.050E+00,  0.050E+00,  0.050E+00,
        4.000E+00,  4.000E+00,  4.000E+00,
        5.000E+00,  5.000E+00,  5.000E+00,
        6.000E+00,  6.000E+00,  6.000E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *df = 0;
        *x = *lambda = *cdf = 0.0E+00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *lambda = lambda_vec[*n_data-1];
        *df = df_vec[*n_data-1];
        *cdf = cdf_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _chi_square_cdf_values ( void * data)
/******************************************************************************/
//  Purpose:
//    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
//  Discussion:
//    The value of CHI_CDF ( DF, X ) can be evaluated in Mathematica by
//    commands like:
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[ChiSquareDistribution[DF], X ]
//  Modified:
//    11 June 2004
//  Author:
//    John Burkardt
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, int *A, the parameter of the function.
//    Output, double *X, the argument of the function.
//    Output, double *FX, the value of the function.
{
	const _2pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * a = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 21

    dim_typ a_vec[N_MAX] =
    {
        1,  2,  1,  2,
        1,  2,  3,  4,
        1,  2,  3,  4,
        5,  3,  3,  3,
        3,  3, 10, 10,
        10
    };
    ityp fx_vec[N_MAX] =
    {
        0.0796557E+00, 0.00498752E+00, 0.112463E+00,    0.00995017E+00,
        0.472911E+00,  0.181269E+00,   0.0597575E+00,   0.0175231E+00,
        0.682689E+00,  0.393469E+00,   0.198748E+00,    0.090204E+00,
        0.0374342E+00, 0.427593E+00,   0.608375E+00,    0.738536E+00,
        0.828203E+00,  0.88839E+00,    0.000172116E+00, 0.00365985E+00,
        0.0185759E+00
    };
    ityp x_vec[N_MAX] =
    {
        0.01E+00, 0.01E+00, 0.02E+00, 0.02E+00,
        0.40E+00, 0.40E+00, 0.40E+00, 0.40E+00,
        1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
        1.00E+00, 2.00E+00, 3.00E+00, 4.00E+00,
        5.00E+00, 6.00E+00, 1.00E+00, 2.00E+00,
        3.00E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *a = 0;
        *x = *fx = 0.0E+00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }
    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _erf_values ( void * data)
/******************************************************************************/
//  Purpose:
//    ERF_VALUES returns some values of the ERF or "error" function.
//  Definition:
//    ERF(X) = ( 2 / sqrt ( M_PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
//  Modified:
//    31 May 2004
//  Author:
//    John Burkardt
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, double *X, the argument of the function.
//    Output, double *FX, the value of the function.
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 21

    ityp fx_vec[N_MAX] =
    {
        0.0000000000E+00, 0.1124629160E+00, 0.2227025892E+00, 0.3286267595E+00,
        0.4283923550E+00, 0.5204998778E+00, 0.6038560908E+00, 0.6778011938E+00,
        0.7421009647E+00, 0.7969082124E+00, 0.8427007929E+00, 0.8802050696E+00,
        0.9103139782E+00, 0.9340079449E+00, 0.9522851198E+00, 0.9661051465E+00,
        0.9763483833E+00, 0.9837904586E+00, 0.9890905016E+00, 0.9927904292E+00,
        0.9953222650E+00
    };
    ityp x_vec[N_MAX] =
    {
        0.0E+00, 0.1E+00, 0.2E+00, 0.3E+00,
        0.4E+00, 0.5E+00, 0.6E+00, 0.7E+00,
        0.8E+00, 0.9E+00, 1.0E+00, 1.1E+00,
        1.2E+00, 1.3E+00, 1.4E+00, 1.5E+00,
        1.6E+00, 1.7E+00, 1.8E+00, 1.9E+00,
        2.0E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *x = *fx = 0.0E+00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }
    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _f_cdf_values ( void * data)
/******************************************************************************/
//  Purpose:
//    F_CDF_VALUES returns some values of the F CDF test function.
//  Discussion:
//    The value of F_CDF ( DFN, DFD, X ) can be evaluated in Mathematica by
//    commands like:
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[FRatioDistribution[ DFN, DFD ], X ]
//  Modified:
//    11 June 2004
//  Author:
//    John Burkardt
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, int *A, int *B, the parameters of the function.
//    Output, double *X, the argument of the function.
//    Output, double *FX, the value of the function.
{
	const _3pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * a = s_data->a1;
	dim_typ * b = s_data->a2;
	ityp * x = s_data->a3;
	ityp * fx = s_data->a4;
	
    # define N_MAX 20

    dim_typ a_vec[N_MAX] =
    {
        1, 1, 5, 1,
        2, 4, 1, 6,
        8, 1, 3, 6,
        1, 1, 1, 1,
        2, 3, 4, 5
    };
    dim_typ b_vec[N_MAX] =
    {
        1,  5,  1,  5,
        10, 20,  5,  6,
        16,  5, 10, 12,
        5,  5,  5,  5,
        5,  5,  5,  5
    };
    ityp fx_vec[N_MAX] =
    {
        0.500000E+00, 0.499971E+00, 0.499603E+00, 0.749699E+00,
        0.750466E+00, 0.751416E+00, 0.899987E+00, 0.899713E+00,
        0.900285E+00, 0.950025E+00, 0.950057E+00, 0.950193E+00,
        0.975013E+00, 0.990002E+00, 0.994998E+00, 0.999000E+00,
        0.568799E+00, 0.535145E+00, 0.514343E+00, 0.500000E+00
    };
    ityp x_vec[N_MAX] =
    {
        1.00E+00,  0.528E+00, 1.89E+00,  1.69E+00,
        1.60E+00,  1.47E+00,  4.06E+00,  3.05E+00,
        2.09E+00,  6.61E+00,  3.71E+00,  3.00E+00,
        10.01E+00, 16.26E+00, 22.78E+00, 47.18E+00,
        1.00E+00,  1.00E+00,  1.00E+00,  1.00E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *a = *b = 0;
        *x = *fx = 0.0E+00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }
    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _f_noncentral_cdf_values ( void * data)
/******************************************************************************/
//  Purpose:
//    F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
//  Discussion:
//    The value of NONCENTRAL_F_CDF ( DFN, DFD, LAMDA, X ) can be evaluated
//    in Mathematica by commands like:
//      Needs["Statistics`ContinuousDistributions`"]
//      CDF[NoncentralFRatioDistribution[ DFN, DFD, LAMBDA ], X ]
//  Modified:
//    12 June 2004
//  Author:
//    John Burkardt
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, int *A, int *B, double *LAMBDA, the
//    parameters of the function.
//    Output, double *X, the argument of the function.
//    Output, double *FX, the value of the function.
{
	const _3pdt3pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * a = s_data->a1;
	dim_typ * b = s_data->a2;
	ityp * lambda = s_data->a3;
	ityp * x = s_data->a4;
	ityp * fx =  s_data->a5;
	
	
    # define N_MAX 22

    dim_typ a_vec[N_MAX] =
    {
        1,  1,  1,  1,
        1,  1,  1,  1,
        1,  1,  2,  2,
        3,  3,  4,  4,
        5,  5,  6,  6,
        8, 16
    };
    dim_typ b_vec[N_MAX] =
    {
        1,  5,  5,  5,
        5,  5,  5,  5,
        5,  5,  5, 10,
        5,  5,  5,  5,
        1,  5,  6, 12,
        16,  8
    };
    ityp fx_vec[N_MAX] =
    {
        0.500000E+00, 0.636783E+00, 0.584092E+00, 0.323443E+00,
        0.450119E+00, 0.607888E+00, 0.705928E+00, 0.772178E+00,
        0.819105E+00, 0.317035E+00, 0.432722E+00, 0.450270E+00,
        0.426188E+00, 0.337744E+00, 0.422911E+00, 0.692767E+00,
        0.363217E+00, 0.421005E+00, 0.426667E+00, 0.446402E+00,
        0.844589E+00, 0.816368E+00
    };
    ityp lambda_vec[N_MAX] =
    {
        0.00E+00,  0.000E+00, 0.25E+00,  1.00E+00,
        1.00E+00,  1.00E+00,  1.00E+00,  1.00E+00,
        1.00E+00,  2.00E+00,  1.00E+00,  1.00E+00,
        1.00E+00,  2.00E+00,  1.00E+00,  1.00E+00,
        0.00E+00,  1.00E+00,  1.00E+00,  1.00E+00,
        1.00E+00,  1.00E+00
    };
    ityp x_vec[N_MAX] =
    {
        1.00E+00,  1.00E+00, 1.00E+00,  0.50E+00,
        1.00E+00,  2.00E+00, 3.00E+00,  4.00E+00,
        5.00E+00,  1.00E+00, 1.00E+00,  1.00E+00,
        1.00E+00,  1.00E+00, 1.00E+00,  2.00E+00,
        1.00E+00,  1.00E+00, 1.00E+00,  1.00E+00,
        2.00E+00,  2.00E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *a = *b = 0;
        *x = *fx = *lambda = 0.0E+00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *b = b_vec[*n_data-1];
        *lambda = lambda_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _negative_binomial_cdf_values ( void * data)
/******************************************************************************/
//  Purpose:
//    NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
//  Discussion:
//    Assume that a coin has a probability P of coming up heads on
//    any one trial.  Suppose that we plan to flip the coin until we
//    achieve a total of S heads.  If we let F represent the number of
//    tails that occur in this process, then the value of F satisfies
//    a negative binomial PDF:
//      PDF(F,S,P) = Choose ( F from F+S-1 ) * P**S * (1-P)**F
//    The negative binomial CDF is the probability that there are F or
//    fewer failures upon the attainment of the S-th success.  Thus,
//      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
//  Modified:
//    07 June 2004
//  Author:
//    John Burkardt
//  Reference:
//    F C Powell,
//    Statistical Tables for Sociology, Biology and Physical Sciences,
//    Cambridge University Press, 1982.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, int *F, the maximum number of failures.
//    Output, int *S, the number of successes.
//    Output, double *P, the probability of a success on one trial.
//    Output, double *CDF, the probability of at most F failures before the
//    S-th success.
{
	const _3pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * f = s_data->a1;
	dim_typ * s = s_data->a2;
	ityp * p = s_data->a3;
	ityp * cdf = s_data->a4;
	
    # define N_MAX 27

    ityp cdf_vec[N_MAX] =
    {
        0.6367, 0.3633, 0.1445,
        0.5000, 0.2266, 0.0625,
        0.3438, 0.1094, 0.0156,
        0.1792, 0.0410, 0.0041,
        0.0705, 0.0109, 0.0007,
        0.9862, 0.9150, 0.7472,
        0.8499, 0.5497, 0.2662,
        0.6513, 0.2639, 0.0702,
        1.0000, 0.0199, 0.0001
    };
    dim_typ f_vec[N_MAX] =
    {
        4,  3,  2,
        3,  2,  1,
        2,  1,  0,
        2,  1,  0,
        2,  1,  0,
        11, 10,  9,
        17, 16, 15,
        9,  8,  7,
        2,  1,  0
    };
    ityp p_vec[N_MAX] =
    {
        0.50, 0.50, 0.50,
        0.50, 0.50, 0.50,
        0.50, 0.50, 0.50,
        0.40, 0.40, 0.40,
        0.30, 0.30, 0.30,
        0.30, 0.30, 0.30,
        0.10, 0.10, 0.10,
        0.10, 0.10, 0.10,
        0.01, 0.01, 0.01
    };
    dim_typ s_vec[N_MAX] =
    {
        4, 5, 6,
        4, 5, 6,
        4, 5, 6,
        4, 5, 6,
        4, 5, 6,
        1, 2, 3,
        1, 2, 3,
        1, 2, 3,
        0, 1, 2
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *f = *s = 0;
        *p = *cdf = 0.0E+00;
    }
    else
    {
        *f = f_vec[*n_data-1];
        *s = s_vec[*n_data-1];
        *p = p_vec[*n_data-1];
        *cdf = cdf_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *   _normal_cdf_values ( void * data)
/******************************************************************************/
//  Purpose:
//    NORMAL_CDF_VALUES returns some values of the Normal CDF.
//  Modified:
//    31 May 2004
//  Author:
//    John Burkardt
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, double *X, the argument of the function.
//    Output double *FX, the value of the function.
{
	const pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	ityp * x = s_data->a1;
	ityp * fx = s_data->a2;
	
    # define N_MAX 13

    ityp fx_vec[N_MAX] =
    {
        0.500000000000000E+00, 0.539827837277029E+00, 0.579259709439103E+00,
        0.617911422188953E+00, 0.655421741610324E+00, 0.691462461274013E+00,
        0.725746882249927E+00, 0.758036347776927E+00, 0.788144601416604E+00,
        0.815939874653241E+00, 0.841344746068543E+00, 0.933192798731142E+00,
        0.977249868051821E+00
    };
    ityp x_vec[N_MAX] =
    {
        0.00E+00, 0.10E+00, 0.20E+00,
        0.30E+00, 0.40E+00, 0.50E+00,
        0.60E+00, 0.70E+00, 0.80E+00,
        0.90E+00, 1.00E+00, 1.50E+00,
        2.00E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
	    *n_data = 0;
	    *x = *fx = 0.0E+00;
    }
    else
    {
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _poisson_cdf_values ( void * data)
/******************************************************************************/
//  Purpose:
//    POISSON_CDF_VALUES returns some values of the Poisson CDF.
//  Discussion:
//    CDF(X)(A) is the probability of at most X successes in unit time,
//    given that the expected mean number of successes is A.
//  Modified:
//    31 May 2004
//  Author:
//    John Burkardt
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 653-658.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, double *A, the parameter of the function.
//    Output, int *X, the argument of the function.
//    Output, double *FX, the value of the function.
{
	const _2pdt2pit * const s_data = data;
	
	dim_typ * n_data = s_data->a0;
	dim_typ * x = s_data->a1;
	ityp * a = s_data->a2;
	ityp * fx =  s_data->a3;
	
    # define N_MAX 21

    ityp a_vec[N_MAX] =
    {
        0.02E+00, 0.10E+00, 0.10E+00, 0.50E+00,
        0.50E+00, 0.50E+00, 1.00E+00, 1.00E+00,
        1.00E+00, 1.00E+00, 2.00E+00, 2.00E+00,
        2.00E+00, 2.00E+00, 5.00E+00, 5.00E+00,
        5.00E+00, 5.00E+00, 5.00E+00, 5.00E+00,
        5.00E+00
    };
    ityp fx_vec[N_MAX] =
    {
        0.980E+00, 0.905E+00, 0.995E+00, 0.607E+00,
        0.910E+00, 0.986E+00, 0.368E+00, 0.736E+00,
        0.920E+00, 0.981E+00, 0.135E+00, 0.406E+00,
        0.677E+00, 0.857E+00, 0.007E+00, 0.040E+00,
        0.125E+00, 0.265E+00, 0.441E+00, 0.616E+00,
        0.762E+00
    };
    dim_typ x_vec[N_MAX] =
    {
        0, 0, 1, 0,
        1, 2, 0, 1,
        2, 3, 0, 1,
        2, 3, 0, 1,
        2, 3, 4, 5,
        6
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *x = 0;
        *fx = *a = 0.0E+00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }
    return NULL;
    # undef N_MAX
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _student_cdf_values ( void * data)
/******************************************************************************/
//  Purpose:
//    STUDENT_CDF_VALUES returns some values of the Student CDF.
//  Modified:
//    31 May 2004
//  Author:
//    John Burkardt
//  Reference:
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//  Parameters:
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//    Output, int *A, the parameter of the function.
//    Output, double *X, the argument of the function.
//    Output, double *FX, the value of the function.
{
	const _2pdt2pit * const s_data = data;
	dim_typ * n_data = s_data->a0;
	dim_typ * a = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 13

    dim_typ a_vec[N_MAX] =
    {
        1, 2, 3, 4,
        5, 2, 5, 2,
        5, 2, 3, 4,
        5
    };
    ityp fx_vec[N_MAX] =
    {
        0.60E+00, 0.60E+00, 0.60E+00, 0.60E+00,
        0.60E+00, 0.75E+00, 0.75E+00, 0.95E+00,
        0.95E+00, 0.99E+00, 0.99E+00, 0.99E+00,
        0.99E+00
    };
    ityp x_vec[N_MAX] =
    {
        0.325E+00, 0.289E+00, 0.277E+00, 0.271E+00,
        0.267E+00, 0.816E+00, 0.727E+00, 2.920E+00,
        2.015E+00, 6.965E+00, 4.541E+00, 3.747E+00,
        3.365E+00
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = *a = *fx = 0;
        *x = *fx = 0.0E+00;
        *fx = 0.0E+00;
    }
    else
    {
        *a = a_vec[*n_data-1];
        *x = x_vec[*n_data-1];
        *fx = fx_vec[*n_data-1];
    }

    return NULL;
    # undef N_MAX
}

#endif
