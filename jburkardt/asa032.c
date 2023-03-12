#ifndef __DISABLEDEEP_ASA032

#include "../dutils.h"

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _gamma_inc (void *data)
/******************************************************************************/
//  Purpose:
//    GAMMA_INC evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
//  Discussion:
//    This is certified spaghetti code.
//  Author:
//    Alfred H Morris, Jr,
//    Naval Surface Weapons Center,
//    Dahlgren, Virginia.
//  Parameters:
//    Input, double *A, *X, the arguments of the incomplete
//    gamma ratio.  A and X must be nonnegative.  A and X cannot
//    both be zero.
//    Output, double *ANS, *QANS.  On normal output,
//    ANS = P(A,X) and QANS = Q(A,X).  However, ANS is set to 2 if
//    A or X is negative, or both are 0, or when the answer is
//    computationally indeterminate because A is extremely large
//    and X is very close to A.
//    Input, int *IND, indicates the accuracy request:
//    0, as much accuracy as possible.
//    1, to within 1 unit of the 6-th significant digit,
//    otherwise, to within 1 unit of the 3rd significant digit.
{
	const pdt4pit * const s_data = data;
	
	dim_typ * ind = s_data->a0;
	ityp * a = s_data->a1;
	ityp * x = s_data->a2;
	ityp * ans = s_data->a3;
	ityp * qans = s_data->a4;
	
	
    static ityp alog10 = 2.30258509299405e0;
    static ityp d10 = -.185185185185185e-02;
    static ityp d20 = .413359788359788e-02;
    static ityp d30 = .649434156378601e-03;
    static ityp d40 = -.861888290916712e-03;
    static ityp d50 = -.336798553366358e-03;
    static ityp d60 = .531307936463992e-03;
    static ityp d70 = .344367606892378e-03;
    static ityp rt2pin = .398942280401433e0;
    static ityp rtpi = 1.77245385090552e0;
    static ityp third = .333333333333333e0;
    static ityp acc0[3] =
    {
        5.e-15,
        5.e-7,
        5.e-4
    };
    static ityp big[3] =
    {
        20.0e0,
        14.0e0,
        10.0e0
    };
    static ityp d0[13] =
    {
        .833333333333333e-01,
        -.148148148148148e-01,
        .115740740740741e-02,
        .352733686067019e-03,
        -.178755144032922e-03,
        .391926317852244e-04,
        -.218544851067999e-05,
        -.185406221071516e-05,
        .829671134095309e-06,
        -.176659527368261e-06,
        .670785354340150e-08,
        .102618097842403e-07,
        -.438203601845335e-08
    };
    static ityp d1[12] =
    {
        -.347222222222222e-02,
        .264550264550265e-02,
        -.990226337448560e-03,
        .205761316872428e-03,
        -.401877572016461e-06,
        -.180985503344900e-04,
        .764916091608111e-05,
        -.161209008945634e-05
        ,.464712780280743e-08,
        .137863344691572e-06,
        -.575254560351770e-07,
        .119516285997781e-07
        };
    static ityp d2[10] =
    {
        -.268132716049383e-02,
        .771604938271605e-03,
        .200938786008230e-05,
        -.107366532263652e-03,
        .529234488291201e-04,
        -.127606351886187e-04,
        .342357873409614e-07,
        .137219573090629e-05,
        -.629899213838006e-06,
        .142806142060642e-06
    };
    static ityp d3[8] =
    {
        .229472093621399e-03,
        -.469189494395256e-03,
        .267720632062839e-03,
        -.756180167188398e-04,
        -.239650511386730e-06,
        .110826541153473e-04,
        -.567495282699160e-05,
        .142309007324359e-05
    };
    static ityp d4[6] =
    {
        .784039221720067e-03,
        -.299072480303190e-03,
        -.146384525788434e-05,
        .664149821546512e-04,
        -.396836504717943e-04,
        .113757269706784e-04
    };
    static ityp d5[4] =
    {
        -.697281375836586e-04,
        .277275324495939e-03,
        -.199325705161888e-03,
        .679778047793721e-04
    };
    static ityp d6[2] =
    {
        -.592166437353694e-03,
        .270878209671804e-03
    };
    static ityp e00[3] =
    {
        .25e-3,
        .25e-1,
        .14e0
    };
    static ityp x00[3] =
    {
        31.0e0,
        17.0e0,
        9.7e0
    };
    static dim_typ K1 = 1;
    static dim_typ K2 = 0;
    static ityp a2n,a2nm1,acc,am0,amn,an,an0,apn,b2n,b2nm1,c,c0,c1,c2,c3,c4,c5,c6,cma,e,e0,g,h,j,l,r,rta,rtx,s,sum,t,t1,tol,twoa,u,w,x0,y,z;
    static dim_typ i,iop,m,MAX,n;
    static ityp wk[20],T3;
    static dim_typ T4,T5;
    static ityp T6,T7;

    //
    //  E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
    //  NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
    //
    e = dpmpar(&K1);
    if(*a < 0.0e0 || *x < 0.0e0) goto S430;
    if(*a == 0.0e0 && *x == 0.0e0) goto S430;
    if(*a**x == 0.0e0) goto S420;
    iop = *ind+1;
    if(iop != 1 && iop != 2) iop = 3;
    acc = fifdmax1(acc0[iop-1],e);
    e0 = e00[iop-1];
    x0 = x00[iop-1];
    //
    //  SELECT THE APPROPRIATE ALGORITHM
    //
    if(*a >= 1.0e0) goto S10;
    if(*a == 0.5e0) goto S390;
    if(*x < 1.1e0) goto S160;
    t1 = *a*log(*x)-*x;
    u = *a*exp(t1);
    if(u == 0.0e0) goto S380;
    r = u*(1.0e0+gam1(a));
    goto S250;
    S10:
        if(*a >= big[iop-1]) goto S30;
        if(*a > *x || *x >= x0) goto S20;
        twoa = *a+*a;
        m = fifidint(twoa);
        if(twoa != (double)m) goto S20;
        i = m/2;
        if(*a == (double)i) goto S210;
        goto S220;
    S20:
        t1 = *a*log(*x)-*x;
        r = exp(t1)/ gamma_x(a);
        goto S40;
    S30:
        l = *x/ *a;
        if(l == 0.0e0) goto S370;
        s = 0.5e0+(0.5e0-l);
        z = rlog(&l);
        if(z >= 700.0e0/ *a) goto S410;
        y = *a*z;
        rta = sqrt(*a);
        if(fabs(s) <= e0/rta) goto S330;
        if(fabs(s) <= 0.4e0) goto S270;
        t = pow(1.0e0/ *a,2.0);
        t1 = (((0.75e0*t-1.0e0)*t+3.5e0)*t-105.0e0)/(*a*1260.0e0);
        t1 -= y;
        r = rt2pin*rta*exp(t1);
    S40:
        if(r == 0.0e0) goto S420;
        if(*x <= fifdmax1(*a,alog10)) goto S50;
        if(*x < x0) goto S250;
        goto S100;
    S50:
        //  TAYLOR SERIES FOR P/R
        apn = *a+1.0e0;
        t = *x/apn;
        wk[0] = t;
        for ( n = 2; n <= 20; n++ )
        {
            apn += 1.0e0;
            t *= (*x/apn);
            if(t <= 1.e-3) goto S70;
            wk[n-1] = t;
        }
        n = 20;
    S70:
        sum = t;
        tol = 0.5e0*acc;
    S80:
        apn += 1.0e0;
        t *= (*x/apn);
        sum += t;
        if(t > tol) goto S80;
        MAX = n-1;
        for ( m = 1; m <= MAX; ++m)
        {
            -- n;
            sum += wk[n-1];
        }
        *ans = r/ *a*(1.0e0+sum);
        *qans = 0.5e0+(0.5e0-*ans);
        return NULL;
    S100:
        //  ASYMPTOTIC EXPANSION
        amn = *a-1.0e0;
        t = amn/ *x;
        wk[0] = t;
        for ( n = 2; n <= 20; ++n )
        {
            amn -= 1.0e0;
            t *= (amn/ *x);
            if(fabs(t) <= 1.e-3) goto S120;
            wk[n-1] = t;
        }
        n = 20;
        S120:
        sum = t;
    S130:
        if(fabs(t) <= acc) goto S140;
        amn -= 1.0e0;
        t *= (amn/ *x);
        sum += t;
        goto S130;
    S140:
        MAX = n-1;
        for ( m = 1; m <= MAX; ++m )
        {
            -- n;
            sum += wk[n-1];
        }
        *qans = r/ *x*(1.0e0+sum);
        *ans = 0.5e0+(0.5e0-*qans);
        return NULL;
    S160:
        //  TAYLOR SERIES FOR P(A,X)/X**A
        an = 3.0e0;
        c = *x;
        sum = *x/(*a+3.0e0);
        tol = 3.0e0*acc/(*a+1.0e0);
    S170:
        an += 1.0e0;
        c = -(c*(*x/an));
        t = c/(*a+an);
        sum += t;
        if(fabs(t) > tol) goto S170;
        j = *a**x*((sum/6.0e0-0.5e0/(*a+2.0e0))**x+1.0e0/(*a+1.0e0));
        z = *a*log(*x);
        h = gam1(a);
        g = 1.0e0+h;
        if(*x < 0.25e0) goto S180;
        if(*a < *x/2.59e0) goto S200;
        goto S190;
    S180:
        if(z > -.13394e0) goto S200;
    S190:
        w = exp(z);
        *ans = w*g*(0.5e0+(0.5e0-j));
        *qans = 0.5e0+(0.5e0-*ans);
        return NULL;
    S200:
        l = rexp(&z);
        w = 0.5e0+(0.5e0+l);
        *qans = (w*j-l)*g-h;
        if(*qans < 0.0e0) goto S380;
        *ans = 0.5e0+(0.5e0-*qans);
        return NULL;
    S210:
        //  FINITE SUMS FOR Q WHEN A .GE. 1 AND 2*A IS AN INTEGER
        sum = exp(-*x);
        t = sum;
        n = 1;
        c = 0.0e0;
        goto S230;
    S220:
        rtx = sqrt(*x);
        sum = error_fc ( &K2, &rtx );
        t = exp(-*x)/(rtpi*rtx);
        n = 0;
        c = -0.5e0;
    S230:
        if(n == i) goto S240;
        ++ n;
        c += 1.0e0;
        t = *x*t/c;
        sum += t;
        goto S230;
    S240:
        *qans = sum;
        *ans = 0.5e0+(0.5e0-*qans);
        return NULL;
    S250:
        //  CONTINUED FRACTION EXPANSION
        tol = fifdmax1(5.0e0*e,acc);
        a2nm1 = a2n = 1.0e0;
        b2nm1 = *x;
        b2n = *x+(1.0e0-*a);
        c = 1.0e0;
        S260:
        a2nm1 = *x*a2n+c*a2nm1;
        b2nm1 = *x*b2n+c*b2nm1;
        am0 = a2nm1/b2nm1;
        c += 1.0e0;
        cma = c-*a;
        a2n = a2nm1+cma*a2n;
        b2n = b2nm1+cma*b2n;
        an0 = a2n/b2n;
        if(fabs(an0-am0) >= tol*an0) goto S260;
        *qans = r*an0;
        *ans = 0.5e0+(0.5e0-*qans);
        return NULL;
    S270:
        //  GENERAL TEMME EXPANSION
        if(fabs(s) <= 2.0e0*e && *a*e*e > 3.28e-3) goto S430;
        c = exp(-y);
        T3 = sqrt(y);
        w = 0.5e0 * error_fc ( &K1, &T3 );
        u = 1.0e0/ *a;
        z = sqrt(z+z);
        if(l < 1.0e0) z = -z;
        T4 = iop-2;
        if(T4 < 0) goto S280;
        else if(T4 == 0) goto S290;
        else  goto S300;
        S280:
        if(fabs(s) <= 1.e-3) goto S340;
        c0 = ((((((((((((d0[12]*z+d0[11])*z+d0[10])*z+d0[9])*z+d0[8])*z+d0[7])*z+d0[
        6])*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
        c1 = (((((((((((d1[11]*z+d1[10])*z+d1[9])*z+d1[8])*z+d1[7])*z+d1[6])*z+d1[5]
        )*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
        c2 = (((((((((d2[9]*z+d2[8])*z+d2[7])*z+d2[6])*z+d2[5])*z+d2[4])*z+d2[3])*z+
        d2[2])*z+d2[1])*z+d2[0])*z+d20;
        c3 = (((((((d3[7]*z+d3[6])*z+d3[5])*z+d3[4])*z+d3[3])*z+d3[2])*z+d3[1])*z+
        d3[0])*z+d30;
        c4 = (((((d4[5]*z+d4[4])*z+d4[3])*z+d4[2])*z+d4[1])*z+d4[0])*z+d40;
        c5 = (((d5[3]*z+d5[2])*z+d5[1])*z+d5[0])*z+d50;
        c6 = (d6[1]*z+d6[0])*z+d60;
        t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
        goto S310;
    S290:
        c0 = (((((d0[5]*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
        c1 = (((d1[3]*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
        c2 = d2[0]*z+d20;
        t = (c2*u+c1)*u+c0;
        goto S310;
    S300:
        t = ((d0[2]*z+d0[1])*z+d0[0])*z-third;
    S310:
        if(l < 1.0e0) goto S320;
        *qans = c*(w+rt2pin*t/rta);
        *ans = 0.5e0+(0.5e0-*qans);
        return NULL;
    S320:
        *ans = c*(w-rt2pin*t/rta);
        *qans = 0.5e0+(0.5e0-*ans);
        return NULL;
    S330:
        //  TEMME EXPANSION FOR L = 1
        if(*a*e*e > 3.28e-3) goto S430;
        c = 0.5e0+(0.5e0-y);
        w = (0.5e0-sqrt(y)*(0.5e0+(0.5e0-y/3.0e0))/rtpi)/c;
        u = 1.0e0/ *a;
        z = sqrt(z+z);
        if(l < 1.0e0) z = -z;
        T5 = iop-2;
        if(T5 < 0) goto S340;
        else if(T5 == 0) goto S350;
        else  goto S360;
    S340:
        c0 = ((((((d0[6]*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-
        third;
        c1 = (((((d1[5]*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
        c2 = ((((d2[4]*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z+d20;
        c3 = (((d3[3]*z+d3[2])*z+d3[1])*z+d3[0])*z+d30;
        c4 = (d4[1]*z+d4[0])*z+d40;
        c5 = (d5[1]*z+d5[0])*z+d50;
        c6 = d6[0]*z+d60;
        t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
        goto S310;
    S350:
        c0 = (d0[1]*z+d0[0])*z-third;
        c1 = d1[0]*z+d10;
        t = (d20*u+c1)*u+c0;
        goto S310;
    S360:
        t = d0[0]*z-third;
        goto S310;
    S370:
        //  SPECIAL CASES
        *ans = 0.0e0;
        *qans = 1.0e0;
        return NULL;
    S380:
        *ans = 1.0e0;
        *qans = 0.0e0;
        return NULL;
    S390:
        if(*x >= 0.25e0) goto S400;
        T6 = sqrt(*x);
        *ans = error_f ( &T6 );
        *qans = 0.5e0+(0.5e0-*ans);
        return NULL;
    S400:
        T7 = sqrt(*x);
        *qans = error_fc ( &K2, &T7 );
        *ans = 0.5e0+(0.5e0-*qans);
        return NULL;
    S410:
        if(fabs(s) <= 2.0e0*e) goto S430;
    S420:
        if(*x <= *a) goto S370;
        goto S380;
    S430:
        //  ERROR RETURN

        *ans = 2.0e0;
    return NULL;
}

/******************************************************************************/
__MATHSUITE __JBURKARDT  void *  _gamma_inc_values ( void *data)
/******************************************************************************/
/*
  Purpose:
    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
  Discussion:
    The (normalized) incomplete Gamma function P(A,X) is defined as:
      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
    With this definition, for all A and X,
      0 <= PN(A,X) <= 1
    and
      PN(A,INFINITY) = 1.0
    In Mathematica, the function can be evaluated by:
      1 - GammaRegularized[A,X]
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    20 November 2004
  Author:
    John Burkardt
  Reference:
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  Parameters:
    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.
    Output, double *A, the parameter of the function.
    Output, double *X, the argument of the function.
    Output, double *FX, the value of the function.
*/
{
	const _2pdt2pit * const s_data = data; 
	dim_typ * n_data = s_data->a0;
	dim_typ * a = s_data->a1;
	ityp * x = s_data->a2;
	ityp * fx = s_data->a3;
	
    # define N_MAX 20

    ityp a_vec[N_MAX] =
    {
        0.10E+00,
        0.10E+00,
        0.10E+00,
        0.50E+00,
        0.50E+00,
        0.50E+00,
        0.10E+01,
        0.10E+01,
        0.10E+01,
        0.11E+01,
        0.11E+01,
        0.11E+01,
        0.20E+01,
        0.20E+01,
        0.20E+01,
        0.60E+01,
        0.60E+01,
        0.11E+02,
        0.26E+02,
        0.41E+02
    };

    ityp fx_vec[N_MAX] =
    {
        0.7382350532339351E+00,
        0.9083579897300343E+00,
        0.9886559833621947E+00,
        0.3014646416966613E+00,
        0.7793286380801532E+00,
        0.9918490284064973E+00,
        0.9516258196404043E-01,
        0.6321205588285577E+00,
        0.9932620530009145E+00,
        0.7205974576054322E-01,
        0.5891809618706485E+00,
        0.9915368159845525E+00,
        0.1018582711118352E-01,
        0.4421745996289254E+00,
        0.9927049442755639E+00,
        0.4202103819530612E-01,
        0.9796589705830716E+00,
        0.9226039842296429E+00,
        0.4470785799755852E+00,
        0.7444549220718699E+00
    };

    ityp x_vec[N_MAX] =
    {
        0.30E-01,
        0.30E+00,
        0.15E+01,
        0.75E-01,
        0.75E+00,
        0.35E+01,
        0.10E+00,
        0.10E+01,
        0.50E+01,
        0.10E+00,
        0.10E+01,
        0.50E+01,
        0.15E+00,
        0.15E+01,
        0.70E+01,
        0.25E+01,
        0.12E+02,
        0.16E+02,
        0.25E+02,
        0.45E+02
    };

    ++ *n_data;

    if ( N_MAX < *n_data )
    {
        *n_data = 0;
        *a = *x = *fx = 0.00;
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
