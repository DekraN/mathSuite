
// ext_math.h
// as External Math Library
// 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be included
exclusively by main.c, mengine.c, programs.c, algebra.c and settings.c project files of my suite program!
I do not assume any responsibilities about the use with any other code-scripts.
*/

// START

#include "dutils.h"
// #include "ExprEval/exprincl.h"

#define LAST_TYPERANGE DOMAIN_LDBL
#define MAX_TYPERANGE LAST_TYPERANGE+1


// This part, except the ext_math
// structured var initialization and definitions,
// is completely taken from BINARYAD.C
/// More informations about the Source Code, which here I adapted and modfied, at:
/// http://www.Planet-Source-Code.com/vb/scripts/ShowCode.asp?txtCodeId=13861&lngWId=3

#define SIZE (sizeof(int)<<3)
// FUNCTIONS DECLARATIONS:

 void  toBinary(int,char *);
 int  toDecimal(char *);
 void  rComplement(char *);
 void  binaryAdd(char *,char *,char *);

__MATHSUITE static ityp   math_sum(mpfr_t, mpfr_t);
__MATHSUITE static ityp   math_sub(mpfr_t, mpfr_t);
__MATHSUITE static double complex   math_csum(register double complex, register double complex);
__MATHSUITE static double complex   math_csub(register double complex, register double complex);

const struct ext_math_type ext_math =
{
    {
        IDENTIFIER_SIN,
        IDENTIFIER_COS,
        IDENTIFIER_SINH,
        IDENTIFIER_COSH,
        IDENTIFIER_CSC,
        IDENTIFIER_CSCH,
        IDENTIFIER_SEC,
        IDENTIFIER_SECH,
        IDENTIFIER_ASIN,
        IDENTIFIER_ASINH,
        IDENTIFIER_ACOS,
        IDENTIFIER_ACOSH,
        IDENTIFIER_ACSC,
        IDENTIFIER_ACSCH,
        IDENTIFIER_ASEC,
        IDENTIFIER_ASECH,
        IDENTIFIER_TAN,
        IDENTIFIER_TANH,
        IDENTIFIER_ATAN,
        IDENTIFIER_ATANH,
        IDENTIFIER_COT,
        IDENTIFIER_COTH,
        IDENTIFIER_ACOT,
        IDENTIFIER_ACOTH,
        IDENTIFIER_HSIN,
        IDENTIFIER_HSINH,
        IDENTIFIER_QSIN,
        IDENTIFIER_QSINH,
        IDENTIFIER_HCOS,
        IDENTIFIER_HCOSH,
        IDENTIFIER_QCOS,
        IDENTIFIER_QCOSH,
        IDENTIFIER_HSEC,
        IDENTIFIER_HSECH,
        IDENTIFIER_QSEC,
        IDENTIFIER_QSECH,
        IDENTIFIER_HCSC,
        IDENTIFIER_HCSCH,
        IDENTIFIER_QCSC,
        IDENTIFIER_QCSCH,
        IDENTIFIER_HTAN,
        IDENTIFIER_HTANH,
        IDENTIFIER_QTAN,
        IDENTIFIER_QTANH,
        IDENTIFIER_HCOT,
        IDENTIFIER_HCOTH,
        IDENTIFIER_QCOT,
        IDENTIFIER_QCOTH,
        IDENTIFIER_VSIN,
        IDENTIFIER_VSINH,
        IDENTIFIER_CVSIN,
        IDENTIFIER_CVSINH,
        IDENTIFIER_VCOS,
        IDENTIFIER_VCOSH,
        IDENTIFIER_CVCOS,
        IDENTIFIER_CVCOSH,
        IDENTIFIER_HVSIN,
        IDENTIFIER_HVSINH,
        IDENTIFIER_HCVSIN,
        IDENTIFIER_HCVSINH,
        IDENTIFIER_QVSIN,
        IDENTIFIER_QVSINH,
        IDENTIFIER_QCVSIN,
        IDENTIFIER_QCVSINH,
        IDENTIFIER_HVCOS,
        IDENTIFIER_HVCOSH,
        IDENTIFIER_HCVCOS,
        IDENTIFIER_HCVCOSH,
        IDENTIFIER_QVCOS,
        IDENTIFIER_QVCOSH,
        IDENTIFIER_QCVCOS,
        IDENTIFIER_QCVCOSH,
        IDENTIFIER_ESEC,
        IDENTIFIER_ESECH,
        IDENTIFIER_ECSC,
        IDENTIFIER_ECSCH,
        IDENTIFIER_HESEC,
        IDENTIFIER_HESECH,
        IDENTIFIER_HECSC,
        IDENTIFIER_HECSCH,
        IDENTIFIER_QESEC,
        IDENTIFIER_QESECH,
        IDENTIFIER_QECSC,
        IDENTIFIER_QECSCH,
        IDENTIFIER_SINC,
        IDENTIFIER_SINCH,
        IDENTIFIER_HSINC,
        IDENTIFIER_HSINCH,
        IDENTIFIER_QSINC,
        IDENTIFIER_QSINCH,
        IDENTIFIER_COSC,
        IDENTIFIER_COSCH,
        IDENTIFIER_HCOSC,
        IDENTIFIER_HCOSCH,
        IDENTIFIER_QCOSC,
        IDENTIFIER_QCOSCH,
        IDENTIFIER_SECC,
        IDENTIFIER_SECCH,
        IDENTIFIER_HSECC,
        IDENTIFIER_HSECCH,
        IDENTIFIER_QSEC,
        IDENTIFIER_QSECH,
        IDENTIFIER_CSCC,
        IDENTIFIER_CSCCH,
        IDENTIFIER_HCSCC,
        IDENTIFIER_HCSCCH,
        IDENTIFIER_QCSCC,
        IDENTIFIER_QCSCCH,
        IDENTIFIER_TANC,
        IDENTIFIER_TANCH,
        IDENTIFIER_HTANC,
        IDENTIFIER_HTANCH,
        IDENTIFIER_QTANC,
        IDENTIFIER_QTANCH,
        IDENTIFIER_COTC,
        IDENTIFIER_COTCH,
        IDENTIFIER_HCOTC,
        IDENTIFIER_HCOTCH,
        IDENTIFIER_QCOT,
        IDENTIFIER_QCOTH,
        IDENTIFIER_LN,
        IDENTIFIER_LOG,
        IDENTIFIER_LOG2,
        IDENTIFIER_LNC,
        IDENTIFIER_LOGC,
        IDENTIFIER_LOG2C,
        IDENTIFIER_LOG1P,
        IDENTIFIER_LOG1PC,
        IDENTIFIER_EXP,
        IDENTIFIER_EXPC,
        IDENTIFIER_EXP10,
        IDENTIFIER_EXP10C,
        IDENTIFIER_EXP2,
        IDENTIFIER_EXP2C,
        IDENTIFIER_ARMONICSUM,
        IDENTIFIER_FIBONACCI,
        IDENTIFIER_FACTORIAL,
        IDENTIFIER_DOUBLEFACTORIAL,
        IDENTIFIER_NPRIMENUMBER,
        IDENTIFIER_PRIMORIAL,
        IDENTIFIER_FIRSTNPRIMENUMBERSSUM,
        IDENTIFIER_FIBONACCIAL,
        IDENTIFIER_FIBONACCISUM,
        IDENTIFIER_FACTORIALSUM,
        IDENTIFIER_DOUBLEFACTORIALSUM,
        IDENTIFIER_FIRSTNNUMBERSSUM,
        IDENTIFIER_FLOOR,
        IDENTIFIER_CEIL,
        IDENTIFIER_MPFRDEG,
        IDENTIFIER_MPFRRAD
    },
    {
        mswrap_sin,
        mswrap_cos,
        mswrap_sinh,
        mswrap_cosh,
        mswrap_sec,
        mswrap_sech,
        mswrap_csc,
        mswrap_csch,
        mswrap_asin,
        mswrap_asinh,
        mswrap_acos,
        mswrap_acosh,
        mpfr_asec,
        mpfr_asech,
        mpfr_acsc,
        mpfr_acsch,
        mswrap_tan,
        mswrap_tanh,
        mswrap_atan,
        mswrap_atanh,
        mswrap_cot,
        mswrap_coth,
        mpfr_acot,
        mpfr_acoth,
        mpfr_hsin,
        mpfr_hsinh,
        mpfr_qsin,
        mpfr_qsinh,
        mpfr_hcos,
        mpfr_hcosh,
        mpfr_qcos,
        mpfr_qcosh,
        mpfr_hsec,
        mpfr_hsech,
        mpfr_qsec,
        mpfr_qsech,
        mpfr_hcsc,
        mpfr_hcsch,
        mpfr_qcsc,
        mpfr_qcsch,
        mpfr_htan,
        mpfr_htanh,
        mpfr_qtan,
        mpfr_qtanh,
        mpfr_hcot,
        mpfr_hcoth,
        mpfr_qcot,
        mpfr_qcoth,
        mpfr_vsin,
        mpfr_vsinh,
        mpfr_cvsin,
        mpfr_cvsinh,
        mpfr_vcos,
        mpfr_vcosh,
        mpfr_cvcos,
        mpfr_cvcosh,
        mpfr_hvsin,
        mpfr_hvsinh,
        mpfr_hcvsin,
        mpfr_hcvsinh,
        mpfr_qvsin,
        mpfr_qvsinh,
        mpfr_qcvsin,
        mpfr_qcvsinh,
        mpfr_hvcos,
        mpfr_hvcosh,
        mpfr_hcvcos,
        mpfr_hcvcosh,
        mpfr_qvcos,
        mpfr_qvcosh,
        mpfr_qcvcos,
        mpfr_qcvcosh,
        mpfr_esec,
        mpfr_esech,
        mpfr_ecsc,
        mpfr_ecsch,
        mpfr_hesec,
        mpfr_hesech,
        mpfr_hecsc,
        mpfr_hecsch,
        mpfr_qesec,
        mpfr_qesech,
        mpfr_qecsc,
        mpfr_qecsch,
        mpfr_sinc,
        mpfr_sinch,
        mpfr_hsinc,
        mpfr_hsinch,
        mpfr_qsinc,
        mpfr_qsinch,
        mpfr_cosc,
        mpfr_cosch,
        mpfr_hcosc,
        mpfr_hcosch,
        mpfr_qcosc,
        mpfr_qcosch,
        mpfr_secc,
        mpfr_secch,
        mpfr_hsecc,
        mpfr_hsecch,
        mpfr_qsecc,
        mpfr_qsecch,
        mpfr_cscc,
        mpfr_cscch,
        mpfr_hcscc,
        mpfr_hcscch,
        mpfr_qcscc,
        mpfr_qcscch,
        mpfr_tanc,
        mpfr_tanch,
        mpfr_htanc,
        mpfr_htanch,
        mpfr_qtanc,
        mpfr_qtanch,
        mpfr_cotc,
        mpfr_cotch,
        mpfr_hcotc,
        mpfr_hcotch,
        mpfr_qcotc,
        mpfr_qcotch,
        mswrap_log,
        mswrap_log10,
        mswrap_log2,
        mpfr_logc,
        mpfr_log10c,
        mpfr_log2c,
        mswrap_log1p,
        mpfr_log1pc,
        mswrap_exp,
        mpfr_expc,
        mswrap_exp10,
        mpfr_exp10c,
        mswrap_exp2,
        mpfr_exp2c,
        mpfr_asum,
        mpfr_fibo,
        mpfr_fact,
        mpfr_dfact,
        N_prime_Number,
        mpfr_primr,
        mpfr_fpnsum,
        mpfr_fibnc,
        mpfr_fsum,
        mpfr_fasum,
        mpfr_sfasum,
        mpfr_fnnsum,
        mswrap_floor,
        mswrap_ceil,
        mpfr_deg,
        mpfr_rad
    },
    {
        sin,
        cos,
        sinh,
        cosh,
        sec,
        sech,
        csc,
        csch,
        asin,
        asinh,
        acos,
        acosh,
        asec,
        asech,
        acsc,
        acsch,
        tan,
        tanh,
        atan,
        atanh,
        cot,
        coth,
        acot,
        acoth,
        hsin,
        hsinh,
        qsin,
        qsinh,
        hcos,
        hcosh,
        qcos,
        qcosh,
        hsec,
        hsech,
        qsec,
        qsech,
        hcsc,
        hcsch,
        qcsc,
        qcsch,
        htan,
        htanh,
        qtan,
        qtanh,
        hcot,
        hcoth,
        qcot,
        qcoth,
        vsin,
        vsinh,
        cvsin,
        cvsinh,
        vcos,
        vcosh,
        cvcos,
        cvcosh,
        hvsin,
        hvsinh,
        hcvsin,
        hcvsinh,
        qvsin,
        qvsinh,
        qcvsin,
        qcvsinh,
        hvcos,
        hvcosh,
        hcvcos,
        hcvcosh,
        qvcos,
        qvcosh,
        qcvcos,
        qcvcosh,
        esec,
        esech,
        ecsc,
        ecsch,
        hesec,
        hesech,
        hecsc,
        hecsch,
        qesec,
        qesech,
        qecsc,
        qecsch,
        sinc,
        sinch,
        hsinc,
        hsinch,
        qsinc,
        qsinch,
        cosc,
        cosch,
        hcosc,
        hcosch,
        qcosc,
        qcosch,
        secc,
        secch,
        hsecc,
        hsecch,
        qsecc,
        qsecch,
        cscc,
        cscch,
        hcscc,
        hcscch,
        qcscc,
        qcscch,
        tanc,
        tanch,
        htanc,
        htanch,
        qtanc,
        qtanch,
        cotc,
        cotch,
        hcotc,
        hcotch,
        qcotc,
        qcotch,
        log,
        log10,
        log2,
        logc,
        log10c,
        log2c,
        log1p,
        log1pc,
        exp,
        expc,
        exp10,
        exp10c,
        exp2,
        exp2c,
        floor,
        ceil,
        deg,
        rad
    },
        {
        "Sunday",
        "Monday",
        "Tuesday",
        "Wednesday",
        "Thursday",
        "Friday",
        "Saturday"
    },
    {
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December"
    },
    {
        0,
        3,
        2,
        5,
        0,
        3,
        5,
        1,
        4,
        6,
        2,
        4
    },
    {
        "M",
        "MM",
        "MMM"
    },
    {
        {
            "C",
            "CC",
            "CCC",
            "CD",
            "D",
            "DX",
            "DXX",
            "DXXX",
            "CM"
        },
        {
            "X",
            "XX",
            "XXX",
            "XL",
            "L",
            "LX",
            "LXX",
            "LXXX",
            "XC"
        },
        {
             "I",
            "II",
            "III",
            "IV",
            "V",
            "VI",
            "VII",
            "VIII",
            "IX"
        }
    }
};

// COMPLEX FUNCTIONS DEFINITIONS

 void  toBinary(int a,char *c)
    {
    char *s=c;
    int r;
    while(a!=0)
        {
        r=a%2;
        *c='0'+r;
        ++c;
        a*=0.5;
        }

    while((c-s)!=(SIZE))
        {
        *c='0';
        ++c;
        }
    *c='\0';
    strrev(s);
    }

 int  toDecimal(char *c)
    {
    int d=0,count=SIZE-1;
    while(count!=-1)
        {
        if(*(c+count)=='1')
            {
            d+=pow(2,SIZE-count-1);
            }
        count--;
        }
    return d;
    }

 void  rComplement(char *c)
    {

    int count=strlen(c)-1;
    while((*(c+count))!='1')
        {
        count--;
        }
    count--;
    while(count!=-1)
        {
        if((*(c+count))=='1')
            *(c+count)='0';
        else
            *(c+count)='1';
        count--;
        }
    }

 void  binaryAdd(char *first,char *second,char *sum)
    {
    int count=SIZE-1,temp;
    char carry='0';
    *(sum+SIZE)='\0';
    while(count!=-1)
        {
        temp=(*(first+count)-'0')+(*(second+count)-'0')+carry-'0';
        if(temp<2)
            {
            *(sum+count)=temp+'0';
            carry='0';
            }
        else if(temp==2)
            {
            *(sum+count)='0';
            carry='1';
            }
        else
            {
            *(sum+count)='1';
            carry='1';
            }
        count--;
        }
    }

/// thanks to UNKNOWN
/// visit: http://www.Planet-Source-Code.com/vb/scripts/ShowCode.asp?txtCodeId=13861&lngWId=3
/// for further informations about the Original Code from which I taken the Core Part.

__MATHSUITE char * const   binaryAlgSum(ityp a, ityp b, bool mode)
{
    int num1, num2;
    char *c, *d, *e;

    num1 = (int) a;
    num2 = (int) b;

    c = (char *) malloc(SIZE+1);
    d = (char *) malloc(SIZE+1);
    e = (char *) malloc(SIZE+1);

    errMem(c, NULL);
    errMem(c, NULL);
    errMem(c, NULL);

    if(mode)
    {
        if(num1 >= 0)
            toBinary(num1,c);
        else
        {
            toBinary(abs(num1),c);
            rComplement(c);
        }

        if(b >= 0)
            toBinary(num2,d);
        else
        {
            toBinary(abs(num2),d);
            if(b!=-1)
                rComplement(d);
        }

        binaryAdd(c,d,e);
        msprintf(COLOR_USER, "\t %s\t %d\n",c,num1);
        msprintf(COLOR_USER, "\t+%s\t+%d =\n",d,num2);

        if(*(e)=='1')
        {
            msprintf(COLOR_USER, "\t %s\t",e);
            rComplement(e);
            msprintf(COLOR_USER, "-%d",toDecimal(e));
        }
        else
        {
            msprintf(COLOR_USER, "\t %s\t ",e);
            msprintf(COLOR_USER, "%d",toDecimal(e));
        }
    }
    else
    {
        if(a>=0)
            toBinary(a,c);
        else
        {
            toBinary(abs(a),c);
            rComplement(c);
        }
        if(b>=0)
        {
            toBinary(b,d);
            if(b!=0)
                rComplement(d);
        }
        else
            toBinary(abs(b),d);

        binaryAdd(c,d,e);
        msprintf(COLOR_USER, "\t %s\t %d\n",c,a);
        if(b<0)
            msprintf(COLOR_USER, "\t+%s\t+%d\n ",d,abs(b));
        else
            msprintf(COLOR_USER, "\t+%s\t-%d\n ",d,b);

        if(*(e)=='1')
        {
            msprintf(COLOR_USER, "\t %s\t",e);
            rComplement(e);
            msprintf(COLOR_USER, "-%d",toDecimal(e));
        }
        else
        {
            msprintf(COLOR_USER, "\t %s\t ",e);
            msprintf(COLOR_USER, "%d",toDecimal(e));
        }
    }

    PRINTN();
    free(d);
    free(e);
    return c;
}

__MATHSUITE char * const   binNumComp(ityp a)
{
    int num;
    num = (int) a;
    char *c = (char *) malloc(SIZE+1);
    char bin_num[sizeof(c)];
    errMem(c, NULL);
    strcpy(c, NULL_CHAR);
    toBinary(a, c);
    strcpy(bin_num, c);
    rComplement(c);
    msprintf(COLOR_USER, "\nInserted DECIMAL NUMBER: %lld;\nIts BINARY Conversion is: %s.\nIts TWO's COMPLEMENT is: %s.\n", num, bin_num, c);
    toBinary(~num, c);
    msprintf(COLOR_USER, "Decimal ONE'S COMPLEMENT: %d;\nBinary ONE'S COMPLEMENT is: %s.\n\n", ~num, c);
    return c;
}

__MATHSUITE inline void  ___cabs(mpfr_t rop, mpfr_t *restrict cpx, const register sel_typ dim)
{
	const register sel_typ algebra_units = exp2(dim);
	mpfr_set_ui(rop, 0, MPFR_RNDN);
	mpfr_t tmp;
	mpfr_init(tmp);
	for(sel_typ i=0; i<algebra_units; ++i)
	{
		mpfr_exp2(tmp, cpx[i], MPFR_RNDN);
		mpfr_add(rop, rop, tmp, MPFR_RNDN);
	}
	mpfr_sqrt(rop, rop, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  _complexAdd(mpfr_t *restrict cpx, mpfr_t complexRes[static MAX_COMPLEX_UNITS])
{
	const register double complex result = (mpfr_get_d(*(cpx + REAL_PART), MPFR_RNDN) + (mpfr_get_d(*(cpx + IMAG_PART), MPFR_RNDN))*I) + (mpfr_get_d(*(cpx + MAX_COMPLEX_UNITS + REAL_PART), MPFR_RNDN) + (mpfr_get_d(*(cpx + MAX_COMPLEX_UNITS + IMAG_PART), MPFR_RNDN))*I);
	#pragma omp parallel sections num_threads(MAX_COMPLEX_UNITS)
	{
		#pragma omp section
		mpfr_set_d(complexRes[REAL_PART], creal(result), MPFR_RNDN);
		#pragma omp section
    	mpfr_set_d(complexRes[IMAG_PART], cimag(result), MPFR_RNDN);
	}
    return;
}

__MATHSUITE inline void  _complexSub(mpfr_t *restrict cpx, mpfr_t complexRes[static MAX_COMPLEX_UNITS])
{
	const register double complex result = (mpfr_get_d(*(cpx + REAL_PART), MPFR_RNDN) + (mpfr_get_d(*(cpx + IMAG_PART), MPFR_RNDN))*I) - (mpfr_get_d(*(cpx + MAX_COMPLEX_UNITS + REAL_PART), MPFR_RNDN) + (mpfr_get_d(*(cpx + MAX_COMPLEX_UNITS + IMAG_PART), MPFR_RNDN))*I);
	#pragma omp parallel sections num_threads(MAX_COMPLEX_UNITS)
	{
		#pragma omp section
		mpfr_set_d(complexRes[REAL_PART], creal(result), MPFR_RNDN);
    	#pragma omp section
		mpfr_set_d(complexRes[IMAG_PART], cimag(result), MPFR_RNDN);
	}
    return;
}

__MATHSUITE inline void  _complexMul(mpfr_t *restrict cpx, mpfr_t complexRes[static MAX_COMPLEX_UNITS])
{
	const register double complex result = (mpfr_get_d(*(cpx + REAL_PART), MPFR_RNDN) + (mpfr_get_d(*(cpx + IMAG_PART), MPFR_RNDN))*I) * (mpfr_get_d(*(cpx + MAX_COMPLEX_UNITS + REAL_PART), MPFR_RNDN) + (mpfr_get_d(*(cpx + MAX_COMPLEX_UNITS + IMAG_PART), MPFR_RNDN))*I);
	#pragma omp parallel sections num_threads(MAX_COMPLEX_UNITS)
	{
		#pragma omp section
		mpfr_set_d(complexRes[REAL_PART], creal(result), MPFR_RNDN);
	    #pragma omp section
		mpfr_set_d(complexRes[IMAG_PART], cimag(result), MPFR_RNDN);
	}
    return;
}

__MATHSUITE inline void  _complexDiv(mpfr_t *restrict cpx, mpfr_t complexRes[static MAX_COMPLEX_UNITS])
{
	const register double complex result = (mpfr_get_d(*(cpx + REAL_PART), MPFR_RNDN) + (mpfr_get_d(*(cpx + IMAG_PART), MPFR_RNDN))*I) / (mpfr_get_d(*(cpx + MAX_COMPLEX_UNITS + REAL_PART), MPFR_RNDN) + (mpfr_get_d(*(cpx + MAX_COMPLEX_UNITS + IMAG_PART), MPFR_RNDN))*I);
	#pragma omp parallel sections num_threads(MAX_COMPLEX_UNITS)
	{
		#pragma omp section
		mpfr_set_d(complexRes[REAL_PART], creal(result), MPFR_RNDN);
	    #pragma omp section
	    mpfr_set_d(complexRes[IMAG_PART], cimag(result), MPFR_RNDN);
	}
    return;
}

__MATHSUITE inline void  _quaternionsAdd(mpfr_t *restrict quaternions, mpfr_t quaternionsRes[static MAX_QUATERNIONS_UNITS])
{
	#pragma omp parallel sections num_threads(MAX_QUATERNIONS_UNITS)
	{
		#pragma omp section
		mpfr_add(quaternionsRes[QUATERNIONS_REALPART], *(quaternions + QUATERNIONS_REALPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(quaternionsRes[QUATERNIONS_IPART], *(quaternions + QUATERNIONS_IPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(quaternionsRes[QUATERNIONS_JPART], *(quaternions + QUATERNIONS_JPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART), MPFR_RNDN);
	    #pragma omp section
	    mpfr_add(quaternionsRes[QUATERNIONS_KPART], *(quaternions + QUATERNIONS_KPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART), MPFR_RNDN);
	}
    return;
}


__MATHSUITE inline void  _quaternionsSub(mpfr_t *restrict quaternions, mpfr_t quaternionsRes[static MAX_QUATERNIONS_UNITS])
{
	#pragma omp parallel sections num_threads(MAX_QUATERNIONS_UNITS)
	{
		#pragma omp section
		mpfr_sub(quaternionsRes[QUATERNIONS_REALPART], *(quaternions + QUATERNIONS_REALPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(quaternionsRes[QUATERNIONS_IPART], *(quaternions + QUATERNIONS_IPART),  *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(quaternionsRes[QUATERNIONS_JPART], *(quaternions + QUATERNIONS_JPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(quaternionsRes[QUATERNIONS_KPART], *(quaternions + QUATERNIONS_KPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART), MPFR_RNDN);
	}
    return;
}

__MATHSUITE void  _quaternionsMul(mpfr_t *restrict quaternions, mpfr_t quaternionsRes[static MAX_QUATERNIONS_UNITS])
{
	mpfr_t tmp, tmp2, tmp3, tmp4;
	mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
	#pragma omp parallel sections num_threads(MAX_QUATERNIONS_UNITS)
	{
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
			mpfr_mul(tmp, *(quaternions + QUATERNIONS_REALPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART), MPFR_RNDN);
			mpfr_mul(tmp2, *(quaternions + QUATERNIONS_IPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART), MPFR_RNDN);
			mpfr_mul(tmp3, *(quaternions + QUATERNIONS_JPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART), MPFR_RNDN);
			mpfr_mul(tmp4, *(quaternions + QUATERNIONS_KPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART), MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_REALPART], tmp, tmp2, MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_REALPART], quaternionsRes[QUATERNIONS_REALPART], tmp3, MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_REALPART], quaternionsRes[QUATERNIONS_REALPART], tmp4, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
			mpfr_mul(tmp, *(quaternions + QUATERNIONS_REALPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART), MPFR_RNDN);
			mpfr_mul(tmp2, *(quaternions + QUATERNIONS_IPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART), MPFR_RNDN);
			mpfr_mul(tmp3, *(quaternions + QUATERNIONS_JPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART), MPFR_RNDN);
			mpfr_mul(tmp4, *(quaternions + QUATERNIONS_KPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART), MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_IPART], tmp, tmp2, MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_IPART], quaternionsRes[QUATERNIONS_IPART], tmp3, MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_IPART], quaternionsRes[QUATERNIONS_IPART], tmp4, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
			mpfr_mul(tmp, *(quaternions + QUATERNIONS_REALPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART), MPFR_RNDN);
			mpfr_mul(tmp2, *(quaternions + QUATERNIONS_JPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART), MPFR_RNDN);
			mpfr_mul(tmp3, *(quaternions + QUATERNIONS_KPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART), MPFR_RNDN);
			mpfr_mul(tmp4, *(quaternions + QUATERNIONS_IPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART), MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_JPART], tmp, tmp2, MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_JPART], quaternionsRes[QUATERNIONS_JPART], tmp3, MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_JPART], quaternionsRes[QUATERNIONS_JPART], tmp4, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
			mpfr_mul(tmp, *(quaternions + QUATERNIONS_REALPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART), MPFR_RNDN);
			mpfr_mul(tmp2, *(quaternions + QUATERNIONS_KPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART), MPFR_RNDN);
			mpfr_mul(tmp3, *(quaternions + QUATERNIONS_IPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART), MPFR_RNDN);
			mpfr_mul(tmp4, *(quaternions + QUATERNIONS_JPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART), MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_KPART], tmp, tmp2, MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_KPART], quaternionsRes[QUATERNIONS_KPART], tmp3, MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_KPART], quaternionsRes[QUATERNIONS_KPART], tmp4, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
		}
	}
    return;
}

__MATHSUITE void  _quaternionsDiv(mpfr_t *restrict quaternions, mpfr_t quaternionsRes[static MAX_QUATERNIONS_UNITS])
{
	mpfr_t squared_qnorm;
	mpfr_init(squared_qnorm);
	_qabs(squared_qnorm, quaternions + MAX_QUATERNIONS_UNITS);
	mpfr_exp2(squared_qnorm, squared_qnorm, MPFR_RNDN);
	
	#pragma omp parallel sections num_threads(MAX_QUATERNIONS_UNITS)
	{
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
			mpfr_mul(tmp, *(quaternions + QUATERNIONS_REALPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART), MPFR_RNDN);
			mpfr_mul(tmp2, *(quaternions + QUATERNIONS_IPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART), MPFR_RNDN);
			mpfr_mul(tmp3, *(quaternions + QUATERNIONS_JPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART), MPFR_RNDN);
			mpfr_mul(tmp4, *(quaternions + QUATERNIONS_KPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART), MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_REALPART], tmp, tmp2, MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_REALPART], quaternionsRes[QUATERNIONS_REALPART], tmp3, MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_REALPART], quaternionsRes[QUATERNIONS_REALPART], tmp4, MPFR_RNDN);
			mpfr_div(quaternionsRes[QUATERNIONS_REALPART], quaternionsRes[QUATERNIONS_REALPART], squared_qnorm, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
			mpfr_mul(tmp, *(quaternions + QUATERNIONS_REALPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART), MPFR_RNDN);
			mpfr_mul(tmp2, *(quaternions + QUATERNIONS_IPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART), MPFR_RNDN);
			mpfr_mul(tmp3, *(quaternions + QUATERNIONS_JPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART), MPFR_RNDN);
			mpfr_mul(tmp4, *(quaternions + QUATERNIONS_KPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART), MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_IPART], tmp, tmp2, MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_IPART], quaternionsRes[QUATERNIONS_IPART], tmp3, MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_IPART], quaternionsRes[QUATERNIONS_IPART], tmp4, MPFR_RNDN);
			mpfr_div(quaternionsRes[QUATERNIONS_IPART], quaternionsRes[QUATERNIONS_IPART], squared_qnorm, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
			mpfr_mul(tmp, *(quaternions + QUATERNIONS_REALPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART), MPFR_RNDN);
			mpfr_mul(tmp2, *(quaternions + QUATERNIONS_JPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART), MPFR_RNDN);
			mpfr_mul(tmp3, *(quaternions + QUATERNIONS_KPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART), MPFR_RNDN);
			mpfr_mul(tmp4, *(quaternions + QUATERNIONS_IPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART), MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_JPART], tmp, tmp2, MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_JPART], quaternionsRes[QUATERNIONS_JPART], tmp3, MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_JPART], quaternionsRes[QUATERNIONS_JPART], tmp4, MPFR_RNDN);
			mpfr_div(quaternionsRes[QUATERNIONS_JPART], quaternionsRes[QUATERNIONS_JPART], squared_qnorm, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
			mpfr_mul(tmp, *(quaternions + QUATERNIONS_REALPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART), MPFR_RNDN);
			mpfr_mul(tmp2, *(quaternions + QUATERNIONS_KPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART), MPFR_RNDN);
			mpfr_mul(tmp3, *(quaternions + QUATERNIONS_IPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART), MPFR_RNDN);
			mpfr_mul(tmp4, *(quaternions + QUATERNIONS_JPART), *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART), MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_KPART], tmp, tmp2, MPFR_RNDN);
			mpfr_sub(quaternionsRes[QUATERNIONS_KPART], quaternionsRes[QUATERNIONS_KPART], tmp3, MPFR_RNDN);
			mpfr_add(quaternionsRes[QUATERNIONS_KPART], quaternionsRes[QUATERNIONS_KPART], tmp4, MPFR_RNDN);
			mpfr_div(quaternionsRes[QUATERNIONS_KPART], quaternionsRes[QUATERNIONS_KPART], squared_qnorm, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
		}
	}

	mpfr_clear(squared_qnorm); 
    return;
}

__MATHSUITE void  _octonionsAdd(mpfr_t *restrict octonions, mpfr_t octonionsRes[static MAX_OCTONIONS_UNITS])
{
	#pragma omp parallel sections num_threads(MAX_OCTONIONS_UNITS)
	{
		#pragma omp section
		mpfr_add(octonionsRes[OCTONIONS_REALPART], *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(octonionsRes[OCTONIONS_E1PART], *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	    #pragma omp section
	    mpfr_add(octonionsRes[OCTONIONS_E2PART], *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(octonionsRes[OCTONIONS_E3PART], *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(octonionsRes[OCTONIONS_E4PART], *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(octonionsRes[OCTONIONS_E5PART], *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(octonionsRes[OCTONIONS_E6PART], *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(octonionsRes[OCTONIONS_E7PART], *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	}
    return;
}

__MATHSUITE void  _octonionsSub(mpfr_t *restrict octonions, mpfr_t octonionsRes[static MAX_OCTONIONS_UNITS])
{
	#pragma omp parallel sections num_threads(MAX_OCTONIONS_UNITS)
	{
		#pragma omp section
		mpfr_sub(octonionsRes[OCTONIONS_REALPART], *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(octonionsRes[OCTONIONS_E1PART], *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(octonionsRes[OCTONIONS_E2PART], *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(octonionsRes[OCTONIONS_E3PART], *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(octonionsRes[OCTONIONS_E4PART], *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(octonionsRes[OCTONIONS_E5PART], *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(octonionsRes[OCTONIONS_E6PART], *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(octonionsRes[OCTONIONS_E7PART], *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	}
    return;
}

__MATHSUITE void  _octonionsMul(mpfr_t *restrict octonions, mpfr_t octonionsRes[static MAX_OCTONIONS_UNITS])
{
	#pragma omp parallel sections num_threads(MAX_OCTONIONS_UNITS)
	{
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
			mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
			mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
			mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
			mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
			mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
			mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
			mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
			mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
			mpfr_sub(octonionsRes[OCTONIONS_REALPART], tmp, tmp2, MPFR_RNDN);
			mpfr_sub(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp3, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp4, MPFR_RNDN);
			mpfr_sub(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp5, MPFR_RNDN);
			mpfr_sub(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp6, MPFR_RNDN);
			mpfr_sub(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp7, MPFR_RNDN);
			mpfr_sub(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp8, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
		}
		
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E1PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp3, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp4, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp5, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp6, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp7, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp8, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E2PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp3, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp4, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp5, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp6, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp7, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp8, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E3PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp3, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp4, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp5, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp6, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp7, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp8, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E4PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp3, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp4, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp5, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp6, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp7, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp8, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E5PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp3, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp4, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp5, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp6, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp7, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp8, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E6PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp3, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp4, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp5, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp6, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp7, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp8, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
		}
	        
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E7PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp3, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp4, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp5, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp6, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp7, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp8, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	    }
	}
    return;
}

__MATHSUITE void  _octonionsDiv(mpfr_t *restrict octonions, mpfr_t octonionsRes[static MAX_OCTONIONS_UNITS])
{
    mpfr_t squared_onorm;
	mpfr_init(squared_onorm);
	_oabs(squared_onorm, octonions + MAX_OCTONIONS_UNITS);
	mpfr_exp2(squared_onorm, squared_onorm, MPFR_RNDN);

	#pragma omp parallel sections num_threads(MAX_OCTONIONS_UNITS)
	{
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
		    mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
		    mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
		    mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
		    mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
		    mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
		    mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
		    mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
		    mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
		    mpfr_add(octonionsRes[OCTONIONS_REALPART], tmp, tmp2, MPFR_RNDN);
		    mpfr_add(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp3, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp4, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp5, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp6, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp7, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], tmp8, MPFR_RNDN);
	        mpfr_div(octonionsRes[OCTONIONS_REALPART], octonionsRes[OCTONIONS_REALPART], squared_onorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E1PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp3, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp4, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp5, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp6, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp7, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], tmp8, MPFR_RNDN);
	        mpfr_div(octonionsRes[OCTONIONS_E1PART], octonionsRes[OCTONIONS_E1PART], squared_onorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E2PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp3, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp4, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp5, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp6, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp7, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], tmp8, MPFR_RNDN);
	        mpfr_div(octonionsRes[OCTONIONS_E2PART], octonionsRes[OCTONIONS_E2PART], squared_onorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E3PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp3, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp4, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp5, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp6, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp7, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], tmp8, MPFR_RNDN);
	        mpfr_div(octonionsRes[OCTONIONS_E3PART], octonionsRes[OCTONIONS_E3PART], squared_onorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
		}
		
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E4PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp3, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp4, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp5, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp6, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp7, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], tmp8, MPFR_RNDN);
	        mpfr_div(octonionsRes[OCTONIONS_E4PART], octonionsRes[OCTONIONS_E4PART], squared_onorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E5PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp3, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp4, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp5, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp6, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp7, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], tmp8, MPFR_RNDN);
	        mpfr_div(octonionsRes[OCTONIONS_E5PART], octonionsRes[OCTONIONS_E5PART], squared_onorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E6PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART), MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E6PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp3, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp4, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp5, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp6, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp7, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], tmp8, MPFR_RNDN);
	        mpfr_div(octonionsRes[OCTONIONS_E6PART], octonionsRes[OCTONIONS_E6PART], squared_onorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	        mpfr_mul(tmp, *(octonions + OCTONIONS_REALPART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(octonions + OCTONIONS_E1PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(octonions + OCTONIONS_E2PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(octonions + OCTONIONS_E3PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(octonions + OCTONIONS_E4PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(octonions + OCTONIONS_E5PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(octonions + OCTONIONS_E7PART), *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART), MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E7PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp3, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp4, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp5, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp6, MPFR_RNDN);
	        mpfr_add(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp7, MPFR_RNDN);
	        mpfr_sub(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], tmp8, MPFR_RNDN);
	        mpfr_div(octonionsRes[OCTONIONS_E7PART], octonionsRes[OCTONIONS_E7PART], squared_onorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
		}
	}
	mpfr_clear(squared_onorm);
    return;
}

__MATHSUITE void  _sedenionsAdd(mpfr_t *restrict sedenions, mpfr_t sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	#pragma omp parallel sections num_threads(MAX_SEDENIONS_UNITS)
	{
		#pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_REALPART], *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E1PART], *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E2PART], *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E3PART], *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E4PART], *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E5PART], *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E6PART], *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E7PART], *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E8PART], *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E9PART], *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E10PART], *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E11PART], *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E12PART], *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E13PART], *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E14PART], *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_add(sedenionsRes[SEDENIONS_E15PART], *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	}
    return;
}

__MATHSUITE void  _sedenionsSub(mpfr_t *restrict sedenions, mpfr_t sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	#pragma omp parallel sections num_threads(MAX_SEDENIONS_UNITS)
	{
		#pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_REALPART], *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E1PART], *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E2PART], *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E3PART], *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E4PART], *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E5PART], *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E6PART], *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E7PART], *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E8PART], *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E9PART], *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E10PART], *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E11PART], *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E12PART], *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E13PART], *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E14PART], *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	    #pragma omp section
		mpfr_sub(sedenionsRes[SEDENIONS_E15PART], *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	}
    return;
}

__MATHSUITE void  _sedenionsMul(mpfr_t *restrict sedenions, mpfr_t sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	#pragma omp parallel sections num_threads(MAX_SEDENIONS_UNITS)
	{
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
		}

        #pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
			mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	    }

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp5,  *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp16, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp16, MPFR_RNDN);
	    	mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}
	}
    return;
}

__MATHSUITE void  _sedenionsDiv(mpfr_t *restrict sedenions, mpfr_t sedenionsRes[static MAX_SEDENIONS_UNITS])
{
    mpfr_t squared_snorm;
	mpfr_init(squared_snorm);
	_sabs(squared_snorm, sedenions + MAX_SEDENIONS_UNITS);
	mpfr_exp2(squared_snorm, squared_snorm, MPFR_RNDN);

    mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
	mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 

	#pragma omp parallel sections num_threads(MAX_SEDENIONS_UNITS)
	{
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_REALPART], sedenionsRes[SEDENIONS_REALPART], squared_snorm, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}
		
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E1PART], sedenionsRes[SEDENIONS_E1PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E2PART], sedenionsRes[SEDENIONS_E2PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E3PART], sedenionsRes[SEDENIONS_E3PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E4PART], sedenionsRes[SEDENIONS_E4PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E5PART], sedenionsRes[SEDENIONS_E5PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E6PART], sedenionsRes[SEDENIONS_E6PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E7PART], sedenionsRes[SEDENIONS_E7PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}
	
		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp8, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E8PART], sedenionsRes[SEDENIONS_E8PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E9PART], sedenionsRes[SEDENIONS_E9PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E10PART], sedenionsRes[SEDENIONS_E10PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp4, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp12, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E11PART], sedenionsRes[SEDENIONS_E11PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp5,  *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp6, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E12PART], sedenionsRes[SEDENIONS_E12PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp7, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp10, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp15, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E13PART], sedenionsRes[SEDENIONS_E13PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp3, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp5, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp9, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp11, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp13, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp14, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E14PART], sedenionsRes[SEDENIONS_E14PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}

		#pragma omp section
		{
			mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
			mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	        mpfr_mul(tmp, *(sedenions + SEDENIONS_REALPART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART), MPFR_RNDN);
	        mpfr_mul(tmp2, *(sedenions + SEDENIONS_E1PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART), MPFR_RNDN);
	        mpfr_mul(tmp3, *(sedenions + SEDENIONS_E2PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART), MPFR_RNDN);
	        mpfr_mul(tmp4, *(sedenions + SEDENIONS_E3PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART), MPFR_RNDN);
	        mpfr_mul(tmp5, *(sedenions + SEDENIONS_E4PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART), MPFR_RNDN);
	        mpfr_mul(tmp6, *(sedenions + SEDENIONS_E5PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART), MPFR_RNDN);
	        mpfr_mul(tmp7, *(sedenions + SEDENIONS_E6PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART), MPFR_RNDN);
	        mpfr_mul(tmp8, *(sedenions + SEDENIONS_E7PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART), MPFR_RNDN);
	        mpfr_mul(tmp9, *(sedenions + SEDENIONS_E8PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp10, *(sedenions + SEDENIONS_E9PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART), MPFR_RNDN);
	        mpfr_mul(tmp11, *(sedenions + SEDENIONS_E10PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART), MPFR_RNDN);
	        mpfr_mul(tmp12, *(sedenions + SEDENIONS_E11PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART), MPFR_RNDN);
	        mpfr_mul(tmp13, *(sedenions + SEDENIONS_E12PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART), MPFR_RNDN);
	        mpfr_mul(tmp14, *(sedenions + SEDENIONS_E13PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART), MPFR_RNDN);
	        mpfr_mul(tmp15, *(sedenions + SEDENIONS_E14PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART), MPFR_RNDN);
	        mpfr_mul(tmp16, *(sedenions + SEDENIONS_E15PART), *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART), MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], tmp, tmp2, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp3, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp4, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp5, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp6, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp7, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp8, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp9, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp10, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp11, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp12, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp13, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp14, MPFR_RNDN);
	        mpfr_add(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp15, MPFR_RNDN);
	        mpfr_sub(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], tmp16, MPFR_RNDN);
	        mpfr_div(sedenionsRes[SEDENIONS_E15PART], sedenionsRes[SEDENIONS_E15PART], squared_snorm, MPFR_RNDN);
	        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
		}
	}

    mpfr_clear(squared_snorm);
    return;
}

__MATHSUITE bool  _secondGradeEquationSolver(mpfr_t *restrict abc, mpfr_t root[static 2])
{
    mpfr_t dscr; // storage for equation's discriminant
    mpfr_t tmp;
    mpfr_inits(tmp, dscr, NULL); 
	mpfr_mul(dscr, abc[COEFF_B], abc[COEFF_B], MPFR_RNDN);
    mpfr_mul(tmp, abc[COEFF_A], abc[COEFF_C], MPFR_RNDN);
    mpfr_mul_ui(tmp, tmp, 4, MPFR_RNDN);
    mpfr_sub(dscr, dscr, tmp, MPFR_RNDN);

    if(mpfLZero(dscr))
    {
        msprintf(COLOR_SYSTEM, "\nThese Equation has Imaginary ROOTS.\n\n");
        mpfr_clears(tmp, dscr, NULL); 
        return false;
    }

    mpfr_mul_ui(tmp, abc[COEFF_A], 2, MPFR_RNDN);
    mpfr_t tmp2;
    mpfr_init(tmp2);
    mpfr_sqrt(tmp2, dscr, MPFR_RNDN);
    mpfr_t tmp3;
    mpfr_init(tmp3);
    mpfr_neg(tmp3, abc[COEFF_B], MPFR_RNDN);
    mpfr_add(root[ROOT_X1], tmp3, tmp2, MPFR_RNDN);
    mpfr_div(root[ROOT_X1], root[ROOT_X1], tmp, MPFR_RNDN);
    mpfr_sub(root[ROOT_X2], tmp3, tmp2, MPFR_RNDN);
    mpfr_div(root[ROOT_X2], root[ROOT_X2], tmp, MPFR_RNDN);

    mpfr_clears(tmp, tmp2, tmp3, dscr, NULL); 
    return true;
}

// char * getDayName(sel_typ day);
// sel_typ getDayNumber(sel_typ dd, sel_typ mm, uint64_t yy);
//returns the name of the day
/// thanks to Bibek Subedi for these two functions,
/// which I renamed, modified and adapted to this program. Link at:
/// http://www.programming-technique.blogspot.it/2013/02/cc-program-to-ditermine-day-of-week.html
 inline const char * const  getDayName(sel_typ dd)
{
    return ext_math.days_week_names[dd];
}

/// by Me (DekraN)
 inline const char * const  getMonthName(sel_typ mm)
{
    return ext_math.months_names[mm-1];
}

// returns the number of the day correspondent to the given INPUT date
 inline const sel_typ  getDayNumber(sel_typ dd, sel_typ mm, uint64_t yyyy)
{
    return ((yyyy -= mm < 3) + yyyy/4 - yyyy/100 + yyyy/400 + ext_math.getdaynum_map[mm-1] + dd) % MAX_DAYSWEEK;
}

/// thanks to Bibek Subedi for this function,
/// which I renamed, modified and adapted to this program. Link at:
/// http://programming-technique.blogspot.it/2011/08/conversion-of-decimal-to-roman-using-c.html
 inline void  getRomanNumber(dim_typ number, char string[static ROMAN_NUMBER_STRING])
{
    dim_typ i, j, k;
    // k = 0;

    for(i=1000,k=0; i >= 1; i /= 10,++k)
        if((j=((number/i)%(k ? 10 : 1))))
            strcat(string, k ? ext_math.romn_map[k-1][j-1] : ext_math.romn_thousand_map[j-1]);

    return;
}

 void  mpow(mpfr_t rop, mpfr_t b, mpfr_t exp)
{
   if (mpfLZero(exp) < 0)
   {
      mpfr_pow_si(b, b, -1, MPFR_RNDN);
      mpfr_mul_si(exp, exp, -1, MPFR_RNDN);
   }
   mpow2(rop, b, exp);
   return;
}

 void  mpow2(mpfr_t rop, mpfr_t b, mpfr_t exp)
{
    if (mpfr_cmp_ui(exp, 2) > 0)
    {
        mpf_t tmp;
        mpfr_init_set_ui(tmp, 2, MPFR_RNDN);
        mpfr_remainder(tmp, exp, tmp, MPFR_RNDN);
        if(mpfr_zero_p(tmp) == 0)
        {
            mpfr_t tmp2, tmp3;
            mpfr_inits(tmp2, tmp3, NULL); 
            mpfr_mul(tmp2, b, b, MPFR_RNDN);
            mpfr_div_ui(tmp3, exp, 2, MPFR_RNDN);
            mpfr_init(rop);
            mpow(rop, tmp2, tmp3);
            if(mpfr_zero_p(tmp) == 0)
            	mpfr_mul(rop, b, rop, MPFR_RNDN);
            mpfr_clears(tmp, tmp2, tmp3, NULL); 
            return;
        }
        mpfr_clear(tmp);
    }
    else if (mpfr_cmp_ui(exp, 2) == 0)
    {
        mpfr_mul(rop, b, b, MPFR_RNDN);
        return;
    }
   else if (mpfr_cmp_ui(exp, 1) == 0)
   {
      mpfr_set(rop, b, MPFR_RNDN);
      return;
   }

    mpfr_set_ui(rop, 1, MPFR_RNDN);
    return; // exp == 0
}

/// special thanks to Bibek Subedi for this function,
/// which I renamed, adapted and modified to this program. Link at:
/// http://programming-technique.blogspot.it/2013/05/pascals-triangle-using-c-program.html
__MATHSUITE  inline void  getPascalTriangle(mpfr_t rows)
{
	
	mpfr_t tmp;
    mpfr_t i, j;
    char buf[MAX_BUFSIZ];
    
    mpfr_inits(i, j, tmp, NULL); 
    PRINTL();

    for(mpfr_init_set_ui(i,0,MPFR_RNDN); mpfr_less_p(i, rows); mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        for (mpfr_init_set_ui(j,0,MPFR_RNDN); mpfr_less_p(rows, i); mpfr_add_ui(j, j, 1, MPFR_RNDN))
            PRINTSPACE();

        for(mpfr_init_set_ui(j,1,MPFR_RNDN); mpfr_less_p(j,i); mpfr_add_ui(j, j, 1, MPFR_RNDN))
        {
        	comb(tmp, i, j);
            mpfr_sprintf(buf, "%Rf ", tmp);
            msprintf(COLOR_USER, buf);
        }
        PRINTN();
    }

	mpfr_clears(i, j, tmp, NULL); 
    PRINTL();
    return;
}

/// PRIME NUMBERS FUNCTIONS
/// thanks to Bibek Subedi for these functions,
/// that I renamed, modified and adapted to this program both each. Link at:
/// http://programming-technique.blogspot.it/2012/02/fastest-way-of-calculating-prime-number.html

/// HISP isPrime bool checking function
/// (High Iteration - Slow Process)
 inline bool   isPrimeHISP(register mpfr_t n)
{
    mpfr_t i, result;
    mpfr_init(result);
    mpfr_init_set_ui(i, 2, MPFR_RNDN);
    do
    {
        mpfr_remainder(result, n, i, MPFR_RNDN);
        if(mpfr_zero_p(result))
        {
            mpfr_clears(i, result, NULL); 
            return false;
        }
        mpfr_add_ui(i, i, 1, MPFR_RNDN);
    }
    while(mpfr_cmp(n, i));
    mpfr_clear(result);
    if(mpfr_cmp(n, i) == 0)
    {
        mpfr_clear(i);
        return true;
    }
    mpfr_clear(i);
    return false;
}

/// LIFP isPrime bool checking function
/// (Low Iterations, Fast process)
  bool   isPrimeLIFP(register mpfr_t n)
{
    mpfr_t i, tmp;
    mpfr_inits(i, tmp, NULL); 
    mpfr_set_ui(tmp, 2, MPFR_RNDN);
    // this because in my program The user won't be able to insert
    // a n minor of 3.
    // if(n == 2) return false; // true;
    mpfr_remainder(i, n, tmp, MPFR_RNDN);
    if(mpfr_zero_p(i))
    {
        mpfr_clears(i, tmp, NULL); 
        return false;
    }
    /*
    # pragma omp parallel \
        shared ( n ) \
        private ( i )

    # pragma omp for
    */
    mpfr_t tmp2;
    mpfr_init(tmp2);
    for(mpfr_set_ui(i, 2, MPFR_RNDN), mpfr_mul(tmp, i, i, MPFR_RNDN); mpfr_cmp(tmp, n) <= 0; mpfr_add_ui(i, i, 2, MPFR_RNDN), mpfr_mul(tmp, i, i, MPFR_RNDN))
    {
        mpfr_remainder(tmp2, n, i, MPFR_RNDN);
        if(mpfr_zero_p(tmp2))
        {
            mpfr_clears(i, tmp, tmp2, NULL); 
            return false;
        }
    }
    mpfr_clears(i, tmp, tmp2, NULL); 
    return true;
}

__MATHSUITE void   prime_N_Number(register mpfr_t m, register mpfr_t n)
{
    mpfr_t i;
    mpfr_t dim;
	char buf[MAX_BUFSIZ];
	
    mpfr_init_set_ui(dim, 0, MPFR_RNDN);
    // tab = malloc(sizeof)
    SHOWPAUSEMESSAGE();
    PRINTL();

    bool (* const prime_check_func)(register mpfr_t) = lazy_exec ? isPrimeHISP : isPrimeLIFP;
    //
    PRINTL();

    if(isnSett(lazy_exec) && mpfr_cmp_ui(m,2) == 0)
    {
    	mpfr_add_ui(dim, dim, 1, MPFR_RNDN);
        msprintf(COLOR_USER, "- 1: 2;\n");
        mpfr_add_ui(m, m, 1, MPFR_RNDN);
    }

    for(mpfr_init_set(i, m, MPFR_RNDN); mpfr_cmp(i, n) <= 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
        if(prime_check_func(i) || mpfr_cmp_ui(m,2) == 0)
        {
        	mpfr_add_ui(dim, dim, 1, MPFR_RNDN);
            mpfr_sprintf(buf, "- %Rf: %Rf;\n", dim, i);
            msprintf(COLOR_USER, buf);
            if(catchPause())
            {
                mpfr_clears(i, dim, NULL); 
                return;
            }
    	}

    PRINTL();
    msprintf(COLOR_USER, "\n%llu PRIME NUMBERS has been correctly printed.\n\n", dim);
    mpfr_clears(i, dim, NULL); 
    return;
}

// N_prime_Number
__MATHSUITE void   N_prime_Number(mpfr_t rop, register mpfr_t n)
{
    mpfr_t i, dim;
    mpfr_t isn;
    bool (*prime_check_func)(register mpfr_t) = lazy_exec ? isPrimeHISP : isPrimeLIFP;

    mpfr_init_set(isn, n, MPFR_RNDN);
    mpfr_init_set_ui(dim, 0, MPFR_RNDN);

    if(isnSett(lazy_exec))
    {
        if(mpfr_cmp_ui(n,1) == 0)
        {
            mpfr_set_ui(rop, 2, MPFR_RNDN);
            mpfr_clears(isn, dim, NULL); 
            return;
        }
        mpfr_sub_ui(isn, isn, 1, MPFR_RNDN);
    }

    if(mpfGZero(n))
        for(mpfr_init_set_ui(i, 2, MPFR_RNDN); ; mpfr_add_ui(i, i, 1, MPFR_RNDN))
            if(prime_check_func(i))
            {
                mpfr_add_ui(dim, dim, 1, MPFR_RNDN);
                if(mpfr_cmp(dim, isn) == 0)
                {
                    mpfr_set(rop, i, MPFR_RNDN);
                    mpfr_clears(i, isn, dim, NULL); 
                    return;
                }
            }

    mpfr_set_d(rop, INVALIDRETURNVALUE_NPRIMENUMBER, MPFR_RNDN);
    mpfr_clears(isn, i, NULL); 
    return;// basic returning 0 in case of error (obviously not a pn)
}

// Fibonaccial
__MATHSUITE inline void  mpfr_fibnc(mpfr_t rop, register mpfr_t n)
{
    mpfr_t c, i;
    mpfr_init(c);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN), mpfr_init_set(c, i, MPFR_RNDN); mpfr_cmp(i, n) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_fibo(rop, i);
        mpfr_mul(c, c, rop, MPFR_RNDN);
    }

    mpfr_set(rop, c, MPFR_RNDN);
    mpfr_clears(i, c, NULL); 
    return;
}

// N_primorial_Number
__MATHSUITE inline void  mpfr_primr(mpfr_t rop, register mpfr_t n)
{
    mpfr_t c, i;
    mpfr_init_set_ui(c, 1, MPFR_RNDN);

    for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, n) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        N_prime_Number(rop, i);
        mpfr_mul(c, c, rop, MPFR_RNDN);
    }

    mpfr_set(rop, c, MPFR_RNDN);
    mpfr_clears(i, c, NULL); 
    return;
}

// First N Prime Number Sum
__MATHSUITE inline void  mpfr_fpnsum(mpfr_t rop, register mpfr_t n)
{
    mpfr_t c, i;
    mpfr_init_set_ui(c, 0, MPFR_RNDN);

    for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, n) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        N_prime_Number(rop, i);
        mpfr_add(c, c, rop, MPFR_RNDN);
    }

    mpfr_set(rop, c, MPFR_RNDN);
    mpfr_clears(i, c, NULL); 
    return;
}

 //compares if the ityp f1 is equal with f2 and returns 1 if true and 0 if false
 static bool compare_double(register mpfr_t f1, register mpfr_t f2)
 {
    mpfr_t tmp, tmp2;
    mpfr_inits(tmp, tmp2, NULL); 
    mpfr_sub_d(tmp, f1, COMPAREDOUBLE_PRECISION, MPFR_RNDN);
    mpfr_add_d(tmp2, f1, COMPAREDOUBLE_PRECISION, MPFR_RNDN);
    if(mpfr_cmp(tmp, f2) < 0 && mpfr_cmp(tmp2, f2) > 0)
    {
        mpfr_clears(tmp, tmp2, NULL); 
        return true;
    }
	mpfr_clears(tmp, tmp2, NULL); 
	return false;
 }

__MATHSUITE bool   trigonometric_domain(register mpfr_t x, register mpfr_t y, register mpfr_t z)
{
    mpfr_t i;
    mpfr_t tmp;
    mpfr_init(tmp);
    for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, x) <= 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_mul(tmp, i, z, MPFR_RNDN);
        mpfr_add(tmp, tmp, y, MPFR_RNDN);
    	if(compare_double(x,tmp))
        {
            mpfr_clears(i, tmp, NULL); 
    		return false;
        }
    }

    mpfr_clears(i, tmp, NULL); 
    return true;
}

__MATHSUITE inline void  stirling(mpfr_t rop, register mpfr_t n)
{
    mpfr_t tmp, tmp2;
    mpfr_inits(tmp, tmp2, NULL); 
    mpfr_div_d(tmp, n, M_E, MPFR_RNDN);
    mpfr_pow(tmp, tmp, n, MPFR_RNDN);
    mpfr_set_ui(rop, 2, MPFR_RNDN);
    mpfr_mul(rop, rop, n, MPFR_RNDN);
    mpfr_const_pi(tmp2, MPFR_RNDN);
    mpfr_mul(rop, rop, tmp2, MPFR_RNDN);
    mpfr_sqrt(rop, rop, MPFR_RNDN);
    mpfr_mul(rop, rop, tmp, MPFR_RNDN);
    mpfr_clears(tmp, tmp2, NULL); 
    return;
}

__MATHSUITE void  mpfr_fibo(mpfr_t rop, register mpfr_t n)
{
    if(mpfLeZero(n))
        return;

    mpfr_t tmp;
    mpfr_init_set_ui(tmp, 0, MPFR_RNDN);

    if(!access(sysMem)[FUNCTION_FIBONACCI].current_max_index)
    {
    	if(!matrixAlloc(&access(sysMem)[FUNCTION_FIBONACCI].memoizer, (dim_typ2){2,1}))
    	{
    		mpfr_set_d(rop, INVALIDRETURNVALUE_FIBONACCI, MPFR_RNDN);
    		return;
    	}
        mpfr_init_set_ui(access(sysMem)[FUNCTION_FIBONACCI].memoizer[0], 0, MPFR_RNDN);
        mpfr_init_set_ui(access(sysMem)[FUNCTION_FIBONACCI].memoizer[1], 1, MPFR_RNDN);
    }

    mpfr_t prev, curr, next, i;

    mpfr_init(i);
    mpfr_init_set(prev, tmp, MPFR_RNDN);
    mpfr_init_set(next, tmp, MPFR_RNDN);
    mpfr_init_set_ui(curr, 1, MPFR_RNDN);

    if(mpfr_cmp_ui(n, access(sysMem)[FUNCTION_FIBONACCI].current_max_index) < 0)
        mpfr_set(curr, access(sysMem)[FUNCTION_FIBONACCI].memoizer[mpfr_get_ui(n, MPFR_RNDN)], MPFR_RNDN);
    else
    {
        if(access(sysMem)[FUNCTION_FIBONACCI].current_max_index < access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI])
        {

            if(access(sysMem)[FUNCTION_FIBONACCI].current_max_index)
            {
                mpfr_set(prev, access(sysMem)[FUNCTION_FIBONACCI].memoizer[access(sysMem)[FUNCTION_FIBONACCI].current_max_index -1], MPFR_RNDN);
                mpfr_set(curr, access(sysMem)[FUNCTION_FIBONACCI].memoizer[access(sysMem)[FUNCTION_FIBONACCI].current_max_index], MPFR_RNDN);
            }

            for (mpfr_set_ui(i, (access(sysMem)[FUNCTION_FIBONACCI].current_max_index ? access(sysMem)[FUNCTION_FIBONACCI].current_max_index+1 : 2), MPFR_RNDN); mpfr_cmp(i, n) <= 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
            {
                mpfr_add(next, curr, prev, MPFR_RNDN);
                mpfr_set(prev, curr, MPFR_RNDN);

                access(sysMem)[FUNCTION_FIBONACCI].memoizer = realloc(access(sysMem)[FUNCTION_FIBONACCI].memoizer, sizeof(mpfr_t)*(mpfr_get_ui(i, MPFR_RNDN)+1));
				errMem(access(sysMem)[FUNCTION_FIBONACCI].memoizer, (void) mpfr_set_d(rop, INVALIDRETURNVALUE_FIBONACCI, MPFR_RNDN));
                mpfr_set(curr, next, MPFR_RNDN);
                mpfr_init_set(access(sysMem)[FUNCTION_FIBONACCI].memoizer[mpfr_get_ui(i, MPFR_RNDN)], curr, MPFR_RNDN);
            }
            access(sysMem)[FUNCTION_FIBONACCI].current_max_index = mpfr_get_ui(n, MPFR_RNDN);
        }
        else
            for (mpfr_set_ui(i, 2, MPFR_RNDN); mpfr_cmp(i, n) <= 0; )
            {
                mpfr_add(next, curr, prev, MPFR_RNDN);
                mpfr_set(prev, curr, MPFR_RNDN);
                mpfr_set(curr, next, MPFR_RNDN);
                mpfr_add_ui(i, i, 1, MPFR_RNDN);
            }
    }

    if(mpfr_cmp_ui(curr, ULONG_MAX))
        mpfr_set_d(curr, INVALIDRETURNVALUE_FIBONACCI, MPFR_RNDN);

    mpfr_set(rop, curr, MPFR_RNDN);
    mpfr_clears(i, next, prev, tmp, curr, NULL); 
    return;
}

__MATHSUITE inline void  mpfr_fact(mpfr_t rop, register mpfr_t n)
{
	mpfr_fac_ui(rop, mpfr_get_ui(n, MPFR_RNDN), MPFR_RNDN);
    return;
}

/*
    if(n < 0.00) return 0.00;

    const uint64_t num = (uint64_t) n;
    uint64_t i, res;

    if(!accessINVALIDRETURNVALUE_EVEN_DOUBLEFACTORIAL(sysMem)[FUNCTION_FACTORIAL].current_max_index)
    {
        access(sysMem)[FUNCTION_FACTORIAL].memoizer = malloc(sizeof(ityp)<<1);
        errMem(access(sysMem)[FUNCTION_FACTORIAL].memoizer, INVALIDRETURNVALUE_FACTORIAL);
        access(sysMem)[FUNCTION_FACTORIAL].memoizer[0] = access(sysMem)[FUNCTION_FACTORIAL].memoizer[1] = 1;
    }

    res = 1;

    if(num < access(sysMem)[FUNCTION_FACTORIAL].current_max_index)
        res = access(sysMem)[FUNCTION_FACTORIAL].memoizer[num];
    else
    {
        // this cycle redundance has been introduced in order
        // to avoid overhead problems.
        if(access(sysMem)[FUNCTION_FACTORIAL].current_max_index < access(curLayout)->max_memoizable_indices[FUNCTION_FACTORIAL])
        {
            if(access(sysMem)[FUNCTION_FACTORIAL].current_max_index)
                res = access(sysMem)[FUNCTION_FACTORIAL].memoizer[access(sysMem)[FUNCTION_FACTORIAL].current_max_index];

            for(i=(access(sysMem)[FUNCTION_FACTORIAL].current_max_index ? access(sysMem)[FUNCTION_FACTORIAL].current_max_index+1 : 2); i<=num; ++i)
            {
                access(sysMem)[FUNCTION_FACTORIAL].memoizer = realloc(access(sysMem)[FUNCTION_FACTORIAL].memoizer, sizeof(ityp)*(i+1));
                errMem(access(sysMem)[FUNCTION_FACTORIAL].memoizer, res);
            	access(sysMem)[FUNCTION_FACTORIAL].memoizer[i] = i<access(curLayout)->min_stirling_number ? (res *= i) : (res=stirling(i));
            }

            access(sysMem)[FUNCTION_FACTORIAL].current_max_index = num;
        }
        else
        {
			if(num<access(curLayout)->min_stirling_number)
	            for(i=1; ++i<=num; )
	                res *= i;
	        else
	        	res = stirling(num);
		}
    }

    return res < ULLONG_MAX ? res : INVALIDRETURNVALUE_FACTORIAL;
}
*/

__MATHSUITE inline void  mpfr_dfact(mpfr_t rop, register mpfr_t n)
{
	mpfr_set_ui(rop, 2, MPFR_RNDN);
    mpfr_remainder(rop, n, rop, MPFR_RNDN);
    if(mpfr_zero_p(rop) == 0)
    {
        sfact_odd(rop, n);
        return;
    }
    sfact_even(rop, n);
    return;
}

__MATHSUITE void  sfact_even(mpfr_t rop, register mpfr_t n)
{
	if(mpfLZero(n))
        return;

    mpfr_t tmp;
    mpfr_init_set_ui(tmp, 0, MPFR_RNDN);

    if(!access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].current_max_index)
    {
    	if(!matrixAlloc(&access(sysMem)[FUNCTION_FIBONACCI].memoizer, (dim_typ2){2,1}))
    	{
    		mpfr_set_d(rop, INVALIDRETURNVALUE_EVEN_DOUBLEFACTORIAL, MPFR_RNDN);
    		return;
    	}
        mpfr_init_set_ui(access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].memoizer[0], 1, MPFR_RNDN);
        mpfr_init_set_ui(access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].memoizer[1], 1, MPFR_RNDN);
    }

    mpfr_t i, res;
    mpfr_init(i);
    mpfr_init_set_ui(res, 1, MPFR_RNDN);

    if(mpfr_cmp_ui(n, access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].current_max_index) < 0)
        mpfr_set(res, access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].memoizer[mpfr_get_ui(n, MPFR_RNDN)], MPFR_RNDN);
    else
    {
    	mpfr_t tmp2;
        mpfr_init(tmp2);
        // this cycle redundance has been introduced in order
        // to avoid overhead problems.
        if(access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].current_max_index < access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_DOUBLEFACTORIAL])
        {
            if(access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].current_max_index)
                mpfr_set(res, access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].memoizer[access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].current_max_index], MPFR_RNDN);

            for(mpfr_set_ui(i,(access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].current_max_index ? access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].current_max_index+2 : 2), MPFR_RNDN); mpfr_cmp(i, n) <= 0; mpfr_add_ui(i, i, 2, MPFR_RNDN))
            {
                access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].memoizer = realloc(access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].memoizer, sizeof(mpfr_t)*(mpfr_get_ui(i, MPFR_RNDN)+1));
                errMem(access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].memoizer, (void) mpfr_set_d(rop, INVALIDRETURNVALUE_EVEN_DOUBLEFACTORIAL, MPFR_RNDN));
                if(access(curLayout)->min_stirling_number%2 || mpfr_cmp_ui(n, (access(curLayout)->min_stirling_number<<1)) < 0)
                    mpfr_mul(res, res, i, MPFR_RNDN);
                else
                {
                    mpfr_exp2(res, n, MPFR_RNDN);
                    stirling(tmp2, n);
                    mpfr_mul(res, res, tmp2, MPFR_RNDN);
                    mpfr_clear(tmp2);
                }
                mpfr_init_set(access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].memoizer[mpfr_get_ui(i, MPFR_RNDN)], res, MPFR_RNDN);
            }
            access(sysMem)[FUNCTION_EVEN_DOUBLEFACTORIAL].current_max_index = mpfr_get_ui(n, MPFR_RNDN);
        }
        else
    	{
    		if(access(curLayout)->min_stirling_number%2 || mpfr_cmp_ui(n, access(curLayout)->min_stirling_number<<1) < 0)
	            for(mpfr_set_ui(i, 2, MPFR_RNDN); mpfr_cmp(i, n) <= 0; mpfr_add_ui(i, i, 2, MPFR_RNDN))
	                mpfr_mul(res, res, i, MPFR_RNDN);
	        else
            {
                mpfr_exp2(res, n, MPFR_RNDN);
                stirling(tmp2, n);
                mpfr_mul(res, res, tmp2, MPFR_RNDN);
                mpfr_clear(tmp2);
            }
	    }
    }

    if(mpfr_cmp_ui(res, ULONG_MAX) >= 0)
        mpfr_set_d(res, INVALIDRETURNVALUE_EVEN_DOUBLEFACTORIAL, MPFR_RNDN);

    mpfr_set(rop, res, MPFR_RNDN);
    mpfr_clears(i, tmp, res, NULL); 
    return;
}

__MATHSUITE void  sfact_odd(mpfr_t rop, register mpfr_t n)
{
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, 0, MPFR_RNDN);

    if(mpfLZero(n))
        return;

    if(!access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].current_max_index)
    {
    	if(!matrixAlloc(&access(sysMem)[FUNCTION_FIBONACCI].memoizer, (dim_typ2){3,1}))
    	{
    		mpfr_set_d(rop, INVALIDRETURNVALUE_ODD_DOUBLEFACTORIAL, MPFR_RNDN);
    		return;
    	}
        mpfr_init_set_ui(access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].memoizer[0], 1, MPFR_RNDN);
        mpfr_init_set_ui(access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].memoizer[1], 1, MPFR_RNDN);
        mpfr_init_set_ui(access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].memoizer[2], 2, MPFR_RNDN);
    }

    mpfr_t i, res;
    mpfr_init(i);
    mpfr_init_set_ui(res, 1, MPFR_RNDN);

    if(mpfr_cmp_ui(n, access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].current_max_index) < 0)
        mpfr_set(res, access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].memoizer[mpfr_get_ui(n, MPFR_RNDN)], MPFR_RNDN);
    else
    {
        // this cycle redundance has been introduced in order
        // to avoid overhead problems.
        if(access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].current_max_index < access(curLayout)->max_memoizable_indices[FUNCTION_ODD_DOUBLEFACTORIAL])
        {
            if(access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].current_max_index)
                mpfr_set(res, access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].memoizer[access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].current_max_index], MPFR_RNDN);

            for(mpfr_set_ui(i,(access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].current_max_index ? access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].current_max_index+2 : 3), MPFR_RNDN); mpfr_cmp(i, n) <= 0; mpfr_add_ui(i, i, 2, MPFR_RNDN))
            {
                access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].memoizer = realloc(access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].memoizer, sizeof(mpfr_t)*(mpfr_get_ui(i, MPFR_RNDN)+1));
                errMem(access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].memoizer, (void) mpfr_set_d(rop, INVALIDRETURNVALUE_EVEN_DOUBLEFACTORIAL, MPFR_RNDN));
                mpfr_mul(res, res, i, MPFR_RNDN);
                mpfr_init_set(access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].memoizer[mpfr_get_ui(i, MPFR_RNDN)], res, MPFR_RNDN);
            }

            access(sysMem)[FUNCTION_ODD_DOUBLEFACTORIAL].current_max_index = mpfr_get_ui(n, MPFR_RNDN);
        }
        else
            for(mpfr_set_ui(i, 1, MPFR_RNDN); mpfr_cmp(i, n) <= 0; mpfr_add_ui(i, i, 2, MPFR_RNDN))
                mpfr_mul(res, res, i, MPFR_RNDN);
    }

    if(mpfr_cmp_ui(res, ULONG_MAX) >= 0)
        mpfr_set_d(res, INVALIDRETURNVALUE_ODD_DOUBLEFACTORIAL, MPFR_RNDN);

    mpfr_set(rop, res, MPFR_RNDN);
    mpfr_clears(i, tmp, res, NULL); 
    return;
}

 inline void  perm(mpfr_t rop, register mpfr_t n)
{
    mpfr_fac_ui(rop, mpfr_get_ui(n, MPFR_RNDN), MPFR_RNDN);
    return;
}

 inline void  perm_rep(mpfr_t rop, register mpfr_t n, register mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t tmp;
    mpfr_init(tmp);
    mpfr_fac_ui(rop, mpfr_get_ui(n, MPFR_RNDN), MPFR_RNDN);
    mpfr_product(tmp, dim, PRODUCTORY_MUL, vector);
    mpfr_div(rop, rop, tmp, MPFR_RNDN);
    mpfr_clear(tmp);
    return;
}

 inline void  kperm(mpfr_t rop, register mpfr_t n, register mpfr_t k)
{
    mpfr_t tmp, tmp2, tmp3;
    mpfr_init(tmp);
    mpfr_sub(tmp, n, k, MPFR_RNDN);
    mpfr_fac_ui(tmp2, mpfr_get_ui(n, MPFR_RNDN), MPFR_RNDN);
    mpfr_fac_ui(tmp3, mpfr_get_ui(tmp, MPFR_RNDN), MPFR_RNDN);
    mpfr_div(tmp, tmp2, tmp3, MPFR_RNDN);
    mpfr_set(rop, tmp, MPFR_RNDN);
    mpfr_clears(tmp, tmp2, tmp3, NULL); 
    return;
}

 inline void  kperm_rep(mpfr_t rop, register mpfr_t n, register mpfr_t k)
{
    mpfr_pow(rop, n, k, MPFR_RNDN);
    return;
}

 inline void  comb(mpfr_t rop, register mpfr_t n, register mpfr_t r)
{
    mpfr_t tmp, tmp2, tmp3;
    mpfr_inits(tmp, tmp2, tmp3, NULL); 
    mpfr_sub(tmp, n, r, MPFR_RNDN);
    mpfr_fac_ui(tmp2, mpfr_get_ui(tmp, MPFR_RNDN), MPFR_RNDN);
    mpfr_fac_ui(tmp3, mpfr_get_ui(r, MPFR_RNDN), MPFR_RNDN);
    mpfr_mul(rop, tmp2, tmp3, MPFR_RNDN);
    mpfr_fac_ui(tmp, mpfr_get_ui(n, MPFR_RNDN), MPFR_RNDN);
    mpfr_div(rop, tmp, rop, MPFR_RNDN);
    mpfr_clears(tmp, tmp2, tmp3, NULL); 
    return;
}

 inline void  comb_rep(mpfr_t rop, register mpfr_t n, register mpfr_t r)
{
    mpfr_t tmp, tmp2, tmp3, tmp4;
    mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
    mpfr_sub_ui(tmp, n, 1, MPFR_RNDN);
    mpfr_add(tmp2, n, r, MPFR_RNDN);
    mpfr_sub_ui(tmp2, tmp2, 1, MPFR_RNDN);
    mpfr_fac_ui(tmp3, mpfr_get_ui(tmp2, MPFR_RNDN), MPFR_RNDN);
    mpfr_fac_ui(tmp4, mpfr_get_ui(tmp, MPFR_RNDN), MPFR_RNDN);
    mpfr_div(rop, tmp3, tmp4, MPFR_RNDN);
    mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
    return;
}

// Somma SUCCESSIONE GEOMETRICA
__MATHSUITE inline void  mpfr_gsum(mpfr_t rop, register mpfr_t a, register mpfr_t exponent)
{
    mpfr_t tmp, tmp2, tmp3;
    mpfr_inits(tmp, tmp2, tmp3, NULL); 
    mpfr_add_ui(tmp, exponent, 1, MPFR_RNDN);
    mpfr_set(tmp2, a, MPFR_RNDN);
    mpfr_neg(tmp2, tmp2, MPFR_RNDN);
    mpfr_add_ui(tmp2, tmp2, 1, MPFR_RNDN);
    mpow2(tmp3, a, tmp);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div(rop, rop, tmp2, MPFR_RNDN);
    mpfr_clears(tmp, tmp2, tmp3, NULL); 
    return;
}

// Somma SUCCESSIONE ARMONICA GENERALIZZATA
__MATHSUITE void  mpfr_gasum(mpfr_t rop, register mpfr_t a, register mpfr_t exponent)
{
    mpfr_t c, i, tmp;
    mpfr_init_set_ui(c, 0, MPFR_RNDN);
    mpfr_init(tmp);

    for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, a) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_add_ui(tmp, i, 1, MPFR_RNDN);
        mpow2(rop, tmp, exponent);
        mpfr_pow_si(tmp, rop, -1, MPFR_RNDN);
        mpfr_add(c, c, tmp, MPFR_RNDN);
    }

    mpfr_set(rop, c, MPFR_RNDN);
    mpfr_clears(c, i, NULL); 
    return;
}

// Somma SUCCESSIONE ARMONICA
__MATHSUITE inline void  mpfr_asum(mpfr_t rop, register mpfr_t a)
{
    mpfr_t tmp;
    mpfr_init_set_ui(tmp, 1, MPFR_RNDN);
    mpfr_gasum(rop, a, tmp);
    mpfr_clear(tmp);
    return;
}

// Somma SUCCESSIONE FIBONACCI
__MATHSUITE void  mpfr_fsum(mpfr_t rop, register mpfr_t a)
{
    mpfr_t c, i;
    mpfr_init_set_ui(c, 0, MPFR_RNDN);

    for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, a) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_fibo(rop, i);
        mpfr_add(c, c, rop, MPFR_RNDN);
    }

    mpfr_set(rop, c, MPFR_RNDN);
    mpfr_clears(c, i, NULL); 
    return;
}

// Somma SUCCESSIONE FATTORIALE
__MATHSUITE void  mpfr_fasum(mpfr_t rop, register mpfr_t a)
{
    mpfr_t c, i;
    mpfr_init_set_ui(c, 0, MPFR_RNDN);

    for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, a) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_fac_ui(rop, mpfr_get_ui(i, MPFR_RNDN), MPFR_RNDN);
        mpfr_add(c, c, rop, MPFR_RNDN);
    }

    mpfr_set(rop, c, MPFR_RNDN);
    mpfr_clears(c, i, NULL); 
    return;
}

// Somma SUCCESSIONE SEMIFATTORIALE
__MATHSUITE void  mpfr_sfasum(mpfr_t rop, register mpfr_t a)
{
    mpfr_t c, i;
    mpfr_init_set_ui(c, 0, MPFR_RNDN);

    for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, a) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_dfact(rop, i);
        mpfr_add(c, c, rop, MPFR_RNDN);
    }

    mpfr_set(rop, c, MPFR_RNDN);
    mpfr_clears(c, i, NULL); 
    return;
}

// Somma PRIMI N NUMERI NATURALI
__MATHSUITE inline void  mpfr_fnnsum(mpfr_t rop, register mpfr_t a)
{
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_mul(rop, rop, rop, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return ;
}

// The Answer To Life The Universe And Everything
__MATHSUITE inline void   answer_to_life_the_universe_and_everything(mpfr_t rop)
{
	mpfr_set_ui(rop, ANSWER_TO_LIFE_THE_UNIVERSE_AND_EVERYTHING, MPFR_RNDN);
	return;
}

// SOMMATORIA
__MATHSUITE inline ityp  summation(uint64_t dim, bool mode, ityp vector[static dim])
{
    ityp res = 0.00;
    const register sel_typ mode_binder = 1-(mode<<1);

	#pragma omp parallel for
    for(uint64_t i=0; i<dim; ++i)
        res += (vector[i])*mode_binder;

    return res;
}

// PRODUTTORIA
__MATHSUITE inline ityp  productory(uint64_t dim, bool mode, ityp vector[static dim])
{
    ityp res = 1.00;
    const register sel_typ mode_binder = 1-(mode<<1);

	#pragma omp parallel for
    for(uint64_t i=0; i<dim; ++i)
        res *= (vector[i])*mode_binder;

    return res;
}

// SOMMATORIA
__MATHSUITE void  mpfr_summation(mpfr_t rop, mpfr_t dim, bool mode, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t i, tmp, mode_binder;
    mpfr_set_ui(rop, 0, MPFR_RNDN);
    mpfr_inits(mode_binder, tmp, NULL); 
    mpfr_set_ui(mode_binder, 2*mode, MPFR_RNDN);
    mpfr_neg(mode_binder, mode_binder, MPFR_RNDN);
    mpfr_add_ui(mode_binder, mode_binder, 1, MPFR_RNDN);

    for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_pow(tmp, vector[mpfr_get_ui(i, MPFR_RNDN)], mode_binder, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_clears(i, tmp, NULL); 
    return;
}

// PRODUTTORIA
__MATHSUITE void  mpfr_product(mpfr_t rop, mpfr_t dim, bool mode, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t i, tmp, mode_binder;
    mpfr_set_ui(rop, 1, MPFR_RNDN);
    mpfr_inits(mode_binder, tmp, NULL); 
    mpfr_set_ui(mode_binder, 2*mode, MPFR_RNDN);
    mpfr_neg(mode_binder, mode_binder, MPFR_RNDN);
    mpfr_add_ui(mode_binder, mode_binder, 1, MPFR_RNDN);

    for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_pow(tmp, vector[mpfr_get_ui(i, MPFR_RNDN)], mode_binder, MPFR_RNDN);
        mpfr_mul(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_clears(i, tmp, NULL); 
    return;
}

__MATHSUITE inline ityp  mbase_mean(const register dim_typ dim, ityp vector[static dim])
{
    return (summation(dim, SUMMATION_SUM, vector)/dim);
}

__MATHSUITE inline void  math_mean(mpfr_t rop, mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_summation(rop, dim, SUMMATION_SUM, vector);
    mpfr_div(rop, rop, dim, MPFR_RNDN);
    return;
}

__MATHSUITE void  math_mode(mpfr_t rop, mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
	mpfr_t vmoda[mpfr_get_ui(dim, MPFR_RNDN)];
	mpfr_t moda;
	mpfr_t i, j;
	mpfr_init(moda);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
		for(mpfr_init_set_ui(vmoda[mpfr_get_ui(i, MPFR_RNDN)], 0, MPFR_RNDN), mpfr_init_set_ui(j, 0, MPFR_RNDN); mpfr_cmp(j, i) && mpfr_cmp(j, dim) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
			if(mpfr_cmp(vector[mpfr_get_ui(j, MPFR_RNDN)], vector[mpfr_get_ui(i, MPFR_RNDN)]) == 0)
                mpfr_add_ui(vmoda[mpfr_get_ui(i, MPFR_RNDN)], vmoda[mpfr_get_ui(i, MPFR_RNDN)], 1, MPFR_RNDN);
	MPFR_MINMAX(NULL, j, dim, vmoda, moda, NULL); // adjust this
    for(mpfr_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
        mpfr_clear(vmoda[mpfr_get_ui(i, MPFR_RNDN)]);
    mpfr_set(rop, vector[mpfr_get_ui(moda, MPFR_RNDN)], MPFR_RNDN);   
	mpfr_clears(i, j, moda, NULL); 
	return;
}

__MATHSUITE void  math_variance(mpfr_t rop, mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t i, tmp;
    mpfr_set_ui(rop, 0, MPFR_RNDN);
    mpfr_t media;
    mpfr_init(media);
    math_mean(media, dim, vector);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_sub(tmp, vector[mpfr_get_ui(i, MPFR_RNDN)], media, MPFR_RNDN);
        mpfr_exp2(tmp, tmp, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_div(rop, rop, dim, MPFR_RNDN);
    mpfr_clears(media, i, tmp, NULL); 
	return;
}

__MATHSUITE void  math_variance2(mpfr_t rop, mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t i, tmp;
    mpfr_set_ui(rop, 0, MPFR_RNDN);
    mpfr_t media;
    mpfr_init(media);
    math_mean(media, dim, vector);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_sub(tmp, vector[mpfr_get_ui(i, MPFR_RNDN)], media, MPFR_RNDN);
        mpfr_exp2(tmp, tmp, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }


    mpfr_sub_ui(tmp, dim, 1, MPFR_RNDN);
    mpfr_div(rop, rop, tmp, MPFR_RNDN);
    mpfr_clears(media, i, tmp, NULL); 
	return;
}

__MATHSUITE void  math_covariance(mpfr_t rop, mpfr_t dim, mpfr_t vector1[static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t vector2[static mpfr_get_ui(dim, MPFR_RNDN)])
{
	mpfr_t i, tmp, tmp2, media[2];
	mpfr_inits(media[FIRST_NUMBER], media[SECOND_NUMBER], i, tmp, tmp2, NULL); 
	math_mean(media[FIRST_NUMBER], dim, vector1);
    math_mean(media[SECOND_NUMBER], dim, vector2);
    mpfr_set_ui(rop, 0, MPFR_RNDN);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_sub(tmp, vector1[mpfr_get_ui(i, MPFR_RNDN)], media[FIRST_NUMBER], MPFR_RNDN);
        mpfr_sub(tmp2, vector2[mpfr_get_ui(i, MPFR_RNDN)], media[SECOND_NUMBER], MPFR_RNDN);
        mpfr_mul(tmp, tmp, tmp2, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_div(rop, rop, dim, MPFR_RNDN);
    mpfr_clears(i, tmp, tmp2, media[FIRST_NUMBER], media[SECOND_NUMBER], NULL); 
	return;
}

__MATHSUITE void  math_covariance2(mpfr_t rop, mpfr_t dim, mpfr_t vector1[static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t vector2[static mpfr_get_ui(dim, MPFR_RNDN)])
{
	mpfr_t i, tmp, tmp2, media[2];
	mpfr_inits(media[FIRST_NUMBER], media[SECOND_NUMBER], i, tmp, tmp2, NULL); 
	math_mean(media[FIRST_NUMBER], dim, vector1);
    math_mean(media[SECOND_NUMBER], dim, vector2);
    mpfr_set_ui(rop, 0, MPFR_RNDN);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_sub(tmp, vector1[mpfr_get_ui(i, MPFR_RNDN)], media[FIRST_NUMBER], MPFR_RNDN);
        mpfr_sub(tmp2, vector2[mpfr_get_ui(i, MPFR_RNDN)], media[SECOND_NUMBER], MPFR_RNDN);
        mpfr_mul(tmp, tmp, tmp2, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_sub_ui(tmp, dim, 1, MPFR_RNDN);
    mpfr_div(rop, rop, tmp, MPFR_RNDN);
    mpfr_clears(i, tmp, tmp2, media[FIRST_NUMBER], media[SECOND_NUMBER], NULL); 
	return;
}

__MATHSUITE void  math_stdcod(mpfr_t rop, mpfr_t dim, mpfr_t vector1[static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t vector2[static mpfr_get_ui(dim, MPFR_RNDN)])
{
	mpfr_t i, tmp, tmp2, media[2];
	mpfr_inits(media[FIRST_NUMBER], media[SECOND_NUMBER], i, tmp, tmp2, NULL); 
	math_mean(media[FIRST_NUMBER], dim, vector1);
    math_mean(media[SECOND_NUMBER], dim, vector2);
    mpfr_set_ui(rop, 0, MPFR_RNDN);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_sub(tmp, vector1[mpfr_get_ui(i, MPFR_RNDN)], media[FIRST_NUMBER], MPFR_RNDN);
        mpfr_sub(tmp2, vector2[mpfr_get_ui(i, MPFR_RNDN)], media[SECOND_NUMBER], MPFR_RNDN);
        mpfr_mul(tmp, tmp, tmp2, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_div(rop, rop, dim, MPFR_RNDN);
    mpfr_sqrt(rop, rop, MPFR_RNDN);
    mpfr_clears(i, tmp, tmp2, NULL); 
	return;
}

__MATHSUITE void  math_stdcod2(mpfr_t rop, mpfr_t dim, mpfr_t vector1[static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t vector2[static mpfr_get_ui(dim, MPFR_RNDN)])
{
	mpfr_t i, tmp, tmp2, media[2];
	mpfr_inits(media[FIRST_NUMBER], media[SECOND_NUMBER], i, tmp, tmp2, NULL); 
	math_mean(media[FIRST_NUMBER], dim, vector1);
    math_mean(media[SECOND_NUMBER], dim, vector2);
    mpfr_set_ui(rop, 0, MPFR_RNDN);

    for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_sub(tmp, vector1[mpfr_get_ui(i, MPFR_RNDN)], media[FIRST_NUMBER], MPFR_RNDN);
        mpfr_sub(tmp2, vector2[mpfr_get_ui(i, MPFR_RNDN)], media[SECOND_NUMBER], MPFR_RNDN);
        mpfr_mul(tmp, tmp, tmp2, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_sub_ui(tmp, dim, 1, MPFR_RNDN);
    mpfr_div(rop, rop, tmp, MPFR_RNDN);
    mpfr_sqrt(rop, rop, MPFR_RNDN);
    mpfr_clears(i, tmp, tmp2, media[FIRST_NUMBER], media[SECOND_NUMBER], NULL); 
    return;
}

__MATHSUITE void  math_stddev(mpfr_t rop, mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t i, tmp;
    mpfr_set_ui(rop, 0, MPFR_RNDN);
    mpfr_t media;
    mpfr_init(media);
    math_mean(media, dim, vector);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_sub(tmp, vector[mpfr_get_ui(i, MPFR_RNDN)], media, MPFR_RNDN);
        mpfr_exp2(tmp, tmp, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_div(rop, rop, dim, MPFR_RNDN);
    mpfr_sqrt(rop, rop, MPFR_RNDN);
    mpfr_clears(i, tmp, media, NULL); 
	return;
}

__MATHSUITE void  math_stddev2(mpfr_t rop, mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t i, tmp;
    mpfr_set_ui(rop, 0, MPFR_RNDN);
    mpfr_t media;
    mpfr_init(media);
    math_mean(media, dim, vector);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_sub(tmp, vector[mpfr_get_ui(i, MPFR_RNDN)], media, MPFR_RNDN);
        mpfr_exp2(tmp, tmp, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_sub_ui(tmp, dim, 1, MPFR_RNDN);
    mpfr_div(rop, rop, tmp, MPFR_RNDN);
    mpfr_sqrt(rop, rop, MPFR_RNDN);
    mpfr_clears(i, tmp, media, NULL); 
	return;
}

__MATHSUITE void  math_pearson(mpfr_t rop, mpfr_t dim, mpfr_t vector1[static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t vector2[static mpfr_get_ui(dim, MPFR_RNDN)])
{
	mpfr_t i, tmp, tmp2, media[2];
	mpfr_inits(media[FIRST_NUMBER], media[SECOND_NUMBER], i, tmp, tmp2, NULL); 
	math_mean(media[FIRST_NUMBER], dim, vector1);
    math_mean(media[SECOND_NUMBER], dim, vector2);
    mpfr_set_ui(rop, 0, MPFR_RNDN);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_sub(tmp, vector1[mpfr_get_ui(i, MPFR_RNDN)], media[FIRST_NUMBER], MPFR_RNDN);
        mpfr_sub(tmp2, vector2[mpfr_get_ui(i, MPFR_RNDN)], media[SECOND_NUMBER], MPFR_RNDN);
        mpfr_mul(tmp, tmp, tmp2, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_div(rop, rop, dim, MPFR_RNDN);
    math_stddev(tmp, dim, vector1);
    math_stddev(i, dim, vector2);
    mpfr_mul(tmp, tmp, i, MPFR_RNDN);
    mpfr_mul(rop, rop, tmp, MPFR_RNDN);
    mpfr_clears(i, tmp, media[FIRST_NUMBER], media[SECOND_NUMBER], NULL); 
	return;
}

__MATHSUITE void  math_pearson2(mpfr_t rop, mpfr_t dim, mpfr_t vector1[static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t vector2[static mpfr_get_ui(dim, MPFR_RNDN)])
{
	mpfr_t i, tmp, tmp2, media[2];
	mpfr_inits(media[FIRST_NUMBER], media[SECOND_NUMBER], i, tmp, tmp2, NULL); 
	math_mean(media[FIRST_NUMBER], dim, vector1);
    math_mean(media[SECOND_NUMBER], dim, vector2);
    mpfr_set_ui(rop, 0, MPFR_RNDN);

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        mpfr_sub(tmp, vector1[mpfr_get_ui(i, MPFR_RNDN)], media[FIRST_NUMBER], MPFR_RNDN);
        mpfr_sub(tmp2, vector2[mpfr_get_ui(i, MPFR_RNDN)], media[SECOND_NUMBER], MPFR_RNDN);
        mpfr_mul(tmp, tmp, tmp2, MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
    }

    mpfr_sub_ui(tmp, dim, 1, MPFR_RNDN);
    mpfr_div(rop, rop, tmp, MPFR_RNDN);
    math_stddev(tmp, dim, vector1);
    math_stddev(i, dim, vector2);
    mpfr_mul(tmp, tmp, i, MPFR_RNDN);
    mpfr_mul(rop, rop, tmp, MPFR_RNDN);
    mpfr_clears(i, tmp, media[FIRST_NUMBER], media[SECOND_NUMBER], NULL); 
	return;
}

__MATHSUITE inline bool  math_outlier(mpfr_t dim, mpfr_t outlier_idx, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
	return math_outlier2(dim, outlier_idx, access(curLayout)->outlier_constant, vector);
}

__MATHSUITE bool  math_outlier2(mpfr_t dim, mpfr_t outlier_idx, ityp outlier_constant, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    const register dim_typ dimd = mpfr_get_ui(dim, MPFR_RNDN);
    mpfr_t Q1, Q3;

	if(dimd%2)
    	mpfr_init_set(Q1, vector[(uint64_t)((dimd-3)*FIRST_QUARTILE_CONSTANT)], MPFR_RNDN);
    else
    {
    	mpfr_init(Q1);
    	mpfr_add(Q1, vector[(uint64_t)((dimd-4)*FIRST_QUARTILE_CONSTANT)], vector[(uint64_t)(dimd*FIRST_QUARTILE_CONSTANT)], MPFR_RNDN);
    	mpfr_div_ui(Q1, Q1, 2, MPFR_RNDN);
    }

    if(dimd%2)
		mpfr_init_set(Q3, vector[(uint64_t)(((dimd+1)*THIRD_QUARTILE_CONSTANT)-1)], MPFR_RNDN);
	else
	{
	    mpfr_init(Q3);
		mpfr_add(Q3, vector[(uint64_t)((dimd*THIRD_QUARTILE_CONSTANT)-1)], vector[(uint64_t)(dimd*THIRD_QUARTILE_CONSTANT)], MPFR_RNDN);
   		mpfr_div_ui(Q3, Q3, 2, MPFR_RNDN);
	}

	qsort(vector, dimd, sizeof(mpfr_t), mpfr_cmpfunc);
	mpfr_t deviation, tmp, tmp2;
	mpfr_inits(deviation, tmp, tmp2, NULL); 
	mpfr_sub(deviation, Q3, Q1, MPFR_RNDN);
	mpfr_mul_d(deviation, deviation, outlier_constant, MPFR_RNDN);
	mpfr_sub(tmp, Q1, deviation, MPFR_RNDN);
	mpfr_add(tmp2, Q3, deviation, MPFR_RNDN);
	const register dim_typ outidx = mpfr_get_ui(outlier_idx, MPFR_RNDN);
	mpfr_clears(Q1, Q3, tmp, tmp2, deviation, NULL); 
	return(mpfr_cmp(vector[outidx], tmp) < 0 || mpfr_cmp(vector[outidx], tmp2) > 0);
}

__MATHSUITE void  math_geomean(mpfr_t rop, mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t tmp;
    mpfr_init(tmp);
    mpfr_pow_si(rop, dim, -1, MPFR_RNDN);
    mpfr_product(tmp, dim, PRODUCTORY_MUL, vector);
    mpfr_pow(rop, tmp , rop, MPFR_RNDN);
    mpfr_clear(tmp);
    return;
}

__MATHSUITE void  math_armean(mpfr_t rop, mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t i;

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
        mpfr_pow_si(vector[mpfr_get_ui(i, MPFR_RNDN)], vector[mpfr_get_ui(i, MPFR_RNDN)], -1, MPFR_RNDN);

    mpfr_summation(rop, dim, SUMMATION_SUM, vector);
    mpfr_div(rop, dim, rop, MPFR_RNDN);
    mpfr_clear(i);
    return;

}

__MATHSUITE void  math_powmean(mpfr_t rop, mpfr_t dim, mpfr_t power, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t i;

	for(mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp(i, dim) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
        mpow2(vector[mpfr_get_ui(i, MPFR_RNDN)], vector[mpfr_get_ui(i, MPFR_RNDN)], power);

    mpfr_t tmp, tmp2;
    mpfr_inits(tmp, tmp2, NULL); 
    mpfr_pow_si(tmp, power, -1, MPFR_RNDN);
    mpfr_pow_si(rop, dim, -1, MPFR_RNDN);
    mpfr_summation(tmp2, dim, SUMMATION_SUM, vector);
    mpfr_mul(rop, rop, tmp2, MPFR_RNDN);
    mpfr_pow(rop, rop, tmp, MPFR_RNDN);
    mpfr_clears(tmp, i, NULL); 
    return;
}

__MATHSUITE inline void  math_scale(mpfr_t rop, mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)])
{
    mpfr_t tmp;
    mpfr_init(tmp);
    _MPFR_MIN(tmp, dim, vector);
    _MPFR_MAX(rop, dim, vector);
    mpfr_add(rop, rop, tmp, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    mpfr_clear(tmp);
    return;
}

__MATHSUITE void  math_first_quartile(mpfr_t rop, const register dim_typ dim, mpfr_t vector[static dim])
{
    
    if((dim-1)%4 == 0)
    {
    	mpfr_t rop2;
    	mpfr_init(rop2);
    	const register dim_typ idx = (dim-1)>>2;
    	randomizedSelect(rop, vector, 0, dim-1, idx);
    	mpfr_mul_d(rop, rop, FIRST_QUARTILE_CONSTANT, MPFR_RNDN);
    	randomizedSelect(rop2, vector, 0, dim-1, idx+1);
    	mpfr_mul_d(rop2, rop2, THIRD_QUARTILE_CONSTANT, MPFR_RNDN);
    	mpfr_add(rop, rop, rop2, MPFR_RNDN);
    	mpfr_clear(rop2);
    }
    else if((dim-3)%4 == 0)
    {
    	mpfr_t rop2;
    	mpfr_init(rop2);
    	const register dim_typ idx = (dim-3)>>2;
    	randomizedSelect(rop, vector, 0, dim-1, idx+1);
    	mpfr_mul_d(rop, rop, THIRD_QUARTILE_CONSTANT, MPFR_RNDN);
    	randomizedSelect(rop2, vector, 0, dim-1, idx+2);
    	mpfr_mul_d(rop2, rop2, FIRST_QUARTILE_CONSTANT, MPFR_RNDN);
    	mpfr_add(rop, rop, rop2, MPFR_RNDN);
    	mpfr_clear(rop2);
    }
    else // CASO PARI --> split
    {
    	if(dim%4)
			randomizedSelect(rop, vector, 0, dim-1, (dim+2)>>2);
		else
		{
			mpfr_t rop2;
			mpfr_init(rop2);
			randomizedSelect(rop, vector, 0, dim-1, dim>>2);
			randomizedSelect(rop2, vector, 0, dim-1, (dim>>2)+1);
			mpfr_add(rop, rop, rop2, MPFR_RNDN);
			mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
			mpfr_clear(rop2);
		}
	}
    return;
}

__MATHSUITE inline void  math_median(mpfr_t rop, const register dim_typ dim, mpfr_t vector[static dim])
{
    
	if(dim % 2)
		randomizedSelect(rop, vector, 0, dim-1, (dim+1)>>1);
	else
	{
		mpfr_t rop2;
		mpfr_init(rop2);
		randomizedSelect(rop, vector, 0, dim-1, dim>>1);
		randomizedSelect(rop2, vector, 0, dim-1, (dim>>1)+1);
		mpfr_add(rop, rop, rop2, MPFR_RNDN);
		mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
		mpfr_clear(rop2);
	}
	
    return;
}

__MATHSUITE void  math_third_quartile(mpfr_t rop, const register dim_typ dim, mpfr_t vector[static dim])
{
    
    if((dim-1)%4 == 0)
    {
    	mpfr_t rop2;
    	mpfr_init(rop2);
    	const register dim_typ idx = (dim-1)>>2;
    	randomizedSelect(rop, vector, 0, dim-1, 3*idx +1);
    	mpfr_mul_d(rop, rop, THIRD_QUARTILE_CONSTANT, MPFR_RNDN);
    	randomizedSelect(rop2, vector, 0, dim-1, 3*idx +2);
    	mpfr_mul_d(rop2, rop2, FIRST_QUARTILE_CONSTANT, MPFR_RNDN);
    	mpfr_add(rop, rop, rop2, MPFR_RNDN);
    	mpfr_clear(rop2);
    }
    else if((dim-3)%4 == 0)
    {
    	mpfr_t rop2;
    	mpfr_init(rop2);
    	const register dim_typ idx = (dim-3)>>2;
    	randomizedSelect(rop, vector, 0, dim-1, 3*idx + 2);
    	mpfr_mul_d(rop, rop, FIRST_QUARTILE_CONSTANT, MPFR_RNDN);
    	randomizedSelect(rop2, vector, 0, dim-1, 3*idx + 3);
    	mpfr_mul_d(rop2, rop2, THIRD_QUARTILE_CONSTANT, MPFR_RNDN);
    	mpfr_add(rop, rop, rop2, MPFR_RNDN);
    	mpfr_clear(rop2);
    }
    else // CASO PARI --> split
    {
    	if(dim%4)
			randomizedSelect(rop, vector, 0, dim-1, ((3*dim)+2)>>2);
		else
		{
			mpfr_t rop2;
			mpfr_init(rop2);
			randomizedSelect(rop, vector, 0, dim-1, 3*(dim>>2));
			randomizedSelect(rop2, vector, 0, dim-1, 3*(dim>>2)+1);
			mpfr_add(rop, rop, rop2, MPFR_RNDN);
			mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
			mpfr_clear(rop2);
		}
	}
    return;
}

__MATHSUITE dim_typ randomizedSelect(mpfr_t rop, mpfr_t A[], dim_typ p, dim_typ r, dim_typ i)
{
	if (p == r)
	{
		mpfr_set(rop, A[p], MPFR_RNDN);
		return p; 
	}
	    
	const register dim_typ q = randomizedPartition(A, p, r);
	const register dim_typ k = q - p + 1;
	
	if (i == k)
		mpfr_set(rop, A[q], MPFR_RNDN);
	else if (i < k)
	    randomizedSelect(rop, A, p, q-1, i) ;
	else
		randomizedSelect(rop, A, q+1, r, i - k);
	
	return USHRT_MAX;
}

__MATHSUITE dim_typ partition(mpfr_t A[], dim_typ p, dim_typ r)
{
	mpfr_t temp;
	mpfr_init(temp);
	register int i=p, j;

	for(j=p;j<r;++j)
	    if(mpfr_lessequal_p(A[j], A[r]))
	    {
	        mpfr_set(temp, A[i], MPFR_RNDN);
	        mpfr_set(A[i++], A[j], MPFR_RNDN); 
	        mpfr_set(A[j], temp, MPFR_RNDN); 
	    }
	
	mpfr_set(temp, A[i], MPFR_RNDN);
	mpfr_set(A[i], A[r], MPFR_RNDN);
	mpfr_set(A[r], temp, MPFR_RNDN);
	mpfr_clear(temp);
	return i;
}

__MATHSUITE inline dim_typ randomizedPartition(mpfr_t A[], dim_typ p, dim_typ r)
{
	mpfr_t temp;
	register dim_typ j = p + rand()%(r-p+1);
	
	mpfr_init_set(temp, A[r], MPFR_RNDN);
	mpfr_set(A[r], A[j], MPFR_RNDN); 
	mpfr_set(A[j], temp, MPFR_RNDN);
	mpfr_clear(temp);
	return partition(A, p, r);
}

 inline int64_t   powi(register int64_t x, register int64_t y) // calculates x^y
{
    int64_t base, res;
    for(base=x, res=1; y; y >>=1, base *= base)
        if(y&1)
            res *= base;
    return res;
}

 inline int64_t   changeBase(register int n, sel_typ start_base, sel_typ base)
{

    uint64_t h, sum;

    for(h=sum=0; n>0; ++h)
    {
        sum += (n%base)*powi(start_base, h);
        n /= base;
    }
    return sum;
}

 void  math_GCD(mpfr_t rop, mpfr_t a, mpfr_t b) // Euclide's Algorythm
{
    mpfr_t tmp, tmp2;
    mpfr_init_set(tmp, a, MPFR_RNDN);
    mpfr_init_set(rop, b, MPFR_RNDN);
    mpfr_init(tmp2);
    while(mpfr_cmp(a,b))
        if(mpfr_cmp(a,b) > 0)
        {
            mpfr_neg(tmp2, rop, MPFR_RNDN);
            mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
        }
        else
        {
            mpfr_neg(tmp2, tmp, MPFR_RNDN);
            mpfr_add(rop, rop, tmp2, MPFR_RNDN);
        }

    mpfr_clears(tmp, tmp2, NULL); 
    return;
}

// thanks elite.polito.it/files/courses/12BHD/progr/Esercizi-C-v2_01.pdf
// I renamed, modified and adapted to this program this function.
 void  math_lcm(mpfr_t rop, mpfr_t a, mpfr_t b)
{
    mpfr_t massimo, minimo;
    mpfr_t conta, tmp;
    bool fine;
    mpfr_init(tmp);
    mpfr_init_set(massimo, mpfr_cmp(a,b) > 0 ? a : b, MPFR_RNDN);
    mpfr_init_set(minimo, mpfr_cmp(a,b) < 0 ? a : b, MPFR_RNDN);
    mpfr_init_set_ui(conta, 1, MPFR_RNDN);
    mpfr_init_set_ui(rop, 0, MPFR_RNDN);
    fine = false;
    while(!(fine))
    {
        mpfr_mul(rop, conta, massimo, MPFR_RNDN);
        mpfr_remainder(tmp, rop, minimo, MPFR_RNDN);
        if(mpfr_zero_p(tmp)) fine = true;
        else mpfr_add_ui(conta, conta, 1, MPFR_RNDN);
    }
    mpfr_clears(massimo, minimo, conta, tmp, NULL); 
    return;
}

__MATHSUITE inline ityp  exp10(register ityp n)
{
    return pow(10.00, n);
}

__MATHSUITE inline ityp  expc(register ityp n)
{
    return (exp(n)/n);
}

__MATHSUITE inline ityp  exp10c(register ityp n)
{
    return (pow(10.00, n)/n);
}

__MATHSUITE inline ityp  exp2c(register ityp n)
{
    return (exp2(n)/n);
}

 inline ityp  logbN(register ityp b, register ityp N)
{
	return(log(N)/log(b));
}

__MATHSUITE inline ityp  logc(register ityp n)
{
    return (log(n)/n);
}

__MATHSUITE inline ityp  log10c(register ityp n)
{
    return (log10(n)/n);
}

__MATHSUITE inline ityp  log2c(register ityp n)
{
    return (log2(n)/n);
}

__MATHSUITE inline ityp  log1pc(register ityp n)
{
    return (log1p(n)/n);
}

__MATHSUITE inline ityp  log101p(register ityp n)
{
	return (logbN(10.00, 1.00+n));
}

__MATHSUITE inline ityp  log101pc(register ityp n)
{
	return (logbN(10.00, 1.00+n)/n);
}

__MATHSUITE inline ityp  log21p(register ityp n)
{
	return (logbN(2.00, 1.00+n));
}

__MATHSUITE inline ityp  log21pc(register ityp n)
{
	return (logbN(2.00, 1.00+n)/n);
}


// Properly Trigonometric Functions DEFINITIONS
// (by inlining them)

// It converts Radiants r to Degrees
__MATHSUITE inline ityp  deg(register ityp r)
{
    return ((180.0 * r) / M_PI);
}

// It converts Degrees d to Radiants
__MATHSUITE inline ityp  rad(register ityp d)
{
    return ((M_PI * d) / 180.0);
}

__MATHSUITE inline ityp  csc(register ityp n)
{
    return (1/sin(n));
}

__MATHSUITE inline ityp  sec(register ityp n)
{
    return (1/cos(n));
}

__MATHSUITE inline ityp  cot(register ityp n)
{
    return (1/tan(n));
}

__MATHSUITE inline ityp  csch(register ityp n)
{
    return (1/sinh(n));
}

__MATHSUITE inline ityp  sech(register ityp n)
{
    return (1/cosh(n));
}

__MATHSUITE inline ityp  coth(register ityp n)
{
    return (1/tanh(n));
}

__MATHSUITE inline ityp  acsc(register ityp n)
{
    return asin(1/n);
}

__MATHSUITE inline ityp  asec(register ityp n)
{
    return acos(1/n);
}

__MATHSUITE inline ityp  acot(register ityp n)
{
    return atan(1/n);
}

__MATHSUITE inline ityp  acsch(register ityp n)
{
    return asinh(1/n);
}

__MATHSUITE inline ityp  asech(register ityp n)
{
    return acosh(1/n);
}

__MATHSUITE inline ityp  acoth(register ityp n)
{
    return atanh(1/n);
}

__MATHSUITE inline ityp  hsin(register ityp n)
{
    return sin(n)/2.00;
}

__MATHSUITE inline ityp  hsinh(register ityp n)
{
    return sinh(n)/2.00;
}

__MATHSUITE inline ityp  qsin(register ityp n)
{
    return sin(n)/4.00;
}

__MATHSUITE inline ityp  qsinh(register ityp n)
{
    return sinh(n)/4.00;
}

__MATHSUITE inline ityp  hcos(register ityp n)
{
    return cos(n)/2.00;
}

__MATHSUITE inline ityp  hcosh(register ityp n)
{
    return cosh(n)/2.00;
}

__MATHSUITE inline ityp  qcos(register ityp n)
{
    return cos(n)/4.00;
}

__MATHSUITE inline ityp  qcosh(register ityp n)
{
    return cosh(n)/4.00;
}

__MATHSUITE inline ityp  hcsc(register ityp n)
{
    return csc(n)/2.00;
}

__MATHSUITE inline ityp  hcsch(register ityp n)
{
    return csch(n)/2.00;
}

__MATHSUITE inline ityp  qcsc(register ityp n)
{
    return csc(n)/4.00;
}

__MATHSUITE inline ityp  qcsch(register ityp n)
{
    return csch(n)/4.00;
}

__MATHSUITE inline ityp  hsec(register ityp n)
{
    return sec(n)/2.00;
}

__MATHSUITE inline ityp  hsech(register ityp n)
{
    return sech(n)/2.00;
}

__MATHSUITE inline ityp  qsec(register ityp n)
{
    return sec(n)/4.00;
}

__MATHSUITE inline ityp  qsech(register ityp n)
{
    return sech(n)/4.00;
}

__MATHSUITE inline ityp  htan(register ityp n)
{
    return tan(n)/2.00;
}

__MATHSUITE inline ityp  htanh(register ityp n)
{
    return tanh(n)/2.00;
}

__MATHSUITE inline ityp  qtan(register ityp n)
{
    return tan(n)/4.00;
}

__MATHSUITE inline ityp  qtanh(register ityp n)
{
    return tanh(n)/4.00;
}

__MATHSUITE inline ityp  hcot(register ityp n)
{
    return cot(n)/2.00;
}

__MATHSUITE inline ityp  hcoth(register ityp n)
{
    return coth(n)/2.00;
}

__MATHSUITE inline ityp  qcot(register ityp n)
{
    return cot(n)/4.00;
}

__MATHSUITE inline ityp  qcoth(register ityp n)
{
    return coth(n)/4.00;
}

__MATHSUITE inline ityp  vsin(register ityp n)
{
    return 1.00-cos(n);
}

__MATHSUITE inline ityp  cvsin(register ityp n)
{
    return 1.00-sin(n);
}

__MATHSUITE inline ityp  vcos(register ityp n)
{
    return 1.00+cos(n);
}

__MATHSUITE inline ityp  cvcos(register ityp n)
{
    return 1.00+sin(n);
}

__MATHSUITE inline ityp  hvsin(register ityp n)
{
    return (1.00-cos(n))/2.00;
}

__MATHSUITE inline ityp  hcvsin(register ityp n)
{
    return (1.00-sin(n))/2.00;
}

__MATHSUITE inline ityp  hvcos(register ityp n)
{
    return (1.00+cos(n))/2.00;
}

__MATHSUITE inline ityp  hcvcos(register ityp n)
{
    return (1.00+sin(n))/2.00;
}

__MATHSUITE inline ityp  qvsin(register ityp n)
{
    return (1.00-cos(n))/4.00;
}

__MATHSUITE inline ityp  qcvsin(register ityp n)
{
    return (1.00-sin(n))/4.00;
}

__MATHSUITE inline ityp  qvcos(register ityp n)
{
    return (1.00+cos(n))/4.00;
}

__MATHSUITE inline ityp  qcvcos(register ityp n)
{
    return (1.00+sin(n))/4.00;
}

__MATHSUITE inline ityp  vsinh(register ityp n)
{
    return 1.00-cosh(n);
}

__MATHSUITE inline ityp  cvsinh(register ityp n)
{
    return 1.00-sinh(n);
}

__MATHSUITE inline ityp  vcosh(register ityp n)
{
    return 1.00+cosh(n);
}

__MATHSUITE inline ityp  cvcosh(register ityp n)
{
    return 1.00+sinh(n);
}

__MATHSUITE inline ityp  hvsinh(register ityp n)
{
    return (1.00-cosh(n))/2.00;
}

__MATHSUITE inline ityp  hcvsinh(register ityp n)
{
    return (1.00-sinh(n))/2.00;
}

__MATHSUITE inline ityp  hvcosh(register ityp n)
{
    return (1.00+cosh(n))/2.00;
}

__MATHSUITE inline ityp  hcvcosh(register ityp n)
{
    return (1.00+sinh(n))/2.00;
}

__MATHSUITE inline ityp  qvsinh(register ityp n)
{
    return (1.00-cosh(n))/4.00;
}

__MATHSUITE inline ityp  qcvsinh(register ityp n)
{
    return (1.00-sinh(n))/4.00;
}

__MATHSUITE inline ityp  qvcosh(register ityp n)
{
    return (1.00+cosh(n))/4.00;
}

__MATHSUITE inline ityp  qcvcosh(register ityp n)
{
    return (1.00+sinh(n))/4.00;
}

__MATHSUITE inline ityp  esec(register ityp n)
{
    return sec(n)-1.00;
}

__MATHSUITE inline ityp  ecsc(register ityp n)
{
    return csc(n)-1.00;
}

__MATHSUITE inline ityp  esech(register ityp n)
{
    return sech(n)-1.00;
}

__MATHSUITE inline ityp  ecsch(register ityp n)
{
    return csch(n)-1.00;
}

__MATHSUITE inline ityp  hesec(register ityp n)
{
    return (sec(n)-1.00)/2.00;
}

__MATHSUITE inline ityp  hecsc(register ityp n)
{
    return (csc(n)-1.00)/2.00;
}

__MATHSUITE inline ityp  hesech(register ityp n)
{
    return (sech(n)-1.00)/2.00;
}

__MATHSUITE inline ityp  hecsch(register ityp n)
{
    return (csch(n)-1.00)/2.00;
}

__MATHSUITE inline ityp  qesec(register ityp n)
{
    return (sec(n)-1.00)/4.00;
}

__MATHSUITE inline ityp  qecsc(register ityp n)
{
    return (csc(n)-1.00)/4.00;
}

__MATHSUITE inline ityp  qesech(register ityp n)
{
    return (sech(n)-1.00)/4.00;
}

__MATHSUITE inline ityp  qecsch(register ityp n)
{
    return (csch(n)-1.00)/4.00;
}

__MATHSUITE inline ityp  sinc(register ityp n)
{
    return n ? sin(n)/n : 1.00;
}

__MATHSUITE inline ityp  sinch(register ityp n)
{
    return n ? sinh(n)/n : 1.00;
}

__MATHSUITE inline ityp  hsinc(register ityp n)
{
    return n ? sin(n)/(2.00*n) : 0.50;
}

__MATHSUITE inline ityp  hsinch(register ityp n)
{
    return n ? sinh(n)/(2.00*n) : 0.50;
}

__MATHSUITE inline ityp  qsinc(register ityp n)
{
    return n ? sin(n)/(4.00*n) : 0.25;
}

__MATHSUITE inline ityp  qsinch(register ityp n)
{
    return n ? sinh(n)/(4.00*n) : 0.25;
}

__MATHSUITE inline ityp  cosc(register ityp n)
{
    return cos(n)/n;
}

__MATHSUITE inline ityp  cosch(register ityp n)
{
    return cosh(n)/n;
}

__MATHSUITE inline ityp  hcosc(register ityp n)
{
    return cos(n)/(2.00*n);
}

__MATHSUITE inline ityp  hcosch(register ityp n)
{
    return cosh(n)/(2.00*n);
}

__MATHSUITE inline ityp  qcosc(register ityp n)
{
    return cos(n)/(4.00*n);
}

__MATHSUITE inline ityp  qcosch(register ityp n)
{
    return cosh(n)/(4.00*n);
}

__MATHSUITE inline ityp  secc(register ityp n)
{
    return sec(n)/n;
}

__MATHSUITE inline ityp  secch(register ityp n)
{
    return sech(n)/n;
}

__MATHSUITE inline ityp  hsecc(register ityp n)
{
    return sec(n)/(2.00*n);
}

__MATHSUITE inline ityp  hsecch(register ityp n)
{
    return sech(n)/(2.00*n);
}

__MATHSUITE inline ityp  qsecc(register ityp n)
{
    return sec(n)/(4.00*n);
}

__MATHSUITE inline ityp  qsecch(register ityp n)
{
    return sech(n)/(4*n);
}

__MATHSUITE inline ityp  cscc(register ityp n)
{
    return csc(n)/n;
}

__MATHSUITE inline ityp  cscch(register ityp n)
{
    return csch(n)/n;
}

__MATHSUITE inline ityp  hcscc(register ityp n)
{
    return csc(n)/(2.00*n);
}

__MATHSUITE inline ityp  hcscch(register ityp n)
{
    return csch(n)/(2.00*n);
}

__MATHSUITE inline ityp  qcscc(register ityp n)
{
    return csc(n)/(4.00*n);
}

__MATHSUITE inline ityp  qcscch(register ityp n)
{
    return csch(n)/(4.00*n);
}

__MATHSUITE inline ityp  tanc(register ityp n)
{
    return tan(n)/n;
}

__MATHSUITE inline ityp  tanch(register ityp n)
{
    return tanh(n)/n;
}

__MATHSUITE inline ityp  htanc(register ityp n)
{
    return tan(n)/(2.00*n);
}

__MATHSUITE inline ityp  htanch(register ityp n)
{
    return tanh(n)/(2.00*n);
}

__MATHSUITE inline ityp  qtanc(register ityp n)
{
    return tan(n)/(4.00*n);
}

__MATHSUITE inline ityp  qtanch(register ityp n)
{
    return tanh(n)/(4.00*n);
}

__MATHSUITE inline ityp  cotc(register ityp n)
{
    return cot(n)/n;
}

__MATHSUITE inline ityp  cotch(register ityp n)
{
    return coth(n)/n;
}

__MATHSUITE inline ityp  hcotc(register ityp n)
{
    return cot(n)/(2.00*n);
}

__MATHSUITE inline ityp  hcotch(register ityp n)
{
    return coth(n)/(2.00*n);
}

__MATHSUITE inline ityp  qcotc(register ityp n)
{
    return cot(n)/(4.00*n);
}

__MATHSUITE inline ityp  qcotch(register ityp n)
{
    return coth(n)/(4.00*n);
}

__MATHSUITE inline void  mpfr_expc(mpfr_t rop, mpfr_t n)
{
    mpfr_exp(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_exp10c(mpfr_t rop, mpfr_t n)
{
    mpfr_exp10(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_exp2c(mpfr_t rop, mpfr_t n)
{
    mpfr_exp2(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

 inline void  mpfr_logbN(mpfr_t rop, mpfr_t b, mpfr_t N)
{
    mpfr_t tmp;
    mpfr_init(tmp);
    mpfr_log(rop, N, MPFR_RNDN);
    mpfr_log(tmp, b, MPFR_RNDN);
    mpfr_div(rop, rop, tmp, MPFR_RNDN);
    mpfr_clear(tmp);
    return;
}

__MATHSUITE inline void  mpfr_logc(mpfr_t rop, mpfr_t n)
{
    mpfr_log(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_log10c(mpfr_t rop, mpfr_t n)
{
    mpfr_log10(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_log2c(mpfr_t rop, mpfr_t n)
{
    mpfr_log2(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mswrap_log1p(mpfr_t rop, mpfr_t n)
{
    mpfr_log1p(rop, n, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_floor(mpfr_t rop, mpfr_t n)
{
    mpfr_floor(rop, n);
	return;
}

__MATHSUITE inline void  mswrap_ceil(mpfr_t rop, mpfr_t n)
{
    mpfr_ceil(rop, n);
	return;
}

__MATHSUITE inline void  mpfr_log1pc(mpfr_t rop, mpfr_t n)
{
    mpfr_log1p(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_log101p(mpfr_t rop, mpfr_t n)
{
    mpfr_t tmp;
    mpfr_set_ui(rop, 10, MPFR_RNDN);
    mpfr_init_set(tmp, n, MPFR_RNDN);
    mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
    mpfr_logbN(rop, rop, tmp);
    mpfr_clear(tmp);
    return;
}

__MATHSUITE inline void  mpfr_log101pc(mpfr_t rop, mpfr_t n)
{
    mpfr_t tmp;
    mpfr_set_ui(rop, 10, MPFR_RNDN);
    mpfr_init_set(tmp, n, MPFR_RNDN);
    mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
    mpfr_logbN(rop, rop, tmp);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_clears(tmp, NULL); 
    return;
}

__MATHSUITE inline void  mpfr_log21p(mpfr_t rop, mpfr_t n)
{
    mpfr_t tmp;
    mpfr_set_ui(rop, 2, MPFR_RNDN);
    mpfr_init_set(tmp, n, MPFR_RNDN);
    mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
    mpfr_logbN(rop, rop, tmp);
    mpfr_clears(tmp, NULL); 
    return;
}


__MATHSUITE inline void  mpfr_log21pc(mpfr_t rop, mpfr_t n)
{
    mpfr_t tmp;
    mpfr_set_ui(rop, 2, MPFR_RNDN);
    mpfr_init_set(tmp, n, MPFR_RNDN);
    mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
    mpfr_logbN(rop, rop, tmp);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_clears(tmp, NULL); 
    return;
}

 inline void  rootnX(mpfr_t rop, mpfr_t n, mpfr_t X)
{
    mpfr_pow_si(rop, n, -1, MPFR_RNDN);
    mpfr_pow(rop, X, rop, MPFR_RNDN);
    return;
}

__MATHSUITE inline double complex  cexpc(register double complex n)
{
    return (cexp(n)/n);
}

__MATHSUITE inline double complex  cexp10(register double complex n)
{
    return cpow(10.00, n);
}

__MATHSUITE inline double complex  cexp10c(register double complex n)
{
    return (cpow(10.00, n)/n);
}

__MATHSUITE inline double complex  cexp2(register double complex n)
{
	return (cpow(2.00, n));
}

__MATHSUITE inline double complex  cexp2c(register double complex n)
{
    return (cpow(2.00, n)/n);
}

 inline double complex  clogbN(register double complex b, register double complex N)
{
	return(clog(N)/clog(b));
}

__MATHSUITE inline double complex  clogc(register double complex n)
{
    return (clog(n)/n);
}

__MATHSUITE inline double complex  clog10(register double complex n)
{
	return (clogbN(10.00, n));
}

__MATHSUITE inline double complex  clog10c(register double complex n)
{
    return (clogbN(10.00, n)/n);
}

__MATHSUITE inline double complex  clog2(register double complex n)
{
	return (clogbN(2.00, n));
}

__MATHSUITE inline double complex  clog2c(register double complex n)
{
    return (clogbN(2.00, n)/n);
}

__MATHSUITE inline double complex  clog1p(register double complex n)
{
	return (clog(1.00+n));
}

__MATHSUITE inline double complex  clog1pc(register double complex n)
{
    return (clog(1.00+n)/n);
}

__MATHSUITE inline double complex  clog101p(register double complex n)
{
	return (clogbN(10.00, 1.00+n));
}

__MATHSUITE inline double complex  clog101pc(register double complex n)
{
    return (clogbN(10.00, 1.00+n)/n);
}

__MATHSUITE inline double complex  clog21p(register double complex n)
{
	return (clogbN(2.00, 1.00+n));
}

__MATHSUITE inline double complex  clog21pc(register double complex n)
{
    return (clogbN(2.00, 1.00+n)/n);
}

 inline double complex  ccbrt(register double complex n)
{
    return(cpow(n, 1/3));
}

 inline double complex  crootnX(register double complex n, register double complex X)
{
    return(cpow(X, 1/n));
}

// TEMPERATURE GRADES CONVERSIONS Functions DEFINITIONS
// if mode is TRUE then the CONVERSION is DIRECT, elsewhere it is INDIRECT

// Celsius to Fahrehneit
__MATHSUITE void  cel_fah(mpfr_t rop, mpfr_t grad, const bool mode)
{
    if(mode)
    {
        mpfr_mul_d(rop, grad, 1.80, MPFR_RNDN);
        mpfr_add_ui(rop, rop, 32, MPFR_RNDN);
    }
    else
    {
        mpfr_sub_ui(rop, grad, 32, MPFR_RNDN);
        mpfr_mul_d(rop, rop, 0.5555555556, MPFR_RNDN);
    }
    return;
}

// Celsius to Kelvin
__MATHSUITE inline void  cel_kel(mpfr_t rop, mpfr_t grad, const bool mode)
{
    mpfr_add_d(rop, grad, 273.15*(1.00-(mode<<1)), MPFR_RNDN);
    return;
}

// Celsius to Rankine
__MATHSUITE void  cel_rank(mpfr_t rop, mpfr_t grad, const bool mode)
{
    if(mode)
    {
        mpfr_mul_d(rop, grad, 1.80, MPFR_RNDN);
        mpfr_add_ui(rop, rop, 32, MPFR_RNDN);
        mpfr_add_d(rop, rop, 459.67, MPFR_RNDN);
    }
    else
    {
        mpfr_sub_ui(rop, grad, 32, MPFR_RNDN);
        mpfr_sub_d(rop, rop, 459.67, MPFR_RNDN);
        mpfr_div_d(rop, rop, 1.80, MPFR_RNDN);
    }
    return;
}

// Celsius to Réaumur
__MATHSUITE inline void  cel_rea(mpfr_t rop, mpfr_t grad, const bool mode)
{
    mpfr_mul_d(rop, grad, mode ? 0.80 : 1.25, MPFR_RNDN);
    return;
}

// Celsius to Newton
__MATHSUITE inline void  cel_new(mpfr_t rop, mpfr_t grad, const bool mode)
{
    mpfr_mul_d(rop, grad, mode ? (33/100) : (100/33), MPFR_RNDN);
    return;
}

// Celsius to Delisle
__MATHSUITE inline void  cel_del(mpfr_t rop, mpfr_t grad, const bool mode)
{
    mpfr_t tmp;
    mpfr_set(rop, grad, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    if(mode)
    {
        mpfr_add_ui(rop, rop, 100, MPFR_RNDN);
        mpfr_mul_d(rop, rop, 1.50, MPFR_RNDN);
    }
    else
    {
        mpfr_mul_d(rop, rop, 1.50, MPFR_RNDN);
        mpfr_add_ui(rop, rop, 100, MPFR_RNDN);
    }
    return;
}

// Celsius to Rømer
__MATHSUITE inline void  cel_rom(mpfr_t rop, mpfr_t grad, const bool mode)
{
    if(mode)
    {
        mpfr_mul_d(rop, grad, 0.52500, MPFR_RNDN);
        mpfr_add_d(rop, rop, 7.50, MPFR_RNDN);
    }
    else
    {
        mpfr_sub_d(rop, grad, 7.50, MPFR_RNDN);
        mpfr_mul_d(rop, rop, 1.9047619047619047619047619047619, MPFR_RNDN);
    }
    return;
}

// Fahrenheit to Kelvin
__MATHSUITE inline void  fah_kel(mpfr_t rop, mpfr_t grad, const bool mode)
{
    if(mode)
    {
        mpfr_add_d(rop, grad, 459.67, MPFR_RNDN);
        mpfr_div_d(rop, rop, 1.80, MPFR_RNDN);
    }
    else
    {
        mpfr_mul_d(rop, grad, 1.80, MPFR_RNDN);
        mpfr_sub_d(rop, rop, 459.67, MPFR_RNDN);
    }
    return;
}

// Fahrenheit to Rankine
__MATHSUITE inline void  fah_rank(mpfr_t rop, mpfr_t grad, const bool mode)
{
    mpfr_add_d(rop, grad, (459.67*(-1+(mode<<1))), MPFR_RNDN);
    return;
}

// Fahrenheit to Réaumur
__MATHSUITE inline void  fah_rea(mpfr_t rop, mpfr_t grad, const bool mode)
{
    if(mode)
    {
        mpfr_sub_ui(rop, grad, 32, MPFR_RNDN);
        mpfr_div_d(rop, rop, 2.25, MPFR_RNDN);
    }
    else
    {
        mpfr_mul_d(rop, grad, 2.25, MPFR_RNDN);
        mpfr_add_ui(rop, rop, 32, MPFR_RNDN);
    }
    return;
}

// Réaumur to Rankine
__MATHSUITE inline void  rea_rank(mpfr_t rop, mpfr_t grad, const bool mode)
{
    if(mode)
    {
        mpfr_mul_d(rop, grad, 2.25, MPFR_RNDN);
        mpfr_add_d(rop, rop, 491.67, MPFR_RNDN);
    }
    else
    {
        mpfr_sub_d(rop, grad, 491.67, MPFR_RNDN);
        mpfr_div_d(rop, rop, 2.25, MPFR_RNDN);
    }
    return;
}
//


// It converts "INTERNATIONAL EUROPEAN" (m/s) speed measures system
// into (km/h) measures system
__MATHSUITE inline void  speed(mpfr_t rop, mpfr_t value, bool mode)
{
    mpfr_t tmp; 
    mpfr_set_d(rop, 3.60, MPFR_RNDN);
    mpfr_init_set_si(tmp, 1.00-(mode<<1), MPFR_RNDN);
    mpow2(rop, rop, tmp);
    mpfr_mul(rop, rop, value, MPFR_RNDN);
    mpfr_clear(tmp);
    return;
}


// Properly Trigonometric Functions DEFINITIONS
// (by inlining them)

__MATHSUITE inline void  mswrap_sin(mpfr_t rop, register mpfr_t r)
{
	mpfr_sin(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_cos(mpfr_t rop, register mpfr_t r)
{
	mpfr_cos(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_sinh(mpfr_t rop, register mpfr_t r)
{
	mpfr_sinh(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_cosh(mpfr_t rop, register mpfr_t r)
{
	mpfr_cosh(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_sec(mpfr_t rop, register mpfr_t r)
{
	mpfr_sec(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_sech(mpfr_t rop, register mpfr_t r)
{
	mpfr_sech(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_csc(mpfr_t rop, register mpfr_t r)
{
	mpfr_csc(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_csch(mpfr_t rop, register mpfr_t r)
{
	mpfr_csch(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_asin(mpfr_t rop, register mpfr_t r)
{
	mpfr_asin(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_asinh(mpfr_t rop, register mpfr_t r)
{
	mpfr_asinh(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_acos(mpfr_t rop, register mpfr_t r)
{
	mpfr_acos(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_acosh(mpfr_t rop, register mpfr_t r)
{
	mpfr_acosh(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_tan(mpfr_t rop, register mpfr_t r)
{
	mpfr_tan(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_tanh(mpfr_t rop, register mpfr_t r)
{
	mpfr_tanh(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_atan(mpfr_t rop, register mpfr_t r)
{
	mpfr_atan(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_atanh(mpfr_t rop, register mpfr_t r)
{
	mpfr_atanh(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_cot(mpfr_t rop, register mpfr_t r)
{
	mpfr_cot(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_coth(mpfr_t rop, register mpfr_t r)
{
	mpfr_coth(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_log(mpfr_t rop, register mpfr_t r)
{
	mpfr_log(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_log10(mpfr_t rop, register mpfr_t r)
{
	mpfr_log10(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_log2(mpfr_t rop, register mpfr_t r)
{
	mpfr_log2(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_exp(mpfr_t rop, register mpfr_t r)
{
	mpfr_exp(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_exp10(mpfr_t rop, register mpfr_t r)
{
	mpfr_exp10(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  mswrap_exp2(mpfr_t rop, register mpfr_t r)
{
	mpfr_exp2(rop, r, MPFR_RNDN);
	return;
}

__MATHSUITE inline ityp  cotan(const register ityp a)
{
	return 1/tan(a);
}

__MATHSUITE inline ityp  degd(const register ityp a)
{
	return a*(180/M_PI);
}

// It converts Radiants r to Degrees
__MATHSUITE inline void  mpfr_deg(mpfr_t rop, register mpfr_t r)
{
    mpfr_t tmp;
    mpfr_init(tmp);
    mpfr_set(rop, r, MPFR_RNDN);
    mpfr_mul_ui(rop, rop, 180, MPFR_RNDN);
    mpfr_const_pi(tmp, MPFR_RNDN);
    mpfr_div(rop, rop, tmp, MPFR_RNDN);
    mpfr_clear(tmp);
    return;
}

// It converts Degrees d to Radiants
__MATHSUITE inline void  mpfr_rad(mpfr_t rop, mpfr_t d)
{
    mpfr_t tmp;
    mpfr_init(tmp);
    mpfr_set(rop, d, MPFR_RNDN);
    mpfr_const_pi(tmp, MPFR_RNDN);
    mpfr_mul(rop, rop, tmp, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 180, MPFR_RNDN);
    mpfr_clear(tmp);
    return;
}

__MATHSUITE inline void  mpfr_acsc(mpfr_t rop, mpfr_t n)
{
    mpfr_pow_si(rop, n, -1, MPFR_RNDN);
    mpfr_asin(rop, rop, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_asec(mpfr_t rop, mpfr_t n)
{
    mpfr_pow_si(rop, n, -1, MPFR_RNDN);
    mpfr_acos(rop, rop, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_acot(mpfr_t rop, mpfr_t n)
{
    mpfr_pow_si(rop, n, -1, MPFR_RNDN);
    mpfr_atan(rop, rop, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_acsch(mpfr_t rop, mpfr_t n)
{
    mpfr_pow_si(rop, n, -1, MPFR_RNDN);
    mpfr_asinh(rop, rop, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_asech(mpfr_t rop, mpfr_t n)
{
    mpfr_pow_si(rop, n, -1, MPFR_RNDN);
    mpfr_acosh(rop, rop, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_acoth(mpfr_t rop, mpfr_t n)
{
    mpfr_pow_si(rop, n, -1, MPFR_RNDN);
    mpfr_atanh(rop, rop, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hsin(mpfr_t rop, mpfr_t n)
{
    mpfr_sin(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hsinh(mpfr_t rop, mpfr_t n)
{
    mpfr_sinh(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qsin(mpfr_t rop, mpfr_t n)
{
    mpfr_sin(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qsinh(mpfr_t rop, mpfr_t n)
{
    mpfr_sinh(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcos(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcosh(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcos(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcosh(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcsc(mpfr_t rop, mpfr_t n)
{
    mpfr_csc(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcsch(mpfr_t rop, mpfr_t n)
{
    mpfr_csch(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcsc(mpfr_t rop, mpfr_t n)
{
    mpfr_csc(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcsch(mpfr_t rop, mpfr_t n)
{
    mpfr_csch(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hsec(mpfr_t rop, mpfr_t n)
{
    mpfr_sec(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hsech(mpfr_t rop, mpfr_t n)
{
    mpfr_sech(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qsec(mpfr_t rop, mpfr_t n)
{
    mpfr_sec(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qsech(mpfr_t rop, mpfr_t n)
{
    mpfr_sech(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_htan(mpfr_t rop, mpfr_t n)
{

    mpfr_tan(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_htanh(mpfr_t rop, mpfr_t n)
{
    mpfr_tanh(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qtan(mpfr_t rop, mpfr_t n)
{
    mpfr_tan(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qtanh(mpfr_t rop, mpfr_t n)
{
    mpfr_tanh(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcot(mpfr_t rop, mpfr_t n)
{
    mpfr_cot(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcoth(mpfr_t rop, mpfr_t n)
{
    mpfr_coth(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcot(mpfr_t rop, mpfr_t n)
{
    mpfr_cot(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcoth(mpfr_t rop, mpfr_t n)
{
    mpfr_coth(rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_vsin(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_cvsin(mpfr_t rop, mpfr_t n)
{
    mpfr_sin(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_vcos(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_cvcos(mpfr_t rop, mpfr_t n)
{
    mpfr_sin(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hvsin(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcvsin(mpfr_t rop, mpfr_t n)
{
    mpfr_sin(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hvcos(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcvcos(mpfr_t rop, mpfr_t n)
{
    mpfr_sin(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qvsin(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcvsin(mpfr_t rop, mpfr_t n)
{
    mpfr_sin(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qvcos(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcvcos(mpfr_t rop, mpfr_t n)
{
    mpfr_sin(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_vsinh(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_cvsinh(mpfr_t rop, mpfr_t n)
{
    mpfr_sinh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_vcosh(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_cvcosh(mpfr_t rop, mpfr_t n)
{
    mpfr_sinh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hvsinh(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcvsinh(mpfr_t rop, mpfr_t n)
{
    mpfr_sinh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hvcosh(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcvcosh(mpfr_t rop, mpfr_t n)
{
    mpfr_sinh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qvsinh(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcvsinh(mpfr_t rop, mpfr_t n)
{
    mpfr_sinh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qvcosh(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcvcosh(mpfr_t rop, mpfr_t n)
{
    mpfr_sinh(rop, n, MPFR_RNDN);
    mpfr_neg(rop, rop, MPFR_RNDN);
    mpfr_add_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_esec(mpfr_t rop, mpfr_t n)
{
    mpfr_sec(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_ecsc(mpfr_t rop, mpfr_t n)
{
    mpfr_csc(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_esech(mpfr_t rop, mpfr_t n)
{
    mpfr_sech(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_ecsch(mpfr_t rop, mpfr_t n)
{
    mpfr_csch(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hesec(mpfr_t rop, mpfr_t n)
{
    mpfr_sec(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hecsc(mpfr_t rop, mpfr_t n)
{
    mpfr_csc(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hesech(mpfr_t rop, mpfr_t n)
{
    mpfr_sech(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hecsch(mpfr_t rop, mpfr_t n)
{
    mpfr_csch(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qesec(mpfr_t rop, mpfr_t n)
{
    mpfr_sec(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qecsc(mpfr_t rop, mpfr_t n)
{
    mpfr_csc(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qesech(mpfr_t rop, mpfr_t n)
{
    mpfr_sech(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qecsch(mpfr_t rop, mpfr_t n)
{
    mpfr_csch(rop, n, MPFR_RNDN);
    mpfr_sub_ui(rop, rop, 1, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_sinc(mpfr_t rop, mpfr_t n)
{
    if(n)
    {
        mpfr_sin(rop, n, MPFR_RNDN);
        mpfr_div(rop, rop, n, MPFR_RNDN);
    }
    else
        mpfr_set_ui(rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_sinch(mpfr_t rop, mpfr_t n)
{
    if(n)
    {
        mpfr_sinh(rop, n, MPFR_RNDN);
        mpfr_div(rop, rop, n, MPFR_RNDN);
    }
    else
        mpfr_set_ui(rop, 1, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hsinc(mpfr_t rop, mpfr_t n)
{
    if(n)
    {
        mpfr_sin(rop, n, MPFR_RNDN);
        mpfr_div(rop, rop, n, MPFR_RNDN);
        mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    }
    else
        mpfr_set_d(rop, 0.50, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hsinch(mpfr_t rop, mpfr_t n)
{
    if(n)
    {
        mpfr_sinh(rop, n, MPFR_RNDN);
        mpfr_div(rop, rop, n, MPFR_RNDN);
        mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    }
    else
        mpfr_set_d(rop, 0.50, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qsinc(mpfr_t rop, mpfr_t n)
{
    if(n)
    {
        mpfr_sin(rop, n, MPFR_RNDN);
        mpfr_div(rop, rop, n, MPFR_RNDN);
        mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    }
    else
        mpfr_set_d(rop, 0.25, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qsinch(mpfr_t rop, mpfr_t n)
{
    if(n)
    {
        mpfr_sinh(rop, n, MPFR_RNDN);
        mpfr_div(rop, rop, n, MPFR_RNDN);
        mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    }
    else
        mpfr_set_d(rop, 0.25, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_cosc(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_cosch(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcosc(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcosch(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcosc(mpfr_t rop, mpfr_t n)
{
    mpfr_cos(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcosch(mpfr_t rop, mpfr_t n)
{
    mpfr_cosh(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_secc(mpfr_t rop, mpfr_t n)
{
    mpfr_sec(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_secch(mpfr_t rop, mpfr_t n)
{
    mpfr_sech(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hsecc(mpfr_t rop, mpfr_t n)
{
    mpfr_sec(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hsecch(mpfr_t rop, mpfr_t n)
{
    mpfr_sech(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qsecc(mpfr_t rop, mpfr_t n)
{
    mpfr_sec(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qsecch(mpfr_t rop, mpfr_t n)
{
    mpfr_sech(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_cscc(mpfr_t rop, mpfr_t n)
{
    mpfr_csc(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_cscch(mpfr_t rop, mpfr_t n)
{
    mpfr_csch(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcscc(mpfr_t rop, mpfr_t n)
{
    mpfr_csc(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcscch(mpfr_t rop, mpfr_t n)
{
    mpfr_csch(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcscc(mpfr_t rop, mpfr_t n)
{
    mpfr_csc(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcscch(mpfr_t rop, mpfr_t n)
{
    mpfr_csch(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_tanc(mpfr_t rop, mpfr_t n)
{
    mpfr_tan(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_tanch(mpfr_t rop, mpfr_t n)
{
    mpfr_tanh(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_htanc(mpfr_t rop, mpfr_t n)
{
    mpfr_tan(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_htanch(mpfr_t rop, mpfr_t n)
{
    mpfr_tanh(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qtanc(mpfr_t rop, mpfr_t n)
{
    mpfr_tan(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qtanch(mpfr_t rop, mpfr_t n)
{
    mpfr_tanh(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_cotc(mpfr_t rop, mpfr_t n)
{
    mpfr_cot(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_cotch(mpfr_t rop, mpfr_t n)
{
    mpfr_coth(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    return;
}
__MATHSUITE inline void  mpfr_hcotc(mpfr_t rop, mpfr_t n)
{
    mpfr_cot(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_hcotch(mpfr_t rop, mpfr_t n)
{
    mpfr_coth(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 2, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcotc(mpfr_t rop, mpfr_t n)
{
    mpfr_cot(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline void  mpfr_qcotch(mpfr_t rop, mpfr_t n)
{
    mpfr_coth(rop, n, MPFR_RNDN);
    mpfr_div(rop, rop, n, MPFR_RNDN);
    mpfr_div_ui(rop, rop, 4, MPFR_RNDN);
    return;
}

__MATHSUITE inline double complex  ccsc(register double complex n)
{
    return (1/csin(n));
}

__MATHSUITE inline double complex  csec(register double complex n)
{
    return (1/ccos(n));
}

__MATHSUITE inline double complex  ccot(register double complex n)
{
    return (1/ctan(n));
}

__MATHSUITE inline double complex  ccsch(register double complex n)
{
    return (1/csinh(n));
}

__MATHSUITE inline double complex  csech(register double complex n)
{
    return (1/ccosh(n));
}

__MATHSUITE inline double complex  ccoth(register double complex n)
{
    return (1/ctanh(n));
}

__MATHSUITE inline double complex  cacsc(register double complex n)
{
    return casin(1/n);
}

__MATHSUITE inline double complex  casec(register double complex n)
{
    return cacos(1/n);
}

__MATHSUITE inline double complex  cacot(register double complex n)
{
    return catan(1/n);
}

__MATHSUITE inline double complex  cacsch(register double complex n)
{
    return casinh(1/n);
}

__MATHSUITE inline double complex  casech(register double complex n)
{
    return cacosh(1/n);
}

__MATHSUITE inline double complex  cacoth(register double complex n)
{
    return catanh(1/n);
}

__MATHSUITE inline double complex  chsin(register double complex n)
{
    return csin(n)/2.00;
}

__MATHSUITE inline double complex  chsinh(register double complex n)
{
    return csinh(n)/2.00;
}

__MATHSUITE inline double complex  cqsin(register double complex n)
{
    return csin(n)/4.00;
}

__MATHSUITE inline double complex  cqsinh(register double complex n)
{
    return csinh(n)/4.00;
}

__MATHSUITE inline double complex  chcos(register double complex n)
{
    return ccos(n)/2.00;
}

__MATHSUITE inline double complex  chcosh(register double complex n)
{
    return ccosh(n)/2.00;
}

__MATHSUITE inline double complex  cqcos(register double complex n)
{
    return ccos(n)/4.00;
}

__MATHSUITE inline double complex  cqcosh(register double complex n)
{
    return ccosh(n)/4.00;
}

__MATHSUITE inline double complex  chcsc(register double complex n)
{
    return ccsc(n)/2.00;
}

__MATHSUITE inline double complex  chcsch(register double complex n)
{
    return ccsch(n)/2.00;
}

__MATHSUITE inline double complex  cqcsc(register double complex n)
{
    return ccsc(n)/4.00;
}

__MATHSUITE inline double complex  cqcsch(register double complex n)
{
    return ccsch(n)/4.00;
}

__MATHSUITE inline double complex  chsec(register double complex n)
{
    return csec(n)/2.00;
}

__MATHSUITE inline double complex  chsech(register double complex n)
{
    return csech(n)/2.00;
}

__MATHSUITE inline double complex  cqsec(register double complex n)
{
    return csec(n)/4.00;
}

__MATHSUITE inline double complex  cqsech(register double complex n)
{
    return csech(n)/4.00;
}

__MATHSUITE inline double complex  chtan(register double complex n)
{
    return ctan(n)/2.00;
}

__MATHSUITE inline double complex  chtanh(register double complex n)
{
    return ctanh(n)/2.00;
}

__MATHSUITE inline double complex  cqtan(register double complex n)
{
    return ctan(n)/4.00;
}

__MATHSUITE inline double complex  cqtanh(register double complex n)
{
    return ctanh(n)/4.00;
}

__MATHSUITE inline double complex  chcot(register double complex n)
{
    return ccot(n)/2.00;
}

__MATHSUITE inline double complex  chcoth(register double complex n)
{
    return ccoth(n)/2.00;
}

__MATHSUITE inline double complex  cqcot(register double complex n)
{
    return ccot(n)/4.00;
}

__MATHSUITE inline double complex  cqcoth(register double complex n)
{
    return ccoth(n)/4.00;
}

__MATHSUITE inline double complex  cpxvsin(register double complex n)
{
    return 1.00-ccos(n);
}

__MATHSUITE inline double complex  ccvsin(register double complex n)
{
    return 1.00-csin(n);
}

__MATHSUITE inline double complex  cpxvcos(register double complex n)
{
    return 1.00+ccos(n);
}

__MATHSUITE inline double complex  ccvcos(register double complex n)
{
    return 1.00+csin(n);
}

__MATHSUITE inline double complex  chvsin(register double complex n)
{
    return (1.00-ccos(n))/2.00;
}

__MATHSUITE inline double complex  chcvsin(register double complex n)
{
    return (1.00-csin(n))/2.00;
}

__MATHSUITE inline double complex  chvcos(register double complex n)
{
    return (1.00+ccos(n))/2.00;
}

__MATHSUITE inline double complex  chcvcos(register double complex n)
{
    return (1.00+csin(n))/2.00;
}

__MATHSUITE inline double complex  cqvsin(register double complex n)
{
    return (1.00-ccos(n))/4.00;
}

__MATHSUITE inline double complex  cqcvsin(register double complex n)
{
    return (1.00-csin(n))/4.00;
}

__MATHSUITE inline double complex  cqvcos(register double complex n)
{
    return (1.00+ccos(n))/4.00;
}

__MATHSUITE inline double complex  cqcvcos(register double complex n)
{
    return (1.00+csin(n))/4.00;
}

__MATHSUITE inline double complex  cpxvsinh(register double complex n)
{
    return 1.00-ccosh(n);
}

__MATHSUITE inline double complex  ccvsinh(register double complex n)
{
    return 1.00-csinh(n);
}

__MATHSUITE inline double complex  cpxvcosh(register double complex n)
{
    return 1.00+ccosh(n);
}

__MATHSUITE inline double complex  ccvcosh(register double complex n)
{
    return 1.00+csinh(n);
}

__MATHSUITE inline double complex  chvsinh(register double complex n)
{
    return (1.00-ccosh(n))/2.00;
}

__MATHSUITE inline double complex  chcvsinh(register double complex n)
{
    return (1.00-csinh(n))/2.00;
}

__MATHSUITE inline double complex  chvcosh(register double complex n)
{
    return (1.00+ccosh(n))/2.00;
}

__MATHSUITE inline double complex  chcvcosh(register double complex n)
{
    return (1.00+csinh(n))/2.00;
}

__MATHSUITE inline double complex  cqvsinh(register double complex n)
{
    return (1.00-ccosh(n))/4.00;
}

__MATHSUITE inline double complex  cqcvsinh(register double complex n)
{
    return (1.00-csinh(n))/4.00;
}

__MATHSUITE inline double complex  cqvcosh(register double complex n)
{
    return (1.00+ccosh(n))/4.00;
}

__MATHSUITE inline double complex  cqcvcosh(register double complex n)
{
    return (1.00+csinh(n))/4.00;
}

__MATHSUITE inline double complex  cesec(register double complex n)
{
    return csec(n)-1.00;
}

__MATHSUITE inline double complex  cecsc(register double complex n)
{
    return ccsc(n)-1.00;
}

__MATHSUITE inline double complex  cesech(register double complex n)
{
    return csech(n)-1.00;
}

__MATHSUITE inline double complex  cecsch(register double complex n)
{
    return ccsch(n)-1.00;
}

__MATHSUITE inline double complex  chesec(register double complex n)
{
    return (csec(n)-1.00)/2.00;
}

__MATHSUITE inline double complex  checsc(register double complex n)
{
    return (ccsc(n)-1.00)/2.00;
}

__MATHSUITE inline double complex  chesech(register double complex n)
{
    return (csech(n)-1.00)/2.00;
}

__MATHSUITE inline double complex  checsch(register double complex n)
{
    return (ccsch(n)-1.00)/2.00;
}

__MATHSUITE inline double complex  cqesec(register double complex n)
{
    return (csec(n)-1.00)/4.00;
}

__MATHSUITE inline double complex  cqecsc(register double complex n)
{
    return (ccsc(n)-1.00)/4.00;
}

__MATHSUITE inline double complex  cqesech(register double complex n)
{
    return (csech(n)-1.00)/4.00;
}

__MATHSUITE inline double complex  cqecsch(register double complex n)
{
    return (ccsch(n)-1.00)/4.00;
}

__MATHSUITE inline double complex  csinc(register double complex n)
{
    return n ? csin(n)/n : 1.00;
}

__MATHSUITE inline double complex  csinch(register double complex n)
{
    return n ? csinh(n)/n : 1.00;
}

__MATHSUITE inline double complex  chsinc(register double complex n)
{
    return n ? csin(n)/(2.00*n) : 0.50;
}

__MATHSUITE inline double complex  chsinch(register double complex n)
{
    return n ? csinh(n)/(2.00*n) : 0.50;
}

__MATHSUITE inline double complex  cqsinc(register double complex n)
{
    return n ? csin(n)/(4.00*n) : 0.25;
}

__MATHSUITE inline double complex  cqsinch(register double complex n)
{
    return n ? csinh(n)/(4.00*n) : 0.25;
}

__MATHSUITE inline double complex  ccosc(register double complex n)
{
    return ccos(n)/n;
}

__MATHSUITE inline double complex  ccosch(register double complex n)
{
    return ccosh(n)/n;
}

__MATHSUITE inline double complex  chcosc(register double complex n)
{
    return ccos(n)/(2.00*n);
}

__MATHSUITE inline double complex  chcosch(register double complex n)
{
    return ccosh(n)/(2.00*n);
}

__MATHSUITE inline double complex  cqcosc(register double complex n)
{
    return ccos(n)/(4.00*n);
}

__MATHSUITE inline double complex  cqcosch(register double complex n)
{
    return ccosh(n)/(4.00*n);
}

__MATHSUITE inline double complex  csecc(register double complex n)
{
    return csec(n)/n;
}

__MATHSUITE inline double complex  csecch(register double complex n)
{
    return csech(n)/n;
}

__MATHSUITE inline double complex  chsecc(register double complex n)
{
    return csec(n)/(2.00*n);
}

__MATHSUITE inline double complex  chsecch(register double complex n)
{
    return csech(n)/(2.00*n);
}

__MATHSUITE inline double complex  cqsecc(register double complex n)
{
    return csec(n)/(4.00*n);
}

__MATHSUITE inline double complex  cqsecch(register double complex n)
{
    return csech(n)/(4*n);
}

__MATHSUITE inline double complex  ccscc(register double complex n)
{
    return ccsc(n)/n;
}

__MATHSUITE inline double complex  ccscch(register double complex n)
{
    return ccsch(n)/n;
}

__MATHSUITE inline double complex  chcscc(register double complex n)
{
    return ccsc(n)/(2.00*n);
}

__MATHSUITE inline double complex  chcscch(register double complex n)
{
    return ccsch(n)/(2.00*n);
}

__MATHSUITE inline double complex  cqcscc(register double complex n)
{
    return ccsc(n)/(4.00*n);
}

__MATHSUITE inline double complex  cqcscch(register double complex n)
{
    return ccsch(n)/(4.00*n);
}

__MATHSUITE inline double complex  ctanc(register double complex n)
{
    return ctan(n)/n;
}

__MATHSUITE inline double complex  ctanch(register double complex n)
{
    return ctanh(n)/n;
}

__MATHSUITE inline double complex  chtanc(register double complex n)
{
    return ctan(n)/(2.00*n);
}

__MATHSUITE inline double complex  chtanch(register double complex n)
{
    return ctanh(n)/(2.00*n);
}

__MATHSUITE inline double complex  cqtanc(register double complex n)
{
    return ctan(n)/(4.00*n);
}

__MATHSUITE inline double complex  cqtanch(register double complex n)
{
    return ctanh(n)/(4.00*n);
}

__MATHSUITE inline double complex  ccotc(register double complex n)
{
    return ccot(n)/n;
}

__MATHSUITE inline double complex  ccotch(register double complex n)
{
    return ccoth(n)/n;
}

__MATHSUITE inline double complex  chcotc(register double complex n)
{
    return ccot(n)/(2.00*n);
}

__MATHSUITE inline double complex  chcotch(register double complex n)
{
    return ccoth(n)/(2.00*n);
}

__MATHSUITE inline double complex  cqcotc(register double complex n)
{
    return ccot(n)/(4.00*n);
}

__MATHSUITE inline double complex  cqcotch(register double complex n)
{
    return ccoth(n)/(4.00*n);
}


__MATHSUITE inline bool   isEqualMatrix(mpfr_t *matrix1, mpfr_t *matrix2, const register dim_typ dim[static 2])
{
    dim_typ i, j;

    for(i=0; i<dim[ROWS]; ++i)
        for(j=0; j<dim[COLUMNS]; ++j)
            if(mpfr_cmp(*(matrix1 + dim[COLUMNS]*i + j), *(matrix2 + dim[COLUMNS]*i + j)))
                return false;

    return true;

}

/// massive thanks to Bibek Subedi. Link at:
/// http://www.programming-technique.blogspot.it/2011/09/numerical-methods-condition-number-and.html
__MATHSUITE  void   norms(mpfr_t rop, mpfr_t *matrix, dim_typ dim)
{
    dim_typ i, j;
    mpfr_t tmp, tmp2;

    mpfr_init(tmp);
    mpfr_init_set_ui(tmp2, 2, MPFR_RNDN);
    mpfr_set_ui(rop, 0, MPFR_RNDN);
    for(i=0; i < dim; ++i)
        for(j = 0; j < dim; ++j)
        {
            mpow2(tmp, *(matrix + dim*i + j), tmp2);
            mpfr_add(rop, rop, tmp, MPFR_RNDN);
        }

    mpfr_sqrt(rop, rop, MPFR_RNDN);
    mpfr_clears(tmp, tmp2, NULL); 
    return;
}

/// thanks to: apatriarca, a User of Matematicamente forum, who is the author of this
/// pseudo-coded duet of functions which I renamed, unified, modified and adapted to this program.
/// This function calculates the Norm, One or Infinity (depending on mode bool param), of a Square Matrix.
/// http://www.matematicamente.it/forum/viewtopic.php?t=51498&p=371845
__MATHSUITE  void   norm(mpfr_t rop, mpfr_t *matrix, dim_typ dim, bool mode)
{
    dim_typ i, j;
    mpfr_t sum, tmp; //  = 0.00; /* tutte le somme sono >= 0 */

    mpfr_inits(sum, tmp, NULL); 
    mpfr_set_ui(rop, 0, MPFR_RNDN);
    for (j = 0; j < dim; ++j)
    {
        mpfr_set_ui(sum, 0, MPFR_RNDN);
        for (i = 0; i < dim; ++i)
        {
            mpfr_abs(tmp, *(matrix + dim*(mode?j:i) + (mode?i:j)), MPFR_RNDN);
            mpfr_add(sum, sum, tmp, MPFR_RNDN);
        }
        if (mpfr_cmp(sum, rop) > 0)
            mpfr_set(rop, sum, MPFR_RNDN);
    }

    mpfr_clears(sum, tmp, NULL); 
    return;
}

__MATHSUITE  void  newtonDifferenceTable(dim_typ n, mpfr_t y[access(curLayout)->max_newton_difftables_dim][access(curLayout)->max_newton_difftables_dim], bool mode)
{
    dim_typ i, j;
    if(mode)
    {
        for(j=0;++j<n;)
            for(i=0;i<(n-j);++i)
                mpfr_sub(y[i][j], y[i+1][j-1], y[i][j-1], MPFR_RNDN);
    }
    else
    {
        for(j=0;++j<n;)
            for(i=n-1;i>(j-1);--i)
                mpfr_sub(y[i][j], y[i][j-1], y[i-1][j-1], MPFR_RNDN);
    }

    return;
}

/*
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of n points.
   dir = true gives forward transform
   dir = false gives reverse transform
*/
__MATHSUITE  bool   FFT(ityp **x, register int n, const bool dir)
{
	int i,i1,j,k,i2,l,l1,l2;
	ityp c1,c2,tx,ty,t1,t2,u1,u2,z;

	/* Calculate the number of points */
	ityp tmp = log2(n);
	int m = (int) tmp;
	// Not a perfect 2 Power
	// so performing ZeroPadding
	if(tmp != m)
	{
		int oldn = n;
		(*x) = realloc((*x), ((n = exp2(++m))<<1)*sizeof(ityp));
		errMem((*x), FFT_ALLOC_ERROR);
		#pragma omp parallel for
		for(i=oldn; i<n; ++i)
			*((*x) + i) = *((*x) + n + i) = 0.00;
	}

	/* Do the bit reversal */
	i2 = n >> 1;
	j = 0;
	for (i=0;i<n-1;++i)
	{
		if(i < j)
		{
	    	tx = *((*x) + i);
	    	ty = *((*x) + n + i);
	     	*((*x) + i) = *((*x) + j);
	     	*((*x) + n + i) = *((*x) + n + j);
	     	*((*x) + j) = tx;
	    	*((*x) + n + j)= ty;
	 	}
		k = i2;
		while (k <= j)
		{
	    	j -= k;
	    	k >>= 1;
	 	}
		j += k;
	}

	/* Compute the FFT */
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	for(l=0;l<m;++l)
	{
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;
		for(j=0;j<l1;++j)
		{
			for(i=j;i<n;i+=l2)
			{
		        i1 = i + l1;
		        t1 = (u1 * *((*x) + i1)) - (u2 * *((*x) + n + i1));
		        t2 = (u1 * *((*x) + n + i1)) + (u2 * *((*x) + i1));
		        *((*x) + i1) = *((*x) + i) - t1;
		        *((*x) + n + i1) = *((*x) + n + i) - t2;
		        *((*x) + i) += t1;
		        *((*x) + n + i) += t2;
	     	}
		    	z =  u1 * c1 - u2 * c2;
		     	u2 = u1 * c2 + u2 * c1;
		     	u1 = z;
	  	}
	  	c2 = sqrt((1.0 - c1) / 2.0);
	  	if (dir)
	    	c2 *= -1;
	  	c1 = sqrt((1.0 + c1) / 2.0);
	}

	/* Scaling for forward transform */
	if(dir)
		#pragma omp parallel for
		for(i=0;i<n;++i)
	    	*((*x) + i) = ((*((*x) + n + i)) /= n);

	return FFT_SUCCESS;
}

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there inot powers of 2
*/
__MATHSUITE  bool   FFT2D(ityp **c, register int dim[static 2], const bool dir)
{
	int i,j;
	ityp *tmp = NULL;

	// Zero Padding
	ityp xy_tmp[2];
	xy_tmp[ROWS] = log2(dim[ROWS]);
	xy_tmp[COLUMNS] = log2(dim[COLUMNS]);
	int m[2];
	m[ROWS] = (int) xy_tmp[ROWS];
	m[COLUMNS] = (int) xy_tmp[COLUMNS];
	const bool xynp2[2] =
	{
		xy_tmp[ROWS] != m[ROWS],
		xy_tmp[COLUMNS] != m[COLUMNS]
	};
	if(xynp2[ROWS] || xynp2[COLUMNS])
	{
		const int oldnxy[2] =
		{
			dim[ROWS],
			dim[COLUMNS]
		};
		if(xynp2[ROWS])
			dim[ROWS] = exp2(++m[ROWS]);
		if(xynp2[COLUMNS])
			dim[COLUMNS] = exp2(++m[COLUMNS]);
		c[REAL_PART] = realloc(c[REAL_PART], dim[ROWS]*sizeof(ityp));
		c[IMAG_PART] = realloc(c[IMAG_PART], dim[COLUMNS]*sizeof(ityp));
		#pragma omp parallel for
		for(i=0; i<dim[ROWS]; ++i)
			#pragma omp parallel for
			for(j=oldnxy[COLUMNS];j<dim[COLUMNS]; ++j)
				c[REAL_PART][dim[COLUMNS]*i + j] = c[IMAG_PART][dim[COLUMNS]*i + j] = 0.00;

		#pragma omp parallel for
		for(i=oldnxy[ROWS]; i<dim[ROWS]; ++i)
			#pragma omp parallel for
			for(j=0; j<dim[COLUMNS]; ++j)
				c[REAL_PART][dim[COLUMNS]*i + j] = c[IMAG_PART][dim[COLUMNS]*i + j] = 0.00;
	}

	/* Transform the rows */

	tmp = malloc(sizeof(ityp)*(dim[ROWS]<<1));
	errMem(tmp, FFT_ALLOC_ERROR);

	for(j=0;j<dim[COLUMNS];++j)
	{
		#pragma omp parallel for
		for(i=0;i<dim[ROWS];++i)
		{
			*(tmp + i) = c[REAL_PART][dim[COLUMNS]*i + j];
			*(tmp + dim[ROWS] + i) = c[IMAG_PART][dim[COLUMNS]+i + j];
		}
		FFT(&tmp, m[ROWS], dir);
		#pragma omp parallel for
		for(i=0;i<dim[ROWS];++i)
		{
			c[REAL_PART][dim[COLUMNS]*i + j] = *(tmp + i);
			c[IMAG_PART][dim[COLUMNS]*i + j] = *(tmp + dim[ROWS] + i);
		}
	}

	/* Transform the columns */

	tmp = realloc(tmp, sizeof(ityp)*(dim[COLUMNS]<<1));
	errMem(tmp, FFT_ALLOC_ERROR);

	for(i=0;i<dim[ROWS];++i)
	{
		#pragma omp parallel for
		for(j=0;j<dim[COLUMNS];++j)
		{
			*(tmp + j) = c[REAL_PART][dim[COLUMNS]*i + j];
			*(tmp + dim[COLUMNS] + j) = c[IMAG_PART][dim[COLUMNS]*i + j];
		}
		FFT(&tmp, m[COLUMNS], dir);
		#pragma omp parallel for
		for(j=0;j<dim[COLUMNS];++j)
		{
			c[REAL_PART][dim[COLUMNS]*i + j] = *(tmp + j);
			c[IMAG_PART][dim[COLUMNS]*i + j] = *(tmp + dim[COLUMNS] + j);
		}
	}

	free(tmp);
	return FFT_SUCCESS;
}

__MATHSUITE inline void   harris(mpfr_t rop, const register mpfr_t alpha, const register dim_typ dim, mpfr_t *restrict a)
{
    mpfr_t tmp, tmp2, tmp3, tmp4;
    mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
    _matrixTrace(tmp4, a, dim);
    mpfr_pow_ui(tmp3, tmp4, 2, MPFR_RNDN);
    mpfr_mul(tmp2, alpha, tmp3, MPFR_RNDN);
    det(tmp, a, dim, NULL);
    mpfr_sub(rop, tmp, tmp2, MPFR_RNDN),
    mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
    return;
}

__MATHSUITE inline void   harris2(mpfr_t rop, const register mpfr_t alpha, mpfr_t *restrict a)
{
    mpfr_t tmp, tmp2, tmp3, tmp4;
    mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
    _matrixTrace(tmp4, a, 2);
    mpfr_pow_ui(tmp3, tmp4, 2, MPFR_RNDN);
    mpfr_mul(tmp2, alpha, tmp3, MPFR_RNDN);
    det(tmp, a, 2, NULL);
    mpfr_sub(rop, tmp, tmp2, MPFR_RNDN),
    mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
    return;
}

__MATHSUITE inline bool   t_harris(const register mpfr_t alpha, const register mpfr_t thresh, const register dim_typ dim, mpfr_t *restrict a)
{
    mpfr_t tmp, tmp2, tmp3, tmp4;
    mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
    _matrixTrace(tmp4, a, dim);
    mpfr_pow_ui(tmp3, tmp4, 2, MPFR_RNDN);
    mpfr_mul(tmp2, alpha, tmp3, MPFR_RNDN);
    det(tmp, a, dim, NULL);
    mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
    const bool retval = mpfr_cmp(tmp, thresh) >=0;
    mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
    return retval;
}

__MATHSUITE inline bool   t_harris2(const register mpfr_t alpha, const register mpfr_t thresh, mpfr_t *restrict a)
{
    mpfr_t tmp, tmp2, tmp3, tmp4;
    mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
    _matrixTrace(tmp4, a, 2);
    mpfr_pow_ui(tmp3, tmp4, 2, MPFR_RNDN);
    mpfr_mul(tmp2, alpha, tmp3, MPFR_RNDN);
    det(tmp, a, 2, NULL);
    mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
    const bool retval = mpfr_cmp(tmp, thresh) >=0;
    mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
    return retval;
}


__MATHSUITE inline void   eval(mpfr_t rop, mpfr_t *restrict a, const register dim_typ dim, const mpfr_t val)
{
    mpfr_t tmp;
    mpfr_init(tmp);
    mpfr_set(rop, *(a + dim-1), MPFR_RNDN);
	#pragma omp parallel for
	for(dim_typ i=0; i<dim-1; ++i)
	{
	    mpfr_pow_ui(tmp, val, dim-1-i, MPFR_RNDN);
        mpfr_mul(tmp, tmp, *(a + i), MPFR_RNDN);
        mpfr_add(rop, rop, tmp, MPFR_RNDN);
	}
	mpfr_clear(tmp);
	return;
}

__MATHSUITE inline void   deval(mpfr_t rop, mpfr_t *restrict a, const register dim_typ dim, const mpfr_t val)
{
    mpfr_t tmp, tmp2, tmp3;
    mpfr_init_set_ui(rop, 0, MPFR_RNDN);
    mpfr_init_set(tmp, *a, MPFR_RNDN);
    mpfr_inits(tmp2, tmp3, NULL); 
    mpfr_set_ui(*a, 0, MPFR_RNDN);
	for(dim_typ i=1; i<dim; ++i)
	{
	    mpfr_set(tmp2, *(a + i), MPFR_RNDN);
	    mpfr_set(*(a + i), tmp, MPFR_RNDN);
	    mpfr_mul_ui(*(a + i), *(a + i), dim-i, MPFR_RNDN);
        mpfr_pow_ui(tmp3, val, dim-1-i, MPFR_RNDN);
        mpfr_mul(tmp, *(a + i), tmp3, MPFR_RNDN);
	    mpfr_add(rop, rop, tmp, MPFR_RNDN);
        mpfr_set(tmp, tmp2, MPFR_RNDN);
	}
	mpfr_clears(tmp, tmp2, tmp3, NULL); 
	return;
}

__MATHSUITE short    _routhTable(mpfr_t **table, const register dim_typ dim, fsel_typ *nullrow)
{
	dim_typ i, j;
	(*table) = realloc((*table), sizeof(mpfr_t)*dim*dim);
	errMem((*table), ROUTHTABLE_ALLOC_ERROR);

	// cleaning before filling
	#pragma omp parallel for
	for(i=1; i<dim; ++i)
		#pragma omp parallel for
		for(j=0; j<dim; ++j)
            mpfr_set_ui(*((*table) + dim*i + j), 0, MPFR_RNDN);

	const register dim_typ columns = ((dim_typ)((dim*0.50) + 1));

	mpfr_t *tmp = NULL;

	if(!matrixAlloc(&tmp, (dim_typ2){1, columns}))
		return ROUTHTABLE_ALLOC_ERROR;

	for(j=0,i=1; i<dim; i+=2,++j)
        mpfr_set(*(tmp + j), *((*table) + i), MPFR_RNDN);
		// *((*table) + dim + j) = *((*table) + i);

	// works!!!
	for(j=1,i=2; i<dim; i+=2,++j)
		mpfr_set(*((*table) + j), *((*table) + i), MPFR_RNDN);

	// cleaning eventual residuals
	#pragma omp parallel for
	for(i=j; i<dim; ++i)
		mpfr_set_ui(*((*table) + i), 0, MPFR_RNDN);


	(*table) = realloc((*table), sizeof(mpfr_t)*dim*columns);
	errMem((*table), ROUTHTABLE_ALLOC_ERROR);

	#pragma omp parallel for
	for(i=0; i<columns; ++i)
		mpfr_set(*((*table) + columns + i), *(tmp + i), MPFR_RNDN);

	tmp = realloc(tmp, 4);
	errMem(tmp, ROUTHTABLE_ALLOC_ERROR);

	mpfr_t pivot;
	dim_typ q;

	mpfr_t mpftmp;

	for(i=2; i<dim; ++i)
	{
	    if(mpfr_zero_p(*((*table) + columns*(i-1))))
			for(j=1; j<columns; ++j)
                if(!mpfr_zero_p(*((*table) + columns*(i-1) + j)))
				{
					mpfr_set_d(*((*table) + columns*(i-1)), ROUTHTABLE_EPSILON, MPFR_RNDN); // 0.0000000001;
					break;
				}
				else if(j == columns-1)
				{
					if(nullrow)
						(*nullrow) = i-1;

					q = dim-i+1;
					for(dim_typ k=0; k<columns-1; ++k)
					{
					    mpfr_set(*((*table) + columns*(i-1) + k), *((*table) + columns*(i-2) + k), MPFR_RNDN);
					    mpfr_mul_ui(*((*table) + columns*(i-1) + k), *((*table) + columns*(i-1) + k), q, MPFR_RNDN);
						q -= 2;
						// shifting derivative of the [i-2] polynom into the substaying nullrow
						// by indexing the polynoms terms with the k variable
					}
				}

        mpfr_neg(pivot, pivot, MPFR_RNDN);
        mpfr_set(pivot, (*((*table) + columns*(i-1))), MPFR_RNDN);
        mpfr_pow_si(pivot, pivot, -1, MPFR_RNDN);
		mpfr_set(*tmp, *((*table) + columns*(i-2)), MPFR_RNDN);
		mpfr_set(*(tmp + 2), *((*table)  + columns*(i-1)), MPFR_RNDN);

		mpfr_init(mpftmp);
		for(j=0; j<columns-1; ++j)
        {
            mpfr_set(*(tmp + 1), *((*table) + columns*(i-2) + j+1), MPFR_RNDN);
            mpfr_set(*(tmp + 2 +1), *((*table) + columns*(i-1) + j+1), MPFR_RNDN);
            det(mpftmp, tmp, 2, NULL);
            mpfr_mul(*((*table) + columns*i + j),  mpftmp, pivot, MPFR_RNDN);
		}
	}

    mpfr_init(mpftmp);
	matrixFree(&tmp, ((dim_typ2){2, 2}));

	register short permanences = 0;

	#pragma omp parallel for
	for(i=0; i<dim-1; ++i)
    {
        mpfr_set(mpftmp, *((*table) + columns*i), MPFR_RNDN);
        mpfr_mul(mpftmp, mpftmp, *((*table) + columns*(i+1)), MPFR_RNDN);
        if(mpfGZero(mpftmp))
			++ permanences;
    }

    mpfr_clear(mpftmp);
	return permanences;
}

__MATHSUITE sel_typ    _juryTable(mpfr_t **table, const register dim_typ dim)
{
	dim_typ i, j;
	const register dim_typ mpfr_deg = dim-1;
	const register dim_typ rows = (mpfr_deg<<1)-3;

	(*table) = realloc((*table), sizeof(mpfr_t)*dim*rows);
	errMem((*table), JURYTABLE_ALLOC_ERROR);

	// cleaning before filling
	#pragma omp parallel for
	for(i=1; i<rows; ++i)
		for(j=0; j<dim; ++j)
			mpfr_set_ui(*((*table) + dim*i + j), 0, MPFR_RNDN);

	#pragma omp parallel for
	for(i=0; i<dim; ++i)
		mpfr_set(*((*table) + dim + i), *((*table) + i), MPFR_RNDN);

	#pragma omp parallel for
	for(i=0; i<dim; ++i)
		mpfr_set(*((*table) + i), *((*table) + (dim<<1) -1-i), MPFR_RNDN);

	mpfr_t *tmp = NULL;

	if(!matrixAlloc(&tmp, (dim_typ2){2, 2}))
		return JURYTABLE_ALLOC_ERROR;

	sel_typ jtest=0;
	mpfr_t mpftmp, mpftmp2;

	mpfr_inits(mpftmp, mpftmp2, NULL); 
	mpfr_abs(mpftmp, *(*table), MPFR_RNDN);

    if(mpfr_cmp(mpftmp, *((*table) + dim)) < 0)
		++ jtest;

	dim_typ k;
	// unrolled optimized loop
	for(k=1,i=3; i<rows; i+=2,++k)
	{
	    mpfr_set(*tmp, *((*table) + dim*(i-3)), MPFR_RNDN);
		mpfr_set(*(tmp + 2), *((*table) + dim*(i-2)), MPFR_RNDN);
		#pragma omp parallel for
		for(j=0; j<dim-k; ++j)
		{
			mpfr_set(*(tmp + 1), *((*table) + dim*(i-3) + j+1), MPFR_RNDN);
			mpfr_set(*(tmp + 2 +1), *((*table) + dim*(i-2) + j+1), MPFR_RNDN);
			mpfr_clear(mpftmp);
			det(mpftmp, tmp, 2, NULL);
			mpfr_set(*((*table) + (dim*i) + j), mpftmp, MPFR_RNDN);
			// if(i != rows-1)
            mpfr_set(*((*table) + (dim*(i-1)) + dim-k-j-1), *((*table) + (dim*i) + j), MPFR_RNDN);
		}

		mpfr_abs(mpftmp, *((*table) + (dim*(i-1))), MPFR_RNDN);
		mpfr_abs(mpftmp2, *((*table) + dim*i), MPFR_RNDN);
		if(mpfr_greater_p(mpftmp, mpftmp2))
			++ jtest;
	}


	mpfr_set(*tmp, *((*table) + dim*(rows-3)), MPFR_RNDN);
	mpfr_set(*(tmp + 2), *((*table) + dim*(rows-2)), MPFR_RNDN);
	#pragma omp parallel for
	for(i=0; i<3; ++i)
	{
		mpfr_set(*(tmp + 1), *((*table) + dim*(rows-3) + i+1), MPFR_RNDN);
		mpfr_set(*(tmp + 2 +1), *((*table) + dim*(rows-2) + i+1), MPFR_RNDN);
        det(mpftmp, tmp, 2, NULL);
		mpfr_set(*((*table) + dim*(rows-1) + 2-i), mpftmp, MPFR_RNDN);
	}

    mpfr_abs(mpftmp, *((*table) + dim*(rows-1)), MPFR_RNDN);
    mpfr_abs(mpftmp2, *((*table) + dim*(rows-1) + 2), MPFR_RNDN);
	if(mpfr_cmp(mpftmp, mpftmp2) > 0)
		++ jtest;

	matrixFree(&tmp, ((dim_typ2){2, 2}));
	mpfr_set_ui(mpftmp, 1, MPFR_RNDN);
	eval(mpftmp, *table, dim, mpftmp);
	mpfr_set_si(mpftmp2, -1, MPFR_RNDN);
	eval(mpftmp2, *table, dim, mpftmp2);
	const ityp binder = mpfr_get_d(mpftmp2, MPFR_RNDN)*(1-(((mpfr_deg%2)!=0)<<1));
	const bool sec_test = mpfGZero(mpftmp) && binder > 0;
	mpfr_clears(mpftmp, mpftmp2, NULL); 
	return jtest == mpfr_deg-1 && sec_test ? JURYTABLE_SATISFIED : JURYTABLE_NOTSATISFIED;
}

__MATHSUITE  sel_typ    _matrixEigenValues(mpfr_t *restrict m, mpfr_t *restrict l, mpfr_t *restrict vc, const register dim_typ n)
{
	mpfr_t *m2 = NULL;
	if(!matrixAlloc(&m2, (dim_typ2){n,n}))
		return EIGVALUES_ALLOC_ERROR;
	{
		register dim_typ i=n;
		while(i--)
		{
			register dim_typ j=n;
			while(j--)
			{
				mpfr_set(*(m2 + n*i + j), *(m + n*i + j), MPFR_RNDN);
				mpfr_set_ui(*(vc + n*i + j), i==j, MPFR_RNDN);
			}
		}
	}
	register dim_typ cnt = 0;
	while(true)
	{
		if(++cnt >= access(curLayout)->max_eigvalues_iterations)
			return EIGVALUES_INFEVS_ERROR;
		mpfr_t mod, q;
		mpfr_init(q);
		mpfr_init_set_ui(mod, 0, MPFR_RNDN);
		register dim_typ i=0, j=0;
		{
			register dim_typ k=n;
			while(k--)
			{
				register dim_typ m=n;
				while((--m)>k)
				{
					mpfr_abs(q, *(m2 + n*k + m), MPFR_RNDN);
  					if(mpfr_cmp(q,mod) > 0)
					{
					    mpfr_set(mod, q, MPFR_RNDN);
						i=k;
						j=m;
					}
				}
			}
		}
		// mpfr_clears(q);
		if(mpfr_cmp_d(mod, EIGENVALUES_PREC) < 0) break;
		mpfr_t th, tmp;
		mpfr_init(tmp);
		mpfr_init_set(th, *(m2 + n*i + j), MPFR_RNDN);
		mpfr_sub(tmp, *(m2 + n*i + i), *(m2 + n*j + j), MPFR_RNDN);
		mpfr_div(th, th, tmp, MPFR_RNDN);
		mpfr_mul_ui(th, th, 2, MPFR_RNDN);
		mpfr_atan(th, th, MPFR_RNDN);
		mpfr_div_ui(th, th, 2, MPFR_RNDN);
		{
		    mpfr_t c, s;
		    mpfr_inits(c, s, NULL); 
		    mpfr_cos(c, th, MPFR_RNDN);
		    mpfr_sin(s, th, MPFR_RNDN);
			inline void twst(mpfr_t *restrict m)
			{
			    mpfr_t t;
			    mpfr_init(t);
				register dim_typ k=n;
				while(k--)
				{
					mpfr_mul(t, *(m + n*i + k), c, MPFR_RNDN);
					mpfr_mul(q, *(m+ n*j + k), s, MPFR_RNDN);
					mpfr_add(t, t, q, MPFR_RNDN);
					mpfr_mul(q, *(m + n*i + k), s, MPFR_RNDN);
					mpfr_mul(tmp, *(m + n*j + k), c, MPFR_RNDN);
					mpfr_neg(q, q, MPFR_RNDN);
					mpfr_add(*(m + n*j + k), q, tmp, MPFR_RNDN);
					mpfr_set(*(m + n*i + k), t, MPFR_RNDN);
				}
				mpfr_clear(t);
			}
			{
				register dim_typ k=n;
				mpfr_t t;
				mpfr_init(t);
				while(k--)
				{
				    mpfr_mul(t, *(m2 + n*k + i), c, MPFR_RNDN);
				    mpfr_mul(q, *(m2 + n*k + j), s, MPFR_RNDN);
				    mpfr_add(t, t, q, MPFR_RNDN);
					mpfr_mul(*(m2 + n*k + j), *(m2 + n*k + i), s, MPFR_RNDN);
					mpfr_neg(*(m2 + n*k + j), *(m2 + n*k + j), MPFR_RNDN);
					mpfr_mul(q, *(m2 + n*k + j), c, MPFR_RNDN);
					mpfr_add(*(m2 + n*k + j), *(m2 + n*k + j), q, MPFR_RNDN);
					mpfr_set(*(m2 + n*k + i), t, MPFR_RNDN);
				}
				mpfr_clear(t);
			}
			twst(m2);
			twst(vc);
			mpfr_clears(c, s, NULL); 
		}
		mpfr_clears(q, tmp, th, mod, NULL); 
	}
	{
		register dim_typ j=n;
		while(j--) mpfr_set(l[j], *(m2 + n*j + j), MPFR_RNDN);
	}
	return EIGVALUES_FOUNDEVS;
}

__MATHSUITE inline void    _matrixSub(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_sum, const register dim_typ dim[static 2], const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
    	{
    		idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + j;
    		idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*i + j;
    		idx[MATRIX_SUM] = pitch[MATRIX_SUM]*i + j;
        	mpfr_sub(*((*matrix_sum) + idx[MATRIX_SUB]), *((*matrix1) + idx[FIRST_MATRIX]), *((*matrix2) + idx[SECOND_MATRIX]), MPFR_RNDN);
    	}

    return;
}

__MATHSUITE void    _matrixCSub(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_sum, const register dim_typ dim[static 2], const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
		{
			idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + j;
    		idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*i + j;
    		idx[MATRIX_SUM] = pitch[MATRIX_SUM]*i + j;
			const register double complex tmpres = (mpfr_get_d(matrix1[REAL_PART][idx[FIRST_MATRIX]], MPFR_RNDN) + mpfr_get_d(matrix1[IMAG_PART][idx[FIRST_MATRIX]], MPFR_RNDN)*I) - (mpfr_get_d(matrix2[REAL_PART][idx[SECOND_MATRIX]], MPFR_RNDN)+ mpfr_get_d(matrix2[IMAG_PART][idx[SECOND_MATRIX]], MPFR_RNDN)*I);
        	#pragma omp parallel sections num_threads(MAX_COMPLEX_UNITS)
       		{
       			#pragma omp section
	            mpfr_set_d(matrix_sum[REAL_PART][idx[MATRIX_SUB]], creal(tmpres), MPFR_RNDN);
	            #pragma omp section
	            mpfr_set_d(matrix_sum[IMAG_PART][idx[MATRIX_SUB]], cimag(tmpres), MPFR_RNDN);
	    	}
		}

    return;
}

__MATHSUITE void    _matrixQSub(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_sum, const register dim_typ dim[static 2], const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];


	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
    	{
    		idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + j;
    		idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*i + j;
    		idx[MATRIX_SUM] = pitch[MATRIX_SUM]*i + j;
        	#pragma omp parallel sections num_threads(MAX_QUATERNIONS_UNITS)
        	{
        		#pragma omp section
	            mpfr_sub(matrix_sum[QUATERNIONS_REALPART][idx[MATRIX_SUB]], matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	            #pragma omp section
				mpfr_sub(matrix_sum[QUATERNIONS_IPART][idx[MATRIX_SUB]], matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	            #pragma omp section
				mpfr_sub(matrix_sum[QUATERNIONS_JPART][idx[MATRIX_SUB]], matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	            #pragma omp section
				mpfr_sub(matrix_sum[QUATERNIONS_KPART][idx[MATRIX_SUB]], matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	    	}
		}

    return;
}

__MATHSUITE void    _matrixOSub(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_sum, const register dim_typ dim[static 2], const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
    	{
        	idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + j;
    		idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*i + j;
    		idx[MATRIX_SUM] = pitch[MATRIX_SUM]*i + j;
        	#pragma omp parallel sections num_threads(MAX_OCTONIONS_UNITS)
			{
				#pragma omp section
				mpfr_sub(matrix_sum[OCTONIONS_REALPART][idx[MATRIX_SUB]], matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
				#pragma omp section
				mpfr_sub(matrix_sum[OCTONIONS_E1PART][idx[MATRIX_SUB]], matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				#pragma omp section
				mpfr_sub(matrix_sum[OCTONIONS_E2PART][idx[MATRIX_SUB]], matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				#pragma omp section
				mpfr_sub(matrix_sum[OCTONIONS_E3PART][idx[MATRIX_SUB]], matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				#pragma omp section
				mpfr_sub(matrix_sum[OCTONIONS_E4PART][idx[MATRIX_SUB]], matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				#pragma omp section
				mpfr_sub(matrix_sum[OCTONIONS_E5PART][idx[MATRIX_SUB]], matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				#pragma omp section
				mpfr_sub(matrix_sum[OCTONIONS_E6PART][idx[MATRIX_SUB]], matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				#pragma omp section
				mpfr_sub(matrix_sum[OCTONIONS_E7PART][idx[MATRIX_SUB]], matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
			}
		}

    return;
}

__MATHSUITE void    _matrixSSub(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_sum, const register dim_typ dim[static 2], const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
        {
	        idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + j;
    		idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*i + j;
    		idx[MATRIX_SUM] = pitch[MATRIX_SUM]*i + j;
	        #pragma omp parallel sections num_threads(MAX_SEDENIONS_UNITS)
			{
				#pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_REALPART][idx[MATRIX_SUB]], matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E1PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E2PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E3PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E4PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E5PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E6PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E7PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E8PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E9PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E10PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E11PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E12PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E13PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E14PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_sub(matrix_sum[SEDENIONS_E15PART][idx[MATRIX_SUB]], matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
			}
		}

    return;
}

__MATHSUITE inline void    _matrixAdd(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_sum, const register dim_typ dim[static 2], const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
    	{
    		idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + j;
    		idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*i + j;
    		idx[MATRIX_SUM] = pitch[MATRIX_SUM]*i + j;
        	mpfr_add(*((*matrix_sum) + idx[MATRIX_SUM]), *((*matrix1) + idx[FIRST_MATRIX]), *((*matrix2) + idx[SECOND_MATRIX]), MPFR_RNDN);
    	}

    return;
}

__MATHSUITE void    _matrixCAdd(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_sum, const register dim_typ dim[static 2], const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
		{
			idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + j;
    		idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*i + j;
    		idx[MATRIX_SUM] = pitch[MATRIX_SUM]*i + j;
			const register double complex tmpres = (mpfr_get_d(matrix1[REAL_PART][idx[FIRST_MATRIX]], MPFR_RNDN) + mpfr_get_d(matrix1[IMAG_PART][idx[FIRST_MATRIX]], MPFR_RNDN)*I) + (mpfr_get_d(matrix2[REAL_PART][idx[SECOND_MATRIX]], MPFR_RNDN)+ mpfr_get_d(matrix2[IMAG_PART][idx[SECOND_MATRIX]], MPFR_RNDN)*I);
        	#pragma omp parallel sections num_threads(MAX_COMPLEX_UNITS)
       		{
       			#pragma omp section
	            mpfr_set_d(matrix_sum[REAL_PART][idx[MATRIX_SUM]], creal(tmpres), MPFR_RNDN);
	            #pragma omp section
	            mpfr_set_d(matrix_sum[IMAG_PART][idx[MATRIX_SUM]], cimag(tmpres), MPFR_RNDN);
	    	}
		}

    return;
}

__MATHSUITE void    _matrixQAdd(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_sum, const register dim_typ dim[static 2], const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];


	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
    	{
    		idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + j;
    		idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*i + j;
    		idx[MATRIX_SUM] = pitch[MATRIX_SUM]*i + j;
        	#pragma omp parallel sections num_threads(MAX_QUATERNIONS_UNITS)
        	{
        		#pragma omp section
	            mpfr_add(matrix_sum[QUATERNIONS_REALPART][idx[MATRIX_SUM]], matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], MPFR_RNDN);
	            #pragma omp section
				mpfr_add(matrix_sum[QUATERNIONS_IPART][idx[MATRIX_SUM]], matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_IPART][idx[FIRST_MATRIX]], MPFR_RNDN);
	            #pragma omp section
				mpfr_add(matrix_sum[QUATERNIONS_JPART][idx[MATRIX_SUM]], matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_JPART][idx[FIRST_MATRIX]], MPFR_RNDN);
	            #pragma omp section
				mpfr_add(matrix_sum[QUATERNIONS_KPART][idx[MATRIX_SUM]], matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_KPART][idx[FIRST_MATRIX]], MPFR_RNDN);
	    	}
		}

    return;
}

__MATHSUITE void    _matrixOAdd(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_sum, const register dim_typ dim[static 2], const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
        {
	        idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + j;
    		idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*i + j;
    		idx[MATRIX_SUM] = pitch[MATRIX_SUM]*i + j;
		    #pragma omp parallel sections num_threads(MAX_OCTONIONS_UNITS)
			{
				#pragma omp section
		        mpfr_add(matrix_sum[OCTONIONS_REALPART][idx[MATRIX_SUM]], matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[OCTONIONS_E1PART][idx[MATRIX_SUM]], matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[OCTONIONS_E2PART][idx[MATRIX_SUM]], matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[OCTONIONS_E3PART][idx[MATRIX_SUM]], matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[OCTONIONS_E4PART][idx[MATRIX_SUM]], matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[OCTONIONS_E5PART][idx[MATRIX_SUM]], matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[OCTONIONS_E6PART][idx[MATRIX_SUM]], matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[OCTONIONS_E7PART][idx[MATRIX_SUM]], matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
			}
		}

    return;
}

__MATHSUITE void    _matrixSAdd(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_sum, const register dim_typ dim[static 2], const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
        {
        	idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + j;
    		idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*i + j;
    		idx[MATRIX_SUM] = pitch[MATRIX_SUM]*i + j;
	        #pragma omp parallel sections num_threads(MAX_SEDENIONS_UNITS)
			{
				#pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_REALPART][idx[MATRIX_SUM]], matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E1PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E2PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E3PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E4PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E5PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E6PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E7PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E8PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E9PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E10PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E11PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E12PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E13PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E14PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
		        #pragma omp section
				mpfr_add(matrix_sum[SEDENIONS_E15PART][idx[MATRIX_SUM]], matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
			}
		}

    return;
}

 inline void   squareOSMM(void (* const prodFunc)(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static MAX_MATRICES], const register dim_typ pitch[static MAX_MATRICES]), mpfr_t **A, mpfr_t **B, mpfr_t **C, const register dim_typ dim, const register dim_typ pitch[static MAX_MATRICES])
{
	dim_typ i, j, k;
	const register dim_typ BLOCK_SIZE = access(curLayout)->block_size;
	for (i = 0; i < dim; i += BLOCK_SIZE)
		for (j = 0; j < dim; j += BLOCK_SIZE)
    		for (k = 0; k < dim; k += BLOCK_SIZE)
				/* Correct block dimensions if block "goes off edge of" the matrix */
				/* Perform individual block osmm */
				prodFunc(A + i + k*dim, B + k + j*dim, C + i + j*dim, (dim_typ3){MIN (BLOCK_SIZE, dim-i), MIN (BLOCK_SIZE, dim-j), MIN (BLOCK_SIZE, dim-k)}, pitch);
}

//
// Volker Strassen algorithm for matrix multiplication.
// Theoretical Runtime is O(N^2.807).
// Assume NxN matrices where N is a power of two.
// Algorithm:
//   Matrices X and Y are split into four smaller
//  (N/2)x(N/2) matrices as follows:
//          _    _          _   _
//     X = | A  B |    Y = | E F |
//         | C  D |        | G H |
//          -    -          -   -
//   Then we build the following 7 matrices (requiring
//   seven (N/2)x(N/2) matrix multiplications -- this is
//   where the 2.807 = log2(7) improvement comes from):
//     P0 = A*(F - H);
//     P1 = (A + B)*H
//     P2 = (C + D)*E
//     P3 = D*(G - E);
//     P4 = (A + D)*(E + H)
//     P5 = (B - D)*(G + H)
//     P6 = (A - C)*(E + F)
//   The final result is
//        _                                            _
//   Z = | (P3 + P4) + (P5 - P1)   P0 + P1              |
//       | P2 + P3              (P0 + P4) - (P2 + P6)|
//        -                                            -
//

__MATHSUITE  void   mmult_fast(const register dim_typ N, const register dim_typ pitch[static MAX_MATRICES], const register sel_typ algunits, mpfr_t **X, mpfr_t **Y,  mpfr_t **Z,
	void (* const matrixMulFunc)(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static MAX_MATRICES], const register dim_typ [static MAX_MATRICES]),
	void (* const matrixAddFunc)(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]),
	void (* const matrixSubFunc)(mpfr_t **, mpfr_t **, mpfr_t **, const register dim_typ [static 2], const register dim_typ [static MAX_MATRICES]))
{
	//
	// Recursive base case.
	// If matrices are 16x16 or smaller we just use
	// the conventional algorithm.
	// At what size we should switch will vary based
	// on hardware platform.
	//
	if (N <= access(curLayout)->min_strassen_dim)
	{
		if (N <= access(curLayout)->min_osmm_dim)
			matrixMulFunc(X,Y,Z,(dim_typ3){N,N,N}, pitch);
		else
			squareOSMM(matrixMulFunc, X,Y,Z, N, pitch);  // matrixMulFunc(X,Y,Z,(dim_typ3){N,N,N});
		return;
	}

	dim_typ i, j;
	const register dim_typ n = N/2;      // size of sub-matrices
	const register dim_typ sz = n*n*sizeof(mpfr_t); // )n*n*sizeof(ityp*);
	const register dim_typ sz2 = sizeof(mpfr_t*)*algunits;

	mpfr_t **A = X;    // A-D matrices embedded in X
	
	mpfr_t *B[algunits];
	mpfr_t *C[algunits];
	mpfr_t *D[algunits];

	mpfr_t **E = Y;    // E-H matrices embeded in Y
	mpfr_t *F[algunits];
	mpfr_t *G[algunits];
	mpfr_t *H[algunits];
	
	//for (i = 0; i < 7; i++)
		//P[i] = malloc(sz);
		
	mpfr_t *P[7][algunits];
	
	mpfr_t *Z0[algunits];
	mpfr_t *Z1[algunits];
	mpfr_t *Z2[algunits];
	mpfr_t **T;
	mpfr_t **U;
	
	dim_typ2 dim =
	{
		n,
		n
	};

	#pragma omp parallel for num_threads(algunits)
	for(i=0; i<algunits; ++i)
	{
		B[i] = X[i] + n;
		C[i] = X[i] + n*pitch[FIRST_MATRIX];
		D[i] = C[i] + n;
		F[i] = Y[i] + n;
		G[i] = Y[i] + n*pitch[SECOND_MATRIX];
		H[i] = G[i] + n;
		Z0[i] = Z[i] + n*pitch[THIRD_MATRIX];
		Z1[i] = Z[i] + n;
		Z2[i] = Z[i]+n*(pitch[THIRD_MATRIX]+1);
		#pragma omp parallel for num_threads(7) 
		for(j=0; j<7; ++j)
			matrixAlloc(&P[j][i], dim);
	}

	#pragma omp parallel sections num_threads(7) private(i, T, U)
	{	
		#pragma omp section
		{
			T = malloc(sz2);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
				matrixAlloc(&T[i], dim);
			// P0 = A*(F - H);
			matrixSubFunc(F,H,T,dim,(dim_typ3){pitch[SECOND_MATRIX], pitch[SECOND_MATRIX], n});
			mmult_fast(n,(dim_typ3){pitch[FIRST_MATRIX], n, n},algunits,A,T,P[0],matrixMulFunc,matrixAddFunc,matrixSubFunc);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
				matrixFree(&T[i], dim);
			free(T);
		}
		
		#pragma omp section
		{
			T = malloc(sz2);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
				matrixAlloc(&T[i], dim);
			// P1 = (A + B)*H
			matrixAddFunc(A,B,T,dim,(dim_typ3){pitch[FIRST_MATRIX], pitch[FIRST_MATRIX], n});
			mmult_fast(n,(dim_typ3){n, pitch[SECOND_MATRIX], n},algunits,T,H,P[1],matrixMulFunc,matrixAddFunc,matrixSubFunc);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
				matrixFree(&T[i], dim);
			free(T);
		}
		
		#pragma omp section
		{
			T = malloc(sz2);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
				matrixAlloc(&T[i], dim);
			// P2 = (C + D)*E
			matrixAddFunc(C,D,T,dim,(dim_typ3){pitch[FIRST_MATRIX], pitch[FIRST_MATRIX], n});
			mmult_fast(n, (dim_typ3){n, pitch[SECOND_MATRIX], n},algunits, T, E, P[2],matrixMulFunc,matrixAddFunc,matrixSubFunc);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
				matrixFree(&T[i], dim);
			free(T);
		}
		
		#pragma omp section
		{
			T = malloc(sz2);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
				matrixAlloc(&T[i], dim);
			// P3 = D*(G - E);
			matrixSubFunc(G,E,T,dim,(dim_typ3){pitch[SECOND_MATRIX], pitch[SECOND_MATRIX], n});
			mmult_fast(n,(dim_typ3){pitch[FIRST_MATRIX], n, n},algunits,D,T,P[3],matrixMulFunc,matrixAddFunc,matrixSubFunc);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
				matrixFree(&T[i], dim);
			free(T);
		}
		
		#pragma omp section
		{
			T = malloc(sz2);
			U = malloc(sz2);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
			{
				matrixAlloc(&T[i], dim);
				matrixAlloc(&U[i], dim);
			}
			// P4 = (A + D)*(E + H)
			matrixAddFunc(A,D,T,dim,(dim_typ3){pitch[FIRST_MATRIX], pitch[FIRST_MATRIX], n});
			matrixAddFunc(E,H,U,dim,(dim_typ3){pitch[SECOND_MATRIX], pitch[SECOND_MATRIX], n});
			mmult_fast(n,(dim_typ3){n,n,n},algunits,T,U,P[4],matrixMulFunc,matrixAddFunc,matrixSubFunc);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
			{
				matrixFree(&U[i], dim);
				matrixFree(&T[i], dim);
			}
			free(U);
			free(T);
		}
		
		#pragma omp section
		{
			T = malloc(sz2);
			U = malloc(sz2);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
			{
				matrixAlloc(&T[i], dim);
				matrixAlloc(&U[i], dim);
			}
			// P5 = (B - D)*(G + H)
			matrixSubFunc(B,D,T,dim,(dim_typ3){pitch[FIRST_MATRIX], pitch[FIRST_MATRIX], n});
			matrixAddFunc(G,H,U,dim,(dim_typ3){pitch[SECOND_MATRIX], pitch[SECOND_MATRIX], n});
			mmult_fast(n,(dim_typ3){n,n,n},algunits,T,U,P[5],matrixMulFunc,matrixAddFunc,matrixSubFunc);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
			{
				matrixFree(&U[i], dim);
				matrixFree(&T[i], dim);
			}
			free(U);
			free(T);
		}
		
		#pragma omp section
		{
			T = malloc(sz2);
			U = malloc(sz2);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
			{
				matrixAlloc(&T[i], dim);
				matrixAlloc(&U[i], dim);
			}
			// P6 = (A - C)*(E + F)
			matrixSubFunc(A,C,T,dim,(dim_typ3){pitch[FIRST_MATRIX], pitch[FIRST_MATRIX], n});
			matrixAddFunc(E,F,U,dim,(dim_typ3){pitch[SECOND_MATRIX], pitch[SECOND_MATRIX], n});
			mmult_fast(n,(dim_typ3){n,n,n},algunits,T,U,P[6],matrixMulFunc,matrixAddFunc,matrixSubFunc);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
			{
				matrixFree(&U[i], dim);
				matrixFree(&T[i], dim);
			}
			free(U);
			free(T);
		}
	}
	
	#pragma omp parallel sections num_threads(4) private(i, T, U)
	{
		#pragma omp section
		{
			T = malloc(sz2);
			U = malloc(sz2);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
			{
				matrixAlloc(&T[i], dim);
				matrixAlloc(&U[i], dim);
			}
			// Z upper left = (P3 + P4) + (P5 - P1)
			matrixAddFunc(P[4], P[3], T, dim, (dim_typ3){n,n,n});
			matrixSubFunc(P[5], P[1], U, dim, (dim_typ3){n,n,n});
			matrixAddFunc(T,U,Z,dim, (dim_typ3){n,n,pitch[THIRD_MATRIX]});
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
			{
				matrixFree(&U[i], dim);
				matrixFree(&T[i], dim);
			}
			free(U);
			free(T);
		}
		
		#pragma omp section
			// Z lower left = P2 + P3
			matrixAddFunc(P[2],P[3], Z0, dim, (dim_typ3){n,n,pitch[THIRD_MATRIX]});
		
		#pragma omp section
			// Z upper right = P0 + P1
			matrixAddFunc(P[0], P[1], Z1, dim, (dim_typ3){n,n,pitch[THIRD_MATRIX]});
		
		#pragma omp section
		{
			T = malloc(sz2);
			U = malloc(sz2);
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
			{
				matrixAlloc(&T[i], dim);
				matrixAlloc(&U[i], dim);
			}
			// Z lower right = (P0 + P4) - (P2 + P6)
			matrixAddFunc(P[0], P[4], T, dim, (dim_typ3){n,n,n});
			matrixAddFunc(P[2], P[6], U, dim, (dim_typ3){n,n,n});
			matrixSubFunc(T,U,Z2, dim, (dim_typ3){n,n,pitch[THIRD_MATRIX]});
			#pragma omp parallel for num_threads(algunits)
			for(i=0; i<algunits; ++i)
			{
				matrixFree(&U[i], dim);
				matrixFree(&T[i], dim);
			}
			free(U);
			free(T);
		}
		
		
	}
	
	#pragma omp parallel for num_threads(7)
		for(i=0; i<7; ++i)
			#pragma omp parallel for num_threads(algunits)
			for(j=0; j<algunits; ++j)
				matrixFree(&P[i][j], dim);
				
	return;
}

__MATHSUITE inline void    _matrixMultiplication(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_product, const register dim_typ dim[static MAX_MATRICES], const register dim_typ pitch[static MAX_MATRICES]) // dim_typ righe, dim_typ colonne, dim_typ colonne2)
{
    dim_typ i, j, k;
    register dim_typ idx[MAX_MATRICES];
    mpfr_t tmp;

    mpfr_init(tmp);
    for(i=0; i<dim[ROWS]; ++i)
        for(k=0; k<dim[COLUMNS]; ++k)
       		for(j=0; j<dim[COLUMNS2]; ++j)
            {
            	#pragma omp parallel sections num_threads(MAX_MATRICES)
        		{
        			#pragma omp section
	        		idx[FIRST_MATRIX] = pitch[FIRST_MATRIX]*i + k;
	        		#pragma omp section
	            	idx[SECOND_MATRIX] = pitch[SECOND_MATRIX]*k + j;
	            	#pragma omp section
	            	idx[MATRIX_PRODUCT] = pitch[MATRIX_PRODUCT]*i + j;
	            }
	            
                mpfr_mul(tmp, *((*matrix1) + idx[FIRST_MATRIX]), *((*matrix2) + idx[SECOND_MATRIX]), MPFR_RNDN);
                mpfr_add(*((*matrix_product) + idx[MATRIX_PRODUCT]), *((*matrix_product) + idx[MATRIX_PRODUCT]), tmp, MPFR_RNDN);
            }

    mpfr_clear(tmp);
    return;
}

__MATHSUITE void    _matrixCMultiplication(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_product, const register dim_typ dim[static MAX_MATRICES], const register dim_typ pitch[static MAX_MATRICES])
{
    dim_typ i, j, k, x;
    dim_typ idx[MAX_MATRICES];

    for(i=0; i<dim[ROWS]; ++i)
        for(k=0; k<dim[COLUMNS]; ++k)
        	for(j=0; j<dim[COLUMNS2]; ++j)
        	{
        		#pragma omp parallel for num_threads(MAX_MATRICES)
        		for(x=0;x<MAX_MATRICES; ++x)
	        		idx[x] = pitch[x]*i + k;
            	const register double complex tmpres = (mpfr_get_d(matrix1[REAL_PART][idx[FIRST_MATRIX]], MPFR_RNDN) + mpfr_get_d(matrix1[IMAG_PART][idx[FIRST_MATRIX]], MPFR_RNDN)*I) * (mpfr_get_d(matrix2[REAL_PART][idx[SECOND_MATRIX]], MPFR_RNDN) + mpfr_get_d(matrix2[IMAG_PART][idx[SECOND_MATRIX]], MPFR_RNDN)*I);
				#pragma omp parallel sections num_threads(MAX_COMPLEX_UNITS) 
				{
					#pragma omp section
					mpfr_add_d(matrix_product[REAL_PART][idx[MATRIX_PRODUCT]], matrix_product[REAL_PART][idx[MATRIX_PRODUCT]], creal(tmpres), MPFR_RNDN);
		            #pragma omp section
					mpfr_add_d(matrix_product[IMAG_PART][idx[MATRIX_PRODUCT]], matrix_product[IMAG_PART][idx[MATRIX_PRODUCT]], cimag(tmpres), MPFR_RNDN);
	     		}
		   	}

    return;
}

__MATHSUITE void    _matrixQMultiplication(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_product, const register dim_typ dim[static MAX_MATRICES], const register dim_typ pitch[static MAX_MATRICES])
{
    dim_typ i, j, k, x;
    dim_typ idx[MAX_MATRICES];

    for(i=0; i<dim[ROWS]; ++i)
        for(k=0; k<dim[COLUMNS]; ++k)
        	for(j=0; j<dim[COLUMNS2]; ++j)
        	{
        		#pragma omp parallel for num_threads(MAX_MATRICES)
        		for(x=0;x<MAX_MATRICES; ++x)
	        		idx[x] = pitch[x]*i + k;
            	#pragma omp parallel sections num_threads(MAX_QUATERNIONS_UNITS)
            	{
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
	                    mpfr_mul(tmp, matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                    mpfr_mul(tmp2, matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                    mpfr_mul(tmp3, matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                    mpfr_mul(tmp4, matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(matrix_product[QUATERNIONS_REALPART][idx[MATRIX_PRODUCT]], matrix_product[QUATERNIONS_REALPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
					}
					
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
		                mpfr_mul(tmp, matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_mul(tmp2, matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_mul(tmp3, matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_mul(tmp4, matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
		                mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
		                mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
		                mpfr_add(matrix_product[QUATERNIONS_IPART][idx[MATRIX_PRODUCT]], matrix_product[QUATERNIONS_IPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
					}
	                	
	                #pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
						mpfr_mul(tmp, matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_mul(tmp2, matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_mul(tmp3, matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_mul(tmp4, matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
		                mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
		                mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
		                mpfr_add(matrix_product[QUATERNIONS_JPART][idx[MATRIX_PRODUCT]], matrix_product[QUATERNIONS_JPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
		                mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
					}

					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
		                mpfr_mul(tmp, matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_mul(tmp2, matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_mul(tmp3, matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_mul(tmp4, matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]], MPFR_RNDN);
		                mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
		                mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
		                mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
		                mpfr_add(matrix_product[QUATERNIONS_KPART][idx[MATRIX_PRODUCT]], matrix_product[QUATERNIONS_KPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
		                mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
					}
	            }
	    	}

    return;
}

__MATHSUITE void    _matrixOMultiplication(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_product, const register dim_typ dim[static MAX_MATRICES], const register dim_typ pitch[static MAX_MATRICES])
{
    dim_typ i, j, k, x;
    mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
    dim_typ idx[MAX_MATRICES];

    mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
    for(i=0; i<dim[ROWS]; ++i)
        for(k=0; k<dim[COLUMNS]; ++k)
        	for(j=0; j<dim[COLUMNS2]; ++j)
        	{
				#pragma omp parallel for num_threads(MAX_MATRICES)
        		for(x=0;x<MAX_MATRICES; ++x)
	        		idx[x] = pitch[x]*i + k;
	        	#pragma omp parallel sections num_threads(MAX_OCTONIONS_UNITS)
				{
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					    mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(matrix_product[OCTONIONS_REALPART][idx[MATRIX_PRODUCT]], matrix_product[OCTONIONS_REALPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					}
			        		
			        #pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
						mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(matrix_product[OCTONIONS_E1PART][idx[MATRIX_PRODUCT]], matrix_product[OCTONIONS_E1PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                    mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					}

					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
				        mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(matrix_product[OCTONIONS_E2PART][idx[MATRIX_PRODUCT]], matrix_product[OCTONIONS_E2PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					}
			       
			       	#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					    mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(matrix_product[OCTONIONS_E3PART][idx[MATRIX_PRODUCT]], matrix_product[OCTONIONS_E3PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					}
			        
			        #pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
						mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(matrix_product[OCTONIONS_E4PART][idx[MATRIX_PRODUCT]], matrix_product[OCTONIONS_E4PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					}
			        
			        #pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
						mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(matrix_product[OCTONIONS_E5PART][idx[MATRIX_PRODUCT]], matrix_product[OCTONIONS_E5PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					}
			        
			        #pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
						mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(matrix_product[OCTONIONS_E6PART][idx[MATRIX_PRODUCT]], matrix_product[OCTONIONS_E6PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                    mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					}

					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
				        mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(matrix_product[OCTONIONS_E7PART][idx[MATRIX_PRODUCT]], matrix_product[OCTONIONS_E7PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                    mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					}
				}
			}

    return;
}

__MATHSUITE void    _matrixSMultiplication(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_product, const register dim_typ dim[static 3], const register dim_typ pitch[static MAX_MATRICES])
{
    dim_typ i, j, k, x;
	dim_typ idx[MAX_MATRICES];

    for(i=0; i<dim[ROWS]; ++i)
        for(k=0; k<dim[COLUMNS]; ++k)
        	for(j=0; j<dim[COLUMNS2]; ++j)
        	{
				#pragma omp parallel for num_threads(MAX_MATRICES)
        		for(x=0;x<MAX_MATRICES; ++x)
	        		idx[x] = pitch[x]*i + k;
        		#pragma omp parallel sections num_threads(MAX_SEDENIONS_UNITS)
				{
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
					    mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                    mpfr_add(matrix_product[SEDENIONS_REALPART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_REALPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
			        
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
						mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                    mpfr_add(matrix_product[SEDENIONS_E1PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E1PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
					
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
				        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                    mpfr_add(matrix_product[SEDENIONS_E2PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E2PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                    mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
		
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
						mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
				        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                    mpfr_add(matrix_product[SEDENIONS_E3PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E3PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                    mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
	
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
						mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                    mpfr_add(matrix_product[SEDENIONS_E4PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E4PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                    mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
	
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
				        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
				        mpfr_add(matrix_product[SEDENIONS_E5PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E5PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
				        
				    #pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
						mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
				        mpfr_add(matrix_product[SEDENIONS_E6PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E6PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
				        
				    #pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
						mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
				        mpfr_add(matrix_product[SEDENIONS_E7PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E7PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
				        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
	
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
				        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
				        mpfr_add(matrix_product[SEDENIONS_E8PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E8PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
				        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
		
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
				        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
				        mpfr_add(matrix_product[SEDENIONS_E9PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E9PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
				        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
	
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
				        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                    mpfr_add(matrix_product[SEDENIONS_E10PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E10PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                    mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
	
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
				        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
				        mpfr_add(matrix_product[SEDENIONS_E11PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E11PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
						mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
				        
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
						mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                    mpfr_add(matrix_product[SEDENIONS_E12PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E12PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                    mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
		
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
				        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                    mpfr_add(matrix_product[SEDENIONS_E13PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E13PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                    mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
				
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
				        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                    mpfr_add(matrix_product[SEDENIONS_E14PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E14PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                    mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
	
					#pragma omp section
					{
						mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
				        mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
				        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
					    mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
				        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
					    mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
					    mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                    mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                    mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
				        mpfr_add(matrix_product[SEDENIONS_E15PART][idx[MATRIX_PRODUCT]], matrix_product[SEDENIONS_E15PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
				        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
					}
				}
			}
    return;
}

__MATHSUITE inline void    _matrixKProduct(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_product, register dim_typ dim[static 2][2])
{
    dim_typ i, j, k, l;
    const register dim_typ pdim = dim[FIRST_MATRIX][COLUMNS]*dim[SECOND_MATRIX][COLUMNS];

	#pragma omp parallel for
    for(i=0; i<dim[FIRST_MATRIX][ROWS]; ++i)
    	#pragma omp parallel for
        for(k=0; k<dim[SECOND_MATRIX][ROWS]; ++k)
			#pragma omp parallel for
        	for(j=0; j<dim[FIRST_MATRIX][COLUMNS]; ++j)
            	#pragma omp parallel for
                for(l=0; l<dim[SECOND_MATRIX][COLUMNS]; ++l)
                	mpfr_mul(*((*matrix_product) + pdim*(k+(i*dim[SECOND_MATRIX][ROWS])) + l+(j*dim[SECOND_MATRIX][COLUMNS])), *((*matrix1) + dim[FIRST_MATRIX][COLUMNS]*i + j), *((*matrix2) + dim[SECOND_MATRIX][COLUMNS]*k + l), MPFR_RNDN);

    return;
}

__MATHSUITE void    _matrixKCProduct(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_product, register dim_typ dim[static 2][2])
{
    dim_typ i, j, k, l;
    register dim_typ idx[MAX_MATRICES];
    const register dim_typ pdim = dim[FIRST_MATRIX][COLUMNS]*dim[SECOND_MATRIX][COLUMNS];

	#pragma omp parallel for
    for(i=0; i<dim[FIRST_MATRIX][ROWS]; ++i)
    	#pragma omp parallel for
        	for(k=0; k<dim[SECOND_MATRIX][ROWS]; ++k)
    		#pragma omp parallel for
        	for(j=0; j<dim[FIRST_MATRIX][COLUMNS]; ++j)
            	#pragma omp parallel for
                for(l=0; l<dim[SECOND_MATRIX][COLUMNS]; ++l)
            	{
            		idx[FIRST_MATRIX] = dim[FIRST_MATRIX][COLUMNS]*i + j;
                	idx[SECOND_MATRIX] = dim[SECOND_MATRIX][COLUMNS]*k + l;
                	idx[MATRIX_PRODUCT] = pdim*(k+(i*dim[SECOND_MATRIX][ROWS])) + l+(j*dim[SECOND_MATRIX][COLUMNS]);
            		const register double complex tmpres = (mpfr_get_d(matrix1[REAL_PART][idx[FIRST_MATRIX]], MPFR_RNDN) + mpfr_get_d(matrix1[IMAG_PART][idx[FIRST_MATRIX]], MPFR_RNDN)*I) * (mpfr_get_d(matrix2[REAL_PART][idx[SECOND_MATRIX]], MPFR_RNDN) + mpfr_get_d(matrix2[IMAG_PART][idx[SECOND_MATRIX]], MPFR_RNDN)*I);
                    #pragma omp parallel sections num_threads(MAX_COMPLEX_UNITS)
                    {
                    	#pragma omp section
						mpfr_set_d(matrix_product[REAL_PART][idx[MATRIX_PRODUCT]], creal(tmpres), MPFR_RNDN);
	                    #pragma omp section
						mpfr_set_d(matrix_product[IMAG_PART][idx[MATRIX_PRODUCT]], cimag(tmpres), MPFR_RNDN);
	                }
	        	}

    return;
}

__MATHSUITE void    _matrixKQProduct(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_product, register dim_typ dim[static 2][2])
{
    dim_typ i, j, k, l;
    register dim_typ idx[MAX_MATRICES];
    const register dim_typ pdim = dim[FIRST_MATRIX][COLUMNS]*dim[SECOND_MATRIX][COLUMNS];

	#pragma omp parallel for
    for(i=0; i<dim[FIRST_MATRIX][ROWS]; ++i)
    	#pragma omp parallel for
        for(k=0; k<dim[SECOND_MATRIX][ROWS]; ++k)
    		#pragma omp parallel for
        	for(j=0; j<dim[FIRST_MATRIX][COLUMNS]; ++j)
            	#pragma omp parallel for
                for(l=0; l<dim[SECOND_MATRIX][COLUMNS]; ++l)
                {
                	idx[FIRST_MATRIX] = dim[FIRST_MATRIX][COLUMNS]*i + j;
                	idx[SECOND_MATRIX] = dim[SECOND_MATRIX][COLUMNS]*k + l;
                	idx[MATRIX_PRODUCT] = pdim*(k+(i*dim[SECOND_MATRIX][ROWS])) + l+(j*dim[SECOND_MATRIX][COLUMNS]);
                	#pragma omp parallel sections num_threads(MAX_QUATERNIONS_UNITS)
                	{
						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
		                    mpfr_mul(tmp, matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_set(matrix_product[QUATERNIONS_REALPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
	                	}
	
						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
	                        mpfr_mul(tmp, matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_set(matrix_product[QUATERNIONS_IPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL); 
	                        mpfr_mul(tmp, matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_set(matrix_product[QUATERNIONS_JPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, NULL);
	                        mpfr_mul(tmp, matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]], matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_set(matrix_product[QUATERNIONS_KPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, NULL); 
	                	}
	    			}
	    		}
 
    return;
}

__MATHSUITE void    _matrixKOProduct(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_product, register dim_typ dim[static 2][2])
{
    dim_typ i, j, k, l;
    register dim_typ idx[MAX_MATRICES];
    const register dim_typ pdim = dim[FIRST_MATRIX][COLUMNS]*dim[SECOND_MATRIX][COLUMNS];

    #pragma omp parallel for
    for(i=0; i<dim[FIRST_MATRIX][ROWS]; ++i)
    	#pragma omp parallel for
        for(k=0; k<dim[SECOND_MATRIX][ROWS]; ++k)
    		#pragma omp parallel for
       	 	for(j=0; j<dim[FIRST_MATRIX][COLUMNS]; ++j)
            	#pragma omp parallel for
                for(l=0; l<dim[SECOND_MATRIX][COLUMNS]; ++l)
                {
					idx[FIRST_MATRIX] = dim[FIRST_MATRIX][COLUMNS]*i + j;
					idx[SECOND_MATRIX] = dim[SECOND_MATRIX][COLUMNS]*k + l;
					idx[MATRIX_PRODUCT] = pdim*(k+(i*dim[SECOND_MATRIX][ROWS])) + l+(j*dim[SECOND_MATRIX][COLUMNS]);
                	#pragma omp parallel sections num_threads(MAX_OCTONIONS_UNITS)
					{
						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL); 
					        mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_set(matrix_product[OCTONIONS_REALPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                        mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_set(matrix_product[OCTONIONS_E1PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                        mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_set(matrix_product[OCTONIONS_E2PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                        mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_set(matrix_product[OCTONIONS_E3PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                	}
	
						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                        mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_set(matrix_product[OCTONIONS_E4PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                        mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_set(matrix_product[OCTONIONS_E5PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                	}
	
						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                        mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_set(matrix_product[OCTONIONS_E6PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                        mpfr_mul(tmp, matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]], matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_set(matrix_product[OCTONIONS_E7PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, NULL);
	                	}
					}
				}
    return;
}

__MATHSUITE void    _matrixKSProduct(mpfr_t **matrix1, mpfr_t **matrix2, mpfr_t **matrix_product, register dim_typ dim[static 2][2])
{

    dim_typ i, j, k, l;
    register dim_typ idx[MAX_MATRICES];
    const register dim_typ pdim = dim[FIRST_MATRIX][COLUMNS]*dim[SECOND_MATRIX][COLUMNS];

	#pragma omp parallel for
    for(i=0; i<dim[FIRST_MATRIX][ROWS]; ++i)
    	#pragma omp parallel for
        for(k=0; k<dim[SECOND_MATRIX][ROWS]; ++k)
    		#pragma omp parallel for
        	for(j=0; j<dim[FIRST_MATRIX][COLUMNS]; ++j)
            	#pragma omp parallel for
                for(l=0; l<dim[SECOND_MATRIX][COLUMNS]; ++l)
                {
                	idx[FIRST_MATRIX] = dim[FIRST_MATRIX][COLUMNS]*i + j;
					idx[SECOND_MATRIX] = dim[SECOND_MATRIX][COLUMNS]*k + l;
					idx[MATRIX_PRODUCT] = pdim*(k+(i*dim[SECOND_MATRIX][ROWS])) + l+(j*dim[SECOND_MATRIX][COLUMNS]);
                	#pragma omp parallel sections num_threads(MAX_SEDENIONS_UNITS)
					{
						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_REALPART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E1PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E2PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E3PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E4PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E5PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E6PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E7PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
                       	 	mpfr_set(matrix_product[SEDENIONS_E8PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
                       	 	mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E9PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E10PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E11PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E12PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E13PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E14PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}

						#pragma omp section
						{
							mpfr_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16;
							mpfr_inits(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL);
	                        mpfr_mul(tmp, matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp2, matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp3, matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp4, matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp5, matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp6, matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp7, matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp8, matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp9, matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp10, matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp11, matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp12, matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp13, matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp14, matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp15, matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_mul(tmp16, matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]], matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]], MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp3, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp4, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp5, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp6, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp7, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp8, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp9, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp10, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp11, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp12, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp13, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp14, MPFR_RNDN);
	                        mpfr_sub(tmp, tmp, tmp15, MPFR_RNDN);
	                        mpfr_add(tmp, tmp, tmp16, MPFR_RNDN);
	                        mpfr_set(matrix_product[SEDENIONS_E15PART][idx[MATRIX_PRODUCT]], tmp, MPFR_RNDN);
	                        mpfr_clears(tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, NULL); 
	                	}
					}
				}
    return;
}

// ENDAlgebra

// END
