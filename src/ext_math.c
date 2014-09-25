// ext_math.h
// as External Math Library
// 16/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be included
exclusively by main.c, geometry.c, programs.c, algebra.c and settings.c project files of my suite program!
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

__MSUTIL_ void __export toBinary(int,char *);
__MSUTIL_ int __export toDecimal(char *);
__MSUTIL_ void __export rComplement(char *);
__MSUTIL_ void __export binaryAdd(char *,char *,char *);

__MSNATIVE_ static ityp __system __export math_sum(register ityp, register ityp);
__MSNATIVE_ static ityp __system __export math_sub(register ityp, register ityp);
__MSNATIVE_ static double complex __system __export math_csum(register double complex, register double complex);
__MSNATIVE_ static double complex __system __export math_csub(register double complex, register double complex);

const struct ext_math_type ext_math =
{
    {
        {
            false,
            true
        },
        {
        	0,
        	UCHAR_MAX
        },
        {
            SCHAR_MIN,
            SCHAR_MAX
        },
        {
            INT_MIN,
            INT_MAX
        },
        {
            SHRT_MIN,
            SHRT_MAX
        },
        {
            LONG_MIN,
            LONG_MAX,
        },
        {
            LLONG_MIN,
            LLONG_MAX
        },
        {
            0,
            UINT_MAX
        },
        {
            0,
            USHRT_MAX
        },
        {
            0,
            ULONG_MAX
        },
        {
            0,
            ULLONG_MAX
        },
        {
            -FLT_MAX,
            FLT_MAX
        },
        {
            -DBL_MAX,
            DBL_MAX
        },
        {
            -LDBL_MAX,
            LDBL_MAX
        }
    },
    {
        IDENTIFIER_SINANDSINH,
        IDENTIFIER_COSANDCOSH,
        IDENTIFIER_SINANDSINH"h",
        IDENTIFIER_COSANDCOSH"h",
        IDENTIFIER_CSCANDCSCH,
        IDENTIFIER_CSCANDCSCH"h",
        IDENTIFIER_SECANDSECH,
        IDENTIFIER_SECANDSECH"h",
        IDENTIFIER_ASINANDASINH,
        IDENTIFIER_ASINANDASINH"h",
        IDENTIFIER_ACOSANDACOSH,
        IDENTIFIER_ACOSANDACOSH"h",
        IDENTIFIER_ACSCANDACSCH,
        IDENTIFIER_ACSCANDACSCH"h",
        IDENTIFIER_ASECANDASECH,
        IDENTIFIER_ASECANDASECH"h",
        IDENTIFIER_TANANDTANH,
        IDENTIFIER_TANANDTANH"h",
        IDENTIFIER_ATANANDATANH,
        IDENTIFIER_ATANANDATANH"h",
        IDENTIFIER_COTANDCOTH,
        IDENTIFIER_COTANDCOTH"h",
        IDENTIFIER_ACOTANDACOTH,
        IDENTIFIER_ACOTANDACOTH"h",
        IDENTIFIER_HSINANDHSINH,
        IDENTIFIER_HSINANDHSINH"h",
        IDENTIFIER_QSINANDQSINH,
        IDENTIFIER_QSINANDQSINH"h",
        IDENTIFIER_HCOSANDHCOSH,
        IDENTIFIER_HCOSANDHCOSH"h",
        IDENTIFIER_QCOSANDQCOSH,
        IDENTIFIER_QCOSANDQCOSH"h",
        IDENTIFIER_HSECANDHSECH,
        IDENTIFIER_HSECANDHSECH"h",
        IDENTIFIER_QSECANDQSECH,
        IDENTIFIER_QSECANDQSECH"h",
        IDENTIFIER_HCSCANDHCSCH,
        IDENTIFIER_HCSCANDHCSCH"h",
        IDENTIFIER_QCSCANDQCSCH,
        IDENTIFIER_QCSCANDQCSCH"h",
        IDENTIFIER_HTANANDHTANH,
        IDENTIFIER_HTANANDHTANH"h",
        IDENTIFIER_QTANANDQTANH,
        IDENTIFIER_QTANANDQTANH"h",
        IDENTIFIER_HCOTANDHCOTH,
        IDENTIFIER_HCOTANDHCOTH"h",
        IDENTIFIER_QCOTANDQCOTH,
        IDENTIFIER_QCOTANDQCOTH"h",
        IDENTIFIER_VSINANDVSINH,
        IDENTIFIER_VSINANDVSINH"h",
        IDENTIFIER_CVSINANDCVSINH,
        IDENTIFIER_CVSINANDCVSINH"h",
        IDENTIFIER_VCOSANDVCOSH,
        IDENTIFIER_VCOSANDVCOSH"h",
        IDENTIFIER_CVCOSANDCVCOSH,
        IDENTIFIER_CVCOSANDCVCOSH"h",
        IDENTIFIER_HVSINANDHVSINH,
        IDENTIFIER_HVSINANDHVSINH"h",
        IDENTIFIER_HCVSINANDHCVSINH,
        IDENTIFIER_HCVSINANDHCVSINH"h",
        IDENTIFIER_QVSINANDQVSINH,
        IDENTIFIER_QVSINANDQVSINH"h",
        IDENTIFIER_QCVSINANDQCVSINH,
        IDENTIFIER_QCVSINANDQCVSINH"h",
        IDENTIFIER_HVCOSANDHVCOSH,
        IDENTIFIER_HVCOSANDHVCOSH"h",
        IDENTIFIER_HCVCOSANDHCVCOSH,
        IDENTIFIER_HCVCOSANDHCVCOSH"h",
        IDENTIFIER_QVCOSANDQVCOSH,
        IDENTIFIER_QVCOSANDQVCOSH"h",
        IDENTIFIER_QCVCOSANDQCVCOSH,
        IDENTIFIER_QCVCOSANDQCVCOSH"h",
        IDENTIFIER_ESECANDESECH,
        IDENTIFIER_ESECANDESECH"h",
        IDENTIFIER_ECSCANDECSCH,
        IDENTIFIER_ECSCANDECSCH"h",
        IDENTIFIER_HESECANDHESECH,
        IDENTIFIER_HESECANDHESECH"h",
        IDENTIFIER_HECSCANDHECSCH,
        IDENTIFIER_HECSCANDHECSCH"h",
        IDENTIFIER_QESECANDQESECH,
        IDENTIFIER_QESECANDQESECH"h",
        IDENTIFIER_QECSCANDQECSCH,
        IDENTIFIER_QECSCANDQECSCH"h",
        IDENTIFIER_SINCANDSINCH,
        IDENTIFIER_SINCANDSINCH"h",
        IDENTIFIER_HSINCANDHSINCH,
        IDENTIFIER_HSINCANDHSINCH"h",
        IDENTIFIER_QSINCANDQSINCH,
        IDENTIFIER_QSINCANDQSINCH"h",
        IDENTIFIER_COSCANDCOSCH,
        IDENTIFIER_COSCANDCOSCH"h",
        IDENTIFIER_HCOSCANDHCOSCH,
        IDENTIFIER_HCOSCANDHCOSCH"h",
        IDENTIFIER_QCOSCANDQCOSCH,
        IDENTIFIER_QCOSCANDQCOSCH"h",
        IDENTIFIER_SECCANDSECCH,
        IDENTIFIER_SECCANDSECCH"h",
        IDENTIFIER_HSECCANDHSECCH,
        IDENTIFIER_HSECCANDHSECCH"h",
        IDENTIFIER_QSECANDQSECH,
        IDENTIFIER_QSECANDQSECH"h",
        IDENTIFIER_CSCCANDCSCCH,
        IDENTIFIER_CSCCANDCSCCH"h",
        IDENTIFIER_HCSCCANDHCSCCH,
        IDENTIFIER_HCSCCANDHCSCCH"h",
        IDENTIFIER_QCSCCANDQCSCCH,
        IDENTIFIER_QCSCCANDQCSCCH"h",
        IDENTIFIER_TANCANDTANCH,
        IDENTIFIER_TANCANDTANCH"h",
        IDENTIFIER_HTANCANDHTANCH,
        IDENTIFIER_HTANCANDHTANCH"h",
        IDENTIFIER_QTANCANDQTANCH,
        IDENTIFIER_QTANCANDQTANCH"h",
        IDENTIFIER_COTCANDCOTCH,
        IDENTIFIER_COTCANDCOTCH"h",
        IDENTIFIER_HCOTCANDHCOTCH,
        IDENTIFIER_HCOTCANDHCOTCH"h",
        IDENTIFIER_QCOTANDQCOTH,
        IDENTIFIER_QCOTANDQCOTH"h",
        IDENTIFIER_LOGARITMO,
        "log10",
        IDENTIFIER_LOGARITMO2,
        IDENTIFIER_LOGARITMOC,
        "log10c",
        IDENTIFIER_LOGARITMO2C,
        IDENTIFIER_LOGARITMO1P,
        IDENTIFIER_LOGARITMO1PC,
        IDENTIFIER_EXPANDEXPC,
        IDENTIFIER_EXPANDEXPC"c",
        IDENTIFIER_EXP10ANDEXP10C,
        IDENTIFIER_EXP10ANDEXP10C"c",
        IDENTIFIER_EXP2ANDEXP2C,
        IDENTIFIER_EXP2ANDEXP2C"c",
        IDENTIFIER_SOMMASUCCESSIONEARMONICA,
        IDENTIFIER_FIBONACCI,
        IDENTIFIER_FATTORIALE,
        IDENTIFIER_SEMIFATTORIALE,
        IDENTIFIER_NESIMONUMEROPRIMO,
        IDENTIFIER_PRIMORIALE,
        IDENTIFIER_SOMMAPRIMINNUMERIPRIMI,
        IDENTIFIER_FIBONACCIALE,
        IDENTIFIER_SOMMASUCCESSIONEFIBONACCI,
        IDENTIFIER_SOMMASUCCESSIONEFATTORIALE,
        IDENTIFIER_SOMMASUCCESSIONESEMIFATTORIALE,
        IDENTIFIER_SOMMAPRIMINNUMERI,
        "floor",
        "ceil",
        "deg",
        "rad",
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
        asum,
        fibo,
        fact,
        sfact,
        N_prime_Number,
        primr,
        fpnsum,
        fibnc,
        fsum,
        fasum,
        sfasum,
        fnnsum,
        floor,
        ceil,
        deg,
        rad,
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

__MSUTIL_ void __export toBinary(int a,char *c)
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

__MSUTIL_ int __export toDecimal(char *c)
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

__MSUTIL_ void __export rComplement(char *c)
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

__MSUTIL_ void __export binaryAdd(char *first,char *second,char *sum)
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

__MSNATIVE_ char * const __system __export binaryAlgSum(ityp a, ityp b, bool mode)
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
        printf2(COLOR_USER, "\t %s\t %d\n",c,num1);
        printf2(COLOR_USER, "\t+%s\t+%d =\n",d,num2);

        if(*(e)=='1')
        {
            printf2(COLOR_USER, "\t %s\t",e);
            rComplement(e);
            printf2(COLOR_USER, "-%d",toDecimal(e));
        }
        else
        {
            printf2(COLOR_USER, "\t %s\t ",e);
            printf2(COLOR_USER, "%d",toDecimal(e));
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
        printf2(COLOR_USER, "\t %s\t %d\n",c,a);
        if(b<0)
            printf2(COLOR_USER, "\t+%s\t+%d\n ",d,abs(b));
        else
            printf2(COLOR_USER, "\t+%s\t-%d\n ",d,b);

        if(*(e)=='1')
        {
            printf2(COLOR_USER, "\t %s\t",e);
            rComplement(e);
            printf2(COLOR_USER, "-%d",toDecimal(e));
        }
        else
        {
            printf2(COLOR_USER, "\t %s\t ",e);
            printf2(COLOR_USER, "%d",toDecimal(e));
        }
    }

    PRINTN();

    // free(c);
    free(d);
    free(e);

    return c;
}

__MSNATIVE_ char * const __system __export binNumComp(ityp a)
{
    int num;
    num = (int) a;
    char *c;

    c = (char *) malloc(SIZE+1);
    errMem(c, NULL);
    strcpy(c, NULL_CHAR);

    toBinary(a, c);

    char bin_num[sizeof(c)];
    strcpy(bin_num, c);

    rComplement(c);

    printf2(COLOR_USER, "\nInserted DECIMAL NUMBER: %lld;\nIts BINARY Conversion is: %s.\n\
Il suo COMPLEMENTO a DUE e': %s.\n", num, bin_num, c);

    toBinary(~num, c);

    printf2(COLOR_USER, "Decimal ONE'S COMPLEMENT: %d;\nBinary ONE'S COMPLEMENT is: %s.\n\n", ~num, c);

    // free(c);

    return c;
}

__MSNATIVE_ inline ityp __export ___cabs(ityp *restrict cpx, const register sel_typ dim)
{
	const register sel_typ algebra_units = exp2(dim);
	ityp accumulate = 0.00;
	for(sel_typ i=0; i<algebra_units; ++i)
		accumulate += exp2(*(cpx + i));
	return sqrt(accumulate);
}

__MSNATIVE_ inline void __export _complexAdd(ityp *restrict cpx, ityp complexRes[static MAX_COMPLEX_UNITS])
{
	const register double complex result = (*(cpx + REAL_PART) + (*(cpx + IMAG_PART))*I) + (*(cpx + MAX_COMPLEX_UNITS + REAL_PART) + (*(cpx + MAX_COMPLEX_UNITS + IMAG_PART))*I);
	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
	{
		complexRes[REAL_PART] = creal(result);
    	complexRes[IMAG_PART] = cimag(result);
	}
    return;
}

__MSNATIVE_ inline void __export _complexSub(ityp *restrict cpx, ityp complexRes[static MAX_COMPLEX_UNITS])
{
	const register double complex result = (*(cpx + REAL_PART) + (*(cpx + IMAG_PART))*I) - (*(cpx + MAX_COMPLEX_UNITS + REAL_PART) + (*(cpx + MAX_COMPLEX_UNITS + IMAG_PART))*I);
	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
	{
		complexRes[REAL_PART] = creal(result);
    	complexRes[IMAG_PART] = cimag(result);
	}
    return;
}

__MSNATIVE_ inline void __export _complexMul(ityp *restrict cpx, ityp complexRes[static MAX_COMPLEX_UNITS])
{
	const register double complex result = (*(cpx + REAL_PART) + (*(cpx + IMAG_PART))*I) * (*(cpx + MAX_COMPLEX_UNITS + REAL_PART) + (*(cpx + MAX_COMPLEX_UNITS + IMAG_PART))*I);
	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
	{
		complexRes[REAL_PART] = creal(result);
	    complexRes[IMAG_PART] = cimag(result);
	}
    return;
}

__MSNATIVE_ inline void __export _complexDiv(ityp *restrict cpx, ityp complexRes[static MAX_COMPLEX_UNITS])
{
	const register double complex result = (*(cpx + REAL_PART) + (*(cpx + IMAG_PART))*I) / (*(cpx + MAX_COMPLEX_UNITS + REAL_PART) + (*(cpx + MAX_COMPLEX_UNITS + IMAG_PART))*I);
	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
	{
		complexRes[REAL_PART] = creal(result);
	    complexRes[IMAG_PART] = cimag(result);
	}
    return;
}

__MSNATIVE_ inline void __export _quaternionsAdd(ityp *restrict quaternions, ityp quaternionsRes[static MAX_QUATERNIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
	{
		quaternionsRes[QUATERNIONS_REALPART] = *(quaternions + QUATERNIONS_REALPART) + *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART);
	    quaternionsRes[QUATERNIONS_IPART] = *(quaternions + QUATERNIONS_IPART) + *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART);
	    quaternionsRes[QUATERNIONS_JPART] = *(quaternions + QUATERNIONS_JPART) + *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART);
	    quaternionsRes[QUATERNIONS_KPART] = *(quaternions + QUATERNIONS_KPART) + *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART);
	}
    return;
}


__MSNATIVE_ inline void __export _quaternionsSub(ityp *restrict quaternions, ityp quaternionsRes[static MAX_QUATERNIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
	{
		quaternionsRes[QUATERNIONS_REALPART] = *(quaternions + QUATERNIONS_REALPART) - *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART);
	    quaternionsRes[QUATERNIONS_IPART] = *(quaternions + QUATERNIONS_IPART) - *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART);
	    quaternionsRes[QUATERNIONS_JPART] = *(quaternions + QUATERNIONS_JPART) - *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART);
	    quaternionsRes[QUATERNIONS_KPART] = *(quaternions + QUATERNIONS_KPART) - *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART);
	}
    return;
}

__MSNATIVE_ void __export _quaternionsMul(ityp *restrict quaternions, ityp quaternionsRes[static MAX_QUATERNIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
	{
		
		quaternionsRes[QUATERNIONS_REALPART] = *(quaternions + QUATERNIONS_REALPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART) -
											   *(quaternions + QUATERNIONS_IPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART) -
											   *(quaternions + QUATERNIONS_JPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART) -
											   *(quaternions + QUATERNIONS_KPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART);
									
		quaternionsRes[QUATERNIONS_IPART] = *(quaternions + QUATERNIONS_REALPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART) +
											*(quaternions + QUATERNIONS_IPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART) +
											*(quaternions + QUATERNIONS_JPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART) -	
											*(quaternions + QUATERNIONS_KPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART);
											
											
		quaternionsRes[QUATERNIONS_JPART] = *(quaternions + QUATERNIONS_REALPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART) +
											*(quaternions + QUATERNIONS_JPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART) +
											*(quaternions + QUATERNIONS_KPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART) -
											*(quaternions + QUATERNIONS_IPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART);
											
		quaternionsRes[QUATERNIONS_KPART] = *(quaternions + QUATERNIONS_REALPART)* *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART) +
											*(quaternions + QUATERNIONS_KPART)* *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART) +
											*(quaternions + QUATERNIONS_IPART)* *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART) -
											*(quaternions + QUATERNIONS_JPART)* *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART);
	}

    return;
}

__MSNATIVE_ void __export _quaternionsDiv(ityp *restrict quaternions, ityp quaternionsRes[static MAX_QUATERNIONS_UNITS])
{
	const register ityp squared_qnorm = exp2(_qabs(quaternions + MAX_QUATERNIONS_UNITS));
	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
	{
		quaternionsRes[QUATERNIONS_REALPART] = (*(quaternions + QUATERNIONS_REALPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART) +
											   *(quaternions + QUATERNIONS_IPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART) +
											   *(quaternions + QUATERNIONS_JPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART) +
											   *(quaternions + QUATERNIONS_KPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART)) / squared_qnorm;
									
		quaternionsRes[QUATERNIONS_IPART] = (*(quaternions + QUATERNIONS_REALPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART) -
											*(quaternions + QUATERNIONS_IPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART) -
											*(quaternions + QUATERNIONS_JPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART) +	
											*(quaternions + QUATERNIONS_KPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART)) / squared_qnorm;
											
											
		quaternionsRes[QUATERNIONS_JPART] = (*(quaternions + QUATERNIONS_REALPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART) -
											*(quaternions + QUATERNIONS_JPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART) -
											*(quaternions + QUATERNIONS_KPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART) +
											*(quaternions + QUATERNIONS_IPART) * *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART)) / squared_qnorm;
											
		quaternionsRes[QUATERNIONS_KPART] = (*(quaternions + QUATERNIONS_REALPART)* *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_KPART) -
											*(quaternions + QUATERNIONS_KPART)* *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_REALPART) -
											*(quaternions + QUATERNIONS_IPART)* *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_JPART) +
											*(quaternions + QUATERNIONS_JPART)* *(quaternions + MAX_QUATERNIONS_UNITS + QUATERNIONS_IPART)) / squared_qnorm;
	}

    return;
}

__MSNATIVE_ void __export _octonionsAdd(ityp *restrict octonions, ityp octonionsRes[static MAX_OCTONIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
	{
		octonionsRes[OCTONIONS_REALPART] = *(octonions + OCTONIONS_REALPART) + *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART);
	    octonionsRes[OCTONIONS_E1PART] = *(octonions + OCTONIONS_E1PART) + *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART);
	    octonionsRes[OCTONIONS_E2PART] = *(octonions + OCTONIONS_E2PART) + *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART);
	    octonionsRes[OCTONIONS_E3PART] = *(octonions + OCTONIONS_E3PART) + *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART);
	    octonionsRes[OCTONIONS_E4PART] = *(octonions + OCTONIONS_E4PART) + *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART);
	    octonionsRes[OCTONIONS_E5PART] = *(octonions + OCTONIONS_E5PART) + *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART);
	    octonionsRes[OCTONIONS_E6PART] = *(octonions + OCTONIONS_E6PART) + *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART);
	    octonionsRes[OCTONIONS_E7PART] = *(octonions + OCTONIONS_E7PART) + *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART);
	}
    return;
}

__MSNATIVE_ void __export _octonionsSub(ityp *restrict octonions, ityp octonionsRes[static MAX_OCTONIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
	{
		octonionsRes[OCTONIONS_REALPART] = *(octonions + OCTONIONS_REALPART) - *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART);
	    octonionsRes[OCTONIONS_E1PART] = *(octonions + OCTONIONS_E1PART) - *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART);
	    octonionsRes[OCTONIONS_E2PART] = *(octonions + OCTONIONS_E2PART) - *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART);
	    octonionsRes[OCTONIONS_E3PART] = *(octonions + OCTONIONS_E3PART) - *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART);
	    octonionsRes[OCTONIONS_E4PART] = *(octonions + OCTONIONS_E4PART) - *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART);
	    octonionsRes[OCTONIONS_E5PART] = *(octonions + OCTONIONS_E5PART) - *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART);
	    octonionsRes[OCTONIONS_E6PART] = *(octonions + OCTONIONS_E6PART) - *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART);
	    octonionsRes[OCTONIONS_E7PART] = *(octonions + OCTONIONS_E7PART) - *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART);
	}
    return;
}

__MSNATIVE_ void __export _octonionsMul(ityp *restrict octonions, ityp octonionsRes[static MAX_OCTONIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
	{
		octonionsRes[OCTONIONS_REALPART] = *(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) -
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) -
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) -
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) -
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) -
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) -
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) -
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART);
										   
		octonionsRes[OCTONIONS_E1PART] = *(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) +
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) +
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) +
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) -
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) +
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) -
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) -
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART);
										   
		octonionsRes[OCTONIONS_E2PART] = *(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) -
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) +
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) +
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) -
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) -
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) +
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) -
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART);
										   
		octonionsRes[OCTONIONS_E3PART] = *(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) -
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) -
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) +
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) +
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) +
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) -
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) +
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART);								   

	    octonionsRes[OCTONIONS_E4PART] = *(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) +
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) -
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) -
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) +
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) +
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) +
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) -
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART);								   


	    octonionsRes[OCTONIONS_E5PART] = *(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) -
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) +
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) -
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) -
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) +
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) +
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) +
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART);								   


	    octonionsRes[OCTONIONS_E6PART] = *(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) +
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) -
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) +
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) -
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) -
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) +
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) +
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART);								   


	    octonionsRes[OCTONIONS_E7PART] = *(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) +
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) +
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) -
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) +
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) -
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) -
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) +
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART);								   
	}
    return;
}

__MSNATIVE_ void __export _octonionsDiv(ityp *restrict octonions, ityp octonionsRes[static MAX_OCTONIONS_UNITS])
{
		const register ityp squared_onorm = exp2(_oabs(octonions + MAX_OCTONIONS_UNITS));
	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
	{
		octonionsRes[OCTONIONS_REALPART] = (*(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) +
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) +
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) +
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) +
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) +
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) +
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) +
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART)) / squared_onorm;
										   
		octonionsRes[OCTONIONS_E1PART] = (*(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) -
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) -
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) -
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) +
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) -
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) +
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) +
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART)) / squared_onorm;
										   
		octonionsRes[OCTONIONS_E2PART] = (*(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) +
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) -
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) -
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) +
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) +
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) -
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) +
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART)) / squared_onorm;
										   
		octonionsRes[OCTONIONS_E3PART] = (*(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) +
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) +
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) -
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) -
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) -
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) +
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) -
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART)) / squared_onorm;								   

	    octonionsRes[OCTONIONS_E4PART] = (*(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) -
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) +
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) +
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) -
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) -
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) -
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) +
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART)) / squared_onorm;							   


	    octonionsRes[OCTONIONS_E5PART] = (*(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) +
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) -
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) +
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) +
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) -
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) -
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) -
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART)) / squared_onorm;							   


	    octonionsRes[OCTONIONS_E6PART] = (*(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) -
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) +
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) -
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) +
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) +
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) -
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART) -
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART)) / squared_onorm;								   


	    octonionsRes[OCTONIONS_E7PART] = (*(octonions + OCTONIONS_REALPART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E7PART) -
										   *(octonions + OCTONIONS_E1PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E3PART) -
										   *(octonions + OCTONIONS_E2PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E6PART) +
										   *(octonions + OCTONIONS_E3PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E1PART) -
										   *(octonions + OCTONIONS_E4PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E5PART) +
										   *(octonions + OCTONIONS_E5PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E4PART) +
										   *(octonions + OCTONIONS_E6PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_E2PART) -
										   *(octonions + OCTONIONS_E7PART) * *(octonions + MAX_OCTONIONS_UNITS + OCTONIONS_REALPART)) / squared_onorm;
	}
    return;
}

__MSNATIVE_ void __export _sedenionsAdd(ityp *restrict sedenions, ityp sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
	{
		sedenionsRes[SEDENIONS_REALPART] = *(sedenions + SEDENIONS_REALPART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART);
	    sedenionsRes[SEDENIONS_E1PART] = *(sedenions + SEDENIONS_E1PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART);
	    sedenionsRes[SEDENIONS_E2PART] = *(sedenions + SEDENIONS_E2PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART);
	    sedenionsRes[SEDENIONS_E3PART] = *(sedenions + SEDENIONS_E3PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART);
	    sedenionsRes[SEDENIONS_E4PART] = *(sedenions + SEDENIONS_E4PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART);
	    sedenionsRes[SEDENIONS_E5PART] = *(sedenions + SEDENIONS_E5PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART);
	    sedenionsRes[SEDENIONS_E6PART] = *(sedenions + SEDENIONS_E6PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART);
	    sedenionsRes[SEDENIONS_E7PART] = *(sedenions + SEDENIONS_E7PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART);
	    sedenionsRes[SEDENIONS_E8PART] = *(sedenions + SEDENIONS_E8PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART);
	    sedenionsRes[SEDENIONS_E9PART] = *(sedenions + SEDENIONS_E9PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART);
	    sedenionsRes[SEDENIONS_E10PART] = *(sedenions + SEDENIONS_E10PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART);
	    sedenionsRes[SEDENIONS_E11PART] = *(sedenions + SEDENIONS_E11PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART);
	    sedenionsRes[SEDENIONS_E12PART] = *(sedenions + SEDENIONS_E12PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART);
	    sedenionsRes[SEDENIONS_E13PART] = *(sedenions + SEDENIONS_E13PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART);
	    sedenionsRes[SEDENIONS_E14PART] = *(sedenions + SEDENIONS_E14PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART);
	    sedenionsRes[SEDENIONS_E15PART] = *(sedenions + SEDENIONS_E15PART) + *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART);
	}
    return;
}

__MSNATIVE_ void __export _sedenionsSub(ityp *restrict sedenions, ityp sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
	{
		sedenionsRes[SEDENIONS_REALPART] = *(sedenions + SEDENIONS_REALPART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART);
	    sedenionsRes[SEDENIONS_E1PART] = *(sedenions + SEDENIONS_E1PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART);
	    sedenionsRes[SEDENIONS_E2PART] = *(sedenions + SEDENIONS_E2PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART);
	    sedenionsRes[SEDENIONS_E3PART] = *(sedenions + SEDENIONS_E3PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART);
	    sedenionsRes[SEDENIONS_E4PART] = *(sedenions + SEDENIONS_E4PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART);
	    sedenionsRes[SEDENIONS_E5PART] = *(sedenions + SEDENIONS_E5PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART);
	    sedenionsRes[SEDENIONS_E6PART] = *(sedenions + SEDENIONS_E6PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART);
	    sedenionsRes[SEDENIONS_E7PART] = *(sedenions + SEDENIONS_E7PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART);
	    sedenionsRes[SEDENIONS_E8PART] = *(sedenions + SEDENIONS_E8PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART);
	    sedenionsRes[SEDENIONS_E9PART] = *(sedenions + SEDENIONS_E9PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART);
	    sedenionsRes[SEDENIONS_E10PART] = *(sedenions + SEDENIONS_E10PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART);
	    sedenionsRes[SEDENIONS_E11PART] = *(sedenions + SEDENIONS_E11PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART);
	    sedenionsRes[SEDENIONS_E12PART] = *(sedenions + SEDENIONS_E12PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART);
	    sedenionsRes[SEDENIONS_E13PART] = *(sedenions + SEDENIONS_E13PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART);
	    sedenionsRes[SEDENIONS_E14PART] = *(sedenions + SEDENIONS_E14PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART);
	    sedenionsRes[SEDENIONS_E15PART] = *(sedenions + SEDENIONS_E15PART) - *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART);
	}
    return;
}

__MSNATIVE_ void __export _sedenionsMul(ityp *restrict sedenions, ityp sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
	{
		sedenionsRes[SEDENIONS_REALPART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART);
										   
		sedenionsRes[SEDENIONS_E1PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART);							

	    sedenionsRes[SEDENIONS_E2PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART);	

	    sedenionsRes[SEDENIONS_E3PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART);	

	    sedenionsRes[SEDENIONS_E4PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART);	

	    sedenionsRes[SEDENIONS_E5PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART);	

	    sedenionsRes[SEDENIONS_E6PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART);	

	    sedenionsRes[SEDENIONS_E7PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART);
										   
		sedenionsRes[SEDENIONS_E8PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART);	

	    sedenionsRes[SEDENIONS_E9PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART);	

	    sedenionsRes[SEDENIONS_E10PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART);	

	    sedenionsRes[SEDENIONS_E11PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART);	

	    sedenionsRes[SEDENIONS_E12PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART);	

	    sedenionsRes[SEDENIONS_E13PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART);	

	    sedenionsRes[SEDENIONS_E14PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART);	

	    sedenionsRes[SEDENIONS_E15PART] = *(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART);	
	}
    return;
}

__MSNATIVE_ void __export _sedenionsDiv(ityp *restrict sedenions, ityp sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	const register ityp squared_snorm = exp2(_sabs(sedenions + MAX_SEDENIONS_UNITS));
	#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
	{
		sedenionsRes[SEDENIONS_REALPART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART)) / squared_snorm;
										   
		sedenionsRes[SEDENIONS_E1PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART)) / squared_snorm;							

	    sedenionsRes[SEDENIONS_E2PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E3PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E4PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E5PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E6PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E7PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART)) / squared_snorm;
										   
		sedenionsRes[SEDENIONS_E8PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E9PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E10PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E11PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E12PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E13PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) +
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) +
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E14PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) +
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E7PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART)) / squared_snorm;	

	    sedenionsRes[SEDENIONS_E15PART] = (*(sedenions + SEDENIONS_REALPART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E15PART) -
										   *(sedenions + SEDENIONS_E1PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E14PART) +
										   *(sedenions + SEDENIONS_E2PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E13PART) +
										   *(sedenions + SEDENIONS_E3PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E12PART) -
										   *(sedenions + SEDENIONS_E4PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E11PART) -
										   *(sedenions + SEDENIONS_E5PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E10PART) -
										   *(sedenions + SEDENIONS_E6PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E9PART) -
										   *(sedenions + SEDENIONS_E7PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E8PART) +
										   *(sedenions + SEDENIONS_E8PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E9PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E6PART) +
										   *(sedenions + SEDENIONS_E10PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E5PART) +
										   *(sedenions + SEDENIONS_E11PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E4PART) -
										   *(sedenions + SEDENIONS_E12PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E3PART) -
										   *(sedenions + SEDENIONS_E13PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E2PART) +
										   *(sedenions + SEDENIONS_E14PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_E1PART) -
										   *(sedenions + SEDENIONS_E15PART) * *(sedenions + MAX_SEDENIONS_UNITS + SEDENIONS_REALPART)) / squared_snorm;
	}
    return;
}

__MSNATIVE_ bool __export _secondGradeEquationSolver(ityp *restrict abc, ityp root[static MAX_DIMENSIONS])
{
    ityp dscr; // storage for equation's discriminant

    dscr = (abc[COEFF_B]*abc[COEFF_B]) - 4*abc[COEFF_A]*abc[COEFF_C];

    if(dscr < 0)
    {
        printf2(COLOR_SYSTEM, "\nThese Equation has Imaginary ROOTS.\n\n");
        return false;
    }

    root[ROOT_X1] = (-abc[COEFF_B] + sqrt(dscr))/(2.00*abc[COEFF_A]);
    root[ROOT_X2] = (-abc[COEFF_B] - sqrt(dscr))/(2.00*abc[COEFF_A]);

    return true;
}

/// thanks to Bibek Subedi for this function,
/// which I renamed, modified and adapted to this program. Link at:
/// http://programming-technique.blogspot.it/2011/08/conversion-of-decimal-to-roman-using-c.html
__MSUTIL_ inline void __export getRomanNumber(dim_typ number, char string[static ROMAN_NUMBER_STRING])
{
    dim_typ i, j, k;
    // k = 0;

    for(i=1000,k=0; i >= 1; i /= 10,++k)
        if((j=((number/i)%(k ? 10 : 1))))
            strcat(string, k ? ext_math.romn_map[k-1][j-1] : ext_math.romn_thousand_map[j-1]);

    return;
}

__MSUTIL_ ityp __export mpow(ityp b, int64_t exp)
{
   if (exp < 0)
   {
      b = 1/b;
      exp *= -1;
   }

   return mpow2(b,exp);
}

__MSUTIL_ ityp __export mpow2(ityp b, int64_t exp)
{
   if (exp > 2)
       return exp%2 ? b*mpow(b*b,(int64_t)exp*0.5) : mpow(b*b,exp*0.5);
   else if (2 == exp)
      return b*b;
   else if (1 == exp)
      return b;

   return 1.00; // exp == 0
}

/// special thanks to Bibek Subedi for this function,
/// which I renamed, adapted and modified to this program. Link at:
/// http://programming-technique.blogspot.it/2013/05/pascals-triangle-using-c-program.html
__MSNATIVE_ __MSUTIL_ inline void __export getPascalTriangle(uint64_t rows)
{
    uint64_t i, j;
    PRINTL();

    for(i = 0; i < rows; ++i)
    {
        for (j = 0; j < rows - i; ++j)
            PRINTSPACE();
        for(j = 0; ++j <= i; )
            printf2(COLOR_USER, "%d ", comb(i, j));
        PRINTN();
    }

    PRINTL();
    return;
}

/// PRIME NUMBERS FUNCTIONS
/// thanks to Bibek Subedi for these functions,
/// that I renamed, modified and adapted to this program both each. Link at:
/// http://programming-technique.blogspot.it/2012/02/fastest-way-of-calculating-prime-number.html

/// HISP isPrime bool checking function
/// (High Iteration - Slow Process)
__MSUTIL_ inline bool __system __export isPrimeHISP(register uint64_t n)
{
    uint64_t i = 2;
    uint64_t result;
    do
    {
        if(!(result = n%i)) break;
        ++ i;
    }
    while(n != i);
    if(n == i) return true;
    return false;
}

/// LIFP isPrime bool checking function
/// (Low Iterations, Fast process)
__MSUTIL_ inline bool __system __export isPrimeLIFP(register uint64_t n)
{
    uint64_t i;
    // this because in my program The user won't be able to insert
    // a n minor of 3.
    // if(n == 2) return false; // true;
    if(!(n%2)) return false;
    for(i = 2; i*i<=n; i+=2)
        if(!(n%i)) return false;
    return true;
}

__MSNATIVE_ void __system __export prime_N_Number(register uint64_t m, register uint64_t n)
{
    uint64_t i;
    uint64_t dim;

    dim = 0;
    // tab = malloc(sizeof)
    SHOWPAUSEMESSAGE();
    PRINTL();

    bool (* const prime_check_func)(register uint64_t) = lazy_exec ? isPrimeHISP : isPrimeLIFP;
    //
    PRINTL();

    if(isnSett(lazy_exec) && m==2)
    {
        printf2(COLOR_USER, "- 1: 2;\n", ++dim, i);
        ++ m;
    }

    for(i=m; i <= n; ++i)
        if(prime_check_func(i) || m==2)
        {
            printf2(COLOR_USER, "- %llu: %llu;\n", ++dim, i);
            if(catchPause()) return;
        }

    PRINTL();
    printf2(COLOR_USER, "\n%llu PRIME NUMBERS has been correctly printed.\n\n", dim);

    return;
}

// N_prime_Number
__MSNATIVE_ inline ityp __system __export N_prime_Number(register ityp n)
{
    uint64_t i, dim;
    uint64_t isn;
    bool (*prime_check_func)(register uint64_t) = lazy_exec ? isPrimeHISP : isPrimeLIFP;

    isn = (uint64_t) n;
    dim = 0;

    if(isnSett(lazy_exec))
    {
        if(n == 1)
            return 2.00;
        -- isn;
    }

    if(n>0)
        for(i=2; ; ++i)
            if(prime_check_func(i) && ++dim == isn) return i;

    return INVALIDRETURNVALUE_NPRIMENUMBER; // basic returning 0 in case of error (obviously not a pn)
}

// Fibonaccial
__MSNATIVE_ inline ityp __export fibnc(register ityp n)
{
    ityp i, c;
        
	for(c=i=0.00; i<n; ++i)
		c *= fibo(i);

    return c;
}

// N_primorial_Number
__MSNATIVE_ inline ityp __export primr(register ityp n)
{
    ityp c = 0.00;

    for(uint64_t i=0; i<n; ++i)
        c *= N_prime_Number(i);

    return c;
}

// First N Prime Number Sum
__MSNATIVE_ inline ityp __export fpnsum(register ityp n)
{
    ityp c = 0.00;

    for(uint64_t i=0; i<n; ++i)
        c += N_prime_Number(i);

    return c;
}

 //compares if the float f1 is equal with f2 and returns 1 if true and 0 if false
 static inline bool compare_double(register ityp f1, register ityp f2)
 {
	return (((f1 - COMPAREDOUBLE_PRECISION) < f2) && ((f1 + COMPAREDOUBLE_PRECISION) > f2));
 }

__MSNATIVE_ inline bool __system __export trigonometric_domain(register ityp x, register ityp y, register ityp z)
{
    for(uint64_t i=0; i<=x; ++i)
    	if(compare_double(x,y+(i*z)))
    		return false;
    return true;
}

__MSNATIVE_ inline ityp __export stirling(register ityp n)
{
	return ((uint64_t)(sqrt(M_PI*(((uint64_t)n)<<1))*pow((n/M_E),n)));
}

__MSNATIVE_ ityp __export fibo(register ityp n)
{
    if(n <= 0.00) return 0.00;

    const uint64_t num = (uint64_t) n;

    if(!access(sysMem)[FUNCTION_FIBONACCI].current_max_index)
    {
        access(sysMem)[FUNCTION_FIBONACCI].memoizer = malloc(sizeof(ityp)<<1);
        errMem(access(sysMem)[FUNCTION_FIBONACCI].memoizer, INVALIDRETURNVALUE_FIBONACCI);
        access(sysMem)[FUNCTION_FIBONACCI].memoizer[0] = 0;
        access(sysMem)[FUNCTION_FIBONACCI].memoizer[1] = 1;
    }

    uint64_t prev, curr, next, i;

    prev = next = 0;
    curr = 1;

    if(num < access(sysMem)[FUNCTION_FIBONACCI].current_max_index)
        curr = access(sysMem)[FUNCTION_FIBONACCI].memoizer[num];
    else
    {
        if(access(sysMem)[FUNCTION_FIBONACCI].current_max_index < access(curLayout)->max_memoizable_indices[FUNCTION_FIBONACCI])
        {

            if(access(sysMem)[FUNCTION_FIBONACCI].current_max_index)
            {
                prev = access(sysMem)[FUNCTION_FIBONACCI].memoizer[access(sysMem)[FUNCTION_FIBONACCI].current_max_index -1];
                curr = access(sysMem)[FUNCTION_FIBONACCI].memoizer[access(sysMem)[FUNCTION_FIBONACCI].current_max_index];
            }

            for (i=(access(sysMem)[FUNCTION_FIBONACCI].current_max_index ? access(sysMem)[FUNCTION_FIBONACCI].current_max_index+1 : 2); i <= num; ++i)
            {
                next = curr+prev;
                prev = curr;
                access(sysMem)[FUNCTION_FIBONACCI].memoizer = realloc(access(sysMem)[FUNCTION_FIBONACCI].memoizer, sizeof(ityp)*(i+1));
                errMem(access(sysMem)[FUNCTION_FIBONACCI].memoizer, -1.00);
                access(sysMem)[FUNCTION_FIBONACCI].memoizer[i] = (curr = next);
            }
            access(sysMem)[FUNCTION_FIBONACCI].current_max_index = num;
        }
        else
            for (i = 1; ++i <= num; )
            {
                next = curr+prev;
                prev = curr;
                curr = next;
            }
    }

    return curr < ULLONG_MAX ? curr : INVALIDRETURNVALUE_FIBONACCI;
}


__MSNATIVE_ ityp __export fact(register ityp n)
{
    if(n < 0.00) return 0.00;

    const uint64_t num = (uint64_t) n;
    uint64_t i, res;

    if(!access(sysMem)[FUNCTION_FATTORIALE].current_max_index)
    {
        access(sysMem)[FUNCTION_FATTORIALE].memoizer = malloc(sizeof(ityp)<<1);
        errMem(access(sysMem)[FUNCTION_FATTORIALE].memoizer, INVALIDRETURNVALUE_FATTORIALE);
        access(sysMem)[FUNCTION_FATTORIALE].memoizer[0] = access(sysMem)[FUNCTION_FATTORIALE].memoizer[1] = 1;
    }

    res = 1;

    if(num < access(sysMem)[FUNCTION_FATTORIALE].current_max_index)
        res = access(sysMem)[FUNCTION_FATTORIALE].memoizer[num];
    else
    {
        // this cycle redundance has been introduced in order
        // to avoid overhead problems.
        if(access(sysMem)[FUNCTION_FATTORIALE].current_max_index < access(curLayout)->max_memoizable_indices[FUNCTION_FATTORIALE])
        {
            if(access(sysMem)[FUNCTION_FATTORIALE].current_max_index)
                res = access(sysMem)[FUNCTION_FATTORIALE].memoizer[access(sysMem)[FUNCTION_FATTORIALE].current_max_index];

            for(i=(access(sysMem)[FUNCTION_FATTORIALE].current_max_index ? access(sysMem)[FUNCTION_FATTORIALE].current_max_index+1 : 2); i<=num; ++i)
            {
                access(sysMem)[FUNCTION_FATTORIALE].memoizer = realloc(access(sysMem)[FUNCTION_FATTORIALE].memoizer, sizeof(ityp)*(i+1));
                errMem(access(sysMem)[FUNCTION_FATTORIALE].memoizer, res);
            	access(sysMem)[FUNCTION_FATTORIALE].memoizer[i] = i<access(curLayout)->min_stirling_number ?  (res *= i) : (res=stirling(i));
            }

            access(sysMem)[FUNCTION_FATTORIALE].current_max_index = num;
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

    return res < ULLONG_MAX ? res : INVALIDRETURNVALUE_FATTORIALE;
}

__MSNATIVE_ inline ityp __export sfact(register ityp n)
{
	return (((uint64_t)n)%2) ? sfact_odd(n) : sfact_even(n);
}

__MSNATIVE_ ityp __export sfact_even(register ityp n)
{
    if(n < 0.00) return 0.00;

    const uint64_t num = (uint64_t) n;
    uint64_t i, res;

    if(!access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index)
    {
        access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer = malloc(sizeof(ityp)<<1);
        errMem(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer, INVALIDRETURNVALUE_EVEN_SFATTORIALE);
        access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer[0] = access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer[1] = 1;
    }

    res = 1;

    if(num < access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index)
        res = access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer[num];
    else
    {
        // this cycle redundance has been introduced in order
        // to avoid overhead problems.
        if(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index < access(curLayout)->max_memoizable_indices[FUNCTION_EVEN_SFATTORIALE])
        {
            if(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index)
                res = access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer[access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index];

            for(i=(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index ? access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index+2 : 2); i<=num; i+=2)
            {
                access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer = realloc(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer, sizeof(ityp)*(i+1));
                errMem(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer, res);
                access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer[i] = (access(curLayout)->min_stirling_number%2 || num < (access(curLayout)->min_stirling_number<<1)) ? (res *= i) : (res = exp2(num)*stirling(num));
            }

            access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index = num;
        }
        else
    	{
    		if(access(curLayout)->min_stirling_number%2 || (num < access(curLayout)->min_stirling_number<<1))
	            for(i=2; i<=num; i+=2)
	                res *= i;
	        else
	    		res = exp2(num)*stirling(num);
	    }
    }

    return res < ULLONG_MAX ? res : INVALIDRETURNVALUE_EVEN_SFATTORIALE;
}

__MSNATIVE_ ityp __export sfact_odd(register ityp n)
{
    if(n < 0.00) return 0.00;

    const uint64_t num = (uint64_t) n;
    uint64_t i, res;

    if(!access(sysMem)[FUNCTION_ODD_SFATTORIALE].current_max_index)
    {
        access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer = malloc(sizeof(ityp)*MAX_ABSTRACT_DIMENSIONS);
        errMem(access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer, INVALIDRETURNVALUE_ODD_SFATTORIALE);
        access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer[0] = access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer[1] = 1;
        access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer[2] = 2;
    }

    res = 1;

    if(num < access(sysMem)[FUNCTION_ODD_SFATTORIALE].current_max_index)
        res = access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer[num];
    else
    {
        // this cycle redundance has been introduced in order
        // to avoid overhead problems.
        if(access(sysMem)[FUNCTION_ODD_SFATTORIALE].current_max_index < access(curLayout)->max_memoizable_indices[FUNCTION_ODD_SFATTORIALE])
        {
            if(access(sysMem)[FUNCTION_ODD_SFATTORIALE].current_max_index)
                res = access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer[access(sysMem)[FUNCTION_ODD_SFATTORIALE].current_max_index];

            for(i=(access(sysMem)[FUNCTION_ODD_SFATTORIALE].current_max_index ? access(sysMem)[FUNCTION_ODD_SFATTORIALE].current_max_index+2 : 3); i<=num; i+=2)
            {
                access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer = realloc(access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer, sizeof(ityp)*(i+1));
                errMem(access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer, res);
                access(sysMem)[FUNCTION_ODD_SFATTORIALE].memoizer[i] = (res *= i);
            }

            access(sysMem)[FUNCTION_ODD_SFATTORIALE].current_max_index = num;
        }
        else
            for(i=1; i<=num; i+=2)
                res *= i;
    }

    return res < ULLONG_MAX ? res : INVALIDRETURNVALUE_ODD_SFATTORIALE;
}

__MSUTIL_ inline ityp __export perm(register ityp n)
{
	return fact(n);
}

__MSUTIL_ inline ityp __export perm_rep(register ityp n, register uint64_t dim, ityp vector[static dim])
{
	 return (fact(n)/productory(dim, PRODUCTORY_MUL, vector));
}

__MSUTIL_ inline ityp __export kperm(register ityp n, register ityp k)
{
    return (fact(n)/(fact(n-k)));
}

__MSUTIL_ inline ityp __export kperm_rep(register ityp n, register ityp k)
{
	return pow(n, k);
}

__MSUTIL_ inline ityp __export comb(register ityp n, register ityp r)
{
    return (fact(n)/(fact(n-r)*fact(r)));
}

__MSUTIL_ inline ityp __export comb_rep(register ityp n, register ityp r)
{
	return (fact(n+r-1))/(fact(n-1)); 
}

// Somma SUCCESSIONE GEOMETRICA
__MSNATIVE_ inline ityp __export gsum(register ityp a, register int64_t exponent)
{
    return (1-(mpow2(a,exponent+1))/(1-a));
}

// Somma SUCCESSIONE ARMONICA GENERALIZZATA
__MSNATIVE_ inline ityp __export gasum(register uint64_t a, register uint64_t exponent)
{
    ityp c = 0.00;

	#pragma omp parallel for
    for(uint64_t i=0; i<a; ++i)
        c += (1/(mpow2(i+1, exponent)));

    return c;
}

// Somma SUCCESSIONE ARMONICA
__MSNATIVE_ inline ityp __export asum(register ityp a)
{
    return gasum((uint64_t) a, 1);
}

// Somma SUCCESSIONE FIBONACCI
__MSNATIVE_ inline ityp __export fsum(register ityp a)
{
    ityp c = 0.00;
    const uint64_t num = (uint64_t) a;

	#pragma omp parallel for
    for(uint64_t i=0; i<num; ++i)
        c += fibo(i);

    return c;
}

// Somma SUCCESSIONE FATTORIALE
__MSNATIVE_ inline ityp __export fasum(register ityp a)
{
    ityp c = 0.00;
    const uint64_t num = (uint64_t) a;

	#pragma omp parallel for
    for(uint64_t i=0; i<num; ++i)
        c += fact(i);

    return c;
}

// Somma SUCCESSIONE SEMIFATTORIALE
__MSNATIVE_ inline ityp __export sfasum(register ityp a)
{
    ityp c = 0.00;
    const uint64_t num = (uint64_t) a;

	#pragma omp parallel for
    for(uint64_t i=0; i<num; ++i)
        c += sfact(i);

    return c;
}

// Somma PRIMI N NUMERI NATURALI
__MSNATIVE_ inline ityp __export fnnsum(register ityp a)
{
    const uint64_t num = (uint64_t) a;
    return ((num*(num+1))*0.5);
}


// SOMMATORIA
__MSNATIVE_ inline ityp __export summation(uint64_t dim, bool mode, ityp vector[static dim])
{
    ityp res = 0.00;
    const register sel_typ mode_binder = 1-(mode<<1);

	#pragma omp parallel for
    for(uint64_t i=0; i<dim; ++i)
        res += (vector[i])*mode_binder;

    return res;
}

// PRODUTTORIA
__MSNATIVE_ inline ityp __export productory(uint64_t dim, bool mode, ityp vector[static dim])
{
    ityp res = 1.00;
    const register sel_typ mode_binder = 1-(mode<<1);

	#pragma omp parallel for
    for(uint64_t i=0; i<dim; ++i)
        res *= pow(vector[i], mode_binder);

    return res;
}

__MSNATIVE_ inline ityp __export math_media(uint64_t dim, ityp vector[static dim])
{
    return (summation(dim, SUMMATION_SUM, vector)/dim);
}

__MSNATIVE_ inline ityp __export math_mode(uint64_t dim, ityp vector[static dim])
{
	ityp vmoda[dim];
	dim_typ moda;
	register uint64_t i, j;
	#pragma omp parallel for
	for(i=0; i<dim; ++i)
		for(j=vmoda[i]=0; j!=i && j<dim; ++j)
			if(vector[j] == vector[i])
				++ vmoda[i];
	(void) MINMAX(dim, vmoda, MAX_MODE, &moda);
	return vector[moda];
}

__MSNATIVE_ inline ityp __export math_variance(uint64_t dim, ityp vector[static dim])
{
	const register ityp media = math_media(dim, vector);
	ityp register res = 0.00;
	
	#pragma omp parallel for
	for(uint64_t i=0; i<dim; ++i)
		res += exp2(vector[i] - media);
		
	return res/dim;
}

__MSNATIVE_ inline ityp __export math_covariance(uint64_t dim, ityp vector1[static dim], ityp vector2[static dim])
{
	const register ityp media[MAX_DIMENSIONS] =
	{
		math_media(dim, vector1),
		math_media(dim, vector2)
	};
	
	ityp register res = 0.00;
	
	#pragma omp parallel for
	for(uint64_t i=0; i<dim; ++i)
		res += (vector1[i]-media[FIRST_NUMBER])*(vector2[i]-media[SECOND_NUMBER]);
		
	return res/dim;
}

__MSNATIVE_ inline ityp __export math_stddev(uint64_t dim, ityp vector[static dim])
{
	const register ityp media = math_media(dim, vector);
	ityp register res = 0.00;
	
	#pragma omp parallel for
	for(uint64_t i=0; i<dim; ++i)
		res += exp2(vector[i] - media);
		
	return sqrt(res/dim);
}

__MSNATIVE_ inline bool __export math_outlier(uint64_t dim, uint64_t outlier_idx, ityp vector[static dim])
{
	return math_outlier2(dim, outlier_idx, access(curLayout)->outlier_constant, vector);
}

__MSNATIVE_ bool __export math_outlier2(uint64_t dim, uint64_t outlier_idx, float outlier_constant, ityp vector[static dim])
{
	qsort(vector, dim, sizeof(ityp), cmpfunc);
	const register ityp Q1 = ((dim%2) ? vector[(uint64_t)((dim-3)*FIRST_QUARTILE_CONSTANT)] : ((vector[(uint64_t)((dim-4)*FIRST_QUARTILE_CONSTANT)]+vector[(uint64_t)(dim*FIRST_QUARTILE_CONSTANT)])*0.5));
	const register ityp Q3 = ((dim%2) ? vector[(uint64_t)(((dim+1)*THIRD_QUARTILE_CONSTANT)-1)] : ((vector[(uint64_t)((dim*THIRD_QUARTILE_CONSTANT)-1)]+vector[(uint64_t)(dim*THIRD_QUARTILE_CONSTANT)])*0.5));
	const register ityp deviation = (Q3-Q1)*outlier_constant;
	return(vector[outlier_idx] < Q1-deviation || vector[outlier_idx] > Q3+deviation);
}

__MSNATIVE_ inline ityp __export math_geomedia(uint64_t dim, ityp vector[static dim])
{
    return (pow(productory(dim, PRODUCTORY_MUL, vector), (1/dim)));
}

__MSNATIVE_ inline ityp __export math_armedia(uint64_t dim, ityp vector[static dim])
{
    #pragma omp parallel for
    for(uint64_t i=0; i<dim; ++i)
        vector[i] = (1/vector[i]);

    return (dim/(summation(dim, SUMMATION_SUM, vector)));

}

__MSNATIVE_ inline ityp __export math_powmedia(uint64_t dim, uint64_t power, ityp vector[static dim])
{
	#pragma omp parallel for
    for(uint64_t i=0; i<dim; ++i)
        vector[i] = mpow2(vector[i], power);

    return (pow(summation(dim, SUMMATION_SUM, vector), (1/power)));
}

__MSNATIVE_ inline ityp __export math_scale(uint64_t dim, ityp vector[static dim])
{
    return ((MIN(dim, vector)+MAX(dim, vector))*0.5);
}

__MSNATIVE_ inline ityp __export math_first_quartile(uint64_t dim, ityp vector[static dim])
{
    qsort(vector, dim, sizeof(ityp), cmpfunc);
    return ((dim%2) ? vector[(uint64_t)((dim-3)*FIRST_QUARTILE_CONSTANT)] : ((vector[(uint64_t)((dim-4)*FIRST_QUARTILE_CONSTANT)]+vector[(uint64_t)(dim*FIRST_QUARTILE_CONSTANT)])*0.5));
}

__MSNATIVE_ inline ityp __export math_mediana(uint64_t dim, ityp vector[static dim])
{
    qsort(vector, dim, sizeof(ityp), cmpfunc);
    return ((dim%2) ? vector[(uint64_t)((dim-1)*SECOND_QUARTILE_CONSTANT)] : ((vector[(uint64_t)((dim-2)*SECOND_QUARTILE_CONSTANT)]+vector[(uint64_t)(dim*SECOND_QUARTILE_CONSTANT)])*0.5));
}

__MSNATIVE_ inline ityp __export math_third_quartile(uint64_t dim, ityp vector[static dim])
{
    qsort(vector, dim, sizeof(ityp), cmpfunc);
    return ((dim%2) ? vector[(uint64_t)(((dim+1)*THIRD_QUARTILE_CONSTANT)-1)] : ((vector[(uint64_t)((dim*THIRD_QUARTILE_CONSTANT)-1)]+vector[(uint64_t)(dim*THIRD_QUARTILE_CONSTANT)])*0.5));
}

__MSUTIL_ inline int64_t _MS__private __export powi(register int64_t x, register int64_t y) // calculates x^y
{
    int64_t base, res;
    for(base=x, res=1; y; y >>=1, base *= base)
        if(y&1)
            res *= base;
    return res;
}

__MSUTIL_ inline int64_t __system __export changeBase(register int n, sel_typ start_base, sel_typ base)
{

    uint64_t h, sum;

    for(h=sum=0; n>0; ++h)
    {
        sum += (n%base)*powi(start_base, h);
        n /= base;
    }
    return sum;
}

__MSUTIL_ inline uint64_t __export math_MCD(register uint64_t a, register uint64_t b) // Euclide's Algorythm
{
    while(a != b)
        if(a > b) a += (-b);
        else b += (-a);
    return b;
}

// thanks elite.polito.it/files/courses/12BHD/progr/Esercizi-C-v2_01.pdf
// I renamed, modified and adapted to this program this function.
__MSUTIL_ inline uint64_t __export math_mcm(register uint64_t a, register uint64_t b)
{
    uint64_t massimo, minimo;
    uint64_t mcm, conta;
    bool fine;
    massimo = a > b ? a : b;
    minimo = a < b ? a : b;
    conta = 1;
    mcm = 0;
    fine = false;
    while(!(fine))
    {
        mcm = conta * massimo;
        if(mcm % minimo == 0) fine = true;
        else ++ conta;
    }
    return mcm;
}

__MSNATIVE_ inline ityp __export exp10(register ityp n)
{
    return pow(10.00, n);
}

__MSNATIVE_ inline ityp __export expc(register ityp n)
{
    return (exp(n)/n);
}

__MSNATIVE_ inline ityp __export exp10c(register ityp n)
{
    return (pow(10.00, n)/n);
}

__MSNATIVE_ inline ityp __export exp2c(register ityp n)
{
    return (exp2(n)/n);
}

__MSUTIL_ inline ityp __export logbN(register ityp b, register ityp N)
{
	return(log(N)/log(b));
}

__MSNATIVE_ inline ityp __export logc(register ityp n)
{
    return (log(n)/n);
}

__MSNATIVE_ inline ityp __export log10c(register ityp n)
{
    return (log10(n)/n);
}

__MSNATIVE_ inline ityp __export log2c(register ityp n)
{
    return (log2(n)/n);
}

__MSNATIVE_ inline ityp __export log1pc(register ityp n)
{
    return (log1p(n)/n);
}

__MSNATIVE_ inline ityp __export log101p(register ityp n)
{
	return (logbN(10.00, 1.00+n));
}

__MSNATIVE_ inline ityp __export log101pc(register ityp n)
{
	return (logbN(10.00, 1.00+n)/n);
}

__MSNATIVE_ inline ityp __export log21p(register ityp n)
{
	return (logbN(2.00, 1.00+n));
}

__MSNATIVE_ inline ityp __export log21pc(register ityp n)
{
	return (logbN(2.00, 1.00+n)/n);
}

__MSUTIL_ inline ityp __export rootnX(register ityp n, register ityp X)
{
    return(pow(X, 1/n));
}

__MSNATIVE_ inline double complex __export cexpc(register double complex n)
{
    return (cexp(n)/n);
}

__MSNATIVE_ inline double complex __export cexp10(register double complex n)
{
    return cpow(10.00, n);
}

__MSNATIVE_ inline double complex __export cexp10c(register double complex n)
{
    return (cpow(10.00, n)/n);
}

__MSNATIVE_ inline double complex __export cexp2(register double complex n)
{
	return (cpow(2.00, n));
}

__MSNATIVE_ inline double complex __export cexp2c(register double complex n)
{
    return (cpow(2.00, n)/n);
}

__MSUTIL_ inline double complex __export clogbN(register double complex b, register double complex N)
{
	return(clog(N)/clog(b));
}

__MSNATIVE_ inline double complex __export clogc(register double complex n)
{
    return (clog(n)/n);
}

__MSNATIVE_ inline double complex __export clog10(register double complex n)
{
	return (clogbN(10.00, n));
}

__MSNATIVE_ inline double complex __export clog10c(register double complex n)
{
    return (clogbN(10.00, n)/n);
}

__MSNATIVE_ inline double complex __export clog2(register double complex n)
{
	return (clogbN(2.00, n));
}

__MSNATIVE_ inline double complex __export clog2c(register double complex n)
{
    return (clogbN(2.00, n)/n);
}

__MSNATIVE_ inline double complex __export clog1p(register double complex n)
{
	return (clog(1.00+n));
}

__MSNATIVE_ inline double complex __export clog1pc(register double complex n)
{
    return (clog(1.00+n)/n);
}

__MSNATIVE_ inline double complex __export clog101p(register double complex n)
{
	return (clogbN(10.00, 1.00+n));
}

__MSNATIVE_ inline double complex __export clog101pc(register double complex n)
{
    return (clogbN(10.00, 1.00+n)/n);
}

__MSNATIVE_ inline double complex __export clog21p(register double complex n)
{
	return (clogbN(2.00, 1.00+n));
}

__MSNATIVE_ inline double complex __export clog21pc(register double complex n)
{
    return (clogbN(2.00, 1.00+n)/n);
}

__MSUTIL_ inline double complex __export ccbrt(register double complex n)
{
    return(cpow(n, 1/3));
}

__MSUTIL_ inline double complex __export crootnX(register double complex n, register double complex X)
{
    return(cpow(X, 1/n));
}


// Properly Trigonometric Functions DEFINITIONS
// (by inlining them)

// It converts Radiants r to Degrees
__MSNATIVE_ inline ityp __export deg(register ityp r)
{
    return ((180.0 * r) / M_PI);
}

// It converts Degrees d to Radiants
__MSNATIVE_ inline ityp __export rad(register ityp d)
{
    return ((M_PI * d) / 180.0);
}

__MSNATIVE_ inline ityp __export csc(register ityp n)
{
    return (1/sin(n));
}

__MSNATIVE_ inline ityp __export sec(register ityp n)
{
    return (1/cos(n));
}

__MSNATIVE_ inline ityp __export cot(register ityp n)
{
    return (1/tan(n));
}

__MSNATIVE_ inline ityp __export csch(register ityp n)
{
    return (1/sinh(n));
}

__MSNATIVE_ inline ityp __export sech(register ityp n)
{
    return (1/cosh(n));
}

__MSNATIVE_ inline ityp __export coth(register ityp n)
{
    return (1/tanh(n));
}

__MSNATIVE_ inline ityp __export acsc(register ityp n)
{
    return asin(1/n);
}

__MSNATIVE_ inline ityp __export asec(register ityp n)
{
    return acos(1/n);
}

__MSNATIVE_ inline ityp __export acot(register ityp n)
{
    return atan(1/n);
}

__MSNATIVE_ inline ityp __export acsch(register ityp n)
{
    return asinh(1/n);
}

__MSNATIVE_ inline ityp __export asech(register ityp n)
{
    return acosh(1/n);
}

__MSNATIVE_ inline ityp __export acoth(register ityp n)
{
    return atanh(1/n);
}

__MSNATIVE_ inline ityp __export hsin(register ityp n)
{
    return sin(n)/2.00;
}

__MSNATIVE_ inline ityp __export hsinh(register ityp n)
{
    return sinh(n)/2.00;
}

__MSNATIVE_ inline ityp __export qsin(register ityp n)
{
    return sin(n)/4.00;
}

__MSNATIVE_ inline ityp __export qsinh(register ityp n)
{
    return sinh(n)/4.00;
}

__MSNATIVE_ inline ityp __export hcos(register ityp n)
{
    return cos(n)/2.00;
}

__MSNATIVE_ inline ityp __export hcosh(register ityp n)
{
    return cosh(n)/2.00;
}

__MSNATIVE_ inline ityp __export qcos(register ityp n)
{
    return cos(n)/4.00;
}

__MSNATIVE_ inline ityp __export qcosh(register ityp n)
{
    return cosh(n)/4.00;
}

__MSNATIVE_ inline ityp __export hcsc(register ityp n)
{
    return csc(n)/2.00;
}

__MSNATIVE_ inline ityp __export hcsch(register ityp n)
{
    return csch(n)/2.00;
}

__MSNATIVE_ inline ityp __export qcsc(register ityp n)
{
    return csc(n)/4.00;
}

__MSNATIVE_ inline ityp __export qcsch(register ityp n)
{
    return csch(n)/4.00;
}

__MSNATIVE_ inline ityp __export hsec(register ityp n)
{
    return sec(n)/2.00;
}

__MSNATIVE_ inline ityp __export hsech(register ityp n)
{
    return sech(n)/2.00;
}

__MSNATIVE_ inline ityp __export qsec(register ityp n)
{
    return sec(n)/4.00;
}

__MSNATIVE_ inline ityp __export qsech(register ityp n)
{
    return sech(n)/4.00;
}

__MSNATIVE_ inline ityp __export htan(register ityp n)
{
    return tan(n)/2.00;
}

__MSNATIVE_ inline ityp __export htanh(register ityp n)
{
    return tanh(n)/2.00;
}

__MSNATIVE_ inline ityp __export qtan(register ityp n)
{
    return tan(n)/4.00;
}

__MSNATIVE_ inline ityp __export qtanh(register ityp n)
{
    return tanh(n)/4.00;
}

__MSNATIVE_ inline ityp __export hcot(register ityp n)
{
    return cot(n)/2.00;
}

__MSNATIVE_ inline ityp __export hcoth(register ityp n)
{
    return coth(n)/2.00;
}

__MSNATIVE_ inline ityp __export qcot(register ityp n)
{
    return cot(n)/4.00;
}

__MSNATIVE_ inline ityp __export qcoth(register ityp n)
{
    return coth(n)/4.00;
}

__MSNATIVE_ inline ityp __export vsin(register ityp n)
{
    return 1.00-cos(n);
}

__MSNATIVE_ inline ityp __export cvsin(register ityp n)
{
    return 1.00-sin(n);
}

__MSNATIVE_ inline ityp __export vcos(register ityp n)
{
    return 1.00+cos(n);
}

__MSNATIVE_ inline ityp __export cvcos(register ityp n)
{
    return 1.00+sin(n);
}

__MSNATIVE_ inline ityp __export hvsin(register ityp n)
{
    return (1.00-cos(n))/2.00;
}

__MSNATIVE_ inline ityp __export hcvsin(register ityp n)
{
    return (1.00-sin(n))/2.00;
}

__MSNATIVE_ inline ityp __export hvcos(register ityp n)
{
    return (1.00+cos(n))/2.00;
}

__MSNATIVE_ inline ityp __export hcvcos(register ityp n)
{
    return (1.00+sin(n))/2.00;
}

__MSNATIVE_ inline ityp __export qvsin(register ityp n)
{
    return (1.00-cos(n))/4.00;
}

__MSNATIVE_ inline ityp __export qcvsin(register ityp n)
{
    return (1.00-sin(n))/4.00;
}

__MSNATIVE_ inline ityp __export qvcos(register ityp n)
{
    return (1.00+cos(n))/4.00;
}

__MSNATIVE_ inline ityp __export qcvcos(register ityp n)
{
    return (1.00+sin(n))/4.00;
}

__MSNATIVE_ inline ityp __export vsinh(register ityp n)
{
    return 1.00-cosh(n);
}

__MSNATIVE_ inline ityp __export cvsinh(register ityp n)
{
    return 1.00-sinh(n);
}

__MSNATIVE_ inline ityp __export vcosh(register ityp n)
{
    return 1.00+cosh(n);
}

__MSNATIVE_ inline ityp __export cvcosh(register ityp n)
{
    return 1.00+sinh(n);
}

__MSNATIVE_ inline ityp __export hvsinh(register ityp n)
{
    return (1.00-cosh(n))/2.00;
}

__MSNATIVE_ inline ityp __export hcvsinh(register ityp n)
{
    return (1.00-sinh(n))/2.00;
}

__MSNATIVE_ inline ityp __export hvcosh(register ityp n)
{
    return (1.00+cosh(n))/2.00;
}

__MSNATIVE_ inline ityp __export hcvcosh(register ityp n)
{
    return (1.00+sinh(n))/2.00;
}

__MSNATIVE_ inline ityp __export qvsinh(register ityp n)
{
    return (1.00-cosh(n))/4.00;
}

__MSNATIVE_ inline ityp __export qcvsinh(register ityp n)
{
    return (1.00-sinh(n))/4.00;
}

__MSNATIVE_ inline ityp __export qvcosh(register ityp n)
{
    return (1.00+cosh(n))/4.00;
}

__MSNATIVE_ inline ityp __export qcvcosh(register ityp n)
{
    return (1.00+sinh(n))/4.00;
}

__MSNATIVE_ inline ityp __export esec(register ityp n)
{
    return sec(n)-1.00;
}

__MSNATIVE_ inline ityp __export ecsc(register ityp n)
{
    return csc(n)-1.00;
}

__MSNATIVE_ inline ityp __export esech(register ityp n)
{
    return sech(n)-1.00;
}

__MSNATIVE_ inline ityp __export ecsch(register ityp n)
{
    return csch(n)-1.00;
}

__MSNATIVE_ inline ityp __export hesec(register ityp n)
{
    return (sec(n)-1.00)/2.00;
}

__MSNATIVE_ inline ityp __export hecsc(register ityp n)
{
    return (csc(n)-1.00)/2.00;
}

__MSNATIVE_ inline ityp __export hesech(register ityp n)
{
    return (sech(n)-1.00)/2.00;
}

__MSNATIVE_ inline ityp __export hecsch(register ityp n)
{
    return (csch(n)-1.00)/2.00;
}

__MSNATIVE_ inline ityp __export qesec(register ityp n)
{
    return (sec(n)-1.00)/4.00;
}

__MSNATIVE_ inline ityp __export qecsc(register ityp n)
{
    return (csc(n)-1.00)/4.00;
}

__MSNATIVE_ inline ityp __export qesech(register ityp n)
{
    return (sech(n)-1.00)/4.00;
}

__MSNATIVE_ inline ityp __export qecsch(register ityp n)
{
    return (csch(n)-1.00)/4.00;
}

__MSNATIVE_ inline ityp __export sinc(register ityp n)
{
    return n ? sin(n)/n : 1.00;
}

__MSNATIVE_ inline ityp __export sinch(register ityp n)
{
    return n ? sinh(n)/n : 1.00;
}

__MSNATIVE_ inline ityp __export hsinc(register ityp n)
{
    return n ? sin(n)/(2.00*n) : 0.50;
}

__MSNATIVE_ inline ityp __export hsinch(register ityp n)
{
    return n ? sinh(n)/(2.00*n) : 0.50;
}

__MSNATIVE_ inline ityp __export qsinc(register ityp n)
{
    return n ? sin(n)/(4.00*n) : 0.25;
}

__MSNATIVE_ inline ityp __export qsinch(register ityp n)
{
    return n ? sinh(n)/(4.00*n) : 0.25;
}

__MSNATIVE_ inline ityp __export cosc(register ityp n)
{
    return cos(n)/n;
}

__MSNATIVE_ inline ityp __export cosch(register ityp n)
{
    return cosh(n)/n;
}

__MSNATIVE_ inline ityp __export hcosc(register ityp n)
{
    return cos(n)/(2.00*n);
}

__MSNATIVE_ inline ityp __export hcosch(register ityp n)
{
    return cosh(n)/(2.00*n);
}

__MSNATIVE_ inline ityp __export qcosc(register ityp n)
{
    return cos(n)/(4.00*n);
}

__MSNATIVE_ inline ityp __export qcosch(register ityp n)
{
    return cosh(n)/(4.00*n);
}

__MSNATIVE_ inline ityp __export secc(register ityp n)
{
    return sec(n)/n;
}

__MSNATIVE_ inline ityp __export secch(register ityp n)
{
    return sech(n)/n;
}

__MSNATIVE_ inline ityp __export hsecc(register ityp n)
{
    return sec(n)/(2.00*n);
}

__MSNATIVE_ inline ityp __export hsecch(register ityp n)
{
    return sech(n)/(2.00*n);
}

__MSNATIVE_ inline ityp __export qsecc(register ityp n)
{
    return sec(n)/(4.00*n);
}

__MSNATIVE_ inline ityp __export qsecch(register ityp n)
{
    return sech(n)/(4*n);
}

__MSNATIVE_ inline ityp __export cscc(register ityp n)
{
    return csc(n)/n;
}

__MSNATIVE_ inline ityp __export cscch(register ityp n)
{
    return csch(n)/n;
}

__MSNATIVE_ inline ityp __export hcscc(register ityp n)
{
    return csc(n)/(2.00*n);
}

__MSNATIVE_ inline ityp __export hcscch(register ityp n)
{
    return csch(n)/(2.00*n);
}

__MSNATIVE_ inline ityp __export qcscc(register ityp n)
{
    return csc(n)/(4.00*n);
}

__MSNATIVE_ inline ityp __export qcscch(register ityp n)
{
    return csch(n)/(4.00*n);
}

__MSNATIVE_ inline ityp __export tanc(register ityp n)
{
    return tan(n)/n;
}

__MSNATIVE_ inline ityp __export tanch(register ityp n)
{
    return tanh(n)/n;
}

__MSNATIVE_ inline ityp __export htanc(register ityp n)
{
    return tan(n)/(2.00*n);
}

__MSNATIVE_ inline ityp __export htanch(register ityp n)
{
    return tanh(n)/(2.00*n);
}

__MSNATIVE_ inline ityp __export qtanc(register ityp n)
{
    return tan(n)/(4.00*n);
}

__MSNATIVE_ inline ityp __export qtanch(register ityp n)
{
    return tanh(n)/(4.00*n);
}

__MSNATIVE_ inline ityp __export cotc(register ityp n)
{
    return cot(n)/n;
}

__MSNATIVE_ inline ityp __export cotch(register ityp n)
{
    return coth(n)/n;
}

__MSNATIVE_ inline ityp __export hcotc(register ityp n)
{
    return cot(n)/(2.00*n);
}

__MSNATIVE_ inline ityp __export hcotch(register ityp n)
{
    return coth(n)/(2.00*n);
}

__MSNATIVE_ inline ityp __export qcotc(register ityp n)
{
    return cot(n)/(4.00*n);
}

__MSNATIVE_ inline ityp __export qcotch(register ityp n)
{
    return coth(n)/(4.00*n);
}

__MSNATIVE_ inline double complex __export ccsc(register double complex n)
{
    return (1/csin(n));
}

__MSNATIVE_ inline double complex __export csec(register double complex n)
{
    return (1/ccos(n));
}

__MSNATIVE_ inline double complex __export ccot(register double complex n)
{
    return (1/ctan(n));
}

__MSNATIVE_ inline double complex __export ccsch(register double complex n)
{
    return (1/csinh(n));
}

__MSNATIVE_ inline double complex __export csech(register double complex n)
{
    return (1/ccosh(n));
}

__MSNATIVE_ inline double complex __export ccoth(register double complex n)
{
    return (1/ctanh(n));
}

__MSNATIVE_ inline double complex __export cacsc(register double complex n)
{
    return casin(1/n);
}

__MSNATIVE_ inline double complex __export casec(register double complex n)
{
    return cacos(1/n);
}

__MSNATIVE_ inline double complex __export cacot(register double complex n)
{
    return catan(1/n);
}

__MSNATIVE_ inline double complex __export cacsch(register double complex n)
{
    return casinh(1/n);
}

__MSNATIVE_ inline double complex __export casech(register double complex n)
{
    return cacosh(1/n);
}

__MSNATIVE_ inline double complex __export cacoth(register double complex n)
{
    return catanh(1/n);
}

__MSNATIVE_ inline double complex __export chsin(register double complex n)
{
    return csin(n)/2.00;
}

__MSNATIVE_ inline double complex __export chsinh(register double complex n)
{
    return csinh(n)/2.00;
}

__MSNATIVE_ inline double complex __export cqsin(register double complex n)
{
    return csin(n)/4.00;
}

__MSNATIVE_ inline double complex __export cqsinh(register double complex n)
{
    return csinh(n)/4.00;
}

__MSNATIVE_ inline double complex __export chcos(register double complex n)
{
    return ccos(n)/2.00;
}

__MSNATIVE_ inline double complex __export chcosh(register double complex n)
{
    return ccosh(n)/2.00;
}

__MSNATIVE_ inline double complex __export cqcos(register double complex n)
{
    return ccos(n)/4.00;
}

__MSNATIVE_ inline double complex __export cqcosh(register double complex n)
{
    return ccosh(n)/4.00;
}

__MSNATIVE_ inline double complex __export chcsc(register double complex n)
{
    return ccsc(n)/2.00;
}

__MSNATIVE_ inline double complex __export chcsch(register double complex n)
{
    return ccsch(n)/2.00;
}

__MSNATIVE_ inline double complex __export cqcsc(register double complex n)
{
    return ccsc(n)/4.00;
}

__MSNATIVE_ inline double complex __export cqcsch(register double complex n)
{
    return ccsch(n)/4.00;
}

__MSNATIVE_ inline double complex __export chsec(register double complex n)
{
    return csec(n)/2.00;
}

__MSNATIVE_ inline double complex __export chsech(register double complex n)
{
    return csech(n)/2.00;
}

__MSNATIVE_ inline double complex __export cqsec(register double complex n)
{
    return csec(n)/4.00;
}

__MSNATIVE_ inline double complex __export cqsech(register double complex n)
{
    return csech(n)/4.00;
}

__MSNATIVE_ inline double complex __export chtan(register double complex n)
{
    return ctan(n)/2.00;
}

__MSNATIVE_ inline double complex __export chtanh(register double complex n)
{
    return ctanh(n)/2.00;
}

__MSNATIVE_ inline double complex __export cqtan(register double complex n)
{
    return ctan(n)/4.00;
}

__MSNATIVE_ inline double complex __export cqtanh(register double complex n)
{
    return ctanh(n)/4.00;
}

__MSNATIVE_ inline double complex __export chcot(register double complex n)
{
    return ccot(n)/2.00;
}

__MSNATIVE_ inline double complex __export chcoth(register double complex n)
{
    return ccoth(n)/2.00;
}

__MSNATIVE_ inline double complex __export cqcot(register double complex n)
{
    return ccot(n)/4.00;
}

__MSNATIVE_ inline double complex __export cqcoth(register double complex n)
{
    return ccoth(n)/4.00;
}

__MSNATIVE_ inline double complex __export cpxvsin(register double complex n)
{
    return 1.00-ccos(n);
}

__MSNATIVE_ inline double complex __export ccvsin(register double complex n)
{
    return 1.00-csin(n);
}

__MSNATIVE_ inline double complex __export cpxvcos(register double complex n)
{
    return 1.00+ccos(n);
}

__MSNATIVE_ inline double complex __export ccvcos(register double complex n)
{
    return 1.00+csin(n);
}

__MSNATIVE_ inline double complex __export chvsin(register double complex n)
{
    return (1.00-ccos(n))/2.00;
}

__MSNATIVE_ inline double complex __export chcvsin(register double complex n)
{
    return (1.00-csin(n))/2.00;
}

__MSNATIVE_ inline double complex __export chvcos(register double complex n)
{
    return (1.00+ccos(n))/2.00;
}

__MSNATIVE_ inline double complex __export chcvcos(register double complex n)
{
    return (1.00+csin(n))/2.00;
}

__MSNATIVE_ inline double complex __export cqvsin(register double complex n)
{
    return (1.00-ccos(n))/4.00;
}

__MSNATIVE_ inline double complex __export cqcvsin(register double complex n)
{
    return (1.00-csin(n))/4.00;
}

__MSNATIVE_ inline double complex __export cqvcos(register double complex n)
{
    return (1.00+ccos(n))/4.00;
}

__MSNATIVE_ inline double complex __export cqcvcos(register double complex n)
{
    return (1.00+csin(n))/4.00;
}

__MSNATIVE_ inline double complex __export cpxvsinh(register double complex n)
{
    return 1.00-ccosh(n);
}

__MSNATIVE_ inline double complex __export ccvsinh(register double complex n)
{
    return 1.00-csinh(n);
}

__MSNATIVE_ inline double complex __export cpxvcosh(register double complex n)
{
    return 1.00+ccosh(n);
}

__MSNATIVE_ inline double complex __export ccvcosh(register double complex n)
{
    return 1.00+csinh(n);
}

__MSNATIVE_ inline double complex __export chvsinh(register double complex n)
{
    return (1.00-ccosh(n))/2.00;
}

__MSNATIVE_ inline double complex __export chcvsinh(register double complex n)
{
    return (1.00-csinh(n))/2.00;
}

__MSNATIVE_ inline double complex __export chvcosh(register double complex n)
{
    return (1.00+ccosh(n))/2.00;
}

__MSNATIVE_ inline double complex __export chcvcosh(register double complex n)
{
    return (1.00+csinh(n))/2.00;
}

__MSNATIVE_ inline double complex __export cqvsinh(register double complex n)
{
    return (1.00-ccosh(n))/4.00;
}

__MSNATIVE_ inline double complex __export cqcvsinh(register double complex n)
{
    return (1.00-csinh(n))/4.00;
}

__MSNATIVE_ inline double complex __export cqvcosh(register double complex n)
{
    return (1.00+ccosh(n))/4.00;
}

__MSNATIVE_ inline double complex __export cqcvcosh(register double complex n)
{
    return (1.00+csinh(n))/4.00;
}

__MSNATIVE_ inline double complex __export cesec(register double complex n)
{
    return csec(n)-1.00;
}

__MSNATIVE_ inline double complex __export cecsc(register double complex n)
{
    return ccsc(n)-1.00;
}

__MSNATIVE_ inline double complex __export cesech(register double complex n)
{
    return csech(n)-1.00;
}

__MSNATIVE_ inline double complex __export cecsch(register double complex n)
{
    return ccsch(n)-1.00;
}

__MSNATIVE_ inline double complex __export chesec(register double complex n)
{
    return (csec(n)-1.00)/2.00;
}

__MSNATIVE_ inline double complex __export checsc(register double complex n)
{
    return (ccsc(n)-1.00)/2.00;
}

__MSNATIVE_ inline double complex __export chesech(register double complex n)
{
    return (csech(n)-1.00)/2.00;
}

__MSNATIVE_ inline double complex __export checsch(register double complex n)
{
    return (ccsch(n)-1.00)/2.00;
}

__MSNATIVE_ inline double complex __export cqesec(register double complex n)
{
    return (csec(n)-1.00)/4.00;
}

__MSNATIVE_ inline double complex __export cqecsc(register double complex n)
{
    return (ccsc(n)-1.00)/4.00;
}

__MSNATIVE_ inline double complex __export cqesech(register double complex n)
{
    return (csech(n)-1.00)/4.00;
}

__MSNATIVE_ inline double complex __export cqecsch(register double complex n)
{
    return (ccsch(n)-1.00)/4.00;
}

__MSNATIVE_ inline double complex __export csinc(register double complex n)
{
    return n ? csin(n)/n : 1.00;
}

__MSNATIVE_ inline double complex __export csinch(register double complex n)
{
    return n ? csinh(n)/n : 1.00;
}

__MSNATIVE_ inline double complex __export chsinc(register double complex n)
{
    return n ? csin(n)/(2.00*n) : 0.50;
}

__MSNATIVE_ inline double complex __export chsinch(register double complex n)
{
    return n ? csinh(n)/(2.00*n) : 0.50;
}

__MSNATIVE_ inline double complex __export cqsinc(register double complex n)
{
    return n ? csin(n)/(4.00*n) : 0.25;
}

__MSNATIVE_ inline double complex __export cqsinch(register double complex n)
{
    return n ? csinh(n)/(4.00*n) : 0.25;
}

__MSNATIVE_ inline double complex __export ccosc(register double complex n)
{
    return ccos(n)/n;
}

__MSNATIVE_ inline double complex __export ccosch(register double complex n)
{
    return ccosh(n)/n;
}

__MSNATIVE_ inline double complex __export chcosc(register double complex n)
{
    return ccos(n)/(2.00*n);
}

__MSNATIVE_ inline double complex __export chcosch(register double complex n)
{
    return ccosh(n)/(2.00*n);
}

__MSNATIVE_ inline double complex __export cqcosc(register double complex n)
{
    return ccos(n)/(4.00*n);
}

__MSNATIVE_ inline double complex __export cqcosch(register double complex n)
{
    return ccosh(n)/(4.00*n);
}

__MSNATIVE_ inline double complex __export csecc(register double complex n)
{
    return csec(n)/n;
}

__MSNATIVE_ inline double complex __export csecch(register double complex n)
{
    return csech(n)/n;
}

__MSNATIVE_ inline double complex __export chsecc(register double complex n)
{
    return csec(n)/(2.00*n);
}

__MSNATIVE_ inline double complex __export chsecch(register double complex n)
{
    return csech(n)/(2.00*n);
}

__MSNATIVE_ inline double complex __export cqsecc(register double complex n)
{
    return csec(n)/(4.00*n);
}

__MSNATIVE_ inline double complex __export cqsecch(register double complex n)
{
    return csech(n)/(4*n);
}

__MSNATIVE_ inline double complex __export ccscc(register double complex n)
{
    return ccsc(n)/n;
}

__MSNATIVE_ inline double complex __export ccscch(register double complex n)
{
    return ccsch(n)/n;
}

__MSNATIVE_ inline double complex __export chcscc(register double complex n)
{
    return ccsc(n)/(2.00*n);
}

__MSNATIVE_ inline double complex __export chcscch(register double complex n)
{
    return ccsch(n)/(2.00*n);
}

__MSNATIVE_ inline double complex __export cqcscc(register double complex n)
{
    return ccsc(n)/(4.00*n);
}

__MSNATIVE_ inline double complex __export cqcscch(register double complex n)
{
    return ccsch(n)/(4.00*n);
}

__MSNATIVE_ inline double complex __export ctanc(register double complex n)
{
    return ctan(n)/n;
}

__MSNATIVE_ inline double complex __export ctanch(register double complex n)
{
    return ctanh(n)/n;
}

__MSNATIVE_ inline double complex __export chtanc(register double complex n)
{
    return ctan(n)/(2.00*n);
}

__MSNATIVE_ inline double complex __export chtanch(register double complex n)
{
    return ctanh(n)/(2.00*n);
}

__MSNATIVE_ inline double complex __export cqtanc(register double complex n)
{
    return ctan(n)/(4.00*n);
}

__MSNATIVE_ inline double complex __export cqtanch(register double complex n)
{
    return ctanh(n)/(4.00*n);
}

__MSNATIVE_ inline double complex __export ccotc(register double complex n)
{
    return ccot(n)/n;
}

__MSNATIVE_ inline double complex __export ccotch(register double complex n)
{
    return ccoth(n)/n;
}

__MSNATIVE_ inline double complex __export chcotc(register double complex n)
{
    return ccot(n)/(2.00*n);
}

__MSNATIVE_ inline double complex __export chcotch(register double complex n)
{
    return ccoth(n)/(2.00*n);
}

__MSNATIVE_ inline double complex __export cqcotc(register double complex n)
{
    return ccot(n)/(4.00*n);
}

__MSNATIVE_ inline double complex __export cqcotch(register double complex n)
{
    return ccoth(n)/(4.00*n);
}


__MSNATIVE_ inline bool __system __export isEqualMatrix(ityp *matrix1, ityp *matrix2, const register dim_typ dim[static MAX_DIMENSIONS])
{
    dim_typ i, j;

    for(i=0; i<dim[ROWS]; ++i)
        for(j=0; j<dim[COLUMNS]; ++j)
            if(*(matrix1 + dim[COLUMNS]*i + j) != *(matrix2 + dim[COLUMNS]*i + j))
                return false;

    return true;

}

/// massive thanks to Bibek Subedi. Link at:
/// http://www.programming-technique.blogspot.it/2011/09/numerical-methods-condition-number-and.html
__MSNATIVE_ __MSUTIL_ ityp __system __export norms(ityp *matrix, dim_typ dim)
{
    dim_typ i, j;
    ityp sum;

    for(i=sum=0; i < dim; ++i)
        for(j = 0; j < dim; ++j)
            sum += mpow2(*(matrix + dim*i + j), 2);

    return sqrt(sum);
}

/// thanks to: apatriarca, a User of Matematicamente forum, who is the author of this
/// pseudo-coded duet of functions which I renamed, unified, modified and adapted to this program.
/// This function calculates the Norm, One or Infinity (depending on mode bool param), of a Square Matrix.
/// http://www.matematicamente.it/forum/viewtopic.php?t=51498&p=371845
__MSNATIVE_ __MSUTIL_ ityp __system __export norm(ityp *matrix, dim_typ dim, bool mode)
{
    dim_typ i, j;
    ityp max; //  = 0.00; /* tutte le somme sono >= 0 */

    for (j = max = 0; j < dim; ++j)
    {
        ityp sum = 0.00;
        for (i = 0; i < dim; ++i)
            sum += fabs(*(matrix + dim*(mode?j:i) + (mode?i:j)));
        if (sum > max)
            max = sum;
    }
    return max;
}

__MSNATIVE_ __MSUTIL_ void __export newtonDifferenceTable(dim_typ n, ityp y[access(curLayout)->max_newton_difftables_dim][access(curLayout)->max_newton_difftables_dim], bool mode)
{
    dim_typ i, j;
    if(mode)
    {
        for(j=0;++j<n;)
            for(i=0;i<(n-j);++i)
                y[i][j] = y[i+1][j-1] - y[i][j-1];
    }
    else
    {
        for(j=0;++j<n;)
            for(i=n-1;i>(j-1);--i)
                y[i][j] = y[i][j-1] - y[i-1][j-1];
    }

    return;
}

// Algebra-Math Properly FUNCTIONS DEFINITIONS
//


// INTRODUCED THESE FUNCTIONS
// in order to Eliminate overheading problem
// into _matrixMultiplication and _matrixAdd functions,
// generated by the Inverse Operations Attribute Check
__MSNATIVE_ static inline ityp __system __export math_sum(register ityp a, register ityp b)
{
    return (a+b);
}

__MSNATIVE_ static inline ityp __system __export math_sub(register ityp a, register ityp b)
{
    return (a-b);
}

__MSNATIVE_ static inline double complex __system __export math_csum(register double complex a, register double complex b)
{
    return (a+b);
}

__MSNATIVE_ static inline double complex __system __export math_csub(register double complex a, register double complex b)
{
    return (a-b);
}

__MSNATIVE_ inline ityp __system __export math_mul(register ityp a, register ityp b)
{
    return (a*b);
}

__MSNATIVE_ inline ityp __system __export math_div(register ityp a, register ityp b)
{
    return (a/b);
}

// It works only with PL Problems with non-negative variables
__MSNATIVE_ sel_typ _MS__private __system __export _simplexMethod(ityp **tableau, ityp **bfs, const register dim_typ dim[static MAX_DIMENSIONS], ityp *constraint_types, bool mode)
{

    dim_typ i, j;
    const dim_typ dim_minus1[MAX_DIMENSIONS] =
    {
        dim[ROWS]-1,
        dim[COLUMNS]-1
    };

    const dim_typ vardims = dim[COLUMNS]+dim_minus1[ROWS];
    
    (*tableau) = realloc((*tableau), sizeof(ityp)*dim[ROWS]*vardims);
    errMem((*tableau), SIMPLEXMETHOD_ALLOC_ERROR);

    /// assert = false;

    for(i=0; i<dim[ROWS]; ++i)
    {
        const bool assert = i == dim_minus1[ROWS];
        *((*tableau)+ vardims*i + vardims-1) = assert ? 0.00 : *((*tableau) + vardims*i + dim_minus1[COLUMNS]);
        #pragma omp parallel for
        for(j=dim_minus1[COLUMNS]; j<vardims-1; ++j)
        	*((*tableau) + vardims*i + j) = assert ? 0.00 : ((1-(2*constraint_types[i]))*(i+dim_minus1[COLUMNS] == j));
    }

    (*bfs) = realloc((*bfs), sizeof(ityp)*(vardims-1));
    errMem((*bfs), SIMPLEXMETHOD_ALLOC_ERROR);


    // INITIALIZING BFS PHASE....
    #pragma omp parallel for
    for(i=0; i<dim_minus1[COLUMNS]; ++i)
    	*((*bfs) + i) = 0.00;

    j=0;
    for(i=dim_minus1[COLUMNS]; i<vardims-1; ++i)
    	*((*bfs) + i) = *((*tableau) + (vardims*(j++)) + vardims-1);


    ityp func_vector[vardims];
    ityp min_ratios_vector[dim_minus1[ROWS]];
    dim_typ itemcheck;

    bool (* const _cmpfunc)(const register ityp, const register ityp) = MAX_PROBLEM ? max_cmpfunc : min_cmpfunc;

    dim_typ bestvaltax_idx, leaving_var_idx;
    bestvaltax_idx = leaving_var_idx = vardims;
    bool another_iteration = false;

    for(dim_typ control=0; ; )
    {
        if(++control >= access(curLayout)->max_simplex_iterations)
            return SIMPLEXMETHOD_FARBFS_ERROR;

        // vector reversing operation...
        for(i=another_iteration=0; i<vardims; ++i)
            if(_cmpfunc((func_vector[i]=*((*tableau) + vardims*dim_minus1[ROWS] + i)), 0.00))
                another_iteration = true;

        if(!another_iteration)
            break;

        (void) MINMAX(vardims, func_vector, mode, &bestvaltax_idx);

        for(i=itemcheck=0; i<dim_minus1[ROWS]; ++i)
        {
            if(*((*tableau) + vardims*i + bestvaltax_idx) <= 0)
            {
                if(++itemcheck == dim_minus1[ROWS])
                    return SIMPLEXMETHOD_INFBFS_ERROR;
                min_ratios_vector[i] = MAX_VAL;
            }
            else
                min_ratios_vector[i] = *((*tableau) + vardims*i + vardims-1)/ *((*tableau) + vardims*i + bestvaltax_idx); // ((*tableau)[i][vardims-1])/((*tableau)[i][bestvaltax_idx]);
        }

        (void) MINMAX(dim_minus1[ROWS], min_ratios_vector, MIN_MODE, &leaving_var_idx);

        const ityp pivot = *((*tableau) + vardims*leaving_var_idx + bestvaltax_idx);

        // gauss-jordan's method application to the tableau, inclusive of function coefficients.
        for(i=0; i<vardims; ++i)
            *((*tableau) + vardims*leaving_var_idx + i) /= pivot;
            
        ityp zero_binder = 0.00;

        for(i=0; i<dim[ROWS]; ++i)
        {
            if(i == leaving_var_idx || ! *((*tableau) + vardims*i + bestvaltax_idx))
                continue;
		
			*((*tableau) + vardims*i + bestvaltax_idx); 
			*((*tableau) + vardims*leaving_var_idx + bestvaltax_idx);
			
            zero_binder = (1-(((*((*tableau) + vardims*i + bestvaltax_idx)* *((*tableau) + vardims*leaving_var_idx + bestvaltax_idx))>0))<<1)*(*((*tableau) + vardims*i + bestvaltax_idx))/(*((*tableau) + vardims*leaving_var_idx + bestvaltax_idx));

			for(j=0; j<vardims; ++j)
				*((*tableau) + vardims*i + j) += *((*tableau) + vardims*leaving_var_idx + j)*zero_binder;
        }

        j=0;
        *((*bfs) + bestvaltax_idx) = *((*tableau) + vardims*leaving_var_idx + vardims-1);
        for(i=dim_minus1[COLUMNS]; i<vardims-1; ++i)
            if(i == dim_minus1[COLUMNS]+leaving_var_idx)
            	*((*bfs) + i) = 0.00;
            else if(*((*bfs) + i))
            	*((*bfs) + i) = *((*tableau) + (vardims*(j++)) + vardims-1);
    }

    return SIMPLEXMETHOD_FOUNDBFS;
}

__MSNATIVE_ __MSUTIL_ sel_typ _MS__private __system __export _matrixEigenValues(ityp *restrict m, ityp *restrict l, ityp *restrict vc, const register dim_typ n)
{
	ityp *m2 = NULL;
	if(!matrixAlloc(&m2, (dim_typ2){n,n}))
		return EIGVALUES_ALLOC_ERROR; 
	{
		register dim_typ i=n;
		while(i--)
		{
			register dim_typ j=n;
			while(j--)
			{
				*(m2 + n*i + j) = *(m + n*i + j);
				*(vc + n*i + j) = i==j;
			}
		}
	}
	register dim_typ cnt = 0;
	while(true)
	{
		if(++cnt >= access(curLayout)->max_eigvalues_iterations)
			return EIGVALUES_INFEVS_ERROR;
		ityp mod = 0;
		register dim_typ i=0, j=0;
		{
			register dim_typ k=n;
			while(k--)
			{
				register dim_typ m=n;
				while((--m)>k)
				{
					ityp q = fabs(*(m2 + n*k + m));
  					if(q > mod)
					{
						mod=q; 
						i=k;
						j=m;
					}
				}
			}
		}
		if(mod < EIGENVALUES_PREC) break;
		ityp th = 0.5*atan(MAX_DIMENSIONS* *(m2 + n*i + j)/(*(m2 + n*i + i) - *(m2 + n*j + j)));
		{
			ityp c = cos(th), s = sin(th);
			inline void twst(ityp *restrict m)
			{
				register dim_typ k=n;
				while(k--)
				{
					ityp t = (c* *(m + n*i + k)) + (s* *(m+ n*j + k));
   					*(m + n*j + k) = -(s* *(m + n*i + k))+(c* *(m + n*j + k));
					*(m + n*i + k) = t;
				}
			}
			{
				register dim_typ k=n;
				while(k--)
				{
					ityp t = (c* *(m2 + n*k + i)) + (s* *(m2 + n*k + j));
					*(m2 + n*k + j) = -(s* *(m2 + n*k + i)) + (c* *(m2 + n*k + j));
					*(m2 + n*k + i) = t;
				}
			}
			twst(m2);
			twst(vc);
		}
	}
	{
		register dim_typ j=n;
		while(j--) l[j] = *(m2 + n*j + j);
	}
	return EIGVALUES_FOUNDEVS;		
}

__MSNATIVE_ inline void _MS__private __system __export _matrixSub(ityp **matrix1, ityp **matrix2, ityp **matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	register dim_typ idx;
    
	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
    	{
        	idx = dim[COLUMNS]*i + j;
        	*((*matrix_sum) + idx) = *((*matrix1) + idx) - *((*matrix2) + idx);
    	}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixCSub(ityp **matrix1, ityp **matrix2, ityp **matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	register dim_typ idx;
	
	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
		{
			idx = dim[COLUMNS]*i + j;
			const register double complex tmpres = (matrix1[REAL_PART][idx] + matrix1[IMAG_PART][idx]*I) - (matrix2[REAL_PART][idx] + matrix2[IMAG_PART][idx]*I);
        	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
       		{
	            matrix_sum[REAL_PART][idx] = creal(tmpres);
	            matrix_sum[IMAG_PART][idx] = cimag(tmpres);
	    	}
		}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixQSub(ityp **matrix1, ityp **matrix2, ityp **matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	register dim_typ idx;
    

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
    	{
    		idx = dim[COLUMNS]*i + j;
        	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
        	{
	            matrix_sum[QUATERNIONS_REALPART][idx] = matrix1[QUATERNIONS_REALPART][idx] - matrix2[QUATERNIONS_REALPART][idx];
	            matrix_sum[QUATERNIONS_IPART][idx] = matrix1[QUATERNIONS_IPART][idx] - matrix2[QUATERNIONS_IPART][idx];
	            matrix_sum[QUATERNIONS_JPART][idx] = matrix1[QUATERNIONS_JPART][idx] - matrix2[QUATERNIONS_JPART][idx];
	            matrix_sum[QUATERNIONS_KPART][idx] = matrix1[QUATERNIONS_KPART][idx] - matrix2[QUATERNIONS_KPART][idx];
	    	}
		}
		
    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixOSub(ityp **matrix1, ityp **matrix2, ityp **matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	register dim_typ idx;

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
    	{
        	idx = dim[COLUMNS]*i + j;
        	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
			{
				matrix_sum[OCTONIONS_REALPART][idx] = matrix1[OCTONIONS_REALPART][idx] + matrix2[OCTONIONS_REALPART][idx];
				matrix_sum[OCTONIONS_E1PART][idx] = matrix1[OCTONIONS_E1PART][idx] + matrix2[OCTONIONS_E1PART][idx];
				matrix_sum[OCTONIONS_E2PART][idx] = matrix1[OCTONIONS_E2PART][idx] + matrix2[OCTONIONS_E2PART][idx];
				matrix_sum[OCTONIONS_E3PART][idx] = matrix1[OCTONIONS_E3PART][idx] + matrix2[OCTONIONS_E3PART][idx];
				matrix_sum[OCTONIONS_E4PART][idx] = matrix1[OCTONIONS_E4PART][idx] + matrix2[OCTONIONS_E4PART][idx];
				matrix_sum[OCTONIONS_E5PART][idx] = matrix1[OCTONIONS_E5PART][idx] + matrix2[OCTONIONS_E5PART][idx];
				matrix_sum[OCTONIONS_E6PART][idx] = matrix1[OCTONIONS_E6PART][idx] + matrix2[OCTONIONS_E6PART][idx];
				matrix_sum[OCTONIONS_E7PART][idx] = matrix1[OCTONIONS_E7PART][idx] + matrix2[OCTONIONS_E7PART][idx];
			}
		}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixSSub(ityp **matrix1, ityp **matrix2, ityp **matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	register dim_typ idx;

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
        {
	        idx = dim[COLUMNS]*i + j;
	        #pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
			{
				matrix_sum[SEDENIONS_REALPART][idx] = matrix1[SEDENIONS_REALPART][idx] + matrix2[SEDENIONS_REALPART][idx];
		        matrix_sum[SEDENIONS_E1PART][idx] = matrix1[SEDENIONS_E1PART][idx] + matrix2[SEDENIONS_E1PART][idx];
		        matrix_sum[SEDENIONS_E2PART][idx] = matrix1[SEDENIONS_E2PART][idx] + matrix2[SEDENIONS_E2PART][idx];
		        matrix_sum[SEDENIONS_E3PART][idx] = matrix1[SEDENIONS_E3PART][idx] + matrix2[SEDENIONS_E3PART][idx];
		        matrix_sum[SEDENIONS_E4PART][idx] = matrix1[SEDENIONS_E4PART][idx] + matrix2[SEDENIONS_E4PART][idx];
		        matrix_sum[SEDENIONS_E5PART][idx] = matrix1[SEDENIONS_E5PART][idx] + matrix2[SEDENIONS_E5PART][idx];
		        matrix_sum[SEDENIONS_E6PART][idx] = matrix1[SEDENIONS_E6PART][idx] + matrix2[SEDENIONS_E6PART][idx];
		        matrix_sum[SEDENIONS_E7PART][idx] = matrix1[SEDENIONS_E7PART][idx] + matrix2[SEDENIONS_E7PART][idx];
		        matrix_sum[SEDENIONS_E8PART][idx] = matrix1[SEDENIONS_E8PART][idx] + matrix2[SEDENIONS_E8PART][idx];
		        matrix_sum[SEDENIONS_E9PART][idx] = matrix1[SEDENIONS_E9PART][idx] + matrix2[SEDENIONS_E9PART][idx];
		        matrix_sum[SEDENIONS_E10PART][idx] = matrix1[SEDENIONS_E10PART][idx] + matrix2[SEDENIONS_E10PART][idx];
		        matrix_sum[SEDENIONS_E11PART][idx] = matrix1[SEDENIONS_E11PART][idx] + matrix2[SEDENIONS_E11PART][idx];
		        matrix_sum[SEDENIONS_E12PART][idx] = matrix1[SEDENIONS_E12PART][idx] + matrix2[SEDENIONS_E12PART][idx];
		        matrix_sum[SEDENIONS_E13PART][idx] = matrix1[SEDENIONS_E13PART][idx] + matrix2[SEDENIONS_E13PART][idx];
		        matrix_sum[SEDENIONS_E14PART][idx] = matrix1[SEDENIONS_E14PART][idx] + matrix2[SEDENIONS_E14PART][idx];
		        matrix_sum[SEDENIONS_E15PART][idx] = matrix1[SEDENIONS_E15PART][idx] + matrix2[SEDENIONS_E15PART][idx];
			}
		}

    return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixAdd(ityp **matrix1, ityp **matrix2, ityp **matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	register dim_typ idx;
    
	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
    	{
        	idx = dim[COLUMNS]*i + j;
        	*((*matrix_sum) + idx) = *((*matrix1) + idx) + *((*matrix2) + idx);
    	}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixCAdd(ityp **matrix1, ityp **matrix2, ityp **matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	register dim_typ idx;
	
	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
		{
			idx = dim[COLUMNS]*i + j;
			const register double complex tmpres = (matrix1[REAL_PART][idx] + matrix1[IMAG_PART][idx]*I) + (matrix2[REAL_PART][idx] + matrix2[IMAG_PART][idx]*I);
        	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
       		{
	            matrix_sum[REAL_PART][idx] = creal(tmpres);
	            matrix_sum[IMAG_PART][idx] = cimag(tmpres);
	    	}
		}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixQAdd(ityp **matrix1, ityp **matrix2, ityp **matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	register dim_typ idx;
    

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
    	{
    		idx = dim[COLUMNS]*i + j;
        	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
        	{
	            matrix_sum[QUATERNIONS_REALPART][idx] = matrix1[QUATERNIONS_REALPART][idx] + matrix2[QUATERNIONS_REALPART][idx];
	            matrix_sum[QUATERNIONS_IPART][idx] = matrix1[QUATERNIONS_IPART][idx] + matrix2[QUATERNIONS_IPART][idx];
	            matrix_sum[QUATERNIONS_JPART][idx] = matrix1[QUATERNIONS_JPART][idx] + matrix2[QUATERNIONS_JPART][idx];
	            matrix_sum[QUATERNIONS_KPART][idx] = matrix1[QUATERNIONS_KPART][idx] + matrix2[QUATERNIONS_KPART][idx];
	    	}
		}
		
    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixOAdd(ityp **matrix1, ityp **matrix2, ityp **matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	register dim_typ idx;
	
	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
        {
	        idx = dim[COLUMNS]*i + j;
		    #pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
			{
		        matrix_sum[OCTONIONS_REALPART][idx] = matrix1[OCTONIONS_REALPART][idx] + matrix2[OCTONIONS_REALPART][idx];
		        matrix_sum[OCTONIONS_E1PART][idx] = matrix1[OCTONIONS_E1PART][idx] + matrix2[OCTONIONS_E1PART][idx];
		        matrix_sum[OCTONIONS_E2PART][idx] = matrix1[OCTONIONS_E2PART][idx] + matrix2[OCTONIONS_E2PART][idx];
		        matrix_sum[OCTONIONS_E3PART][idx] = matrix1[OCTONIONS_E3PART][idx] + matrix2[OCTONIONS_E3PART][idx];
		        matrix_sum[OCTONIONS_E4PART][idx] = matrix1[OCTONIONS_E4PART][idx] + matrix2[OCTONIONS_E4PART][idx];
		        matrix_sum[OCTONIONS_E5PART][idx] = matrix1[OCTONIONS_E5PART][idx] + matrix2[OCTONIONS_E5PART][idx];
		        matrix_sum[OCTONIONS_E6PART][idx] = matrix1[OCTONIONS_E6PART][idx] + matrix2[OCTONIONS_E6PART][idx];
		        matrix_sum[OCTONIONS_E7PART][idx] = matrix1[OCTONIONS_E7PART][idx] + matrix2[OCTONIONS_E7PART][idx];
			}
		}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixSAdd(ityp **matrix1, ityp **matrix2, ityp **matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	register dim_typ idx;

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
        {
        	idx = dim[COLUMNS]*i + j;
	        #pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
			{
				matrix_sum[SEDENIONS_REALPART][idx] = matrix1[SEDENIONS_REALPART][idx] + matrix2[SEDENIONS_REALPART][idx];
		        matrix_sum[SEDENIONS_E1PART][idx] = matrix1[SEDENIONS_E1PART][idx] + matrix2[SEDENIONS_E1PART][idx];
		        matrix_sum[SEDENIONS_E2PART][idx] = matrix1[SEDENIONS_E2PART][idx] + matrix2[SEDENIONS_E2PART][idx];
		        matrix_sum[SEDENIONS_E3PART][idx] = matrix1[SEDENIONS_E3PART][idx] + matrix2[SEDENIONS_E3PART][idx];
		        matrix_sum[SEDENIONS_E4PART][idx] = matrix1[SEDENIONS_E4PART][idx] + matrix2[SEDENIONS_E4PART][idx];
		        matrix_sum[SEDENIONS_E5PART][idx] = matrix1[SEDENIONS_E5PART][idx] + matrix2[SEDENIONS_E5PART][idx];
		        matrix_sum[SEDENIONS_E6PART][idx] = matrix1[SEDENIONS_E6PART][idx] + matrix2[SEDENIONS_E6PART][idx];
		        matrix_sum[SEDENIONS_E7PART][idx] = matrix1[SEDENIONS_E7PART][idx] + matrix2[SEDENIONS_E7PART][idx];
		        matrix_sum[SEDENIONS_E8PART][idx] = matrix1[SEDENIONS_E8PART][idx] + matrix2[SEDENIONS_E8PART][idx];
		        matrix_sum[SEDENIONS_E9PART][idx] = matrix1[SEDENIONS_E9PART][idx] + matrix2[SEDENIONS_E9PART][idx];
		        matrix_sum[SEDENIONS_E10PART][idx] = matrix1[SEDENIONS_E10PART][idx] + matrix2[SEDENIONS_E10PART][idx];
		        matrix_sum[SEDENIONS_E11PART][idx] = matrix1[SEDENIONS_E11PART][idx] + matrix2[SEDENIONS_E11PART][idx];
		        matrix_sum[SEDENIONS_E12PART][idx] = matrix1[SEDENIONS_E12PART][idx] + matrix2[SEDENIONS_E12PART][idx];
		        matrix_sum[SEDENIONS_E13PART][idx] = matrix1[SEDENIONS_E13PART][idx] + matrix2[SEDENIONS_E13PART][idx];
		        matrix_sum[SEDENIONS_E14PART][idx] = matrix1[SEDENIONS_E14PART][idx] + matrix2[SEDENIONS_E14PART][idx];
		        matrix_sum[SEDENIONS_E15PART][idx] = matrix1[SEDENIONS_E15PART][idx] + matrix2[SEDENIONS_E15PART][idx];
			}
		}

    return;
}

__MSNATIVE_ inline void __system __call_OSMM(ityp **matrix1, ityp **matrix2, ityp **matrix_product, const register dim_typ dim [static MAX_MATRICES], const register sel_typ algunits,
	void (* const prodFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_MATRICES]),
	void (* const sumFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_DIMENSIONS]),
	void (* const subFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_DIMENSIONS]))
{
	squareOSMM(prodFunc, matrix1, matrix2, matrix_product, dim[ROWS]);
	return;
}

__MSNATIVE_ inline void __system __call_STRASSENMM(ityp **matrix1, ityp **matrix2, ityp **matrix_product, const register dim_typ dim [static MAX_MATRICES], const register sel_typ algunits,
	void (* const prodFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_MATRICES]),
	void (* const sumFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_DIMENSIONS]),
	void (* const subFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_DIMENSIONS]))
{
	mmult_fast(dim[ROWS], algunits, matrix1, matrix2, matrix_product, prodFunc, sumFunc, subFunc);
	return;
}


__MSNATIVE_ inline void __system __call_NORMALMM(ityp **matrix1, ityp **matrix2, ityp **matrix_product, const register dim_typ dim [static MAX_MATRICES], const register sel_typ algunits,
	void (* const prodFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_MATRICES]),
	void (* const sumFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_DIMENSIONS]),
	void (* const subFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_DIMENSIONS]))
{
	prodFunc(matrix1, matrix2, matrix_product, dim);
	return;
}

__MSUTIL_ inline void __system __export squareOSMM(void (* const prodFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_MATRICES]), ityp **A, ityp **B, ityp **C, const register dim_typ dim)
{
	dim_typ i, j, k;
	const register dim_typ BLOCK_SIZE = access(curLayout)->block_size;
	for (i = 0; i < dim; i += BLOCK_SIZE)
		for (j = 0; j < dim; j += BLOCK_SIZE)
    		for (k = 0; k < dim; k += BLOCK_SIZE)
				/* Correct block dimensions if block "goes off edge of" the matrix */
				/* Perform individual block osmm */
				prodFunc(A + i + k*dim, B + k + j*dim, C + i + j*dim, (dim_typ3){min (BLOCK_SIZE, dim-i), min (BLOCK_SIZE, dim-j), min (BLOCK_SIZE, dim-k)});
}

// #define STRASSEN_BASECASEDIM

//
// Volker Strassen algorithm for matrix multiplication.
// Theoretical Runtime is O(N^2.807).
// Assume NxN matrices where N is a power of two.
// Algorithm:
//   Matrices X and Y are split into four smaller
//   (N/2)x(N/2) matrices as follows:
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
//       | P2 + P3                 (P0 + P4) - (P2 + P6)|
//        -                                            -
//
__MSNATIVE_ __MSUTIL_ void __system __export mmult_fast(const register dim_typ N, const register sel_typ algunits, ityp **X, ityp **Y,  ityp **Z,
	void (* const matrixMulFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_MATRICES]), 
	void (* const matrixAddFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_DIMENSIONS]),
	void (* const matrixSubFunc)(ityp **, ityp **, ityp **, const register dim_typ [static MAX_DIMENSIONS]))
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
		matrixMulFunc(X,Y,Z,(dim_typ3){N,N,N});
		return;
	}
	
	const register dim_typ n = N/2;      // size of sub-matrices
	
	ityp **A = X;    // A-D matrices embedded in X
	ityp **B = X + n;
	ityp **C = X + n*N;
	ityp **D = C + n;
	
	ityp **E = Y;    // E-H matrices embeded in Y
	ityp **F = Y + n;
	ityp **G = Y + n*N;
	ityp **H = G + n;
	
	ityp **P[7];   // allocate temp matrices off heap
	const register dim_typ sz = n*n*sizeof(ityp); // )n*n*sizeof(ityp*);
	const register dim_typ sz2 = sizeof(ityp*)*algunits;
	dim_typ i, j;
	//for (i = 0; i < 7; i++)
		//P[i] = malloc(sz);
		
	ityp **T = malloc(sz2); // malloc(sz);
	ityp **U = malloc(sz2); // malloc(sz);
	
	#pragma omp parallel for num_threads(algunits)
	for(i=0; i<algunits; ++i)
	{
		T[i] = malloc(sz);
		U[i] = malloc(sz);
		P[i] = malloc(sz2);
		for(j=6; j>=0; --j)
			P[i][j] = malloc(sz);
	}
	
	dim_typ2 dim =
	{
		n,
		n
	};
	
	// P0 = A*(F - H);
	matrixSubFunc(F,H,T,dim);
	mmult_fast(n,algunits,A,T,P[0],matrixMulFunc,matrixAddFunc,matrixSubFunc);
	
	// P1 = (A + B)*H
	matrixAddFunc(A,B,T,dim);
	mmult_fast(n,algunits,T,H,P[1],matrixMulFunc,matrixAddFunc,matrixSubFunc);
	
	// P2 = (C + D)*E
	matrixAddFunc(C,D,T,dim);
	mmult_fast(n, algunits, T, E, P[2],matrixMulFunc,matrixAddFunc,matrixSubFunc);
	
	// P3 = D*(G - E);
	matrixSubFunc(G,E,T,dim);
	mmult_fast(n,algunits,D,T,P[3],matrixMulFunc,matrixAddFunc,matrixSubFunc);
	
	// P4 = (A + D)*(E + H)
	matrixAddFunc(A,D,T,dim);
	matrixAddFunc(E,H,U,dim);
	mmult_fast(n,algunits,T,U,P[4],matrixMulFunc,matrixAddFunc,matrixSubFunc);
	
	// P5 = (B - D)*(G + H)
	matrixSubFunc(B,D,T,dim);
	matrixAddFunc(G,H,U,dim);
	mmult_fast(n,algunits,T,U,P[5],matrixMulFunc,matrixAddFunc,matrixSubFunc);
	
	// P6 = (A - C)*(E + F)
	matrixSubFunc(A,C,T,dim);
	matrixAddFunc(E,F,U,dim);
	mmult_fast(n,algunits,T,U,P[6],matrixMulFunc,matrixAddFunc,matrixSubFunc);
	
	// Z upper left = (P3 + P4) + (P5 - P1)
	matrixAddFunc(P[4], P[3], T, dim);
	matrixSubFunc(P[5], P[1], U, dim);
	matrixAddFunc(T,U,Z,dim);
	
	// Z lower left = P2 + P3
	matrixAddFunc(P[2],P[3], Z + n*N, dim);
	
	// Z upper right = P0 + P1
	matrixAddFunc(P[0], P[1], Z + n, dim);
	
	// Z lower right = (P0 + P4) - (P2 + P6)
	matrixAddFunc(P[0], P[4], T, dim);
	matrixAddFunc(P[2], P[6], U, dim);
	matrixSubFunc(T,U,Z+n*(N+1), dim);
	

	#pragma omp parallel for num_threads(algunits)
	for(i=0; i<algunits; ++i)
	{
		free(U[i]);
		free(T[i]);
		for (j = 6; j >= 0; --j)
			free(P[i][j]);
		free(P[i]);
	}
	
	free(U);  // deallocate temp matrices
	free(T);
	
	return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixMultiplication(ityp **matrix1, ityp **matrix2, ityp **matrix_product, const register dim_typ dim[static MAX_MATRICES]) // dim_typ righe, dim_typ colonne, dim_typ colonne2)
{
    dim_typ i, j, k;

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(k=0; k<dim[COLUMNS]; ++k)
        	#pragma omp parallel for
       		for(j=0; j<dim[COLUMNS2]; ++j)
            	*((*matrix_product) + dim[COLUMNS2]*i + j) += (*((*matrix1) + dim[COLUMNS]*i + k) * *((*matrix2) + dim[COLUMNS2]*k + j));

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixCMultiplication(ityp **matrix1, ityp **matrix2, ityp **matrix_product, const register dim_typ dim[static MAX_MATRICES])
{
    dim_typ i, j, k;
    register dim_typ idx[MAX_MATRICES];
    
	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(k=0; k<dim[COLUMNS]; ++k)
    		#pragma omp parallel for
        	for(j=0; j<dim[COLUMNS2]; ++j)
        	{
        		idx[FIRST_MATRIX] = dim[COLUMNS]*i + k;
            	idx[SECOND_MATRIX] = dim[COLUMNS2]*k + j;
            	idx[MATRIX_PRODUCT] = dim[COLUMNS2]*i + j;
            	const register double complex tmpres = ((matrix1[REAL_PART][idx[FIRST_MATRIX]] + matrix1[IMAG_PART][idx[FIRST_MATRIX]]*I) * (matrix2[REAL_PART][idx[SECOND_MATRIX]] + matrix2[IMAG_PART][idx[SECOND_MATRIX]]*I));
				matrix_product[REAL_PART][idx[MATRIX_PRODUCT]] += creal(tmpres);
	            matrix_product[IMAG_PART][idx[MATRIX_PRODUCT]] += cimag(tmpres);
	       	}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixQMultiplication(ityp **matrix1, ityp **matrix2, ityp **matrix_product, const register dim_typ dim[static MAX_MATRICES])
{
    dim_typ i, j, k;
    register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(k=0; k<dim[COLUMNS]; ++k)
    		#pragma omp parallel for
        	for(j=0; j<dim[COLUMNS2]; ++j)
        	{
        		idx[FIRST_MATRIX] = dim[COLUMNS]*i + k;
            	idx[SECOND_MATRIX] = dim[COLUMNS2]*k + j;
            	idx[MATRIX_PRODUCT] = dim[COLUMNS2]*i + j;
            	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
            	{
            		
            		
	                matrix_product[QUATERNIONS_REALPART][idx[MATRIX_PRODUCT]] += matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]] -
	                                                              matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]] -
	                                                              matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]] -
	                                                              matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]];
	                matrix_product[QUATERNIONS_IPART][idx[MATRIX_PRODUCT]] += matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]] +
	                                                        	matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]] +
	                                                           	matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]] -
	                                                           	matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]];
	                matrix_product[QUATERNIONS_JPART][idx[MATRIX_PRODUCT]] += matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]] +
	                                                           matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]] +
	                                                           matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]] -
	                                                           matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]];
	                matrix_product[QUATERNIONS_KPART][idx[MATRIX_PRODUCT]] += matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]] +
	                                                           matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]] +
	                                                           matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]] -
	                                                           matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]];
	            }
	    	}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixOMultiplication(ityp **matrix1, ityp **matrix2, ityp **matrix_product, const register dim_typ dim[static MAX_MATRICES])
{
    dim_typ i, j, k;
    register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
   		#pragma omp parallel for
        for(k=0; k<dim[COLUMNS]; ++k)
    		#pragma omp parallel for
        	for(j=0; j<dim[COLUMNS2]; ++j)
        	{
				idx[FIRST_MATRIX] = dim[COLUMNS]*i + k;
				idx[SECOND_MATRIX] = dim[COLUMNS2]*k + j;
				idx[MATRIX_PRODUCT] = dim[COLUMNS2]*i + j;
	        	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
				{
					matrix_product[OCTONIONS_REALPART][idx[MATRIX_PRODUCT]] +=  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] -
			                                                             matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]];
			
			        matrix_product[OCTONIONS_E1PART][idx[MATRIX_PRODUCT]] +=  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]];
			
			        matrix_product[OCTONIONS_E2PART][idx[MATRIX_PRODUCT]] +=  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]];
			
			        matrix_product[OCTONIONS_E3PART][idx[MATRIX_PRODUCT]] +=  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]];
			
			        matrix_product[OCTONIONS_E4PART][idx[MATRIX_PRODUCT]] +=  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]];
			
			        matrix_product[OCTONIONS_E5PART][idx[MATRIX_PRODUCT]] +=  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]];
			
			        matrix_product[OCTONIONS_E6PART][idx[MATRIX_PRODUCT]] +=  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]];
			
			        matrix_product[OCTONIONS_E7PART][idx[MATRIX_PRODUCT]] +=  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]];
				}
			}
            	
    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixSMultiplication(ityp **matrix1, ityp **matrix2, ityp **matrix_product, const register dim_typ dim[static MAX_ABSTRACT_DIMENSIONS])
{
    dim_typ i, j, k;
	register dim_typ idx[MAX_MATRICES];

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(k=0; k<dim[COLUMNS]; ++k)
			#pragma omp parallel for
        	for(j=0; j<dim[COLUMNS2]; ++j)
        	{
				idx[FIRST_MATRIX] = dim[COLUMNS]*i + k;
				idx[SECOND_MATRIX] = dim[COLUMNS2]*k + j;
				idx[MATRIX_PRODUCT] = dim[COLUMNS2]*i + j;
        		#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
				{
					matrix_product[SEDENIONS_REALPART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                             matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E1PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E2PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E3PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E4PART][idx[MATRIX_PRODUCT]] +=   matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E5PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]];
			                                                   
			        matrix_product[SEDENIONS_E6PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]];
			 
			        matrix_product[SEDENIONS_E7PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E8PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E9PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                   matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
			                                                   matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E10PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E11PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E12PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E13PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E14PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]];
			
			        matrix_product[SEDENIONS_E15PART][idx[MATRIX_PRODUCT]] +=  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
			                                                    matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
			                                                    matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]];
				}
			}
		
    return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixKProduct(ityp **matrix1, ityp **matrix2, ityp **matrix_product, register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS])
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
                	*((*matrix_product) + pdim*(k+(i*dim[SECOND_MATRIX][ROWS])) + l+(j*dim[SECOND_MATRIX][COLUMNS])) = (*((*matrix1) + dim[FIRST_MATRIX][COLUMNS]*i + j) * *((*matrix2) + dim[SECOND_MATRIX][COLUMNS]*k + l));

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixKCProduct(ityp **matrix1, ityp **matrix2, ityp **matrix_product, register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS])
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
                	const register double complex tmpres = ((matrix1[REAL_PART][idx[FIRST_MATRIX]] + matrix1[IMAG_PART][idx[FIRST_MATRIX]]*I)*(matrix2[REAL_PART][idx[SECOND_MATRIX]] + matrix2[IMAG_PART][idx[SECOND_MATRIX]]*I));
                    matrix_product[REAL_PART][idx[MATRIX_PRODUCT]] = creal(tmpres);
                    matrix_product[IMAG_PART][idx[MATRIX_PRODUCT]] = cimag(tmpres);
	        	}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixKQProduct(ityp **matrix1, ityp **matrix2, ityp **matrix_product, register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS])
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
                	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
                	{
	                    
	                    matrix_product[QUATERNIONS_REALPART][idx[MATRIX_PRODUCT]] = matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]] -
	                                                                                                                               matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]] -
	                                                                                                                               matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]] -
	                                                                                                                               matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]];

	                    matrix_product[QUATERNIONS_IPART][idx[MATRIX_PRODUCT]] = matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]] +
	                                                                                                                            matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]] +
	                                                                                                                            matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]] -
	                                                                                                                            matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]];

	                    matrix_product[QUATERNIONS_JPART][idx[MATRIX_PRODUCT]] = matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]] +
	                                                                                                                            matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]] +
	                                                                                                                            matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]] -
	                                                                                                                            matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]];

	                    matrix_product[QUATERNIONS_KPART][idx[MATRIX_PRODUCT]] = matrix1[QUATERNIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_KPART][idx[SECOND_MATRIX]] +
	                                                                                                                            matrix1[QUATERNIONS_KPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_REALPART][idx[SECOND_MATRIX]] +
	                                                                                                                            matrix1[QUATERNIONS_IPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_JPART][idx[SECOND_MATRIX]] -
	                                                                                                                            matrix1[QUATERNIONS_JPART][idx[FIRST_MATRIX]] * matrix2[QUATERNIONS_IPART][idx[SECOND_MATRIX]];
	    			}
	    		}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixKOProduct(ityp **matrix1, ityp **matrix2, ityp **matrix_product, register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS])
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
                	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
					{
				        matrix_product[OCTONIONS_REALPART][idx[MATRIX_PRODUCT]] = matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]];
				
				        matrix_product[OCTONIONS_E1PART][idx[MATRIX_PRODUCT]] =  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]];
				
				        matrix_product[OCTONIONS_E2PART][idx[MATRIX_PRODUCT]] =  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]];
				
				        matrix_product[OCTONIONS_E3PART][idx[MATRIX_PRODUCT]] =  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]];
				
				        matrix_product[OCTONIONS_E4PART][idx[MATRIX_PRODUCT]] =  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]];
				
				        matrix_product[OCTONIONS_E5PART][idx[MATRIX_PRODUCT]] =  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]];
				
				        matrix_product[OCTONIONS_E6PART][idx[MATRIX_PRODUCT]] =  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]];
				
				        matrix_product[OCTONIONS_E7PART][idx[MATRIX_PRODUCT]] =  matrix1[OCTONIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E6PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E4PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[OCTONIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_E2PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[OCTONIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[OCTONIONS_REALPART][idx[SECOND_MATRIX]];
					}
				}
				
    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixKSProduct(ityp **matrix1, ityp **matrix2, ityp **matrix_product, register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS])
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
                	#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
					{
				        matrix_product[SEDENIONS_REALPART][idx[MATRIX_PRODUCT]] = matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                   matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E1PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E2PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E3PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E4PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E5PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E6PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E7PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E8PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E9PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                 matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
				                                                                                 matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E10PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E11PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E12PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E13PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E14PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E7PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_REALPART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]];
				
				        matrix_product[SEDENIONS_E15PART][idx[MATRIX_PRODUCT]] =  matrix1[SEDENIONS_REALPART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E15PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E1PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E2PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E13PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E3PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E12PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E4PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E11PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E5PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E10PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E6PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E9PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E7PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E8PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E8PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E9PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E6PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E10PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E5PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E11PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E4PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E12PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E3PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E13PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E2PART][idx[SECOND_MATRIX]] -
				                                                                                  matrix1[SEDENIONS_E14PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E1PART][idx[SECOND_MATRIX]] +
				                                                                                  matrix1[SEDENIONS_E15PART][idx[FIRST_MATRIX]] * matrix2[SEDENIONS_E14PART][idx[SECOND_MATRIX]];
					}
				}
    return;
}

// ENDAlgebra

// END
