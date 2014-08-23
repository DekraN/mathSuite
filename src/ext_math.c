// ext_math.h
// as External Math Library
// 23/08/2014 Marco Chiarelli aka DekraN
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
            SCHAR_MAX,
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
            ULLONG_MAX
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

__MSNATIVE_ inline void __export _complexSum(ityp **complex, ityp complexRes[static MAX_COMPLEX_UNITS])
{
	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
	{
		complexRes[REAL_PART] = complex[FIRST_NUMBER][REAL_PART] + complex[SECOND_NUMBER][REAL_PART];
    	complexRes[IMAG_PART] = complex[FIRST_NUMBER][IMAG_PART] + complex[SECOND_NUMBER][IMAG_PART];
	}
    return;
}

__MSNATIVE_ inline void __export _complexProd(ityp **complex, ityp complexRes[static MAX_COMPLEX_UNITS])
{
	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
	{
		complexRes[REAL_PART] =((complex[FIRST_NUMBER][REAL_PART]*complex[SECOND_NUMBER][REAL_PART])-(complex[FIRST_NUMBER][IMAG_PART]*complex[SECOND_NUMBER][IMAG_PART]));
	    complexRes[IMAG_PART] =((complex[FIRST_NUMBER][REAL_PART]*complex[SECOND_NUMBER][IMAG_PART])+(complex[FIRST_NUMBER][IMAG_PART]*complex[SECOND_NUMBER][REAL_PART]));
	}
    return;
}

__MSNATIVE_ inline void __export _quaternionsSum(ityp **quaternions, ityp quaternionsRes[static MAX_QUATERNIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
	{
		quaternionsRes[QUATERNIONS_REALPART] = quaternions[FIRST_NUMBER][REAL_PART] + quaternions[SECOND_NUMBER][REAL_PART];
	    quaternionsRes[QUATERNIONS_IPART] = quaternions[FIRST_NUMBER][QUATERNIONS_IPART] + quaternions[SECOND_NUMBER][QUATERNIONS_IPART];
	    quaternionsRes[QUATERNIONS_JPART] = quaternions[FIRST_NUMBER][QUATERNIONS_JPART] + quaternions[SECOND_NUMBER][QUATERNIONS_JPART];
	    quaternionsRes[QUATERNIONS_KPART] = quaternions[FIRST_NUMBER][QUATERNIONS_KPART] + quaternions[SECOND_NUMBER][QUATERNIONS_KPART];
	}
    return;
}

__MSNATIVE_ inline void __export _quaternionsProduct(ityp **quaternions, ityp quaternionsRes[static MAX_QUATERNIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
	{
		quaternionsRes[QUATERNIONS_REALPART] = (quaternions[FIRST_NUMBER][QUATERNIONS_REALPART]*quaternions[SECOND_NUMBER][QUATERNIONS_REALPART]) -
	                                            (quaternions[FIRST_NUMBER][QUATERNIONS_IPART]*quaternions[SECOND_NUMBER][QUATERNIONS_IPART]) -
	                                            (quaternions[FIRST_NUMBER][QUATERNIONS_JPART]*quaternions[SECOND_NUMBER][QUATERNIONS_JPART]) -
	                                            (quaternions[FIRST_NUMBER][QUATERNIONS_KPART]*quaternions[SECOND_NUMBER][QUATERNIONS_KPART]);
	    quaternionsRes[QUATERNIONS_IPART] = (quaternions[FIRST_NUMBER][QUATERNIONS_REALPART]*quaternions[SECOND_NUMBER][QUATERNIONS_IPART]) +
	                                        (quaternions[FIRST_NUMBER][QUATERNIONS_IPART]*quaternions[SECOND_NUMBER][QUATERNIONS_REALPART]) +
	                                        (quaternions[FIRST_NUMBER][QUATERNIONS_JPART]*quaternions[SECOND_NUMBER][QUATERNIONS_KPART]) -
	                                        (quaternions[FIRST_NUMBER][QUATERNIONS_KPART]*quaternions[SECOND_NUMBER][QUATERNIONS_JPART]);
	    quaternionsRes[QUATERNIONS_JPART] = (quaternions[FIRST_NUMBER][QUATERNIONS_REALPART]*quaternions[SECOND_NUMBER][QUATERNIONS_JPART]) +
	                                        (quaternions[FIRST_NUMBER][QUATERNIONS_JPART]*quaternions[SECOND_NUMBER][QUATERNIONS_REALPART]) +
	                                        (quaternions[FIRST_NUMBER][QUATERNIONS_KPART]*quaternions[SECOND_NUMBER][QUATERNIONS_IPART]) -
	                                        (quaternions[FIRST_NUMBER][QUATERNIONS_IPART]*quaternions[SECOND_NUMBER][QUATERNIONS_KPART]);
	    quaternionsRes[QUATERNIONS_KPART] = (quaternions[FIRST_NUMBER][QUATERNIONS_REALPART]*quaternions[SECOND_NUMBER][QUATERNIONS_KPART]) +
	                                        (quaternions[FIRST_NUMBER][QUATERNIONS_KPART]*quaternions[SECOND_NUMBER][QUATERNIONS_REALPART]) +
	                                        (quaternions[FIRST_NUMBER][QUATERNIONS_IPART]*quaternions[SECOND_NUMBER][QUATERNIONS_JPART]) -
	                                        (quaternions[FIRST_NUMBER][QUATERNIONS_JPART]*quaternions[SECOND_NUMBER][QUATERNIONS_IPART]);
	}

    return;
}

__MSNATIVE_ void __export _octonionsEMTSum(ityp **octonions, ityp octonionsRes[static MAX_OCTONIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
	{
		octonionsRes[OCTONIONS_REALPART] = octonions[FIRST_NUMBER][OCTONIONS_REALPART] + octonions[SECOND_NUMBER][OCTONIONS_REALPART];
	    octonionsRes[OCTONIONS_E1PART] = octonions[FIRST_NUMBER][OCTONIONS_E1PART] + octonions[SECOND_NUMBER][OCTONIONS_E1PART];
	    octonionsRes[OCTONIONS_E2PART] = octonions[FIRST_NUMBER][OCTONIONS_E2PART] + octonions[SECOND_NUMBER][OCTONIONS_E2PART];
	    octonionsRes[OCTONIONS_E3PART] = octonions[FIRST_NUMBER][OCTONIONS_E3PART] + octonions[SECOND_NUMBER][OCTONIONS_E3PART];
	    octonionsRes[OCTONIONS_E4PART] = octonions[FIRST_NUMBER][OCTONIONS_E4PART] + octonions[SECOND_NUMBER][OCTONIONS_E4PART];
	    octonionsRes[OCTONIONS_E5PART] = octonions[FIRST_NUMBER][OCTONIONS_E5PART] + octonions[SECOND_NUMBER][OCTONIONS_E5PART];
	    octonionsRes[OCTONIONS_E6PART] = octonions[FIRST_NUMBER][OCTONIONS_E6PART] + octonions[SECOND_NUMBER][OCTONIONS_E6PART];
	    octonionsRes[OCTONIONS_E7PART] = octonions[FIRST_NUMBER][OCTONIONS_E7PART] + octonions[SECOND_NUMBER][OCTONIONS_E7PART];
	}
	return;
}

__MSNATIVE_ void __export _octonionsMTSum(ityp **octonions, ityp octonionsRes[static MAX_OCTONIONS_UNITS])
{
	#pragma omp parallel
	{
		octonionsRes[OCTONIONS_REALPART] = octonions[FIRST_NUMBER][OCTONIONS_REALPART] + octonions[SECOND_NUMBER][OCTONIONS_REALPART];
	    octonionsRes[OCTONIONS_E1PART] = octonions[FIRST_NUMBER][OCTONIONS_E1PART] + octonions[SECOND_NUMBER][OCTONIONS_E1PART];
	    octonionsRes[OCTONIONS_E2PART] = octonions[FIRST_NUMBER][OCTONIONS_E2PART] + octonions[SECOND_NUMBER][OCTONIONS_E2PART];
	    octonionsRes[OCTONIONS_E3PART] = octonions[FIRST_NUMBER][OCTONIONS_E3PART] + octonions[SECOND_NUMBER][OCTONIONS_E3PART];
	    octonionsRes[OCTONIONS_E4PART] = octonions[FIRST_NUMBER][OCTONIONS_E4PART] + octonions[SECOND_NUMBER][OCTONIONS_E4PART];
	    octonionsRes[OCTONIONS_E5PART] = octonions[FIRST_NUMBER][OCTONIONS_E5PART] + octonions[SECOND_NUMBER][OCTONIONS_E5PART];
	    octonionsRes[OCTONIONS_E6PART] = octonions[FIRST_NUMBER][OCTONIONS_E6PART] + octonions[SECOND_NUMBER][OCTONIONS_E6PART];
	    octonionsRes[OCTONIONS_E7PART] = octonions[FIRST_NUMBER][OCTONIONS_E7PART] + octonions[SECOND_NUMBER][OCTONIONS_E7PART];
	}
	return;
}

__MSNATIVE_ inline void __export _octonionsSum(ityp **octonions, ityp octonionsRes[static MAX_OCTONIONS_UNITS])
{
	void (* const omp_func)(ityp **, ityp [static MAX_OCTONIONS_UNITS]) = isSett(BOOLS_EXTENSIVEMULTITHREADING) ? _octonionsEMTSum : _octonionsMTSum;
	omp_func(octonions, octonionsRes);
    return;
}

__MSNATIVE_ void __export _octonionsEMTProduct(ityp **octonions, ityp octonionsRes[static MAX_OCTONIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
	{
		octonionsRes[OCTONIONS_REALPART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]);

	    octonionsRes[OCTONIONS_E1PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]);


	    octonionsRes[OCTONIONS_E2PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]);

	    octonionsRes[OCTONIONS_E3PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]);

	    octonionsRes[OCTONIONS_E4PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]);

	    octonionsRes[OCTONIONS_E5PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]);

	    octonionsRes[OCTONIONS_E6PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]);

	    octonionsRes[OCTONIONS_E7PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]);
	}
	return;
}

__MSNATIVE_ void __export _octonionsMTProduct(ityp **octonions, ityp octonionsRes[static MAX_OCTONIONS_UNITS])
{
	#pragma omp parallel
	{
		octonionsRes[OCTONIONS_REALPART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) -
                                       (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]);

	    octonionsRes[OCTONIONS_E1PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]);


	    octonionsRes[OCTONIONS_E2PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]);

	    octonionsRes[OCTONIONS_E3PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]);

	    octonionsRes[OCTONIONS_E4PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]);

	    octonionsRes[OCTONIONS_E5PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]);

	    octonionsRes[OCTONIONS_E6PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]);

	    octonionsRes[OCTONIONS_E7PART] = (octonions[FIRST_NUMBER][OCTONIONS_REALPART]*octonions[SECOND_NUMBER][OCTONIONS_E7PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E1PART]*octonions[SECOND_NUMBER][OCTONIONS_E3PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E2PART]*octonions[SECOND_NUMBER][OCTONIONS_E6PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E3PART]*octonions[SECOND_NUMBER][OCTONIONS_E1PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E4PART]*octonions[SECOND_NUMBER][OCTONIONS_E5PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E5PART]*octonions[SECOND_NUMBER][OCTONIONS_E4PART]) -
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E6PART]*octonions[SECOND_NUMBER][OCTONIONS_E2PART]) +
	                                     (octonions[FIRST_NUMBER][OCTONIONS_E7PART]*octonions[SECOND_NUMBER][OCTONIONS_REALPART]);
	}
	return;
}

__MSNATIVE_ inline void __export _octonionsProduct(ityp **octonions, ityp octonionsRes[static MAX_OCTONIONS_UNITS])
{
	void (* const omp_func)(ityp **, ityp [static MAX_OCTONIONS_UNITS]) = isSett(BOOLS_EXTENSIVEMULTITHREADING) ? _octonionsEMTProduct : _octonionsMTProduct;
	omp_func(octonions, octonionsRes);
    return;
}

__MSNATIVE_ void __export _sedenionsEMTSum(ityp **sedenions, ityp sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
	{
		sedenionsRes[SEDENIONS_REALPART] = sedenions[FIRST_NUMBER][SEDENIONS_REALPART] + sedenions[SECOND_NUMBER][SEDENIONS_REALPART];
	    sedenionsRes[SEDENIONS_E1PART] = sedenions[FIRST_NUMBER][SEDENIONS_E1PART] + sedenions[SECOND_NUMBER][SEDENIONS_E1PART];
	    sedenionsRes[SEDENIONS_E2PART] = sedenions[FIRST_NUMBER][SEDENIONS_E2PART] + sedenions[SECOND_NUMBER][SEDENIONS_E2PART];
	    sedenionsRes[SEDENIONS_E3PART] = sedenions[FIRST_NUMBER][SEDENIONS_E3PART] + sedenions[SECOND_NUMBER][SEDENIONS_E3PART];
	    sedenionsRes[SEDENIONS_E4PART] = sedenions[FIRST_NUMBER][SEDENIONS_E4PART] + sedenions[SECOND_NUMBER][SEDENIONS_E4PART];
	    sedenionsRes[SEDENIONS_E5PART] = sedenions[FIRST_NUMBER][SEDENIONS_E5PART] + sedenions[SECOND_NUMBER][SEDENIONS_E5PART];
	    sedenionsRes[SEDENIONS_E6PART] = sedenions[FIRST_NUMBER][SEDENIONS_E6PART] + sedenions[SECOND_NUMBER][SEDENIONS_E6PART];
	    sedenionsRes[SEDENIONS_E7PART] = sedenions[FIRST_NUMBER][SEDENIONS_E7PART] + sedenions[SECOND_NUMBER][SEDENIONS_E7PART];
	    sedenionsRes[SEDENIONS_E8PART] = sedenions[FIRST_NUMBER][SEDENIONS_E8PART] + sedenions[SECOND_NUMBER][SEDENIONS_E8PART];
	    sedenionsRes[SEDENIONS_E9PART] = sedenions[FIRST_NUMBER][SEDENIONS_E9PART] + sedenions[SECOND_NUMBER][SEDENIONS_E9PART];
	    sedenionsRes[SEDENIONS_E10PART] = sedenions[FIRST_NUMBER][SEDENIONS_E10PART] + sedenions[SECOND_NUMBER][SEDENIONS_E10PART];
	    sedenionsRes[SEDENIONS_E11PART] = sedenions[FIRST_NUMBER][SEDENIONS_E11PART] + sedenions[SECOND_NUMBER][SEDENIONS_E11PART];
	    sedenionsRes[SEDENIONS_E12PART] = sedenions[FIRST_NUMBER][SEDENIONS_E12PART] + sedenions[SECOND_NUMBER][SEDENIONS_E12PART];
	    sedenionsRes[SEDENIONS_E13PART] = sedenions[FIRST_NUMBER][SEDENIONS_E13PART] + sedenions[SECOND_NUMBER][SEDENIONS_E13PART];
	    sedenionsRes[SEDENIONS_E14PART] = sedenions[FIRST_NUMBER][SEDENIONS_E14PART] + sedenions[SECOND_NUMBER][SEDENIONS_E14PART];
	    sedenionsRes[SEDENIONS_E15PART] = sedenions[FIRST_NUMBER][SEDENIONS_E15PART] + sedenions[SECOND_NUMBER][SEDENIONS_E15PART];
	}
	return;
}

__MSNATIVE_ void __export _sedenionsMTSum(ityp **sedenions, ityp sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	#pragma omp parallel
	{
		sedenionsRes[SEDENIONS_REALPART] = sedenions[FIRST_NUMBER][SEDENIONS_REALPART] + sedenions[SECOND_NUMBER][SEDENIONS_REALPART];
	    sedenionsRes[SEDENIONS_E1PART] = sedenions[FIRST_NUMBER][SEDENIONS_E1PART] + sedenions[SECOND_NUMBER][SEDENIONS_E1PART];
	    sedenionsRes[SEDENIONS_E2PART] = sedenions[FIRST_NUMBER][SEDENIONS_E2PART] + sedenions[SECOND_NUMBER][SEDENIONS_E2PART];
	    sedenionsRes[SEDENIONS_E3PART] = sedenions[FIRST_NUMBER][SEDENIONS_E3PART] + sedenions[SECOND_NUMBER][SEDENIONS_E3PART];
	    sedenionsRes[SEDENIONS_E4PART] = sedenions[FIRST_NUMBER][SEDENIONS_E4PART] + sedenions[SECOND_NUMBER][SEDENIONS_E4PART];
	    sedenionsRes[SEDENIONS_E5PART] = sedenions[FIRST_NUMBER][SEDENIONS_E5PART] + sedenions[SECOND_NUMBER][SEDENIONS_E5PART];
	    sedenionsRes[SEDENIONS_E6PART] = sedenions[FIRST_NUMBER][SEDENIONS_E6PART] + sedenions[SECOND_NUMBER][SEDENIONS_E6PART];
	    sedenionsRes[SEDENIONS_E7PART] = sedenions[FIRST_NUMBER][SEDENIONS_E7PART] + sedenions[SECOND_NUMBER][SEDENIONS_E7PART];
	    sedenionsRes[SEDENIONS_E8PART] = sedenions[FIRST_NUMBER][SEDENIONS_E8PART] + sedenions[SECOND_NUMBER][SEDENIONS_E8PART];
	    sedenionsRes[SEDENIONS_E9PART] = sedenions[FIRST_NUMBER][SEDENIONS_E9PART] + sedenions[SECOND_NUMBER][SEDENIONS_E9PART];
	    sedenionsRes[SEDENIONS_E10PART] = sedenions[FIRST_NUMBER][SEDENIONS_E10PART] + sedenions[SECOND_NUMBER][SEDENIONS_E10PART];
	    sedenionsRes[SEDENIONS_E11PART] = sedenions[FIRST_NUMBER][SEDENIONS_E11PART] + sedenions[SECOND_NUMBER][SEDENIONS_E11PART];
	    sedenionsRes[SEDENIONS_E12PART] = sedenions[FIRST_NUMBER][SEDENIONS_E12PART] + sedenions[SECOND_NUMBER][SEDENIONS_E12PART];
	    sedenionsRes[SEDENIONS_E13PART] = sedenions[FIRST_NUMBER][SEDENIONS_E13PART] + sedenions[SECOND_NUMBER][SEDENIONS_E13PART];
	    sedenionsRes[SEDENIONS_E14PART] = sedenions[FIRST_NUMBER][SEDENIONS_E14PART] + sedenions[SECOND_NUMBER][SEDENIONS_E14PART];
	    sedenionsRes[SEDENIONS_E15PART] = sedenions[FIRST_NUMBER][SEDENIONS_E15PART] + sedenions[SECOND_NUMBER][SEDENIONS_E15PART];
	}
	return;
}

__MSNATIVE_ inline void __export _sedenionsSum(ityp **sedenions, ityp sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	void (* const omp_func)(ityp **, ityp [static MAX_SEDENIONS_UNITS]) = isSett(BOOLS_EXTENSIVEMULTITHREADING) ? _sedenionsEMTSum : _sedenionsMTSum;
	omp_func(sedenions, sedenionsRes);
    return;
}

__MSNATIVE_ void __export _sedenionsEMTProduct(ityp **sedenions, ityp sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
	{
		sedenionsRes[SEDENIONS_REALPART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]);

	    sedenionsRes[SEDENIONS_E1PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]);

	    sedenionsRes[SEDENIONS_E2PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]);

	    sedenionsRes[SEDENIONS_E3PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]);

	    sedenionsRes[SEDENIONS_E4PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]);

	    sedenionsRes[SEDENIONS_E5PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]);

	    sedenionsRes[SEDENIONS_E6PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]);



	    sedenionsRes[SEDENIONS_E7PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]);


	    sedenionsRes[SEDENIONS_E8PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]);

	    sedenionsRes[SEDENIONS_E9PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]);

	    sedenionsRes[SEDENIONS_E10PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]);

	    sedenionsRes[SEDENIONS_E11PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]);

	    sedenionsRes[SEDENIONS_E12PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]);

	    sedenionsRes[SEDENIONS_E13PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]);

	    sedenionsRes[SEDENIONS_E14PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]);

	    sedenionsRes[SEDENIONS_E15PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]);
	}

	return;
}

__MSNATIVE_ void __export _sedenionsMTProduct(ityp **sedenions, ityp sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	#pragma omp parallel
	{
		sedenionsRes[SEDENIONS_REALPART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
                                       (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]);

	    sedenionsRes[SEDENIONS_E1PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]);

	    sedenionsRes[SEDENIONS_E2PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]);

	    sedenionsRes[SEDENIONS_E3PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]);

	    sedenionsRes[SEDENIONS_E4PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]);

	    sedenionsRes[SEDENIONS_E5PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]);

	    sedenionsRes[SEDENIONS_E6PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]);



	    sedenionsRes[SEDENIONS_E7PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]);


	    sedenionsRes[SEDENIONS_E8PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]);

	    sedenionsRes[SEDENIONS_E9PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                     (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]);

	    sedenionsRes[SEDENIONS_E10PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]);

	    sedenionsRes[SEDENIONS_E11PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]);

	    sedenionsRes[SEDENIONS_E12PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]);

	    sedenionsRes[SEDENIONS_E13PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]);

	    sedenionsRes[SEDENIONS_E14PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E7PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]);

	    sedenionsRes[SEDENIONS_E15PART] = (sedenions[FIRST_NUMBER][SEDENIONS_REALPART]*sedenions[SECOND_NUMBER][SEDENIONS_E15PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E1PART]*sedenions[SECOND_NUMBER][SEDENIONS_E14PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E2PART]*sedenions[SECOND_NUMBER][SEDENIONS_E13PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E3PART]*sedenions[SECOND_NUMBER][SEDENIONS_E12PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E4PART]*sedenions[SECOND_NUMBER][SEDENIONS_E11PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E5PART]*sedenions[SECOND_NUMBER][SEDENIONS_E10PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E6PART]*sedenions[SECOND_NUMBER][SEDENIONS_E9PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E7PART]*sedenions[SECOND_NUMBER][SEDENIONS_E8PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E8PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E9PART]*sedenions[SECOND_NUMBER][SEDENIONS_E6PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E10PART]*sedenions[SECOND_NUMBER][SEDENIONS_E5PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E11PART]*sedenions[SECOND_NUMBER][SEDENIONS_E4PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E12PART]*sedenions[SECOND_NUMBER][SEDENIONS_E3PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E13PART]*sedenions[SECOND_NUMBER][SEDENIONS_E2PART]) -
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E14PART]*sedenions[SECOND_NUMBER][SEDENIONS_E1PART]) +
	                                      (sedenions[FIRST_NUMBER][SEDENIONS_E15PART]*sedenions[SECOND_NUMBER][SEDENIONS_REALPART]);
	}

	return;
}

__MSNATIVE_ inline void __export _sedenionsProduct(ityp **sedenions, ityp sedenionsRes[static MAX_SEDENIONS_UNITS])
{
	void (* const omp_func)(ityp **, ityp [static MAX_SEDENIONS_UNITS]) = isSett(BOOLS_EXTENSIVEMULTITHREADING) ? _sedenionsEMTProduct : _sedenionsMTProduct;
	omp_func(sedenions, sedenionsRes);
    return;
}

__MSNATIVE_ bool __export _secondGradeEquationSolver(ityp *abc, ityp root[static MAX_DIMENSIONS])
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

// char * getDayName(sel_typ day);
// sel_typ getDayNumber(sel_typ dd, sel_typ mm, uint64_t yy);
//returns the name of the day
/// thanks to Bibek Subedi for these two functions,
/// which I renamed, modified and adapted to this program. Link at:
/// http://www.programming-technique.blogspot.it/2013/02/cc-program-to-ditermine-day-of-week.html
__MSUTIL_ inline const char * const __export getDayName(sel_typ dd)
{
    return ext_math.days_week_names[dd];
}

/// by Me (DekraN)
__MSUTIL_ inline const char * const __export getMonthName(sel_typ mm)
{
    return ext_math.months_names[mm-1];
}

// returns the number of the day correspondent to the given INPUT date
__MSUTIL_ inline const sel_typ __export getDayNumber(sel_typ dd, sel_typ mm, uint64_t yyyy)
{
    return ((yyyy -= mm < 3) + yyyy/4 - yyyy/100 + yyyy/400 + ext_math.getdaynum_map[mm-1] + dd) % MAX_DAYSWEEK;
}
//

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
__MSNATIVE_ __MSUTIL_ inline void __export getPascalTriangle(uint64_t raws)
{
    uint64_t i, j;
    PRINTL();

    for(i = 0; i < raws; ++i)
    {
        for (j = 0; j < raws - i; ++j)
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

__MSNATIVE_ inline bool __system __export trigonometric_domain(register ityp x, register ityp y, register ityp z)
{
    for(uint64_t i=0; i<=x; ++i)
        if(x == y + (i*z)) return false;
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
        if(access(sysMem)[FUNCTION_FIBONACCI].current_max_index < access(curLayout)->max_memoizable_indeces[FUNCTION_FIBONACCI])
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
        if(access(sysMem)[FUNCTION_FATTORIALE].current_max_index < access(curLayout)->max_memoizable_indeces[FUNCTION_FATTORIALE])
        {
            if(access(sysMem)[FUNCTION_FATTORIALE].current_max_index)
                res = access(sysMem)[FUNCTION_FATTORIALE].memoizer[access(sysMem)[FUNCTION_FATTORIALE].current_max_index];

            for(i=(access(sysMem)[FUNCTION_FATTORIALE].current_max_index ? access(sysMem)[FUNCTION_FATTORIALE].current_max_index+1 : 2); i<=num; ++i)
            {
                access(sysMem)[FUNCTION_FATTORIALE].memoizer = realloc(access(sysMem)[FUNCTION_FATTORIALE].memoizer, sizeof(ityp)*(i+1));
                errMem(access(sysMem)[FUNCTION_FATTORIALE].memoizer, res);
            	access(sysMem)[FUNCTION_FATTORIALE].memoizer[i] = i<access(curLayout)->min_stirlingrequires_number ?  (res *= i) : (res=stirling(i));
            }

            access(sysMem)[FUNCTION_FATTORIALE].current_max_index = num;
        }
        else
        {
			if(num<access(curLayout)->min_stirlingrequires_number)
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
        if(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index < access(curLayout)->max_memoizable_indeces[FUNCTION_EVEN_SFATTORIALE])
        {
            if(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index)
                res = access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer[access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index];

            for(i=(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index ? access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index+2 : 2); i<=num; i+=2)
            {
                access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer = realloc(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer, sizeof(ityp)*(i+1));
                errMem(access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer, res);
                access(sysMem)[FUNCTION_EVEN_SFATTORIALE].memoizer[i] = (access(curLayout)->min_stirlingrequires_number%2 || num < (access(curLayout)->min_stirlingrequires_number<<1)) ? (res *= i) : (res = exp2(num)*stirling(num));
            }

            access(sysMem)[FUNCTION_EVEN_SFATTORIALE].current_max_index = num;
        }
        else
    	{
    		if(access(curLayout)->min_stirlingrequires_number%2 || (num < access(curLayout)->min_stirlingrequires_number<<1))
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
        if(access(sysMem)[FUNCTION_ODD_SFATTORIALE].current_max_index < access(curLayout)->max_memoizable_indeces[FUNCTION_ODD_SFATTORIALE])
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

/*
__MSUTIL_ inline ityp __export sfact(register ityp n)
{
    return (n <= 1.00) ? 1.00 : n*sfact(n-2.00);
}
*/

__MSUTIL_ inline ityp __export perm(register ityp n, register ityp r)
{
    return (fact(n)/(fact(n-r)));
}

__MSUTIL_ inline ityp __export comb(register ityp n, register ityp r)
{
    return (fact(n)/(fact(n-r)*fact(r)));
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
    ityp res = 0.00;
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
	return math_outlier2(dim, outlier_idx, OUTLIER_CONSTANT, vector);
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
    return (pow(productory(dim, false, vector), (1/dim)));
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

__MSUTIL_ inline ityp __export rootnX(register ityp n, register ityp X)
{
    return(pow(X, 1/n));
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

__MSNATIVE_ inline bool __system __export isEqualMatrix(ityp **matrix1, ityp **matrix2, const register dim_typ dim[static MAX_DIMENSIONS])
{
    dim_typ i, j;

    for(i=0; i<dim[RAWS]; ++i)
        for(j=0; j<dim[COLUMNS]; ++j)
            if(matrix1[i][j] != matrix2[i][j])
                return false;

    return true;

}

/// massive thanks to Bibek Subedi. Link at:
/// http://www.programming-technique.blogspot.it/2011/09/numerical-methods-condition-number-and.html
__MSNATIVE_ __MSUTIL_ ityp __system __export norms(ityp **matrix, dim_typ dim)
{
    dim_typ i, j;
    ityp sum;

    for(i=sum=0; i < dim; ++i)
        for(j = 0; j < dim; ++j)
            sum += mpow2(matrix[i][j], 2);

    return sqrt(sum);
}

/// thanks to: apatriarca, a User of Matematicamente forum, who is the author of this
/// pseudo-coded duet of functions which I renamed, unified, modified and adapted to this program.
/// This function calculates the Norm, One or Infinity (depending on mode bool param), of a Square Matrix.
/// http://www.matematicamente.it/forum/viewtopic.php?t=51498&p=371845
__MSNATIVE_ __MSUTIL_ ityp __system __export norm(ityp **matrix, dim_typ dim, bool mode)
{
    dim_typ i, j;
    ityp max; //  = 0.00; /* tutte le somme sono >= 0 */

    for (j = max = 0; j < dim; ++j)
    {
        ityp sum = 0.00;
        for (i = 0; i < dim; ++i)
            sum += fabs(matrix[mode ? j:i][mode ? i:j]);
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
// into _matrixProduct and _matrixSum functions,
// generated by the Inverse Operations Attribute Check
__MSNATIVE_ inline ityp __system __export math_sum(register ityp a, register ityp b)
{
    return (a+b);
}

__MSNATIVE_ inline ityp __system __export math_sub(register ityp a, register ityp b)
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
__MSNATIVE_ sel_typ _MS__private __system __export _simplexMethod(ityp ***tableau, ityp ***bfs, const register dim_typ dim[static MAX_DIMENSIONS], const bool constraint_types[static dim[RAWS]-1], bool mode)
{

    dim_typ i, j;
    const dim_typ dim_minus1[MAX_DIMENSIONS] =
    {
        dim[RAWS]-1,
        dim[COLUMNS]-1
    };

    const dim_typ vardims = dim[COLUMNS]+dim_minus1[RAWS];

    printf("\nreached x1\n");

    /// assert = false;

    for(i=0; i<dim[RAWS]; ++i)
    {
        (*tableau)[i] = realloc((*tableau)[i], sizeof(ityp)*vardims);
        errMem((*tableau)[i], SIMPLEXMETHOD_ALLOC_ERROR);
        const bool assert = i == dim_minus1[RAWS];
        (*tableau)[i][vardims-1] = assert ? 0.00 : (*tableau)[i][dim_minus1[COLUMNS]];
        #pragma omp parallel for
        for(j=dim_minus1[COLUMNS]; j<vardims-1; ++j)
            (*tableau)[i][j] = assert ? 0.00 : ((1-(2*constraint_types[i]))*(i+dim_minus1[COLUMNS] == j));
    }

    (*bfs)[RAWS] = realloc((*bfs)[RAWS], sizeof(ityp)*(vardims-1));
    errMem(bfs, SIMPLEXMETHOD_ALLOC_ERROR);


    // INITIALIZING BFS PHASE....
    #pragma omp parallel for
    for(i=0; i<dim_minus1[COLUMNS]; ++i)
        (*bfs)[RAWS][i] = 0.00;

    j=0;
    for(i=dim_minus1[COLUMNS]; i<vardims-1; ++i)
        (*bfs)[RAWS][i] = (*tableau)[j++][vardims-1];


    ityp func_vector[vardims];
    ityp min_ratios_vector[dim_minus1[RAWS]];
    dim_typ itemcheck;

    bool (* const _cmpfunc)(const register ityp, const register ityp) = MAX_PROBLEM ? max_cmpfunc : min_cmpfunc;

    dim_typ bestvaltax_idx, leaving_var_idx;
    bestvaltax_idx = leaving_var_idx = vardims;
    bool another_iteration = false;

    for(dim_typ control=0; ; )
    {
        if(control++ > access(curLayout)->max_simplex_iterations)
            return SIMPLEXMETHOD_FARBFS_ERROR;

        // vector reversing operation...
        for(i=another_iteration=0; i<vardims; ++i)
            if(_cmpfunc((func_vector[i] = (*tableau)[dim_minus1[RAWS]][i]), 0.00))
                another_iteration = true;

        if(!another_iteration)
            break;

        MINMAX(vardims, func_vector, mode, &bestvaltax_idx);

        for(i=itemcheck=0; i<dim_minus1[RAWS]; ++i)
        {
            if((*tableau)[i][bestvaltax_idx] <= 0)
            {
                if(++itemcheck == dim_minus1[RAWS])
                    return SIMPLEXMETHOD_INFBFS_ERROR;
                min_ratios_vector[i] = MAX_VAL;
            }
            else
                min_ratios_vector[i] = ((*tableau)[i][vardims-1])/((*tableau)[i][bestvaltax_idx]);
        }

        MINMAX(dim_minus1[RAWS], min_ratios_vector, MIN_MODE, &leaving_var_idx);

        const ityp pivot = (*tableau)[leaving_var_idx][bestvaltax_idx];

        // gauss-jordan's method application to the tableau, inclusive of function coefficients.
        for(i=0; i<vardims; ++i)
            (*tableau)[leaving_var_idx][i] /= pivot;

        ityp zero_binder = 0.00;

        for(i=0; i<dim[RAWS]; ++i)
        {
            if(i == leaving_var_idx || !(*tableau)[i][bestvaltax_idx])
                continue;

            zero_binder = (1-2*(((*tableau)[i][bestvaltax_idx]*(*tableau)[leaving_var_idx][bestvaltax_idx])>0))*((*tableau)[i][bestvaltax_idx])/((*tableau)[leaving_var_idx][bestvaltax_idx]);

			for(j=0; j<vardims; ++j)
                (*tableau)[i][j] += (*tableau)[leaving_var_idx][j]*zero_binder;
        }

        j=0;
        (*bfs)[RAWS][bestvaltax_idx] = (*tableau)[leaving_var_idx][vardims-1];
        for(i=dim_minus1[COLUMNS]; i<vardims-1; ++i)
            if(i == dim_minus1[COLUMNS]+leaving_var_idx)
                (*bfs)[RAWS][i] = 0.00;
            else if((*bfs)[RAWS][i])
                (*bfs)[RAWS][i] = (*tableau)[j++][vardims-1];

    }

    return SIMPLEXMETHOD_FOUNDBFS;
}

__MSNATIVE_ _MS__private void __system __export _matrixOEMTSum(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_sum, const register dim_typ i, const register dim_typ j)
{
	ityp (* const sum_func)(register ityp, register ityp) = INVERSE_OPS ? math_sub : math_sum;
	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
	{
        matrix_sum[OCTONIONS_REALPART][i][j] = sum_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_REALPART][i][j]);
        matrix_sum[OCTONIONS_E1PART][i][j] = sum_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E1PART][i][j]);
        matrix_sum[OCTONIONS_E2PART][i][j] = sum_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E2PART][i][j]);
        matrix_sum[OCTONIONS_E3PART][i][j] = sum_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E3PART][i][j]);
        matrix_sum[OCTONIONS_E4PART][i][j] = sum_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E4PART][i][j]);
        matrix_sum[OCTONIONS_E5PART][i][j] = sum_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E5PART][i][j]);
        matrix_sum[OCTONIONS_E6PART][i][j] = sum_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E6PART][i][j]);
        matrix_sum[OCTONIONS_E7PART][i][j] = sum_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E7PART][i][j]);
	}
	return;
}

__MSNATIVE_ _MS__private void __system __export _matrixOMTSum(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_sum, const register dim_typ i, const register dim_typ j)
{
	ityp (* const sum_func)(register ityp, register ityp) = INVERSE_OPS ? math_sub : math_sum;
	#pragma omp parallel
	{
        matrix_sum[OCTONIONS_REALPART][i][j] = sum_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_REALPART][i][j]);
        matrix_sum[OCTONIONS_E1PART][i][j] = sum_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E1PART][i][j]);
        matrix_sum[OCTONIONS_E2PART][i][j] = sum_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E2PART][i][j]);
        matrix_sum[OCTONIONS_E3PART][i][j] = sum_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E3PART][i][j]);
        matrix_sum[OCTONIONS_E4PART][i][j] = sum_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E4PART][i][j]);
        matrix_sum[OCTONIONS_E5PART][i][j] = sum_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E5PART][i][j]);
        matrix_sum[OCTONIONS_E6PART][i][j] = sum_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E6PART][i][j]);
        matrix_sum[OCTONIONS_E7PART][i][j] = sum_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E7PART][i][j]);
	}
	return;
}

__MSNATIVE_ _MS__private void __system __export _matrixSEMTSum(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_sum, const register dim_typ i, const register dim_typ j)
{
	ityp (* const sum_func)(register ityp, register ityp) = INVERSE_OPS ? math_sub : math_sum;
	#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
	{
		matrix_sum[SEDENIONS_REALPART][i][j] = sum_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_REALPART][i][j]);
        matrix_sum[SEDENIONS_E1PART][i][j] = sum_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E1PART][i][j]);
        matrix_sum[SEDENIONS_E2PART][i][j] = sum_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E2PART][i][j]);
        matrix_sum[SEDENIONS_E3PART][i][j] = sum_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E3PART][i][j]);
        matrix_sum[SEDENIONS_E4PART][i][j] = sum_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E4PART][i][j]);
        matrix_sum[SEDENIONS_E5PART][i][j] = sum_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E5PART][i][j]);
        matrix_sum[SEDENIONS_E6PART][i][j] = sum_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E6PART][i][j]);
        matrix_sum[SEDENIONS_E7PART][i][j] = sum_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E7PART][i][j]);
        matrix_sum[SEDENIONS_E8PART][i][j] = sum_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E8PART][i][j]);
        matrix_sum[SEDENIONS_E9PART][i][j] = sum_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E9PART][i][j]);
        matrix_sum[SEDENIONS_E10PART][i][j] = sum_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E10PART][i][j]);
        matrix_sum[SEDENIONS_E11PART][i][j] = sum_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E11PART][i][j]);
        matrix_sum[SEDENIONS_E12PART][i][j] = sum_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E12PART][i][j]);
        matrix_sum[SEDENIONS_E13PART][i][j] = sum_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E13PART][i][j]);
        matrix_sum[SEDENIONS_E14PART][i][j] = sum_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E14PART][i][j]);
        matrix_sum[SEDENIONS_E15PART][i][j] = sum_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E15PART][i][j]);
	}

	return;
}

__MSNATIVE_ _MS__private void __system __export _matrixSMTSum(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_sum, const register dim_typ i, const register dim_typ j)
{
	ityp (* const sum_func)(register ityp, register ityp) = INVERSE_OPS ? math_sub : math_sum;
	#pragma omp parallel
	{
		matrix_sum[SEDENIONS_REALPART][i][j] = sum_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_REALPART][i][j]);
        matrix_sum[SEDENIONS_E1PART][i][j] = sum_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E1PART][i][j]);
        matrix_sum[SEDENIONS_E2PART][i][j] = sum_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E2PART][i][j]);
        matrix_sum[SEDENIONS_E3PART][i][j] = sum_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E3PART][i][j]);
        matrix_sum[SEDENIONS_E4PART][i][j] = sum_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E4PART][i][j]);
        matrix_sum[SEDENIONS_E5PART][i][j] = sum_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E5PART][i][j]);
        matrix_sum[SEDENIONS_E6PART][i][j] = sum_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E6PART][i][j]);
        matrix_sum[SEDENIONS_E7PART][i][j] = sum_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E7PART][i][j]);
        matrix_sum[SEDENIONS_E8PART][i][j] = sum_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E8PART][i][j]);
        matrix_sum[SEDENIONS_E9PART][i][j] = sum_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E9PART][i][j]);
        matrix_sum[SEDENIONS_E10PART][i][j] = sum_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E10PART][i][j]);
        matrix_sum[SEDENIONS_E11PART][i][j] = sum_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E11PART][i][j]);
        matrix_sum[SEDENIONS_E12PART][i][j] = sum_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E12PART][i][j]);
        matrix_sum[SEDENIONS_E13PART][i][j] = sum_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E13PART][i][j]);
        matrix_sum[SEDENIONS_E14PART][i][j] = sum_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E14PART][i][j]);
        matrix_sum[SEDENIONS_E15PART][i][j] = sum_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E15PART][i][j]);
	}

	return;
}

__MSNATIVE_ _MS__private void __system __export _matrixOEMTProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ i, const register dim_typ j, const register dim_typ k)
{
	ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;
	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
	{
		matrix_product[OCTONIONS_REALPART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_REALPART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E7PART][k][j]);

        matrix_product[OCTONIONS_E1PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E3PART][k][j]);

        matrix_product[OCTONIONS_E2PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E6PART][k][j]);

        matrix_product[OCTONIONS_E3PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E1PART][k][j]);

        matrix_product[OCTONIONS_E4PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E5PART][k][j]);

        matrix_product[OCTONIONS_E5PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E4PART][k][j]);

        matrix_product[OCTONIONS_E6PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E2PART][k][j]);

        matrix_product[OCTONIONS_E7PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_REALPART][k][j]);
	}
	return;
}

__MSNATIVE_ _MS__private void __system __export _matrixOMTProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ i, const register dim_typ j, const register dim_typ k)
{
	ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;
	#pragma omp parallel
	{
		matrix_product[OCTONIONS_REALPART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_REALPART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) -
                                                            mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E7PART][k][j]);

        matrix_product[OCTONIONS_E1PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E3PART][k][j]);

        matrix_product[OCTONIONS_E2PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E6PART][k][j]);

        matrix_product[OCTONIONS_E3PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E1PART][k][j]);

        matrix_product[OCTONIONS_E4PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E5PART][k][j]);

        matrix_product[OCTONIONS_E5PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E4PART][k][j]);

        matrix_product[OCTONIONS_E6PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_E2PART][k][j]);

        matrix_product[OCTONIONS_E7PART][i][j] += mul_func(matrix1[OCTONIONS_REALPART][i][k], matrix2[OCTONIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E1PART][i][k], matrix2[OCTONIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E2PART][i][k], matrix2[OCTONIONS_E6PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E3PART][i][k], matrix2[OCTONIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E4PART][i][k], matrix2[OCTONIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E5PART][i][k], matrix2[OCTONIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[OCTONIONS_E6PART][i][k], matrix2[OCTONIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[OCTONIONS_E7PART][i][k], matrix2[OCTONIONS_REALPART][k][j]);
	}
	return;
}

__MSNATIVE_ _MS__private void __system __export _matrixSEMTProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ i, const register dim_typ j, const register dim_typ k)
{
	ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;
	#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
	{
		matrix_product[SEDENIONS_REALPART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E15PART][k][j]);

        matrix_product[SEDENIONS_E1PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E14PART][k][j]);

        matrix_product[SEDENIONS_E2PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E13PART][k][j]);

        matrix_product[SEDENIONS_E3PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E12PART][k][j]);

        matrix_product[SEDENIONS_E4PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E11PART][k][j]);

        matrix_product[SEDENIONS_E5PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E10PART][k][j]);

        matrix_product[SEDENIONS_E6PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E15PART][k][j]);

        matrix_product[SEDENIONS_E7PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E8PART][k][j]);

        matrix_product[SEDENIONS_E8PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E7PART][k][j]);

        matrix_product[SEDENIONS_E9PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E9PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E6PART][k][j]);

        matrix_product[SEDENIONS_E10PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E10PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E5PART][k][j]);

        matrix_product[SEDENIONS_E11PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E4PART][k][j]);

        matrix_product[SEDENIONS_E12PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E3PART][k][j]);

        matrix_product[SEDENIONS_E13PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E2PART][k][j]);

        matrix_product[SEDENIONS_E14PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E1PART][k][j]);

        matrix_product[SEDENIONS_E15PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_REALPART][k][j]);
	}
}

__MSNATIVE_ _MS__private void __system __export _matrixSMTProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ i, const register dim_typ j, const register dim_typ k)
{
	ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;
	#pragma omp parallel
	{
		matrix_product[SEDENIONS_REALPART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                            mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E15PART][k][j]);

        matrix_product[SEDENIONS_E1PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E14PART][k][j]);

        matrix_product[SEDENIONS_E2PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E13PART][k][j]);

        matrix_product[SEDENIONS_E3PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E12PART][k][j]);

        matrix_product[SEDENIONS_E4PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E11PART][k][j]);

        matrix_product[SEDENIONS_E5PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E10PART][k][j]);

        matrix_product[SEDENIONS_E6PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E15PART][k][j]);

        matrix_product[SEDENIONS_E7PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E8PART][k][j]);

        matrix_product[SEDENIONS_E8PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E7PART][k][j]);

        matrix_product[SEDENIONS_E9PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E9PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                  mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                  mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E6PART][k][j]);

        matrix_product[SEDENIONS_E10PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E10PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E5PART][k][j]);

        matrix_product[SEDENIONS_E11PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E4PART][k][j]);

        matrix_product[SEDENIONS_E12PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E3PART][k][j]);

        matrix_product[SEDENIONS_E13PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E2PART][k][j]);

        matrix_product[SEDENIONS_E14PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E15PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E7PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_REALPART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_E1PART][k][j]);

        matrix_product[SEDENIONS_E15PART][i][j] += mul_func(matrix1[SEDENIONS_REALPART][i][k], matrix2[SEDENIONS_E15PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E1PART][i][k], matrix2[SEDENIONS_E14PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E2PART][i][k], matrix2[SEDENIONS_E13PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E3PART][i][k], matrix2[SEDENIONS_E12PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E4PART][i][k], matrix2[SEDENIONS_E11PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E5PART][i][k], matrix2[SEDENIONS_E10PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E6PART][i][k], matrix2[SEDENIONS_E9PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E7PART][i][k], matrix2[SEDENIONS_E8PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E8PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E9PART][i][k], matrix2[SEDENIONS_E6PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E10PART][i][k], matrix2[SEDENIONS_E5PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E11PART][i][k], matrix2[SEDENIONS_E4PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E12PART][i][k], matrix2[SEDENIONS_E3PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E13PART][i][k], matrix2[SEDENIONS_E2PART][k][j]) -
                                                   mul_func(matrix1[SEDENIONS_E14PART][i][k], matrix2[SEDENIONS_E1PART][k][j]) +
                                                   mul_func(matrix1[SEDENIONS_E15PART][i][k], matrix2[SEDENIONS_REALPART][k][j]);
	}
}

__MSNATIVE_ _MS__private void __system __export _matrixKOEMTProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS], const register dim_typ i, const register dim_typ j, const register dim_typ k, const register dim_typ l)
{
	ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;
	const register dim_typ dim_cache[MAX_DIMENSIONS] =
	{
		k+(i*dim[SECOND_MATRIX][RAWS]),
		l+(j*dim[SECOND_MATRIX][COLUMNS])
	};
	#pragma omp parallel num_threads(MAX_OCTONIONS_UNITS)
	{
        matrix_product[OCTONIONS_REALPART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_REALPART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E7PART][k][l]);

        matrix_product[OCTONIONS_E1PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E3PART][k][l]);

        matrix_product[OCTONIONS_E2PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E1PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E6PART][k][l]);

        matrix_product[OCTONIONS_E3PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E1PART][k][l]);

        matrix_product[OCTONIONS_E4PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E5PART][k][l]);

        matrix_product[OCTONIONS_E5PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E3PART][k][l]);

        matrix_product[OCTONIONS_E6PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E2PART][k][l]);

        matrix_product[OCTONIONS_E7PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_REALPART][k][l]);
	}
	return;
}

__MSNATIVE_ _MS__private void __system __export _matrixKOMTProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS], const register dim_typ i, const register dim_typ j, const register dim_typ k, const register dim_typ l)
{
	ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;
	const register dim_typ dim_cache[MAX_DIMENSIONS] =
	{
		k+(i*dim[SECOND_MATRIX][RAWS]),
		l+(j*dim[SECOND_MATRIX][COLUMNS])
	};
	#pragma omp parallel
	{
        matrix_product[OCTONIONS_REALPART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_REALPART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) -
                                                                                  mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E7PART][k][l]);

        matrix_product[OCTONIONS_E1PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E3PART][k][l]);

        matrix_product[OCTONIONS_E2PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E1PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E6PART][k][l]);

        matrix_product[OCTONIONS_E3PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E1PART][k][l]);

        matrix_product[OCTONIONS_E4PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E5PART][k][l]);

        matrix_product[OCTONIONS_E5PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E3PART][k][l]);

        matrix_product[OCTONIONS_E6PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_E2PART][k][l]);

        matrix_product[OCTONIONS_E7PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[OCTONIONS_REALPART][i][j], matrix2[OCTONIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E1PART][i][j], matrix2[OCTONIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E2PART][i][j], matrix2[OCTONIONS_E6PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E3PART][i][j], matrix2[OCTONIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E4PART][i][j], matrix2[OCTONIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E5PART][i][j], matrix2[OCTONIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[OCTONIONS_E6PART][i][j], matrix2[OCTONIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[OCTONIONS_E7PART][i][j], matrix2[OCTONIONS_REALPART][k][l]);
	}
	return;
}

__MSNATIVE_ _MS__private void __system __export _matrixKSEMTProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS], const register dim_typ i, const register dim_typ j, const register dim_typ k, const register dim_typ l)
{
	ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;
	const register dim_typ dim_cache[MAX_DIMENSIONS] =
	{
		k+(i*dim[SECOND_MATRIX][RAWS]),
		l+(j*dim[SECOND_MATRIX][COLUMNS])
	};
	#pragma omp parallel num_threads(MAX_SEDENIONS_UNITS)
	{
        matrix_product[SEDENIONS_REALPART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E15PART][k][l]);

        matrix_product[SEDENIONS_E1PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E14PART][k][l]);

        matrix_product[SEDENIONS_E2PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E13PART][k][l]);

        matrix_product[SEDENIONS_E3PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E12PART][k][l]);

        matrix_product[SEDENIONS_E4PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E11PART][k][l]);

        matrix_product[SEDENIONS_E5PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E10PART][k][l]);

        matrix_product[SEDENIONS_E6PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E9PART][k][l]);

        matrix_product[SEDENIONS_E7PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E8PART][k][l]);

        matrix_product[SEDENIONS_E8PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E7PART][k][l]);

        matrix_product[SEDENIONS_E9PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E9PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E6PART][k][l]);

        matrix_product[SEDENIONS_E10PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E10PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E5PART][k][l]);

        matrix_product[SEDENIONS_E11PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E4PART][k][l]);

        matrix_product[SEDENIONS_E12PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E3PART][k][l]);

        matrix_product[SEDENIONS_E13PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E2PART][k][l]);

        matrix_product[SEDENIONS_E14PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E1PART][k][l]);

        matrix_product[SEDENIONS_E15PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E14PART][k][l]);
	}

	return;
}

__MSNATIVE_ _MS__private void __system __export _matrixKSMTProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS], const register dim_typ i, const register dim_typ j, const register dim_typ k, const register dim_typ l)
{
	ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;
	const register dim_typ dim_cache[MAX_DIMENSIONS] =
	{
		k+(i*dim[SECOND_MATRIX][RAWS]),
		l+(j*dim[SECOND_MATRIX][COLUMNS])
	};
	#pragma omp parallel
	{
        matrix_product[SEDENIONS_REALPART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                  mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E15PART][k][l]);

        matrix_product[SEDENIONS_E1PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E14PART][k][l]);

        matrix_product[SEDENIONS_E2PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E13PART][k][l]);

        matrix_product[SEDENIONS_E3PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E12PART][k][l]);

        matrix_product[SEDENIONS_E4PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E11PART][k][l]);

        matrix_product[SEDENIONS_E5PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E10PART][k][l]);

        matrix_product[SEDENIONS_E6PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E9PART][k][l]);

        matrix_product[SEDENIONS_E7PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E8PART][k][l]);

        matrix_product[SEDENIONS_E8PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E7PART][k][l]);

        matrix_product[SEDENIONS_E9PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E9PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E6PART][k][l]);

        matrix_product[SEDENIONS_E10PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E10PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E5PART][k][l]);

        matrix_product[SEDENIONS_E11PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E4PART][k][l]);

        matrix_product[SEDENIONS_E12PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E3PART][k][l]);

        matrix_product[SEDENIONS_E13PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E2PART][k][l]);

        matrix_product[SEDENIONS_E14PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E15PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E7PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_REALPART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E1PART][k][l]);

        matrix_product[SEDENIONS_E15PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[SEDENIONS_REALPART][i][j], matrix2[SEDENIONS_E15PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E1PART][i][j], matrix2[SEDENIONS_E14PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E2PART][i][j], matrix2[SEDENIONS_E13PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E3PART][i][j], matrix2[SEDENIONS_E12PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E4PART][i][j], matrix2[SEDENIONS_E11PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E5PART][i][j], matrix2[SEDENIONS_E10PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E6PART][i][j], matrix2[SEDENIONS_E9PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E7PART][i][j], matrix2[SEDENIONS_E8PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E8PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E9PART][i][j], matrix2[SEDENIONS_E6PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E10PART][i][j], matrix2[SEDENIONS_E5PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E11PART][i][j], matrix2[SEDENIONS_E4PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E12PART][i][j], matrix2[SEDENIONS_E3PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E13PART][i][j], matrix2[SEDENIONS_E2PART][k][l]) -
                                                                                 mul_func(matrix1[SEDENIONS_E14PART][i][j], matrix2[SEDENIONS_E1PART][k][l]) +
                                                                                 mul_func(matrix1[SEDENIONS_E15PART][i][j], matrix2[SEDENIONS_E14PART][k][l]);
	}

	return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixSum(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
    ityp (* const sum_func)(register ityp, register ityp) = INVERSE_OPS ? math_sub : math_sum;

	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
            (*matrix_sum)[i][j] = sum_func((*matrix1)[i][j], (*matrix2)[i][j]);

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixCSum(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
    ityp (* const sum_func)(register ityp, register ityp) = INVERSE_OPS ? math_sub : math_sum;

	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
        	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
       		{
	            matrix_sum[REAL_PART][i][j] = sum_func(matrix1[REAL_PART][i][j], matrix2[REAL_PART][i][j]);
	            matrix_sum[IMAG_PART][i][j] = sum_func(matrix1[IMAG_PART][i][j], matrix2[IMAG_PART][i][j]);
	    	}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixQSum(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
    ityp (* const sum_func)(register ityp, register ityp) = INVERSE_OPS ? math_sub : math_sum;

	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
        	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
        	{
	            matrix_sum[QUATERNIONS_REALPART][i][j] = sum_func(matrix1[QUATERNIONS_REALPART][i][j], matrix2[QUATERNIONS_REALPART][i][j]);
	            matrix_sum[QUATERNIONS_IPART][i][j] = sum_func(matrix1[QUATERNIONS_IPART][i][j], matrix2[QUATERNIONS_IPART][i][j]);
	            matrix_sum[QUATERNIONS_JPART][i][j] = sum_func(matrix1[QUATERNIONS_JPART][i][j], matrix2[QUATERNIONS_JPART][i][j]);
	            matrix_sum[QUATERNIONS_KPART][i][j] = sum_func(matrix1[QUATERNIONS_KPART][i][j], matrix2[QUATERNIONS_KPART][i][j]);
	    	}

    return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixOSum(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
    void (* const omp_func)(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ) = isSett(BOOLS_EXTENSIVEMULTITHREADING) ? _matrixOEMTSum : _matrixOMTSum;

	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
        	omp_func(matrix1, matrix2, matrix_sum, i, j);

    return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixSSum(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_sum, const register dim_typ dim[static MAX_DIMENSIONS])
{
	dim_typ i, j;
	void (* const omp_func)(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ) = isSett(BOOLS_EXTENSIVEMULTITHREADING) ? _matrixSEMTSum : _matrixSMTSum;


	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
        	omp_func(matrix1, matrix2, matrix_sum, i, j);

    return;
}


__MSNATIVE_ inline void _MS__private __system __export _matrixProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_ABSTRACT_DIMENSIONS]) // dim_typ righe, dim_typ colonne, dim_typ colonne2)
{
    dim_typ i, j, k;
    ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;

	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS2]; ++j)
        	#pragma omp parallel for
            for(k=(*matrix_product)[i][j]=0; k<dim[COLUMNS]; ++k)
                (*matrix_product)[i][j] += mul_func((*matrix1)[i][k], (*matrix2)[k][j]); /// INVERSE_OPS ? (matrix1[i][k] / matrix2[k][j]):(matrix1[i][k] * matrix2[k][j]);

                        //
                        // wasn't able to check OVERFLOW here, for technical and formatting reasons.

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixCProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_ABSTRACT_DIMENSIONS])
{
    dim_typ i, j, k;
    ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;

	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS2]; ++j)
        	#pragma omp parallel for
            for(k=matrix_product[REAL_PART][i][j]=matrix_product[IMAG_PART][i][j]=0; k<dim[COLUMNS]; ++k)
            	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
            	{
	                matrix_product[REAL_PART][i][j] += mul_func(matrix1[REAL_PART][i][k], matrix2[REAL_PART][k][j]) - mul_func(matrix1[IMAG_PART][i][k], matrix2[IMAG_PART][k][j]);
	                matrix_product[IMAG_PART][i][j] += mul_func(matrix1[REAL_PART][i][k], matrix2[IMAG_PART][k][j]) + mul_func(matrix1[IMAG_PART][i][k], matrix2[REAL_PART][k][j]);
	       		}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixQProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_ABSTRACT_DIMENSIONS])
{
    dim_typ i, j, k;
    ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;

	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS2]; ++j)
        	#pragma omp parallel for
            for(k=matrix_product[QUATERNIONS_REALPART][i][j]=
                matrix_product[QUATERNIONS_IPART][i][j]   =
                matrix_product[QUATERNIONS_JPART][i][j]   =
                matrix_product[QUATERNIONS_KPART][i][j]=0; k<dim[COLUMNS]; ++k)
            	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
            	{
	                matrix_product[QUATERNIONS_REALPART][i][j] += mul_func(matrix1[QUATERNIONS_REALPART][i][k], matrix2[QUATERNIONS_REALPART][k][j]) -
	                                                              mul_func(matrix1[QUATERNIONS_IPART][i][k], matrix2[QUATERNIONS_IPART][k][j]) -
	                                                              mul_func(matrix1[QUATERNIONS_JPART][i][k], matrix2[QUATERNIONS_JPART][k][j]) -
	                                                              mul_func(matrix1[QUATERNIONS_KPART][i][k], matrix2[QUATERNIONS_KPART][k][j]);
	                matrix_product[QUATERNIONS_IPART][i][j] += mul_func(matrix1[QUATERNIONS_REALPART][i][k], matrix2[QUATERNIONS_IPART][k][j]) +
	                                                           mul_func(matrix1[QUATERNIONS_IPART][i][k], matrix2[QUATERNIONS_REALPART][k][j]) +
	                                                           mul_func(matrix1[QUATERNIONS_JPART][i][k], matrix2[QUATERNIONS_KPART][k][j]) -
	                                                           mul_func(matrix1[QUATERNIONS_KPART][i][k], matrix2[QUATERNIONS_JPART][k][j]);
	                matrix_product[QUATERNIONS_JPART][i][j] += mul_func(matrix1[QUATERNIONS_REALPART][i][k], matrix2[QUATERNIONS_JPART][k][j]) +
	                                                           mul_func(matrix1[QUATERNIONS_JPART][i][k], matrix2[QUATERNIONS_REALPART][k][j]) +
	                                                           mul_func(matrix1[QUATERNIONS_KPART][i][k], matrix2[QUATERNIONS_IPART][k][j]) -
	                                                           mul_func(matrix1[QUATERNIONS_IPART][i][k], matrix2[QUATERNIONS_KPART][k][j]);
	                matrix_product[QUATERNIONS_KPART][i][j] += mul_func(matrix1[QUATERNIONS_REALPART][i][k], matrix2[QUATERNIONS_KPART][k][j]) +
	                                                           mul_func(matrix1[QUATERNIONS_KPART][i][k], matrix2[QUATERNIONS_REALPART][k][j]) +
	                                                           mul_func(matrix1[QUATERNIONS_IPART][i][k], matrix2[QUATERNIONS_JPART][k][j]) -
	                                                           mul_func(matrix1[QUATERNIONS_JPART][i][k], matrix2[QUATERNIONS_IPART][k][j]);
	            }

    return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixOProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_ABSTRACT_DIMENSIONS])
{
    dim_typ i, j, k;
	void (* const omp_func)(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ, const register dim_typ) = isSett(BOOLS_EXTENSIVEMULTITHREADING) ? _matrixOEMTProduct : _matrixOMTProduct;

	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS2]; ++j)
        	#pragma omp parallel for
            for(k=matrix_product[OCTONIONS_REALPART][i][j]=
                matrix_product[OCTONIONS_E1PART][i][j]  =
                matrix_product[OCTONIONS_E2PART][i][j]  =
                matrix_product[OCTONIONS_E3PART][i][j]  =
                matrix_product[OCTONIONS_E4PART][i][j]  =
                matrix_product[OCTONIONS_E5PART][i][j]  =
                matrix_product[OCTONIONS_E6PART][i][j]  =
                matrix_product[OCTONIONS_E7PART][i][j]=0; k<dim[COLUMNS]; ++k)
            	omp_func(matrix1, matrix2, matrix_product, i, j, k);
    return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixSProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_ABSTRACT_DIMENSIONS])
{
    dim_typ i, j, k;
	void (* const omp_func)(ityp ***, ityp ***, ityp ***, const register dim_typ, const register dim_typ, const register dim_typ) = isSett(BOOLS_EXTENSIVEMULTITHREADING) ? _matrixSEMTProduct : _matrixSMTProduct;

	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS2]; ++j)
        	#pragma omp parallel for
            for(k=matrix_product[SEDENIONS_REALPART][i][j]=
                matrix_product[SEDENIONS_E1PART][i][j]  =
                matrix_product[SEDENIONS_E2PART][i][j]  =
                matrix_product[SEDENIONS_E3PART][i][j]  =
                matrix_product[SEDENIONS_E4PART][i][j]  =
                matrix_product[SEDENIONS_E5PART][i][j]  =
                matrix_product[SEDENIONS_E6PART][i][j]  =
                matrix_product[SEDENIONS_E7PART][i][j]  =
                matrix_product[SEDENIONS_E8PART][i][j]  =
                matrix_product[SEDENIONS_E9PART][i][j]  =
                matrix_product[SEDENIONS_E10PART][i][j]  =
                matrix_product[SEDENIONS_E11PART][i][j]  =
                matrix_product[SEDENIONS_E12PART][i][j]  =
                matrix_product[SEDENIONS_E13PART][i][j]  =
                matrix_product[SEDENIONS_E14PART][i][j]  =
                matrix_product[SEDENIONS_E15PART][i][j]=0; k<dim[COLUMNS]; ++k)
            	omp_func(matrix1, matrix2, matrix_product, i, j, k);
    return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixKProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS])
{
    dim_typ i, j, k, l;
    ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;

	#pragma omp parallel for
    for(i=0; i<dim[FIRST_MATRIX][RAWS]; ++i)
	#pragma omp parallel for
        for(j=0; j<dim[FIRST_MATRIX][COLUMNS]; ++j)
        	#pragma omp parallel for
            for(k=0; k<dim[SECOND_MATRIX][RAWS]; ++k)
            	#pragma omp parallel for
                for(l=0; l<dim[SECOND_MATRIX][COLUMNS]; ++l)
                    (*matrix_product)[k+(i*dim[SECOND_MATRIX][RAWS])][l+(j*dim[SECOND_MATRIX][COLUMNS])] = mul_func((*matrix1)[i][j], (*matrix2)[k][l]);

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixKCProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS])
{
    dim_typ i, j, k, l;
    dim_typ dim_cache[MAX_DIMENSIONS];
    ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;

	#pragma omp parallel for
    for(i=0; i<dim[FIRST_MATRIX][RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[FIRST_MATRIX][COLUMNS]; ++j)
        	#pragma omp parallel for
            for(k=0; k<dim[SECOND_MATRIX][RAWS]; ++k)
            	#pragma omp parallel for
                for(l=0; l<dim[SECOND_MATRIX][COLUMNS]; ++l)
                	#pragma omp parallel num_threads(MAX_COMPLEX_UNITS)
                	{
	                    dim_cache[RAWS] = k+(i*dim[SECOND_MATRIX][RAWS]);
	                    dim_cache[COLUMNS] = l+(j*dim[SECOND_MATRIX][COLUMNS]);
	                    matrix_product[REAL_PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[REAL_PART][i][j], matrix2[REAL_PART][k][l]) - mul_func(matrix1[IMAG_PART][i][j], matrix2[IMAG_PART][k][l]);
	                    matrix_product[IMAG_PART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[REAL_PART][i][j], matrix2[IMAG_PART][k][l]) - mul_func(matrix1[IMAG_PART][i][j], matrix2[REAL_PART][k][l]);
	            	}

    return;
}

__MSNATIVE_ void _MS__private __system __export _matrixKQProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS])
{
    dim_typ i, j, k, l;
    dim_typ dim_cache[MAX_DIMENSIONS];
    ityp (* const mul_func)(register ityp, register ityp) = INVERSE_OPS ? math_div : math_mul;

	#pragma omp parallel for
    for(i=0; i<dim[FIRST_MATRIX][RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[FIRST_MATRIX][COLUMNS]; ++j)
        	#pragma omp parallel for
            for(k=0; k<dim[SECOND_MATRIX][RAWS]; ++k)
            	#pragma omp parallel for
                for(l=0; l<dim[SECOND_MATRIX][COLUMNS]; ++l)
                	#pragma omp parallel num_threads(MAX_QUATERNIONS_UNITS)
                	{
	                    dim_cache[RAWS] = k+(i*dim[SECOND_MATRIX][RAWS]);
	                    dim_cache[COLUMNS] = l+(j*dim[SECOND_MATRIX][COLUMNS]);
	                    matrix_product[QUATERNIONS_REALPART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[QUATERNIONS_REALPART][i][j], matrix2[QUATERNIONS_REALPART][k][l]) -
	                                                                                                                              mul_func(matrix1[QUATERNIONS_IPART][i][j], matrix2[QUATERNIONS_IPART][k][l]) -
	                                                                                                                              mul_func(matrix1[QUATERNIONS_JPART][i][j], matrix2[QUATERNIONS_JPART][k][l]) -
	                                                                                                                              mul_func(matrix1[QUATERNIONS_KPART][i][j], matrix2[QUATERNIONS_KPART][k][l]);

	                    matrix_product[QUATERNIONS_IPART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[QUATERNIONS_REALPART][i][j], matrix2[QUATERNIONS_IPART][k][l]) +
	                                                                                                                           mul_func(matrix1[QUATERNIONS_IPART][i][j], matrix2[QUATERNIONS_REALPART][k][l])+
	                                                                                                                           mul_func(matrix1[QUATERNIONS_JPART][i][j], matrix2[QUATERNIONS_KPART][k][l]) -
	                                                                                                                           mul_func(matrix1[QUATERNIONS_KPART][i][j], matrix2[QUATERNIONS_JPART][k][l]);

	                    matrix_product[QUATERNIONS_JPART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[QUATERNIONS_REALPART][i][j], matrix2[QUATERNIONS_JPART][k][l]) +
	                                                                                                                           mul_func(matrix1[QUATERNIONS_JPART][i][j], matrix2[QUATERNIONS_REALPART][k][l])+
	                                                                                                                           mul_func(matrix1[QUATERNIONS_KPART][i][j], matrix2[QUATERNIONS_IPART][k][l]) -
	                                                                                                                           mul_func(matrix1[QUATERNIONS_IPART][i][j], matrix2[QUATERNIONS_KPART][k][l]);

	                    matrix_product[QUATERNIONS_KPART][dim_cache[RAWS]][dim_cache[COLUMNS]] = mul_func(matrix1[QUATERNIONS_REALPART][i][j], matrix2[QUATERNIONS_KPART][k][l]) +
	                                                                                                                           mul_func(matrix1[QUATERNIONS_KPART][i][j], matrix2[QUATERNIONS_REALPART][k][l])+
	                                                                                                                           mul_func(matrix1[QUATERNIONS_IPART][i][j], matrix2[QUATERNIONS_JPART][k][l]) -
	                                                                                                                           mul_func(matrix1[QUATERNIONS_JPART][i][j], matrix2[QUATERNIONS_IPART][k][l]);
	    			}

    return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixKOProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS])
{
    dim_typ i, j, k, l;
	void (* const omp_func)(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS], const register dim_typ, const register dim_typ, const register dim_typ, const register dim_typ) = isSett(BOOLS_EXTENSIVEMULTITHREADING) ? _matrixKOEMTProduct : _matrixKOMTProduct;

    #pragma omp parallel for
    for(i=0; i<dim[FIRST_MATRIX][RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[FIRST_MATRIX][COLUMNS]; ++j)
        	#pragma omp parallel for
            for(k=0; k<dim[SECOND_MATRIX][RAWS]; ++k)
            	#pragma omp parallel for
                for(l=0; l<dim[SECOND_MATRIX][COLUMNS]; ++l)
                	omp_func(matrix1, matrix2, matrix_product, dim, i, j, k, l);
    return;
}

__MSNATIVE_ inline void _MS__private __system __export _matrixKSProduct(ityp ***matrix1, ityp ***matrix2, ityp ***matrix_product, const register dim_typ dim[static MAX_DIMENSIONS][MAX_DIMENSIONS])
{

    dim_typ i, j, k, l;
	void (* const omp_func)(ityp ***, ityp ***, ityp ***, const register dim_typ [static MAX_DIMENSIONS][MAX_DIMENSIONS], const register dim_typ, const register dim_typ, const register dim_typ, const register dim_typ) = isSett(BOOLS_EXTENSIVEMULTITHREADING) ? _matrixKSEMTProduct : _matrixKSMTProduct;

	#pragma omp parallel for
    for(i=0; i<dim[FIRST_MATRIX][RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[FIRST_MATRIX][COLUMNS]; ++j)
        	#pragma omp parallel for
            for(k=0; k<dim[SECOND_MATRIX][RAWS]; ++k)
            	#pragma omp parallel for
                for(l=0; l<dim[SECOND_MATRIX][COLUMNS]; ++l)
                	omp_func(matrix1, matrix2, matrix_product, dim, i, j, k, l);
    return;
}

// ENDAlgebra

// END
