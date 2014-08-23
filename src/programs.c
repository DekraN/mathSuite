// programs.c 23/08/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.
*/

#include "dutils.h" // DA RENDERE VISIBILE       SIA AL COMPILATORE CHE AL LINKER
// #include "ExprEval/exprincl.h" // In order To use Redefined MATH_ constants


static void printOpersIdentifiers(void);

struct operations operazioni[MAX_OPERATIONS];

static const struct operations default_operazioni[MAX_OPERATIONS] =
{
    {
        "Exit from Program",
        NULL_CHAR,
        {
            EXIT_CMD
        },
        DOMAIN_INT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Sum",
        NULL_CHAR,
        {
            OPERATOR_ADDIZIONE,
            IDENTIFIER_ADDIZIONE
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Subtraction",
        NULL_CHAR,
        {
            OPERATOR_SOTTRAZIONE,
            IDENTIFIER_SOTTRAZIONE
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Multiplication",
        NULL_CHAR,
        {
            OPERATOR_MOLTIPLICAZIONE,
            IDENTIFIER_MOLTIPLICAZIONE
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Division",
        NULL_CHAR,
        {
            OPERATOR_DIVISIONE,
            IDENTIFIER_DIVISIONE
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Rest",
        NULL_CHAR,
        {
            OPERATOR_RESTO,
            IDENTIFIER_RESTO,
            IDENTIFIER_ALIAS_RESTO
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Binary Sum",
        NULL_CHAR,
        {
            OPERATOR_ADDIZIONEBINARIA,
            IDENTIFIER_ADDIZIONEBINARIA
        },
        DOMAIN_INT,
        DOMAIN_INT,
        binary_function
    },
    {
        "Binary Subtraction",
        NULL_CHAR,
        {
            OPERATOR_SOTTRAZIONEBINARIA,
            IDENTIFIER_SOTTRAZIONEBINARIA
        },
        DOMAIN_INT,
        DOMAIN_INT,
        binary_function
    },
    {
        "Complement",
        NULL_CHAR,
        {
            OPERATOR_COMPLEMENT,
            IDENTIFIER_COMPLEMENT
        },
        DOMAIN_INT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Exponentiation",
        NULL_CHAR,
        {
            OPERATOR_ELEVAMENTOAPOTENZA,
            IDENTIFIER_ELEVAMENTOAPOTENZA
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Exponential (Cardinal)",
        "(a == 0 -> exp, 1 -> expc)",
        {
            OPERATOR_EXPANDEXPC,
            IDENTIFIER_EXPANDEXPC,
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function,
    },
    {
        "10 Base Exponential (Cardinal)",
        "(a == 0 -> exp10, 1 -> exp10c)",
        {
            OPERATOR_EXP10ANDEXP10C,
            OPERATOR_ALIAS_EXP10ANDEXP10C,
            IDENTIFIER_EXP10ANDEXP10C
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "2 Base Exponential (Cardinal)",
        "(a == 0 -> exp2, 1 -> exp2c)",
        {
            OPERATOR_EXP2ANDEXP2C,
            OPERATOR_ALIAS_EXP2ANDEXP2C,
            IDENTIFIER_EXP2ANDEXP2C
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "N Root",
        NULL_CHAR,
        {
            IDENTIFIER_RADICENESIMA
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Logarithm",
        "(a == 0 -> ln, 1 -> log)",
        {
            IDENTIFIER_LOGARITMO,
            IDENTIFIER_ALIAS_LOGARITMO,
            IDENTIFIER_ALIAS2_LOGARITMO,
            IDENTIFIER_ALIAS3_LOGARITMO
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Base 2 Logarithm",
        NULL_CHAR,
        {
            IDENTIFIER_LOGARITMO2,
            IDENTIFIER_ALIAS_LOGARITMO2,
            IDENTIFIER_ALIAS2_LOGARITMO2,
            IDENTIFIER_ALIAS3_LOGARITMO2
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "N Base Logarithm",
        "(b as N)",
        {
            IDENTIFIER_LOGARITMOBN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Cardinal Logarithm",
        "(a == 0 -> lnc, 1 -> logc)",
        {
            IDENTIFIER_LOGARITMOC
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Base 2 Cardinal Logarithm",
        NULL_CHAR,
        {
            IDENTIFIER_LOGARITMO2C
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Plus 1 Logarithm",
        NULL_CHAR,
        {
            IDENTIFIER_LOGARITMO1P
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Plus 1 Cardinal Logarithm",
        NULL_CHAR,
        {
            IDENTIFIER_LOGARITMO1PC
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Bit Counter",
        NULL_CHAR,
        {
            IDENTIFIER_BITCOUNTER,
            IDENTIFIER_ALIAS_BITCOUNTER,
            IDENTIFIER_ALIAS2_BITCOUNTER,
            IDENTIFIER_ALIAS3_BITCOUNTER,
            IDENTIFIER_ALIAS4_BITCOUNTER
        },
        DOMAIN_LONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Unsigned Bit Counter",
        NULL_CHAR,
        {
            IDENTIFIER_UBITCOUNTER,
            IDENTIFIER_ALIAS_UBITCOUNTER,
            IDENTIFIER_ALIAS2_UBITCOUNTER,
            IDENTIFIER_ALIAS3_UBITCOUNTER,
            IDENTIFIER_ALIAS4_UBITCOUNTER
        },
        DOMAIN_ULONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Show Prog Version",
        NULL_CHAR,
        {
            IDENTIFIER_VERSION,
            IDENTIFIER_ALIAS_VERSION,
            IDENTIFIER_ALIAS2_VERSION,
            IDENTIFIER_ALIAS3_VERSION
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Change Precision",
        NULL_CHAR,
        {
            IDENTIFIER_PREC
        },
        DOMAIN_USHRT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Change Stabilizer Factor",
        NULL_CHAR,
        {
            IDENTIFIER_STABFACT,
            IDENTIFIER_ALIAS_STABFACT,
            IDENTIFIER_ALIAS2_STABFACT
        },
        DOMAIN_USHRT,
        DOMAIN_NULL,
        unary_function
    },
    {
    	"Change Min Stirling-Requires Number",
    	NULL_CHAR,
    	{
    		IDENTIFIER_MINSRNUMBER,
    		IDENTIFIER_ALIAS_MINSRNUMBER,
    		IDENTIFIER_ALIAS2_MINSRNUMBER,
    		IDENTIFIER_ALIAS3_MINSRNUMBER
    	},
    	DOMAIN_USHRT,
    	DOMAIN_NULL,
    	unary_function
    },
    {
        "Change Algebra",
        NULL_CHAR,
        {
            IDENTIFIER_ALGEBRA,
            IDENTIFIER_ALIAS_ALGEBRA
        },
        DOMAIN_USHRT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Change Random Seed",
        NULL_CHAR,
        {
            IDENTIFIER_RSEED,
            IDENTIFIER_ALIAS_RSEED,
            IDENTIFIER_ALIAS2_RSEED
        },
        DOMAIN_USHRT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Change Fibonacci MMI",
        NULL_CHAR,
        {
            IDENTIFIER_MMIFIBO,
            IDENTIFIER_ALIAS_MMIFIBO
        },
        DOMAIN_USHRT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Change Factorial MMI",
        NULL_CHAR,
        {
            IDENTIFIER_MMIFACT,
            IDENTIFIER_ALIAS_MMIFACT
        },
        DOMAIN_USHRT,
        DOMAIN_NULL,
        unary_function
    },
    {
    	"Change Even SemiFactorial MMI",
    	NULL_CHAR,
    	{
    		IDENTIFIER_MMIEVENSFACT,
    		IDENTIFIER_ALIAS_MMIEVENSFACT
    	},
    	DOMAIN_USHRT,
    	DOMAIN_NULL,
    	unary_function
    },
    {
    	"Change Odd SemiFactorial MMI",
    	NULL_CHAR,
    	{
    		IDENTIFIER_MMIODDSFACT,
    		IDENTIFIER_ALIAS_MMIODDSFACT
    	},
    	DOMAIN_USHRT,
    	DOMAIN_NULL,
    	unary_function
    },
    {
        "Transform Angles",
        "(a == 0 -> Degrees to Radiant, 1 -> Radiant to Degrees)",
        {
            IDENTIFIER_TRASFORMAANGOLI
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Sin (Hyperbolic)",
        "(a == 0 -> sin, 1 -> sinh)",
        {
            IDENTIFIER_SINANDSINH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Cos (Hyperbolic)",
        "(a == 0 -> cos, 1 -> cosh)",
        {
            IDENTIFIER_COSANDCOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Tan (Hyperbolic)",
        "(a == 0 -> tan, 1 -> tanh)",
        {
            IDENTIFIER_TANANDTANH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Csc (Hyperbolic)",
        "(a == 0 -> csc, 1 -> csch)",
        {
            IDENTIFIER_CSCANDCSCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Sec (Hyperbolic)",
        "(a == 0 -> sec, 1 -> sech)",
        {
            IDENTIFIER_SECANDSECH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Cot (Hyperbolic)",
        "(a == 0 -> cot, 1 -> coth)",
        {
            IDENTIFIER_COTANDCOTH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Hsin (Hyperbolic)",
        "(a == 0 -> hsin, 1 -> hsinh)",
        {
            IDENTIFIER_HSINANDHSINH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qsin (Hyperbolic)",
        "(a == 0 -> qsin, 1 -> qsinh)",
        {
            IDENTIFIER_QSINANDQSINH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Hcos (Hyperbolic)",
        "(a == 0 -> hcos, 1 -> hcosh)",
        {
            IDENTIFIER_HCOSANDHCOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qcos (Hyperbolic)",
        "(a == 0 -> qcos, 1 -> qcosh)",
        {
            IDENTIFIER_QCOSANDQCOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Hsec (Hyperbolic)",
        "(a == 0 -> hsec, 1 -> hsech)",
        {
            IDENTIFIER_HSECANDHSECH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qsec (Hyperbolic)",
        "(a == 0 -> qsec, 1 -> qsech)",
        {
            IDENTIFIER_QSECANDQSECH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Hcsc (Hyperbolic)",
        "(a == 0 -> hcsc, 1 -> hcsch)",
        {
            IDENTIFIER_HCSCANDHCSCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qcsc (Hyperbolic)",
        "(a == 0 -> qcsc, 1 -> qcsch)",
        {
            IDENTIFIER_QCSCANDQCSCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Htan (Hyperbolic)",
        "(a == 0 -> htan, 1 -> htanh)",
        {
            IDENTIFIER_HTANANDHTANH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qtan (Hyperbolic)",
        "(a == 0 -> qtan, 1 -> qtanh)",
        {
            IDENTIFIER_QTANANDQTANH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Hcot (Hyperbolic)",
        "(a == 0 -> hcot, 1 -> hcoth)",
        {
            IDENTIFIER_HCOTANDHCOTH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qcot (Hyperbolic)",
        "(a == 0 -> qcot, 1 -> qcoth)",
        {
            IDENTIFIER_QCOTANDQCOTH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Vsin (Hyperbolic)",
        "(a == 0 -> vsin, 1 -> vsinh)",
        {
            IDENTIFIER_VSINANDVSINH,
            IDENTIFIER_ALIAS_VSINANDVSINH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "CVsin (Hyperbolic)",
        "(a == 0 -> cvsin, 1 -> cvsinh)",
        {
            IDENTIFIER_CVSINANDCVSINH,
            IDENTIFIER_ALIAS_CVSINANDCVSINH,
            IDENTIFIER_ALIAS2_CVSINANDCVSINH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Vcos (Hyperbolic)",
        "(a == 0 -> vcos, 1 -> vcosh)",
        {
            IDENTIFIER_VCOSANDVCOSH,
            IDENTIFIER_ALIAS_VCOSANDVCOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "CVcos (Hyperbolic)",
        "(a == 0 -> cvcos, 1 -> cvcosh)",
        {
            IDENTIFIER_CVCOSANDCVCOSH,
            IDENTIFIER_ALIAS_CVCOSANDCVCOSH,
            IDENTIFIER_ALIAS2_CVCOSANDCVCOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "HVsin (Hyperbolic)",
        "(a == 0 -> hvsin, 1 -> hvsinh)",
        {
            IDENTIFIER_HVSINANDHVSINH,
            IDENTIFIER_ALIAS_HVSINANDHVSINH,
            IDENTIFIER_ALIAS2_HVSINANDHVSINH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "HCVsin (Hyperbolic)",
        "(a == 0 -> hcvsin, 1 -> hcvsinh)",
        {
            IDENTIFIER_HCVSINANDHCVSINH,
            IDENTIFIER_ALIAS_HCVSINANDHCVSINH,
            IDENTIFIER_ALIAS2_HCVSINANDHCVSINH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "QVsin (Hyperbolic)",
        "(a == 0 -> qvsin, 1 -> qvsinh)",
        {
            IDENTIFIER_QVSINANDQVSINH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "QCVsin (Hyperbolic)",
        "(a == 0 -> qcvsin, 1 -> qcvsinh)",
        {
            IDENTIFIER_QCOSANDQCOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "HVcos (Hyperbolic)",
        "(a == 0 -> hvcos, 1 -> hvcosh)",
        {
            IDENTIFIER_HVCOSANDHVCOSH,
            IDENTIFIER_ALIAS_HVCOSANDHVCOSH,
            IDENTIFIER_ALIAS2_HVCOSANDHVCOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "HCVcos (Hyperbolic)",
        "(a == 0 -> hcvcos, 1 -> hcvcosh)",
        {
            IDENTIFIER_HCVCOSANDHCVCOSH,
            IDENTIFIER_ALIAS_HCVCOSANDHCVCOSH,
            IDENTIFIER_ALIAS2_HCVCOSANDHCVCOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "QVcos (Hyperbolic)",
        "(a == 0 -> qvcos, 1 -> qvcosh)",
        {
            IDENTIFIER_QVCOSANDQVCOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "QCVcos (Hyperbolic)",
        "(a == 0 -> qcvcos, 1 -> qcvcosh)",
        {
            IDENTIFIER_QCVCOSANDQCVCOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Esec (Hyperbolic)",
        "(a == 0 -> esec, 1 -> esech)",
        {
            IDENTIFIER_ESECANDESECH,
            IDENTIFIER_ALIAS_ESECANDESECH,
            IDENTIFIER_ALIAS2_ESECANDESECH,
            IDENTIFIER_ALIAS3_ESECANDESECH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Ecsc (Hyperbolic)",
        "(a == 0 -> esec, 1 -> esech)",
        {
            IDENTIFIER_ECSCANDECSCH,
            IDENTIFIER_ALIAS_ECSCANDECSCH,
            IDENTIFIER_ALIAS2_ECSCANDECSCH,
            IDENTIFIER_ALIAS3_ECSCANDECSCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "HEsec (Hyperbolic)",
        "(a == 0 -> hesec, 1 -> hesech)",
        {
            IDENTIFIER_HESECANDHESECH,
            IDENTIFIER_ALIAS_HESECANDHESECH,
            IDENTIFIER_ALIAS2_HESECANDHESECH,
            IDENTIFIER_ALIAS3_HESECANDHESECH,
            IDENTIFIER_ALIAS4_HESECANDHESECH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "HEcsc (Hyperbolic)",
        "(a == 0 -> hecsc, 1 -> hesech)",
        {
            IDENTIFIER_HECSCANDHECSCH,
            IDENTIFIER_ALIAS_HECSCANDHECSCH,
            IDENTIFIER_ALIAS2_HECSCANDHECSCH,
            IDENTIFIER_ALIAS3_HECSCANDHECSCH,
            IDENTIFIER_ALIAS4_HECSCANDHECSCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "QEsec (Hyperbolic)",
        "(a == 0 -> qesec, 1 -> qesech)",
        {
            IDENTIFIER_QESECANDQESECH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "QEcsc (Hyperbolic)",
        "(a == 0 -> qecsc, 1 -> qecsch)",
        {
            IDENTIFIER_QECSCANDQECSCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "sinc (Hyperbolic)",
        "(a == 0 -> sinc, 1 -> sinch)",
        {
            IDENTIFIER_SINCANDSINCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Hsinc (Hyperbolic)",
        "(a == 0 -> hsinc, 1 -> hsinch)",
        {
            IDENTIFIER_HSINCANDHSINCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qsinc (Hyperbolic)",
        "(a == 0 -> qsinc, 1 -> qsinch)",
        {
            IDENTIFIER_QSINCANDQSINCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "cosc (Hyperbolic)",
        "(a == 0 -> cosc, 1 -> cosch)",
        {
            IDENTIFIER_COSCANDCOSCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Hcosc (Hyperbolic)",
        "(a == 0 -> hcosc, 1 -> hcosch)",
        {
            IDENTIFIER_HCOSCANDHCOSCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qcosc (Hyperbolic)",
        "(a == 0 -> qcosc, 1 -> qcosch)",
        {
            IDENTIFIER_QCOSCANDQCOSCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "secc (Hyperbolic)",
        "(a == 0 -> secc, 1 -> secch)",
        {
            IDENTIFIER_SECCANDSECCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Hsecc (Hyperbolic)",
        "(a == 0 -> hsecc, 1 -> hsecch)",
        {
            IDENTIFIER_HSECCANDHSECCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qsecc (Hyperbolic)",
        "(a == 0 -> qsecc, 1 -> qsecch)",
        {
            IDENTIFIER_QSECCANDQSECCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "cscc (Hyperbolic)",
        "(a == 0 -> cscc, 1 -> cscch)",
        {
            IDENTIFIER_CSCCANDCSCCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Hcscc (Hyperbolic)",
        "(a == 0 -> hcscc, 1 -> hcscch)",
        {
            IDENTIFIER_HCSCCANDHCSCCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qcscc (Hyperbolic)",
        "(a == 0 -> qcscc, 1 -> qcscch)",
        {
            IDENTIFIER_QCSCCANDQCSCCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "tanc (Hyperbolic)",
        "(a == 0 -> tanc, 1 -> tanhc)",
        {
            IDENTIFIER_TANCANDTANCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Htanc (Hyperbolic)",
        "(a == 0 -> htanc, 1 -> htanch)",
        {
            IDENTIFIER_HTANCANDHTANCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qtanc (Hyperbolic)",
        "(a == 0 -> qtanc, 1 -> qtanch)",
        {
            IDENTIFIER_QTANCANDQTANCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "cotc (Hyperbolic)",
        "(a == 0 -> cotc, 1 -> cotch)",
        {
            IDENTIFIER_COTCANDCOTCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Hcotc (Hyperbolic)",
        "(a == 0 -> hcotc, 1 -> hcotch)",
        {
            IDENTIFIER_HCOTCANDHCOTCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Qcotc (Hyperbolic)",
        "(a == 0 -> qcotc, 1 -> qcotch)",
        {
            IDENTIFIER_QCOTCANDQCOTCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Asin (Hyperbolic)",
        "(a == 0 -> asin, 1 -> asinh)",
        {
            IDENTIFIER_ASINANDASINH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Acos (Hyperbolic)",
        "(a == 0 -> acos, 1 -> acosh)",
        {
            IDENTIFIER_ACOSANDACOSH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Atan (Hyperbolic)",
        "(a == 0 -> atan, 1 -> atanh)",
        {
            IDENTIFIER_ATANANDATANH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "y/x atan2 in: ]-M_PI,M_PI]",
        "(a as y, b as x)",
        {
            IDENTIFIER_ATAN2
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Acsc (Hyperbolic)",
        "(a == 0 -> acsc, 1 -> acsch)",
        {
            IDENTIFIER_ACSCANDACSCH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Asec (Hyperbolic)",
        "(a == 0 -> asec, 1 -> asech)",
        {
            IDENTIFIER_ASECANDASECH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Acot (Hyperbolic)",
        "(a == 0 -> acot, 1 -> acoth)",
        {
            IDENTIFIER_ACOTANDACOTH
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Greatest Common Divisor",
        "("IDENTIFIER_ALIAS_MCD"(a,b))",
        {
            IDENTIFIER_ALIAS_MCD,
            IDENTIFIER_MCD
        },
        DOMAIN_ULLONG,
        DOMAIN_ULLONG,
        binary_function
    },
    {
        "Least Common Multiple",
        "("IDENTIFIER_ALIAS_MCM"(a,b))",
        {
            IDENTIFIER_ALIAS_MCM,
            IDENTIFIER_MCM
        },
        DOMAIN_ULLONG,
        DOMAIN_ULLONG,
        binary_function
    },
    {
        "Approximation",
        "(a == 0 -> floor, 1 ceil)",
        {
            IDENTIFIER_APPROSSIMAZIONE
        },
        DOMAIN_BOOL,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Geometric Succession Sum",
        "(a = (base > 1), n = exp)",
        {
            IDENTIFIER_SOMMASUCCESSIONEGEOMETRICA
        },
        DOMAIN_DEFAULT,
        DOMAIN_ULLONG,
        binary_function
    },
    {
        "Armonic Succession Sum",
        "(a = k -> sum = (1/k))",
        {
            IDENTIFIER_SOMMASUCCESSIONEARMONICA
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Generalized Armonic Succession Sum",
        "(a = (base > 0), n = exp)",
        {
            IDENTIFIER_SOMMASUCCESSIONEARMONICAGEN
        },
        DOMAIN_ULLONG,
        DOMAIN_ULLONG,
        binary_function
    },
    {
        "Fibonacci Succession Sum",
        NULL_CHAR,
        {
            IDENTIFIER_SOMMASUCCESSIONEFIBONACCI
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Factorial Succession Sum",
        NULL_CHAR,
        {
            IDENTIFIER_SOMMASUCCESSIONEFATTORIALE
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "SemiFactorial Succession Sum",
        NULL_CHAR,
        {
            IDENTIFIER_SOMMASUCCESSIONESEMIFATTORIALE
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Summation",
        "(a != 0 -> Increase Sum, 0 -> Reset Sum)",
        {
            IDENTIFIER_SOMMATORIA
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Productory",
        "(a != 0 -> Increase Sum, 0 -> Reset Sum)",
        {
            IDENTIFIER_PRODUTTORIA
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Media",
        "(b == 1 -> Enter another Element, 0 -> Stop Entering)",
        {
            IDENTIFIER_MEDIA,
            IDENTIFIER_ALIAS_MEDIA,
            IDENTIFIER_ALIAS2_MEDIA,
            IDENTIFIER_ALIAS3_MEDIA
        },
        DOMAIN_DEFAULT,
        DOMAIN_BOOL,
        binary_function
    },
    {
    	"Variance",
    	"(b == 1 -> Enter another Element, 0 -> Stop Entering)",
    	{
    		IDENTIFIER_VARIANCE,
    		IDENTIFIER_ALIAS_VARIANCE
    	},
    	DOMAIN_DEFAULT,
    	DOMAIN_BOOL,
    	binary_function
    },
    {
    	"Standard Deviation",
    	"(b == 1 -> Enter another Element, 0 -> Stop Entering)",
    	{
    		IDENTIFIER_STDDEV,
    		IDENTIFIER_ALIAS_STDDEV,
    		IDENTIFIER_ALIAS2_STDDEV
    	},
    	DOMAIN_DEFAULT,
    	DOMAIN_BOOL,
    	binary_function
    },
    {
    	"Outlier Test",
    	"(b = Outlier Index, Enter different value to Stop Entering)",
    	{
    		IDENTIFIER_OUTLIER,
    		IDENTIFIER_ALIAS_OUTLIER,
    		IDENTIFIER_ALIAS2_OUTLIER
    	},
    	DOMAIN_DEFAULT,
    	DOMAIN_LLONG,
    	binary_function
    },
    {
        "Geometric Media",
        "(a != 0 -> Increasing Media, 0 -> Stop Media)",
        {
            IDENTIFIER_MEDIAGEOMETRICA
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Armonic Media",
        "(a != 0 -> Increasing Media, 0 -> Stop Media)",
        {
            IDENTIFIER_MEDIAARMONICA
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Power Media",
        "(b != 0 -> Increasing Media, Both Two's 0 -> Stop Media, b = Media Root != 0)",
        {
            IDENTIFIER_MEDIAPOTENZA
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Mid Value",
        "(a != 0 -> Increasing Media, 0 -> Stop Media)",
        {
            IDENTIFIER_VALORECENTRALE
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
    	"First Quartile",
    	"(a != 0 -> Enter another Element, 0 -> Stop Entering)",
    	{
    		IDENTIFIER_FIRSTQUARTILE,
    		IDENTIFIER_ALIAS_FIRSTQUARTILE,
    		IDENTIFIER_ALIAS2_FIRSTQUARTILE
    	},
    	DOMAIN_DEFAULT,
    	DOMAIN_BOOL,
    	unary_function
    },
    {
        "Median",
        "(a != 0 -> Entering another Element, 0 -> Stop Entering)",
        {
            IDENTIFIER_MEDIANA,
            IDENTIFIER_ALIAS_MEDIANA,
            IDENTIFIER_ALIAS2_MEDIANA
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
    	"Third Quartile",
    	"(a != 0 -> Enter another Element, 0 -> Stop Entering)",
    	{
    		IDENTIFIER_THIRDQUARTILE,
    		IDENTIFIER_ALIAS_THIRDQUARTILE,
    		IDENTIFIER_ALIAS2_THIRDQUARTILE
    	},
    	DOMAIN_DEFAULT,
    	DOMAIN_NULL,
    	unary_function
    },
    {
        "First N Number Sum",
        NULL_CHAR,
        {
            IDENTIFIER_SOMMAPRIMINNUMERI
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Square Root",
        NULL_CHAR,
        {
            IDENTIFIER_RADICEQUADRATA
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Cubic Root",
        NULL_CHAR,
        {
            IDENTIFIER_RADICECUBICA
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Factorial",
        "(b == 0 -> Shows Only that number, 1 -> Shows all the succession till it reaches Result)",
        {
            OPERATOR_FATTORIALE,
            IDENTIFIER_FATTORIALE,
            IDENTIFIER_ALIAS_FATTORIALE,
            IDENTIFIER_ALIAS2_FATTORIALE,
            IDENTIFIER_ALIAS3_FATTORIALE
        },
        DOMAIN_DEFAULT,
        DOMAIN_BOOL,
        binary_function
    },
    {
        "SemiFactorial",
        "(b == 0 -> Shows Only that number, 1 -> Shows all the succession till it reaches Result)",
        {
            OPERATOR_SEMIFATTORIALE,
            IDENTIFIER_SEMIFATTORIALE,
            IDENTIFIER_ALIAS_SEMIFATTORIALE,
            IDENTIFIER_ALIAS2_SEMIFATTORIALE,
            IDENTIFIER_ALIAS3_SEMIFATTORIALE
        },
        DOMAIN_DEFAULT,
        DOMAIN_BOOL,
        binary_function
    },
    {
    	"Stirling's Approximation",
    	"(b == 0 -> Shows Only that number, 1 -> Shows all the succession till it reaches Result)",
        {
            IDENTIFIER_STIRLING,
            IDENTIFIER_ALIAS_STIRLING,
            IDENTIFIER_ALIAS2_STIRLING,
            IDENTIFIER_ALIAS3_STIRLING
            
        },
        DOMAIN_DEFAULT,
        DOMAIN_BOOL,
        binary_function
	},
    {
        "Fibonacci Succession",
        "(b == 0 -> Shows Only that number, 1 -> Shows all the succession till it reaches Result)",
        {
            OPERATOR_FIBONACCI,
            IDENTIFIER_FIBONACCI,
            IDENTIFIER_ALIAS_FIBONACCI,
            IDENTIFIER_ALIAS2_FIBONACCI
        },
        DOMAIN_DEFAULT,
        DOMAIN_BOOL,
        binary_function
    },
    {
        "Permutations",
        NULL_CHAR,
        {
            IDENTIFIER_PERMUTATION
        },
        DOMAIN_ULLONG,
        DOMAIN_ULLONG,
        binary_function
    },
    {
        "Combinations",
        NULL_CHAR,
        {
            IDENTIFIER_COMBINATION
        },
        DOMAIN_ULLONG,
        DOMAIN_ULLONG,
        binary_function
    },
    {
        "Pascal's Triangle",
        NULL_CHAR,
        {
            IDENTIFIER_PASCALTRIANGLE
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Roman Numbers Conversion",
        NULL_CHAR,
        {
            IDENTIFIER_NUMERIROMANI
        },
        DOMAIN_INT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "First N Prime Numbers",
        NULL_CHAR,
        {
            IDENTIFIER_PRIMINNUMERIPRIMI
        },
        DOMAIN_ULLONG,
        DOMAIN_ULLONG,
        binary_function // unary_function
    },
    {
        "N-th Prime Number",
        NULL_CHAR,
        {
            IDENTIFIER_NESIMONUMEROPRIMO
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Primorial",
        NULL_CHAR,
        {
            OPERATOR_PRIMORIALE,
            IDENTIFIER_PRIMORIALE
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "First N Prime Number Sum",
        NULL_CHAR,
        {
            IDENTIFIER_SOMMAPRIMINNUMERIPRIMI
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Fibonaccial",
        NULL_CHAR,
        {
            IDENTIFIER_FIBONACCIALE
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Base Conversion",
        "(Switch a base: from 10 to b)",
        {
            IDENTIFIER_CAMBIAMENTODIBASE
        },
        DOMAIN_ULLONG,
        DOMAIN_UCHAR,
        binary_function
    },
    {
        "Random Matrix Generation",
        "(Generate an [aXb] Random Matrix)",
        {
            IDENTIFIER_GENERATOREMATRICIRANDOM
        },
        DOMAIN_USHRT,
        DOMAIN_USHRT,
        binary_function
    },
    {
        "Calculator Informations",
        NULL_CHAR,
        {
            OPERATOR_NOTANOPERATION,
            IDENTIFIER_NOTANOPERATION
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    }
};

__MSSHELL_WRAPPER_ void calcolatoreDiBase(const register sel_typ argc, char ** argv)
{

    if(isSett(BOOLS_BASECALCPARSING))
        // requires(NULL, yes to show result, yes to show varlist, yes to show difftime, yes to save results, not to reset lists
        requires(argc ? argv[0] : NULL, argc ? NULL : "Enter an Expression", "Result is", PARSER_SHOWRESULT | PARSER_SHOWVARLIST | PARSER_SHOWDIFFTIME | PARSER_SAVERESULT);
    else
    {

        #define viewInsertedValue                                                                               \
            printf2(COLOR_SYSTEM, "\nValue: ");                                                                               \
            printf2(COLOR_CREDITS, OUTPUT_CONVERSION_FORMAT, a);                                                                \
            printf2(COLOR_SYSTEM, " inserted in [%llu] position of elaborating Media Vector...\n\n", accumulate);    \
            return ;


        sel_typ tmp;
        char identifier[MAX_IDENTIFIER_LENGTH];

        printf2(COLOR_CREDITS, "\nEnter Operations INPUT as expected:\n");
        printf2(COLOR_SYSTEM, "a identifier b\n\n");

        ityp n1, n2;

        // Inizializzazione VARIABILI di CALCOLO
        // e di riconoscimento OPERAZIONE

        n1 = n2 = 0.00;
        strcpy(identifier, OPERATOR_NOTANOPERATION);

        CLEARBUFFER();

        if(argc)
        {
            if((tmp = sscanf(argv[0], "%lf %s %lf", &n1, identifier, &n2)) != 3 && tmp != 2)
            {
                printErr(5, "Invalid Format or non-existent Identifier");
                return ;
            }
        }
        else
            while((tmp = scanf("%lf %s %lf", &n1, identifier, &n2)) != 3 && tmp != 2)
            {
                CLEARBUFFER();
                printErr(5, "Invalid Format or non-existent Identifier");
            }

        register ityp a, b;

        a = n1;
        b = n2;

        // if(isDomainForbidden(a, INPUT) || isDomainForbidden(b, INPUT)) return;


        register ityp c;
        sel_typ output_ID[MAX_DIMENSIONS] =
        {
            BASECALC_INFORMAZIONI,
            0
        };

        c = 0;
        CLEARBUFFER();

        doesExistOperIdentifier(identifier, output_ID);

        const sel_typ oprID = output_ID[OPERATION_ID];

        const TypeRange minmax[MAX_DIMENSIONS] =
        {
            {
                ext_math.types_range[operazioni[oprID].domA].min,
                ext_math.types_range[operazioni[oprID].domA].max
            },
            {
                ext_math.types_range[operazioni[oprID].domB].min,
                ext_math.types_range[operazioni[oprID].domB].max
            }
        };

        if(a < minmax[LEFT_OPR].min || a > minmax[LEFT_OPR].max ||
        (operazioni[oprID].bin_or_unary && (b < minmax[RIGHT_OPR].min || b > minmax[RIGHT_OPR].max)))
        {
            printf2(COLOR_ERROR, "\nINPUT OVERFLOW ERROR.\n");
            printErr(33, "One of two values is incorrect. You have to enter:\na between %.*G and %.*G, and b between %.*G and %.*G",
                            DEFAULT_PRECISION, minmax[LEFT_OPR].min, DEFAULT_PRECISION, minmax[LEFT_OPR].max, DEFAULT_PRECISION, minmax[RIGHT_OPR].min, DEFAULT_PRECISION, minmax[RIGHT_OPR].max);
            return ;
        }

        const bool dg_ent = isSett(BOOLS_DEGREESENTERING);

        if(dg_ent && oprID > BASECALC_TRASFORMAANGOLI && oprID < BASECALC_ASINANDASINH && a == 0) b = rad(b);

        switch(oprID)
        {

            #include "bcilfs.h"

            default:
                printErr(1, "Non-existent Operation");
                return ;
        }

        #undef viewInsertedValue


        if(dg_ent && oprID > BASECALC_QECSCANDQECSCH && oprID < BASECALC_MCM && a == 0) c = deg(c);

        CLEARBUFFER();

        if(isDomainForbidden(c,OUTPUT)) return ;

        printf2(COLOR_USER, "\n\n%s Operation RESULT is: ", operazioni[oprID].name);
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, c);
        printf2(COLOR_USER, ".\n");
    }

    return ;
}

__MSSHELL_WRAPPER_ void __apnt calcolatoreAvanzato(const register sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_ADVCALC_PROGS, adv_calc,
                        main_menu[MAIN_ADVANCEDCALCULATOR].name,
                        #if MAX_ADVCALC_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
                            BY_NUMBERS
                        #else
                            BY_CHARS
                        #endif // MAX_SETTINGS
                        );
    return ;
}

__MSSHELL_WRAPPER_ void __apnt mhssManager(const register sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_MHSSMANAGER_PROGS, mhss_manager,
                        main_menu[MAIN_MHSSMANAGER].name,
                        #if MAX_MHSSMANAGER_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
                            BY_NUMBERS
                        #else
                            BY_CHARS
                        #endif // MAX_SETTINGS
                        );
    return ;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ void _MS__private __system __export operationsGroupMenu(dim_typ dim, sprog programs[static dim], const char access_point[static INFO_STRING], bool list)
{
    for( ;; )
    {

        printf("\nSelect desired Program:\n");
        PRINTL();
        // PRINTN();

        dim_typ i;
        dim_typ tmp;

        // tmp = PROGRAM_BUSY;

        CLEARBUFFER();

        if(list == BY_NUMBERS)
        {
            for(i=0; i<dim; ++i)
                printf("- %hu: %s;\n", i, programs[i].name);

            printf("- %hu: Clear SCREEN;\n", i);
            printf("- %hu: PROGRAM Informations;\n", i+1);
            printf("- %hu: Exit from PROGRAM.\n\n", i+2);
            PRINTL();


            ityp tmp2;

            while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp2 = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
                    scanf(INPUT_CONVERSION_FORMAT, &tmp2) != 1) || tmp2 != (tmp = (dim_typ)tmp2) || tmp < 0 || tmp > dim+2)
                printErr(1, "Invalid PROGRAM Mode");

        }
        else
        {
            // printf(programs[MAIN_ALGEBRAOPERATIONS].name);
            for(i=0; i<dim; ++i)
                printf("- %c: %s;\n", i+'A', programs[i].name);
            printf("- %c: Clear SCREEN;\n", i+'A');
            printf("- %c: PROGRAM Informations;\n", i+1+'A');
            printf("- %c: Exit from PROGRAM.\n\n", i+2+'A');
            PRINTL();

            sel_typ tmp2;

            do
                tmp2 = toupper(getch());
            while(tmp2 < 'A' || tmp2 > dim+2+'A');

            // tmp = tmp2-65;
            // tmp = (toupper(tmp2)-65);

            tmp = tmp2-'A';
        }

        __pmode__ = tmp;


        CLEARBUFFER();
        PRINT2N();


        if(tmp == dim)
        {
            pulisciSchermo;
            printf2(COLOR_USER, "\nSCREEN has been correctly cleaned.\n\n");
            PRINTL();
            operationsGroupMenu(dim, programs, access_point, list);
            return;
        }

        if(tmp == dim+1)
        {
            if(!strcmp(access_point, NULL_CHAR))
                progInfo(WITH_DESCRIPTION);
            else
            {
                char low_name[FILENAME_MAX] = NULL_CHAR;
                char tmpname[MAX_PATH_LENGTH] = DESCRIPTIONS_FOLDER;
                strfnm(access_point, low_name);
                strcat(tmpname, low_name);
                PRINT2N();
                printf2(COLOR_CREDITS, "\t\t       %s\n", access_point);

                if(!readFile(tmpname))
                    printErr(2, "Non-existent DESCRIPTION File in Directory "DESCRIPTIONS_FOLDER":\n%s", low_name);

            }

            continue;
        }

        if(tmp == dim+2) break;



        sprog prog_chosen = programs[tmp];


        // PRINTING PROGRAM NAME
        char str[INFO_STRING];

        strcpy(str, prog_chosen.name);
        toupper_s(str);
        printf2(COLOR_CREDITS, str);
        PRINTN();

        if(isSett(BOOLS_SHOWDESCRIPTION))
        {
            if(prog_chosen.program_function == calcolatoreDiBase)
                if(isSett(BOOLS_BASECALCPARSING))
                {
                    char low_name[FILENAME_MAX];
                    strfnm(prog_chosen.name, low_name);

                    if(file_exists(low_name))
                    {
                        #if WINOS
                            sprintf(str, "notepad %s", low_name);
                            system(str);
                        #else
                            readFile(low_name);
                        #endif
                    }
                    else
                    {
                        printErr(2, "Unable to open File:\n%s", low_name);

                        printf("\nDo you want to Reach Expression Help Page\n"); // : %s\n", EXPREVAL_TMPL);
                        printf("for more informations about EXPREVAL basilar inline functions?\n");
                        printf("[Y (Yes) / N (Not)]\n");

                        char selection;

                        do selection = toupper(getch());
                        while(selection != 'Y' && selection != 'N');

                        if(selection == 'Y')
                        {
                            sprintf(str, "START %s", EXPREVAL_TMPL);
                            system(str);
                        }

                        PRINTN();
                        CLEARBUFFER();

                    }

                }
                else
                    printOpersIdentifiers();
            else
            {
                char low_name[FILENAME_MAX];
                char tmpname[MAX_PATH_LENGTH] = DESCRIPTIONS_FOLDER;
                strfnm(prog_chosen.name, low_name);
                strcat(tmpname, low_name);
                printf2(COLOR_CREDITS, "\t\t       %s\n", prog_chosen.name);

                if(!readFile(tmpname))
                    printErr(2, "Non-existent DESCRIPTION File in Directory "DESCRIPTIONS_FOLDER":\n%s", low_name);

                printf2(COLOR_SYSTEM, ".\nCMDNAME: %s ; USAGE: %s .\n", prog_chosen.cmdname, prog_chosen.usage);
                PRINTL();
            }
        }

        volatile bool rep_check;

        do
        {
            // printf("\n_____________________________________________________\n\n");

            time_t t1;
            const bool asrt = isSett(BOOLS_SHOWEXECTIME);

            CLEARBUFFER();

            if(asrt)
                t1 = time(NULL);


            prog_chosen.program_function(0, NULL); // RICHIAMA LA FUNZIONE O METODO DEL PROGRAMMA SELEZIONATO
            // avendo a disposizione l'indirizzo della funzione corrispondente al subprogram scelto.
            if(asrt)
                printf("\nAverage Execution Time: %.*f.\n\n", DEFAULT_PRECISION, difftime(time(NULL), t1));

            if(isSett(BOOLS_RESETLISTS)) refreshExprEvalLists();

            if((rep_check = (!prog_chosen.isFather) && isSett(BOOLS_PROGREPEATCHECK) && (!prog_chosen.automatic)))
            {
                PRINTL();
                printf2(COLOR_CREDITS, "Press any key to repeat\nor press %c to go Back to Main Menu.\n", access(exit_char));
            }
            CLEARBUFFER();
        }
        while(rep_check && getch() != access(exit_char));
    }
    
    CLEARBUFFER();
    PRINTL();

    // #undef prog_chosen

    return;
}

// New Families Access Points
//
__MSSHELL_WRAPPER_ void __apnt algebraOperations(const register sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_ALGEBRA_OPERATIONS,
                        alg_operations, main_menu[MAIN_ALGEBRAOPERATIONS].name,
                        #if MAX_ALGEBRA_OPERATIONS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
                            BY_NUMBERS
                        #else
                            BY_CHARS
                        #endif // MAX_ALGEBRA_OPERATIONS
                        );
    return;
}

__MSSHELL_WRAPPER_ void __apnt changeProgramSettings(const register sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_SETTINGS, change_settings,
                        main_menu[MAIN_CHANGESETTINGS].name,
                        #if MAX_SETTINGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
                            BY_NUMBERS
                        #else
                            BY_CHARS
                        #endif // MAX_SETTINGS
                        );
    return;
}

//

__MSSHELL_WRAPPER_ __MSNATIVE_ void __system progInfo(sel_typ skip)
{
    if(skip)
    {
        FILE *fp = NULL;

        if((fp = fopen(DESCRIPTIONS_FOLDER"suite.NFO", "r")) == NULL)
            printErr(2, "Unable to open File containing\nProgram Description");
        else
        {
            printf2(COLOR_CREDITS, "\t          "PROG__NAME" V"PROG__VERSION" by ");
            printf2(access(colors)[MEMBER_COLORAUTHOR], PROG__AUTHOR);
            printf2(COLOR_CREDITS, ".\n\n\t   //________________________________________________________\\\\\
\n\nLAST UPDATE DATE: "PROG__LASTUPDATEDATE"\n");

            if(!readFile(DESCRIPTIONS_FOLDER"suite.NFO"))
                printErr(2, "Non-existent DESCRIPTION File in Directory "DESCRIPTIONS_FOLDER":\nsuite."DEFAULT_HELP_FILE_EXTENSION);

            PRINTL();
        }
    }
    else
    {
        printf2(COLOR_CREDITS, PROG__NAME" V"PROG__VERSION" by ");
        printf2(access(colors)[MEMBER_COLORAUTHOR], PROG__AUTHOR);
        printf2(COLOR_CREDITS, ".\n");
        PRINTL();
    }

    return;
}

__MSNATIVE_ static void _MS__private __system printOpersIdentifiers(void)
{
    printf2(COLOR_CREDITS, "\nTo use with expected INPUT format: a identifier b\n");
    printf2(COLOR_CREDITS, "Available Operations Lists:\n");

    SHOWPAUSEMESSAGE();
    PRINTL();

    dim_typ i, j;

    for(i=0; i<MAX_OPERATIONS; ++i)
    {
        printf("{\n\t- %s:", operazioni[i].name);
        for(j=0; j<MAX_ALIAS && strcmp(operazioni[i].identifiers[j], NULL_CHAR); ++j)
            printf("\n\t\ta %s %c %s\n", operazioni[i].identifiers[j], operazioni[i].bin_or_unary ? 'b':'0', operazioni[i].description);
        printf("}\n");
        if(catchPause()) return;
    }

    return;
}

__MSNATIVE_ inline bool _MS__private __system doesExistOperIdentifier(const char identifier[static MAX_IDENTIFIER_LENGTH], sel_typ output_ID[static MAX_DIMENSIONS])
{
    sel_typ i, j;
    for(i=0; i<MAX_OPERATIONS; ++i)
        for(j=0; j<MAX_ALIAS; ++j)
            if(!strcmp(identifier, operazioni[i].identifiers[j]))
            {
                output_ID[OPERATION_ID] = i;
                output_ID[ALIAS_ID] = j;
                return true;
            }

    return false;
}

__MSUTIL_ inline void __system __export freeExprEvalLists()
{
    exprValListFree(access(const_list));
    exprFuncListFree(access(func_list));
    return;
}

__MSUTIL_ void __system __export refreshExprEvalVarList(dim_typ which_env)
{
    volatile int err;
    jmp_buf jumper;

    time_t t1;
    double diff;

    diff = 0.00;

    exprType * const tmp = malloc(sizeof(exprType));

    errMem(tmp, VSPACE);

    tmp->var_list = INIT_VARLIST;
    ///    tmp->const_list = INIT_CONSTLIST;
    tmp->e_ANS = INIT_ANS;
    tmp->global_var = DEFAULT_GLOBALVAL;

    /* Create varlist */
    err = exprValListCreate(&(tmp->var_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Var List Creation Error\n");
        longjmp(jumper, err);
    }

    /* Init variable list */
    // err = exprValListAddAddress((*vlist), "global", &(suite.exprVars.global_var));
    err = exprValListAddAddress(tmp->var_list, "global", &(tmp->global_var));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Var List Init Error\n");
        longjmp(jumper, err);
    }

    // err = exprValListAdd((*vlist), DEFAULT_ENVS_ANSVALNAME, 0.0);
    err = exprValListAdd(tmp->var_list, DEFAULT_ENVS_ANSVALNAME, 0.00);
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Error adding variable \'%s\'\n", DEFAULT_ENVS_ANSVALNAME);
        longjmp(jumper, err);
    }

    // exprValListGetAddress((*vlist), DEFAULT_ENVS_ANSVALNAME, &(suite.exprVars.e_res));
    exprValListGetAddress(tmp->var_list, DEFAULT_ENVS_ANSVALNAME, &(tmp->e_ANS));
    if(tmp->e_ANS == NULL)
    {
        sprint("Unable to get address of \'%s\'\n", DEFAULT_ENVS_ANSVALNAME);
        longjmp(jumper, EXPR_ERROR_UNKNOWN);
    }

    FILE *fp = NULL;
    const bool assert = getItemsListNo(ENVS) != STARTING_ENVSNO;
    const bool asrt = isSett(BOOLS_SHOWDIFFTIME);

    if(assert && (fp = checkForFHErrors(listNo(which_env, ENVS)->path, "r")) == NULL)
        return;

    exprObj *e = INIT_OBJLIST;
    int start, end;

    if(assert)
    {
        // char c;
        char str[MAX_BUFSIZ]; // char str[MIN_STRING];

        // strcpy(str, NULL_CHAR);

        while(fgets(str, MAX_FILE_LINES, fp) != NULL)
        {
            // exprObj *e = INIT_OBJLIST;
            // int start, end;
            // err = setjmp(jumper);


            //if(e) // if(err && e)
            //  exprFree(e);

            /*err = */
            // exprCreate(&e, suite.exprVars.func_list, (*vlist), suite.exprVars.const_list, NULL, 0);
            err = exprCreate(&e, access(func_list), tmp->var_list, access(const_list), NULL, 0);
            if(err != EXPR_ERROR_NOERROR)
            {
                sprint("Expr Creation Error.\n");
                exprFree(e);
                continue;
            }


            err = exprParse(e, str);
            if(err != EXPR_ERROR_NOERROR)
            {
                exprGetErrorPosition(e, &start, &end);
                sprint("Parse Error (%d,%d).\n", start, end);
                exprFree(e);
                continue;
            }

            if(asrt)
                t1 = time(NULL);


            ityp val;
            err = exprEval(e, &val);

            if(err != EXPR_ERROR_NOERROR)
            {
                sprint("Eval Error: %d.\n", err);
                exprFree(e);
                continue;
            }

            if(asrt)
                diff += difftime(time(NULL), t1);

            exprFree(e);
        }

        fclose(fp);

        if(asrt)
        {
            PRINTL();
            printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", DEFAULT_PRECISION, (EXPRTYPE)diff);
            PRINTL();
        }
    }

    listNo(which_env, ENVS)->data = tmp;

    return;
}

__MSUTIL_ void __system __export refreshExprEvalLists()
{
    freeExprEvalLists();


    access(func_list) = INIT_FUNCLIST;
    access(const_list) = INIT_CONSTLIST;
    /* Set error buffer */

    volatile int err;
    jmp_buf jumper;

    err = setjmp(jumper);

    if(err)
        {
        /* Free stuff */

        if(access(func_list))
            exprFuncListFree(access(func_list));

        if(err != -1)
            sprint("Error: %d\n", err);

        return;
        }

    err = exprFuncListCreate(&(access(func_list)));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Func List Creation Error\n");
        longjmp(jumper, 1);
    }

    /* Init funclist */
    err = exprFuncListInit(access(func_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Error initializing internal functions\n");
        longjmp(jumper, err);
    }

    /* Create constlist */
    err = exprValListCreate(&access(const_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Const List Creation Error\n");
        longjmp(jumper, err);
    }

    /* Init constlist */
    err = exprValListInit(access(const_list));
    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Error initializing internal constants\n");
        longjmp(jumper, err);
    }

    return;
}

__MSNATIVE_ inline void __system __export setCurrentMatrix(dim_typ which_mat)
{
    if(getItemsListNo(MATRICES) != STARTING_MATNO && !extractMat(which_mat))
    {
        matrixObj * const tmp = ((matrixObj *)(listNo(which_mat, MATRICES)->data));
        matrixFree(&(tmp->matrix), tmp->dim[RAWS]);
    }

    return;
}


__MSSHELL_WRAPPER_ __MSNATIVE_ void __system setDefaults()
{

    static bool once_executed = false;
    access(mode) = PROGRAM_BUSY;

    if(once_executed)
        resetProgramSettings(access(curLayout), listNo(access(lists)[LAYOUTS].cur_item, LAYOUTS)->path);

    once_executed = true;

    access(exit_char) = EXIT_CHAR;

    randomize;
    access(random_seed) = starting_random_seed;

    dim_typ i;

	#pragma omp parallel for
    for(i=0; i<MAX_OPERATIONS; ++i)
        operazioni[i] = default_operazioni[i];

    refreshExprEvalLists();

    if(access(lmpMatrix))
    {
        if(access(lmpMatrix)->matrix)
            matrixFree(&(access(lmpMatrix)->matrix), access(lmpMatrix)->dim[RAWS]);
        free(access(lmpMatrix));
    }
    
    access(lmpMatrix) = malloc(sizeof(matrixObj));
    errMem(access(lmpMatrix), VSPACE);

	resetLmpMatrix();

    return;
}
