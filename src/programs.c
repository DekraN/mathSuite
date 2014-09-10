// programs.c 10/09/2014 Marco Chiarelli aka DekraN
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
        "Base 10 Plus 1 Logarithm",
        NULL_CHAR,
        {
            IDENTIFIER_LOGARITMO101P
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Base 10 Plus 1 Cardinal Logarithm",
        NULL_CHAR,
        {
            IDENTIFIER_LOGARITMO101PC
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Base 2 Plus 1 Logarithm",
        NULL_CHAR,
        {
            IDENTIFIER_LOGARITMO21P
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Base 2 Plus 1 Cardinal Logarithm",
        NULL_CHAR,
        {
            IDENTIFIER_LOGARITMO21PC
        },
        DOMAIN_DEFAULT,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Complex Exponential",
        "(a = RPART, b = IPART)",
        {
            OPERATOR_CEXP,
            IDENTIFIER_CEXP,
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function,
    },
    {
        "Complex Cardinal Exponential",
        "(a = RPART, b = IPART)",
        {
            OPERATOR_CEXP"c",
            IDENTIFIER_CEXP"c",
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function,
    },
    {
        "Complex 10 Base Exponential",
        "(a = RPART, b = IPART)",
        {
            OPERATOR_CEXP10,
            OPERATOR_ALIAS_CEXP10,
            IDENTIFIER_CEXP10
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex 10 Base Cardinal Exponential",
        "(a = RPART, b = IPART)",
        {
            OPERATOR_CEXP10"c",
            OPERATOR_ALIAS_CEXP10"c",
            IDENTIFIER_CEXP10"c"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex 2 Base Exponential",
        "(a = RPART, b = IPART)",
        {
            OPERATOR_CEXP2,
            OPERATOR_ALIAS_CEXP2,
            IDENTIFIER_CEXP2
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex 2 Base Cardinal Exponential",
        "(a = RPART, b = IPART)",
        {
            OPERATOR_CEXP2"c",
            OPERATOR_ALIAS_CEXP2"c",
            IDENTIFIER_CEXP2"c"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Square Root",
        NULL_CHAR,
        {
            IDENTIFIER_CSQRT
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Cubic Root",
        NULL_CHAR,
        {
            IDENTIFIER_CCBRT
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG,
            IDENTIFIER_ALIAS_CLOG,
            IDENTIFIER_ALIAS2_CLOG,
            IDENTIFIER_ALIAS3_CLOG
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Cardinal Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG"c",
            IDENTIFIER_ALIAS_CLOG"c",
            IDENTIFIER_ALIAS2_CLOG"c",
            IDENTIFIER_ALIAS3_CLOG"c"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Base 10 Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG10,
            IDENTIFIER_ALIAS_CLOG10,
            IDENTIFIER_ALIAS2_CLOG10,
            IDENTIFIER_ALIAS3_CLOG10
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Base 10 Cardinal Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG10"c",
            IDENTIFIER_ALIAS_CLOG10"c",
            IDENTIFIER_ALIAS2_CLOG10"c",
            IDENTIFIER_ALIAS3_CLOG10"c"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Base 2 Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG2,
            IDENTIFIER_ALIAS_CLOG2,
            IDENTIFIER_ALIAS2_CLOG2,
            IDENTIFIER_ALIAS3_CLOG2
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Base 2 Cardinal Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG2"c",
            IDENTIFIER_ALIAS_CLOG2"c",
            IDENTIFIER_ALIAS2_CLOG2"c",
            IDENTIFIER_ALIAS3_CLOG2"c"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Plus 1 Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG1P
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Cardinal Plus 1 Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG1P"c"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Base 10 Plus 1 Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG101P
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Base 10 Cardinal Plus 1 Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG101P"c"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Base 2 Plus 1 Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG21P
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Base 2 Cardinal Plus 1 Logarithm",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CLOG21P"c"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
    	"Complex Argument"
    	"(a = RPART, b = IPART)",
    	{
    		IDENTIFIER_CARG
    	},
    	DOMAIN_DEFAULT,
    	DOMAIN_DEFAULT,
    	binary_function

    },
    {
    	"Complex Absolute Value"
    	"(a = RPART, b = IPART)",
    	{
    		IDENTIFIER_CABS
    	},
    	DOMAIN_DEFAULT,
    	DOMAIN_DEFAULT,
    	binary_function
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
		"Change Exit Char",
		NULL_CHAR,
		{
			IDENTIFIER_EXITCHAR,
			IDENTIFIER_ALIAS_EXITCHAR,
			IDENTIFIER_ALIAS2_EXITCHAR,
			IDENTIFIER_ALIAS3_EXITCHAR
		},
		DOMAIN_SCHAR,
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
    	"Change Min Stirling Number",
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
    	"Change Outlier Constant",
    	NULL_CHAR,
    	{
    		IDENTIFIER_OUTLIERCONST,
    		IDENTIFIER_ALIAS_OUTLIERCONST,
    		IDENTIFIER_ALIAS2_OUTLIERCONST
    	},
    	DOMAIN_FLT,
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
        "Complex Sin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CSIN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Sin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CSIN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Cos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Cos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Tan",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CTAN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Tan",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CTAN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Csc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCSC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Csc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCSC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Sec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CSEC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Sec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CSEC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Cot",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCOT
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Cot",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCOT"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHSIN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Hsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHSIN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQSIN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQSIN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Hcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hsec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHSEC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hypernolic Hsec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHSEC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qsec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQSEC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qsec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQSEC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hcsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCSC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Hcsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCSC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qcsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCSC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qcsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCSC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Htan",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHTAN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Htan",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHTAN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qtan",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQTAN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qtan",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQTAN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hcot",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCOT
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Hcot",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCOT"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qcot",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCOT
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qcot",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCOT"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Vsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CVSIN,
            IDENTIFIER_ALIAS_CVSIN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Vsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CVSIN"h",
            IDENTIFIER_ALIAS_CVSIN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex CVsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCVSIN,
            IDENTIFIER_ALIAS_CCVSIN,
            IDENTIFIER_ALIAS2_CCVSIN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic CVsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCVSIN"h",
            IDENTIFIER_ALIAS_CCVSIN"h",
            IDENTIFIER_ALIAS2_CCVSIN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Vcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CVCOS,
            IDENTIFIER_ALIAS_CVCOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Vcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CVCOS"h",
            IDENTIFIER_ALIAS_CVCOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex CVcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCVCOS,
            IDENTIFIER_ALIAS_CCVCOS,
            IDENTIFIER_ALIAS2_CCVCOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic CVcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCVCOS"h",
            IDENTIFIER_ALIAS_CCVCOS"h",
            IDENTIFIER_ALIAS2_CCVCOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex HVsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHVSIN,
            IDENTIFIER_ALIAS_CHVSIN,
            IDENTIFIER_ALIAS2_CHVSIN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic HVsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHVSIN"h",
            IDENTIFIER_ALIAS_CHVSIN"h",
            IDENTIFIER_ALIAS2_CHVSIN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex HCVsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCVSIN,
            IDENTIFIER_ALIAS_CHCVSIN,
            IDENTIFIER_ALIAS2_CHCVSIN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic HCVsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCVSIN"h",
            IDENTIFIER_ALIAS_CHCVSIN"h",
            IDENTIFIER_ALIAS2_CHCVSIN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex QVsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQVSIN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic QVsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQVSIN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex QCVsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic QCVsin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex HVcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHVCOS,
            IDENTIFIER_ALIAS_CHVCOS,
            IDENTIFIER_ALIAS2_CHVCOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic HVcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHVCOS"h",
            IDENTIFIER_ALIAS_CHVCOS"h",
            IDENTIFIER_ALIAS2_CHVCOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex HCVcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCVCOS,
            IDENTIFIER_ALIAS_CHCVCOS,
            IDENTIFIER_ALIAS2_CHCVCOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic HCVcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCVCOS"h",
            IDENTIFIER_ALIAS_CHCVCOS"h",
            IDENTIFIER_ALIAS2_CHCVCOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex QVcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQVCOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic QVcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQVCOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex QCVcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCVCOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic QCVcos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCVCOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Esec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CESEC,
            IDENTIFIER_ALIAS_CESEC,
            IDENTIFIER_ALIAS2_CESEC,
            IDENTIFIER_ALIAS3_CESEC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Esec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CESEC"h",
            IDENTIFIER_ALIAS_CESEC"h",
            IDENTIFIER_ALIAS2_CESEC"h",
            IDENTIFIER_ALIAS3_CESEC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Ecsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CECSC,
            IDENTIFIER_ALIAS_CECSC,
            IDENTIFIER_ALIAS2_CECSC,
            IDENTIFIER_ALIAS3_CECSC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Ecsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CECSC"h",
            IDENTIFIER_ALIAS_CECSC"h",
            IDENTIFIER_ALIAS2_CECSC"h",
            IDENTIFIER_ALIAS3_CECSC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex HEsec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHESEC,
            IDENTIFIER_ALIAS_CHESEC,
            IDENTIFIER_ALIAS2_CHESEC,
            IDENTIFIER_ALIAS3_CHESEC,
            IDENTIFIER_ALIAS4_CHESEC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic HEsec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHESEC"h",
            IDENTIFIER_ALIAS_CHESEC"h",
            IDENTIFIER_ALIAS2_CHESEC"h",
            IDENTIFIER_ALIAS3_CHESEC"h",
            IDENTIFIER_ALIAS4_CHESEC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex HEcsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHECSC,
            IDENTIFIER_ALIAS_CHECSC,
            IDENTIFIER_ALIAS2_CHECSC,
            IDENTIFIER_ALIAS3_CHECSC,
            IDENTIFIER_ALIAS4_CHECSC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic HEcsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHECSC"h",
            IDENTIFIER_ALIAS_CHECSC"h",
            IDENTIFIER_ALIAS2_CHECSC"h",
            IDENTIFIER_ALIAS3_CHECSC"h",
            IDENTIFIER_ALIAS4_CHECSC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex QEsec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQESEC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic QEsec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQESEC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex QEcsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQECSC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic QEcsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQECSC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex sinc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CSINC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic sinc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CSINC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hsinc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHSINC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Hsinc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHSINC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qsinc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQSINC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qsinc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQSINC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex cosc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCOSC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic cosc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCOSC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hcosc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCOSC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Hcosc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCOSC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qcosc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCOSC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qcosc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCOSC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex secc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CSECC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic secc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CSECC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hsecc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHSECC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Hsecc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHSECC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qsecc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQSECC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qsecc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQSECC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex cscc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCSCC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic cscc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCSCC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hcscc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCSCC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Hcscc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCSCC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qcscc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCSCC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qcscc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCSCC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex tanc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CTANC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic tanc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CTANC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Htanc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHTANC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Htanc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHTANC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qtanc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQTANC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qtanc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQTANC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex cotc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCOTC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic cotc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CCOTC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hcotc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCOTC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Hcotc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CHCOTC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Qcotc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCOTC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Qcotc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CQCOTC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Asin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CASIN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Asin",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CASIN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Acos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CACOS
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Acos",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CACOS"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Atan",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CATAN
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Atan",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CATAN"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Acsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CACSC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Acsc",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CACSC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Asec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CASEC
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Asec",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CASEC"h"
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Acot",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CACOT
        },
        DOMAIN_DEFAULT,
        DOMAIN_DEFAULT,
        binary_function
    },
    {
        "Complex Hyperbolic Acot",
        "(a = RPART, b = IPART)",
        {
            IDENTIFIER_CACOT"h"
        },
        DOMAIN_DEFAULT,
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
    	"(b = Outlier Index, Enter different value from the 1st to Stop Entering)",
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
    	"Map",
    	"(b = FID, Enter different value from the 1st to Stop Entering)",
    	{
    		IDENTIFIER_MAP,
    		IDENTIFIER_ALIAS_MAP,
    		IDENTIFIER_ALIAS2_MAP,
    		IDENTIFIER_ALIAS3_MAP,
    	},
    	DOMAIN_DEFAULT,
		DOMAIN_UCHAR,
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
            IDENTIFIER_PERMUTATIONS
        },
        DOMAIN_ULLONG,
        DOMAIN_NULL,
        unary_function
    },
    {
        "Permutations with Repetitions",
    	"(Enter b as the i-th Item of the Repetition Vector. Enter 0 to Stop)",
        {
            IDENTIFIER_PERMUTATIONSREP,
            IDENTIFIER_ALIAS_PERMUTATIONSREP,
    		IDENTIFIER_ALIAS2_PERMUTATIONSREP,
    		IDENTIFIER_ALIAS3_PERMUTATIONSREP
        },
        DOMAIN_ULLONG,
        DOMAIN_ULLONG,
        binary_function
    },
    {
        "K-Permutations",
        NULL_CHAR,
        {
            IDENTIFIER_KPERMUTATIONS
        },
        DOMAIN_ULLONG,
        DOMAIN_ULLONG,
        binary_function
    },
    {
    	"K-Permutations with Repetitions",
		NULL_CHAR,
    	{
    		IDENTIFIER_KPERMUTATIONSREP,
    		IDENTIFIER_ALIAS_KPERMUTATIONSREP,
    		IDENTIFIER_ALIAS2_KPERMUTATIONSREP,
    		IDENTIFIER_ALIAS3_KPERMUTATIONSREP
    	},
    	DOMAIN_ULLONG,
    	DOMAIN_ULLONG,
    	binary_function
    },
    {
        "Combinations",
        NULL_CHAR,
        {
            IDENTIFIER_COMBINATIONS
        },
        DOMAIN_ULLONG,
        DOMAIN_ULLONG,
        binary_function
    },
    {
        "Combinations with Repetitions",
        NULL_CHAR,
        {
            IDENTIFIER_COMBINATIONSREP,
            IDENTIFIER_ALIAS_COMBINATIONSREP,
            IDENTIFIER_ALIAS2_COMBINATIONSREP,
            IDENTIFIER_ALIAS3_COMBINATIONSREP
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

__MSNATIVE_ __MSSHELL_WRAPPER_ static void viewComplexResult(const bool checkdg, const char oprname[static DINFO_STRING], register double complex result);

__MSNATIVE_ __MSSHELL_WRAPPER_ static inline void __system viewComplexResult(const bool checkdg, const char oprname[static DINFO_STRING], register double complex result)
{																				 								 		
	if(checkdg && isSett(BOOLS_DEGREESENTERING)) result = deg(creal(result)) + deg(cimag(result))*I;				   
	printf2(COLOR_USER, "\nRESULT of %s operation is: ", oprname);	
	printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, creal(result));					
	printf2(COLOR_USER, " + ");														
	printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, cimag(result));					
	printf2(COLOR_USER, "*i;\n\n");													
	return;
}

__MSSHELL_WRAPPER_ void basicCalculator(const register sel_typ argc, char ** argv)
{

    if(isSett(BOOLS_BASECALCPARSING))
        // requires(NULL, yes to show result, yes to show varlist, yes to show difftime, yes to save results, not to reset lists
        requires(argc ? argv[0] : NULL, argc ? NULL : "Enter an Expression", "Result is", PARSER_SHOWRESULT | PARSER_SHOWVARLIST | PARSER_SHOWDIFFTIME | PARSER_SAVERESULT);
    else
    {

        #define viewInsertedValue                                                                                \
            printf2(COLOR_SYSTEM, "\nValue: ");                                                                  \
            printf2(COLOR_CREDITS, OUTPUT_CONVERSION_FORMAT, a);                                                 \
            printf2(COLOR_SYSTEM, " inserted in [%llu] position of elaborating Media Vector...\n\n", accumulate);\
            return;

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
        fsel_typ output_ID[MAX_DIMENSIONS] =
        {
            MAX_OPERATIONS,
            0
        };

        c = 0;
        CLEARBUFFER();

        doesExistOperIdentifier(identifier, output_ID);

        const fsel_typ oprID = output_ID[OPERATION_ID];

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
    	const bool complex_entries = ((oprID >= BCALC_CSIN && oprID <= BCALC_CACOTH) || (oprID >= BCALC_CEXP && oprID <= BCALC_CLOG21PC) || oprID == BCALC_CSQRT || oprID == BCALC_CCBRT);

		if(dg_ent)
		{
			if(complex_entries)
				a = rad(a);
       	 	if((complex_entries || (oprID > BCALC_TRASFORMAANGOLI && oprID < BCALC_ASINANDASINH)) && a == 0)
				b = rad(b);
		}
		
        switch(oprID)
        {

            #include "bcilfs.h"

            default:
                printErr(1, "Non-existent Operation");
                return ;
        }

        #undef viewInsertedValue


        if(dg_ent && oprID > BCALC_QECSCANDQECSCH && oprID < BCALC_MCM && a == 0) c = deg(c);

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

__MSSHELL_WRAPPER_ void __apnt mssManager(const register sel_typ argc, char ** argv)
{
    operationsGroupMenu(MAX_MSSMANAGER_PROGS, mss_manager,
                        main_menu[MAIN_MSSMANAGER].name,
                        #if MAX_MSSMANAGER_PROGS > MAX_CASEINSENSITIVE_CHARS_ALPHABET
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
            if(prog_chosen.program_function == basicCalculator)
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

        bool rep_check;

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
                printf2(COLOR_CREDITS, "Press any key to repeat\nor press %c to go Back to Main Menu.\n", access(curLayout)->exit_char);
            }
            CLEARBUFFER();
        }
        while(rep_check && getch() != access(curLayout)->exit_char);
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

__MSNATIVE_ inline bool _MS__private __system doesExistOperIdentifier(const char identifier[static MAX_IDENTIFIER_LENGTH], fsel_typ output_ID[static MAX_DIMENSIONS])
{

    for(fsel_typ i=0; i<MAX_OPERATIONS; ++i)
        for(fsel_typ j=0; j<MAX_ALIAS; ++j)
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
    int err;
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

    int err;
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
        matrixFree(&(tmp->matrix));
    }

    return;
}


__MSSHELL_WRAPPER_ __MSNATIVE_ void __system setDefaults()
{

    static bool once_executed = false;
    access(mode) = PROGRAM_BUSY;

    if(once_executed)
        resetProgramSettings(access(curLayout), listNo(access(lists)[LAYOUTS].cur_item, LAYOUTS)->path);
    else
    	once_executed = true;
    	
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
            matrixFree(&(access(lmpMatrix)->matrix));
        free(access(lmpMatrix));
    }
    
    access(lmpMatrix) = malloc(sizeof(matrixObj));
    errMem(access(lmpMatrix), VSPACE);

	resetLmpMatrix();

    return;
}
