/*
    File: exprincl.h
    Auth: Brian Allen Vanderburg II
    Date: Thursday, April 24, 2003
    Desc: Includes, macros, etc needed by this library

    This file is part of ExprEval.
*/

#ifndef __BAVII_EXPRINCL_H
#define __BAVII_EXPRINCL_H


/* Includes and macros and whatnot for building the library */

/* Memory routines.  memory.h for VC++, mem.h for BC++ */
#ifdef __TURBOC__
#include <mem.h>
#else
#include <memory.h>
#endif

/* Memory allocation */
#include <malloc.h>

/* String routines */
#include <string.h>

/* Character manipulation routines */
#include <ctype.h>

/* Standard routines */
#include <stdlib.h>

/* Math routines */
#include <math.h>

/* Time */
#include <time.h>


/* Math constants.  VC++ does not seem to have these */
#ifndef M_E
#define M_E 2.7182818284590452354
#endif

#ifndef M_LOG2E
#define M_LOG2E 1.4426950408889634074
#endif

#ifndef M_LOG10E
#define M_LOG10E 0.43429448190325182765
#endif

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_4
#define M_PI_4 0.78539816339744830962
#endif

#ifndef M_1_PI
#define M_1_PI 0.31830988618379067154
#endif

#ifndef M_2_PI
#define M_2_PI 0.63661977236758134308
#endif

#ifndef M_1_SQRTPI
#define M_1_SQRTPI 0.56418958354776
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_1_SQRT2
#define M_1_SQRT2 0.70710678118654752440
#endif


// Physics Constants
#ifndef F_G
#define F_G 6.67259E-11
#endif

#ifndef F_GA
#define F_GA 9.80665
#endif

#ifndef F_LA
#define F_LA 1.6300000
#endif

#ifndef F_MRS
#define F_MRS 3.7200000
#endif

#ifndef F_KB
#define F_KB 1.3806503E-23
#endif

#ifndef F_R
#define F_R 8.31400000
#endif

#ifndef F_C
#define F_C 299792458
#endif

#ifndef F_EPS0
#define F_EPS0 8.854187817E-12
#endif

#ifndef F_MI0
#define F_MI0 M_PI*4E-7
#endif

#ifndef F_P
#define F_P 6.62606876E-34
#endif

#ifndef F_EC
#define F_EC 1.602176462E-19
#endif

#ifndef F_ESM
#define F_ESM 9.10938188E-31
#endif

#ifndef F_PSM
#define F_PSM 1.67262158E-27
#endif

#ifndef F_NSM
#define F_NSM 1.67492716E-27
#endif

#ifndef F_AMU
#define F_AMU 1.66053873E-27
#endif

#ifndef F_AV
#define F_AV 6.02214199E23
#endif

#ifndef F_F
#define F_F 9.648534150000
#endif

#ifndef F_ALPHA
#define F_ALPHA 0.007297352533
#endif

#ifndef F_BHR
#define F_BHR 5.291772083E-11
#endif

#ifndef F_RDB
#define F_RDB 1.0973731568549E7
#endif

#ifndef F_MBHR
#define F_MBHR 9.27400899E-24
#endif

#ifndef F_IV
#define F_IV 22.710981
#endif

#ifndef F_EH
#define F_EH 4.35974381E-18
#endif

#ifndef F_EMM
#define F_EMM -9.28476362E-24
#endif

#ifndef F_EMP
#define F_EMP 1.41060761E-26
#endif

#ifndef F_NM
#define F_NM 5.0507866E-27
#endif

#ifndef F_PGMR
#define F_PGMR 2.6752212800000000
#endif

#ifndef F_SB
#define F_SB 5.670400E-8
#endif

#ifndef F_FR
#define F_FR 3.7417749E-16
#endif

#ifndef F_SR
#define F_SR 0.1438769
#endif

// Informatic Constants
#ifndef M_KIB
#define M_KIB 1024
#endif

#ifndef M_MIB
#define M_MIB 1048576
#endif

#ifndef M_GIB
#define M_GIB 1073741824
#endif

#ifndef M_TIB
#define M_TIB 1099511627776
#endif

#ifndef M_PIB
#define M_PIB 1125899906842624
#endif

#ifndef M_EIB
#define M_EIB 1152921504606846976
#endif


// Time constants
#ifndef T_M
#define T_M 60
#endif

#ifndef T_H
#define T_H 3600
#endif

#ifndef T_D
#define T_D 86400
#endif

#ifndef T_W
#define T_W 604800
#endif

#ifndef T_MT
#define T_MT 18144000
#endif

#ifndef T_M1
#define T_M1 18748800
#endif

#ifndef T_BM
#define T_BM 16934400
#endif

#ifndef T_BM1
#define T_BM1 17539200
#endif

#ifndef T_Y
#define T_Y 217728000
#endif

#ifndef T_Y1
#define T_Y1 224985600
#endif

#ifndef T_BY
#define T_BY 203212800
#endif

#ifndef T_BY1
#define T_BY1 210470400
#endif


#endif /* __BAVII_EXPRINCL_H */
