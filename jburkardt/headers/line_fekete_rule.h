#ifndef WRAPPER_LINE_FEKETE_RULE_H_INCLUDED
#define WRAPPER_LINE_FEKETE_RULE_H_INCLUDED

#define cheby_van1(a,b,c,d,e) FUNCNAME_CHEBVYVAN1(C_DT2ITDTPIT(a,b,c,d,e))
#define legendre_van(a,b,c,d,e) FUNCNAME_LEGENDREVAN(C_DT2ITDTPIT(a,b,c,d,e))
#define line_fekete_chebyshev(a,b,c,d,e,f,g,h) FUNCNAME_LINEFEKETECHEBYSHEV(C_DT2ITDTPITPDT2PIT(a,b,c,d,e,f,g,h))
#define line_fekete_legendre(a,b,c,d,e,f,g,h) FUNCNAME_LINEFEKETELEGENDRE(C_DT2ITDTPITPDT2PIT(a,b,c,d,e,f,g,h))
#define line_fekete_monomial(a,b,c,d,e,f,g,h) FUNCNAME_LINEFEKETEMONOMIAL(C_DT2ITDTPITPDT2PIT(a,b,c,d,e,f,g,h))
#define line_monomial_moments(a,b,c) FUNCNAME_LINEMONOMIALMOMENTS(C_2ITDT(a,b,c))

__MATHSUITE __JBURKARDT void * FUNCNAME_LINEFEKETECHEBYSHEV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEFEKETELEGENDRE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEFEKETEMONOMIAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHEBVYVAN1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LEGENDREVAN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LINEMONOMIALMOMENTS(void *);

#endif
