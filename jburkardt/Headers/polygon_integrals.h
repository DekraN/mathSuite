#ifndef WRAPPER_POLYGON_INTEGRALS_H_INCLUDED
#define WRAPPER_POLYGON_INTEGRALS_H_INCLUDED

#define moment(a,b,c,d,e) R_DBL(FUNCNAME_MOMENT(C_3DT2PIT(a,d,e,b,c)))
#define moment_central(a,b,c,d,e) R_DBL(FUNCNAME_MOMENTCENTRAL(C_3DT2PIT(a,d,e,b,c)))
#define moment_normalized(a,b,c,d,e) R_DBL(FUNCNAME_MOMENTNORMALIZED(C_3DT2PIT(a,d,e,b,c)))

__MATHSUITE __JBURKARDT void * FUNCNAME_MOMENT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MOMENTCENTRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MOMENTNORMALIZED(void *);

#endif
