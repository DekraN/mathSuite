#ifndef WRAPPER_TRIANGLE_INTEGRALS_H_INCLUDED
#define WRAPPER_TRIANGLE_INTEGRALS_H_INCLUDED

#define triangle01_monomial_integral_(a,b) R_DBL(FUNCNAME_TRIANGLE01MONOMIALINTEGRAL_(C_PINT2(a,b)))
#define triangle01_poly_integral(a,b) R_DBL(FUNCNAME_TRIANGLE01POLYINTEGRAL(C_DTPIT(a,b)))
#define poly_product(a,b,c,d) FUNCNAME_POLYPRODUCT(C_2DT2PIT(a,c,b,d))
#define poly_power_linear(a,b,c) FUNCNAME_POLYPOWERLINEAR(C_2DTPIT(a,c,b))
#define rs_to_xy_map(a,b,c,d,e,f,g) FUNCNAME_RSTOXYMAP(C_PPDBL7(a,b,c,d,e,f,g))
#define triangle_area_from_vertex(a) R_DBL(FUNCNAME_TRIANGLEAREAFROMVERTEX(a))
#define triangle_monomial_integral(a,b,c) R_DBL(FUNCNAME_TRIANGLEMONOMIALINTEGRAL(C_2DTPIT(a,b,c)))
#define triangle_poly_integral(a,b,c) R_DBL(FUNCNAME_TRIANGLEPOLYINTEGRAL(C_DT2PIT(a,b,c)))
#define triangle_xy_integral(a,b,c,d,e,f) FUNCNAME_TRIANGLEXYINTEGRAL(C_PDBL6(a,b,c,d,e,f))
#define xy_to_rs_map(a,b,c,d,e,f,g) FUNCNAME_XYTORSMAP(C_PPDBL7(a,b,c,d,e,f,g))

__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLE01MONOMIALINTEGRAL_(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLE01POLYINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RSTOXYMAP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEAREAFROMVERTEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEMONOMIALINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEPOLYINTEGRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_XYTORSMAP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYPRODUCT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POLYPOWERLINEAR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_TRIANGLEXYINTEGRAL(void *);

#endif
