#ifndef WRAPPER_LAGRANGE_APPROX_1D_H_INCLUDED
#define WRAPPER_LAGRANGE_APPROX_1D_H_INCLUDED

#define lagrange_approx_1d(a,b,c,d,e,f) FUNCNAME_LAGRANGEAPPROX1D(C_3DT3PIT(a,b,e,c,d,f))
#define lagrange_basis_1d(a,b,c,d) FUNCNAME_LAGRANGEBASIS1D(C_2DT2PIT(a,c,b,d))

__MATHSUITE __JBURKARDT void * FUNCNAME_LAGRANGEAPPROX1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LAGRANGEBASIS1D(void *);

#endif