#ifndef WRAPPER_LOBATTO_POLYNOMIAL_H_INCLUDED
#define WRAPPER_LOBATTO_POLYNOMIAL_H_INCLUDED

#define lobatto_polynomial_value(a,b) R_DBL(FUNCNAME_LOBATTOPOLYNOMIALVALUE(C_DTIT(a,b)))
#define lobatto_polynomial_derivative(a,b) R_DBL(FUNCNAME_LOBATTOPOLYNOMIALDERIVATIVE(C_DTIT(a,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_LOBATTOPOLYNOMIALVALUE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LOBATTOPOLYNOMIALDERIVATIVE(void *);

#endif