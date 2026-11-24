#ifndef WRAPPER_FEM_BASIS_H_INCLUDED
#define WRAPPER_FEM_BASIS_H_INCLUDED

#define fem_basis_1d(a,b,c) R_DBL(FUNCNAME_FEMBASIS1D(C_2DTIT(a,b,c)))
#define fem_basis_2d(a,b,c,d,e) R_DBL(FUNCNAME_FEMBASIS2D(C_3DT2IT(a,b,c,d,e)))
#define fem_basis_3d(a,b,c,d,e,f,g) R_DBL(FUNCNAME_FEMBASIS3D(C_3DT2IT(a,b,c,d,e,f,g)))
#define fem_basis_md(a,b,c) R_DBL(FUNCNAME_FEMBASISMD(C_DTPITPI(a,c,b)))
#define fem_basis_prism_triangle(a,b,c) R_DBL(FUNCNAME_FEMBASISPRISMTRIANGLE(C_2PDTPIT(a,b,c)))
#define r8_fraction(a,b) R_DBL(FUNCNAME_R8FRACTION(C_PINT2(a,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_FEMBASISMD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FEMBASISPRISMTRIANGLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8FRACTION(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FEMBASIS1D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FEMBASIS2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FEMBASIS3D(void *);

#endif
