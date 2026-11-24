#ifndef WRAPPER_HERMITE_RULE_H_INCLUDED
#define WRAPPER_HERMITE_RULE_H_INCLUDED

#define cdgqf(a,b,c,d,e,f) FUNCNAME_CDGQF(C_2DT2IT2PIT(a,b,c,d,e,f))
#define cgqf(a,b,c,d,e,f,g,h) FUNCNAME_CGQF(C_2DT4IT2PIT(a,b,c,d,e,f,g,h))
#define class_matrix(a,b,c,d,e,f) R_DBL(FUNCNAME_CLASSMATRIX(C_2DT2IT2PIT(a,b,c,d,e,f)))
#define scqf(a,b,c,d,e,f,g,h,i,j,k,l,m) FUNCNAME_SCQF(C_DTPITPIIPITPI2PITDT4IT(a,b,c,d,e,f,g,h,i,j,k,l,m))
#define sgqf(a,b,c,d,e,f) FUNCNAME_SGQF(C_DT2PITIT2PIT(a,b,c,d,e,f))

__MATHSUITE __JBURKARDT void * FUNCNAME_CDGQF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CGQF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CLASSMATRIX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SCQF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SGQF(void *);

#endif
