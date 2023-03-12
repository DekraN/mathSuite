#ifndef WRAPPER_ASA226_H_INCLUDED
#define WRAPPER_ASA226_H_INCLUDED

#define betanc(a,b,c,d) R_DBL(FUNCNAME_BETANC(C_PDBL4(a,b,c,d)))
#define beta_noncentral_cdf_values(a,b,c,d,e,f) FUNCNAME_BETANONCENTRALCDFVALUES((C_DT5PIT(a,b,c,d,e,f))

__MATHSUITE __JBURKARDT void * FUNCNAME_BETANC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETANONCENTRALCDFVALUES(void *);

#endif
