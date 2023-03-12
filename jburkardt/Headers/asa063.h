#ifndef WRAPPER_ASA063_H_INCLUDED
#define WRAPPER_ASA063_H_INCLUDED

#define betain(a,b,c,d) R_DBL(FUNCNAME_BETAIN(C_PDBL4(a,b,c,d)))
#define beta_inc_values(a,b,c,d,e) FUNCNAME_BETAINCVALUES(C_PDT4PIT(a,b,c,d,e))

__MATHSUITE __JBURKARDT void * FUNCNAME_BETAIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETAINCVALUES(void *);

#endif
