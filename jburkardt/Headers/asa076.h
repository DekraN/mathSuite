#ifndef WRAPPER_ASA076_H_INCLUDED
#define WRAPPER_ASA076_H_INCLUDED

#define owen_values(a,b,c,d) FUNCNAME_OWENVALUES(C_PDT3PIT(a,b,c,d))
#define tha(a,b,c,d,e) R_DBL(FUNCNAME_THA(C_PDBL5(a,b,c,d,e)))


__MATHSUITE __JBURKARDT void * FUNCNAME_THA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_OWENVALUES(void *);

#endif
