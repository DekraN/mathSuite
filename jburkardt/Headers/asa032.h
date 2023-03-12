#ifndef WRAPPER_ASA032_H_INCLUDED
#define WRAPPER_ASA032_H_INCLUDED

#define gamma_inc(a,b,c,d,e) FUNCNAME_GAMMAINC(C_PDT4PIT(e,a,b,c,d))
#define gamma_inc_values(a,b,c,d) FUNCNAME_GAMMAINCVALUES(C_2PDT2PIT(a,b,c,d))

__MATHSUITE __JBURKARDT void * FUNCNAME_GAMMAINCVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAMMAINC(void *);

#endif
