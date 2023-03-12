#ifndef WRAPPER_ASA103_H_INCLUDED
#define WRAPPER_ASA103_H_INCLUDED

#define digama(a) R_DBL(FUNCNAME_DIGAMA(C_SDBL(a)))
#define psi_values(a,b,c) FUNCNAME_PSIVALUES(C_PDT2PIT(a,b,c))

__MATHSUITE __JBURKARDT void * FUNCNAME_DIGAMA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PSIVALUES(void *);

#endif
