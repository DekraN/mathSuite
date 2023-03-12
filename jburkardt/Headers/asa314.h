#ifndef WRAPPER_ASA314_H_INCLUDED
#define WRAPPER_ASA314_H_INCLUDED

#define invmod(a,b,c,d,e) R_UCHR(FUNCNAME_INVMOD(C_DT2PIT2PDT(a,b,c,d,e)))
#define msort(a,b,c,d,e,f,g) FUNCNAME_MSORT(C_DT2PIT2PDT2PIT(a,b,c,d,e,f,g))
#define musort(a,b,c,d,e,f,g) FUNCNAME_MUSORT(C_DT2PIT2PDT2PIT(a,b,c,d,e,f,g))


__MATHSUITE __JBURKARDT void * FUNCNAME_INVMOD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MSORT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MUSORT(void *);

#endif
