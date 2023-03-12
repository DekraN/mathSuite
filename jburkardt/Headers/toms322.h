#ifndef WRAPPER_TOMS322_H_INCLUDED
#define WRAPPER_TOMS322_H_INCLUDED

#define fisher(a,b,c) R_DBL(FUNCNAME_FISHER(C_2DTIT(a,b,c)))
#define student(a,b) R_DBL(FUNCNAME_STUDENT(C_DTIT(a,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_STUDENT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FISHER(void *);

#endif
