#ifndef WRAPPER_ISBN_H_INCLUDED
#define WRAPPER_ISBN_H_INCLUDED

#define ch_is_digit(a) R_UCHR(FUNCNAME_CHISDIGIT(C_SCHR(a)))
#define ch_to_digit(a) R_INT(FUNCNAME_CHTODIGIT(C_SCHR(a)))
#define s_to_digits(a,b) FUNCNAME_STODIGITS(C_PCDT(a,b))

__MATHSUITE __JBURKARDT void * FUNCNAME_STODIGITS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHISDIGIT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHTODIGIT(void *);

#endif
