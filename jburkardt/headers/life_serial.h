#ifndef WRAPPER_LIFE_SERIAL_H_INCLUDED
#define WRAPPER_LIFE_SERIAL_H_INCLUDED

#define life_init(a,b,c,d) FUNCNAME_LIFEINIT(C_2DTITPI(b,c,a,d))
#define life_update(a,b,c) FUNCNAME_LIFEUPDATE(C_2DTPI(a,b,c))

__MATHSUITE __JBURKARDT void * FUNCNAME_LIFEUPDATE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LIFEINIT(void *);

#endif

