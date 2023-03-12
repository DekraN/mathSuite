#ifndef WRAPPER_PCE_ODE_HERMITE_H_INCLUDED
#define WRAPPER_PCE_ODE_HERMITE_H_INCLUDED

#define pce_ode_hermite(a,b,c,d,e,f,g,h,i) FUNCNAME_PCEODEHERMITE(C_2ITDTITDT2IT2PIT(a,b,c,d,e,f,g,h,i))

__MATHSUITE __JBURKARDT void * FUNCNAME_PCEODEHERMITE(void *);

#endif
