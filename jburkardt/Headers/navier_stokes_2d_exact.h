#ifndef WRAPPER_NAVIER_STOKES_2D_EXACT_H_INCLUDED
#define WRAPPER_NAVIER_STOKES_2D_EXACT_H_INCLUDED

#define resid_lucas(a,b,c,d,e,f,g,h,i) FUNCNAME_RESIDLUCAS(C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i))
#define resid_poiseuille(a,b,c,d,e,f,g,h,i) FUNCNAME_RESIDPOISEUILLE(C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i))
#define resid_spiral(a,b,c,d,e,f,g,h,i) FUNCNAME_RESIDSPIRAL(C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i))
#define resid_taylor(a,b,c,d,e,f,g,h,i) FUNCNAME_RESIDTAYLOR(C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i))
#define resid_vortex(a,b,c,d,e,f,g,h,i) FUNCNAME_RESIDVORTEX(C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i))
#define rhs_lucas(a,b,c,d,e,f,g,h,i) FUNCNAME_RHSLUCAS(C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i))
#define rhs_poiseuille(a,b,c,d,e,f,g,h,i) FUNCNAME_RHSPOISEUILLE(C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i))
#define rhs_spiral(a,b,c,d,e,f,g,h,i) FUNCNAME_RHSSPIRAL(C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i))
#define rhs_taylor(a,b,c,d,e,f,g,h,i) FUNCNAME_RHSTAYLOR(C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i))
#define rhs_vortex(a,b,c,d,e,f,g,h,i) FUNCNAME_RHSVORTEX(C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i))
#define resid_burgers(a,b,c,d,e,f,g,h,i,j) FUNCNAME_RESIDBURGERS(C_ITDT8PIT(a,b,c,d,e,f,g,h,i,j))
#define resid_ethier(a,b,c,d,e,f,g,h,i,j,k) FUNCNAME_RESIDETHIER(C_2ITDT8PIT(a,b,c,d,e,f,g,h,i,j,k))

__MATHSUITE __JBURKARDT void * FUNCNAME_RESIDLUCAS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RESIDPOISEUILLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RESIDSPIRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RESIDTAYLOR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RESIDVORTEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RHSLUCAS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RHSPOISEUILLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RHSSPIRAL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RHSTAYLOR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RHSVORTEX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RESIDBURGERS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RESIDETHIER(void *);

#endif
