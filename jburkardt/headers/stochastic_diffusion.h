#ifndef WRAPPER_STOCHASTIC_DIFFUSION_H_INCLUDED
#define WRAPPER_STOCHASTIC_DIFFUSION_H_INCLUDED

#define diffusivity_1d_xk(a,b,c,d,e) FUNCNAME_DIFFUSIVITY1DXK(C_IT2DT2PIT(a,b,d,c,e))
#define diffusivity_2d_bnt(a,b,c,d,e) FUNCNAME_DIFFUSIVITY2DBNT(C_DTIT3PIT(c,a,d,e,b))
#define diffusivity_2d_elman(a,b,c,d,e,f,g,h,i) FUNCNAME_DIFFUSIVITY2DELMAN(C_3ITDTPIT2DT2PIT(a,b,c,d,e,f,g,h,i))
#define diffusivity_2d_ntw(a,b,c,d,e,f,g,h) FUNCNAME_DIFFUSIVITY2NTW(C_2ITDTPITDT2PIT(a,b,c,d,e,f,g,h))
#define r8vec_mesh_2d(a,b,c,d,e,f) FUNCNAME_R8VECMESH2D(C_2DT4PIT(a,b,c,d,e,f))
#define theta_solve(a,b,c) FUNCNAME_THETASOLVE(C_2ITDT(a,b,c))

__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFUSIVITY1DXK(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFUSIVITY2DBNT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFUSIVITY2DELMAN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DIFFUSIVITY2NTW(void *);	
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECMESH2D(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_THETASOLVE(void *);

#endif
