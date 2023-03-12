#ifndef WRAPPER_SPHERE_LLQ_GRID_H_INCLUDED
#define WRAPPER_SPHERE_LLQ_GRID_H_INCLUDED

#define sphere_llq_grid_point_count(a,b) R_USHRT(FUNCNAME_SPHERELLQGRIDPOINTCOUNT(C_PUSHRT2(a,b)))
#define sphere_llq_grid_points(a,b,c,d,e) FUNCNAME_SPHERELLQGRIDPOINTS(C_ITPIT3DT(a,b,c,d,e))

__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERELLQGRIDPOINTS(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_SPHERELLQGRIDPOINTCOUNT(void *);

#endif

