#ifndef WRAPPER_MD_H_INCLUDED
#define WRAPPER_MD_H_INCLUDED

#define md_compute(a,b,c,d,e,f,g,h) FUNCNAME_MDCOMPUTE(C_2DT2PITIT3PIT(a,b,c,d,e,f,g,h))
#define md_dist(a,b,c,d) R_DBL(FUNCNAME_MDDIST(C_DT3PIT(a,b,c,d)))
#define md_initialize(a,b,c,d,e,f,g) FUNCNAME_MDINITIALIZE(C_2DTPITPI3PIT(a,b,c,d,e,f,g))
#define md_update(a,b,c,d,e,f,g,h) FUNCNAME_MDUPDATE(C_2DT4PIT2IT(a,b,c,d,e,f,g,h))

__MATHSUITE __JBURKARDT void * FUNCNAME_MDCOMPUTE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MDDIST(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MDINITIALIZE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_MDUPDATE(void *);

#endif
