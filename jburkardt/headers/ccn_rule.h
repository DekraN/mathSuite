#ifndef WRAPPER_CCN_RULE_H_INCLUDED
#define WRAPPER_CCN_RULE_H_INCLUDED

#define ccn_compute_points_new(a,b) FUNCNAME_CCNCOMPUTEPOINTSNEW(C_DTPIT(a,b))
#define rescale(a,b,c,d,e) FUNCNAME_RESCALE(C_DT2IT2PIT(c,a,b,d,e))

__MATHSUITE __JBURKARDT void * FUNCNAME_RESCALE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CCNCOMPUTEPOINTSNEW(void *);

#endif
