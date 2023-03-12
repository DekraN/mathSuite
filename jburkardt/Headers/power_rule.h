#ifndef WRAPPER_POWER_RULE_H_INCLUDED
#define WRAPPER_POWER_RULE_H_INCLUDED

#define power_rule_set(a,b,c,d,e,f,g,h,i) FUNCNAME_POWERRULESET(C_3DT6PIT(a,e,f,b,c,d,g,h,i))
#define power_rule_size(a,b) R_INT(FUNCNAME_POWERRULESIZE(C_PUSHRT2(a,b)))

__MATHSUITE __JBURKARDT void * FUNCNAME_POWERRULESET(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POWERRULESIZE(void *);

#endif
