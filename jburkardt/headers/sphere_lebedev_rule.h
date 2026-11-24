#ifndef WRAPPER_SPHERE_LEBEDEV_RULE_H_INCLUDED
#define WRAPPER_SPHERE_LEBEDEV_RULE_H_INCLUDED

#define available_table(a) R_SHRT(FUNCNAME_AVAILABLETABLE(C_SUSHRT(a)))
#define gen_oh(a,b,c,d,e,f,g,h) R_USHRT(FUNCNAME_GENOH(C_DT3IT4PIT(a,b,c,d,e,f,g,h)))
#define ld_by_order(a,b,c,d,e) FUNCNAME_LDBYORDER(C_DT4PIT(a,b,c,d,e))
#define ld0006(a,b,c,d) FUNCNAME_LD0006(C_PPDBL4(a,b,c,d))
#define ld0014(a,b,c,d) FUNCNAME_LD0014(C_PPDBL4(a,b,c,d))
#define ld0026(a,b,c,d) FUNCNAME_LD0026(C_PPDBL4(a,b,c,d))
#define ld0038(a,b,c,d) FUNCNAME_LD0038(C_PPDBL4(a,b,c,d))
#define ld0050(a,b,c,d) FUNCNAME_LD0050(C_PPDBL4(a,b,c,d))
#define ld0074(a,b,c,d) FUNCNAME_LD0074(C_PPDBL4(a,b,c,d))
#define ld0086(a,b,c,d) FUNCNAME_LD0086(C_PPDBL4(a,b,c,d))
#define ld0110(a,b,c,d) FUNCNAME_LD0110(C_PPDBL4(a,b,c,d))
#define ld0146(a,b,c,d) FUNCNAME_LD0146(C_PPDBL4(a,b,c,d))
#define ld0170(a,b,c,d) FUNCNAME_LD0170(C_PPDBL4(a,b,c,d))
#define ld0194(a,b,c,d) FUNCNAME_LD0194(C_PPDBL4(a,b,c,d))
#define ld0230(a,b,c,d) FUNCNAME_LD0230(C_PPDBL4(a,b,c,d))
#define ld0266(a,b,c,d) FUNCNAME_LD0266(C_PPDBL4(a,b,c,d))
#define ld0302(a,b,c,d) FUNCNAME_LD0302(C_PPDBL4(a,b,c,d))
#define ld0350(a,b,c,d) FUNCNAME_LD0350(C_PPDBL4(a,b,c,d))
#define ld0434(a,b,c,d) FUNCNAME_LD0434(C_PPDBL4(a,b,c,d))
#define ld0590(a,b,c,d) FUNCNAME_LD0590(C_PPDBL4(a,b,c,d))
#define ld0770(a,b,c,d) FUNCNAME_LD0770(C_PPDBL4(a,b,c,d))
#define ld0974(a,b,c,d) FUNCNAME_LD0974(C_PPDBL4(a,b,c,d))
#define ld1202(a,b,c,d) FUNCNAME_LD1202(C_PPDBL4(a,b,c,d))
#define ld1454(a,b,c,d) FUNCNAME_LD1454(C_PPDBL4(a,b,c,d))
#define ld1730(a,b,c,d) FUNCNAME_LD1730(C_PPDBL4(a,b,c,d))
#define ld2030(a,b,c,d) FUNCNAME_LD2030(C_PPDBL4(a,b,c,d))
#define ld2354(a,b,c,d) FUNCNAME_LD2354(C_PPDBL4(a,b,c,d))
#define ld2702(a,b,c,d) FUNCNAME_LD2702(C_PPDBL4(a,b,c,d))
#define ld3074(a,b,c,d) FUNCNAME_LD3074(C_PPDBL4(a,b,c,d))
#define ld3470(a,b,c,d) FUNCNAME_LD3470(C_PPDBL4(a,b,c,d))
#define ld3890(a,b,c,d) FUNCNAME_LD3890(C_PPDBL4(a,b,c,d))
#define ld4334(a,b,c,d) FUNCNAME_LD4334(C_PPDBL4(a,b,c,d))
#define ld4802(a,b,c,d) FUNCNAME_LD4802(C_PPDBL4(a,b,c,d))
#define ld5294(a,b,c,d) FUNCNAME_LD5294(C_PPDBL4(a,b,c,d))
#define ld5810(a,b,c,d) FUNCNAME_LD5810(C_PPDBL4(a,b,c,d))
#define order_table(a) R_USHRT(FUNCNAME_ORDERTABLE(C_SUSHRT(a)))
#define precision_table(a) R_USHRT(FUNCNAME_PRECISIONTABLE(C_SUSHRT(a)))

__MATHSUITE __JBURKARDT void * FUNCNAME_GENOH(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LDBYORDER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0006(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0014(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0026(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0038(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0050(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0074(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0086(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0110(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0146(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0170(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0194(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0230(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0266(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0302(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0350(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0434(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0590(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0770(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD0974(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD1202(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD1454(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD1730(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD2030(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD2354(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD2702(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD3074(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD3470(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD3890(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD4334(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD4802(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD5294(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_LD5810(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_AVAILABLETABLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ORDERTABLE(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PRECISIONTABLE(void *);

#endif
