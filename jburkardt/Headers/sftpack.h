#ifndef WRAPPER_SFTPACK_H_INCLUDED
#define WRAPPER_SFTPACK_H_INCLUDED

#define c8vec_sftb(a,b) FUNCNAME_C8VECSFTB(C_DTPCX(a,b))
#define c8vec_sftf(a,b) FUNCNAME_C8VECSFTF(C_DTPCX(a,b))
#define r8vec_sct(a,b) FUNCNAME_R8VECSCT(C_DTPIT(a,b))
#define r8vec_sftb(a,b,c,d) FUNCNAME_R8VECSFTB(C_DTIT2PIT(a,b,c,d))
#define r8vec_sftf(a,b,c,d,e) FUNCNAME_R8VECSFTF(C_DT4PIT(a,b,c,d,e))
#define r8vec_sht(a,b) FUNCNAME_R8VECSHT(C_DTPIT(a,b))
#define r8vec_sqctb(a,b) FUNCNAME_R8VECSQCTB(C_DTPIT(a,b))
#define r8vec_sqctf(a,b) FUNCNAME_R8VECSQCTF(C_DTPIT(a,b))
#define r8vec_sqstb(a,b) FUNCNAME_R8VECSQSTB(C_DTPIT(a,b))
#define r8vec_sqstf(a,b) FUNCNAME_R8VECSQSTF(C_DTPIT(a,b))
#define r8vec_sst(a,b) FUNCNAME_R8VECSST(C_DTPIT(a,b))

__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSFTF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSCT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSFTB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSHT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSQCTB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSQCTF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSQSTB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSQSTF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_R8VECSST(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8VECSFTB(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_C8VECSFTF(void *);


#endif
