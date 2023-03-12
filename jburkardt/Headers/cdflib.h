#ifndef WRAPPER_CDFLIB_H_INCLUDED
#define WRAPPER_CDFLIB_H_INCLUDED

#define algdiv(a,b) R_DBL(FUNCNAME_ALGDIV(C_PPDBL2(a,b)))
#define alnrel(a) R_DBL(FUNCNAME_ALNREL(a))
#define apser(a,b,c,d) R_DBL(FUNCNAME_APSER(C_PPDBL4(a,b,c,d)))
#define bcorr(a,b) R_DBL(FUNCNAME_BCORR(C_PPDBL2(a,b)))
#define beta(a,b) R_DBL(FUNCNAME_BETA(C_PDBL2(a,b)))
#define beta_asym(a,b,c,d) R_DBL(FUNCNAME_BETAASYM(C_PPDBL4(a,b,c,d)))
#define beta_frac(a,b,c,d,e,f) R_DBL(FUNCNAME_BETAFRAC(C_PPDBL6(a,b,c,d,e,f)))
#define beta_grat(a,b,c,d,e,f,g) FUNCNAME_BETAGRAT(C_PDT6PIT(g,a,b,c,d,e,f))
#define beta_inc(a,b,c,d,e,f,g) FUNCNAME_BETAINC(C_PDT6PIT(g,a,b,c,d,e,f))
#define beta_log(a,b) R_DBL(FUNCNAME_BETALOG(C_PPDBL2(a,b)))
#define beta_pser(a,b,c,d) R_DBL(FUNCNAME_BETAPSER(C_PPDBL4(a,b,c,d)))
#define beta_rcomp(a,b,c,d) R_DBL(FUNCNAME_BETARCOMP(C_PPDBL4(a,b,c,d)))
#define beta_rcomp1(a,b,c,d,e) R_DBL(FUNCNAME_BETARCOMP1(C_PDT4PIT(a,b,c,d,e)))
#define beta_up(a,b,c,d,e,f) R_DBL(FUNCNAME_BETAUP(C_4PITPDTPIT(a,b,c,d,e,f)))
#define cdfbet(a,b,c,d,e,f,g,h,i) FUNCNAME_CDFBET(C_PDT6PITPSPIT(a,b,c,d,e,f,g,h,i))
#define cdfbin(a,b,c,d,e,f,g,h,i) FUNCNAME_CDFBIN(C_PDT6PITPSPIT(a,b,c,d,e,f,g,h,i))
#define cdfchi(a,b,c,d,e,f,g) FUNCNAME_CDFCHI(C_PDT4PITPSPIT(a,b,c,d,e,f,g))
#define cdfchn(a,b,c,d,e,f,g,h) FUNCNAME_CDFCHN(C_PDT5PITPSPIT(a,b,c,d,e,f,g,h))
#define cdff(a,b,c,d,e,f,g,h) FUNCNAME_CDFF(C_PDT5PITPSPIT(a,b,c,d,e,f,g,h))
#define cdffnc(a,b,c,d,e,f,g,h,i) FUNCNAME_CDFFNC(C_PDT6PITPSPIT(a,b,c,d,e,f,g,h,i))
#define cdfgam(a,b,c,d,e,f,g,h) FUNCNAME_CDFGAM(C_PDT5PITPSPIT(a,b,c,d,e,f,g,h))
#define cdfnbn(a,b,c,d,e,f,g,h,i) FUNCNAME_CDFNBN(C_PDT6PITPSPIT(a,b,c,d,e,f,g,h,i))
#define cdfnor(a,b,c,d,e,f,g,h) FUNCNAME_CDFNOR(C_PDT5PITPSPIT(a,b,c,d,e,f,g,h))
#define cdfpoi(a,b,c,d,e,f,g) FUNCNAME_CDFPOI(C_PDT4PITPSPIT(a,b,c,d,e,f,g))
#define cdft(a,b,c,d,e,f,g) FUNCNAME_CDFT(C_PDT4PITPSPIT(a,b,c,d,e,f,g))
#define cumbet(a,b,c,d,e,f) FUNCNAME_CUMBET(C_PPDBL6(a,b,c,d,e,f))
#define cumbin(a,b,c,d,e,f) FUNCNAME_CUMBIN(C_PPDBL6(a,b,c,d,e,f))
#define cumchi(a,b,c,d) FUNCNAME_CUMCHI(C_PPDBL4(a,b,c,d))
#define cumchn(a,b,c,d,e) FUNCNAME_CUMCHN(C_PPDBL5(a,b,c,d,e))
#define cumf(a,b,c,d,e) FUNCNAME_CUMF(C_PPDBL5(a,b,c,d,e))
#define cumfnc(a,b,c,d,e,f) FUNCNAME_CUMFNC(C_PPDBL6(a,b,c,d,e,f))
#define cumgam(a,b,c,d) FUNCNAME_CUMGAM(C_PPDBL4(a,b,c,d))
#define cumnbn(a,b,c,d,e,f) FUNCNAME_CUMNBN(C_PPDBL6(a,b,c,d,e,f))
#define cumnor(a,b,c) FUNCNAME_CUMNOR(C_PPDBL3(a,b,c))
#define cumpoi(a,b,c,d) FUNCNAME_CUMPOI(C_PPDBL4(a,b,c,d))
#define cumt(a,b,c,d) FUNCNAME_CUMT(C_PPDBL4(a,b,c,d))
#define dbetrm(a,b) R_DBL(FUNCNAME_DBETRM(C_PPDBL2(a,b)))
#define dexpm1(a) R_DBL(FUNCNAME_DEXPM1(a))
#define dinvnr(a,b) R_DBL(FUNCNAME_DINVNR(C_PPDBL2(a,b)))
#define dinvr(a,b,c,d,e) FUNCNAME_DINVR(C_PS2PIT2PUL(a,b,c,d,e))
#define dlanor(a) R_DBL(FUNCNAME_DLANOR(a))
#define dpmpar(a) R_DBL(FUNCNAME_DPMPAR(a))
#define dstinv(a,b,c,d,e,f,g) FUNCNAME_DSTINV(C_PPDBL7(a,b,c,d,e,f,g))
#define dstrem(a) R_DBL(FUNCNAME_DSTREM(a))
#define dstzr(a,b,c,d) FUNCNAME_DSTZR(C_PPDBL4(a,b,c,d))
#define dt1(a,b,c) R_DBL(FUNCNAME_DT1(C_PPDBL3(a,b,c)))
#define dzror(a,b,c,d,e,f,g) FUNCNAME_DZROR(C_PS4PIT2PUL(a,b,c,d,e,f,g))
#define E0000(a,b,c,d,e,f,g,h,i,j,k,l,m) FUNCNAME_E0000(C_IPS2PIT2PUL7PIT(a,b,c,d,e,f,g,h,i,j,k,l,m))
#define E0001(a,b,c,d,e,f,g,h,i,j,k,l) FUNCNAME_E0001(C_IPS4PIT2PUL4PIT(a,b,c,d,e,f,g,h,i,j,k,l))
#define error_f(a) R_DBL(FUNCNAME_ERRORF(a))
#define error_fc(a,b) R_DBL(FUNCNAME_ERRORFC(C_PDTPIT(a,b)))
#define esum(a,b) R_DBL(FUNCNAME_ESUM(C_PDTPIT(a,b)))
#define eval_pol(a,b,c) R_DBL(FUNCNAME_EVALPOL(C_PDT2PIT(b,a,c)))
#define exparg(a) R_DBL(FUNCNAME_EXPARG(a))
#define fifdint(a) R_DBL(FUNCNAME_FIFDINT(C_SUSHRT(a)))
#define fifdmax1(a,b) R_DBL(FUNCNAME_FIFDMAX1(C_PDBL2(a,b)))
#define fifdmin1(a,b) R_DBL(FUNCNAME_FIFDMIN1(C_PDBL2(a,b)))
#define fifdsign(a,b) R_DBL(FUNCNAME_FIFDSIGN(C_PDBL2(a,b)))
#define fifidint(a) R_LNG(FUNCNAME_FIFIDINT(C_SDBL(a)))
#define fifmod(a,b) R_LNG(FUNCNAME_FIFMOD(C_PLNG2(a,b)))
#define fpser(a,b,c,d) R_DBL(FUNCNAME_FPSER(C_PPDBL4(a,b,c,d)))
#define gam1(a) R_DBL(FUNCNAME_GAM1(a))
#define gamma_inc_inv(a,b,c,d,e,f) FUNCNAME_GAMMAINCINV(C_4PITPDTPIT(a,b,c,d,f,e))
#define gamma_ln1(a) R_DBL(FUNCNAME_GAMMALN1(a))
#define gamma_log(a) R_DBL(FUNCNAME_GAMMALOG(a))
#define gamma_rat1(a,b,c,d,e,f) FUNCNAME_GAMMARAT1(C_PPDBL6(a,b,c,d,e,f))
#define gamma_x(a) R_DBL(FUNCNAME_GAMMAX(a))
#define gsumln(a,b) R_DBL(FUNCNAME_GSUMLN(C_PPDBL2(a,b)))
#define ipmpar(a) R_INT(FUNCNAME_IPMPAR(a))
#define psi(a) R_DBL(FUNCNAME_PSI(a))
#define rcomp(a,b) R_DBL(FUNCNAME_RCOMP(C_PPDBL2(a,b)))
#define rexp(a) R_DBL(FUNCNAME_REXP(a))
#define rlog(a) R_DBL(FUNCNAME_RLOG(a))
#define rlog1(a) R_DBL(FUNCNAME_RLOG1(a))
#define stvaln(a) R_DBL(FUNCNAME_STVALN(a))
#define binomial_cdf_values(a,b,c,d,e) FUNCNAME_BINOMIALCDFVALUES(C_3PDT2PIT(a,b,d,c,e))
#define chi_noncentral_cdf_values(a,b,c,d,e) FUNCNAME_CHINONCENTRALCDFVALUES(C_PDT2PITPDTPIT(a,b,c,d,e))
#define chi_square_cdf_values(a,b,c,d) FUNCNAME_CHISQUARECDFVALUES(C_2PDT2PIT(a,b,c,d))
#define erf_values(a,b,c) FUNCNAME_ERFCDFVALUES(C_PDT2PIT(a,b,c))
#define f_cdf_values(a,b,c,d,e) FUNCNAME_FCDFVALUES(C_3PDT2PIT(a,b,c,d,e))
#define f_noncentral_cdf_values(a,b,c,d,e,f) FUNCNAME_FNONCENTRALCDFVALUES(C_3PDT3PIT(a,b,c,d,e,f))
#define negative_binomial_cdf_values(a,b,c,d,e) FUNCNAME_NEGATIVEBINOMIALCDFVALUES(C_3PDT2PIT(a,b,c,d,e))
#define normal_cdf_values(a,b,c) FUNCNAME_NORMALCDFVALUES(C_PDT2PIT(a,b,c))
#define poisson_cdf_values(a,b,c,d) FUNCNAME_POISSONCDFVALUES(C_PDTPITPDTPIT(a,c,b,d))
#define student_cdf_values(a,b,c,d) FUNCNAME_STUDENTCDFVALUES(C_2PDT2PIT(a,b,c,d))

__MATHSUITE __JBURKARDT void * FUNCNAME_ALNREL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DLANOR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DPMPAR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSTREM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ERRORF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EXPARG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FIFDINT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FIFIDINT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAM1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAMMAINCINV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAMMALN1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAMMALOG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAMMARAT1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GAMMAX(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_IPMPAR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_PSI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_REXP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RLOG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RLOG1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ALGDIV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BCORR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETA(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETAASYM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETALOG(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DBETRM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DINVNR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DINVR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ERRORFC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ESUM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FIFDMAX1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FIFDMIN1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FIFDSIGN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FIFMOD(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_GSUMLN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_RCOMP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMNOR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DT1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_APSER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETAPSER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETARCOMP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMCHI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMGAM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMPOI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMBIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSTZR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FPSER(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETARCOMP1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMCHN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMFNC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMBET(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CUMNBN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETAFRAC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETAGRAT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETAINC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BETAUP(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFBET(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFBIN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFCHI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFCHN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFF(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFFNC(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFGAM(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFNBN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFNOR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFPOI(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CDFT(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DSTINV(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DEXPM1(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_DZROR(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_E0000(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_E0001(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_EVALPOL(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_STVALN(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_BINOMIALCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHINONCENTRALCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_CHISQUARECDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_ERFCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_FNONCENTRALCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NEGATIVEBINOMIALCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_NORMALCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_POISSONCDFVALUES(void *);
__MATHSUITE __JBURKARDT void * FUNCNAME_STUDENTCDFVALUES(void *);


#endif
