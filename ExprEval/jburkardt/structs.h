// Structure and Macro Definitions
// for EXPREVAL C MapperParser

#define R_INT(a) (*(int *)(a))
#define R_UINT(a) (*(unsigned int *)(a))
#define R_SHRT(a) (*(short *)(a))
#define R_USHRT(a) (*(dim_typ *)(a))
#define R_CHR(a) (*(char *)(a))
#define R_UCHR(a) (*(unsigned char *)(a))
#define R_LNG(a) (*(long *)(a))
#define R_ULNG(a) (*(unsigned long *)(a))
#define R_LLNG(a) (*(long long *)(a))
#define R_ULLNG(a) (*(unsigned long long *)(a))
#define R_FLT(a) (*(float *)(a))
#define R_DBL(a) (*(double *)(a))
#define R_LDBL(a) (*(long double *)(a))

#define C_SINT(a) (&((sint){(a)}))
#define C_SUINT(a) (&((suint){(a)}))
#define C_SSHRT(a) (&((sshrt){(a)}))
#define C_SUSHRT(a) (&((sushrt){(a)}))
#define C_SCHR(a) (&((schr){(a)}))
#define C_SUCHR(a) (&((suchr){(a)}))
#define C_SLNG(a) (&((slng){(a)}))
#define C_SULNG(a) (&((sulng){(a)}))
#define C_SLLNG(a) (&((sllng){(a)}))
#define C_SULLNG(a) (&((sullng){(a)}))
#define C_SFLT(a) (&((sflt){(a)}))
#define C_SDBL(a) (&((sdbl){(a)}))
#define C_SLDBL(a) (&((sldbl){(a)}))
#define C_2DT2PIT(a,b,c,d) (&((_2dt2pit){(a),(b),(c),(d)}))
#define C_ITB(a,b) (&((itb){(a),(b)}))
#define C_PIT2IPITPDT(a,b,c,d,e) (&((pit2ipitpdt){(a),(b),(c),(d),(e)}))
#define C_PITDTPDTPITPDTDTPIT(a,b,c,d,e,f,g) (&((pitdtpdtpitpdtdtpit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2PDT2PIT(a,b,c,d) (&((_2pdt2pit){(a),(b),(c),(d)}))
#define C_PDT5PITPDT2DT(a,b,c,d,e,f,g,h,i) (&((pdt5pitpdt2dt){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_PDT4PIT(a,b,c,d,e) (&((pdt4pit){(a),(b),(c),(d),(e)}))
#define C_PDT3PIT(a,b,c,d) (&((pdt3pit){(a),(b),(c),(d)}))
#define C_PDT2PIT(a,b,c) (&((pdt2pit){(a),(b),(c)}))
#define C_4PDTPITPB(a,b,c,d,e,f) (&((_4pdtpitpb){(a),(b),(c),(d),(e),(f)}))
#define C_B3DT(a,b,c,d,e) (&((b3dt){(a),(b),(c),(d),(e)}))
#define C_5PDTPIT(a,b,c,d,e,f) (&((_5pdtpit){(a),(b),(c),(d),(e),(f)}))
#define C_DTPIT(a,b) (&((dtpit){(a),(b)}))
#define C_2BDTPITPDTPIT(a,b,c,d,e,f) (&((_2bdtpitpdtpit){(a),(b),(c),(d),(e),(f)}))
#define C_2DTPIT(a,b,c) (&((_2dtpit){(a),(b),(c)}))
#define C_DT2PIT2PDT(a,b,c,d,e) (&((dt2pit2pdt){(a),(b),(c),(d),(e)}))
#define C_DT2PIT2PDT2PIT(a,b,c,d,e,f,g) (&((dt2pit2pdt2pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DT2I2PI(a,b,c,d,e) (&((dt2i2pi){(a),(b),(c),(d),(e)}))
#define C_DT3ITPIT(a,b,c,d,e) (&((dt3itpit){(a),(b),(c),(d),(e)}))
#define C_3PDT2PIT(a,b,c,d,e) (&((_3pdt2pit){(a),(b),(c),(d),(e)}))
#define C_DTIT(a,b) (&((dtit){(a),(b)}))
#define C_DTITPITDTPITDT(a,b,c,d,e,f) (&((dtitpitdtpitdt){(a),(b),(c),(d),(e),(f)}))
#define C_DTPITI(a,b,c) (&((dtpiti){(a),(b),(c)}))
#define C_DTPITIPITI2IT(a,b,c,d,e,f,g) (&((dtpitipiti2it){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DTITPITI(a,b,c,d) (&((dtitpiti){(a),(b),(c),(d)}))
#define C_DTPITIPITI(a,b,c,d,e) (&((dtpitipiti){(a),(b),(c),(d),(e)}))
#define C_7ITFITPIT(a,b,c,d,e,f,g,h,i) (&((_7itfitpit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_4PITFITPIT(a,b,c,d,e,f) (&((_4pitfitpit){(a),(b),(c),(d),(e),(f)}))
#define C_4PITFIT(a,b,c,d,e) (&((_4pitfit){(a),(b),(c),(d),(e)}))
#define C_2DT2ITPITPI(a,b,c,d,e,f) (&((_2dt2itpitpi){(a),(b),(c),(d),(e),(f)}))
#define C_DT2PIT(a,b,c) (&((dt2pit){(a),(b),(c)}))
#define C_DT3PIT(a,b,c,d) (&((dt3pit){(a),(b),(c),(d)}))
#define C_DTIT3PIT(a,b,c,d,e) (&((dtit3pit){(a),(b),(c),(d),(e)}))
#define C_2ITDT3PIT(a,b,c,d,e,f) (&((_2itdt3pit){(a),(b),(c),(d),(e),(f)}))
#define C_2DT2PIDT3PI(a,b,c,d,e,f,g,h) (&((_2dt2pidt3pi){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_4PITPDTPIT(a,b,c,d,e,f) (&((_4pitpdtpit){(a),(b),(c),(d),(e),(f)}))
#define C_PDT6PITPSPIT(a,b,c,d,e,f,g,h,i) (&((pdt6pitpspit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_PDT4PITPSPIT(a,b,c,d,e,f,g) (&((pdt4pitpspit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_PDT5PITPSPIT(a,b,c,d,e,f,g,h) (&((pdt5pitpspit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_PS2PIT2PUL(a,b,c,d,e) (&((ps2pit2pul){(a),(b),(c),(d),(e)}))
#define C_IPS2PIT2PUL7PIT(a,b,c,d,e,f,g,h,i,j,k,l,m) (&((ips2pit2pul7pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m)}))
#define C_IPS4PIT2PUL4PIT(a,b,c,d,e,f,g,h,i,j,k,l) (&((ips4pit2pul4pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l)}))
#define C_PS4PIT2PUL(a,b,c,d,e,f,g) (&((ps4pit2pul){(a),(b),(c),(d),(e),(f),(g)}))
#define C_PDTPIT(a,b) (&((pdtpit){(a),(b)}))
#define C_PDT2PITPDTPIT(a,b,c,d,e) (&((pdt2pitpdtpit){(a),(b),(c),(d),(e)}))
#define C_3PDT3PIT(a,b,c,d,e,f) (&((_3pdt3pit){(a),(b),(c),(d),(e),(f)}))
#define C_2IPI(a,b,c) (&((_2ipi){(a),(b),(c)}))
#define C_DTPI(a,b) (&((dtpi){(a),(b)}))
#define C_DTFITPIT(a,b,c) (&((dtfitpit){(a),(b),(c)}))
#define C_DTFIT2ITPIT(a,b,c,d,e) (&((dtfit2itpit){(a),(b),(c),(d),(e)}))
#define C_DTPITDTPIPIT(a,b,c,d,e) (&((dtpitdtpipit){(a),(b),(c),(d),(e)}))
#define C_ITPITDT4PIT(a,b,c,d,e,f,g) (&((itpitdt4pit ){(a),(b),(c),(d),(e),(f),(g)}))
#define C_IT3PIT(a,b,c,d) (&((it3pit){(a),(b),(c),(d)}))
#define C_IT3PITDTPI(a,b,c,d,e,f) (&((it3pitdtpi){(a),(b),(c),(d),(e),(f)}))
#define C_ITPIT2ITPIT(a,b,c,d,e) (&((itpit2itpit){(a),(b),(c),(d),(e)}))
#define C_IT3PITDTPI2PIT(a,b,c,d,e,f,g,h) (&((it3pitdtpi2pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_2ITDTPI2PIT(a,b,c,d,e,f) (&((_2itdtpi2pit){(a),(b),(c),(d),(e),(f)}))
#define C_DT2ITPIPIT(a,b,c,d,e) (&((dt2itpipit){(a),(b),(c),(d),(e)}))
#define C_2DTPI(a,b,c) (&((_2dtpi){(a),(b),(c)}))
#define C_2DTPIPB(a,b,c,d) (&((_2dtpipb){(a),(b),(c),(d)}))
#define C_DTPII(a,b,c) (&((dtpii){(a),(b),(c)}))
#define C_2DT6PI(a,b,c,d,e,f,g,h) (&((_2dt6pi){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_2ITDT(a,b,c) (&((_2itdt){(a),(b),(c)}))
#define C_DTPITIT(a,b,c) (&((dtpitit){(a),(b),(c)}))
#define C_IT2PIT(a,b,c) (&((it2pit){(a),(b),(c)}))
#define C_2IT2PITPDT(a,b,c,d,e) (&((_2it2pitpdt){(a),(b),(c),(d),(e)}))
#define C_DT2ITPIT(a,b,c,d) (&((dt2itpit){(a),(b),(c),(d)}))
#define C_FIT2ITDT(a,b,c,d) (&((fit2itdt){(a),(b),(c),(d)}))
#define C_IPI(a,b) (&((ipi){(a),(b)}))
#define C_FDTDT2PDT(a,b,c,d) (&((fdtdt2pdt){(a),(b),(c),(d)}))
#define C_IT2DT2PIT(a,b,c,d,e) (&((it2dt2pit){(a),(b),(c),(d),(e)}))
#define C_DT2PITDTPITPS(a,b,c,d,e,f) (&((dt2pitdtpitps){(a),(b),(c),(d),(e),(f)}))
#define C_DT2PITDTPS(a,b,c,d,e) (&((dt2pitdtps){(a),(b),(c),(d),(e)}))
#define C_3DT2PIT(a,b,c,d,e) (&((_3dt2pit){(a),(b),(c),(d),(e)}))
#define C_DT2PITPDT2PIT(a,b,c,d,e,f) (&((dt2pitpdt2pit){(a),(b),(c),(d),(e),(f)}))
#define C_DT2PIT2ITPDT2PIT(a,b,c,d,e,f,g,h) (&((dt2pit2itpdt2pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DT4PIT(a,b,c,d,e) (&((dt4pit){(a),(b),(c),(d),(e)}))
#define C_DT2PITI2PIT(a,b,c,d,e,f) (&((dt2piti2pit){(a),(b),(c),(d),(e),(f)}))
#define C_DT2IT2PIT(a,b,c,d,e) (&((dt2it2pit){(a),(b),(c),(d),(e)}))
#define C_DTPITPDT2PIT(a,b,c,d,e) (&((dtpitpdt2pit){(a),(b),(c),(d),(e)}))
#define C_DTPITITPI(a,b,c,d) (&((dtpititpi){(a),(b),(c),(d)}))
#define C_2DT2PITITPI(a,b,c,d,e,f) (&((_2dt2pititpi){(a),(b),(c),(d),(e),(f)}))
#define C_DT4FITPIT(a,b,c,d,e,f) (&((dt4fitpit){(a),(b),(c),(d),(e),(f)}))
#define C_DT3FITPIT(a,b,c,d,e) (&((dt3fitpit){(a),(b),(c),(d),(e)}))
#define C_DT2PITFIT(a,b,c,d) (&((dt2pitfit){(a),(b),(c),(d)}))
#define C_DT4IT2FITPIT(a,b,c,d,e,f,g,h) (&((dt4it2fitpit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DTPITDTFIT3PIT(a,b,c,d,e,f,g) (&((dtpitdtfit3pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2DTPI3PDT(a,b,c,d,e,f) (&((_2dtpi3pdt){(a),(b),(c),(d),(e),(f)}))
#define C_2ITDT2PITDT2PITDT(a,b,c,d,e,f,g,h,i) (&((_2itdt2pitdt2pitdt){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_DTPDT2PITDTPITDTPIT(a,b,c,d,e,f,g,h) (&((dtpdt2pitdtpitdtpit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_6IT2PIT(a,b,c,d,e,f,g,h) (&((_6it2pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_2DTIT(a,b,c) (&((_2dtit){(a),(b),(c)}))
#define C_3DT2IT(a,b,c,d,e) (&((_3dt2it){(a),(b),(c),(d),(e)}))
#define C_4DT3IT(a,b,c,d,e,f,g) (&((_4dt3it){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2PDTPIT(a,b,c) (&((_2pdtpit){(a),(b),(c)}))
#define C_PIT3DTIT2PIPIT(a,b,c,d,e,f,g,h) (&((pit3dtit2pipit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_PIT3DTITPI3PITPIPITI(a,b,c,d,e,f,g,h,i,j,k,l) (&((pit3dtitpi3pitpipiti){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l)}))
#define C_PIT3DTIT2PI3PITPIPITI(a,b,c,d,e,f,g,h,i,j,k,l,m) (&((pit3dtit2pi3pitpipiti){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m)}))
#define C_PIT4DT3PITPIPIT(a,b,c,d,e,f,g,h,i,j) (&((pit4dt3pitpipit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)}))
#define C_IT5PIT(a,b,c,d,e,f) (&((it5pit){(a),(b),(c),(d),(e),(f)}))
#define C_ITPIT4ITPDTPIT(a,b,c,d,e,f,g,h) (&((itpit4itpdtpit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_IT4PIT(a,b,c,d,e) (&((it4pit){(a),(b),(c),(d),(e)}))
#define C_ITPIT2IT(a,b,c,d) (&((itpit2it){(a),(b),(c),(d)}))
#define C_ITPIT2ITPIT(a,b,c,d,e) (&((itpit2itpit){(a),(b),(c),(d),(e)}))
#define C_3IT3PIT(a,b,c,d,e,f) (&((_3it3pit){(a),(b),(c),(d),(e),(f)}))
#define C_CHITDT2PIT(a,b,c,d,e) (&((chitdt2pit){(a),(b),(c),(d),(e)}))
#define C_3DTPIT2PI(a,b,c,d,e,f) (&((_3dtpit2pi){(a),(b),(c),(d),(e),(f)}))
#define C_2PITITDTPI(a,b,c,d,e) (&((_2pititdtpi){(a),(b),(c),(d),(e)}))
#define C_DTPITPDT(a,b,c) (&((dtpitpdt){(a),(b),(c)}))
#define C_ITPITPI(a,b,c) (&((itpitpi){(a),(b),(c)}))
#define C_3DTPITPDTPI3DTPIT2PDT(a,b,c,d,e,f,g,h,i,j,k,l) (&((_3dtpitpdtpi3dtpit2pdt){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l)}))
#define C_3DTPIT6PDT(a,b,c,d,e,f,g,h,i,j,k) (&((_3dtpit6pdt){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)}))
#define C_2ITPIT(a,b,c) (&((_2itpit){(a),(b),(c)}))
#define C_PIT3ITDTPIT(a,b,c,d,e,f) (&((pit3itdtpit){(a),(b),(c),(d),(e),(f)}))
#define C_PIT5ITDTPIT(a,b,c,d,e,f,g,h) (&((pit5itdtpit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_6IT3PIT(a,b,c,d,e,f,g,h,i) (&((_6it3pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_4IT2PIT(a,b,c,d,e,f) (&((_4it2pit){(a),(b),(c),(d),(e),(f)}))
#define C_ITDT3ITPIT(a,b,c,d,e,f) (&((itdt3itpit){(a),(b),(c),(d),(e),(f)}))
#define C_4DTPIT3PI(a,b,c,d,e,f,g,h) (&((_4dtpit3pi){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_3PITPB(a,b,c,d) (&((_3pitpb){(a),(b),(c),(d)}))
#define C_3IT4PIT(a,b,c,d,e,f,g) (&((_3it4pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_4ITPIT(a,b,c,d,e) (&((_4itpit){(a),(b),(c),(d),(e)}))
#define C_6ITPIT(a,b,c,d,e,f,g) (&((_6itpit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_4IT3PIT(a,b,c,d,e,f,g) (&((_4it3pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_4PITPDTPIT(a,b,c,d,e,f) (&((_4pitpdtpit){(a),(b),(c),(d),(e),(f)}))
#define C_6ITPDTPIT(a,b,c,d,e,f,g,h) (&((_6itpdtpit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_8IT3PIT(a,b,c,d,e,f,g,h,i,j,k) (&((_8it3pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)}))
#define C_6IT5PIT(a,b,c,d,e,f,g,h,i,j,k) (&((_6it5pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)}))
#define C_3PITDT3PIT(a,b,c,d,e,f,g) (&((_3pitdt3pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_4PITDT2PITPDT(a,b,c,d,e,f,g,h) (&((_4pitdt2pitpdt){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_10ITPIT(a,b,c,d,e,f,g,h,i,j,k) (&((_10itpit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)}))
#define C_2PIT4IT3PIT(a,b,c,d,e,f,g,h,i) (&((_2pit4it3pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_4ITPITPIPIT(a,b,c,d,e,f,g) (&((_4itpitpipit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2PIT2ITPDTPIT(a,b,c,d,e,f) (&((_2pit2itpdtpit){(a),(b),(c),(d),(e),(f)}))
#define C_PI2PIT(a,b,c) (&((pi2pit){(a),(b),(c)}))
#define C_DTPITPDTPI(a,b,c,d) (&((dtpitpdtpi){(a),(b),(c),(d)}))
#define C_2DT3PIT(a,b,c,d,e) (&((_2dt3pit){(a),(b),(c),(d),(e)}))
#define C_DTPITPI(a,b,c) (&((dtpitpi){(a),(b),(c)}))
#define C_PIT2DTPDTDTPDT(a,b,c,d,e,f) (&((pit2dtpdtdtpdt){(a),(b),(c),(d),(e),(f)}))
#define C_3DTPIT2PDTPIT(a,b,c,d,e,f,g) (&((_3dtpit2pdtpit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2DTPITIT(a,b,c,d) (&((_2dtpitit){(a),(b),(c),(d)}))
#define C_2DT4PIT(a,b,c,d,e,f) (&((_2dt4pit){(a),(b),(c),(d),(e),(f)}))
#define C_DTPITIT2PDT(a,b,c,d,e) (&((dtpitit2pdt){(a),(b),(c),(d),(e)}))
#define C_IT3PI(a,b,c,d) (&((it3pi){(a),(b),(c),(d)}))
#define C_DTPIPDTPB(a,b,c,d) (&((dtpipdtpb){(a),(b),(c),(d)}))
#define C_DT2PIPB(a,b,c,d) (&((dt2pipb){(a),(b),(c),(d)}))
#define C_DT3PDT(a,b,c,d) (&((dt3pdt){(a),(b),(c),(d)}))
#define C_3DT2PDT(a,b,c,d,e) (&((_3dt2pdt){(a),(b),(c),(d),(e)}))
#define C_IT2PITDT(a,b,c,d) (&((it2pitdt){(a),(b),(c),(d)}))
#define C_ITPIT4DTPDTPIT(a,b,c,d,e,f,g,h) (&((itpit4dtpdtpit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_ITPIT3DTPDTPIT(a,b,c,d,e,f,g) (&((itpit3dtpdtpit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_ITPITDTPITDTPIT2IT(a,b,c,d,e,f,g,h) (&((itpitdtpitdtpit2it){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_IT6PIT(a,b,c,d,e,f,g) (&((it6pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DT2PIT3PDT(a,b,c,d,e,f) (&((dt2pit3pdt){(a),(b),(c),(d),(e),(f)}))
#define C_PIPDTPB(a,b,c) (&((pipdtpb){(a),(b),(c)}))
#define C_PITDTPIPIT(a,b,c,d) (&((pitdtpipit){(a),(b),(c),(d)}))
#define C_2PITITCH(a,b,c,d) (&((_2pititch){(a),(b),(c),(d)}))
#define C_2PITPCHIT(a,b,c,d) (&((_2pitpchit){(a),(b),(c),(d)}))
#define C_3PITPBPIT(a,b,c,d,e) (&((_3pitpbpit){(a),(b),(c),(d),(e)}))
#define C_PIT2DTPDTPIT(a,b,c,d,e) (&((pit2dtpdtpit){(a),(b),(c),(d),(e)}))
#define C_PIT3ITPDTPIT(a,b,c,d,e,f) (&((pit3itpdtpit){(a),(b),(c),(d),(e),(f)}))
#define C_2PITPB(a,b,c) (&((_2pitpb){(a),(b),(c)}))
#define C_DT4PDT(a,b,c,d,e) (&((dt4pdt){(a),(b),(c),(d),(e)}))
#define C_4DT4PDT(a,b,c,d,e,f,g,h) (&((_4dt4pdt){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_PDTPI3DTPI(a,b,c,d,e,f) (&((pdtpi3dtpi){(a),(b),(c),(d),(e),(f)}))
#define C_DT7PIT(a,b,c,d,e,f,g,h) (&((dt7pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DT4PITDT3PIT(a,b,c,d,e,f,g,h,i) (&((dt4pitdt3pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_DTPDT(a,b) (&((dtpdt){(a),(b)}))
#define C_3DTPIPDT(a,b,c,d,e) (&((_3dtpipdt){(a),(b),(c),(d),(e)}))
#define C_DTPITPIPDTPITPI(a,b,c,d,e,f) (&((dtpitpipdtpitpi){(a),(b),(c),(d),(e),(f)}))
#define C_2DT2IT2PIT(a,b,c,d,e,f) (&((_2dt2it2pit){(a),(b),(c),(d),(e),(f)}))
#define C_2DT4IT2PIT(a,b,c,d,e,f,g,h) (&((_2dt4it2pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DTPITPIIPITPI2PITDT4IT(a,b,c,d,e,f,g,h,i,j,k,l,m) (&((dtpitpiipitpi2pitdt4it){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m)}))
#define C_DT2PITIT2PIT(a,b,c,d,e,f) (&((dt2pitit2pit){(a),(b),(c),(d),(e),(f)}))
#define C_3DTPI(a,b,c,d) (&((_3dtpi){(a),(b),(c),(d)}))
#define C_2DTPDT2PITPDT(a,b,c,d,e,f) (&((_2dtpdt2pitpdt){(a),(b),(c),(d),(e),(f)}))
#define C_2DTITPITPI(a,b,c,d,e) (&((_2dtitpitpi){(a),(b),(c),(d),(e)}))
#define C_IT4PI(a,b,c,d,e) (&((it4pi){(a),(b),(c),(d),(e)}))
#define C_2I2PI(a,b,c,d) (&((_2i2pi){(a),(b),(c),(d)}))
#define C_2DTPI2PDTPIT(a,b,c,d,e,f) (&((_2dtpi2pdtpit){(a),(b),(c),(d),(e),(f)}))
#define C_I2PI(a,b,c) (&((i2pi){(a),(b),(c)}))
#define C_ITI(a,b) (&((iti){(a),(b)}))
#define C_2DTPDT(a,b,c) (&((_2dtpdt){(a),(b),(c)}))
#define C_PPPI3DT(a,b,c,d) (&((pppi3dt){(a),(b),(c),(d)}))
#define C_PPI2DT(a,b,c) (&((ppi2dt){(a),(b),(c)}))
#define C_2DTPI2DT(a,b,c,d,e) (&((_2dtpi2dt){(a),(b),(c),(d),(e)}))
#define C_2DT2PI(a,b,c,d) (&((_2dt2pi){(a),(b),(c),(d)}))
#define C_2DTPII2PDT(a,b,c,d,e,f) (&((_2dtpii2pdt){(a),(b),(c),(d),(e),(f)}))
#define C_2DTPI2I2PDT(a,b,c,d,e,f,g) (&((_2dtpi2i2pdt){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2DTPIPS(a,b,c,d) (&((_2dtpips){(a),(b),(c),(d)}))
#define C_DTI2PDT(a,b,c,d) (&((dti2pdt){(a),(b),(c),(d)}))
#define C_3I3PDT(a,b,c,d,e,f) (&((_3i3pdt){(a),(b),(c),(d),(e),(f)}))
#define C_3I2F(a,b,c,d,e) (&((_3i2f){(a),(b),(c),(d),(e)}))
#define C_3I2IT(a,b,c,d,e) (&((_3i2it){(a),(b),(c),(d),(e)}))
#define C_2DTPI2PDT(a,b,c,d,e) (&((_2dtpi2pdt){(a),(b),(c),(d),(e)}))
#define C_3DT2PI(a,b,c,d,e) (&((_3dt2pi){(a),(b),(c),(d),(e)}))
#define C_DT2PI(a,b,c) (&((dt2pi){(a),(b),(c)}))
#define C_2DT2I2PI(a,b,c,d,e,f) (&((_2dt2i2pi){(a),(b),(c),(d),(e),(f)}))
#define C_2DT2IPI(a,b,c,d,e) (&((_2dt2ipi){(a),(b),(c),(d),(e)}))
#define C_DTPIPDT(a,b,c) (&((dtpipdt){(a),(b),(c)}))
#define C_DTIPIIPII(a,b,c,d,e,f) (&((dtipiipii){(a),(b),(c),(d),(e),(f)}))
#define C_DTPII2PS(a,b,c,d,e) (&((dtpii2ps){(a),(b),(c),(d),(e)}))
#define C_2DTPI2DTPI(a,b,c,d,e,f) (&((_2dtpi2dtpi){(a),(b),(c),(d),(e),(f)}))
#define C_PDTPI(a,b) (&((pdtpi){(a),(b)}))
#define C_PDTPII(a,b,c) (&((pdtpii){(a),(b),(c)}))
#define C_DT2PII3PI(a,b,c,d,e,f,g) (&((dt2pii3pi){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DTPIPDT2PIPS(a,b,c,d,e,f) (&((dtpipdt2pips){(a),(b),(c),(d),(e),(f)}))
#define C_3PII(a,b,c,d) (&((_3pii){(a),(b),(c),(d)}))
#define C_DT4PI(a,b,c,d,e) (&((dt4pi){(a),(b),(c),(d),(e)}))
#define C_2PIPDT(a,b,c) (&((_2pipdt){(a),(b),(c)}))
#define C_DTPI2I(a,b,c,d) (&((dtpi2i){(a),(b),(c),(d)}))
#define C_2DT3PI(a,b,c,d,e) (&((_2dt3pi){(a),(b),(c),(d),(e)}))
#define C_2DTPIIPI(a,b,c,d,e) (&((_2dtpiipi){(a),(b),(c),(d),(e)}))
#define C_DT3PI(a,b,c,d) (&((dt3pi){(a),(b),(c),(d)}))
#define C_DTPIDTPDT2PI(a,b,c,d,e,f) (&((dtpidtpdt2pi){(a),(b),(c),(d),(e),(f)}))
#define C_DT3PII(a,b,c,d,e) (&((dt3pii){(a),(b),(c),(d),(e)}))
#define C_IDTPB(a,b,c) (&((idtpb){(a),(b),(c)}))
#define C_2DTPIT2DT(a,b,c,d,e) (&((_2dtpit2dt){(a),(b),(c),(d),(e)}))
#define C_DTPITIT2PI(a,b,c,d,e) (&((dtpitit2pi){(a),(b),(c),(d),(e)}))
#define C_PCDT(a,b) (&((pcdt){(a),(b)}))
#define C_2DTITPI(a,b,c,d) (&((_2dtitpi){(a),(b),(c),(d)}))
#define C_2DTPIPIT(a,b,c,d) (&((_2dtpipit){(a),(b),(c),(d)}))
#define C_DTPITI2PIT2PI(a,b,c,d,e,f,g) (&((dtpiti2pit2pi){(a),(b),(c),(d),(e),(f),(g)}))
#define C_PDT2IT(a,b,c) (&((pdt2it){(a),(b),(c)}))
#define C_DT3IT(a,b,c,d) (&((dt3it){(a),(b),(c),(d)}))
#define C_2DT4PI(a,b,c,d,e,f) (&((_2dt4pi){(a),(b),(c),(d),(e),(f)}))
#define C_ITPITPDT(a,b,c) (&((itpitpdt){(a),(b),(c)}))
#define C_DT2ITPI(a,b,c,d) (&((dt2itpi){(a),(b),(c),(d)}))
#define C_2DT3PITDT2PIT(a,b,c,d,e,f,g,h) (&((_2dt3pitdt2pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_9PLI4PIT(a,b,c,d,e,f,g,h,i,j,k,l,m) (&((_9pli4pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m)}))
#define C_DT2PITDT2PIDTPIT(a,b,c,d,e,f,g,h) (&((dt2pitdt2pidtpit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DTPIT3PI(a,b,c,d,e) (&((dtpit3pi){(a),(b),(c),(d),(e)}))
#define C_DT3PIDTPITDT3PI(a,b,c,d,e,f,g,h,i,j) (&((dt3pidtpitdt3pi){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)}))
#define C_2ITDTPITDT6PI(a,b,c,d,e,f,g,h,i,j,k) (&((_2itdtpitdt6pi){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)}))
#define C_2PITPI2PIT(a,b,c,d,e) (&((_2pitpi2pit){(a),(b),(c),(d),(e)}))
#define C_2ITDT3PITPDT(a,b,c,d,e,f,g) (&((_2itdt3pitpdt){(a),(b),(c),(d),(e),(f),(g)}))
#define C_4DT2PIT(a,b,c,d,e,f) (&((_4dt2pit){(a),(b),(c),(d),(e),(f)}))
#define C_3DT3PIT(a,b,c,d,e,f) (&((_3dt3pit){(a),(b),(c),(d),(e),(f)}))
#define C_5DT3PIT(a,b,c,d,e,f,g,h) (&((_5dt3pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_7DT3PIT(a,b,c,d,e,f,g,h,i,j) (&((_7dt3pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)}))
#define C_9DT3PIT(a,b,c,d,e,f,g,h,i,j,k,l) (&((_9dt3pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l)}))
#define C_2DT3PIPB(a,b,c,d,e,f) (&((_2dt3pipb){(a),(b),(c),(d),(e),(f)}))
#define C_3IPI(a,b,c,d) (&((_3ipi){(a),(b),(c),(d)}))
#define C_DTPIPB(a,b,c) (&((dtpipb){(a),(b),(c)}))
#define C_DTPIPB2PDT(a,b,c,d,e) (&((dtpipb2pdt){(a),(b),(c),(d),(e)}))
#define C_3I3PI(a,b,c,d,e,f) (&((_3i3pi){(a),(b),(c),(d),(e),(f)}))
#define C_DTPDT2PIPB(a,b,c,d,e) (&((dtpdt2pipb){(a),(b),(c),(d),(e)}))
#define C_DTPIPDTPI(a,b,c,d) (&((dtpipdtpi){(a),(b),(c),(d)}))
#define C_DT4PIPDT(a,b,c,d,e,f) (&((dt4pipdt){(a),(b),(c),(d),(e),(f)}))
#define C_DT4I2PI(a,b,c,d,e,f,g) (&((dt4i2pi){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2DTPITPI2PDTPITPI(a,b,c,d,e,f,g,h) (&((_2dtpitpi2pdtpitpi){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_2DTPDT2PI(a,b,c,d,e) (&((_2dtpdt2pi){(a),(b),(c),(d),(e)}))
#define C_DT4PIPB(a,b,c,d,e,f) (&((dt4pipb){(a),(b),(c),(d),(e),(f)}))
#define C_DT3PIPITPB(a,b,c,d,e,f) (&((dt3pipitpb){(a),(b),(c),(d),(e),(f)}))
#define C_DTPIT2PIPITITPB(a,b,c,d,e,f,g) (&((dtpit2pipititpb){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DTPI2DTPIPB(a,b,c,d,e,f) (&((dtpi2dtpipb){(a),(b),(c),(d),(e),(f)}))
#define C_DTPIT3PI2ITPB(a,b,c,d,e,f,g,h) (&((dtpit3pi2itpb){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DTPIT2PI2ITPB(a,b,c,d,e,f,g) (&((dtpit2pi2itpb){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DT3PIPB(a,b,c,d,e) (&((dt3pipb){(a),(b),(c),(d),(e)}))
#define C_2DT2ITI3PIT(a,b,c,d,e,f,g,h) (&((_2dt2iti3pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_2DT2ITI4PIT(a,b,c,d,e,f,g,h,i) (&((_2dt2iti4pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DTIT2PIT(a,b,c,d) (&((dtit2pit){(a),(b),(c),(d)}))
#define C_DT2ITDTPIT(a,b,c,d,e) (&((dt2itdtpit){(a),(b),(c),(d),(e)}))
#define C_DT2ITDTPITPDT2PIT(a,b,c,d,e,f,g,h) (&((dt2itdtpitpdt2pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DT2ITDT(a,b,c,d) (&((dt2itdt){(a),(b),(c),(d)}))
#define C_PIT2DTPITPDTDT(a,b,c,d,e,f) (&((pit2dtpitpdtdt){(a),(b),(c),(d),(e),(f)}))
#define C_PIT2DT2PIT2DT4PIT(a,b,c,d,e,f,g,h,i,j,k) (&((pit2dt2pit2dt4pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)}))
#define C_PIT4DTPIT2DT2PITDT(a,b,c,d,e,f,g,h,i,j,k) (&((pit4dtpit2dt2pitdt){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)}))
#define C_PIT4DTPDTPIT(a,b,c,d,e,f,g) (&((pit4dtpdtpit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_PIT4DTPIPIT(a,b,c,d,e,f,g) (&((pit4dtpipit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_PIT4DTPDT(a,b,c,d,e,f) (&((pit4dtpdt){(a),(b),(c),(d),(e),(f)}))
#define C_PIT4DTPDTPITDT(a,b,c,d,e,f,g,h) (&((pit4dtpdtpitdt){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_PIT2DTPI2PITDT(a,b,c,d,e,f,g) (&((pit2dtpi2pitdt){(a),(b),(c),(d),(e),(f),(g)}))
#define C_PIT2DTPDT(a,b,c,d) (&((pit2dtpdt){(a),(b),(c),(d)}))
#define C_PIT3DTPITPIPITDT(a,b,c,d,e,f,g,h) (&((pit3dtpitpipitdt){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_PIT3DT7PITDT(a,b,c,d,e,f,g,h,i,j,k,l) (&((pit3dt7pitdt){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l)}))
#define C_PIT3DT3PITDTPITDTPITDT(a,b,c,d,e,f,g,h,i,j,k,l) (&((pit3dt3pitdtpitdtpitdt){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l)}))
#define C_DTPCXDTPCX(a,b,c,d) (&((dtpcxdtpcx){(a),(b),(c),(d)}))
#define C_2DT3PCX(a,b,c,d,e) (&((_2dt3pcx){(a),(b),(c),(d),(e)}))
#define C_3DT3PCX(a,b,c,d,e,f) (&((_3dt3pcx){(a),(b),(c),(d),(e),(f)}))
#define C_2DTITPCXIT2PCX(a,b,c,d,e,f,g) (&((_2dtitpcxit2pcx){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2DTITPCX(a,b,c,d) (&((_2dtitpcx){(a),(b),(c),(d)}))
#define C_2DTPCX(a,b,c) (&((_2dtpcx){(a),(b),(c)}))
#define C_DTPCX(a,b) (&((dtpcx){(a),(b)}))
#define C_2DT2PITIT3PIT(a,b,c,d,e,f,g,h) (&((_2dt2pitit3pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_2DTPITPI3PIT(a,b,c,d,e,f,g) (&((_2dtpitpi3pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2DT4PIT2IT(a,b,c,d,e,f,g,h) (&((_2dt4pit2it){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_2DT3PPIT(a,b,c,d,e) (&((_2dt3ppit){(a),(b),(c),(d),(e)}))
#define C_2DTPIPDT(a,b,c,d) (&((_2dtpipdt){(a),(b),(c),(d)}))
#define C_DT4IT2FITPDTPIT(a,b,c,d,e,f,g,h,i) (&((dt4it2fitpdtpit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_4ITDT5PIT(a,b,c,d,e,f,g,h,i,j) (&((_4itdt5pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)}))
#define C_2ITDT2PITIT3PIT(a,b,c,d,e,f,g,h,i) (&((_2itdt2pitit3pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_ITDT8PIT(a,b,c,d,e,f,g,h,i,j) (&((itdt8pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)}))
#define C_2ITDT8PIT(a,b,c,d,e,f,g,h,i,j,k) (&((_2itdt8pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)}))
#define C_ITPITIT2PIT2DT2PIT(a,b,c,d,e,f,g,h,i) (&((itpitit2pit2dt2pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_2ITDTITDT2IT2PIT(a,b,c,d,e,f,g,h,i) (&((_2itdtitdt2it2pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_2DT2ITPPIT2DT2PPIT(a,b,c,d,e,f,g,h,i) (&((_2dt2itppit2dt2ppit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_3DTPIT(a,b,c,d) (&((_3dtpit){(a),(b),(c),(d)}))
#define C_BDTPI(a,b,c) (&((bdtpi){(a),(b),(c)}))
#define C_DTS2IT2PIT(a,b,c,d,e,f) (&((dts2it2pit){(a),(b),(c),(d),(e),(f)}))
#define C_3DT2PI2PIT(a,b,c,d,e,f,g) (&((_3dt2pi2pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DT2PITDTITPITPDT(a,b,c,d,e,f,g) (&((dt2pitdtitpitpdt){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DT4IT(a,b,c,d,e) (&((dt4it){(a),(b),(c),(d),(e)}))
#define C_DTPITPB(a,b,c) (&((dtpitpb){(a),(b),(c)}))
#define C_PITPIPITPI(a,b,c,d) (&((pitpipitpi){(a),(b),(c),(d)}))
#define C_ITPIPITPI(a,b,c,d) (&((itpipitpi){(a),(b),(c),(d)}))
#define C_2ITPI(a,b,c) (&((_2itpi){(a),(b),(c)}))
#define C_ITIPI(a,b,c) (&((itipi){(a),(b),(c)}))
#define C_IDTIT(a,b,c) (&((idtit){(a),(b),(c)}))
#define C_3IT2I(a,b,c,d,e) (&((_3it2i){(a),(b),(c),(d),(e)}))
#define C_3ITI(a,b,c,d) (&((_3iti){(a),(b),(c),(d)}))
#define C_I2IT(a,b,c) (&((i2it){(a),(b),(c)}))
#define C_DTPIIPIT(a,b,c,d) (&((dtpiipit){(a),(b),(c),(d)}))
#define C_2DTPIT2PDT(a,b,c,d,e) (&((_2dtpit2pdt){(a),(b),(c),(d),(e)}))
#define C_2DTPIIPIT(a,b,c,d,e) (&((_2dtpiipit){(a),(b),(c),(d),(e)}))
#define C_2DTIPIT(a,b,c,d) (&((_2dtipit){(a),(b),(c),(d)}))
#define C_2DTPITDTIT2PI(a,b,c,d,e,f,g) (&((_2dtpitdtit2pi){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DTIPDTDT2PI(a,b,c,d,e,f) (&((dtipdtdt2pi){(a),(b),(c),(d),(e),(f)}))
#define C_PPIT2DT(a,b,c) (&((ppit2dt){(a),(b),(c)}))
#define C_DTPIDTPITPIPIT(a,b,c,d,e,f) (&((dtpidtpitpipit){(a),(b),(c),(d),(e),(f)}))
#define C_DTPITPIITPDTPITPI(a,b,c,d,e,f,g) (&((dtpitpiitpdtpitpi){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DTPITPIPDTPITPI(a,b,c,d,e,f) (&((dtpitpipdtpitpi){(a),(b),(c),(d),(e),(f)}))
#define C_PDTPITPIIT(a,b,c,d) (&((pdtpitpiit){(a),(b),(c),(d)}))
#define C_DTPITPIIT3PDT(a,b,c,d,e,f,g) (&((dtpitpiit3pdt){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DTPITPI2IT2PI(a,b,c,d,e,f,g) (&((dtpitpi2it2pi){(a),(b),(c),(d),(e),(f),(g)}))
#define C_PDTPITPI(a,b,c) (&((pdtpitpi){(a),(b),(c)}))
#define C_DTPIT2PDT(a,b,c,d) (&((dtpit2pdt){(a),(b),(c),(d)}))
#define C_DTPIT2IT3PIT(a,b,c,d,e,f,g) (&((dtpit2it3pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DTPITDTIT2PI(a,b,c,d,e,f) (&((dtpitdtit2pi){(a),(b),(c),(d),(e),(f)}))
#define C_DTPITITPDT(a,b,c,d) (&((dtpititpdt){(a),(b),(c),(d)}))
#define C_DTPITITDTPDTPITPI(a,b,c,d,e,f,g) (&((dtpititdtpdtpitpi){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DT2PITPDT(a,b,c,d) (&((dt2pitpdt){(a),(b),(c),(d)}))
#define C_DT2PITPDTPI(a,b,c,d,e) (&((dt2pitpdtpi){(a),(b),(c),(d),(e)}))
#define C_DTPB(a,b) (&((dtpb){(a),(b)}))
#define C_DTPIT2PI2IT(a,b,c,d,e,f) (&((dtpit2pi2it){(a),(b),(c),(d),(e),(f)}))
#define C_DT2PITITDTPIT(a,b,c,d,e,f) (&((dt2pititdtpit){(a),(b),(c),(d),(e),(f)}))
#define C_DT3PITITDT2PIT(a,b,c,d,e,f,g,h) (&((dt3pititdt2pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_2DT2PITITDTPIT(a,b,c,d,e,f,g) (&((_2dt2pititdtpit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_I2DT2PIT(a,b,c,d,e) (&((i2dt2pit){(a),(b),(c),(d),(e)}))
#define C_IDT(a,b) (&((idt){(a),(b)}))
#define C_3DTPIPIT(a,b,c,d,e) (&((_3dtpipit){(a),(b),(c),(d),(e)}))
#define C_DT2ITFIT(a,b,c,d) (&((dt2itfit){(a),(b),(c),(d)}))
#define C_DTPITDTPIIPIT(a,b,c,d,e,f) (&((dtpitdtpiipit){(a),(b),(c),(d),(e),(f)}))
#define C_4ITDT2PIT(a,b,c,d,e,f,g) (&((_4itdt2pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DT2PITDT2PIDT2PIT(a,b,c,d,e,f,g,h,i) (&((dt2pitdt2pidt2pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_7DTPDTPIT(a,b,c,d,e,f,g,h,i) (&((_7dtpdtpit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_ITPIT3DT(a,b,c,d,e) (&((itpit3dt){(a),(b),(c),(d),(e)}))
#define C_DT3IT4PIT(a,b,c,d,e,f,g,h) (&((dt3it4pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DT3PITPI(a,b,c,d,e) (&((dt3pitpi){(a),(b),(c),(d),(e)}))
#define C_DTPI2PITPI(a,b,c,d,e) (&((dtpi2pitpi){(a),(b),(c),(d),(e)}))
#define C_3ITDTPIT2DT2PIT(a,b,c,d,e,f,g,h,i) (&((_3itdtpit2dt2pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_2ITDTPITDT2PIT(a,b,c,d,e,f,g) (&((_2itdtpitdt2pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_4IT2FITPI(a,b,c,d,e,f,g) (&((_4it2fitpi){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DT5PIT(a,b,c,d,e,f) (&((dt5pit){(a),(b),(c),(d),(e),(f)}))
#define C_PITDT2IT4PIT(a,b,c,d,e,f,g,h) (&((pitdt2it4pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DTPIPSPIB(a,b,c,d,e) (&((dtpipspib){(a),(b),(c),(d),(e)}))
#define C_2DTPITPI2DTPIPITPI(a,b,c,d,e,f,g,h,i) (&((_2dtpitpi2dtpipitpi){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_DTPIDTPIT2PIDTPIT(a,b,c,d,e,f,g,h) (&((dtpidtpit2pidtpit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DTI2PIT(a,b,c,d) (&((dti2pit){(a),(b),(c),(d)}))
#define C_3ITPI(a,b,c,d) (&((_3itpi){(a),(b),(c),(d)}))
#define C_2DT4IT(a,b,c,d,e,f) (&((_2dt4it){(a),(b),(c),(d),(e),(f)}))
#define C_DTPITPIDTPI2DT2PITDTIT(a,b,c,d,e,f,g,h,i,j,k) (&((dtpitpidtpi2dt2pitdtit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)}))
#define C_2DT2PITDT2PIT(a,b,c,d,e,f,g) (&((_2dt2pitdt2pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2DT2ITDT2PIT(a,b,c,d,e,f,g) (&((_2dt2itdt2pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2PITPI2DTPI4DT2ITDT2IT(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) (&((_2pitpi2dtpi4dt2itdt2it){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o)}))
#define C_2PITPI2DTPIDTPIT3DT2ITDT(a,b,c,d,e,f,g,h,i,j,k,l,m,n) (&((_2pitpi2dtpidtpit3dt2itdt){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n)}))
#define C_DTPITPIDTPI2DT4ITDT(a,b,c,d,e,f,g,h,i,j,k,l) (&((dtpitpidtpi2dt4itdt){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l)}))
#define C_DTPITPIDTPI2DT2ITDT(a,b,c,d,e,f,g,h,i,j) (&((dtpitpidtpi2dt2itdt){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)}))
#define C_DTPITDT4ITDT(a,b,c,d,e,f,g,h) (&((dtpitdt4itdt){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_3DTITPITDT3PIT(a,b,c,d,e,f,g,h,i) (&((_3dtitpitdt3pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_2C3DTITPITDTPITDTITPITDT(a,b,c,d,e,f,g,h,i,j,k,l,m) (&((_2c3dtitpitdtpitdtitpitdt){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m)}))
#define C_3DT6PIT(a,b,c,d,e,f,g,h,i) (&((_3dt6pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_3DT6PI(a,b,c,d,e,f,g,h,i) (&((_3dt6pi){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_DTIPIPB(a,b,c,d) (&((dtipipb){(a),(b),(c),(d)}))
#define C_DTI2PI(a,b,c,d) (&((dti2pi){(a),(b),(c),(d)}))
#define C_DT5PI(a,b,c,d,e,f) (&((dt5pi){(a),(b),(c),(d),(e),(f)}))
#define C_DT6PI(a,b,c,d,e,f,g) (&((dt6pi){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DT2PIPDTPI(a,b,c,d,e) (&((dt2pipdtpi){(a),(b),(c),(d),(e)}))
#define C_DT3PITDT3PIT(a,b,c,d,e,f,g,h) (&((dt3pitdt3pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_2IT2PIT(a,b,c,d) (&((_2it2pit){(a),(b),(c),(d)}))
#define C_ITPIT(a,b) (&((itpit){(a),(b)}))
#define C_DT5PITDT(a,b,c,d,e,f,g) (&((dt5pitdt){(a),(b),(c),(d),(e),(f),(g)}))
#define C_DTPIT2DTPI3PIT(a,b,c,d,e,f,g,h) (&((dtpit2dtpi3pit){(a),(b),(c),(d),(e),(f),(g),(h)}))
#define C_DTPIT2DTPI4PIT(a,b,c,d,e,f,g,h,i) (&((dtpit2dtpi4pit){(a),(b),(c),(d),(e),(f),(g),(h),(i)}))
#define C_DTPIT2DTPI5PIT(a,b,c,d,e,f,g,h,i,j) (&((dtpit2dtpi5pit){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)}))
#define C_3I5PIDT2PI(a,b,c,d,e,f,g,h,i,j,k) (&((_3i5pidt2pi){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)}))
#define C_PIPIT(a,b) (&((pipit){(a),(b)}))
#define C_DTPIT2DTPIPIT(a,b,c,d,e,f) (&((dtpit2dtpipit){(a),(b),(c),(d),(e),(f)}))
#define C_DTPIT2DT2PIPITPIPITPI(a,b,c,d,e,f,g,h,i,j) (&((dtpit2dt2pipitpipitpi){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)}))
#define C_2DTPII(a,b,c,d) (&((_2dtpii){(a),(b),(c),(d)}))
#define C_2DT2PIDTPI(a,b,c,d,e,f) (&((_2dt2pidtpi){(a),(b),(c),(d),(e),(f)}))
#define C_DTPIT2DT2PIPITPI3PIT2PI(a,b,c,d,e,f,g,h,i,j,k,l,m) (&((dtpit2dt2pipitpi3pit2pi){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m)}))
#define C_PDTPI3PIT(a,b,c,d,e) (&((pdtpi3pit){(a),(b),(c),(d),(e)}))
#define C_2PDTPU(a,b,c) (&((_2pdtpu){(a),(b),(c)}))
#define C_3PDTPU(a,b,c,d) (&((_3pdtpu){(a),(b),(c),(d)}))
#define C_2PDT4PIT(a,b,c,d,e,f) (&((_2pdt4pit){(a),(b),(c),(d),(e),(f)}))
#define C_2PDTPI(a,b,c) (&((_2pdtpi){(a),(b),(c)}))
#define C_2PDTPS(a,b,c) (&((_2pdtps){(a),(b),(c)}))
#define C_2PDTPS4PIT(a,b,c,d,e,f,g) (&((_2pdtps4pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2DT2PI3PIT(a,b,c,d,e,f,g) (&((_2dt2pi3pit){(a),(b),(c),(d),(e),(f),(g)}))
#define C_2DT3PDT(a,b,c,d,e) (&((_2dt3pdt){(a),(b),(c),(d),(e)}))
#define C_3DTPI2PDT(a,b,c,d,e,f) (&((_3dtpi2pdt){(a),(b),(c),(d),(e),(f)}))
#define C_SIT(a,b) (&((sit){(a),(b)}))
#define C_DTS2PIT(a,b,c,d) (&((dts2pit){(a),(b),(c),(d)}))
#define C_2IT3PIT(a,b,c,d,e) (&((_2it3pit){(a),(b),(c),(d),(e)}))
#define C_2PDTPIT(a,b,c) (&((_2pdtpit){(a),(b),(c)}))
#define C_PDT6PIT(a,b,c,d,e,f,g) (&((pdt6pit){(a),(b),(c),(d),(e),(f),(g)}))

#define C_PINT2(a,b) ((int []){(a),(b)})
#define C_PINT3(a,b,c) ((int []){(a),(b),(c)})
#define C_PINT4(a,b,c,d) ((int []){(a),(b),(c),(d)})
#define C_PINT5(a,b,c,d,e) ((int []){(a),(b),(c),(d),(e)})
#define C_PINT6(a,b,c,d,e,f) ((int []){(a),(b),(c),(d),(e),(f)})
#define C_PINT7(a,b,c,d,e,f,g) ((int []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PINT8(a,b,c,d,e,f,g,h) ((int []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PINT9(a,b,c,d,e,f,g,h,i) ((int []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PINT10(a,b,c,d,e,f,g,h,i,j) ((int []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPINT2(a,b) ((int *[]){(a),(b)})
#define C_PPINT3(a,b,c) ((int *[]){(a),(b),(c)})
#define C_PPINT4(a,b,c,d) ((int *[]){(a),(b),(c),(d)})
#define C_PPINT5(a,b,c,d,e) ((int *[]){(a),(b),(c),(d),(e)})
#define C_PPINT6(a,b,c,d,e,f) ((int *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPINT7(a,b,c,d,e,f,g) ((int *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPINT8(a,b,c,d,e,f,g,h) ((int *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPINT9(a,b,c,d,e,f,g,h,i) ((int *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPINT10(a,b,c,d,e,f,g,h,i,j) ((int *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PUINT2(a,b) ((unsigned int []){(a),(b)})
#define C_PUINT3(a,b,c) ((unsigned int []){(a),(b),(c)})
#define C_PUINT4(a,b,c,d) ((unsigned int []){(a),(b),(c),(d)})
#define C_PUINT5(a,b,c,d,e) ((unsigned int []){(a),(b),(c),(d),(e)})
#define C_PUINT6(a,b,c,d,e,f) ((unsigned int []){(a),(b),(c),(d),(e),(f)})
#define C_PUINT7(a,b,c,d,e,f,g) ((unsigned int []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PUINT8(a,b,c,d,e,f,g,h) ((unsigned int []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PUINT9(a,b,c,d,e,f,g,h,i) ((unsigned int []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PUINT10(a,b,c,d,e,f,g,h,i,j) ((unsigned int []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPUINT2(a,b) ((unsigned int *[]){(a),(b)})
#define C_PPUINT3(a,b,c) ((unsigned int *[]){(a),(b),(c)})
#define C_PPUINT4(a,b,c,d) ((unsigned int *[]){(a),(b),(c),(d)})
#define C_PPUINT5(a,b,c,d,e) ((unsigned int *[]){(a),(b),(c),(d),(e)})
#define C_PPUINT6(a,b,c,d,e,f) ((unsigned int *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPUINT7(a,b,c,d,e,f,g) ((unsigned int *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPUINT8(a,b,c,d,e,f,g,h) ((unsigned int *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPUINT9(a,b,c,d,e,f,g,h,i) ((unsigned int *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPUINT10(a,b,c,d,e,f,g,h,i,j) ((unsigned int *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PSHRT2(a,b) ((short []){(a),(b)})
#define C_PSHRT3(a,b,c) ((short []){(a),(b),(c)})
#define C_PSHRT4(a,b,c,d) ((short []){(a),(b),(c),(d)})
#define C_PSHRT5(a,b,c,d,e) ((short []){(a),(b),(c),(d),(e)})
#define C_PSHRT6(a,b,c,d,e,f) ((short []){(a),(b),(c),(d),(e),(f)})
#define C_PSHRT7(a,b,c,d,e,f,g) ((short []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PSHRT8(a,b,c,d,e,f,g,h) ((short []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PSHRT9(a,b,c,d,e,f,g,h,i) ((short []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PSHRT10(a,b,c,d,e,f,g,h,i,j) ((short []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPSHRT2(a,b) ((short *[]){(a),(b)})
#define C_PPSHRT3(a,b,c) ((short *[]){(a),(b),(c)})
#define C_PPSHRT4(a,b,c,d) ((short *[]){(a),(b),(c),(d)})
#define C_PPSHRT5(a,b,c,d,e) ((short *[]){(a),(b),(c),(d),(e)})
#define C_PPSHRT6(a,b,c,d,e,f) ((short *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPSHRT7(a,b,c,d,e,f,g) ((short *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPSHRT8(a,b,c,d,e,f,g,h) ((short *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPSHRT9(a,b,c,d,e,f,g,h,i) ((short *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPSHRT10(a,b,c,d,e,f,g,h,i,j) ((short *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PUSHRT2(a,b) ((unsigned short []){(a),(b)})
#define C_PUSHRT3(a,b,c) ((unsigned short []){(a),(b),(c)})
#define C_PUSHRT4(a,b,c,d) ((unsigned short []){(a),(b),(c),(d)})
#define C_PUSHRT5(a,b,c,d,e) ((unsigned short []){(a),(b),(c),(d),(e)})
#define C_PUSHRT6(a,b,c,d,e,f) ((unsigned short []){(a),(b),(c),(d),(e),(f)})
#define C_PUSHRT7(a,b,c,d,e,f,g) ((unsigned short []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PUSHRT8(a,b,c,d,e,f,g,h) ((unsigned short []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PUSHRT9(a,b,c,d,e,f,g,h,i) ((unsigned short []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PUSHRT10(a,b,c,d,e,f,g,h,i,j) ((unsigned short []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPUSHRT2(a,b) ((unsigned short *[]){(a),(b)})
#define C_PPUSHRT3(a,b,c) ((unsigned short *[]){(a),(b),(c)})
#define C_PPUSHRT4(a,b,c,d) ((unsigned short *[]){(a),(b),(c),(d)})
#define C_PPUSHRT5(a,b,c,d,e) ((unsigned short *[]){(a),(b),(c),(d),(e)})
#define C_PPUSHRT6(a,b,c,d,e,f) ((unsigned short *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPUSHRT7(a,b,c,d,e,f,g) ((unsigned short *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPUSHRT8(a,b,c,d,e,f,g,h) ((unsigned short *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPUSHRT9(a,b,c,d,e,f,g,h,i) ((unsigned short *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPUSHRT10(a,b,c,d,e,f,g,h,i,j) ((unsigned short *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PCHR2(a,b) ((char []){(a),(b)})
#define C_PCHR3(a,b,c) ((char []){(a),(b),(c)})
#define C_PCHR4(a,b,c,d) ((char []){(a),(b),(c),(d)})
#define C_PCHR5(a,b,c,d,e) ((char []){(a),(b),(c),(d),(e)})
#define C_PCHR6(a,b,c,d,e,f) ((char []){(a),(b),(c),(d),(e),(f)})
#define C_PCHR7(a,b,c,d,e,f,g) ((char []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PCHR8(a,b,c,d,e,f,g,h) ((char []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PCHR9(a,b,c,d,e,f,g,h,i) ((char []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PCHR10(a,b,c,d,e,f,g,h,i,j) ((char []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPCHR2(a,b) ((char *[]){(a),(b)})
#define C_PPCHR3(a,b,c) ((char *[]){(a),(b),(c)})
#define C_PPCHR4(a,b,c,d) ((char *[]){(a),(b),(c),(d)})
#define C_PPCHR5(a,b,c,d,e) ((char *[]){(a),(b),(c),(d),(e)})
#define C_PPCHR6(a,b,c,d,e,f) ((char *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPCHR7(a,b,c,d,e,f,g) ((char *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPCHR8(a,b,c,d,e,f,g,h) ((char *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPCHR9(a,b,c,d,e,f,g,h,i) ((char *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPCHR10(a,b,c,d,e,f,g,h,i,j) ((char *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PUCHR2(a,b) ((unsigned char []){(a),(b)})
#define C_PUCHR3(a,b,c) ((unsigned char []){(a),(b),(c)})
#define C_PUCHR4(a,b,c,d) ((unsigned char []){(a),(b),(c),(d)})
#define C_PUCHR5(a,b,c,d,e) ((unsigned char []){(a),(b),(c),(d),(e)})
#define C_PUCHR6(a,b,c,d,e,f) ((unsigned char []){(a),(b),(c),(d),(e),(f)})
#define C_PUCHR7(a,b,c,d,e,f,g) ((unsigned char []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PUCHR8(a,b,c,d,e,f,g,h) ((unsigned char []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PUCHR9(a,b,c,d,e,f,g,h,i) ((unsigned char []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PUCHR10(a,b,c,d,e,f,g,h,i,j) ((unsigned char []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPUCHR2(a,b) ((unsigned char *[]){(a),(b)})
#define C_PPUCHR3(a,b,c) ((unsigned char *[]){(a),(b),(c)})
#define C_PPUCHR4(a,b,c,d) ((unsigned char *[]){(a),(b),(c),(d)})
#define C_PPUCHR5(a,b,c,d,e) ((unsigned char *[]){(a),(b),(c),(d),(e)})
#define C_PPUCHR6(a,b,c,d,e,f) ((unsigned char *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPUCHR7(a,b,c,d,e,f,g) ((unsigned char *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPUCHR8(a,b,c,d,e,f,g,h) ((unsigned char *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPUCHR9(a,b,c,d,e,f,g,h,i) ((unsigned char *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPUCHR10(a,b,c,d,e,f,g,h,i,j) ((unsigned char *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PLNG2(a,b) ((long []){(a),(b)})
#define C_PLNG3(a,b,c) ((long []){(a),(b),(c)})
#define C_PLNG4(a,b,c,d) ((long []){(a),(b),(c),(d)})
#define C_PLNG5(a,b,c,d,e) ((long []){(a),(b),(c),(d),(e)})
#define C_PLNG6(a,b,c,d,e,f) ((long []){(a),(b),(c),(d),(e),(f)})
#define C_PLNG7(a,b,c,d,e,f,g) ((long []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PLNG8(a,b,c,d,e,f,g,h) ((long []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PLNG9(a,b,c,d,e,f,g,h,i) ((long []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PLNG10(a,b,c,d,e,f,g,h,i,j) ((long []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPLNG2(a,b) ((long *[]){(a),(b)})
#define C_PPLNG3(a,b,c) ((long *[]){(a),(b),(c)})
#define C_PPLNG4(a,b,c,d) ((long *[]){(a),(b),(c),(d)})
#define C_PPLNG5(a,b,c,d,e) ((long *[]){(a),(b),(c),(d),(e)})
#define C_PPLNG6(a,b,c,d,e,f) ((long *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPLNG7(a,b,c,d,e,f,g) ((long *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPLNG8(a,b,c,d,e,f,g,h) ((long *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPLNG9(a,b,c,d,e,f,g,h,i) ((long *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPLNG10(a,b,c,d,e,f,g,h,i,j) ((long *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PULNG2(a,b) ((unsigned long []){(a),(b)})
#define C_PULNG3(a,b,c) ((unsigned long []){(a),(b),(c)})
#define C_PULNG4(a,b,c,d) ((unsigned long []){(a),(b),(c),(d)})
#define C_PULNG5(a,b,c,d,e) ((unsigned long []){(a),(b),(c),(d),(e)})
#define C_PULNG6(a,b,c,d,e,f) ((unsigned long []){(a),(b),(c),(d),(e),(f)})
#define C_PULNG7(a,b,c,d,e,f,g) ((unsigned long []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PULNG8(a,b,c,d,e,f,g,h) ((unsigned long []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PULNG9(a,b,c,d,e,f,g,h,i) ((unsigned long []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PULNG10(a,b,c,d,e,f,g,h,i,j) ((unsigned long []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPULNG2(a,b) ((unsigned long *[]){(a),(b)})
#define C_PPULNG3(a,b,c) ((unsigned long *[]){(a),(b),(c)})
#define C_PPULNG4(a,b,c,d) ((unsigned long *[]){(a),(b),(c),(d)})
#define C_PPULNG5(a,b,c,d,e) ((unsigned long *[]){(a),(b),(c),(d),(e)})
#define C_PPULNG6(a,b,c,d,e,f) ((unsigned long *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPULNG7(a,b,c,d,e,f,g) ((unsigned long *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPULNG8(a,b,c,d,e,f,g,h) ((unsigned long *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPULNG9(a,b,c,d,e,f,g,h,i) ((unsigned long *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPULNG10(a,b,c,d,e,f,g,h,i,j) ((unsigned long *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PLLNG2(a,b) ((long long []){(a),(b)})
#define C_PLLNG3(a,b,c) ((long long []){(a),(b),(c)})
#define C_PLLNG4(a,b,c,d) ((long long []){(a),(b),(c),(d)})
#define C_PLLNG5(a,b,c,d,e) ((long long []){(a),(b),(c),(d),(e)})
#define C_PLLNG6(a,b,c,d,e,f) ((long long []){(a),(b),(c),(d),(e),(f)})
#define C_PLLNG7(a,b,c,d,e,f,g) ((long long []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PLLNG8(a,b,c,d,e,f,g,h) ((long long []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PLLNG9(a,b,c,d,e,f,g,h,i) ((long long []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PLLNG10(a,b,c,d,e,f,g,h,i,j) ((long long []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPLLNG2(a,b) ((long long *[]){(a),(b)})
#define C_PPLLNG3(a,b,c) ((long long *[]){(a),(b),(c)})
#define C_PPLLNG4(a,b,c,d) ((long long *[]){(a),(b),(c),(d)})
#define C_PPLLNG5(a,b,c,d,e) ((long long *[]){(a),(b),(c),(d),(e)})
#define C_PPLLNG6(a,b,c,d,e,f) ((long long *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPLLNG7(a,b,c,d,e,f,g) ((long long *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPLLNG8(a,b,c,d,e,f,g,h) ((long long *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPLLNG9(a,b,c,d,e,f,g,h,i) ((long long *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPLLNG10(a,b,c,d,e,f,g,h,i,j) ((long long *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PULLNG2(a,b) ((unsigned long long []){(a),(b)})
#define C_PULLNG3(a,b,c) ((unsigned long long []){(a),(b),(c)})
#define C_PULLNG4(a,b,c,d) ((unsigned long long []){(a),(b),(c),(d)})
#define C_PULLNG5(a,b,c,d,e) ((unsigned long long []){(a),(b),(c),(d),(e)})
#define C_PULLNG6(a,b,c,d,e,f) ((unsigned long long []){(a),(b),(c),(d),(e),(f)})
#define C_PULLNG7(a,b,c,d,e,f,g) ((unsigned long long []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PULLNG8(a,b,c,d,e,f,g,h) ((unsigned long long []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PULLNG9(a,b,c,d,e,f,g,h,i) ((unsigned long long []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PULLNG10(a,b,c,d,e,f,g,h,i,j) ((unsigned long long []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPULLNG2(a,b) ((unsigned long long *[]){(a),(b)})
#define C_PPULLNG3(a,b,c) ((unsigned long long *[]){(a),(b),(c)})
#define C_PPULLNG4(a,b,c,d) ((unsigned long long *[]){(a),(b),(c),(d)})
#define C_PPULLNG5(a,b,c,d,e) ((unsigned long long *[]){(a),(b),(c),(d),(e)})
#define C_PPULLNG6(a,b,c,d,e,f) ((unsigned long long *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPULLNG7(a,b,c,d,e,f,g) ((unsigned long long *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPULLNG8(a,b,c,d,e,f,g,h) ((unsigned long long *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPULLNG9(a,b,c,d,e,f,g,h,i) ((unsigned long long *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPULLNG10(a,b,c,d,e,f,g,h,i,j) ((unsigned long long *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PFLT2(a,b) ((float []){(a),(b)})
#define C_PFLT3(a,b,c) ((float []){(a),(b),(c)})
#define C_PFLT4(a,b,c,d) ((float []){(a),(b),(c),(d)})
#define C_PFLT5(a,b,c,d,e) ((float []){(a),(b),(c),(d),(e)})
#define C_PFLT6(a,b,c,d,e,f) ((float []){(a),(b),(c),(d),(e),(f)})
#define C_PFLT7(a,b,c,d,e,f,g) ((float []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PFLT8(a,b,c,d,e,f,g,h) ((float []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PFLT9(a,b,c,d,e,f,g,h,i) ((float []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PFLT10(a,b,c,d,e,f,g,h,i,j) ((float []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPFLT2(a,b) ((float *[]){(a),(b)})
#define C_PPFLT3(a,b,c) ((float *[]){(a),(b),(c)})
#define C_PPFLT4(a,b,c,d) ((float *[]){(a),(b),(c),(d)})
#define C_PPFLT5(a,b,c,d,e) ((float *[]){(a),(b),(c),(d),(e)})
#define C_PPFLT6(a,b,c,d,e,f) ((float *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPFLT7(a,b,c,d,e,f,g) ((float *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPFLT8(a,b,c,d,e,f,g,h) ((float *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPFLT9(a,b,c,d,e,f,g,h,i) ((float *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPFLT10(a,b,c,d,e,f,g,h,i,j) ((float *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PDBL2(a,b) ((double []){(a),(b)})
#define C_PDBL3(a,b,c) ((double []){(a),(b),(c)})
#define C_PDBL4(a,b,c,d) ((double []){(a),(b),(c),(d)})
#define C_PDBL5(a,b,c,d,e) ((double []){(a),(b),(c),(d),(e)})
#define C_PDBL6(a,b,c,d,e,f) ((double []){(a),(b),(c),(d),(e),(f)})
#define C_PDBL7(a,b,c,d,e,f,g) ((double []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PDBL8(a,b,c,d,e,f,g,h) ((double []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PDBL9(a,b,c,d,e,f,g,h,i) ((double []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PDBL10(a,b,c,d,e,f,g,h,i,j) ((double []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})
#define C_PDBL11(a,b,c,d,e,f,g,h,i,j,k) ((double []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k)})
#define C_PDBL12(a,b,c,d,e,f,g,h,i,j,k,l) ((double []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l)})


#define C_PPDBL2(a,b) ((double *[]){(a),(b)})
#define C_PPDBL3(a,b,c) ((double *[]){(a),(b),(c)})
#define C_PPDBL4(a,b,c,d) ((double *[]){(a),(b),(c),(d)})
#define C_PPDBL5(a,b,c,d,e) ((double *[]){(a),(b),(c),(d),(e)})
#define C_PPDBL6(a,b,c,d,e,f) ((double *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPDBL7(a,b,c,d,e,f,g) ((double *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPDBL8(a,b,c,d,e,f,g,h) ((double *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPDBL9(a,b,c,d,e,f,g,h,i) ((double *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPDBL10(a,b,c,d,e,f,g,h,i,j) ((double *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PLDBL2(a,b) ((long double []){(a),(b)})
#define C_PLDBL3(a,b,c) ((long double []){(a),(b),(c)})
#define C_PLDBL4(a,b,c,d) ((long double []){(a),(b),(c),(d)})
#define C_PLDBL5(a,b,c,d,e) ((long double []){(a),(b),(c),(d),(e)})
#define C_PLDBL6(a,b,c,d,e,f) ((long double []){(a),(b),(c),(d),(e),(f)})
#define C_PLDBL7(a,b,c,d,e,f,g) ((long double []){(a),(b),(c),(d),(e),(f),(g)})
#define C_PLDBL8(a,b,c,d,e,f,g,h) ((long double []){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PLDBL9(a,b,c,d,e,f,g,h,i) ((long double []){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PLDBL10(a,b,c,d,e,f,g,h,i,j) ((long double []){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

#define C_PPLDBL2(a,b) ((long double *[]){(a),(b)})
#define C_PPLDBL3(a,b,c) ((long double *[]){(a),(b),(c)})
#define C_PPLDBL4(a,b,c,d) ((long double *[]){(a),(b),(c),(d)})
#define C_PPLDBL5(a,b,c,d,e) ((long double *[]){(a),(b),(c),(d),(e)})
#define C_PPLDBL6(a,b,c,d,e,f) ((long double *[]){(a),(b),(c),(d),(e),(f)})
#define C_PPLDBL7(a,b,c,d,e,f,g) ((long double *[]){(a),(b),(c),(d),(e),(f),(g)})
#define C_PPLDBL8(a,b,c,d,e,f,g,h) ((long double *[]){(a),(b),(c),(d),(e),(f),(g),(h)})
#define C_PPLDBL9(a,b,c,d,e,f,g,h,i) ((long double *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i)})
#define C_PPLDBL10(a,b,c,d,e,f,g,h,i,j) ((long double *[]){(a),(b),(c),(d),(e),(f),(g),(h),(i),(j)})

typedef struct
{
	int a;
} sint;

typedef struct
{
	unsigned a;
} suint;

typedef struct
{
	short a;
} sshrt;

typedef struct
{
	unsigned short a;
} sushrt;

typedef struct
{
	char a;
} schr;

typedef struct
{
	unsigned char a;
} suchr;

typedef struct
{
	long a;
} slng;

typedef struct
{
	unsigned long a;
} sulng;

typedef struct
{
	long long a;
} sllng;

typedef struct
{
	unsigned long long a;
} sullng;

typedef struct
{
	float a;
} sflt;

typedef struct
{
	double a;
} sdbl;

typedef struct
{
	long double a;
} sldbl;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	ityp * a3;
} _2dt2pit;

typedef struct
{
	ityp a0;
	bool a1;
} itb;

typedef struct
{
	ityp * a0;
	int a1;
	int a2;
	ityp * a3;
	dim_typ * a4;
} pit2ipitpdt;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ * a2;
	ityp * a3;
	dim_typ * a4;
	dim_typ a5;
	ityp * a6;
} pitdtpdtpitpdtdtpit;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	ityp * a2;
	ityp * a3;
} _2pdt2pit;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	dim_typ * a6;
	dim_typ a7;
	dim_typ a8;
} pdt5pitpdt2dt;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
} pdt4pit;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
} pdt3pit;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	ityp * a2;
} pdt2pit;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	dim_typ * a2;
	dim_typ * a3;
	ityp * a4;
	bool * a5;
} _4pdtpitpb;

typedef struct
{
	bool a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
} b3dt;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	dim_typ * a2;
	dim_typ * a3;
	dim_typ * a4;
	ityp * a5;	
} _5pdtpit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
} dtpit;

typedef struct
{
	bool a0;
	bool a1;
	dim_typ a2;
	ityp * a3;
	dim_typ * a4;
	ityp * a5;
} _2bdtpitpdtpit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
} _2dtpit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ * a3;
	dim_typ * a4;
} dt2pit2pdt;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ * a3;
	dim_typ * a4;
	ityp * a5;
	ityp * a6;
} dt2pit2pdt2pit;

typedef struct
{
	dim_typ a0;
	int a1;
	int a2;
	int * a3;
	int * a4;
} dt2i2pi;

typedef struct
{
	ityp a0;
	ityp * a1;
} itpit;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp * a4;
} dt3itpit;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	dim_typ * a2;
	ityp * a3;
	ityp * a4;
} _3pdt2pit;

typedef struct
{
	dim_typ a0;
	ityp a1;
} dtit;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp * a2;
	dim_typ a3;
	ityp * a4;
	dim_typ a5;
} dtitpitdtpitdt;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int a2;
} dtpiti;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int a2;
	ityp * a3;
	int a4;
	ityp a5;
	ityp a6;
} dtpitipiti2it;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp * a2;
	int a3;
} dtitpiti;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int a2;
	ityp * a3;
	int a4;
} dtpitipiti;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	ityp a6;
	ityp (* a7)(ityp);
	ityp * a8;
} _7itfitpit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp (* a4)(ityp);
	ityp * a5;
} _4pitfitpit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp (* a4)(ityp);
} _4pitfit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp a2;
	ityp a3;
	ityp * a4;
	int * a5;
} _2dt2itpitpi;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
} dt2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
} dt3pit;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
} dtit3pit;

typedef struct
{
	ityp a0;
	ityp a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
} _2itdt3pit;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
	short * a7;
	ityp * a8;
} pdt6pitpspit;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	short * a5;
	ityp * a6;	
} pdt4pitpspit;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	short * a6;
	ityp * a7;	
} pdt5pitpspit;

typedef struct 
{
	short * a0;
	ityp * a1;
	ityp * a2;
	unsigned long * a3;
	unsigned long * a4;
} ps2pit2pul;

typedef struct
{
	int a0;
	short * a1;
	ityp * a2;
	ityp * a3;
	unsigned long * a4;
	unsigned long * a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
	ityp * a9;
	ityp * a10;
	ityp * a11;
	ityp * a12;
} ips2pit2pul7pit;

typedef struct
{
	int a0;
	short * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	unsigned long * a6;
	unsigned long * a7;
	ityp * a8;
	ityp * a9;
	ityp * a10;
	ityp * a11;
} ips4pit2pul4pit;

typedef struct 
{
	short * a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4; 
	unsigned long * a5;
	unsigned long * a6;
} ps4pit2pul;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
} pdtpit;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	ityp * a2;
	dim_typ * a3;
	ityp * a4;
} pdt2pitpdtpit;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	dim_typ * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
} _3pdt3pit;

typedef struct
{
	int a0;
	int a1;
	int * a2;
} _2ipi;

typedef struct
{
	dim_typ a0;
	int * a1;
} dtpi;

typedef struct
{
	dim_typ a0;
	ityp (* a1)(ityp);
	ityp * a2;
} dtfitpit;

typedef struct
{
	dim_typ a0;
	ityp (* a1)(ityp);
	ityp a2;
	ityp a3;
	ityp * a4;
} dtfit2itpit; 

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ a2;
	int * a3;
	ityp * a4;
} dtpitdtpipit;

typedef struct
{
	ityp a0;
	ityp * a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
} itpitdt4pit;  

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
} it3pit;

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	dim_typ a4;
	int * a5;
} it3pitdtpi;

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp a2;
	ityp a3;
	ityp * a4;
} itpit2itpit;

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	dim_typ a4;
	int * a5;
	ityp * a6;
	ityp * a7;
} it3pitdtpi2pit;

typedef struct
{
	ityp a0;
	ityp a1;
	dim_typ a2;
	int * a3;
	ityp * a4;
	ityp * a5;
} _2itdtpi2pit;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	int * a3;
	ityp * a4;
} dt2itpipit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
} _2dtpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	bool * a3; 
} _2dtpipb;

typedef struct
{
	dim_typ a0;
	int * a1;
	int a2;
} dtpii;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int * a3;
	int * a4;
	int * a5;
	int * a6;
	int * a7;
} _2dt6pi;

typedef struct
{
	ityp a0;
	ityp a1;
	dim_typ a2;
} _2itdt;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp a2;
} dtpitit;

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp * a2;
} it2pit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp * a2;
	ityp * a3;
	dim_typ * a4;
} _2it2pitpdt;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	ityp * a3;
} dt2itpit;

typedef struct
{
	ityp (* a0)(ityp); 
	ityp a1;
	ityp a2;
	dim_typ a3;
} fit2itdt;

typedef struct
{
	int a0;
	int * a1;
} ipi; 

typedef struct
{
	dim_typ (* a0)(dim_typ);
	dim_typ a1;
	dim_typ * a2;
	dim_typ * a3;
} fdtdt2pdt;

typedef struct
{
	ityp a0;
	dim_typ a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
} it2dt2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ a3;
	ityp * a4;
	short * a5;
} dt2pitdtpitps;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ a3;
	short * a4;
} dt2pitdtps;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
} _3dt2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ * a3;
	ityp * a4;
	ityp * a5;
} dt2pitpdt2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp a3;
	ityp a4;
	dim_typ * a5;
	ityp * a6;
	ityp * a7;
} dt2pit2itpdt2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
} dt4pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	int a3;
	ityp * a4;
	ityp * a5;
} dt2piti2pit;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	ityp * a3;
	ityp * a4;
} dt2it2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ * a2;
	ityp * a3;
	ityp * a4;
} dtpitpdt2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp a2;
	int * a3;
} dtpititpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	ityp * a3;
	ityp a4;
	int * a5;
} _2dt2pititpi;

typedef struct
{
	dim_typ a0;
	ityp (* a1)(ityp);
	ityp (* a2)(ityp);
	ityp (* a3)(ityp);
	ityp (* a4)(ityp);
	ityp * a5;
} dt4fitpit;

typedef struct
{
	dim_typ a0;
	ityp (* a1)(ityp);
	ityp (* a2)(ityp);
	ityp (* a3)(ityp);
	ityp * a4;
} dt3fitpit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp (* a3)(ityp);
} dt2pitfit;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp (* a5)(ityp);
	ityp (* a6)(ityp);
	ityp * a7;
} dt4it2fitpit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ a2;
	ityp (* a3)(ityp);
	ityp * a4;
	ityp * a5;
	ityp * a6;
} dtpitdtfit3pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	dim_typ * a3;
	dim_typ * a4;
	dim_typ * a5;
} _2dtpi3pdt;

typedef struct
{
	dim_typ a0;
	dim_typ * a1;
	ityp * a2;
	ityp * a3;
	dim_typ a4;
	ityp * a5;
	dim_typ a6;
	ityp * a7;
} dtpdt2pitdtpitdtpit; 

typedef struct
{
	ityp a0;
	ityp a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
	dim_typ a5;
	ityp * a6;
	ityp * a7;
	dim_typ a8;
} _2itdt2pitdt2pitdt;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	ityp * a6;
	ityp * a7;
} _6it2pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp a2;
} _2dtit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	ityp a3;
	ityp a4;
} _3dt2it;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	ityp a4;
	ityp a5;
	ityp a6; 
} _4dt3it;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	ityp * a2;
} _2pdtpit;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	ityp a4;
	int * a5;
	int * a6;
	ityp * a7;
} pit3dtit2pipit;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	ityp a4;
	int * a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
	int * a9;
	ityp * a10;
	int a11;
} pit3dtitpi3pitpipiti;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	ityp a4;
	int * a5;
	int * a6;
	ityp * a7;
	ityp * a8;
	ityp * a9;
	int * a10;
	ityp * a11;
	int a12;
} pit3dtit2pi3pitpipiti;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
	int * a8;
	ityp * a9;
} pit4dt3pitpipit;

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
} it5pit;

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	dim_typ * a6;
	ityp * a7;
} itpit4itpdtpit;

typedef struct
 {
 	ityp a0;
 	ityp * a1;
 	ityp * a2;
 	dim_typ a3;
} it2pitdt;

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
} it4pit;

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp a2;
	ityp a3; 
} itpit2it;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp * a3;
	ityp * a4;
	ityp * a5; 
} _3it3pit;

typedef struct
{
	char a0;
	ityp a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
} chitdt2pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	ityp * a3;
	int * a4;
	int * a5;
} _3dtpit2pi;

typedef struct
{
	ityp * a0;
	ityp * a1;
	ityp a2;
	dim_typ a3;
	int * a4;
} _2pititdtpi;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ * a2;
} dtpitpdt;

typedef struct
{
	ityp a0;
	ityp * a1;
	int * a2;
} itpitpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	ityp * a3;
	dim_typ * a4;
	int * a5;
	dim_typ a6;
	dim_typ a7;
	dim_typ a8;
	ityp * a9;
	dim_typ * a10;
	dim_typ * a11;
} _3dtpitpdtpi3dtpit2pdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	ityp * a4;
	dim_typ * a5;
	dim_typ * a6;
	dim_typ * a7;
	dim_typ * a8;
	dim_typ * a9;
	dim_typ * a10;
} _3dtpit6pdt;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp * a2;
} _2itpit; 

typedef struct
{
	ityp * a0;
	ityp a1;
	ityp a2;
	ityp a3;
	dim_typ a4;
	ityp * a5;
} pit3itdtpit;

typedef struct
{
	ityp * a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	dim_typ a6;
	ityp * a7;
} pit5itdtpit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
} _6it3pit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp * a4;
	ityp * a5;
} _4it2pit;

typedef struct
{
	ityp a0;
	dim_typ a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp * a5;
} itdt3itpit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	ityp * a4;
	int * a5;
	int * a6;
	int * a7;
} _4dtpit3pi;

typedef struct
{
	ityp * a0;
	ityp * a1;
	ityp * a2;
	bool * a3;
} _3pitpb;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
} _3it4pit; 

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp * a4;
} _4itpit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	ityp * a6;
} _6itpit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
} _4it3pit;

typedef struct
{
	ityp * a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	dim_typ * a4;
	ityp * a5;
} _4pitpdtpit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5; 
	dim_typ * a6;
	ityp * a7;
} _6itpdtpit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	ityp a6;
	ityp a7;
	ityp * a8;
	ityp * a9;
	ityp * a10;
} _8it3pit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
	ityp * a9;
	ityp * a10;	
} _6it5pit;

typedef struct
{
	ityp * a0;
	ityp * a1;
	ityp * a2;
	dim_typ a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
} _3pitdt3pit;

typedef struct
{
	ityp * a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	dim_typ a4;
	ityp * a5;
	ityp * a6;
	dim_typ * a7;
} _4pitdt2pitpdt;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	ityp a6;
	ityp a7;
	ityp a8;
	ityp a9;
	ityp * a10;
} _10itpit; 

typedef struct
{
	ityp * a0;
	ityp * a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
} _2pit4it3pit;

typedef struct
{ 
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp * a4;
	int * a5;
	ityp * a6;
} _4itpitpipit;

typedef struct
{
	ityp * a0;
	ityp * a1;
	ityp a2;
	ityp a3;
	dim_typ * a4;
	ityp * a5;
} _2pit2itpdtpit;

typedef struct
{
	int * a0;
	ityp * a1;
	ityp * a2;
} pi2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ * a2;
	int * a3;
} dtpitpdtpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
} _2dt3pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
} dtpitpi;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ * a3;
	dim_typ a4;
	dim_typ * a5;
} pit2dtpdtdtpdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	ityp * a3;
	dim_typ * a4;
	dim_typ * a5;
	ityp * a6;
} _3dtpit2pdtpit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	ityp a3;
} _2dtpitit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
} _2dt4pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp a2;
	dim_typ * a3;
	dim_typ * a4; 
} dtpitit2pdt;

typedef struct
{
	ityp a0;
	int * a1;
	int * a2;
	int * a3;
} it3pi;

typedef struct
{
	dim_typ a0;
	int * a1;
	dim_typ * a2;
	bool * a3;
} dtpipdtpb;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	bool * a3;
} dt2pipb;

typedef struct
{
	dim_typ a0;
	dim_typ * a1;
	dim_typ * a2;
	dim_typ * a3;
} dt3pdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ * a3;
	dim_typ * a4;
} _3dt2pdt;

typedef struct
{
	ityp a0;
	ityp * a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	dim_typ a5;
	dim_typ * a6;
	ityp * a7;
} itpit4dtpdtpit;

typedef struct
{
	ityp a0;
	ityp * a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	dim_typ * a5;
	ityp * a6;
} itpit3dtpdtpit;

typedef struct
{
	ityp a0;
	ityp * a1;
	dim_typ a2;
	ityp * a3;
	dim_typ a4;
	ityp * a5;
	ityp a6;
	ityp a7;
} itpitdtpitdtpit2it;

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;	
	ityp * a6;
} it6pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ * a3;
	dim_typ * a4;
	dim_typ * a5;
} dt2pit3pdt;

typedef struct
{
	int * a0;
	dim_typ * a1;
	bool * a2;
} pipdtpb;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	int * a2;
	ityp * a3;
} pitdtpipit;

typedef struct
{
	ityp * a0;
	ityp * a1;
	ityp a2;
	char a3;
} _2pititch; 

typedef struct
{
	ityp * a0;
	ityp * a1;
	char * a2;
	ityp a3;
} _2pitpchit; 

typedef struct
{
	ityp * a0;
	ityp * a1;
	ityp * a2;
	bool * a3;
	ityp * a4;
} _3pitpbpit; 

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ * a3;
	ityp * a4;
} pit2dtpdtpit;

typedef struct
{
	ityp * a0;
	ityp a1;
	ityp a2;
	ityp a3;
	dim_typ * a4;
	ityp * a5;
} pit3itpdtpit;

typedef struct
{
	ityp * a0;
	ityp * a1;
	bool * a2;
} _2pitpb;

typedef struct
{
	dim_typ a0;
	dim_typ * a1;
	dim_typ * a2;
	dim_typ * a3;
	dim_typ * a4;
} dt4pdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ * a4;
	dim_typ * a5;
	dim_typ * a6;
	dim_typ * a7;
} _4dt4pdt;

typedef struct
{
	dim_typ * a0;
	int * a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	int * a5;
} pdtpi3dtpi;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
} dt7pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	dim_typ a5; 
	ityp * a6;
	ityp * a7;
	ityp * a8;
} dt4pitdt3pit;

typedef struct
{
	dim_typ a0;
	dim_typ * a1;
} dtpdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	int * a3;
	dim_typ * a4;
} _3dtpipdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp a2;
	ityp a3;
	ityp * a4;
	ityp * a5;
} _2dt2it2pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
	ityp * a6;
	ityp * a7;
} _2dt4it2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	int a3;
	ityp * a4;
	int * a5;
	ityp * a6;
	ityp * a7;
	dim_typ a8;
	ityp a9;
	ityp a10;
	ityp a11;
	ityp a12;
} dtpitpiipitpi2pitdt4it;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp a3;
	ityp * a4;
	ityp * a5;	
} dt2pitit2pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	int * a3;
} _3dtpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ * a2;
	ityp * a3;
	ityp * a4;
	dim_typ * a5;
} _2dtpdt2pitpdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp a2;
	ityp * a3;
	int * a4;
} _2dtitpitpi;

typedef struct
{
	ityp a0;
	int * a1;
	int * a2;
	int * a3;
	int * a4;
} it4pi;

typedef struct
{
	int a0;
	int a1;
	int * a2; 
	int * a3;
} _2i2pi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	dim_typ * a3;
	dim_typ * a4;
	ityp * a5;
} _2dtpi2pdtpit;

typedef struct
{
	int a0;
	int * a1;
	int * a2;
} i2pi;

typedef struct
{
	ityp a0;
	int a1;
} iti;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ * a2;
} _2dtpdt;

typedef struct
{
	int *** a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
} pppi3dt;

typedef struct
{
	int ** a0;
	dim_typ a1;
	dim_typ a2;
} ppi2dt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	dim_typ a3;
	dim_typ a4;
} _2dtpi2dt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int * a3;
} _2dt2pi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int a3;
	dim_typ * a4;
	dim_typ * a5;
} _2dtpii2pdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int a3;
	int a4;
	dim_typ * a5;
	dim_typ * a6;
} _2dtpi2i2pdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	short * a3;
} _2dtpips;

typedef struct
{
	dim_typ a0;
	int a1;
	dim_typ * a2;
	dim_typ * a3;
} dti2pdt;

typedef struct
{
	int a0;
	int a1;
	int a2;
	dim_typ * a3;
	dim_typ * a4;
	dim_typ * a5;
} _3i3pdt;

typedef struct
{
	int a0;
	int a1;
	int a2;
	float a3;
	float a4;
} _3i2f;

typedef struct
{
	int a0;
	int a1;
	int a2;
	ityp a3;
	ityp a4;
} _3i2it;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	dim_typ * a3;
	dim_typ * a4;
} _2dtpi2pdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	int * a3;
	int * a4;
} _3dt2pi;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
} dt2pi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int a2;
	int a3;
	int * a4;
	int * a5;
} _2dt2i2pi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int a2;
	int a3;
	int * a4;
} _2dt2ipi;

typedef struct
{
	dim_typ a0;
	int * a1;
	dim_typ * a2;
} dtpipdt;

typedef struct
{
	dim_typ a0;
	int a1;
	int * a2;
	int a3;
	int * a4;
	int a5;
} dtipiipii;

typedef struct
{
	dim_typ a0;
	int * a1;
	int a2;
	short * a3;
	short * a4;
} dtpii2ps;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	dim_typ a3;
	dim_typ a4;
	int * a5;
} _2dtpi2dtpi;

typedef struct
{
	dim_typ * a0;
	int * a1;
} pdtpi;

typedef struct
{
	dim_typ * a0;
	int * a1;
	int a2;
} pdtpii;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int a3;
	int * a4;
	int * a5;
	int * a6;
} dt2pii3pi;

typedef struct
{
	dim_typ a0;
	int * a1;
	dim_typ * a2;
	int * a3;
	int * a4;
	short * a5;
} dtpipdt2pips;

typedef struct
{
	int * a0;
	int * a1;
	int * a2;
	int a3;
} _3pii;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int * a3;
	int * a4; 
} dt4pi;

typedef struct
{
	int * a0;
	int * a1;
	dim_typ * a2; 
} _2pipdt;

typedef struct
{
	dim_typ a0;
	int * a1; 
	int a2;
	int a3;
} dtpi2i;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int * a3;
	int * a4; 	
} _2dt3pi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int a3;
	int * a4; 	
} _2dtpiipi;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int * a3;
} dt3pi;

typedef struct
{
	dim_typ a0;
	int * a1;
	dim_typ a2;
	dim_typ * a3;
	int * a4;
	int * a5; 
} dtpidtpdt2pi;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int * a3;
	int a4;	
} dt3pii;

typedef struct
{
	int a0;
	dim_typ a1;
	bool * a2; 
} idtpb;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	dim_typ a3;
	dim_typ a4;
} _2dtpit2dt;

typedef struct
{
	dim_typ a0; 
	ityp * a1;
	ityp a2;
	int * a3;
	int * a4;
} dtpitit2pi;

typedef struct
{
	char * a0;
	dim_typ a1;
} pcdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp a2;
	int * a3; 
} _2dtitpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	ityp * a3;  
} _2dtpipit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int a2;
	ityp * a3;
	ityp * a4;
	int * a5;
	int * a6;
} dtpiti2pit2pi;

typedef struct
{
	dim_typ * a0;
	ityp a1;
	ityp a2;
} pdt2it;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	ityp a3;
} dt3it;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int * a3;
	int * a4;
	int * a5;
} _2dt4pi;

typedef struct
{
	ityp a0;
	ityp * a1;
	dim_typ * a2;
} itpitpdt;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	int * a3;
} dt2itpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	dim_typ a5;
	ityp * a6;
	ityp * a7;
} _2dt3pitdt2pit;

typedef struct
{
	long int * a0;
	long int * a1;
	long int * a2;
	long int * a3;
	long int * a4;
	long int * a5;
	long int * a6;
	long int * a7;
	long int * a8;
	ityp * a9;
	ityp * a10;
	ityp * a11;
	ityp * a12;
} _9pli4pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ a3;
	int * a4;
	int * a5;
	dim_typ a6;
	ityp * a7;
} dt2pitdt2pidtpit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	int * a3;
	int * a4;
} dtpit3pi;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int * a3;
	dim_typ a4;
	ityp * a5;
	dim_typ a6;
	int * a7;
	int * a8;
	int * a9;
} dt3pidtpitdt3pi;

typedef struct
{
	ityp a0;
	ityp a1;
	dim_typ a2;
	ityp * a3;
	dim_typ a4;
	int * a5;
	int * a6;
	int * a7;
	int * a8;
	int * a9;
	int * a10;
} _2itdtpitdt6pi;

typedef struct
{
	ityp * a0;
	ityp * a1;
	int * a2;
	ityp * a3;
	ityp * a4;
} _2pitpi2pit;

typedef struct
{
	ityp a0;
	ityp a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	dim_typ * a6;
} _2itdt3pitpdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	ityp * a4;
	ityp * a5;
} _4dt2pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
} _3dt3pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
} _5dt3pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	dim_typ a5;
	dim_typ a6;
	ityp * a7;
	ityp * a8;
	ityp * a9;
} _7dt3pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	dim_typ a5;
	dim_typ a6;
	dim_typ a7;
	dim_typ a8;
	ityp * a9;
	ityp * a10;
	ityp * a11;
} _9dt3pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int * a3;
	int * a4;
	bool * a5;
} _2dt3pipb;

typedef struct
{
	int a0;
	int a1;
	int a2;
	int * a3;
} _3ipi;

typedef struct
{
	dim_typ a0;
	int * a1;
	bool * a2;
} dtpipb;

typedef struct
{
	dim_typ a0;
	int * a1;
	bool * a2;
	dim_typ * a3;
	dim_typ * a4;
} dtpipb2pdt;

typedef struct
{
	int a0;
	int a1;
	int a2;
	int * a3;
	int * a4;
	int * a5;
} _3i3pi;

typedef struct
{
	dim_typ a0;
	dim_typ * a1;
	int * a2;
	int * a3;
	bool * a4;
} dtpdt2pipb;

typedef struct
{
	dim_typ a0;
	int * a1;
	dim_typ * a2;
	int * a3;
} dtpipdtpi;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int * a3;
	int * a4;
	dim_typ * a5;
} dt4pipdt;

typedef struct
{
	dim_typ a0;
	int a1;
	int a2;
	int a3;
	int a4;
	int * a5;
	int * a6;
} dt4i2pi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	int * a3;
	dim_typ * a4;
	dim_typ * a5;
	ityp * a6;
	int * a7;
} _2dtpitpi2pdtpitpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ * a2;
	int * a3;
	int * a4;
} _2dtpdt2pi;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int * a3;
	int * a4;
	bool * a5; 
} dt4pipb;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int * a3;
	ityp * a4;
	bool * a5; 
} dt3pipitpb;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	int * a3;
	ityp * a4;
	ityp a5;
	bool * a6;
} dtpit2pipititpb;

typedef struct
{
	dim_typ a0;
	int * a1;
	dim_typ a2;
	dim_typ a3;
	int * a4;
	bool * a5;
} dtpi2dtpipb;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	int * a3;
	int * a4;
	ityp a5;
	ityp a6;
	bool * a7;
} dtpit3pi2itpb;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	int * a3;
	ityp a4;
	ityp a5;
	bool * a6;
} dtpit2pi2itpb;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int * a3;
	bool * a4;
} dt3pipb;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp a2;
	ityp a3;
	int a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
} _2dt2iti3pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp a2;
	ityp a3;
	int a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
} _2dt2iti4pit;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp * a2;
	ityp * a3;
} dtit2pit;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	dim_typ a3;
	ityp * a4;
} dt2itdtpit;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	dim_typ a3;
	ityp * a4;
	dim_typ * a5;
	ityp * a6;
	ityp * a7;
} dt2itdtpitpdt2pit;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	dim_typ a3;
} dt2itdt;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	ityp * a3;
	dim_typ * a4;
	dim_typ a5;
} pit2dtpitpdtdt;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
	dim_typ a5;
	dim_typ a6;
	ityp * a7;
	ityp * a8;
	ityp * a9;
	ityp * a10;
} pit2dt2pit2dt4pit;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	ityp * a5;
	dim_typ a6;
	dim_typ a7;
	ityp * a8;
	ityp * a9;
	dim_typ a10;
} pit4dtpit2dt2pitdt;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	dim_typ * a5;
	ityp * a6;
} pit4dtpdtpit;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	int * a5;
	ityp * a6;
} pit4dtpipit;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	dim_typ * a5;
} pit4dtpdt;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	dim_typ * a5;
	ityp * a6;
	dim_typ a7;
} pit4dtpdtpitdt;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	int * a3;
	ityp * a4;
	ityp * a5;
	dim_typ a6; 
} pit2dtpi2pitdt;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ * a3;
} pit2dtpdt;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	ityp * a4;
	int * a5;
	ityp * a6;
	dim_typ a7;
} pit3dtpitpipitdt;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
	ityp * a9;
	ityp * a10;
	dim_typ a11;
} pit3dt7pitdt;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	dim_typ a2;
	dim_typ a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
	dim_typ a7;
	ityp * a8;
	dim_typ a9;
	ityp * a10;
	dim_typ a11;
} pit3dt3pitdtpitdtpitdt;

typedef struct
{
	dim_typ a0; 
	double complex * a1;
	dim_typ a2;
	double complex * a3;
} dtpcxdtpcx;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	double complex * a2;
	double complex * a3;
	double complex * a4;
} _2dt3pcx;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	dim_typ a2;
	double complex * a3;
	double complex * a4;
	double complex * a5;
} _3dt3pcx;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	ityp a2;
	double complex * a3;
	ityp a4;
	double complex * a5;
	double complex * a6;
} _2dtitpcxit2pcx;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	ityp a2;
	double complex * a3;
} _2dtitpcx;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	double complex * a2;
} _2dtpcx;

typedef struct
{
	dim_typ a0; 
	double complex * a1;
} dtpcx;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	ityp * a3;
	ityp a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
} _2dt2pitit3pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	ityp * a2;
	int * a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
} _2dtpitpi3pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	ityp a6;
	ityp a7;
} _2dt4pit2it;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	ityp ** a2;
	ityp ** a3;
	ityp ** a4;
} _2dt3ppit;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	int * a2;
	dim_typ * a3;
} _2dtpipdt;

typedef struct
{
	dim_typ a0; 
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp (* a5)(ityp);
	ityp (* a6)(ityp);
	dim_typ * a7;
	ityp * a8;
} dt4it2fitpdtpit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	ityp a3;	
	dim_typ a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
	ityp * a9;
} _4itdt5pit;

typedef struct
{
	ityp a0;
	ityp a1;	
	dim_typ a2;
	ityp * a3;
	ityp * a4;
	ityp a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
} _2itdt2pitit3pit;

typedef struct
{
	ityp a0;	
	dim_typ a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
	ityp * a9;
} itdt8pit;

typedef struct
{
	ityp a0;
	ityp a1;	
	dim_typ a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
	ityp * a9;
	ityp * a10;
} _2itdt8pit;

typedef struct
{
	ityp a0;
	ityp * a1;
	ityp a2;
	ityp * a3;
	ityp * a4;
	dim_typ a5;
	dim_typ a6;
	ityp * a7;
	ityp * a8;
} itpitit2pit2dt2pit;

typedef struct
{
	ityp a0;
	ityp a1;	
	dim_typ a2;
	ityp a3;
	dim_typ a4;
	ityp a5;
	ityp a6;
	ityp * a7;
	ityp * a8;
} _2itdtitdt2it2pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	ityp a2;
	ityp a3;
	ityp ** a4;
	dim_typ a5;
	dim_typ a6;
	ityp ** a7;
	ityp ** a8;
} _2dt2itppit2dt2ppit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	ityp * a3;
} _3dtpit;

typedef struct
{
	bool a0;
	dim_typ a1;
	int * a2;
} bdtpi;

typedef struct
{
	dim_typ a0;
	short a1;
	ityp a2;
	ityp a3;
	ityp * a4;
	ityp * a5;
} dts2it2pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1; 
	dim_typ a2;
	int * a3;
	int * a4;
	ityp * a5;
	ityp * a6;
} _3dt2pi2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ a3;
	ityp a4;
	ityp * a5;
	dim_typ * a6;
} dt2pitdtitpitpdt;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2;
	ityp a3;
	ityp a4;
} dt4it;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	bool * a2;	
} dtpitpb;

typedef struct
{
	ityp * a0;
	int * a1;
	ityp * a2; 
	int * a3;
} pitpipitpi;

typedef struct
{
	ityp a0;
	int * a1;
	ityp * a2; 
	int * a3;
} itpipitpi;

typedef struct
{
	ityp a0;
	ityp a1;	
	int * a2;
} _2itpi;

typedef struct
{
	ityp a0;
	int a1;	
	int * a2;
} itipi;

typedef struct
{
	int a0;
	dim_typ a1;
	ityp a2;
} idtit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	int a3;
	int a4;
} _3it2i;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	int a3;
} _3iti;

typedef struct
{
	int a0;
	ityp a1;
	ityp a2;
} i2it;

typedef struct
{
	dim_typ a0;
	int * a1;
	int a2;
	ityp * a3;
} dtpiipit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	dim_typ * a3;
	dim_typ * a4;
} _2dtpit2pdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int a3;
	ityp * a4;
} _2dtpiipit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int a2;
	ityp * a3;
} _2dtipit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	dim_typ a3;
	ityp a4;
	int * a5;
	int * a6;
} _2dtpitdtit2pi;

typedef struct
{
	dim_typ a0;
	int a1;
	dim_typ * a2;
	dim_typ a3;
	int * a4;
	int * a5;
} dtipdtdt2pi;

typedef struct
{
	ityp ** a0;
	dim_typ a1;
	dim_typ a2;
} ppit2dt;

typedef struct
{
	dim_typ a0;
	int * a1;
	dim_typ a2;
	ityp * a3;
	int * a4;
	ityp * a5;
} dtpidtpitpipit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	ityp a3;
	dim_typ * a4;
	ityp * a5;
	int * a6;
} dtpitpiitpdtpitpi;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	dim_typ * a3;
	ityp * a4;
	int * a5;
} dtpitpipdtpitpi;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	int * a2;
	ityp a3;
} pdtpitpiit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	ityp a3;
	dim_typ * a4;
	dim_typ * a5;
	dim_typ * a6;
} dtpitpiit3pdt;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	ityp a3;
	ityp a4;
	int * a5;
	int * a6;
} dtpitpi2it2pi;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	int * a2;
} pdtpitpi;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ * a2;
	dim_typ * a3;
} dtpit2pdt;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp a2;
	ityp a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
} dtpit2it3pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ a2;
	ityp a3;
	int * a4;
	int * a5;
} dtpitdtit2pi;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp a2;
	dim_typ * a3;
} dtpititpdt;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp a2;
	dim_typ a3;
	dim_typ * a4;
	ityp * a5;
	int * a6;
} dtpititdtpdtpitpi;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ * a3;
} dt2pitpdt;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ * a3;
	int * a4;
} dt2pitpdtpi;

typedef struct
{
	dim_typ a0;
	bool * a1;
} dtpb;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	int * a3;
	ityp a4;
	ityp a5;
} dtpit2pi2it;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp a3;
	dim_typ a4;
	ityp * a5; 
} dt2pititdtpit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp a4;
	dim_typ a5;
	ityp * a6; 
	ityp * a7; 
} dt3pititdt2pit;

typedef struct
{
	dim_typ a0; 
	dim_typ a1;
	ityp * a2;
	ityp * a3;
	ityp a4;
	dim_typ a5;
	ityp * a6;
} _2dt2pititdtpit;

typedef struct
{
	int a0; 
	dim_typ a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
} i2dt2pit;

typedef struct
{
	int a0;
	dim_typ a1;
} idt;

typedef struct
{
	dim_typ a0;
	int * a1;
	dim_typ a2;
	ityp * a3;
	ityp * a4;
} dtpidt2pit;

typedef struct
{
	dim_typ a0; 
	dim_typ a1;
	dim_typ a2; 
	int * a3;
	ityp * a4;
} _3dtpipit;

typedef struct
{
	dim_typ a0; 
	ityp a1;
	ityp a2; 
	ityp (* a3)(ityp);
} dt2itfit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ a2;
	int * a3;
	int a4;
	ityp * a5;
} dtpitdtpiipit;

typedef struct
{
	ityp a0;
	ityp a1; 
	ityp a2;
	ityp a3;
	dim_typ a4;
	ityp * a5;
	ityp * a6;
} _4itdt2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	dim_typ a3;
	int * a4;
	int * a5;
	dim_typ a6;
	ityp * a7;
	ityp * a8;
} dt2pitdt2pidt2pit;

typedef struct
{
	dim_typ a0; 
	dim_typ a1;
	dim_typ a2; 
	dim_typ a3;
	dim_typ a4; 
	dim_typ a5;
	dim_typ a6; 
	dim_typ * a7;
	ityp * a8;
} _7dtpdtpit;

typedef struct
{
	ityp a0;
	ityp * a1;
	dim_typ a2; 
	dim_typ a3;
	dim_typ a4; 
} itpit3dt;

typedef struct
{
	dim_typ a0;
	ityp a1;
	ityp a2; 
	ityp a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
} dt3it4pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2; 
	ityp * a3;
	int * a4;
} dt3pitpi;

typedef struct
{
	dim_typ a0;
	int * a1;
	ityp * a2; 
	ityp * a3;
	int * a4;
} dtpi2pitpi;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2; 
	dim_typ a3;
	ityp * a4;
	dim_typ a5;
	dim_typ a6;
	ityp * a7; 
	ityp * a8;
} _3itdtpit2dt2pit;

typedef struct
{
	ityp a0;
	ityp a1; 
	dim_typ a2;
	ityp * a3;
	dim_typ a4;
	ityp * a5; 
	ityp * a6;
} _2itdtpitdt2pit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2; 
	ityp a3;
	ityp (* a4)(ityp);
	ityp (* a5)(ityp);
	int * a6;
} _4it2fitpi;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2; 
	ityp * a3;
	ityp * a4;
	ityp * a5;
} dt5pit;

typedef struct
{
	ityp * a0;
	dim_typ a1;
	ityp a2;
	ityp a3;
	ityp * a4;
	ityp * a5;
	ityp * a6; 
	ityp * a7;
} pitdt2it4pit;

typedef struct
{
	dim_typ a0;
	int * a1;
	short * a2;
	int * a3;
	bool a4;
} dtpipspib;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	int * a3;
	dim_typ a4;
	dim_typ a5;
	int * a6;
	ityp * a7;
	int * a8;
} _2dtpitpi2dtpipitpi;

typedef struct
{
	dim_typ a0;
	int * a1;
	dim_typ a2;
	ityp * a3;
	int * a4;
	int * a5;
	dim_typ a6;
	ityp * a7;
} dtpidtpit2pidtpit;

typedef struct
{
	dim_typ a0;
	int a1;
	ityp * a2;
	ityp * a3;
} dti2pit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp a2;
	int * a3;
} _3itpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp a2;
	ityp a3;
	ityp a4;
	ityp a5;
} _2dt4it;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	int * a2;
	dim_typ a3;
	int * a4;
	dim_typ a5;
	dim_typ a6;
	ityp * a7;
	ityp * a8;
	dim_typ a9;
	ityp a10;
} dtpitpidtpi2dt2pitdtit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp * a2;
	ityp * a3;
	dim_typ a4;
	ityp * a5;
	ityp * a6;
} _2dt2pitdt2pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	ityp a2;
	ityp a3;
	dim_typ a4;
	ityp * a5;
	ityp * a6;
} _2dt2itdt2pit;

typedef struct
{
	ityp * a0;
	ityp * a1;
	int * a2;
	dim_typ a3;
	dim_typ a4;
	int * a5;
	dim_typ a6;
	dim_typ a7;
	dim_typ a8; 
	dim_typ a9;
	ityp a10;
	ityp a11;
	dim_typ a12;
	ityp a13;
	ityp a14;
} _2pitpi2dtpi4dt2itdt2it;

typedef struct
{
	ityp * a0;
	ityp * a1;
	int * a2;
	dim_typ a3;
	dim_typ a4;
	int * a5;
	dim_typ a6;
	ityp * a7;
	dim_typ a8;
	dim_typ a9; 
	dim_typ a10;
	ityp a11;
	ityp a12;
	dim_typ a13;
} _2pitpi2dtpidtpit3dt2itdt;

typedef struct
{
	dim_typ a0; 
	ityp * a1;
	int * a2;
	dim_typ a3;
	int * a4;
	dim_typ a5;
	dim_typ a6;
	ityp a7;
	ityp a8;
	ityp a9;
	ityp a10;
	dim_typ a11;
} dtpitpidtpi2dt4itdt;

typedef struct
{
	dim_typ a0; 
	ityp * a1;
	int * a2;
	dim_typ a3;
	int * a4;
	dim_typ a5;
	dim_typ a6;
	ityp a7;
	ityp a8;
	dim_typ a9;
} dtpitpidtpi2dt2itdt;

typedef struct
{
	dim_typ a0; 
	ityp * a1;
	dim_typ a2;
	ityp a3;
	ityp a4;
	ityp a5;
	ityp a6;
	dim_typ a7;
} dtpitdt4itdt;

typedef struct
{
	dim_typ a0; 
	dim_typ a1;
	dim_typ a2;
	ityp a3;
	ityp * a4;
	dim_typ a5;
	ityp * a6; 
	ityp * a7;
	ityp * a8;
} _3dtitpitdt3pit;

typedef struct
{
	char a0;
	char a1;
	dim_typ a2;
	dim_typ a3;
	dim_typ a4;
	ityp a5;
	ityp * a6; 
	dim_typ a7;
	ityp * a8;
	dim_typ a9;
	ityp a10;
	ityp * a11;
	dim_typ a12;
} _2c3dtitpitdtpitdtitpitdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	ityp * a3; 
	ityp * a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
} _3dt6pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	dim_typ a2;
	int * a3; 
	int * a4;
	int * a5;
	int * a6;
	int * a7;
	int * a8;
} _3dt6pi;

typedef struct
{
	dim_typ a0;
	int a1;
	int * a2;
	bool * a3;
} dtipipb;

typedef struct
{
	dim_typ a0;
	int a1;
	int * a2;
	int * a3;
} dti2pi;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int * a3;
	int * a4; 
	int * a5;
} dt5pi;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	int * a3;
	int * a4; 
	int * a5;
	int * a6;
} dt6pi;

typedef struct
{
	dim_typ a0;
	int * a1;
	int * a2;
	dim_typ * a3;
	int * a4; 
} dt2pipdtpi;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	dim_typ a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
} dt3pitdt3pit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp * a2;
	ityp * a3;
} _2it2pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	dim_typ a6;
} dt5pitdt;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ a2;
	dim_typ a3;
	int * a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
} dtpit2dtpi3pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ a2;
	dim_typ a3;
	int * a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
} dtpit2dtpi4pit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ a2;
	dim_typ a3;
	int * a4;
	ityp * a5;
	ityp * a6;
	ityp * a7;
	ityp * a8;
	ityp * a9;
} dtpit2dtpi5pit;

typedef struct
{
	int a0;
	int a1;
	int a2;
	int * a3; 
	int * a4;
	int * a5;
	int * a6;
	int * a7;
	dim_typ a8;
	int * a9;
	int * a10;	
} _3i5pidt2pi;

typedef struct
{
	int * a0;
	ityp * a1;
} pipit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ a2;
	dim_typ a3;
	int * a4;
	ityp * a5;
} dtpit2dtpipit;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ a2;
	dim_typ a3;
	int * a4;
	int * a5;
	ityp * a6;
	int * a7;
	ityp * a8;
	int * a9;
} dtpit2dt2pipitpipitpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int a3;
} _2dtpii;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int * a3;
	dim_typ a4;
	int * a5;
} _2dt2pidtpi;

typedef struct
{
	dim_typ a0;
	dim_typ a1;
	int * a2;
	int * a3;
	dim_typ a4;
	int * a5;
	int * a6;
	int * a7;
} _2dt2pidt3pi;

typedef struct
{
	dim_typ a0;
	ityp * a1;
	dim_typ a2;
	dim_typ a3;
	int * a4;
	int * a5;
	ityp * a6;
	int * a7;
	ityp * a8;
	ityp * a9;
	ityp * a10;
	int * a11;
	int * a12;
} dtpit2dt2pipitpi3pit2pi;

typedef struct
{
	dim_typ * a0;
	int * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
} pdtpi3pit;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	unsigned * a2;
} _2pdtpu;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	dim_typ * a2;
	unsigned * a3;
} _3pdtpu;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
} _2pdt4pit;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	int * a2;
} _2pdtpi;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	short * a2;
} _2pdtps;

typedef struct
{
	dim_typ * a0;
	dim_typ * a1;
	short * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
} _2pdtps4pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;	
	int * a2;
	int * a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
} _2dt2pi3pit;

typedef struct
{
	dim_typ a0;
	dim_typ a1;	
	dim_typ * a2;
	dim_typ * a3;
	dim_typ * a4;
} _2dt3pdt;

typedef struct
{
	dim_typ a0;
	dim_typ a1;	
	dim_typ a2;
	int * a3;
	dim_typ * a4;
	dim_typ * a5;
} _3dtpi2pdt;

typedef struct
{
	short a0;
	ityp a1; 
} sit;

typedef struct
{
	dim_typ a0;
	short a1;
	ityp * a2;
	ityp * a3;
} dts2pit;

typedef struct
{
	ityp a0;
	ityp a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
} _2it3pit;

typedef struct
{
	dim_typ * a0;
	ityp * a1;
	ityp * a2;
	ityp * a3;
	ityp * a4;
	ityp * a5;
	ityp * a6;
} pdt6pit;

__MATHSUITE __JBURKARDT  void RE_PINT( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_UINT( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PSHRT( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PUSHRT( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PCHR( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PUCHR( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PLNG( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PULNG( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PLLNG( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PULLNG( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PFLT( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PDBL( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_PLDBL( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_INT( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_UINT( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_SHRT( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_USHRT( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_CHR( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_UCHR( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_LNG( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_ULNG( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_LLNG( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_ULLNG( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_FLT( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_DBL( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT  void RE_LDBL( void *, EXPRTYPE [], exprNode *, mpfr_t);
__MATHSUITE __JBURKARDT void RE_CPLX( void *, EXPRTYPE [], exprNode *, mpfr_t);

__MATHSUITE __JBURKARDT void * H_CPLX (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT3PCX (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DT3PCX (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);

__MATHSUITE __JBURKARDT void * H_PINT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PINT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PUINT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PSHRT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PUSHRT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PCHR (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PUCHR (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PLNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PULNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PLLNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PULLNG (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PFLT (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PDBL (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PLDBL (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);

__MATHSUITE __JBURKARDT  uint64_t RC_MUL( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_ADD( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_ADDADDCONST( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_3MUL( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_4MUL( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_TRIANGLENUMADD( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_MULTRIANGLENUMADD( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_SCALADDCONSTMUL( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_ADDCONSTSCAL( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_SCAL( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_PPMULSCAL( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_PPMULDIV( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_PPPMUL2MULDIV( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_2ADD2CONSTDIV( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_PPMULADD( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_ADDCONST( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_DEREF( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_DEREFSCAL( EXPRTYPE args[], sel_typ opr[]);
__MATHSUITE __JBURKARDT  uint64_t RC_WAVELET( EXPRTYPE args[], sel_typ opr[]);;

__MATHSUITE __JBURKARDT void * H_PINT2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PINT3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PINT4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PINT5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PINT6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PINT7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PINT8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PINT9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PINT10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPINT2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPINT3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPINT4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPINT5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPINT6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPINT7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPINT8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPINT9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPINT10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUINT2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUINT3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUINT4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUINT5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUINT6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUINT7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUINT8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUINT9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUINT10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUINT2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUINT3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUINT4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUINT5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUINT6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUINT7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUINT8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUINT9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUINT10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PSHRT2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PSHRT3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PSHRT4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PSHRT5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PSHRT6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PSHRT7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PSHRT8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PSHRT9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PSHRT10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPSHRT2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPSHRT3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPSHRT4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPSHRT5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPSHRT6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPSHRT7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPSHRT8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPSHRT9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPSHRT10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUSHRT2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUSHRT3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUSHRT4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUSHRT5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUSHRT6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUSHRT7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUSHRT8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUSHRT9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUSHRT10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUSHRT2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUSHRT3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUSHRT4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUSHRT5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUSHRT6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUSHRT7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUSHRT8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUSHRT9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUSHRT10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PCHR2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PCHR3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PCHR4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PCHR5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PCHR6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PCHR7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PCHR8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PCHR9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PCHR10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPCHR2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPCHR3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPCHR4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPCHR5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPCHR6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPCHR7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPCHR8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPCHR9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPCHR10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUCHR2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUCHR3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUCHR4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUCHR5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUCHR6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUCHR7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUCHR8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUCHR9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PUCHR10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUCHR2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUCHR3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUCHR4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUCHR5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUCHR6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUCHR7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUCHR8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUCHR9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPUCHR10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLNG2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLNG3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLNG4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLNG5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLNG6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLNG7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLNG8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLNG9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLNG10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLNG2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLNG3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLNG4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLNG5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLNG6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLNG7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLNG8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLNG9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLNG10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULNG2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULNG3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULNG4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULNG5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULNG6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULNG7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULNG8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULNG9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULNG10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULNG2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULNG3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULNG4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULNG5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULNG6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULNG7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULNG8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULNG9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULNG10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLLNG2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLLNG3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLLNG4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLLNG5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLLNG6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLLNG7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLLNG8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLLNG9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLLNG10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLLNG2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLLNG3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLLNG4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLLNG5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLLNG6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLLNG7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLLNG8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLLNG9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLLNG10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULLNG2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULLNG3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULLNG4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULLNG5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULLNG6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULLNG7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULLNG8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULLNG9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PULLNG10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULLNG2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULLNG3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULLNG4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULLNG5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULLNG6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULLNG7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULLNG8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULLNG9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPULLNG10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PFLT2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PFLT3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PFLT4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PFLT5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PFLT6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PFLT7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PFLT8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PFLT9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PFLT10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPFLT2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPFLT3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPFLT4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPFLT5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPFLT6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPFLT7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPFLT8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPFLT9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPFLT10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDBL2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDBL3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDBL4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDBL5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDBL6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDBL7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDBL8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDBL9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDBL10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPDBL2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPDBL3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPDBL4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPDBL5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPDBL6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPDBL7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPDBL8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPDBL9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPDBL10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDBL12 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PDBL12 (EXPRTYPE val, exprObj *obj, exprNode *nodes, EXPRERRTYPE *err, EXPRTYPE args[]);
__MATHSUITE __JBURKARDT void * H_PLDBL2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLDBL3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLDBL4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLDBL5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLDBL6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLDBL7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLDBL8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLDBL9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PLDBL10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLDBL2 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLDBL3 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLDBL4 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLDBL5 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLDBL6 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLDBL7 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLDBL8 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLDBL9 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPLDBL10 (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);

__MATHSUITE __JBURKARDT void * H_SINT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SUINT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SSHRT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SUSHRT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SCHR (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SUCHR (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SLNG (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SULNG (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SLLNG (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SULLNG (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SFLT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SDBL (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SLDBL (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT2IPITPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PITDTPDTPITPDTDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDT4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDT5PITPDT2DT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4PITPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4PDTPITPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_B3DT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_5PDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2BDTPITPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PIT2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PIT2PDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2I2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3ITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3PDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTITPITDTPITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITIPITI2IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTITPITI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITIPITI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_7ITFITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4PITFITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4PITFIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2ITPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITITDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTIT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITDT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2IT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDT6PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDT6PITPSPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDT4PITPSPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDT5PITPSPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PS2PIT2PUL (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IPS2PIT2PUL7PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IPS4PIT2PUL4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PS4PIT2PUL (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDT2PITPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3PDT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2IPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTFITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTFIT2ITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTIT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITIT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPITDT4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT2PITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT3PITDTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPIT2ITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT3PITDTPI2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITDTPI2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2ITPIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPIPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPII (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT6PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2IT2PITPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2ITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_FIT2ITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4DT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT3PIPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_FDTDT2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT2DT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITDTPITPS (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITDTPS (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITPDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PIT2ITPDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITI2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PITDTPIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2PITITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4FITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3FITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITFIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4IT2FITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITDTFIT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPI3PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITDT2PITDT2PITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_6IT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DT2IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4DT3IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT3DTIT2PIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT3DTITPI3PITPIPITI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT3DTIT2PI3PITPIPITI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT4DT3PITPIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT5PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3IT4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPIT4ITPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT3ITDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPIT2IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3IT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PIT2ITPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_CHITDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DTPIT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PITITDTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT2DTPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DTPITPDTPI3DTPIT2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DTPIT6PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT5ITDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_6IT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4IT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITDT3ITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4DTPIT3PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3PITPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4ITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_6ITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4IT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_6ITPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_8IT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_6IT5PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3PITDT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT5PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4PITDT2PITPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_10ITPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PIT4IT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4ITPITPIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PI2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPDTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT2DTPDTDTPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DTPIT2PDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPITIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITIT2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT3PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIPDTPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PIPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DT2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPIT4DTPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPIT3DTPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPITDTPITDTPIT2IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT6PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PIT3PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4ITDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIPDTPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PITITCH (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PITPCHIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3PITPBPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT3ITPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PITPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4DT4PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDTPI3DTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT7PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4PITDT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DTPIPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPI2PDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT2DTPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPIPDTPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT2DTPITPDTDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2IT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT4IT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPIIPITPI2PITDT4IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPDT2PITPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTITPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IT4PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2I2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_I2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPPI3DT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPI2DT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPI2DT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPII2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPI2I2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPIPS (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTI2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3I3PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3I2F (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3I2IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPI2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2I2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2IPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTIPIIPII (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPII2PS (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT3PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPI2DTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDTPII (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTI2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PII3PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIPDT2PIPS (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3PII (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PIPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIPDTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPI2I (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPIIPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIDTPDT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3PII (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IDTPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPIT2DT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITIT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PCDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITI2PIT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDT2IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT4PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPITPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2ITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT3PITDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_9PLI4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITDT2PIDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT3PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3PIDTPITDT3PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITDTPITDT6PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PITPI2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITDT3PITPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_5DT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_7DT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_9DT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3IPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIPB2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3I3PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPDT2PIPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4PIPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4I2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPITPI2PDTPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPDT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4PIPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3PIPITPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2PIPITITPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPI2DTPIPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT3PI2ITPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2PI2ITPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3PIPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITDTPIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2ITI3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2ITI4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPDT2PITDTPITDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2ITDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2ITDTPITPDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2ITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT2DT2PIT2DT4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT4DTPIT2DT2PITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT4DTPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT4DTPIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT4DTPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT4DTPDTPITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT2DTPI2PITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT3DTPITPIPITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT3DT7PITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2PI3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIT3DT3PITDTPITDTPITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPCXDTPCX (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITDTPITDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTITPCXIT2PCX (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTITPCX (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPCX (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPCX (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2PITIT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPITPI3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT4PIT2IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT3PPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPIPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4IT2FITPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4ITDT5PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITDT2PITIT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITDT8PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITDT8PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPITIT2PIT2DT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITDTITDT2IT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2ITPPIT2DT2PPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_BDTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTS2IT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DT2PI2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2DTPIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITDTITPITPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DT6PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT4IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2PIDT3PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PITPIPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPIPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2ITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITIPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IDTIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3IT2I (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3ITI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_I2IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPIT2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPIIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPITDTIT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTIPDTDT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PPIT2DT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIDTPITPIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPIITPDTPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDTPITPIIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPIIT3PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITPDTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPI2IT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDTPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2IT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2PI2IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITDTIT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITITPDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITITDTPDTPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3PITITDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2PITITDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_I2DT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_IDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DTPIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2ITFIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITDTPIIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2DTPI3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PITDT2PIDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_7DTPDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_ITPIT3DT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3IT4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3PITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPI2PITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3ITDTPIT2DT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_4IT2FITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PITDT2IT4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIPSPIB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2PIDTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPITPI2DTPIPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIDTPIT2PIDTPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTI2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3ITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT4IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPIDTPI2DT2PITDTIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2PITDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT2ITDT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PITPI2DTPI4DT2ITDT2IT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PITPI2DTPIDTPIT3DT2ITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPIDTPI2DT4ITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITPIDTPI2DT2ITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPITDT4ITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DTITPITDT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2C3DTITPITDTPITDTITPITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DT6PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTIPIPB (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT5PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT6PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT2PIPDTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT3PITDT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2IT2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DT5PITDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2DTPI4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2DTPI5PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3I5PIDT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PIPIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2DT2PIPITPIPITPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DTPII (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTPIT2DT2PIPITPI3PIT2PI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_PDTPI3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PDTPU (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3PDTPU (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PDT4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PDTPI (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PDTPS (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2PDTPS4PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2DT3PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_3DTPI2PDT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_SIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_DTS2PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);
__MATHSUITE __JBURKARDT void * H_2IT3PIT (EXPRTYPE, exprObj *, exprNode *, EXPRERRTYPE *, EXPRTYPE []);

