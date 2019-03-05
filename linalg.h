#include<pari/pari.h>

GEN GetFq1(GEN);
GEN Z2Fq(GEN n, GEN T);
long ZX_is0mod(GEN,GEN);
GEN FpXM_red(GEN,GEN);
GEN FpXM_add(GEN,GEN,GEN);
GEN FpXM_sub(GEN,GEN,GEN);
GEN FqV_Fq_mul(GEN,GEN,GEN,GEN);
GEN FqM_Fq_mul(GEN,GEN,GEN,GEN);
GEN ZXM_Z_mul(GEN,GEN);
GEN RandVec_1(GEN A,GEN pe);
GEN RandVec_padic(GEN,GEN,GEN,GEN);
GEN matkerpadic(GEN,GEN,GEN,long);
GEN mateqnpadic(GEN,GEN,GEN,long);
GEN matF(GEN,GEN,GEN,long);
GEN mat2col(GEN);
GEN col2mat(GEN,long,long);
GEN M2ABCD(GEN,GEN);
GEN M2ABCD_1block(GEN,ulong,ulong,GEN);
GEN VecSmallCompl(GEN,ulong);
GEN FqM_MinorCompl(GEN,GEN,GEN);
