#include<pari/pari.h>

long ZX_is0mod(GEN,GEN);
GEN FpXM_red(GEN,GEN);
GEN FpXM_add(GEN,GEN,GEN);
GEN FpXM_sub(GEN,GEN,GEN);
GEN ZXM_add(GEN,GEN);
GEN ZXM_sub(GEN,GEN);
GEN FqV_Fq_mul(GEN,GEN,GEN,GEN);
GEN FqM_Fq_mul(GEN,GEN,GEN,GEN);
GEN ZXM_Z_mul(GEN,GEN);
GEN RandVec_1(GEN A,GEN pe);
GEN RandVec_padic(GEN,GEN,GEN,GEN);
GEN ZpXQMinv(GEN A, GEN T, GEN pe, GEN p, long e);
GEN matkerpadic(GEN A, GEN T, GEN pe, GEN p, long e);
GEN mateqnpadic(GEN A, GEN T, GEN pe, GEN p, long e);
GEN matF(GEN B,GEN T,GEN pe, GEN p, long e);
GEN mat2col(GEN);
GEN col2mat(GEN,ulong,ulong);
GEN M2ABCD(GEN,GEN);
GEN M2ABCD_1block(GEN,ulong,ulong,GEN);
GEN VecSmallCompl(GEN,ulong);
GEN FqM_MinorCompl(GEN,GEN,GEN);
GEN Subspace_normalize(GEN V, GEN I, GEN T, GEN pe, GEN p, long e, long drop);
