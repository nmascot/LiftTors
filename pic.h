#include<pari/pari.h>

#define lgJ 16
/* See pic.c for contents of a Jacobian */

GEN Jgetf(GEN J);
long Jgetg(GEN J);
long Jgetd0(GEN J);
GEN JgetL(GEN J);
GEN JgetT(GEN J);
GEN Jgetp(GEN J);
long Jgete(GEN J);
GEN Jgetpe(GEN J);
GEN JgetFrob(GEN J);
GEN JgetV(GEN J, ulong n);
GEN JgetKV(GEN J, ulong n);
GEN JgetW0(GEN J);
GEN JgetEvalData(GEN J);
GEN JgetZ(GEN J);
GEN JgetFrobCyc(GEN J);
GEN JgetAutsCyc(GEN J);
void JgetTpe(GEN J, GEN* T, GEN* pe, GEN* p, long* e);

GEN PicRed(GEN J, ulong e);

GEN ZpXQ_FrobMat(GEN T, GEN p, long e, GEN pe);
GEN Frob(GEN x, GEN FrobMat, GEN T, GEN pe);

GEN DivMul(GEN f, GEN W, GEN T, GEN pe);
GEN DivAdd(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess);
GEN DivSub(GEN WA, GEN WB, GEN KV, ulong d, GEN T, GEN p, long e, GEN pe, ulong nIGS);
ulong PicMember_val(GEN J, GEN W);
GEN PicNeg(GEN J, GEN W, long flag);
GEN PicChord(GEN J, GEN WA, GEN WB, long flag);
GEN PicAdd(GEN J, GEN WA, GEN WB);
GEN PicSub(GEN J, GEN WA, GEN WB);
GEN PicMul(GEN J, GEN W, GEN n, long flag);
GEN PicFrob(GEN J, GEN W);
GEN PicFrobInv(GEN J, GEN W);
GEN PicFrobPow(GEN J, GEN W, long n);
GEN PicFrobPoly(GEN J, GEN W, GEN F);
GEN PicAut(GEN J, GEN W, ulong n);
ulong PicEq_val(GEN J, GEN WA, GEN WB);
ulong PicIsZero_val(GEN J, GEN W);
GEN PicRand0(GEN J, GEN randseed);

GEN PicChart(GEN J, GEN W, ulong P0, GEN P1);
