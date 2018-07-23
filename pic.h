#include<pari/pari.h>

#define lgJ 12

long Jgetg(GEN J);
long Jgetd0(GEN J);
GEN JgetT(GEN J);
GEN Jgetp(GEN J);
long Jgete(GEN J);
GEN Jgetpe(GEN J);
GEN JgetFrob(GEN J);
GEN JgetV(GEN J);
GEN JgetKV(GEN J);
GEN JgetW0(GEN J);
GEN JgetZ(GEN J);
GEN JgetFrobCyc(GEN J);
void JgetTpe(GEN J, GEN* T, GEN* pe, GEN* p, long* e);

GEN PicRed(GEN J, ulong e);

GEN DivAdd(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess);
GEN DivSub(GEN WA, GEN WB, GEN KV, ulong d, GEN T, GEN p, long e, GEN pe, ulong nIGS);
GEN PicChord(GEN J, GEN WA, GEN WB, long flag);
GEN PicAdd(GEN J, GEN WA, GEN WB);
GEN PicSub(GEN J, GEN WA, GEN WB);
GEN PicNeg(GEN J, GEN W);
GEN PicMul(GEN J, GEN W, GEN n, long flag);
GEN PicFrob(GEN J, GEN W);
GEN PicFrobPoly(GEN J, GEN W, GEN F);
long PicEq(GEN J, GEN WA, GEN WB);
long PicIsZero(GEN J, GEN W);
GEN PicRand0(GEN J);

GEN PicChart(GEN J, GEN W, ulong P0);
