#include<pari/pari.h>

#define lgJ 13

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

GEN DivAdd(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess);
GEN DivSub(GEN WA, GEN WB, GEN KV, ulong d, GEN T, GEN p, long e, GEN pe, ulong nIGS);
GEN PicChord(GEN J, GEN WA, GEN WB, long flag);
GEN PicAdd(GEN J, GEN WA, GEN WB);
GEN PicSub(GEN J, GEN WA, GEN WB);
GEN PicNeg(GEN J, GEN W);
GEN PicMul(GEN J, GEN W, GEN n, long flag);
