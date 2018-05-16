#include<pari/pari.h>

#define lgJ 13

ulong Jgetg(GEN J);
void Jsetg(GEN J, ulong g);
ulong Jgetd0(GEN J);
void Jsetd0(GEN J, ulong d0);
GEN JgetT(GEN J);
GEN Jgetp(GEN J);
long Jgete(GEN J);
void Jsete(GEN J, ulong e);
GEN Jgetpe(GEN J);
GEN JgetFrob(GEN J);
GEN JgetV(GEN J);
GEN JgetKV(GEN J);
GEN JgetW0(GEN J);
GEN JgetZ(GEN J);
GEN JgetFrobCyc(GEN J);

GEN DivAdd(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEP pe, ulong excess);
GEN DivSub(GEN WA, GEN WB, GEN KV, ulong d, GEN T, GEN p, long e, GEN pe, ulong nIGS);
GEN PicChord(GEN J, GEN WA, GEN WB, long flag);
GEN PicAdd(GEN J, GEN WA, GEN WB);
GEN PicSub(GEN J, GEN WA, GEN WB);
GEN PicNeg(GEN J, GEN W);
GEN PicMul(GEN J, GEN W, GEN n, long flag);
