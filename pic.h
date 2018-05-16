#include<pari/pari.h>

#define lgJ 13

ulong Jgetg(GEN J);
void Jsetg(GEN J, ulong g);
ulong Jgetd0(GEN J);
void Jsetd0(GEN J, ulong d0);
GEN JgetT(GEN J);
void JsetT(GEN J, GEN T);
GEN Jgetp(GEN J);
void Jsetp(GEN J, GEN p);
long Jgete(GEN J);
void Jsete(GEN J, ulong e);
GEN Jgetpe(GEN J);
void Jsetpe(GEN J, GEN pe);
GEN JgetFrob(GEN J);
void JsetFrob(GEN J, GEN Frob);
GEN JgetV(GEN J);
void JsetV(GEN J, GEN V);
GEN JgetKV(GEN J);
void JsetKV(GEN J, GEN KV);
GEN JgetW0(GEN J);
void JsetW0(GEN J, GEN W0);
GEN JgetZ(GEN J);
void JsetZ(GEN J, GEN Z);
GEN JgetFrobCyc(GEN J);
void JsetFrobCyc(GEN J, GEN FrobCyc);

GEN DivAdd(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess);
GEN DivSub(GEN WA, GEN WB, GEN KV, ulong d, GEN T, GEN p, long e, GEN pe, ulong nIGS);
GEN PicChord(GEN J, GEN WA, GEN WB, long flag);
GEN PicAdd(GEN J, GEN WA, GEN WB);
GEN PicSub(GEN J, GEN WA, GEN WB);
GEN PicNeg(GEN J, GEN W);
GEN PicMul(GEN J, GEN W, GEN n, long flag);
