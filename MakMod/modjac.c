#include "zn.h"

GEN l2(GEN EN, GEN P, GEN Q, GEN T, GEN pe, GEN p, long e)
{
	P = GetCoef(EN,P);
	Q = GetCoef(EN,Q);
	return ZpXQ_div(ZX_sub(gel(Q,2),gel(P,2)),ZX_sub(gel(Q,1),gel(P,1)),T,pe,p,e);
}

GEN l1(GEN EN, GEN P, GEN Q, GEN T, GEN pe, GEN p, long e)
{
	pari_sp av = avma;
	ulong N,n;
	GEN S,nP;

	/* TODO methode Kamal addchain */
	N = lg(EN)-1;
	S = l2(EN,P,Q,T,pe,p,e);
	nP = P;
	for(n=1;n<N;n++)
	{
		S = ZX_add(S,l2(EN,P,ZC_add(Q,nP),T,pe,p,e));
		nP = ZC_add(nP,P);
	}
	return gerepileupto(av,S);
}
