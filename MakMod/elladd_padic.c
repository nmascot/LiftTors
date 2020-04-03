#include<pari/pari.h>

GEN elladd_padic(GEN a4, GEN P, GEN Q, GEN T, GEN pe, GEN p, long e)
{
	pari_sp av = avma;
	GEN P0,xP,yP,xQ,yQ,dx,dy,l,m,xR,yR,R;
	
	P0 = mkvec(gen_0); /* [0] */
	if(gequal(P,P0))
	{
		avma = av;
		return Q;
	}
	if(gequal(Q,P0))
	{
		avma = av;
		return P;
	}
	xP = gel(P,1);
	yP = gel(P,2);
	xQ = gel(Q,1);
	yQ = gel(Q,2);
	if(gequal(FpX_red(xP,p),FpX_red(xQ,p)))
	{
		/* Dangerous case */
		if(gequal(FpX_red(xP,pe),FpX_red(xQ,pe))==0)
			pari_err(e_IMPL,"case P!=Q but P=Q mod p");
		if(gequal0(FpX_red(ZX_add(yP,yQ),pe)))
			return gerepileupto(av,P0);
		dx = ZX_Z_mul(yP,gen_2);
		dy = Fq_sqr(xP,T,pe);
		dy = ZX_Z_mul(dy,utoi(3));
		dy = gadd(dy,a4);
	}
	else
	{
		dx = ZX_sub(xQ,xP);
		dy = ZX_sub(yQ,yP);
	}
	l = ZpXQ_div(dy,dx,T,pe,p,e);
	m = Fq_mul(l,xP,T,pe);
	m = ZX_sub(yP,m);
	xR = ZX_add(xP,xQ);
	xR = FpX_sub(Fq_sqr(l,T,pe),xR,pe);
	yR = Fq_mul(l,xR,T,pe);
	yR = ZX_add(yR,m);
	yR = FpX_neg(yR,pe);
	R = mkvec2(xR,yR);
	return gerepilecopy(av,R);
}
