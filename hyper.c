#include "pic.h"

GEN HyperRandPt(GEN f, GEN T, GEN p, ulong e, GEN pe)
{
	pari_sp av = avma;
	long vT,dT;
	GEN x,y2,y,P;

	vT = varn(T);
	dT = degree(T);
	for(;;)
	{
		avma = av;
		x = random_FpX(dT-1,vT,p);
		y2 = poleval(f,x);
		y2 = FpQX_red(y2,T,pe);
		if(gequal0(FpQX_red(y2,T,p))) continue;
		y = ZpXQ_sqrt(y2,T,p,e);
		if(y == NULL) continue;
		P = mkvec2(x,y);
		return gerepilecopy(av,P);
	}
}

/* Matrix of values of x^i and x^i*y and the points in Ps */
RReval(GEN Ps, ulong n, ulong d, GEN T, GEN pe)
{
	pari_sp av = avma;
	ulong l,m,i,j;
  GEN R,P,x,y,xpow;

	/* Size of matrix */
	l = 2*n-d/2+3;
	m = lg(Ps);
	R = cgetg(l,t_MAT);
	for(j=1;j<l;j++)
	{
		gel(R,j) = cgetg(m,t_COL);
	}

	for(i=1;i<m;i++)
	{
		P = gel(Ps,i);
		x = gel(P,1);
		y = gel(P,2);
		xpow = FpXQ_powers(x,n,T,pe);
		for(j=1;j<=n+1;j++)
		{
			gcoeff(R,i,j) = gel(xpow,j);
		}
		for(j=1;j<=n-d/2+1;j++)
		{
			gcoeff(R,i,j+n+1) = Fq_mul(gel(xpow,j),y,T,pe);
		}
	}

	return gerepilecopy(av,R);
}

GEN HyperInit(GEN f, GEN p, ulong a, ulong e)
{
	pari_sp av = avma;

}
