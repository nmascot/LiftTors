#include "linalg.h"

GEN DivAdd(GEN WA, GEN WB, unsigned long d, GEN T, GEN p, long e, GEN pe, unsigned long excess)
{
	pari_sp av,av1,av0=avma;
	unsigned long nZ,j,P,r;
	GEN WAB,s,t,st;
	nZ = lg(gel(WA,1));
	WAB = cgetg(d+excess+1,t_MAT);
	while(1)
	{
		av1 = avma;
		for(j=1;j<=d+excess;j++)
		{ 
      av = avma;
			s = RandVec_padic(WA,T,p,pe); /* random fn in WA */
			t = RandVec_padic(WB,T,p,pe); /* random fn in WB */
      st = cgetg(nZ,t_COL); /* Product */
			for(P=1;P<nZ;P++)
			{
				gel(st,P) = Fq_mul(gel(s,P),gel(t,P),T,pe);
			}
			gel(WAB,j) = gerepileupto(av,st);
		}
		r = FqM_rank(WAB,T,p); /* TODO faut-il reduire WAB d'abord? */
		if(r==d)
		{
			if(excess)
			{
				WAB = gerepileupto(av0,matimagepadic(WAB,T,p,e));
			}
			return WA;
		}
		printf("add%lu/%lu",r,d);
		avma = av1;
	}
}

GEN DivSub(GEN WA, GEN WB, GEN KV, unsigned long d, GEN T, GEN p, long e, GEN pe, unsigned long nIGS)
{
	pari_sp av1,av = avma;
	unsigned long nZ,P,nE,E,nV,nB,n,r;
	GEN KB,K,col,s;
	nZ = lg(KV);
	nV = lg(gel(KV,1))-1;
	KB = mateqnpadic(WB,T,p,e);
  nB = lg(gel(KB,1))-1;
	/* Prepare a mat K of size a v stack of KV + nIGS copies of KB */
	/* and copy KV at the top */
	nE = nV + nIGS*nB;
	K = cgetg(nZ,t_MAT);
	for(P=1;P<nZ;P++)
	{
		col = cgetg(nE+1,t_COL);
		for(E=1;E<=nV;E++)
		{
			gel(col,E) = gcoeff(KV,E,P);
		}
		gel(K,P) = col;
	}
	av1 = avma;
	while(1)
	{
		/* nIGS times, take rand s in WA, and stack s.KB down K */
		for(n=1;n<=nIGS;n++)
		{
			s = RandVec_padic(WA,T,p,pe);
			for(E=1;E<=nB;E++)
			{
				for(P=1;P<nZ;P++)
				{
					gcoeff(K,nV+(n-1)*nB+E,P) = Fq_mul(gel(s,P),gcoeff(KB,E,P),T,pe);
				}
			}
		}
		r = lg(FqM_ker(K,T,p)); /* TODO faut-il reduire K d'abord? */
		/* TODO take rand subset of eqns */
		/* TODO write fn for that */
		/* TODO case e==1 */
		if(r==d)
		{
			return gerepileupto(av,matkerpadic(K,T,p,e));
		}
		printf("sub%lu/%lu",r,d);
		avma = av1;
	}
}
