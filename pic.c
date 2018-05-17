#include "linalg.h"
#include "pic.h"

long Jgetg(GEN J) {return itos(gel(J,1));}
long Jgetd0(GEN J) {return itos(gel(J,2));}
GEN JgetT(GEN J) {return gel(J,3);}
GEN Jgetp(GEN J) {return gel(J,4);}
long Jgete(GEN J) {return itos(gel(J,5));}
GEN Jgetpe(GEN J) {return gel(J,6);}
GEN JgetFrob(GEN J) {return gel(J,7);}
GEN JgetV(GEN J) {return gel(J,8);}
GEN JgetKV(GEN J) {return gel(J,9);}
GEN JgetW0(GEN J) {return gel(J,10);}
GEN JgetZ(GEN J) {return gel(J,11);}
GEN JgetFrobCyc(GEN J) {return gel(J,12);}

GEN DivAdd(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess)
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
			return WAB;
		}
		printf("add(%lu/%lu)",r,d);
		avma = av1;
	}
}

GEN DivSub(GEN WA, GEN WB, GEN KV, ulong d, GEN T, GEN p, long e, GEN pe, ulong nIGS)
{
	pari_sp av1,av = avma;
	unsigned long nZ,P,nE,E,nV,nB,n,r;
	GEN KB,K,col,s,res;
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
		res = FqM_ker(K,T,p); /* TODO faut-il reduire K d'abord? */
		/* TODO take rand subset of eqns */
		/* TODO write fn for that */
		r = lg(res)-1;
		if(r==d)
		{
			return gerepileupto(av,e==1?res:matkerpadic(K,T,p,e));
		}
		printf("sub(%lu/%lu)",r,d);
		avma = av1;
	}
}

GEN PicChord(GEN J, GEN WA, GEN WB, long flag)
{
	pari_sp av = avma;
	GEN WAWB,WAB,s,col,sV,WC,res;
	GEN V,KV,T,p,pe;
	ulong nZ,nV,P,j;
	long g,d0,e;

	V = JgetV(J);
	KV = JgetKV(J);
	T = JgetT(J);
	p = Jgetp(J);
	pe = Jgetpe(J);
	g = Jgetg(J);
	d0 = Jgetd0(J);
	e = Jgete(J);

	WAWB = DivAdd(WA,WB,4*d0+1-g,T,p,e,pe,0);
	WAB = DivSub(V,WAWB,KV,d0+1-g,T,p,e,pe,2);
	/* TODO can free some memory here */
	if(flag & 1) s = RandVec_padic(WAB,T,p,pe);
	else s = gel(WAB,1);
	nZ = lg(s);
	nV = lg(V);
	sV = cgetg(nV,t_MAT);
  for(j=1;j<nV;j++)
	{
		col = cgetg(nZ,t_COL);
		for(P=1;P<nZ;P++)
		{
			gel(col,P) = Fq_mul(gel(s,P),gcoeff(V,P,j),T,pe);
		}
		gel(sV,j) = col;
	}
  WC = DivSub(WAB,sV,KV,2*d0+1-g,T,p,e,pe,2);

	if(flag & 2)
	{
		res = cgetg(3,t_VEC);
		gel(res,1) = gcopy(WC);
		gel(res,2) = gcopy(s);
		return gerepileupto(av,res);
	}
	else
	{
		return gerepileupto(av,WC);
	}
}

GEN PicAdd(GEN J, GEN WA, GEN WB)
{
	pari_sp av = avma;
	GEN W;
	W = PicChord(J,WA,WB,0);
	W = PicChord(J,W,JgetW0(J),0);
	return gerepileupto(av,W);
}

GEN PicSub(GEN J, GEN WA, GEN WB)
{
	pari_sp av = avma;
	GEN W;
	W = PicChord(J,WB,JgetW0(J),0);
	W = PicChord(J,W,WA,0);
	return gerepileupto(av,W);
}

GEN PicNeg(GEN J, GEN W) { return PicChord(J,W,JgetW0(J),0); }

GEN PicMul(GEN J, GEN W, GEN n, long flag)
{
	pari_sp av = avma;
	GEN W2;

	if(signe(n)==0) return JgetW0(J);
	if(gequal(n,gen_1)) return gcopy(W);
	flag = (flag & 1);
	if(gequal(n,gen_m1)) return PicChord(J,W,JgetW0(J),flag);
	if(mpodd(n))
	{
		W2 = PicMul(J,W,mpdiv(mpadd(n,gen_1),gen_2),0);
		W2 = PicChord(J,W2,W2,0);
		W2 = PicChord(J,W2,W,flag);
	}
	else
	{
		W2 = PicMul(J,W,mpdiv(n,gen_m2),0);
		W2 = PicChord(J,W2,W2,flag);
	}
	return gerepileupto(av,W2);
}

GEN PicFrob(GEN J, GEN W)
{
	GEN W2,T,pe,Frob,FrobCyc;
	ulong o,i,j,k,c,nW,nZ,nCyc;

	T = JgetT(J);
	pe = Jgetpe(J);
	Frob = JgetFrob(J);
	FrobCyc = JgetFrobCyc(J);
	nW = lg(W);
	nZ = lg(JgetZ(J));
	nCyc = lg(FrobCyc);
	Frob = JgetFrob(J);

	W2 = cgetg(nW,t_MAT);
	for(j=1;j<nW;j++)
	{
		gel(W2,j) = cgetg(nZ,t_COL);
	}

	i = 0;
	for(o=1;o<nCyc;o++)
	{
		c = FrobCyc[o];
		for(k=1;k<c;k++)
		{
			for(j=1;j<nW;j++)
			{
				gcoeff(W2,i+k+1,j) = FpX_FpXQ_eval(gcoeff(W,i+k,j),Frob,T,pe);
			}
		}
		for(j=1;j<nW;j++)
		{
			gcoeff(W2,i+1,j) = FpX_FpXQ_eval(gcoeff(W,i+c,j),Frob,T,pe);
		}
		i += c;
	}
	return W2;
}
