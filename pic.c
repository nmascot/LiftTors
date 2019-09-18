#include "linalg.h"
#include "exp.h"
#include "pic.h"

GEN Jgetf(GEN J) {return gel(J,1);}
long Jgetg(GEN J) {return itos(gel(J,2));}
long Jgetd0(GEN J) {return itos(gel(J,3));}
GEN JgetL(GEN J) {return gel(J,4);}
GEN JgetT(GEN J) {return gel(J,5);}
GEN Jgetp(GEN J) {return gel(J,6);}
long Jgete(GEN J) {return itos(gel(J,7));}
GEN Jgetpe(GEN J) {return gel(J,8);}
GEN JgetFrobMat(GEN J) {return gel(J,9);}
GEN JgetV(GEN J) {return gel(J,10);}
GEN JgetKV(GEN J) {return gel(J,11);}
GEN JgetW0(GEN J) {return gel(J,12);}
GEN JgetZ(GEN J) {return gel(J,13);}
GEN JgetFrobCyc(GEN J) {return gel(J,14);}
GEN JgetV3(GEN J) {return gel(J,15);}
GEN JgetKV3(GEN J) {return gel(J,16);}
GEN JgetEvalData(GEN J) {return gel(J,17);}

void JgetTpe(GEN J, GEN* T, GEN* pe, GEN* p, long* e)
{
	*T = gel(J,5);
	*p = gel(J,6);
	*e = itos(gel(J,7));
	*pe = gel(J,8);
}


GEN PicRed(GEN J, ulong e)
{
	pari_sp av = avma;
	GEN Je,p,pe,U,Ue;
	ulong i,j,n;
	if(Jgete(J)<e) pari_err(e_MISC,"Cannot perform this reduction");
	Je = cgetg(lgJ+1,t_VEC);
	gel(Je,1) = gcopy(Jgetf(J));
	gel(Je,2) = stoi(Jgetg(J));
	gel(Je,3) = stoi(Jgetd0(J));
	gel(Je,4) = gcopy(JgetL(J));
	gel(Je,5) = gcopy(JgetT(J));
	gel(Je,6) = p = gcopy(Jgetp(J));
	gel(Je,7) = utoi(e);
	gel(Je,8) = pe = powiu(p,e);
	gel(Je,9) = FpM_red(JgetFrobMat(J),pe);
	gel(Je,10) = FpXM_red(JgetV(J),pe);
	gel(Je,11) = FpXM_red(JgetKV(J),pe);
	gel(Je,12) = FpXM_red(JgetW0(J),pe);
	gel(Je,13) = FpXT_red(JgetZ(J),pe);
	gel(Je,14) = gcopy(JgetFrobCyc(J));
	gel(Je,15) = FpXM_red(JgetV3(J),pe);
	gel(Je,16) = FpXM_red(JgetKV3(J),pe);
	U = JgetEvalData(J);
	Ue = cgetg(3,t_VEC);
	for(i=1;i<=2;i++)
	{
		n = lg(gel(U,i));
		gel(Ue,i) = cgetg(n,t_VEC);
		for(j=1;j<n;j++)
			gmael(Ue,i,j) = FpXM_red(gmael(U,i,j),pe);
	}
	gel(Je,17) = Ue;
	return gerepileupto(av,Je);
}

GEN DivMul(GEN f, GEN W, GEN T, GEN pe)
{
	ulong nW,nZ,i,j;
	GEN fW,col;
	nW = lg(W);
	nZ = lg(f);
	fW = cgetg(nW,t_MAT);
	for(j=1;j<nW;j++)
	{
		col = cgetg(nZ,t_COL);
		for(i=1;i<nZ;i++)
		{
			gel(col,i) = Fq_mul(gel(f,i),gcoeff(W,i,j),T,pe);
		}
		gel(fW,j) = col;
	}
	return fW;
}

/* Does products s_i * t_j with random i and j */
/* Not uniform enough, do not use. */
/*GEN DivAdd0(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess)
{
	pari_sp av,av1,av0=avma;
  unsigned long nZ,nA,nB,j,P,r;
  GEN WAB,s,t,st;
  nZ = lg(gel(WA,1));
	nA = lg(WA)-1;
	nB = lg(WB)-1;
  WAB = cgetg(d+excess+1,t_MAT);
  while(1)
  {
    av1 = avma;
    for(j=1;j<=d+excess;j++)
    {
      av = avma;
			s = gel(WA,1+random_Fl(nA));
			t = gel(WB,1+random_Fl(nB));
      st = cgetg(nZ,t_COL);
      for(P=1;P<nZ;P++)
      {
        gel(st,P) = Fq_mul(gel(s,P),gel(t,P),T,pe);
      }
      gel(WAB,j) = gerepileupto(av,st);
    }
    r = FqM_rank(WAB,T,p);
    if(r==d)
    {
      if(excess)
      {
        WAB = gerepileupto(av0,matimagepadic(WAB,T,p,e));
      }
      return WAB;
    }
    printf("add0(%lu/%lu)",r,d);
    avma = av1;
  }
}*/

GEN DivAdd(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess)
{ /* Does products s*t, with s=sum n_i s_i, n_i = 0 50%, -1 25%, +1 25%; similarly for t */
	/* Fails from time to time but overall good speedup */
	pari_sp av=avma;
  unsigned long nZ,j,P,r;
  GEN WAB,s,t,st;
  nZ = lg(gel(WA,1));
  while(1)
  {
  	WAB = cgetg(d+excess+1,t_MAT);
    for(j=1;j<=d+excess;j++)
    {
      s = RandVec_1(WA,pe); /* random fn in WA */
      t = RandVec_1(WB,pe); /* random fn in WB */
      st = cgetg(nZ,t_COL); /* Product */
      for(P=1;P<nZ;P++)
      {
        gel(st,P) = Fq_mul(gel(s,P),gel(t,P),T,pe);
      }
      gel(WAB,j) = st;
    }
		WAB = FqM_image(WAB,T,p);
		r = lg(WAB)-1;
    if(r==d)
      return gerepileupto(av,WAB);
    if(DEBUGLEVEL) err_printf("add1(%lu/%lu)",r,d);
		excess++;
    avma = av;
  }
}

/* Does products s*t, s and t chosen at random */
/* Safest but slowest */
/*GEN DivAdd(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess)
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
			s = RandVec_padic(WA,T,p,pe);
			t = RandVec_padic(WB,T,p,pe);
			st = cgetg(nZ,t_COL);
			for(P=1;P<nZ;P++)
			{
				gel(st,P) = Fq_mul(gel(s,P),gel(t,P),T,pe);
			}
			gel(WAB,j) = gerepileupto(av,st);
		}
		r = FqM_rank(WAB,T,p);
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
}*/

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
			s = RandVec_padic(WA,T,p,pe); /* Note: RandVec_1 would be slower here */
			for(E=1;E<=nB;E++)
			{
				for(P=1;P<nZ;P++)
				{
					gcoeff(K,nV+(n-1)*nB+E,P) = Fq_mul(gel(s,P),gcoeff(KB,E,P),T,pe);
				}
			}
		}
		res = matkerpadic(K,T,p,e);
		r = lg(res)-1;
		if(r==d) return gerepileupto(av,res);
		if(DEBUGLEVEL) err_printf("sub(%lu/%lu)",r,d);
		avma = av1;
	}
}

GEN DivSub_IGS(GEN GA, GEN WB, GEN KV, GEN T, GEN p, long e, GEN pe)
{
  pari_sp av = avma;
  ulong nIGS,nZ,P,nE,E,nV,nB,n;
  GEN KB,K,col,res;
	nIGS = lg(GA)-1;
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
  /* nIGS times, take rand s in WA, and stack s.KB down K */
  for(n=1;n<=nIGS;n++)
  {
    for(E=1;E<=nB;E++)
    {
      for(P=1;P<nZ;P++)
      {
        gcoeff(K,nV+(n-1)*nB+E,P) = Fq_mul(gcoeff(GA,P,n),gcoeff(KB,E,P),T,pe);
      }
    }
  }
  res = matkerpadic(K,T,p,e);
  return gerepileupto(av,res);
}

GEN PicNeg(GEN J, GEN W, long flag)
{ /* flag: 1: choose s randomly, 2: also return s */
  pari_sp av = avma;
  GEN s,sV,WN,res;
  GEN V,KV,T,p,pe;
  long g,d0,e;

  V = JgetV(J);
  KV = JgetKV(J);
  JgetTpe(J,&T,&pe,&p,&e);
  g = Jgetg(J);
  d0 = Jgetd0(J);

  /* (s) = -2_0-D-N */
  if(flag & 1) s = RandVec_padic(W,T,p,pe);
  else s = gel(W,1);
  sV = DivMul(s,V,T,pe); /* L(4D_0-D-N) */
  WN = DivSub(W,sV,KV,d0+1-g,T,p,e,pe,2); /* L(2D_0-N) */

  if(flag & 2)
  {
    res = cgetg(3,t_VEC);
    gel(res,1) = gcopy(WN);
    gel(res,2) = gcopy(s);
    return gerepileupto(av,res);
  }
  else
  {
    return gerepileupto(av,WN);
  }
}

long PicMember(GEN J, GEN W)
{
	pari_sp av = avma;
	GEN T,pe,p,w,V,wV,KV,KwV,K,W2;
	long e;
	ulong nW,nZ,nKV,nE;
	ulong P,E,n;
	long res;

	JgetTpe(J,&T,&pe,&p,&e);
	V = JgetV(J);
	KV = JgetKV(J);
	nZ = lg(KV)-1;
	nKV = lg(gel(KV,1))-1;
	nW = lg(W)-1;

	do
		w = RandVec_1(W,pe);
	while(gequal0(w));
	wV = DivMul(w,V,T,pe);
	KwV = mateqnpadic(wV,T,p,e);
  /* Prepare a mat K of size a v stack of KV + nW copies of KW */
  /* and copy KV at the top */
  nE = nKV*(nW+1);
  K = cgetg(nZ+1,t_MAT);
  for(P=1;P<=nZ;P++)
  {
    gel(K,P) = cgetg(nE+1,t_COL);
    for(E=1;E<=nKV;E++)
    {
      gcoeff(K,E,P) = gcoeff(KV,E,P);
			for(n=1;n<=nW;n++)
			{
				gcoeff(K,n*nKV+E,P) = Fq_mul(gcoeff(W,P,n),gcoeff(KwV,E,P),T,pe);
			}
		}
  }
  W2 = matkerpadic_safe(K,T,p,e);
	res = (lg(W2)-1==nW?1:0);
	avma = av;
	return res;
}


GEN PicChord(GEN J, GEN WA, GEN WB, long flag)
{ /* flag: 1: choose s randomly, 2: also return s */
	pari_sp av = avma;
	GEN WAWB,WAB,s,sV,WC,res;
	GEN V,KV,KV3,W0,T,p,pe;
	long g,d0,e;

	V = JgetV(J);
	KV = JgetKV(J);
	KV3 = JgetKV3(J);
	W0 = JgetW0(J);
	JgetTpe(J,&T,&pe,&p,&e);
	g = Jgetg(J);
	d0 = Jgetd0(J);

	WAWB = DivAdd(WA,WB,2*d0+1-g,T,p,e,pe,0);
	WAB = DivSub(W0,WAWB,KV3,d0+1-g,T,p,e,pe,2);
	/* TODO can free some memory here */
	if(flag & 1) s = RandVec_padic(WAB,T,p,pe);
	else s = gel(WAB,1);
	sV = DivMul(s,V,T,pe);
	WC = DivSub(WAB,sV,KV,d0+1-g,T,p,e,pe,2);

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
	W = PicNeg(J,W,0);
	return gerepileupto(av,W);
}

GEN PicSub(GEN J, GEN WA, GEN WB)
{
	pari_sp av = avma;
	GEN W;
	W = PicNeg(J,WB,0);
	W = PicChord(J,W,WA,0);
	return gerepileupto(av,W);
}

GEN PicMul(GEN J, GEN W, GEN n, long flag)
{ /* flag: 2: sign matters, 1: pass to PicChord and PicNeg */
	pari_sp av = avma;
	GEN C,Wlist,WA,WB;
	ulong nC,i;
	long a,b;

	if(gequal0(n)) return gcopy(JgetW0(J));
	if(gequal(n,gen_1)) return gcopy(W);
	C = AddChain(n,flag&2);
	nC = lg(C);
	if(DEBUGLEVEL)
	{
		if(flag&2)
			pari_printf("   PicMul : Mul by %Ps in %lu steps\n",n,nC-2);
		else
			pari_printf("   PicMul : Mul by ±%Ps in %lu steps\n",n,nC-2);
	}
	Wlist = cgetg(nC,t_VEC);
	gel(Wlist,1) = W;
	for(i=2;i<nC;i++)
	{
		a = gmael(C,i,2)[1];
		b = gmael(C,i,2)[2];
		if(a)
		{
			WA = gel(Wlist,a);
			if(b)
			{
				WB = gel(Wlist,b);
				gel(Wlist,i) = PicChord(J,WA,WB,(i==nC-1)&&(flag&1));
			}
			else gel(Wlist,i) = PicNeg(J,WA,(i==nC-1)&&(flag&1));
		}
		else gel(Wlist,i) = b?PicNeg(J,gel(Wlist,b),(i==nC-1)&&(flag&1)):gcopy(JgetW0(J));
		/* Does not respect flag if a==b==0 but this is not supposed to happen */
	}
	return gerepileupto(av,gel(Wlist,nC-1));
}

GEN ZpXQ_FrobMat(GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN F,M,col,Fj;
	long v = gvar(T),d = degpol(T),i,j;
	F = ZpX_Frobenius(T,p,e);
	M = cgetg(d+1,t_MAT);
	col = cgetg(d+1,t_COL);
	gel(col,1) = gen_1;
	for(i=2;i<=d;i++) gel(col,i) = gen_0;
	gel(M,1) = col;
	if(d==1) return gerepileupto(av,M);
	col = cgetg(d+1,t_COL);
	for(i=0;i<d;i++) gel(col,i+1) = polcoef(F,i,v);
	gel(M,2) = col;
	Fj = F;
	for(j=2;j<d;j++)
	{
		Fj = Fq_mul(Fj,F,T,pe);
		col = cgetg(d+1,t_COL);
		for(i=0;i<d;i++) gel(col,i+1) = polcoef(Fj,i,v);
		gel(M,j+1) = col;
	}
	return gerepilecopy(av,M);
}

GEN Frob(GEN x, GEN FrobMat, GEN T, GEN pe)
{
	pari_sp av = avma;
	GEN cx,cy,y;
	long var = gvar(T), d = degpol(T), i;
	cx = cgetg(t_COL,d+1);
	for(i=0;i<d;i++) gel(cx,i+1) = polcoef(x,i,var);
	cy = FpM_FpC_mul(FrobMat,cx,pe);
	y = cgetg(d+2,t_POL);
	setvarn(y,var);
	for(i=1;i<=d;i++) gel(y,i+1) = gel(cy,i);
	y = normalizepol(y);
	return gerepilecopy(av,y);
}

GEN PicFrob(GEN J, GEN W)
{
	GEN W2,T,pe,FrobMat,FrobCyc;
	ulong o,i,j,k,c,nW,nZ,nCyc;

	T = JgetT(J);
	pe = Jgetpe(J);
	FrobMat = JgetFrobMat(J);
	FrobCyc = JgetFrobCyc(J);
	nW = lg(W);
	nZ = lg(JgetZ(J));
	nCyc = lg(FrobCyc);

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
				gcoeff(W2,i+k+1,j) = Frob(gcoeff(W,i+k,j),FrobMat,T,pe);
			}
		}
		for(j=1;j<nW;j++)
		{
			gcoeff(W2,i+1,j) = Frob(gcoeff(W,i+c,j),FrobMat,T,pe);
		}
		i += c;
	}
	return W2;
}

GEN PicFrobPoly(GEN J, GEN W, GEN F)
{
	pari_sp av = avma;
	ulong d,i;
	GEN n,FW,res;

	d = degree(F);
	FW = W;
	n = truecoeff(F,0);
	if(d&1L) n = negi(n);
	res = PicMul(J,W,n,2);
	for(i=1;i<=d;i++)
	{
		FW = PicFrob(J,FW);
		n = truecoeff(F,i);
		if((d+1-i)&1L) n = negi(n);
		res = PicChord(J,res,PicMul(J,FW,n,2),0);
	}
	return gerepileupto(av,res);
}

long PicEq(GEN J, GEN WA, GEN WB)
{
	pari_sp av = avma;
	long e,r;
	GEN s,sWB,KsWB,K,KV,col,T,p,pe;
	ulong P,i,j,nZ,nW,nKV,nKsB,nK;

	JgetTpe(J,&T,&pe,&p,&e);
	KV = JgetKV(J);

	/* Take s in WA */
	s = gel(WA,1);
	nZ = lg(s)-1;
	nW = lg(WA)-1;
	nKV = lg(gel(KV,1))-1;
	nKsB = nZ-nW;
	nK = nKV+nW*nKsB;

	/* Compute s*WB */
	sWB = cgetg(nW+1,t_MAT);
	for(j=1;j<=nW;j++)
	{
		col = cgetg(nZ+1,t_COL);
		for(i=1;i<=nZ;i++)
		{
			gel(col,i) = Fq_mul(gel(s,i),gcoeff(WB,i,j),T,pe);
		}
		gel(sWB,j) = col;
	}

	/* Equations for s*WB */
	KsWB = mateqnpadic(sWB,T,p,e);

	K = cgetg(nZ+1,t_MAT);
	for(j=1;j<=nZ;j++)
	{
		gel(K,j) = cgetg(nK+1,t_COL);
	}

	/* Find { v in V | v*WA c s*WB } */
	/* This space is nontrivial iff. A~B */
	for(j=1;j<=nW;j++)
	{
		for(i=1;i<=nKsB;i++)
		{
			for(P=1;P<=nZ;P++)
			{
				gcoeff(K,(j-1)*nKsB+i,P) = Fq_mul(gcoeff(WA,P,j),gcoeff(KsWB,i,P),T,pe);
			}
		}
	}
	for(i=1;i<=nKV;i++)
	{
		for(P=1;P<=nZ;P++)
		{
			gcoeff(K,nW*nKsB+i,P) = gcoeff(KV,i,P);
		}
	}

	r = lg(matkerpadic_safe(K,T,p,e))-1;

	avma = av;
	return r;
}

long PicIsZero(GEN J, GEN W)
{
	return PicEq(J,W,JgetW0(J));
}

GEN PicChart(GEN J, GEN W, ulong P0, GEN P1) /* /!\ Not Galois-equivariant ! */
{
	pari_sp av = avma;
	ulong d0,g,n1,nZ,nW;
	ulong j,P;
	long e;
	GEN V,KV,T,p,pe;
	GEN K,col,s,sV,U;

	g = Jgetg(J);
	d0 = Jgetd0(J);
	n1 = d0-g;
	V = JgetV(J);
	KV = JgetKV(J);
	nZ = lg(gel(V,1))-1;
	nW = lg(W)-1;
	JgetTpe(J,&T,&pe,&p,&e);

	K = cgetg(nW+1,t_MAT);
	for(j=1;j<=nW;j++)
	{ /* C = sum of pts P0+1, ..., P0+d0-g */
		col = cgetg(n1+1,t_COL);
		for(P=1;P<=n1;P++) gel(col,P) = gcoeff(W,P+P0,j);
		gel(K,j) = col;
	}
	K = matkerpadic(K,T,p,e);
	if(lg(K)!=2)
	{
		pari_printf("Genericity 1 failed in PicChart\n");
		avma = av;
		return NULL;
	}
	s = FqM_FqC_mul(W,gel(K,1),T,pe); /* Generator of L(2D0-D-C) */

	sV = DivMul(s,V,T,pe); /* L(4D0-D-C-E_D), deg E_D = g */
	U = DivSub(W,sV,KV,d0+1-g,T,p,e,pe,2); /* L(2D0-C-E_D) */
	for(j=1;j<=nW;j++) /* Remove zero rows */
	{
		for(P=1;P<=n1;P++) gcoeff(U,P0+P,j) = gcoeff(U,P0+n1+P,j);
		setlg(gel(U,j),nZ-n1+1);
	}
	if(P1)
		U = Subspace_normalize(U,P1,T,pe,p,e,1);
	return gerepilecopy(av,mat2col(U));
}

GEN rand_subset(ulong n, ulong r)
{
	pari_sp av;
	GEN X,S;
	ulong m,i;
	S = cgetg(r+1,t_VECSMALL);
	av = avma;
	X = cgetg(n+1,t_VECSMALL);
	for(i=1;i<=n;i++) X[i] = 1;
	m = 0;
	while(m<r)
	{
		i = random_Fl(n)+1;
		if(X[i])
		{
			X[i] = 0;
			m++;
			S[m] = i;
		}
	}
	avma = av;
	return S;
}

GEN PicRand0(GEN J)
{
	pari_sp av = avma;
	ulong d0,nZ,nV;
	ulong i,j;
	long e;
	GEN T,p,pe,V;
	GEN S,col,K;

	d0 = Jgetd0(J);
	JgetTpe(J,&T,&pe,&p,&e);
	V = JgetV(J);
	nV = lg(V);
	nZ = lg(gel(V,1));

	K = cgetg(nV,t_MAT);
	S = rand_subset(nZ-1,d0);
	for(j=1;j<nV;j++)
	{
		col = cgetg(d0+1,t_COL);
		for(i=1;i<=d0;i++)
		{
			gel(col,i) = gcoeff(V,S[i],j);
		}
		gel(K,j) = col;
	}
	K = matkerpadic(K,T,p,e);
	K = FqM_mul(V,K,T,pe);
	return gerepileupto(av,K);
}
