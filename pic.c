#include "linalg.h"
#include "exp.h"
#include "pic.h"

/*
Structure of a Jacobian:
1: f, equation of curve (0 if not available)
2: g, genus
3: d0, degree of line bundle B
4: If B=O_X(D0), basis of RR spaces L(D0), L(2D0), L(2D0-E1), L(2D0-E2) in terms of x and y (the latter 2 used for EvalData); else []
5: T, poly in Z[t] such that Qq = Qp[t]/T
6: p
7: e, such that we work in J(Zq/p^e)
8: p^e
9: FrobMat, gives matrix of Frob on the power basis 1,t,t²,... of Qq
10: V = [V1,V2,V3] where Vn = H0(n*B). (Having a nice basis of V2 improves the evaluaton map.)
11: KV = [KV1,KV2,KV3], where KVn = equation matrix for Vn
12: W0 = f*V1 for some f in V1, subspace of V2 representing the origin
13: EvalData [U1,U2,I,M]: pair of subspaces Ui of the form V2(-E) with E effective of degree d0-g, used for construction of eval map, then vecsmall I of row indices, and matrix M such that v in V should be taken to M*(v_I) for evaluation  
14: If B=O_X(D0), vector Z of points at which the sections are evaluated; else []
15: FrobCyc, permutation describing the action of Frob on Z

Note: usually, B=O_X(D0), in which case W0=V1=L(D0).
Note: if either of f, one of the RR spaces L(...), or Z are not available, then the p-adic accuracy cannot be increased.
Note: if accuracy is increased, assume that in 13, the block froms by the I-rows of V2 is invertible, and that M is the inverse of that block.
*/

GEN Jgetf(GEN J) {return gel(J,1);}
long Jgetg(GEN J) {return itos(gel(J,2));}
long Jgetd0(GEN J) {return itos(gel(J,3));}
GEN JgetL(GEN J) {return gel(J,4);}
GEN JgetT(GEN J) {return gel(J,5);}
GEN Jgetp(GEN J) {return gel(J,6);}
long Jgete(GEN J) {return itos(gel(J,7));}
GEN Jgetpe(GEN J) {return gel(J,8);}
GEN JgetFrobMat(GEN J) {return gel(J,9);}
GEN JgetV(GEN J, ulong n) {return gmael(J,10,n);}
GEN JgetV_all(GEN J) {return gel(J,10);}
GEN JgetKV(GEN J, ulong n) {return gmael(J,11,n);}
GEN JgetKV_all(GEN J) {return gel(J,11);}
GEN JgetW0(GEN J) {return gel(J,12);}
GEN JgetEvalData(GEN J) {return gel(J,13);}
GEN JgetZ(GEN J) {return gel(J,14);}
GEN JgetFrobCyc(GEN J) {return gel(J,15);}

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
	GEN Je,p,pe;
	ulong i;
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
	gel(Je,9) = FpXM_red(JgetFrobMat(J),pe);
	gel(Je,10) = cgetg(4,t_VEC);
	for(i=1;i<=3;i++) gmael(Je,10,i) = FpXM_red(gmael(J,10,i),pe);
	gel(Je,11) = cgetg(4,t_VEC);
	for(i=1;i<=3;i++) gmael(Je,11,i) = FpXM_red(gmael(J,11,i),pe);
	gel(Je,12) = FpXM_red(JgetW0(J),pe);
	gel(Je,13) = cgetg(5,t_VEC);
	gmael(Je,13,1) = FpXM_red(gmael(J,13,1),pe);
	gmael(Je,13,2) = FpXM_red(gmael(J,13,2),pe);
	gmael(Je,13,3) = gcopy(gmael(J,13,3));
	gmael(Je,13,4) = FpXM_red(gmael(J,13,4),pe);
	gel(Je,14) = FpXT_red(JgetZ(J),pe);
	gel(Je,15) = gcopy(JgetFrobCyc(J));
	return gerepileupto(av,Je);
}

GEN DivMul(GEN f, GEN W, GEN T, GEN pe)
{
	ulong nW,nZ,i,j;
	GEN fW;
	nW = lg(W);
	nZ = lg(f);
	fW = cgetg(nW,t_MAT);
	for(j=1;j<nW;j++)
	{
		gel(fW,j) = cgetg(nZ,t_COL);
		for(i=1;i<nZ;i++)
			gcoeff(fW,i,j) = Fq_mul(gel(f,i),gcoeff(W,i,j),T,pe);
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
{ /* { v in Ker KV | v*WA c WB } */
	pari_sp av1,av = avma;
	unsigned long nZ,P,nE,E,nV,nB,n,r;
	GEN KB,K,s,res;
	nZ = lg(KV);
	nV = lg(gel(KV,1))-1;
	KB = mateqnpadic(WB,T,pe,p,e);
	nB = lg(gel(KB,1))-1;
	/* Prepare a mat K of size a v stack of KV + nIGS copies of KB */
	/* and copy KV at the top */
	nE = nV + nIGS*nB;
	K = cgetg(nZ,t_MAT);
	for(P=1;P<nZ;P++)
	{
		gel(K,P) = cgetg(nE+1,t_COL);
		for(E=1;E<=nV;E++)
			gcoeff(K,E,P) = gcoeff(KV,E,P);
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
					gcoeff(K,nV+(n-1)*nB+E,P) = Fq_mul(gel(s,P),gcoeff(KB,E,P),T,pe);
			}
		}
		res = matkerpadic(K,T,pe,p,e);
		r = lg(res)-1;
		if(r==d) return gerepileupto(av,res);
		if(DEBUGLEVEL) err_printf("sub(%lu/%lu)",r,d);
		avma = av1;
	}
}

GEN DivSub_safe(GEN WA, GEN WB, GEN KV, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
  unsigned long nZ,P,nE,E,nV,nA,nB,n;
  GEN KB,K,res;
  nZ = lg(KV);
  nV = lg(gel(KV,1))-1;
  KB = mateqnpadic(WB,T,pe,p,e);
  nB = lg(gel(KB,1))-1;
	nA = lg(WA)-1;

	/* Prepare a mat K of size a v stack of KV + nA copies of KB */
  /* and copy KV at the top */
  nE = nV + nA*nB;
  K = cgetg(nZ,t_MAT);
  for(P=1;P<nZ;P++)
  {
    gel(K,P) = cgetg(nE+1,t_COL);
    for(E=1;E<=nV;E++)
      gcoeff(K,E,P) = gcoeff(KV,E,P);
  }
  /* Stack a.KB down K for a in basis of WA */
  for(n=1;n<=nA;n++)
  {
    for(E=1;E<=nB;E++)
    {
      for(P=1;P<nZ;P++)
        gcoeff(K,nV+(n-1)*nB+E,P) = Fq_mul(gcoeff(WA,P,n),gcoeff(KB,E,P),T,pe);
    }
  }
  res = matkerpadic_safe(K,T,p,e);
  return gerepileupto(av,res);
}

GEN PicNeg(GEN J, GEN W, long flag)
{ /* flag: 1: choose s randomly, 2: also return s */
  pari_sp av = avma;
  GEN s,sV,WN,res;
  GEN V,KV,T,p,pe;
  long g,d0,e;

  V = JgetV(J,2);
  KV = JgetKV(J,2);
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
	GEN T,pe,p,w,V,wV,KV,W2;
	long e;
	ulong nW;
	long res;

	JgetTpe(J,&T,&pe,&p,&e);
	V = JgetV(J,2);
	KV = JgetKV(J,2);
	nW = lg(W)-1;

	do
		w = RandVec_1(W,pe);
	while(gequal0(w));
	wV = DivMul(w,V,T,pe);

	W2 = DivSub_safe(W,wV,KV,T,p,e,pe);
	res = (lg(W2)-1==nW?1:0);
	avma = av;
	return res;
}


GEN PicChord(GEN J, GEN WA, GEN WB, long flag)
{ /* flag: 1: choose s randomly, 2: also return s */
	pari_sp av = avma;
	GEN WAWB,WAB,s,sV,WC,res;
	GEN V,KV,KV3,V1,T,p,pe;
	long g,d0,e;

	V = JgetV(J,2);
	KV = JgetKV(J,2);
	KV3 = JgetKV(J,3);
	V1 = JgetV(J,1);
	JgetTpe(J,&T,&pe,&p,&e);
	g = Jgetg(J);
	d0 = Jgetd0(J);

	/* L(4D0-A-B) */
	WAWB = DivAdd(WA,WB,2*d0+1-g,T,p,e,pe,0);
	/* L(3D0-A-B) */
	WAB = DivSub(V1,WAWB,KV3,d0+1-g,T,p,e,pe,2);
	if(gc_needed(av,1)) WAB = gerepileupto(av,WAB);
	if(flag & 1) s = RandVec_padic(WAB,T,p,pe);
	else s = gel(WAB,1);
	/* s in WB: (s) = -3D0+A+B+C */
	/* L(5D0-A-B-C) */
	sV = DivMul(s,V,T,pe);
	/* L(2D0-C) */
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
{ /* Matrix of Frob on basis 1,t,t²,... of Z[t]/(T,p^e) */
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
	y[1] = 0;
	setsigne(y,1);
	setvarn(y,var);
	for(i=1;i<=d;i++) gel(y,i+1) = gel(cy,i);
	y = normalizepol(y);
	return gerepilecopy(av,y);
}

GEN PicFrob(GEN J, GEN W)
{
	GEN W2,T,pe,FrobMat,FrobCyc;
	ulong nW,nZ,i,j;

	T = JgetT(J);
	pe = Jgetpe(J);
	FrobMat = JgetFrobMat(J);
	FrobCyc = JgetFrobCyc(J);
	nW = lg(W);
	nZ = lg(FrobCyc);

	W2 = cgetg(nW,t_MAT);
	for(j=1;j<nW;j++)
	{
		gel(W2,j) = cgetg(nZ,t_COL);
		for(i=1;i<nZ;i++)
			gcoeff(W2,FrobCyc[i],j) = Frob(gcoeff(W,i,j),FrobMat,T,pe);
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
	GEN s,sWB,K,KV,T,p,pe;

	JgetTpe(J,&T,&pe,&p,&e);
	KV = JgetKV(J,2);

	/* Take s in WA: (s) = -2D0+A+A1 */
	s = gel(WA,1);
	/* Compute s*WB = L(4D0-A-B-A1) */
	sWB = DivMul(s,WB,T,pe);
	/* Find { v in V | v*WA c s*WB } = L(2D0-B-A1) */
	/* This space is nontrivial iff. A~B */
	K = DivSub_safe(WA,sWB,KV,T,p,e,pe);
	r = lg(K)-1;
	avma = av;
	return r;
}

long PicIsZero(GEN J, GEN W)
{
	pari_sp av = avma;
	long r,e;
	GEN K,V1,KV1,T,p,pe;

	JgetTpe(J,&T,&pe,&p,&e);
	V1 = JgetV(J,1);
	KV1 = JgetKV(J,1);

	K = DivSub_safe(V1,W,KV1,T,p,e,pe);
	r = lg(K)-1;
	avma = av;
	return r;
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
	V = JgetV(J,2);
	KV = JgetKV(J,2);
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
	K = matkerpadic(K,T,pe,p,e);
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
	V = JgetV(J,2);
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
	K = matkerpadic(K,T,pe,p,e);
	K = FqM_mul(V,K,T,pe);
	return gerepileupto(av,K);
}
