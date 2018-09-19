#include "pic.h"
#include "linalg.h"
#include "freyruck.h"

GEN PlaneRegRandPt(GEN f, GEN T, GEN p, long e)
{
	pari_sp av = avma;
	long vT,dT;
  GEN x,fx,y,dfx,dy,P;
	vT = varn(T);
  dT = degree(T);
  for(;;)
  {
    avma = av;
    x = random_FpX(dT,vT,p);
		if(ZX_is0mod(x,p)) continue; /* Want x != 0 */
		fx = poleval(f,x);
    y = polrootsmod(fx,mkvec2(T,p));
		if(lg(y)==1) continue; /* No roots */
		y = gel(y,itos(genrand(stoi(lg(y)-1)))+1);
		dfx = RgX_deriv(fx);
		dy = poleval(dfx,y);
		if(gequal0(dy)) continue; /* Bad for Hensel */
		y = gmodulo(liftall(y),T);
		y = gadd(y,zeropadic(p,e));
		y = padicappr(fx,y);
		y = gel(y,1);
		y = liftall(y);
    P = mkvec2(x,y);
    return gerepilecopy(av,P);
  }
}

/* TODO refactor */

GEN Plane_V(long d, long rmin, long rsmax, GEN Z, GEN T, GEN p, long e, GEN pe)
{
	/* Assumes rsmax - rmin >= d-2 */
	pari_sp av = avma;
	GEN V,x,y,x1,xpow,x1pow,ypow,xr,xrys,Fq1;
	long r,s;
	ulong nV,nZ,j,P;

	nZ = lg(Z);
	nV = d*(rsmax-rmin+1)-(d*(d-1))/2;
	V = cgetg(nV+1,t_MAT);
	for(j=1;j<=nV;j++) gel(V,j) = cgetg(nZ,t_COL);
	xpow = cgetg(rsmax+1,t_VEC);
	ypow = cgetg(d,t_VEC);
	x1pow = NULL;
	if(rmin<0) x1pow = cgetg(1-rmin,t_VEC);
	Fq1 = GetFq1(T);

	for(P=1;P<nZ;P++)
	{
		x = gmael(Z,P,1);
		gel(xpow,1) = x;
		for(j=1;j<rsmax;j++) gel(xpow,j+1) = Fq_mul(gel(xpow,j),x,T,pe); 
		y = gmael(Z,P,2);
		gel(ypow,1) = y;
    for(j=1;j<d;j++) gel(ypow,j+1) = Fq_mul(gel(ypow,j),y,T,pe);
		if(rmin<0)
		{
			x1 = ZpXQ_inv(x,T,p,e);
			gel(x1pow,1) = x1;
			for(j=1;j<-rmin;j++) gel(x1pow,j+1) = Fq_mul(gel(x1pow,j),x1,T,pe);
		}
		j = 0;
		for(s=0;s<d;s++)
		{
			for(r=rmin;r+s<=rsmax;r++)
			{
				j++;
				xr = Fq1;
				if(r<0) xr = gel(x1pow,-r);
				if(r>0) xr = gel(xpow,r);
				xrys = xr;
				if(s) xrys = Fq_mul(xr,gel(ypow,s),T,pe);
				gcoeff(V,P,j) = xrys;
			}
		}
	}
	return gerepilecopy(av,V);
}



GEN Vmat_upto(long d, long n, GEN Z, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN V,J,Fq1,x,y,x1,xpow,x1pow,ypow,xr,ys,xrys;
	long r,s,m;
	ulong g,nZ,P,dm,j;

	/* Initialisation */
	g = (d-1)*(d-2);
	g = g/2;
	nZ = lg(Z);
	V = cgetg(n+1,t_VEC);
	J = cgetg(n+1,t_VECSMALL);
	for(m=1;m<=n;m++)
	{
		dm = m*d*(d-2)+1-g;
		gel(V,m) = cgetg(dm+1,t_MAT);
		for(j=1;j<=dm;j++) gmael(V,m,j) = cgetg(nZ,t_COL);
	}
	Fq1 = GetFq1(T);
	xpow = cgetg(n*(d-3)+1,t_VEC);
	x1pow = cgetg(n+1,t_VEC);
	ypow = cgetg(d,t_VEC);

	for(P=1;P<nZ;P++)
	{
		x = gmael(Z,P,1);
		y = gmael(Z,P,2);
		x1 = ZpXQ_inv(x,T,p,e);
		gel(xpow,1) = x;
		for(m=1;m<n*(d-3);m++) gel(xpow,m+1) = Fq_mul(gel(xpow,m),x,T,pe);
		gel(x1pow,1) = x1;
		for(m=1;m<n;m++) gel(x1pow,m+1) = Fq_mul(gel(x1pow,m),x1,T,pe);
		gel(ypow,1) = y;
		for(m=1;m<d-1;m++) gel(ypow,m+1) = Fq_mul(gel(ypow,m),y,T,pe);
		for(m=1;m<=n;m++) J[m] = 1;
		for(s=0;s<d;s++)
		{
			for(r=-n;r+s<=(d-3)*n;r++)
			{
				xr = Fq1;
				if(r>0) xr = gel(xpow,r);
				if(r<0) xr = gel(x1pow,-r);
				if(s)
				{
					ys = gel(ypow,s);
					xrys = Fq_mul(xr,ys,T,pe);
				}
				else xrys = xr;
				for(m=1;m<=n;m++)
				{
					if(r >= -m && r+s <= m*(d-3))
					{
						gcoeff(gel(V,m),P,J[m]) = xrys;
						J[m]++;
					}
				}
			}
		}
	}

	return gerepilecopy(av,V);
}

GEN LinForm(GEN abc, GEN Z, GEN T, GEN pe)
/* /!\ Not memory-clean */
{
	GEN L;
	ulong nZ,P;

	nZ = lg(Z);
	L = cgetg(nZ,t_COL);
	for(P=1;P<nZ;P++)
	{
		gel(L,P) = ZX_Z_mul(gmael(Z,P,1),gel(abc,1));
		gel(L,P) = ZX_add(gel(L,P),ZX_Z_mul(gmael(Z,P,2),gel(abc,2)));
		gel(L,P) = ZX_Z_add(gel(L,P),gel(abc,3));
	}
	return L;
}

GEN V_HE(long d, GEN abc, GEN E, GEN Z, GEN T, GEN p, long e, GEN pe)
/* L(2D0-H-E) where H = hyper sect ax+by+c=0. If deg E=(g-2) then dim 1. */
{
	pari_sp av = avma;
	GEN VH,EvE;

	VH = Plane_V(d,-2,2*(d-3)-1,Z,T,p,e,pe);
	VH = DivMul(LinForm(abc,Z,T,pe),VH,T,pe);
	if(E)
	{
		EvE = Plane_V(d,-2,2*(d-3)-1,E,T,p,e,pe);
		EvE = DivMul(LinForm(abc,E,T,pe),EvE,T,pe);
		VH = FqM_mul(VH,matkerpadic(EvE,T,p,e),T,pe);
	}
	return gerepileupto(av,VH);
}

GEN PlaneInit(GEN f, GEN p, ulong a, long e)
{
	pari_sp avP,av = avma;
  int newpt;
  ulong d1,d2,df,g,d0,nZ,n,ncyc,i;
  GEN vars,pe,t,T,FrobMat,Z,Zp,P,Pp,Q,FrobCyc,x,y,W0,V,KV,V3,KV3,J;

	vars = variables_vecsmall(f);
	d1 = poldegree(f,vars[1]);
	d2 = poldegree(f,vars[2]);
	if(d1>d2) df = d1;
	else df = d2;
	g = (df-1)*(df-2);
	g = g/2;
	d0 = df*(df-2);
	nZ = 5*d0+1;

	t = varlower("t",vars[2]);
  T = liftint(ffinit(p,a,varn(t)));
  pe = powiu(p,e);
	FrobMat = ZpXQ_FrobMat(T,p,e,pe);

  n = ncyc = 0;
  Z = cgetg(nZ+a,t_VEC);
  Zp = cgetg(nZ+a,t_VEC);
  /* TODO sort Zp -> quasilin complexity */
  FrobCyc = cgetg(nZ+1,t_VECSMALL);
  while(n<nZ)
  {
    avP = avma;
    P = PlaneRegRandPt(f,T,p,e);
    /* Already have it ? */
    Pp = FpXV_red(P,p);
    newpt = 1;
    for(i=1;i<=n;i++)
    {
      if(gequal(Pp,gel(Zp,i)))
      {
        newpt = 0;
        avma = avP;
        break;
      }
    }
    if(newpt == 0) continue;
    ncyc++;
    Q = P;
    i = 0;
    do
    {
      i++;
      n++;
      gel(Z,n) = Q;
      gel(Zp,n) = FpXV_red(Q,p);
      x = Frob(gel(Q,1),FrobMat,T,pe);
      y = Frob(gel(Q,2),FrobMat,T,pe);
      Q = mkvec2(x,y);
    } while(!gequal(Q,P));
  FrobCyc[ncyc] = i;
  }
  setlg(Z,n+1);
  setlg(FrobCyc,ncyc+1);

	V = Vmat_upto(df,3,Z,T,p,e,pe);
	W0 = gel(V,1);
	V3 = gel(V,3);
	V = gel(V,2);
  KV = mateqnpadic(V,T,p,e);
  KV3 = mateqnpadic(V3,T,p,e);

  J = mkvecn(lgJ,f,stoi(g),stoi(d0),T,p,stoi(e),pe,FrobMat,V,KV,W0,Z,FrobCyc,V3,KV3);
	return gerepilecopy(av,J);
}

GEN PlaneZeta(GEN f, ulong p)
{
	pari_sp av = avma, av1;
	ulong d1,d2,df,g,a,pa,i,ix,m,n;
	GEN P,N,t,T,vars,Mf,x,fx,f00,Z,D,L,Pi;

	vars = variables_vecsmall(f);
  d1 = poldegree(f,vars[1]);
  d2 = poldegree(f,vars[2]);
  if(d1>d2) df = d1;
  else df = d2;
  g = (df-1)*(df-2);
  g = g/2;
  t = varlower("t",vars[2]);
	Mf = RgXX_to_RgM(f,df+1);
	f00 = cgetg(df+3,t_POL);
	setsigne(f00,1);
	setvarn(f00,0);
	for(i=0;i<=df;i++)
	{
		gel(f00,i+2) = gcoeff(Mf,i+1,df+1-i);
	}
	P = utoi(p);

	N = cgetg(g+1,t_VECSMALL);
	x = cgetg(g+2,t_POL);
	setvarn(x,varn(t));
	for(i=1;i<=g;i++) gel(x,i+1) = gen_0;
	pa = 1;
	for(a=1;a<=g;a++)
	{
		if(a==1) T = t;
		else T=liftint(ffinit(P,a,varn(t)));
		n = 0;
		pa *= p;
		av1 = avma;
		for(ix=0;ix<pa;ix++)
		{
			avma = av1;
			m = ix;
			for(i=1;i<=a;i++)
			{
				gel(x,i+1) = utoi(m%p);
				m = m/p;
			}
			fx = poleval(f,x);
			n += lg(polrootsmod(fx,mkvec2(T,P)))-1;
		}
		if(itou(gel(f00,df+2))%p == 0) n++;
		n += lg(polrootsmod(f00,mkvec2(T,P)))-1;
		N[a] = n;
	}
	Z = cgetg(g+2,t_SER);
	setsigne(Z,1);
	setvarn(Z,0);
	setvalp(Z,1);
	for(i=1;i<=g;i++) gel(Z,i+1) = gdiv(stoi(N[i]),utoi(i));
	Z = gexp(Z,0);
	D = mkpoln(2,gen_m1,gen_1);
	setvarn(D,0);
	Z = gmul(Z,D);
	D = mkpoln(2,gneg(P),gen_1);
  setvarn(D,0);
  Z = gmul(Z,D);
	L = cgetg(2*g+3,t_POL);
  setvarn(L,0);
	setsigne(L,1);
	gel(L,2*g+2) = gen_1;
	for(i=1;i<=g;i++) gel(L,2*g+2-i) = gel(Z,i+2);
	Pi = gen_1;
	for(i=1;i<=g;i++)
	{
		Pi = muliu(Pi,p);
		gel(L,g+2-i) = mulii(gel(L,g+2+i),Pi);
	}

	return gerepilecopy(av,L);
}

GEN PlaneEval(GEN J, GEN W, GEN abc, GEN E, GEN abc1, GEN abc2)
{
	pari_sp av = avma;
	GEN T,p,pe,V,Z,KV;
	long e;
	GEN FqE,VHE,S1,S2,s2,VHE1,VHE2,u1,u2;
	ulong d,d0,g,nE,i;
	
	JgetTpe(J,&T,&pe,&p,&e);
	d0 = Jgetd0(J);
	g = Jgetg(J);
	V = JgetV(J);
	KV = JgetKV(J);
	Z = JgetZ(J);
	for(d = 2; d0 != d*(d-2); d++) {} 

	nE = lg(E);
	FqE = cgetg(nE,t_VEC);
	for(i=1;i<nE;i++)
	{
		gel(FqE,i) = mkvec2(Z2Fq(gmael(E,i,1),T),Z2Fq(gmael(E,i,2),T));
	}
	VHE = V_HE(d,abc,FqE,Z,T,p,e,pe); /* L(2D0-H-E) */
	S1 = DivAdd(W,VHE,3*d0-d-(g-2)+1-g,T,p,e,pe,0); /* L(4D0-D-H-E) */
	S1 = DivSub(V,S1,KV,1,T,p,e,pe,2); /* L(2D0-D-H-E) */
	S2 = DivMul(gel(S1,1),V,T,pe); /* L(4D0-D-H-E-ED) */
	S2 = DivSub(W,S2,KV,d0+1-g,T,p,e,pe,2); /* L(2D0-H-E-ED) */
	S2 = DivAdd(S2,VHE,4*d0-2*d-2*(g-2)-g+1-g,T,p,e,pe,0); /* L(4D0-2H-2E-ED) */
	S2 = DivSub(V,S2,KV,1,T,p,e,pe,2); /* L(2D0-2H-2E-ED) */
	s2 = gel(S2,1);
	VHE1 = V_HE(d,abc1,NULL,Z,T,p,e,pe);
	u1 = PicNorm(J,s2,VHE1);
	VHE2 = V_HE(d,abc2,NULL,Z,T,p,e,pe);
	u2 = PicNorm(J,s2,VHE2);
	u2 = ZpXQ_inv(u2,T,p,e);
	return gerepileupto(av,Fq_mul(u1,u2,T,pe));
}

GEN PlaneEval0_data(GEN J, GEN abc1, GEN E1, GEN abc2, GEN E2)
{
	pari_sp av=avma;
	GEN T,p,pe,Z,FqE1,FqE2,VHE1,VHE2;
	long e;
	ulong d,d0,nE,i;

	JgetTpe(J,&T,&pe,&p,&e);
  d0 = Jgetd0(J);
  Z = JgetZ(J);
  for(d = 2; d0 != d*(d-2); d++) {}

  nE = lg(E1);
  FqE1 = cgetg(nE,t_VEC);
  FqE2 = cgetg(nE,t_VEC);
  for(i=1;i<nE;i++)
  {
    gel(FqE1,i) = mkvec2(Z2Fq(gmael(E1,i,1),T),Z2Fq(gmael(E1,i,2),T));
    gel(FqE2,i) = mkvec2(Z2Fq(gmael(E2,i,1),T),Z2Fq(gmael(E2,i,2),T));
  }
  VHE1 = V_HE(d,abc1,FqE1,Z,T,p,e,pe); /* L(2D0-H1-E1) */
  VHE2 = V_HE(d,abc2,FqE2,Z,T,p,e,pe); /* L(2D0-H2-E2) */
	return gerepilecopy(av,mkvec2(VHE1,VHE2));
}

GEN PlaneEval0(GEN J, GEN W, GEN VHE)
{
  pari_sp av = avma;
  GEN T,p,pe,V,KV;
  long e;
  GEN S1,S2,s2,K;
  ulong d0,g,nV,i;

  JgetTpe(J,&T,&pe,&p,&e);
  d0 = Jgetd0(J);
  g = Jgetg(J);
  V = JgetV(J);
  KV = JgetKV(J);
  nV = lg(V);

  S1 = DivAdd(W,gel(VHE,1),2*d0+1,T,p,e,pe,0); /* L(4D0-D-H1-E1) */
  S1 = DivSub(V,S1,KV,1,T,p,e,pe,2); /* L(2D0-D-H1-E1) */
  S2 = DivMul(gel(S1,1),V,T,pe); /* L(4D0-D-H1-E1-ED) */
  S2 = DivSub(W,S2,KV,d0+1-g,T,p,e,pe,2); /* L(2D0-H1-E1-ED) */
  S2 = DivAdd(S2,gel(VHE,2),2*d0+1,T,p,e,pe,0); /* L(4D0-H1-E1-H2-E2-ED) */
  S2 = DivSub(V,S2,KV,1,T,p,e,pe,2); /* L(2D0-H1-E1-H2-E2-ED) */
  s2 = gel(S2,1);
  K = cgetg(nV+1,t_MAT);
  for(i=1;i<nV;i++) gel(K,i) = gel(V,i);
  gel(K,nV) = s2;
  K = matkerpadic(K,T,p,e);
	return gerepileupto(av,gel(K,1));
}

/*GEN PolExpId(GEN Z, GEN T, GEN pe)
{
	pari_sp av = avma;
	GEN f,a;
	ulong nZ,i;
	nZ = lg(Z);
	f = cgetg(nZ,t_VEC);
	for(i=1;i<nZ;i++) gel(f,i) = mkpoln(2,gen_1,gel(Z,i));
	f = liftpol(factorback(gmodulo(gmodulo(f,pe),T)));
	a = bestappr(f,NULL);
	return gerepilecopy(av,mkvecn(3,Z,f,a));
}

GEN AllPols0(GEN F, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN F1,f,R,pols;
	ulong nF,lF,npols,n,i,j,k;

	nF = lg(F);
	lF = lg(gel(F,1))-1;
	F1 = cgetg(lF,t_VEC);
	npols = 0;
	for(i=1;i<lF;i++)
	{ 
		npols++;
		gel(F1,i) = cgetg(nF,t_VEC);
		for(j=1;j<nF;j++)
		{
			f = gmael(F,j,i);
			if(ZX_is0mod(f,p))
			{
				gel(F1,i) = NULL;
				npols--;
				break;
			}
			gmael(F1,i,j) = ZpXQ_inv(f,T,p,e);
		}
	}
	npols *= (lF-2);
	pols = cgetg(npols+1,t_VEC);
	R = cgetg(nF,t_VEC);
	n = 0;
	for(i=1;i<lF;i++)
	{
		if(gel(F1,i)==NULL) continue;
		for(j=1;j<lF;j++)
		{
			if(j==i) continue;
			R = cgetg(nF,t_VEC);
			for(k=1;k<nF;k++) gel(R,k) = Fq_mul(gmael(F,k,j),gmael(F1,i,k),T,pe);
			n++;
			gel(pols,n) = PolExpId(R,T,pe);
		}
	}		
  
	return gerepilecopy(av,pols);
}*/
