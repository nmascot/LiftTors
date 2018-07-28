#include "pic.h"
#include "linalg.h"

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
	nV = d*(rsmax-rmin+1)+(d*(d-1))/2;
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
				if(r>0) xr = gel(x1pow,r);
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

GEN LinForm(GEN a, GEN b, GEN c, GEN Z, GEN T, GEN pe)
/* /!\ Not memory-clean */
{
	GEN L;
	ulong nZ,P;

	nZ = lg(Z);
	L = cgetg(nZ,t_COL);
	for(P=1;P<nZ;P++)
	{
		gel(L,P) = ZX_Z_mul(gmael(Z,P,1),a);
		gel(L,P) = ZX_add(gel(L,P),ZX_Z_mul(gmael(Z,P,2),b));
		gel(L,P) = ZX_Z_add(gel(L,P),c);
	}
	return L;
}

GEN sV_HE(long d, GEN a, GEN b, GEN c, GEN E, GEN Z, GEN T, GEN p, long e, GEN pe)
/* L(2D0-H-E) where H = hyper sect ax+by+c=0 and deg E=(g-2) -> dim 1 */
{
	pari_sp av = avma;
	GEN VH,EvE;

	VH = Plane_V(d,-2,2*(d-3)-1,Z,T,p,e,pe);
	VH = DivMul(LinForm(a,b,c,Z,T,pe),VH,T,pe);
	EvE = Plane_V(d,-2,2*(d-3)-1,E,T,p,e,pe);
	EvE = DivMul(LinForm(a,b,c,E,T,pe),EvE,T,pe);
	VH = FqM_mul(VH,matkerpadic(EvE,T,p,e),T,pe);
	return gerepileupto(av,VH);
}

GEN PlaneInit(GEN f, GEN p, ulong a, long e)
{
	pari_sp avP,av = avma;
  int newpt;
  ulong d1,d2,df,g,d0,nZ,n,ncyc,i;
  GEN vars,pe,t,T,Frob,Z,Zp,P,Pp,Q,FrobCyc,x,y,W0,V,KV,V3,KV3,J;

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
  Frob = ZpX_Frobenius(T,p,e);
  pe = powiu(p,e);

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
      x = FpX_FpXQ_eval(gel(Q,1),Frob,T,pe);
      y = FpX_FpXQ_eval(gel(Q,2),Frob,T,pe);
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

  J = mkvecn(lgJ,stoi(g),stoi(d0),T,p,stoi(e),pe,Frob,V,KV,W0,Z,FrobCyc,V3,KV3);
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