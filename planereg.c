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

GEN Vmat_upto(ulong d, ulong n, GEN Z, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN V,J,Fq1,x,y,x1,xpow,x1pow,ypow,xr,ys,xrys;
	long r,s;
	ulong g,nZ,P,m,dm,j;

	/* Initialisation */
	g = (d-1)*(d-2);
	g = g/2;
	nZ = lg(Z);
	V = cgetg(n+1,t_VEC);
	J = cgetg(n+1,t_VECSMALL);
	for(m=1;m<=n;m++)
	{
		J[m] = 1;
		dm = m*d*(d-2);
		gel(V,m) = cgetg(dm+1,t_MAT);
		for(j=1;j<=dm;j++) gmael(V,m,j) = cgetg(nZ,t_COL);
	}
	Fq1 = GetFq1(T);
	xpow = cgetg(n*(d-3)+1,t_VEC);
	x1pow = cgetg(n+1,t_VEC);
	ypow = cgetg(d,t_VEC);

	for(P=1;P<=nZ;P++)
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
		
		for(s=0;s<d;s++)
		{
			for(r=-n;r+s<=(d-3)*n;r++)
			{
				if(r==0) xr = Fq1;
				if(r>0) xr = gel(xpow,r);
				else xr = gel(x1pow,-r);
				if(s)
				{
					ys = gel(y,s);
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
	d0 = df*(df-3);
	nZ = 5*d0+1;

	t = varlower("t",varn(f));
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
