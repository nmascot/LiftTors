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

