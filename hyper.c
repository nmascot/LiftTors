#include "pic.h"
#include "linalg.h"

GEN HyperRandPt(GEN f, GEN T, GEN p, ulong e, GEN pe)
{
	pari_sp av2,av = avma;
	long vT,dT;
	GEN x,y2,y,P;

	vT = varn(T);
	dT = degree(T);
	for(;;)
	{
		avma = av;
		x = random_FpX(dT-1,vT,p);
		av2 = avma;
		y2 = poleval(f,x);
		y2 = FpXQ_red(y2,T,pe);
		if(gequal0(FpXQ_red(y2,T,p))) continue;
		y = ZpXQ_sqrt(y2,T,p,e);
		if(y == NULL) continue;
    y = gerepileupto(av2,y);
		P = mkvec2(x,y);
		return P;
	}
}

/* Matrix of values of x^i and x^i*y and the points in Ps */
GEN RReval(GEN Ps, ulong n, ulong d, GEN T, GEN pe)
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
	pari_sp avP,av = avma;
	int newpt;
  ulong df,g,d0,nZ,n,ncyc,i;
	GEN pe,t,T,Frob,Z,Zp,P,Pp,Q,FrobCyc,x,y,W0,V,KV,J;
  
	df = degree(f);
	/* TODO if(df%2) error0("Polynomial must be of even degree!"); */
	g = df/2-1;
	d0 = df;
	nZ = 6*d0+1;

	t = varlower("t",varn(f));
  T = liftint(ffinit(p,a,varn(t)));
	Frob = ZpX_Frobenius(T,p,e);
	pe = powiu(p,e);
	
	n = ncyc = 0;
	Z = cgetg(nZ+a,t_VEC);
	Zp = cgetg(nZ+a,t_VEC);
	FrobCyc = cgetg(nZ,t_VECSMALL);
	while(n<nZ)
	{
		avP = avma;
		P = HyperRandPt(f,T,p,e,pe);
		/* Already have it ? */
		Pp = FpXV_red(P,p);
	  newpt = 1;
		for(i=1;i<=n;i++)
		{
			if(gequal(Pp,gel(Zp,i)))
			{
				newpt = 0;
				av = avP;
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
			x = gel(Q,1); y = gel(Q,2);
			x = FpXQ_red(poleval(x,Frob),T,pe);
			y = FpXQ_red(poleval(y,Frob),T,pe);
			Q = mkvec2(x,y);
		} while(!gequal(Q,P));
    FrobCyc[ncyc] = i;
	}
	setlg(Z,n+1);
	setlg(FrobCyc,ncyc+1);

	W0 = RReval(Z,d0,df,T,pe);
	V = RReval(Z,3*d0/2,df,T,pe);
	KV = shallowtrans(matkerpadic(shallowtrans(V),T,p,e));

	J = cgetg(lgJ,t_VEC);
	Jsetg(J,g);
	Jsetd0(J,d0);
	JsetT(J,T);
	Jsetp(J,p);
	Jsete(J,e);
	Jsetpe(J,pe);
	JsetFrob(J,Frob);
	JsetV(J,V);
	JsetKV(J,KV);
	JsetW0(J,W0);
	JsetZ(J,Z);
	JsetFrobCyc(J,FrobCyc);

	return gerepilecopy(av,J);
}
