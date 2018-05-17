#include "pic.h"
#include "linalg.h"

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
		y2 = FpXQ_red(y2,T,pe);
		if(gequal0(FpXQ_red(y2,T,p))) continue;
		y = ZpXQ_sqrt(y2,T,p,e);
		if(y == NULL) continue;
		P = mkvec2(x,y);
		return gerepilecopy(av,P);
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

GEN HyperInit(GEN f, GEN p, ulong a, long e)
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
	printf("A");

	t = varlower("t",varn(f));
  T = liftint(ffinit(p,a,varn(t)));
	Frob = ZpX_Frobenius(T,p,e);
	pe = powiu(p,e);
	printf("B");

	
	n = ncyc = 0;
	Z = cgetg(nZ+a,t_VEC);
	Zp = cgetg(nZ+a,t_VEC);
	/* TODO sort Zp -> quasilin complexity */
	FrobCyc = cgetg(nZ+1,t_VECSMALL);
	while(n<nZ)
	{
		printf("Loop, n=%lu\n",n);
		avP = avma;
		P = HyperRandPt(f,T,p,e,pe);
		pari_printf("%Ps\n",P);
		/* Already have it ? */
		Pp = FpXV_red(P,p);
	  newpt = 1;
		for(i=1;i<=n;i++)
		{
			if(gequal(Pp,gel(Zp,i)))
			{
				newpt = 0;
				printf("Notnew");
				avma = avP;
				break;
			}
		}
		printf("End verif\n");
		if(newpt == 0) continue;
		ncyc++;
		Q = P;
		i = 0;
		printf("Starting cycle %lu\n",ncyc);
		do
		{
			i++;
			n++;
			printf("Point %lu, step %lu\n",n,i);
			pari_printf("%Ps\n",Q);
			gel(Z,n) = Q;
			gel(Zp,n) = FpXV_red(Q,p);
			printf("Now apply Frob\n");
			x = FpX_FpXQ_eval(gel(Q,1),Frob,T,pe);
			y = FpX_FpXQ_eval(gel(Q,2),Frob,T,pe);
			Q = mkvec2(x,y);
		} while(!gequal(Q,P));
    FrobCyc[ncyc] = i;
		printf("Cycle %lu has length %lu\n",ncyc,i);
	}
	setlg(Z,n+1);
	setlg(FrobCyc,ncyc+1);
	printf("D");

	W0 = RReval(Z,d0,df,T,pe);
	V = RReval(Z,3*d0/2,df,T,pe);
	KV = shallowtrans(matkerpadic(shallowtrans(V),T,p,e));
	printf("E");

	J = mkvecn(lgJ,stoi(g),stoi(d0),T,p,stoi(e),pe,Frob,V,KV,W0,Z,FrobCyc);
	printf("F");

	return gerepilecopy(av,J);
}
