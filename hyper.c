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
		x = random_FpX(dT,vT,p);
		y2 = poleval(f,x);
		y2 = FpXQ_red(y2,T,pe);
		if(gequal0(FpXQ_red(y2,T,p))) continue;
		y = ZpXQ_sqrt(y2,T,p,e);
		if(y == NULL) continue;
		P = mkvec2(x,y);
		return gerepilecopy(av,P);
	}
}

/* Matrix of values of x^i and x^i*y at the points in Ps */
/* x^i up to i=n, x^i*y up to n-d/2 where d=deg f */
/* This is L(n(infty_+ + infty_-)) */
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
	GEN pe,t,T,Frob,Z,Zp,P,Pp,Q,FrobCyc,x,y,W0,V,KV,V3,KV3,J;
  
	df = degree(f);
	/* TODO if(df%2) error0("Polynomial must be of even degree!"); */
	g = df/2-1;
	d0 = df; /* = 2g+2 */
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
		P = HyperRandPt(f,T,p,e,pe);
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

	W0 = RReval(Z,d0/2,df,T,pe);
	V = RReval(Z,d0,df,T,pe);
	KV = mateqnpadic(V,T,p,e);
	V3 = RReval(Z,3*d0/2,df,T,pe);
	KV3 = mateqnpadic(V3,T,p,e);

	J = mkvecn(lgJ,stoi(g),stoi(d0),T,p,stoi(e),pe,Frob,V,KV,W0,Z,FrobCyc,V3,KV3);
	return gerepilecopy(av,J);
}

GEN HyperPicRand(GEN J,GEN f) /* TODO not generic */
{
	pari_sp av1,av = avma;
	GEN T,p,pe,V;
	long e;
	ulong g,d0,df,i,j;
	GEN E[2];

	JgetTpe(J,&T,&pe,&p,&e);
	V = JgetV(J);
	d0 = Jgetd0(J);
	df = degree(f);
	g = Jgetg(J);

	/* Do twice : */
	for(j=0;j<2;j++)
	{
		av1 = avma;
		do
		{
			avma = av1;
			/* Take d0 random points */
			E[j] = cgetg(d0+1,t_VEC);
			for(i=1;i<=d0;i++)
			{
				gel(E[j],i) = HyperRandPt(f,T,p,e,pe);
			}
			/* Form the corresponding W */
			E[j] = RReval(E[j],d0,df,T,pe);
			E[j] = matkerpadic(E[j],T,p,e);
		} while(lg(E[j])!=d0+1-g+1); /* Check that the random points are independent */
		E[j] = FqM_mul(V,E[j],T,pe);
	}
	/* Return the chord of these two W */
	return gerepileupto(av,PicChord(J,E[0],E[1],0));
}

GEN ordJ(GEN f, GEN p, ulong a) /* Cardinal of Jac(y^2=f(x))(F_q), where q=p^a */
{
	pari_sp av = avma;
	GEN fp,chi,xa1,N;
	ulong i;

	fp = gmodulo(f,p);
	chi = hyperellcharpoly(fp);
	xa1 = cgetg(a+3,t_POL);
	for(i=1;i<a;i++) gel(xa1,i+2) = gen_0;
	gel(xa1,2) = gen_1;
	gel(xa1,a+2) = gen_m1;
	setvarn(xa1,varn(f));
	N = ZX_resultant(chi,xa1);

	return gerepileupto(av,N);
}

GEN HyperPicRandTors(GEN J, GEN f, GEN l, GEN C)
{
	pari_sp av1,av = avma;
	GEN T,p,fp,chi,chirem,xa1,chiC,N,M,W,lW,lv,fa;
	long v,a;
	ulong i;

	T = JgetT(J);
	p = Jgetp(J);
	a = degree(T);

	fp = gmodulo(f,p);
  chi = hyperellcharpoly(fp);
  xa1 = cgetg(a+3,t_POL);
  for(i=1;i<a;i++) gel(xa1,i+2) = gen_0;
  gel(xa1,2) = gen_1;
  gel(xa1,a+2) = gen_m1;
  setvarn(xa1,varn(f));
  N = ZX_resultant(chi,xa1);

	v = Z_pvalrem(N,l,&M); /* N = p^v*M */
	if(v==0) pari_err(e_MISC,"No rational %Ps-torsion",l);
	if(signe(C))
	{
		chiC = FpX_divrem(chi,C,l,&chirem);
		if(signe(chirem)) pari_err(e_MISC,"Incorrect characteristic polynomial");
		av1 = avma;
		if(degree(FpX_gcd(chiC,C,l))) pari_err(e_MISC,"Eigen multiplicity");
		lv = powis(l,v);
		fa = mkvec2(C,chiC);
		fa = polhensellift(chi,fa,l,v);
		chiC = gel(fa,2);
		chiC = FpX_center(chiC,lv,shifti(lv,-1));
		chiC = gerepileupto(av1,chiC);
	}
	else chiC = 0;

	do
	{
		W = HyperPicRand(J,f);
		W = PicMul(J,W,M,0);
		if(signe(C)) W = PicFrobPoly(J,W,chiC);
	} while(PicIsZero(J,W));

	lW = PicMul(J,W,l,0);
	av1 = avma;
	while(PicIsZero(J,lW)==0)
	{
		W = lW;
		av1 = avma;
		lW = PicMul(J,W,l,0);
	}
	
	avma = av1;
	return gerepileupto(av,W);
}

GEN HyperPicEvalData(GEN J)
{
	pari_sp av = avma;
	long e;
	GEN T,p,pe,Fq1,Z,col,U1,U2;
	ulong g,nZ,i,j;

	JgetTpe(J,&T,&pe,&p,&e);
  g = Jgetg(J);
	Z = JgetZ(J);
	nZ = lg(Z);
	Fq1 = GetFq1(T);

	if(g%2) pari_err(e_IMPL,"odd genus");
	U1 = RReval(Z,3*(g/2)+2,2*g+2,T,pe);
	U2 = cgetg(g+2,t_MAT);
	col = cgetg(nZ,t_COL);
	for(i=1;i<nZ;i++) gel(col,i) = Fq1;
	gel(U2,1) = col;
	col = cgetg(nZ,t_COL);
  for(i=1;i<nZ;i++) gel(col,i) = gmael(Z,i,2);
  gel(U2,g+1) = col;
  for(j=2;j<=g;j++)
	{
		col = cgetg(nZ,t_COL);
		for(i=1;i<nZ;i++) gel(col,i) = Fq_mul(gcoeff(U2,i,j-1),gmael(Z,i,1),T,pe);
    gel(U2,j) = col;
  }
	return gerepilecopy(av,mkvec2(U1,U2));
}

GEN HyperPicEval(GEN J, GEN W, GEN U)
{
	pari_sp av = avma;
	long e;
	ulong g,d0;
	GEN T,p,pe,V,KV,W0,U1,U2;
	ulong i;
	GEN WW0,S,s,sV,WE,K,res,u;
	
	JgetTpe(J,&T,&pe,&p,&e);
	g = Jgetg(J);
	d0 = Jgetd0(J);
	V = JgetV(J);
	KV = JgetKV(J);
	W0 = JgetW0(J);
	U1 = gel(U,1);
	U2 = gel(U,2);

	if(g%2) pari_err(e_IMPL,"odd genus");

	WW0 = DivAdd(W,W0,2*d0+1-g,T,p,e,pe,0);
	S = DivSub(U1,WW0,KV,1,T,p,e,pe,g+2); /* TODO true sub */
	sV = DivMul(gel(S,1),V,T,pe);
	WE = DivSub(W,sV,KV,g+3,T,p,e,pe,2);
	K = cgetg(g+1+g+3+1,t_MAT);
	for(i=1;i<=g+1;i++) gel(K,i) = gel(U2,i);
	for(i=1;i<=g+3;i++) gel(K,i+g+1) = gel(WE,i);
	K = matkerpadic(K,T,p,e);
	if(lg(K)>2) pari_err(e_MISC,"Genericity 2 failed");
	s = gel(K,1);
	u = gel(s,g+1);
	if(ZX_is0mod(u,p)) pari_err(e_MISC,"Genericity 3 failed");
	u = ZpXQ_inv(u,T,p,e);
	res = cgetg(g+1,t_VEC);
	for(i=1;i<=g;i++) gel(res,i) = gel(s,i);
	res = FqV_Fq_mul(res,u,T,pe);
	return gerepileupto(av,res);
}

