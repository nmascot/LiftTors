#include "pic.h"
#include "linalg.h"
#include "freyruck.h"

GEN FnSubstMod(GEN F, long var, GEN val, GEN T, GEN pe) /* /!\ Not memory-clean */
{
	GEN valm,res;
	pari_CATCH(e_INV)
	{
		printf("Z");
		res = gsubst(F,var,val);
		return res;
		printf("u");
		res = gmodulo(gmodulo(res,T),pe);
		printf("v");
		res = liftall(res);
		printf("w");
		return res;
	}
	pari_TRY
	{
		valm = gmodulo(gmodulo(val,T),pe);
		res = gsubst(F,var,valm);
		res = liftall(res);
	}
	pari_ENDCATCH
	return res;
}

GEN EvalRatMod(GEN F, long var, GEN x, GEN T, GEN p, long e, GEN pe) /* /!\ Not memory-clean */
{
	GEN N,D;
	if(typ(F)==t_INT) return Z2Fq(F,T);
	if(gvar(F)!=var) pari_err(e_MISC,"Bad var 1");
	if(typ(F)==t_POL)
	{
		N = liftall(poleval(F,gmodulo(gmodulo(x,T),pe)));
		if(typ(N)==t_INT) N=Z2Fq(N,T);
		return N;
	}
	N = liftall(poleval(gel(F,1),gmodulo(gmodulo(x,T),pe)));
	if(typ(N)==t_INT) N=Z2Fq(N,T);
	D = liftall(poleval(gel(F,2),gmodulo(gmodulo(x,T),pe)));
	if(typ(D)==t_INT) D=Z2Fq(D,T);
	N = ZpXQ_div(N,D,T,pe,p,e);
	return N;
}

GEN FnEvalAt(GEN F, GEN P, GEN vars, GEN T, GEN p, long e/*GEN E*/, GEN pe)
/* F=N/D, N,D=R(y)x^n+...+R(y), R(y) rat fracs. Assumes P=(a,b) is s.t. the denom of R(b) is nonzero mod p for all R. */
{
	pari_sp av = avma;
	GEN N,D,Fy;
	long /*e = itos(E),*/d;
	ulong i;
	if(typ(F)==t_INT) return Z2Fq(F,T);
	if(typ(F)==t_RFRAC)
	{
		N = FnEvalAt(gel(F,1),P,vars,T,p,e,pe);
		if(typ(N)==t_INT) N=Z2Fq(N,T);
		D = FnEvalAt(gel(F,2),P,vars,T,p,e,pe);
		if(typ(D)==t_INT) D=Z2Fq(D,T);
		return gerepileupto(av,ZpXQ_div(N,D,T,pe,p,e));
	}
	if(gvar(F)==vars[2]) return liftall(poleval(F,gmodulo(gmodulo(gel(P,2),T),pe)));
	if(gvar(F)!=vars[1]) pari_err(e_MISC,"Bad var 2");
	d = lg(F);
	Fy = cgetg(d,t_POL);
	setsigne(Fy,1);
	setvarn(Fy,vars[1]);
	for(i=2;i<d;i++) gel(Fy,i) = EvalRatMod(gel(F,i),vars[2],gel(P,2),T,p,e,pe);
	F = liftall(poleval(Fy,gmodulo(gmodulo(gel(P,1),T),pe)));
	return gerepilecopy(av,F);
}

GEN FnsEvalAt(GEN Fns, GEN Z, GEN vars, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN A;//E;
	ulong nF,nZ;
	long i,j;//k;
	/*struct pari_mt pt;
	GEN worker,done;
	long pending,workid;*/

	/*E = stoi(e);*/
	nF = lg(Fns);
	nZ = lg(Z);
	A = cgetg(nF,t_MAT);
	for(j=1;j<nF;j++)
	{
		gel(A,j) = cgetg(nZ,t_COL);
		for(i=1;i<nZ;i++)
  	{
			//printf("%ld,%ld\n",i,j);
      gcoeff(A,i,j) = FnEvalAt(gel(Fns,j),gel(Z,i),vars,T,p,e,pe);
  	}
	}
	/* Abandoned parallel version (not useful)
	nF--;nZ--;
	pending = 0;
  worker = strtofunction("FnEvalAt");
  mt_queue_start(&pt,worker);
	for(k=0;k<nF*nZ||pending;k++)
	{
		if(k<nF*nZ)
		{
			i = 1 + (k%nZ);
			j = 1 + (k/nZ);
			printf("%ld,%ld\n",i,j);
			mt_queue_submit(&pt,k,mkvecn(7,gel(Fns,j),gel(Z,i),vars,T,p,E,pe));
		}
		else mt_queue_submit(&pt,k,NULL);
		done = mt_queue_get(&pt,&workid,&pending);
		if(done)
		{
			i = 1 + (k%nZ);
      j = 1 + (k/nZ);
			gcoeff(A,i,j) = done;
		}
	}*/
	return gerepilecopy(av,A);
}

GEN FnsEvalAt_Rescale(GEN Fns, GEN Z, GEN vars, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN F,S,K,f,redo,rF;
	ulong i,j,k,nF,nK;
	F = gcopy(Fns);
	nF = lg(F);
	S = FnsEvalAt(F,Z,vars,T,p,e,pe);
	while(1)
	{
		K = FqM_ker(S,T,p);
		nK = lg(K);
		/* Are the evals (and hence the fns) independent ? */
		if(nK==1)
		{
			printf("Good, no relation\n");
			return gerepileupto(av,S);
		}
		pari_printf("Found %ld relations, eliminating and re-evaluating\n",nK-1);
		/* No. We assume Z def / Q, so K has entries in Fp */
		/* Do elimination and start over */
		redo = cgetg(nK,t_VECSMALL);
		rF = cgetg(nK,t_VEC);
		for(j=1;j<nK;j++)
		{
			/* k = pivot = last nonzero entry of the col (it's a 1) */
			k = 0;
			for(i=1;i<nF;i++)
			{
				if(!gequal0(gcoeff(K,i,j))) k=i;
			}
			redo[j]=k;
			/* Form corresponding lin comb, and div by p */
			f = gel(F,k);
			for(i=1;i<k;i++)
			{
				if(!gequal0(gcoeff(K,i,j)))
				{
					f = gadd(f,gmul(centerlift(gmodulo(gcoeff(K,i,j),p)),gel(F,i)));
				}
			}
			gel(F,k) = gel(rF,j) = gdiv(f,p);
		}
		rF = FnsEvalAt(rF,Z,vars,T,p,e,pe);
		for(j=1;j<nK;j++) gel(S,redo[j]) = gel(rF,j);
	}
}

GEN CurveRandPt(GEN f, GEN T, GEN p, long e, GEN bad)
{
	pari_sp av = avma, av1;
	long vT,dT;
  GEN vars,x,fx,y,badpt,dfx,dy,P;
	vT = varn(T);
  dT = degree(T);
	vars = variables_vecsmall(f);
	av1 = avma;
  for(;;)
  {
    avma = av1;
    x = random_FpX(dT,vT,p);
		if(ZX_is0mod(x,p)) continue; /* Want x != 0 */
		fx = poleval(f,x);
    y = polrootsmod(fx,mkvec2(T,p));
		if(lg(y)==1) continue; /* No roots */
		y = gel(y,itos(genrand(stoi(lg(y)-1)))+1);
		badpt = FnEvalAt(bad,mkvec2(x,liftall(y)),vars,T,p,1,p);
		badpt = Fq_red(badpt,T,p);
		if(gequal0(badpt)) continue; /* Forbidden locus */
		dfx = RgX_deriv(fx);
		dy = poleval(dfx,y);
		if(gequal0(dy)) continue; /* Bad for Hensel */
		/* TODO check if bad */
		y = gmodulo(liftall(y),T);
		y = gadd(y,zeropadic(p,e));
		y = padicappr(fx,y);
		y = gel(y,1);
		y = liftall(y);
    P = mkvec2(x,y);
    return gerepilecopy(av,P);
  }
}

GEN RRInit(GEN f, ulong g, ulong d0, GEN L, GEN L2, GEN bad, GEN p, ulong a, long e)
{
	pari_sp avP,av = avma;
  int newpt;
  ulong nZ,n,ncyc,i;
  GEN vars,pe,t,T,FrobMat,Z,Zp,P,Pp,Q,FrobCyc,x,y,V1,V2,V3,W0,V,KV,KV3,J;

	vars = variables_vecsmall(f);
	nZ = 5*d0+1;

	t = varlower("t",vars[2]);
  T = liftint(ffinit(p,a,varn(t)));
  pe = powiu(p,e);
  FrobMat = ZpXQ_FrobMat(T,p,e,pe);

	printf("PicInit: Finding points\n");
  n = ncyc = 0;
  Z = cgetg(nZ+a,t_VEC);
  Zp = cgetg(nZ+a,t_VEC);
  /* TODO sort Zp -> quasilin complexity */
  FrobCyc = cgetg(nZ+1,t_VECSMALL);
  while(n<nZ)
  {
    avP = avma;
    P = CurveRandPt(f,T,p,e,bad);
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

	printf("PicInit: Evaluating rational functions\n");
	V1 = FnsEvalAt_Rescale(L,Z,vars,T,p,e,pe);
	V2 = FnsEvalAt_Rescale(L2,Z,vars,T,p,e,pe);
	V3 = DivAdd1(V1,V2,3*d0+1-g,T,p,e,pe,0);
	W0 = V1;
	V = V2;
  KV = mateqnpadic(V,T,p,e);
  KV3 = mateqnpadic(V3,T,p,e);

  J = mkvecn(lgJ,f,stoi(g),stoi(d0),T,p,stoi(e),pe,FrobMat,V,KV,W0,Z,FrobCyc,V3,KV3);
	return gerepilecopy(av,J);
}

GEN RREvalInit(GEN J, GEN Li)
{
	pari_sp av = avma;
	GEN Z,T,p,pe,vars,res;
	long e;
	ulong i;
	Z = JgetZ(J);
	JgetTpe(J,&T,&pe,&p,&e);
	vars = variables_vecsmall(Jgetf(J));
	res = cgetg(3,t_VEC);
	for(i=1;i<=2;i++)
	{
		gel(res,i) = FnsEvalAt_Rescale(gel(Li,i),Z,vars,T,p,e,pe);
	}
	return gerepilecopy(av,res);
}

GEN RREval(GEN J, GEN W, GEN Li) /* Li = L(2D0-Ei), deg Ei = d0-g (i=1,2) */
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
	nV = lg(V);
	KV = JgetKV(J);

	S1 = DivAdd1(W,gel(Li,1),2*d0+1,T,p,e,pe,0); /* L(4D0-D-E1) */
	S1 = DivSub(V,S1,KV,1,T,p,e,pe,2); /* L(2D0-D-E1) */
	S2 = DivMul(gel(S1,1),V,T,pe); /* L(4D0-D-E1-ED) */
	S2 = DivSub(W,S2,KV,d0+1-g,T,p,e,pe,2); /* L(2D0-E1-ED) */
	S2 = DivAdd1(S2,gel(Li,2),2*d0+1,T,p,e,pe,0); /* L(4D0-E1-E2-ED) */
	S2 = DivSub(V,S2,KV,1,T,p,e,pe,2); /* L(2D0-E1-E2-ED) */
  s2 = gel(S2,1);
	s2 = gerepilecopy(av,s2);
  K = cgetg(nV+1,t_MAT);
  for(i=1;i<nV;i++) gel(K,i) = gel(V,i);
  gel(K,nV) = s2;
  K = matkerpadic(K,T,p,e);
	K = gel(K,1);
	setlg(K,nV);
	return gerepilecopy(av,K);
}

GEN PolExpId(GEN Z, GEN T, GEN pe) /* bestappr of prod(x-z), z in Z */
{
	pari_sp av = avma;
	GEN f,c,a;
	ulong nZ,i;
	nZ = lg(Z);
	f = cgetg(nZ,t_VEC);
	for(i=1;i<nZ;i++) gel(f,i) = mkpoln(2,gen_1,gmodulo(gmodulo(gel(Z,i),pe),T));
	f = liftpol(factorback(f));
	for(i=2;i<lg(f);i++)
	{
		c = gel(f,i);
		if(degpol(c)>0) pari_err(e_MISC,"Irrational coefficient: %Ps k=%lu",c,i-2);
		if(degpol(c)==-1) c = gen_0;
		else c = gel(c,2);
		gel(f,i) = c;
	}
	a = bestappr(f,NULL);
	return gerepilecopy(av,mkvecn(3,Z,f,a));
}

GEN OnePol(GEN N, GEN D, GEN T, GEN pe)
{
	pari_sp av = avma;
	GEN R,F;
	ulong k,n;
	n = lg(N);
	R = cgetg(n,t_VEC);
  for(k=1;k<n;k++) gel(R,k) = Fq_mul(gel(N,k),gel(D,k),T,pe);
	F = PolExpId(R,T,pe);
	return gerepileupto(av,F);
}

GEN AllPols(GEN F, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN Ft,F1,f,pols;
	ulong nF,lF,npols,n,i,j,m;
	struct pari_mt pt;
	GEN worker,done;
	long pending,workid;

	nF = lg(F); /* Number of vectors */
	lF = lg(gel(F,1))-1; /* Size of each vector */
	Ft = cgetg(lF,t_VEC);
	for(i=1;i<lF;i++)
	{ 
		gel(Ft,i) = cgetg(nF,t_VEC);
		for(j=1;j<nF;j++) gmael(Ft,i,j) = gmael(F,j,i);
	}
	F1 = cgetg(lF,t_VEC);
	npols = 0;
	for(i=1;i<lF;i++) /* Find the i such that the ith coord of the vectors are all invertible */
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
	pending = 0;
  worker = strtofunction("OnePol");
  mt_queue_start(&pt,worker);
	i=m=1;
	j=0;n=0;
	done = NULL;
	for(i=j=m=n=1;i<lF||pending;n++,j++)
	{
		if(j==lg(F1))
		{
			j=1;
			i++;
		}
		if(gel(F1,j)==NULL || i==j) continue;
		mt_queue_submit(&pt,n,i<lF?mkvecn(4,gel(Ft,i),gel(F1,j),T,pe):NULL);
		done = mt_queue_get(&pt,&workid,&pending);
		if(done)
		{
			//printf("Getting poly number %ld\n",m);
			gel(pols,m) = done;
			m++;
		}
	}		
  mt_queue_end(&pt);
	return gerepilecopy(av,pols);
}
