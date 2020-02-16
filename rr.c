#include "pic.h"
#include "linalg.h"
#include "freyruck.h"

/*GEN FnSubstMod(GEN F, long var, GEN val, GEN T, GEN pe)*/ /* /!\ Not memory-clean */
/*{
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
}*/

GEN EvalRatMod(GEN F, long var, GEN x, GEN T, GEN p, long e, GEN pe) /* /!\ Not memory-clean */
{
	GEN N,D;
	if(typ(F)==t_INT) return Z2Fq(F,T);
	if(typ(F)==t_FRAC) return Z2Fq(Fp_div(gel(F,1),gel(F,2),pe),T);
	if(varn(F)==varn(T)) return F;
	if(gvar(F)!=var) pari_err(e_MISC,"Bad var 1 in %Ps",F);
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
	if(typ(F)==t_FRAC) return Z2Fq(Fp_div(gel(F,1),gel(F,2),pe),T);
	if(typ(F)==t_RFRAC)
	{
		N = FnEvalAt(gel(F,1),P,vars,T,p,e,pe);
		if(typ(N)==t_INT) N=Z2Fq(N,T);
		D = FnEvalAt(gel(F,2),P,vars,T,p,e,pe);
		if(typ(D)==t_INT) D=Z2Fq(D,T);
		return gerepileupto(av,ZpXQ_div(N,D,T,pe,p,e));
	}
	if(gvar(F)==vars[2]) return liftall(poleval(F,gmodulo(gmodulo(gel(P,2),T),pe)));
	if(gvar(F)!=vars[1]) pari_err(e_MISC,"Bad var 2 in %Ps",F);
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
			if(DEBUGLEVEL) printf("Good, no relation\n");
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

GEN Fn1_FieldOfDef(GEN f, long var)
{
	GEN W,W1,c;
	long i,d;
	if(typ(f)==t_RFRAC)
	{
		W1 = Fn1_FieldOfDef(gel(f,1),var);
		W = Fn1_FieldOfDef(gel(f,2),var);
		if(W1 && W && !gequal(W1,W)) pari_err(e_MODULUS,"function",W1,W);
		return W1?W1:W;
	}

	if(typ(f)==t_POL)
	{
		W = NULL;
		d = poldegree(f,var);
		for(i=0;i<=d;i++)
		{
			c = polcoef(f,i,var);
			if(typ(c)==t_POLMOD)
			{
				W1=gel(c,1);
				if(W && !gequal(W1,W)) pari_err(e_MODULUS,"function",W1,W);
				W = W1;
			}
		}
		return W;
	}
	if(typ(f)==t_POLMOD) return gel(f,1);
	return NULL;
}

GEN Fn2_FieldOfDef(GEN f, GEN vars)
{
	GEN W,W1,c;
  long i,d;
  if(typ(f)==t_RFRAC)
  {
    W1 = Fn2_FieldOfDef(gel(f,1),vars);
    W = Fn2_FieldOfDef(gel(f,2),vars);
    if(W1 && W && !gequal(W1,W)) pari_err(e_MODULUS,"function",W1,W);
    return W1?W1:W;
  }

  if(typ(f)==t_POL)
  {
    W = NULL;
    d = poldegree(f,vars[1]);
    for(i=0;i<=d;i++)
    {
      c = polcoef(f,i,vars[1]);
			W1 = Fn1_FieldOfDef(c,vars[2]);
      if(W1)
      {
        if(W && !gequal(W1,W)) pari_err(e_MODULUS,"function",W1,W);
        W = W1;
      }
    }
    return W;
  }
  if(typ(f)==t_POLMOD) return gel(f,1);
  return NULL;
}

GEN RRspace_FieldOfDef(GEN L, GEN vars)
{
	GEN W,W1;
	ulong n,i;
	n = lg(L);
	W = NULL;
	for(i=1;i<n;i++)
	{
		W1 = Fn2_FieldOfDef(gel(L,i),vars);
		if(W1)
		{
			if(W && !gequal(W1,W)) pari_err(e_MODULUS,"function",W1,W);
      W = W1;
    }
	}
	return W;
}

GEN RRspaceEval(GEN L, GEN vars, GEN pts, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	long w,dW;
	GEN W,S,s,Li,res;
	ulong i;

	W = RRspace_FieldOfDef(L,vars);
	/* TODO delete
	 * L1 = gel(L,1);
	coeff1 = pollead(numerator(pollead(numerator(L1),vars[1])),vars[2]);
	ty = typ(coeff1);*/
	if(W) /* Algebraic case */
	{
		/* Get field of definition */
		w = gvar(W);
		dW = poldegree(W,w);
		/* Lift */
		L = liftpol(L);
		/* Find embeddings into Qp[t]/(T) */
		S = polrootsmod(W,mkvec2(T,p));
		if(lg(S) < dW+1) /* Should have all embeddings */
			pari_err(e_MISC,"Field of definition of Riemann-Roch space cannot be totally embedded into Q_q");
		res = cgetg(dW+1,t_VEC);
		for(i=1;i<=dW;i++)
		{
			s = gel(S,i); /* Embedding mod p */
			s = liftint(s);
			s = gadd(s,zeropadic(p,e));
			s = padicappr(W,s);
			s = gel(s,1);
			s = liftint(s);
			s = gmodulo(s,pe);
			Li = gsubst(L,w,s);
			Li = liftall(Li);
			/* TODO if there is a rescale, loss of p-adic accuracy? */
			gel(res,i) = FnsEvalAt_Rescale(Li,pts,vars,T,p,e,pe);
		}
		return gerepilecopy(av,res);
	}
	else /* Rational case */
  {
    avma = av;
    return mkvec(FnsEvalAt_Rescale(L,pts,vars,T,p,e,pe));
  }
}

GEN CurveLiftPty(GEN fx, GEN y, GEN T, GEN p, long e)
{
	pari_sp av = avma;
	y = gmodulo(liftall(y),T);
  y = gadd(y,zeropadic(p,e));
  y = padicappr(fx,y);
  y = gel(y,1);
  y = liftall(y);
	return gerepileupto(av,y);
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
		y = CurveLiftPty(fx,y,T,p,e);
    P = mkvec2(x,y);
    return gerepilecopy(av,P);
  }
}

GEN RREvalInit(GEN L, GEN vars, GEN Z, GEN T, GEN p, long e, GEN pe)
{
  pari_sp av = avma;
  GEN res;
  ulong i;
  res = cgetg(3,t_VEC);
  for(i=1;i<=2;i++)
    gel(res,i) = RRspaceEval(gel(L,i+2),vars,Z,T,p,e,pe);
  return gerepilecopy(av,res);
}

GEN RRInit(GEN f, ulong g, ulong d0, GEN L, GEN bad, GEN p, ulong a, long e)
{
	pari_sp avP,av = avma;
  int newpt;
  ulong nZ,n,cyc_top,i;
  GEN vars,pe,t,T,FrobMat,Z,Zp,P,Pp,Q,FrobCyc,x,y,V1,V2,V3,W0,V,KV,U,J;

	vars = variables_vecsmall(f);
	nZ = 5*d0+1;

	t = varlower("t",vars[2]);
  T = liftint(ffinit(p,a,varn(t)));
  pe = powiu(p,e);
  FrobMat = ZpXQ_FrobMat(T,p,e,pe);

	if(DEBUGLEVEL) printf("PicInit: Finding points\n");
  n = 0;
  Z = cgetg(nZ+a,t_VEC);
  Zp = cgetg(nZ+a,t_VEC);
  /* TODO sort Zp -> quasilin complexity */
  FrobCyc = cgetg(nZ+a,t_VECSMALL);
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
    Q = P;
    cyc_top = n+1;
    do
    {
			n++;
			FrobCyc[n] = n+1;
      gel(Z,n) = Q;
      gel(Zp,n) = FpXV_red(Q,p);
      x = Frob(gel(Q,1),FrobMat,T,pe);
      y = Frob(gel(Q,2),FrobMat,T,pe);
      Q = mkvec2(x,y);
    } while(!gequal(Q,P));
		FrobCyc[n] = cyc_top;
  }
  setlg(Z,n+1);
  setlg(FrobCyc,n+1);

	if(DEBUGLEVEL) printf("PicInit: Evaluating rational functions\n");
	V1 = FnsEvalAt_Rescale(gel(L,1),Z,vars,T,p,e,pe);
	V2 = FnsEvalAt_Rescale(gel(L,2),Z,vars,T,p,e,pe);
	V3 = DivAdd(V1,V2,3*d0+1-g,T,p,e,pe,0);
	W0 = V1;
	V = mkvecn(3,V1,V2,V3);
	if(DEBUGLEVEL) printf("PicInit: Computing equation matrices\n");
	KV = cgetg(4,t_VEC);
	for(i=1;i<=3;i++) gel(KV,i) = mateqnpadic(gel(V,i),T,p,e);
	if(DEBUGLEVEL) printf("PicInit: Constructing evaluation maps\n");
	U = RREvalInit(L,vars,Z,T,p,e,pe);
  J = mkvecn(lgJ,f,stoi(g),stoi(d0),L,T,p,stoi(e),pe,FrobMat,V,KV,W0,U,Z,FrobCyc);
	return gerepilecopy(av,J);
}

GEN Jlift(GEN J, ulong e2)
{
	pari_sp av = avma, avZ;
	GEN J2,T,p,pe2,f,vars,L,FrobCyc,FrobMat2;
	long g,d0;
	GEN Z,Z2,V1,V2,V3,W0,V,KV,KV3,P,x,y,fx,U;
	ulong nZ,i,j;
  if(Jgete(J)>=e2)
	{
		pari_warn(warner,"Current accuracy already higher than required in Jlift, not changing anything");
		return gcopy(J);
	}
	f = Jgetf(J);
	if(gequal0(f))
		pari_err(e_MISC,"Cannot increase accuracy for this curve (missing equation)");
	L = JgetL(J);
	if(lg(L)!=5)
    pari_err(e_MISC,"Cannot increase accuracy for this curve (missing RR spaces)");
	Z = JgetZ(J);
	if(lg(Z)==1)
    pari_err(e_MISC,"Cannot increase accuracy for this curve (missing points)");
	T = JgetT(J);
	p = Jgetp(J);
	vars = variables_vecsmall(f);
	g = Jgetg(J);
	d0 = Jgetd0(J);

  J2 = cgetg(lgJ+1,t_VEC);
  gel(J2,1) = f;
  gel(J2,2) = stoi(g);
  gel(J2,3) = stoi(d0);
  gel(J2,4) = L;
  gel(J2,5) = T;
  gel(J2,6) = p;
  gel(J2,7) = utoi(e2);
  gel(J2,8) = pe2 = powiu(p,e2);
  gel(J2,9) = FrobMat2 = ZpXQ_FrobMat(T,p,e2,pe2);
	
	nZ = lg(Z);
	FrobCyc = JgetFrobCyc(J);
	avZ = avma;
	Z2 = cgetg(nZ,t_VEC);
	/* Need to lift while respecting Frob orbits */
	for(i=1;i<nZ;i++) /* Mark points as not done */
		gel(Z2,i) = NULL;
	for(i=1;i<nZ;i++)
	{
		if(gel(Z2,i)) continue; /* This point is part of a Frob orbit we have already treated */
		j = i; /* Start of a new orbit */
		P = gel(Z,i); /* Arbitrarily lift this point only */
    x = gel(P,1);
    y = gel(P,2);
    fx = poleval(f,x);
    y = CurveLiftPty(fx,y,T,p,e2);
		while(1)
		{
			gel(Z2,j) = mkvec2(x,y);
			j = FrobCyc[j];
			if(j==i) break; /* End of orbit */
			x = Frob(x,FrobMat2,T,pe2);
      y = Frob(y,FrobMat2,T,pe2);
		}
	}
	Z2 = gerepilecopy(avZ,Z2);
	V1 = FnsEvalAt_Rescale(gel(L,1),Z2,vars,T,p,e2,pe2);
  V2 = FnsEvalAt_Rescale(gel(L,2),Z2,vars,T,p,e2,pe2);
  V3 = DivAdd(V1,V2,3*d0+1-g,T,p,e2,pe2,0);
  W0 = V1; /* TODO can it happen that W0 != V1 even though all data is present? */
  V = V2;
  KV = mateqnpadic(V,T,p,e2);
  KV3 = mateqnpadic(V3,T,p,e2);
	U = cgetg(3,t_VEC);
  for(i=1;i<=2;i++)
    gel(U,i) = RRspaceEval(gel(L,i+2),vars,Z2,T,p,e2,pe2);
  gel(J2,10) = V1;
  gel(J2,11) = V;
  gel(J2,12) = V3;
  gel(J2,13) = KV;
  gel(J2,14) = KV3;
  gel(J2,15) = W0;
	gel(J2,16) = U;
  gel(J2,17) = Z2;
  gel(J2,18) = JgetFrobCyc(J);
  return gerepilecopy(av,J2);
}

GEN RREval(GEN J, GEN W)
{
	pari_sp av = avma, av1, av2;
	GEN T,p,pe,V,KV,U,res;
	long e;
	ulong n1,n2,i1,i2;
	GEN S1,S2,s2,K;
	ulong d0,g,nV,i;
	
	JgetTpe(J,&T,&pe,&p,&e);
	d0 = Jgetd0(J);
	g = Jgetg(J);
	V = JgetV(J,2);
	nV = lg(V);
	KV = JgetKV(J,2);
	U = JgetEvalData(J); /* L(2D0-Ei), deg Ei = d0-g (i=1,2), repeated for each embedding into Qq */
	n1 = lg(gel(U,1)); /* Deg of E1 / Q */
	n2 = lg(gel(U,2)); /* Deg of E2 / Q */

	res = cgetg(n1,t_MAT);
	for(i1=1;i1<n1;i1++)
	{
		gel(res,i1) = cgetg(n2,t_COL);
		av1 = avma;
		S1 = DivAdd(W,gmael(U,1,i1),2*d0+1,T,p,e,pe,0); /* L(4D0-D-E1) */
		S1 = DivSub(V,S1,KV,1,T,p,e,pe,2); /* L(2D0-D-E1), generically 1-dimensional */
		S1 = gel(S1,1); /* Generator */
		S1 = gerepileupto(av1,S1);
		for(i2=1;i2<n2;i2++)
		{
			av2 = avma;
			S2 = DivMul(S1,V,T,pe); /* L(4D0-D-E1-ED) */
			S2 = DivSub(W,S2,KV,d0+1-g,T,p,e,pe,2); /* L(2D0-E1-ED) */
			S2 = DivAdd(S2,gmael(U,2,i2),2*d0+1,T,p,e,pe,0); /* L(4D0-E1-E2-ED) */
			S2 = gerepileupto(av2,S2);
			S2 = DivSub(V,S2,KV,1,T,p,e,pe,2); /* L(2D0-E1-E2-ED), generically 1-dimensional */
  		s2 = gel(S2,1); /* Generator */
			s2 = gerepileupto(av2,s2);
			/* get coords of s2 w.r.t. V */
  		K = cgetg(nV+1,t_MAT);
  		for(i=1;i<nV;i++) gel(K,i) = gel(V,i);
  		gel(K,nV) = s2;
  		K = matkerpadic(K,T,p,e);
			K = gel(K,1);
			setlg(K,nV);
			gcoeff(res,i2,i1)=gerepilecopy(av2,K);
		}
	}
	return gerepilecopy(av,res);
}

GEN PolExpId(GEN Z, GEN T, GEN pe) /* bestappr of prod(x-z), z in Z */
{
	pari_sp av = avma;
	GEN f,a;
	f = FqV_roots_to_pol(Z,T,pe,0);
	if(poldegree(f,varn(T))>0) pari_err(e_MISC,"Irrational coefficient: %Ps",f);
	f = simplify_shallow(f);
	f = gmodulo(f,pe);
	a = bestappr(f,NULL);
	return gerepilecopy(av,mkvecn(3,Z,f,a));
}

GEN OnePol(GEN N, GEN D, GEN T, GEN pe)
{ /* Actually returns a vector of n1*n2 pols (all elem. symm. fns) */
	pari_sp av = avma;
	GEN R,Z,F;
	ulong k,n,i1,i2,i;
	long n1,n2,n12;
	n = lg(N);
	RgM_dimensions(gel(N,1),&n2,&n1);
	n12 = n1*n2;
	R = cgetg(n,t_VEC);
	Z = cgetg(n12+1,t_VEC);
  for(k=1;k<n;k++)
	{
		i=1;
		for(i1=1;i1<=n1;i1++)
		{
			for(i2=1;i2<=n2;i2++)
			{
				gel(Z,i) = Fq_mul(gmael3(N,k,i1,i2),gmael3(D,k,i1,i2),T,pe);
				i++;
			}
		}
		gel(R,k) = FqV_roots_to_pol(Z,T,pe,0);
	}
	F = cgetg(n12+1,t_VEC);
	Z = cgetg(n,t_VEC);
	for(i=0;i<n12;i++)
	{
		for(k=1;k<n;k++)
			gel(Z,k) = polcoef(gel(R,k),i,0);
		gel(F,i+1) = PolExpId(Z,T,pe);
	}
	return gerepileupto(av,F);
}

GEN AllPols(GEN F, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN Ft,F1,f,pols;
	ulong nF,lF,npols,n,i,j,i1,i2,m,k;
	long n1,n2,n12;
	struct pari_mt pt;
	GEN worker,done;
	long pending,workid;

	nF = lg(F); /* Number of vectors */
	RgM_dimensions(gel(F,1),&n2,&n1);
	n12 = n1*n2;
	lF = lg(gmael3(F,1,1,1))-1; /* Size of each vector */
	/* F = list of nF-1 matrices of size n2*n1 of vectors of size lF */
	Ft = cgetg(lF,t_VEC);
	/* Ft[j,i,i1,i2] = F[i,i1,i2,j] */
	for(j=1;j<lF;j++)
	{
		gel(Ft,j) = cgetg(nF,t_VEC);
		for(i=1;i<nF;i++)
		{
			gmael(Ft,j,i) = cgetg(n1+1,t_MAT);
			for(i1=1;i1<=n1;i1++)
			{
				gmael3(Ft,j,i,i1) = cgetg(n2+1,t_COL);
				for(i2=1;i2<=n2;i2++)
					gmael4(Ft,j,i,i1,i2) = gmael4(F,i,i1,i2,j);
			}
		}
	}
	F1 = cgetg(lF,t_VEC);
	npols = 0;
	for(i=1;i<lF;i++) /* Find the i such that the ith coord of all the vectors in all the matrices are invertible, and store there inverses */
	{
		npols++;
		gel(F1,i) = cgetg(nF,t_VEC);
		for(j=1;j<nF;j++)
		{
			gmael(F1,i,j) = cgetg(n1+1,t_MAT);
			for(i1=1;i1<=n1;i1++)
			{
				gmael3(F1,i,j,i1) = cgetg(n2+1,t_COL);
				for(i2=1;i2<=n2;i2++)
				{
					f = gmael4(F,j,i1,i2,i);
					if(ZX_is0mod(f,p))
					{
						gel(F1,i) = NULL;
						npols--;
						i2=n2+1;i1=n1+1;j=nF;
					}
					else
						gmael4(F1,i,j,i1,i2) = ZpXQ_inv(f,T,p,e); /* TODO do it later in case we give up this i */
				}
			}
		}
	}
	npols *= (lF-2)*n12;
	pols = cgetg(npols+1,t_VEC);
	pending = 0;
  worker = strtofunction("OnePol");
  mt_queue_start(&pt,worker);
	done = NULL;
	for(i=j=m=n=1;i<lF||pending;n++,j++)
	{
		if(j==lg(F1))
		{
			j=1;
			i++;
		}
		if(gel(F1,j)==NULL || i==j) continue; /* Skip if denom=0 or if denom=num */
		mt_queue_submit(&pt,n,i<lF?mkvecn(4,gel(Ft,i),gel(F1,j),T,pe):NULL);
		done = mt_queue_get(&pt,&workid,&pending);
		if(done)
		{
			for(k=1;k<=n12;k++)
			{
				gel(pols,m) = gel(done,k);
				m++;
			}
		}
	}		
  mt_queue_end(&pt);
	return gerepilecopy(av,pols);
}
