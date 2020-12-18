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
	if(typ(F)==t_INT) return scalarpol(F,varn(T));
	if(typ(F)==t_FRAC) return scalarpol(Fp_div(gel(F,1),gel(F,2),pe),varn(T));
	if(varn(F)==varn(T)) return F;
	if(gvar(F)!=var) pari_err(e_MISC,"Bad var 1 in %Ps",F);
	if(typ(F)==t_POL)
	{
		N = liftall(poleval(F,gmodulo(gmodulo(x,T),pe)));
		if(typ(N)==t_INT) N=scalarpol(N,varn(T));
		return N;
	}
	N = liftall(poleval(gel(F,1),gmodulo(gmodulo(x,T),pe)));
	if(typ(N)==t_INT) N=scalarpol(N,varn(T));
	D = liftall(poleval(gel(F,2),gmodulo(gmodulo(x,T),pe)));
	if(typ(D)==t_INT) D=scalarpol(D,varn(T));
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
	if(typ(F)==t_INT) return scalarpol(F,varn(T));
	if(typ(F)==t_FRAC) return scalarpol(Fp_div(gel(F,1),gel(F,2),pe),varn(T));
	if(typ(F)==t_RFRAC)
	{
		N = FnEvalAt(gel(F,1),P,vars,T,p,e,pe);
		if(typ(N)==t_INT) N=scalarpol(N,varn(T));
		D = FnEvalAt(gel(F,2),P,vars,T,p,e,pe);
		if(typ(D)==t_INT) D=scalarpol(D,varn(T));
		return gerepileupto(av,ZpXQ_div(N,D,T,pe,p,e));
	}
	if(gvar(F)==vars[2]) return liftall(poleval(F,gmodulo(gmodulo(gel(P,2),T),pe)));
	if(gvar(F)!=vars[1]) pari_err(e_MISC,"Bad var 2 in %Ps",F);
	d = lg(F);
	Fy = cgetg(d,t_POL);
	Fy[1] = 0;
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

GEN RREvalInit(GEN L, GEN vars, GEN Z, GEN V2, GEN T, GEN p, long e, GEN pe)
{
  pari_sp av = avma;
  GEN res,I,M;
  ulong i,j,nV2;
  res = cgetg(5,t_VEC);
  for(i=1;i<=2;i++) /* TODO parallelise */
    gel(res,i) = RRspaceEval(gel(L,i+2),vars,Z,T,p,e,pe);
	nV2 = lg(V2);
	I = FqM_indexrank(V2,T,p);
	I = gel(I,1); /* Rows of V2 forming invertible block */
	gel(res,3) = I;
	/* That invertible block */
	M = cgetg(nV2,t_MAT);
	for(j=1;j<nV2;j++)
	{
		gel(M,j) = cgetg(nV2,t_COL);
		for(i=1;i<nV2;i++)
			gcoeff(M,i,j) = gcoeff(V2,I[i],j);
	}
	gel(res,4) = ZpXQMinv(M,T,pe,p,e);
  return gerepilecopy(av,res);
}

ulong FindMod(GEN P, GEN Z, ulong n, GEN p, int check)
{ /* Finds 1<=i<=n such that P = Z[i], else returns 0 */
  pari_sp av = avma;
  ulong i;
  GEN Pp;
  if(n==0) return 0; /* Empty list */
  Pp = FpXV_red(P,p);
  for(i=1;i<=n;i++)
  {
    if(gequal(Pp,FpXV_red(gel(Z,i),p)))
    {
      if(check && !gequal(P,gel(Z,i))) pari_err(e_MISC,"Points agree mod p but not mod pe");
      avma = av;
      return i;
    }
  }
  avma = av;
  return 0;
}

GEN ApplyAut(GEN aut, GEN P, GEN vars, GEN T, GEN pe, GEN p, long e)
{ /* aut = [X,Y,Z] function of x,y. Return [X/Z,Y/Z]. */
  pari_sp av = avma;
  GEN a,b,c,Q;
  a = FnEvalAt(gel(aut,1),P,vars,T,p,e,pe);
  b = FnEvalAt(gel(aut,2),P,vars,T,p,e,pe);
  c = FnEvalAt(gel(aut,3),P,vars,T,p,e,pe);
  c = ZpXQ_inv(c,T,p,e);
  Q = cgetg(3,t_VEC);
  gel(Q,1) = FpXQ_mul(a,c,T,pe);
  gel(Q,2) = FpXQ_mul(b,c,T,pe);
  return gerepileupto(av,Q);
}

GEN ApplyFrob(GEN P, GEN FrobMat, GEN T, GEN pe)
{
  pari_sp av = avma;
  GEN Q;
  Q = cgetg(3,t_VEC);
  gel(Q,1) = Frob(gel(P,1),FrobMat,T,pe);
  gel(Q,2) = Frob(gel(P,2),FrobMat,T,pe);
  return gerepileupto(av,Q);
}

GEN VecExtend(GEN V)
{ /* Doubles length, leaves second half uninitialised */
  ulong n = lg(V),i;
  GEN V2;
  V2 = cgetg(2*n,t_VEC);
  for(i=1;i<n;i++) gel(V2,i) = gel(V,i);
  return V2;
}

GEN VecSmallExtend(GEN V)
{ /* Doubles length, leaves second half uninitialised */
  ulong n = lg(V),i;
  GEN V2;
  V2 = cgetg(2*n,t_VECSMALL);
  for(i=1;i<n;i++) V2[i] = V[i];
  return V2;
}

GEN AutFrobClosure(GEN P, GEN Auts, GEN vars, GEN FrobMat, GEN T, GEN pe, GEN p, long e)
{ /* Orbit of point P under Frob and Auts, and induced permutations */
  pari_sp av = avma;
  GEN OP,sFrob,sAuts,res;
  ulong nAuts,nO,nmax;
  ulong i,j,k,m,n;

  OP = cgetg(2,t_VEC); /* Orbit of P, will grow as needed */
  gel(OP,1) = P;
  nO = 1; /* Size of orbit */
  sFrob = cgetg(2,t_VECSMALL); /* Perm nduced by Frob, will grow as needed */
  sFrob[1] = 0;
  nAuts = lg(Auts);
  sAuts = cgetg(nAuts,t_VEC); /* Perms induced by Auts, will grow as needed */
  for(i=1;i<nAuts;i++)
  {
    gel(sAuts,i) = cgetg(2,t_VECSMALL);
    gel(sAuts,i)[1] = 0;
  }
	/* A 0 in a permutation means we do not know yet the image by the permutation */
  nmax = 2; /* Current size of vectors. If size nO of orbit reaches this, the vectors must grow! */

  for(n=1;n<=nO;n++)
  {
    /* Do we know what happens to OP[n] for all auts? */
    for(i=0;i<nAuts;i++) /* i=0: Frob. i>0: Auts[i] */
    {
      if(gel(i?gel(sAuts,i):sFrob,n)==0)
      { /* We do not know, let's find out. */
        m = n; /* Index of the point we are following */
        P = gel(OP,m);
        for(;;)
        {
					/* Apply aut to P */
          if(i) P = ApplyAut(gel(Auts,i),P,vars,T,pe,p,e);
          else P = ApplyFrob(P,FrobMat,T,pe);
					/* Is the result a point we already know? */
          k = FindMod(P,OP,nO,p,1);
          if(k)
          { /* We're back to a pt we know, stop search */
            (i?gel(sAuts,i):sFrob)[m] = k;
            break;
          }
          /* This is a new pt. Add it to orbit and create placeholders for its transfos */
          nO++;
          (i?gel(sAuts,i):sFrob)[m] = nO;
          if(nO==nmax) /* Must extend all vectors */
          {
            nmax*=2;
            OP = VecExtend(OP);
            sFrob = VecSmallExtend(sFrob);
            for(j=1;j<nAuts;j++) gel(sAuts,j) = VecSmallExtend(gel(sAuts,j));
          }
          gel(OP,nO)=P;
          m = nO;
          sFrob[nO]=0;
          for(j=1;j<nAuts;j++) gel(sAuts,j)[nO]=0;
        }
      }
    }
  }
  setlg(OP,nO+1);
  setlg(sFrob,nO+1);
  for(i=1;i<nAuts;i++) setlg(gel(sAuts,i),nO+1);
  res = mkvecn(3,OP,sFrob,sAuts);
  return gerepilecopy(av,res);
}

GEN RRInit(GEN f, GEN Auts, ulong g, ulong d0, GEN L, GEN bad, GEN p, ulong a, long e)
{
	pari_sp av = avma;
	long t;
  ulong nZ,nAuts,n,nOP,m,i;
  GEN vars,pe,T,FrobMat,Z,P,FrobCyc,AutsCyc,OP,V1,V2,V3,W0,V,KV,U,J;
	struct pari_mt pt;
	GEN worker,done,E;
	long workid,pending,k;

	vars = variables_vecsmall(f);
	nZ = 5*d0+1; /* min required #pts */

	t = fetch_var();
	name_var(t,"t");
  T = liftint(ffinit(p,a,t));
  pe = powis(p,e);
  FrobMat = ZpXQ_FrobMat(T,p,e,pe);

	if(DEBUGLEVEL) printf("PicInit: Finding points\n");
  n = 0; /* current #pts */
	Z = cgetg(1,t_VEC); /* list of pts */
	/* Initialise empty cycles */
  FrobCyc = cgetg(1,t_VECSMALL);
	nAuts = lg(Auts);
  AutsCyc = cgetg(nAuts,t_VEC);
  for(i=1;i<nAuts;i++)
    gel(AutsCyc,i) = cgetg(1,t_VECSMALL);
	/* Loop until we have enough pts */
  while(n<nZ)
  {
    /* Get new point */
    P = CurveRandPt(f,T,p,e,bad);
    /* Is it new mod p ? */
    if(FindMod(P,Z,n,p,0)) continue;
    if(DEBUGLEVEL) printf("Got new pt\n");
    /* Compute closure under Frob and Auts */
    OP = AutFrobClosure(P,Auts,vars,FrobMat,T,pe,p,e);
    nOP = lg(gel(OP,1))-1; /* # new pts */
    if(DEBUGLEVEL) printf("Got closure of size %lu\n",nOP);
    /* Add new pts */
    Z = gconcat(Z,gel(OP,1));
    /* Shift permutation describing Frob and Auts */
    for(m=1;m<=nOP;m++)
    {
      gel(OP,2)[m] += n;
      for(i=1;i<nAuts;i++)
        gmael(OP,3,i)[m] += n;
    }
    /* Add these permutaton data */
    FrobCyc = gconcat(FrobCyc,gel(OP,2));
    for(i=1;i<nAuts;i++)
      gel(AutsCyc,i) = gconcat(gel(AutsCyc,1),gmael(OP,3,i));
    /* Update # pts */
    n += nOP;
  }

	if(DEBUGLEVEL) printf("PicInit: Evaluating rational functions\n");
	V1 = FnsEvalAt_Rescale(gel(L,1),Z,vars,T,p,e,pe);
	V2 = FnsEvalAt_Rescale(gel(L,2),Z,vars,T,p,e,pe);
	V3 = DivAdd(V1,V2,3*d0+1-g,T,p,e,pe,0);
	W0 = V1;
	V = mkvecn(3,V1,V2,V3);
	if(DEBUGLEVEL) printf("PicInit: Computing equation matrices\n");
	KV = cgetg(4,t_VEC);
	E = stoi(e);
	worker = strtofunction("mateqnpadic");
  mt_queue_start_lim(&pt,worker,3);
  for(k=1;k<=3||pending;k++)
  {
    mt_queue_submit(&pt,k,k<=3?mkvecn(5,gel(V,k),T,pe,p,E):NULL);
    done = mt_queue_get(&pt,&workid,&pending);
    if(done) gel(KV,workid) = done;
  }
  mt_queue_end(&pt);
	if(DEBUGLEVEL) printf("PicInit: Constructing evaluation maps\n");
	U = RREvalInit(L,vars,Z,V2,T,p,e,pe);
  J = mkvecn(lgJ,f,stoi(g),stoi(d0),L,T,p,stoi(e),pe,FrobMat,V,KV,W0,U,Z,FrobCyc,AutsCyc);
	return gerepilecopy(av,J);
}

GEN Jlift(GEN J, ulong e2)
{
	pari_sp av = avma, avZ;
	GEN J2,T,p,pe2,f,vars,L,FrobCyc,FrobMat2;
	long g,d0;
	GEN Z,Z2,V,W0,KV,P,x,y,fx,U,V2,I,M;
	ulong nZ,i,j,nV2;
	struct pari_mt pt;
  GEN worker,done,E2;
  long pending,k,workid;
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
  gel(J2,7) = E2 = utoi(e2);
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
	V = cgetg(4,t_VEC);
	KV = cgetg(4,t_VEC);
	U = cgetg(5,t_VEC);
	/* TODO why does that parallel version not work?
	printf("RR\n");
	worker = strtofunction("RRspaceEval");
	mt_queue_start(&pt,worker);
  for(k=1;k<=4||pending;k++)
  {
		mt_queue_submit(&pt,k,k<=4?mkvecn(7,gel(L,k),vars,Z2,T,p,E2,pe2):NULL);
		done = mt_queue_get(&pt,&workid,&pending);
		if(done)
		{
			if(workid<=2) gel(V,workid) = gel(done,1);
			else gel(U,workid-2) = done;
		}
	}
	mt_queue_end(&pt);
	printf("End RR\n"); */
	gel(V,1) = gel(RRspaceEval(gel(L,1),vars,Z2,T,p,e2,pe2),1);
	V2 = gel(V,2) = gel(RRspaceEval(gel(L,2),vars,Z2,T,p,e2,pe2),1);
	gel(U,1) = RRspaceEval(gel(L,3),vars,Z2,T,p,e2,pe2);
	gel(U,2) = RRspaceEval(gel(L,4),vars,Z2,T,p,e2,pe2);
  gel(V,3) = DivAdd(gel(V,1),gel(V,2),3*d0+1-g,T,p,e2,pe2,0);
  W0 = gel(V,1); /* TODO can it happen that W0 != V1 even though all data is present? */
	I = gel(U,3) = gmael(J,13,3);
	nV2 = lg(V2);
	M = cgetg(nV2,t_MAT);
  for(j=1;j<nV2;j++)
  {
    gel(M,j) = cgetg(nV2,t_COL);
    for(i=1;i<nV2;i++)
      gcoeff(M,i,j) = gcoeff(V2,I[i],j);
  }
  gel(U,4) = ZpXQMinv(M,T,pe2,p,e2);
	worker = strtofunction("mateqnpadic");
  mt_queue_start_lim(&pt,worker,3);
  for(k=1;k<=3||pending;k++)
  {
    mt_queue_submit(&pt,k,k<=3?mkvecn(5,gel(V,k),T,pe2,p,E2):NULL);
    done = mt_queue_get(&pt,&workid,&pending);
    if(done) gel(KV,workid) = done;
  }
  mt_queue_end(&pt);
  gel(J2,10) = V;
  gel(J2,11) = KV;
  gel(J2,12) = W0;
  gel(J2,13) = U;
  gel(J2,14) = Z2;
  gel(J2,15) = JgetFrobCyc(J);
  gel(J2,16) = JgetAutsCyc(J);
  return gerepilecopy(av,J2);
}

GEN RREval(GEN J, GEN W)
{
	pari_sp av = avma, av1, av2;
	GEN T,p,pe,V,KV,U,res,resi1;
	long e;
	ulong n1,n2,i1,i2;
	GEN S1,S2,I,M,s2,s2I,K;
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
	I = gel(U,3); /* Row indices to look at to ID an elt of V */
	M = gel(U,4); /* Matrix to apply to the I-entries */

	res = cgetg(n1,t_MAT);
	for(i1=1;i1<n1;i1++)
	{
		av1 = avma;
		S1 = DivAdd(W,gmael(U,1,i1),2*d0+1,T,p,e,pe,0); /* L(4D0-D-E1) */
		S1 = DivSub(V,S1,KV,1,T,p,e,pe,2); /* L(2D0-D-E1), generically 1-dimensional */
		S1 = gerepileupto(av1,gel(S1,1));
		S1 = DivMul(S1,V,T,pe); /* L(4D0-D-E1-ED) */
		S1 = gerepileupto(av1,S1);
		S1 = DivSub(W,S1,KV,d0+1-g,T,p,e,pe,2); /* L(2D0-E1-ED) */
		S1 = gerepileupto(av1,S1);
		resi1 = cgetg(n2,t_COL);
		for(i2=1;i2<n2;i2++)
		{
			av2 = avma;
			S2 = DivAdd(S1,gmael(U,2,i2),2*d0+1,T,p,e,pe,0); /* L(4D0-E1-E2-ED) */
			S2 = gerepileupto(av2,S2);
			S2 = DivSub(V,S2,KV,1,T,p,e,pe,2); /* L(2D0-E1-E2-ED), generically 1-dimensional */
  		s2 = gel(S2,1); /* Generator */
			s2 = gerepileupto(av2,s2);
			/* get coords of s2 w.r.t. V */
			s2I = cgetg(nV,t_COL);
			for(i=1;i<nV;i++)
				gel(s2I,i) = gel(s2,I[i]);
  		K = FqM_FqC_mul(M,s2I,T,pe);
			gel(resi1,i2) = gerepileupto(av2,K);
		}
		gel(res,i1) = gerepileupto(av1,resi1);
	}
	return gerepileupto(av,res);
}

GEN RREval_worker(GEN W, GEN J)
{
	return RREval(J,W);
}
