#include "pic.h"
#include "linalg.h"
#include "freyruck.h"

GEN FnSubstMod(GEN F, long var, GEN val, GEN T, GEN pe) /* /!\ Not memory-clean */
{
	GEN valm,res;
	pari_CATCH(e_INV)
	{
		//printf("Z");
		res = gsubst(F,var,val);
		return res;
		/*printf("u");
		res = gmodulo(gmodulo(res,T),pe);
		printf("v");
		res = liftall(res);
		printf("w");
		return res;*/
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

GEN FnEvalAt(GEN F, GEN P, GEN vars, GEN T, GEN p, long e, GEN pe) /* /!\ Not memory-clean */
{
	GEN N,D,C;
	//printf(" a");
	if(typ(F)==t_RFRAC)
	{
		//printf("b");
		N = FnEvalAt(gel(F,1),P,vars,T,p,e,pe);
		if(typ(N)==t_INT) N=Z2Fq(N,T);
		D = FnEvalAt(gel(F,2),P,vars,T,p,e,pe);
		if(typ(D)==t_INT) D=Z2Fq(D,T);
		//printf("c\n");
		return ZpXQ_div(N,D,T,pe,p,e);
	}
	//printf("d");
	C = content(F); // TODO wrt y
	F = gdiv(F,C);
	F = FnSubstMod(F,vars[1],gel(P,1),T,pe);
	F = gmul(F,C);
	//printf("e");
	if(typ(F)==t_RFRAC)
	{
		//printf("f");
		N = FnSubstMod(gel(F,1),vars[2],gel(P,2),T,pe);
		if(typ(N)==t_INT) N=Z2Fq(N,T);
		D = FnSubstMod(gel(F,2),vars[2],gel(P,2),T,pe);
		if(typ(D)==t_INT) D=Z2Fq(D,T);
		//printf("g");
		/*pari_printf(" %Ps %Ps ",N,D);*/
		F = ZpXQ_div(N,D,T,pe,p,e);
		//printf("h");
	}
	else{//printf("i");
	F = liftall(gsubst(F,vars[2],gmodulo(gmodulo(gel(P,2),T),pe)));
	}
	//printf("j\n");
	return F;
}

GEN FnsEvalAt(GEN Fns, GEN Z, GEN vars, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN A,col;
	ulong i,j,nF,nZ;
	nF = lg(Fns);
	nZ = lg(Z);
	A = cgetg(nF,t_MAT);
	for(j=1;j<nF;j++)
	{
		printf("Eval fn %lu out of %lu\n",j,nF-1);
		col = cgetg(nZ,t_COL);
		for(i=1;i<nZ;i++)
		{
			gel(col,i) = Fq_red(FnEvalAt(gel(Fns,j),gel(Z,i),vars,T,p,e,pe),T,pe);
		}
		gel(A,j) = col;
	}
	return gerepilecopy(av,A);
}

GEN FnsEvalAt_Rescale(GEN Fns, GEN Z, GEN vars, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN F,S,K,f;
	ulong i,j,k,nF;
	F = gcopy(Fns);
	nF = lg(F);
	while(1)
	{
		/*Eval the fns */
		S = FnsEvalAt(F,Z,vars,T,p,e,pe);
		K = FqM_ker(S,T,p);
		/* Are the evals (and hence the fns) independent ? */
		if(lg(K)==1)
		{
			printf("Good, no relation\n");
			return gerepileupto(av,S);
		}
		pari_printf("Found %ld relations, eliminating and restarting\n",lg(K)-1);
		/* No. We assume Z def / Q, so K has entries in Fp */
		/* Do elimination and start over */
		for(j=1;j<lg(K);j++)
		{
			/* k = pivot = last nonzero entry of the col (it's a 1) */
			k = 0;
			for(i=1;i<nF;i++)
			{
				if(!gequal0(gcoeff(K,i,j))) k=i;
			}
			/* Form corresponding lin comb, and div by p */
			f = gel(F,k);
			for(i=1;i<k;i++)
			{
				if(!gequal0(gcoeff(K,i,j)))
				{
					f = gadd(f,gmul(gcoeff(K,i,j),gel(F,i)));
				}
			}
			gel(F,k) = gdiv(f,p);
		}
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

GEN RRInit(GEN f, ulong g, ulong d0, GEN L, GEN bad, GEN p, ulong a, long e)
{
	pari_sp avP,av = avma;
  int newpt;
  ulong nZ,n,ncyc,i;
  GEN vars,pe,t,T,Frob,Z,Zp,P,Pp,Q,FrobCyc,x,y,V1,V2,V3,W0,V,KV,KV3,J;

	vars = variables_vecsmall(f);
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
      x = FpX_FpXQ_eval(gel(Q,1),Frob,T,pe);
      y = FpX_FpXQ_eval(gel(Q,2),Frob,T,pe);
      Q = mkvec2(x,y);
    } while(!gequal(Q,P));
  FrobCyc[ncyc] = i;
  }
  setlg(Z,n+1);
  setlg(FrobCyc,ncyc+1);
	//pari_printf("T=%Ps\nZ=%Ps\n",T,Z);

	V1 = FnsEvalAt_Rescale(L,Z,vars,T,p,e,pe);
	//printf("Rank V1: %ld, expected %ld\n",FqM_rank(V1,T,p),lg(V1)-1);
	//pari_printf("%Ps\n",FqM_ker(V1,T,p));
	V2 = DivAdd1(V1,V1,2*d0+1-g,T,p,e,pe,0);
	V3 = DivAdd1(V1,V2,3*d0+1-g,T,p,e,pe,0);
	W0 = V1;
	V = V2;
  KV = mateqnpadic(V,T,p,e);
  KV3 = mateqnpadic(V3,T,p,e);

  J = mkvecn(lgJ,f,stoi(g),stoi(d0),T,p,stoi(e),pe,Frob,V,KV,W0,Z,FrobCyc,V3,KV3);
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

	S1 = DivAdd(W,gel(Li,1),2*d0+1,T,p,e,pe,0); /* L(4D0-D-E1) */
	S1 = DivSub(V,S1,KV,1,T,p,e,pe,2); /* L(2D0-D-E1) */
	S2 = DivMul(gel(S1,1),V,T,pe); /* L(4D0-D-E1-ED) */
	S2 = DivSub(W,S2,KV,d0+1-g,T,p,e,pe,2); /* L(2D0-E1-ED) */
	S2 = DivAdd(S2,gel(Li,2),2*d0+1,T,p,e,pe,0); /* L(4D0-E1-E2-ED) */
	S2 = DivSub(V,S2,KV,1,T,p,e,pe,2); /* L(2D0-E1-E2-ED) */
  s2 = gel(S2,1);
  K = cgetg(nV+1,t_MAT);
  for(i=1;i<nV;i++) gel(K,i) = gel(V,i);
  gel(K,nV) = s2;
  K = matkerpadic(K,T,p,e);
	return gerepileupto(av,gel(K,1));
}

GEN PolExpId(GEN Z, GEN T, GEN pe) /* bestappr of prod(x-z), z in Z */
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
	/*printf("Inverting\n");*/
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
	/*printf("Getting %lu pols\n",npols);*/
	pols = cgetg(npols+1,t_VEC);
	R = cgetg(nF,t_VEC);
	/* TODO parallelise */
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
			/*printf("%lu ",n);*/
			gel(pols,n) = PolExpId(R,T,pe);
		}
	}		
  
	return gerepilecopy(av,pols);
}
