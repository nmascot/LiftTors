#include "linalg.h"
#include "pic.h"

GEN PicDeflate(GEN J, GEN W, ulong nIGS)
{ /* Finds nIGS elts w1, .. , wn in W s.t. (W) = min (wi) */
	pari_sp av,av1;
	GEN V,T,pe,p,IGS,IV,wV;
	ulong g,nV,nW,r,i,j,k;
	long e;

	V = JgetV(J);
	JgetTpe(J,&T,&pe,&p,&e);
	g = Jgetg(J);
	nV = lg(V)-1;
	nW = lg(W)-1;
	
	IGS = cgetg(nIGS+1,t_VEC);
	av = avma;
	while(1)
	{
		avma = av;
		for(i=1;i<=nIGS;i++) gel(IGS,i) = RandVec_padic(W,T,p,pe);
		av1 = avma;
		IV = cgetg(nIGS*nV+1,t_MAT);
		k=1;
		for(i=1;i<=nIGS;i++)
		{
			wV = DivMul(gel(IGS,i),V,T,p);
			for(j=1;j<=nV;j++)
			{
				gel(IV,k) = gel(wV,j);
				k++;
			}
		}
		r = FqM_rank(IV,T,p);
		if(r==nV+nW+g-1)
		{
			avma = av1;
			return IGS;
		}
		printf("IGS[%lu,%lu]\n",r,nV+nW+g-1);
	}
}

GEN PicDeflate_U(GEN J, GEN W, ulong nIGS)
{
	pari_sp av = avma;
	GEN V,T,p;
	long e;
	ulong nV;
	GEN GW,K,U;
	ulong j;

	V = JgetV(J);
	T = JgetT(J);
	p = Jgetp(J);
	e = Jgete(J);
	nV = lg(V);

	GW = PicDeflate(J,W,nIGS); /* IGS of W */
  K = cgetg(nV+nIGS,t_MAT);
  for(j=1;j<nV;j++) gel(K,j) = gel(V,j);
  for(j=1;j<=nIGS;j++) gel(K,nV-1+j) = gel(GW,j);
  U = matkerpadic(K,T,p,e); /* Coords of IGS of W // basis of V */
  for(j=1;j<=nIGS;j++) setlg(gel(U,j),nV);
  return gerepilecopy(av,U);
}

GEN PicInflate_U(GEN J, GEN U) /* Takes IGS given by coords // V */
{
	pari_sp av = avma;
	GEN T,pe,p;
	long e;
	GEN V,KV,GWV,wV,W;
	ulong d0,g;
	ulong nU,nV;
	ulong i,j,k;

	JgetTpe(J,&T,&pe,&p,&e);
	V = JgetV(J);
	KV = JgetKV(J);
	nU = lg(U)-1;
	nV = lg(V)-1;
	d0 = Jgetd0(J);
	g = Jgetg(J);

	GWV = cgetg(nU*nV+1,t_MAT); /* w*V for w in GW */
  k = 1;
  for(i=1;i<=nU;i++)
  {
    wV = DivMul(FqM_FqC_mul(V,gel(U,i),T,pe),V,T,pe); /* TODO useful to precompte V[i]*V ? */
    for(j=1;j<=nV;j++)
    {
      gel(GWV,k) = gel(wV,j);
      k++;
    }
  }
  GWV = FqM_image(GWV,T,p);
  W = DivSub(V,GWV,KV,d0+1-g,T,p,e,pe,3); /* TODO pass precomputed IGS of V */
	return gerepileupto(av,W);
}

GEN PicLift_RandLift_U(GEN U, GEN U0, GEN KM, GEN T, GEN p, GEN pe1, GEN pe21, long e21)
{
	pari_sp av;

	GEN K,red,newU;
	ulong nU,nU0,nV;
	ulong i,j,k,m;

	nU = lg(U);
	nU0 = lg(U0);
	nV = lg(gel(U,1));

	av=avma;
  do
  {
    avma=av;
    K = RandVec_1(KM,pe21);
    red = gel(K,lg(K)-1);
  } while(ZX_is0mod(red,p));
  red = ZpXQ_inv(red,T,p,e21);
  setlg(K,lg(K)-1);
  K = FqC_Fq_mul(K,red,T,pe21);
  newU = cgetg(nU,t_VEC);
  k = 1;
  for(i=1;i<nU;i++)
  {
		gel(newU,i) = cgetg(lg(gel(U,i)),t_COL);
		/* Correct U[i] */
    for(j=1;j<nU0;j++)
    { /* Add the proper multiple of U0[j] to it */
      for(m=1;m<nV;m++)
        gmael(newU,i,m) =
					ZX_add(gmael(U,i,m),ZX_Z_mul(Fq_mul(gel(K,k),gcoeff(U0,m,j),T,pe21),pe1));
      k++;
    }
  }
	return gerepileupto(av,newU);
}

GEN PicLiftTors_Chart_worker(GEN J, GEN l, GEN U, GEN U0, GEN KM, GEN pe1, GEN pe21, long e21, GEN c0, ulong k0, ulong P0, GEN randseed)
{
  pari_sp av = avma;
	GEN T,p,pe2;
	long e2;
	GEN W,c;
  ulong nc,i;
  setrand(randseed);
	JgetTpe(J,&T,&pe2,&p,&e2);
  nc = lg(c0)-1;

	/* Get a random lift */
	U = PicLift_RandLift_U(U,U0,KM,T,p,pe1,pe21,e21);
  W = PicInflate_U(J,U);
  /* Mul by l, get coordinates, and compare them to those of W0 */
  c = PicChart(J,PicMul(J,W,l,0),P0);
  c = FqV_Fq_mul(c,ZpXQ_inv(gel(c,k0),T,p,e2),T,pe2); /* Normalize proj coords */
  for(i=1;i<=nc;i++) /* The coords are c0 mod pe1 -> divide */
    gel(c,i) = ZX_Z_divexact(FpX_sub(gel(c,i),gel(c0,i),pe2),pe1);
  return gerepilecopy(av,mkvecn(3,U,W,c));
}

GEN PicLiftTors(GEN J, GEN W, long eini, GEN l)
{
  pari_sp av=avma,av1,av2;
	GEN T,p,V;
  long efin,e1,e2,e21;
  GEN pefin,pe1,pe21,pe2;
  GEN J1,J2;
  GEN sW,Vs,U0,V0;
  GEN K,U;
	GEN GWV,wV;
  GEN ABCD,uv,Ainv,CAinv,AinvB,rho;
  GEN drho,abcd;
  GEN KM,red;
  ulong g,nZ,nW;
  ulong nGW=2,nV,d0,nc,k0=0,P0=0;
	GEN c0=NULL;
	int P0_tested=0;
	GEN Clifts,Ulifts,Wlifts,Ktors;
	struct pari_mt pt;
  GEN randseed,args,worker,done;
  long pending,workid;
  ulong r,i,j,k,n;

	JgetTpe(J,&T,&pefin,&p,&efin);
	g = Jgetg(J);
  d0 = Jgetd0(J);
  V = JgetV(J);
  /*GV = JgetGV(J2);*/
  nV = lg(V)-1;
  nZ = lg(gel(V,1))-1;
  nW = lg(W)-1;

  sW = gel(FqM_indexrank(W,T,p),1); /* rows s.t. this block is invertible, # = nW, we won't change them */
  Vs = cgetg(nV+1,t_MAT); /* V with only the rows in sW */
  for(j=1;j<=nV;j++)
  {
    gel(Vs,j) = cgetg(nW+1,t_COL);
    for(i=1;i<=nW;i++) gcoeff(Vs,i,j) = gcoeff(V,sW[i],j);
  }
	U0 = matkerpadic(Vs,T,p,efin); /* TODO half prec ? */ /* # = nV-nW = d0 */
	U0 = gerepileupto(av,U0);
	V0 = cgetg(d0+1,t_VEC);
	for(i=1;i<=d0;i++) gel(V0,i) = DivMul(FqM_FqC_mul(V,gel(U0,i),T,pefin),V,T,pefin); /* s*V for s in subspace of V whose rows in sW are 0 */ /* TODO can half precision here */
	/* TODO parallel */

	e1 = eini;
	pe1 = powiu(p,e1);
	av1 = avma; /* Use to collect garbage at each iteration */
	J1 = PicRed(J,e1);
	U = PicDeflate_U(J1,W,nGW); /* IGS of W1 // basis of V */

	r = 3*d0+1-g; /* Wanted rank of GWV */
	/* Main loop */
  while(1)
  {
    e2 = 2*e1;
    if(e2>efin) e2 = efin;
		e21 = e2-e1;
		pe21 = e21==e1 ? pe1 : powiu(p,e21);
    pari_printf("Lifting from prec O(%Ps^%lu) to O(%Ps^%lu)\n",p,e1,p,e2);
    J2 = e2<efin ? PicRed(J,e2) : J;
		pe2 = Jgetpe(J2);
		/* START LIFTING */
  	GWV = cgetg(nGW*nV+1,t_MAT); /* w*V for w in GW */
		/* We need it to have rk r, it is already the case mod pe1 */
  	k = 1;
  	for(i=1;i<=nGW;i++)
  	{
    	wV = DivMul(FqM_FqC_mul(V,gel(U,i),T,pe2),V,T,pe2);
    	for(j=1;j<=nV;j++)
    	{
      	gel(GWV,k) = gel(wV,j);
      	k++;
    	}
  	}
  	
  	uv = FqM_MinorCompl(GWV,T,p); /* How to split GWV */
  	ABCD = M2ABCD(GWV,uv); /* Splitting */
  	Ainv = ZpXQM_inv(gel(ABCD,1),T,p,e2);
  	CAinv = FqM_mul(gel(ABCD,3),Ainv,T,pe2);
  	AinvB = FqM_mul(Ainv,gel(ABCD,2),T,pe2);
  	rho = FqM_mul(CAinv,gel(ABCD,2),T,pe2);
  	rho = FpXM_sub(gel(ABCD,4),rho,pe2); /* size nZ-r,nGW*nV-r(=dW if nGW-2); tests if rk = r */
  	for(i=1;i<=nZ-r;i++)
  	{
    	for(j=1;j<=nGW*nV-r;j++) gcoeff(rho,i,j) = ZX_Z_divexact(gcoeff(rho,i,j),pe1);
    }
		
  	/* Now deform the w in GW by p^e1*V0. Actually we deform U1 by p^e1*U0. */
  	/* TODO do not deform U1[1]? */
  	/* TODO swap loops and parallelise */
  	K = cgetg(nGW*d0+2,t_MAT);
  	k = 1;
  	for(i=1;i<=nGW;i++)
  	{
    	for(j=1;j<=d0;j++)
    	{ /* GW[i] shifts by p^e1*V0[j] */
      	abcd = M2ABCD_1block(gel(V0,j),0,(i-1)*nV,uv); /* Split */
      	drho = FpXM_sub(FqM_mul(gel(abcd,1),AinvB,T,pe21),gel(abcd,2),pe21); /* aA^-1B-b */
      	drho = FqM_mul(CAinv,drho,T,pe21); /* CA^-1aA^-1B - CA^-1b */
      	drho = FpXM_add(gel(abcd,4),drho,pe21); /* d + CA^-1aA^-1B - CA^-1b */
      	drho = FpXM_sub(drho,FqM_mul(gel(abcd,3),AinvB,T,pe21),pe21); /* d + CA^-1aA^-1B - CA^-1b - cA^-1B */
      	gel(K,k) = mat2col(drho);
      	k++;
    	}
  	}
  	gel(K,nGW*d0+1) = mat2col(rho);
  	/* Find a random solution to the inhomogeneous system */
  	KM = matkerpadic(K,T,p,e21);
		printf("dim ker lift: %ld\n",lg(KM)-1);
		U = PicLift_RandLift_U(U,U0,KM,T,p,pe1,pe21,e21);
		e1 = e2;
    pe1 = pe2;
    printf("END LOOP\n");
		if(e2==efin)
		{
			W = PicInflate_U(J2,U);
			return gerepileupto(av,W);
		}
		continue;
		if(cmpii(pe21,powiu(l,g+1))<=0)
  	{
    	printf("Lift by mul\n");
			U = PicLift_RandLift_U(U,U0,KM,T,p,pe1,pe21,e21);
			printf("A");
			W = PicInflate_U(J2,U);
			printf("B");
    	W = PicMul(J2,W,pe21,0); /* Make it l-tors */
			printf("C");
			if(e2==efin) /* Already done ? */
				return gerepileupto(av,W);
			printf("D\n");
			U = PicDeflate_U(J2,W,nGW); /* Update U -> ready for new iteration */
		}
		else
		{
			printf("Lift by chart\n");
			Clifts = cgetg(g+2,t_MAT);
			Ulifts = cgetg(g+2,t_VEC);
			Wlifts = cgetg(g+2,t_VEC);
  		worker = strtofunction("PicLiftTors_worker");
  		av2 = avma;
  		while(1)
  		{
				avma = av2;
    		if(c0==NULL) /* Compute coords of 0 if not already done */
    		{
      		/* Find coords of 0 */
      		c0 = PicChart(J,JgetW0(J),P0);
      		nc = lg(c0)-1;
      		/* Find index to dehomogenise */
      		red = NULL;
      		for(k0=1;k0<=nc;k0++)
      		{
        		red = gel(c0,k0);
        		if(!ZX_is0mod(red,p)) break;
      		}
      		c0 = FqV_Fq_mul(c0,ZpXQ_inv(red,T,p,e2),T,pe2); /* Dehomogenize */
					c0 = gerepileupto(av2,c0);
					av2 = avma;
    		}
    		/* Find g+1 lifts in parallel */
    		pending = 0;
    		mt_queue_start(&pt,worker);
    		for(i=1;i<=g+1||pending;i++)
    		{
      		if(i<=g+1)
      		{
        		randseed = utoi(pari_rand());
        		args = mkvecn(12,J2,l,U,U0,KM,pe1,pe21,stoi(e21),c0,utoi(k0),utoi(P0),randseed);
        		mt_queue_submit(&pt,i,args);
      		}
      		else mt_queue_submit(&pt,i,NULL);
      		done = mt_queue_get(&pt,&workid,&pending);
      		if(done)
      		{
        		gel(Ulifts,workid) = gel(done,1);
        		gel(Wlifts,workid) = gel(done,2);
        		gel(Clifts,workid) = gel(done,3);
      		}
    		}
    		mt_queue_end(&pt);
				Ktors = matkerpadic(Clifts,T,p,e21); /* Find comb with coord = 0 */
    		n = lg(Ktors)-1;
    		if(n>1)
    		{ /* l-tors is étale, so this can only happen if CHart is not diffeo - > change chart */
      		printf("Dim ker tors = %ld, changing charts\n",n);
      		P0++; /* New chart */
					P0_tested = 0;
					c0 = NULL; /* Coords of 0 must be recomputed */
      		continue; /* Try again with this new chart */
    		}
    		Ktors = gel(Ktors,1);
    		red = gel(Ktors,1);
    		for(i=2;i<=g+1;i++)
    		{
      		red = FpX_add(red,gel(Ktors,i),pe2);
    		}
    		if(ZX_is0mod(red,p)) /* TODO can this happen ? why, or why not ? */
      		printf("Sum of Ktors is zero!\n"); continue;
    		Ktors = FqC_Fq_mul(Ktors,ZpXQ_inv(red,T,p,e2),T,pe2); /* Normalise so that sum = 1 */
				U = W = NULL; /* If done, return updated W; else update U. */
				/* But first check if really l-tors, as the chart might not be injective ! */
    		if(P0_tested == 0)
				{
					for(i=1;i<=g+1;i++) gel(Wlifts,i) = FqM_Fq_mul(gel(Wlifts,i),gel(Ktors,i),T,pe21);
	    		W = gel(Wlifts,1);
  	  		for(i=2;i<=g+1;i++) W = FpXM_add(W,gel(Wlifts,i),pe2);
      		if(!PicIsZero(J2,PicMul(J2,W,l,0)))
      		{
        		printf("Not actually l-torsion!!! Changing charts\n");
        		P0++;
        		c0 = NULL;
        		continue;
      		}
					P0_tested = 1;
    		}
				if(e2 == efin)
				{ /* Alredy done ? */
					if(W == NULL) /* Update W, if notalready done, and return it */
					{
					}
					return gerepileupto(av,W);
				}
				else
				{ /* Update U -> ready for next iteration */
					for(i=1;i<=g+1;i++) gel(Ulifts,i) = FqM_Fq_mul(gel(Ulifts,i),gel(Ktors,i),T,pe21);
    			U = gel(Ulifts,1);
    			for(i=2;i<=g+1;i++) U = FpXM_add(U,gel(Wlifts,i),pe2);
				}
			}
		}
		/* END LIFTING */
    e1 = e2;
		pe1 = pe2;
		printf("END LOOP\n");
		/*if(c0)
			gerepileall(av1,3,&U,&pe1,&c0);
		else
			gerepileall(av1,2,&U,&pe1);*/
  }
}

