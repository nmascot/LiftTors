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

GEN PicInflateV(GEN J, GEN U) /* Takes IGS given by coords // V */
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

GEN PicLift_RandLiftU(GEN U, GEN U0, GEN KM, GEN T, GEN p, GEN pe1, GEN pe21, long e21)
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

  newU = cgetg(nU+1,t_VEC);
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

GEN PicLift_XP(GEN J, GEN W, long eini)
{
  pari_sp av=avma,av1;
	GEN T,p,V;
  long efin,e1,e2,e21;
  GEN pefin,pe1,pe21,pe2;
  GEN J1,J2;
  GEN sW,Vs,U0,V0;
  GEN GW,K,U;
	GEN GWV,wV;
  GEN ABCD,uv,Ainv,CAinv,AinvB,rho;
  GEN drho,abcd;
  GEN KM,red;
  ulong g,nZ,nW;
  ulong nGW=2,nV,d0;
  ulong r,i,j,k,m;

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
	av1 = avma;
	J1 = PicRed(J,e1);
	GW = PicDeflate(J1,W,nGW); /* IGS of W1 */
	K = cgetg(nV+nGW+1,t_MAT);
	/* TODO reduce V first for performance */
	/* TODO check other places */
	for(j=1;j<=nV;j++) gel(K,j) = gel(V,j);
	for(j=1;j<=nGW;j++) gel(K,nV+j) = gel(GW,j);
	U = matkerpadic(K,T,p,e1); /* Coords of IGS of W1 // basis of V */
	for(j=1;j<=nGW;j++) setlg(gel(U,j),nV+1); /* U is what we will lift */
	U = gerepilecopy(av1,U);

	/* P0 = NULL; */
	r = 3*d0+1-g; /* Wanted rank of GWV */
	/* Main loop */
  while(e1<efin)
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
		if(cmpii(pe21,powiu(l,g+1))<=0)
  	{
    	printf("Lift by mul\n");
			U = PicLift_RandLiftU(U,U0,KM,T,p,pe1,pe21,e21);
			W = PicInflateV(J2,U);
    	W = PicMul(J2,W,pe21,0);
		}
		else
		{
			TODO;
		}
  	av1=avma;
  	do
  	{
    	avma=av1;
    	K = RandVec_1(KM,pe21);
    	red = gel(K,1+d0*nGW);
  	} while(ZX_is0mod(red,p));
	  red = ZpXQ_inv(red,T,p,e21);
	  setlg(K,nGW*d0+1);
	  K = FqC_Fq_mul(K,red,T,pe21);
  	/* Correct U */
  	k = 1;
  	for(i=1;i<=nGW;i++)
  	{ /* Correct U[i] */
    	for(j=1;j<=d0;j++)
    	{ /* Add the proper multiple of U0[j] to it */
      	for(m=1;m<=nV;m++)
        	gcoeff(U,m,i) = ZX_add(gcoeff(U,m,i),ZX_Z_mul(Fq_mul(gel(K,k),gcoeff(U0,m,j),T,pe21),pe1));
      	k++;
    	}
  	}
		/* END LIFTING */
		/* TODO memory management */
    /*P0 = gel(W,2);
    if(P0)
    {
      W = gerepileupto(av,W);
      P0 = gel(W,2);
      W = gel(W,1);
    }
    else
    {
      W = gerepileupto(av,gel(W,1));
    }*/
    e1 = e2;
		pe1 = pe2;
  }

	W = PicInflateV(J,U);
	return gerepileupto(av,W);
}








GEN PicLift_worker(ulong nW, ulong nZ, ulong nKV, GEN KwVk, GEN VFlist, GEN uv, GEN AinvB, GEN CAinv, GEN T, GEN pe21)
{
	pari_sp av = avma;
	GEN M,dK,dKi,ABCD;
	ulong j,P,i,e;

	M = cgetg(nW+1,t_MAT);
	for(j=1;j<=nW;j++) gel(M,j) = NULL;
  for(j=1;j<=nW;j++)
  {
		if(j==1)
		{
			dK = cgetg(nZ+1,t_MAT);
			for(P=1;P<=nZ;P++) gel(dK,P) = cgetg(nKV*(nW-1)+1,t_COL);
			for(i=2;i<=nW;i++)
			{
				dKi = FqM_mul(KwVk,gel(VFlist,i),T,pe21);
				for(P=1;P<=nZ;P++)
				{
					for(e=1;e<=nKV;e++)
					{
						gcoeff(dK,(i-2)*nKV+e,P) = gcoeff(dKi,e,P);
					}
				}
			}
			ABCD = M2ABCD_1block(dK,nKV,0,uv);
		}
		else
		{
			ABCD = M2ABCD_1block(KwVk,(j-1)*nKV,0,uv);
		}
		dK = FpXM_sub(FqM_mul(gel(ABCD,1),AinvB,T,pe21),gel(ABCD,2),pe21);
		dK = FqM_mul(CAinv,dK,T,pe21);
		dK = FpXM_add(gel(ABCD,4),dK,pe21);
		dK = FpXM_sub(dK,FqM_mul(gel(ABCD,3),AinvB,T,pe21),pe21);
		gel(M,j) = mat2col(dK);
  }
	return gerepileupto(av,M);
}

/*GEN PicLift_RandLift(GEN W1, GEN KM, GEN V0, ulong nZ, ulong nW, ulong d0, GEN T, GEN p, long e21, GEN pe1, GEN pe2, GEN pe21)
{
	pari_sp av = avma;
	GEN K,red,W;
	ulong n,j,k,P;
  do
  {
    K = RandVec_1(KM,pe21);
    red = gel(K,1+d0*nW);
  } while(ZX_is0mod(red,p));
  red = ZpXQ_inv(red,T,p,e21);
  setlg(K,d0*nW+1);
  K = FqV_Fq_mul(K,red,T,pe21);

  n = 0;
  W = zeromatcopy(nZ,nW);
  for(k=1;k<=d0;k++)
  {
    for(j=1;j<=nW;j++)
    {
      n++;
      for(P=1;P<=nZ;P++)
      {
        gcoeff(W,P,j) = FpX_add(gcoeff(W,P,j),Fq_mul(gel(K,n),gcoeff(V0,P,k),T,pe21),pe21);
      }
    }
  }
  W = FpXM_add(W1,ZXM_Z_mul(W,pe1),pe2);
	return gerepileupto(av,W);
}*/

GEN PicLiftTors_worker(GEN J, GEN W1, GEN l, GEN KM, GEN c0, GEN V0, ulong d0, ulong nW, ulong nZ, ulong k0, GEN T, GEN p, long e2, GEN pe2, long e21, GEN pe21, GEN pe1, GEN randseed, ulong P0)
{
	pari_sp av = avma;
	GEN W,c;
	ulong nc,i;
	setrand(randseed);
	nc = lg(c0)-1;
	/* Find a random solution to the inhomogeneous system */
	W = PicLift_RandLift(W1,KM,V0,nZ,nW,d0,T,p,e21,pe1,pe2,pe21);
  /* Mul by l, get coordinates, and compare them to those of W0 */
  c = PicChart(J,PicMul(J,W,l,0),P0);
  c = FqV_Fq_mul(c,ZpXQ_inv(gel(c,k0),T,p,e2),T,pe2);
  for(i=1;i<=nc;i++)
  {
    gel(c,i) = ZX_Z_divexact(FpX_sub(gel(c,i),gel(c0,i),pe2),pe1);
  }
	return gerepilecopy(av,mkvec2(W,c));
}

GEN PicLiftTors_2(GEN J2, GEN W1, long e1, GEN l, GEN P0_hint)
{
	pari_sp av1,av2,av=avma;
	GEN V,KV,W0,T,p,pe1,pe2,pe21;
	GEN col,K,wV,KwV,Ainv,AinvB,CAinv,rho,ABCD,uv;
	GEN F,VF,VFlist,sW,cW,Vs,V0,KwVlist,M,KM,worker,done;
	GEN c0,Wlifts,W,red,Ktors;
	GEN NW,NZ,NKV,D0,K0,E2,E21,randseed,args;
	ulong g,d0,nW,nV,nKV,nZ,nc,r;
	long e2,e21,pending,workid;
	ulong e,i,j,k,P,k0,n;
	GEN P0;
	int todo_c0;

	struct pari_mt pt;

	JgetTpe(J2,&T,&pe2,&p,&e2);
	if(e2<=e1)
	{
		return FpXM_red(W1,pe2);
	}
	if(P0_hint) P0 = P0_hint;
	else P0 = gen_0;

	V = JgetV(J2);
	KV = JgetKV(J2);
	W0 = JgetW0(J2);
	g = Jgetg(J2);
	d0 = Jgetd0(J2);
	nW = lg(W1)-1;
	nV = lg(V)-1;
	nKV = lg(gel(KV,1))-1;
	nZ = lg(gel(V,1))-1;
	pe1 = powiu(p,e1);
	e21 = e2-e1;
	pe21=powiu(p,e21);
	E2 = stoi(e2);
	E21 = stoi(e21);
	NW = utoi(nW);
	NKV = utoi(nKV);
	NZ = utoi(nZ);
	D0 = utoi(d0);

	/* Lift W1 as a subspace of V */
	av1 = avma;
	K = cgetg(nV+nW+1,t_MAT);
	for(j=1;j<=nV;j++) gel(K,j) = gel(V,j);
	for(j=1;j<=nW;j++) gel(K,nV+j) = gel(W1,j);
	K = matkerpadic(K,T,p,e1);
  for(j=1;j<=nW;j++) setlg(gel(K,j),nV+1);
	W1 = FqM_mul(V,K,T,pe2);
	W1 = gerepileupto(av1,W1);

	/* Write matrix K */
	wV = cgetg(nV+1,t_MAT);
	for(j=1;j<=nV;j++)
	{
		col = cgetg(nZ+1,t_COL);
		for(P=1;P<=nZ;P++)
		{
			gel(col,P) = Fq_mul(gcoeff(W1,P,1),gcoeff(V,P,j),T,pe2);
		}
		gel(wV,j) = col;
	}
	KwV = mateqnpadic(wV,T,p,e2);

	av1 = avma;
	K = cgetg(nZ+1,t_MAT);
	for(P=1;P<=nZ;P++)
	{
		col = cgetg(nW*nKV+1,t_COL);
		for(e=1;e<=nKV;e++)
		{
			gel(col,e) = gcoeff(KV,e,P);
			for(j=2;j<=nW;j++)
			{
				gel(col,(j-1)*nKV+e) = Fq_mul(gcoeff(KwV,e,P),gcoeff(W1,P,j),T,pe2);
			}
			gel(K,P) = col;
		}
	}
	r = nZ-(d0+1-g);
	uv = FqM_MinorCompl(K,T,p);
	ABCD = M2ABCD(K,uv);
	Ainv = ZpXQM_inv(gel(ABCD,1),T,p,e2);
	CAinv = FqM_mul(gel(ABCD,3),Ainv,T,pe2);
	AinvB = FqM_mul(Ainv,gel(ABCD,2),T,pe2);
	rho = FqM_mul(CAinv,gel(ABCD,2),T,pe2);
	rho = FpXM_sub(gel(ABCD,4),rho,pe2); /* size nW*nKV-r,(d0+1-g=nZ-r) */
	for(i=1;i<=nW*nKV-r;i++)
	{
		for(j=1;j<=nZ-r;j++)
		{
			gcoeff(rho,i,j) = ZX_Z_divexact(gcoeff(rho,i,j),pe1);
		}
	}

	F = matF(wV,T,p,e21);
	/* Negate */
	for(j=1;j<=nZ;j++)
	{
		for(i=1;i<=nV;i++)
		{
			gcoeff(F,i,j) = FpX_neg(gcoeff(F,i,j),pe21);
		}
	}
	VF = FqM_mul(V,F,T,pe21);
	VFlist = cgetg(nW+1,t_VEC); /* [j] = VF*diag(W[j]) */
	for(j=2;j<=nW;j++)
	{
		gel(VFlist,j) = cgetg(nZ+1,t_MAT);
		for(P=1;P<=nZ;P++)
		{
			gmael(VFlist,j,P) = FqC_Fq_mul(gel(VF,P),gcoeff(W1,P,j),T,pe21);
		}
	}
	gel(VFlist,1) = gen_0;
	/* Now find defs that leave a minor of W fixed */
	sW = gel(FqM_indexrank(W1,T,p),1); /* # = nW */
	cW = VecSmallCompl(sW,nZ);
	av1 = avma;
	Vs = cgetg(nV+1,t_MAT);
	for(j=1;j<=nV;j++)
	{
		col = cgetg(1+nW,t_COL);
		for(i=1;i<=nW;i++)
		{
			gel(col,i) = gcoeff(V,sW[i],j);
		}
		gel(Vs,j) = col;
	}
	V0 = FqM_mul(V,matkerpadic(Vs,T,p,e21),T,pe21); /* subspace of V whose rows in sW are 0 */
	V0 = gerepileupto(av1,V0); /* # = nV-nW = d0 */
	KwVlist = cgetg(d0+1,t_VEC); /* [i] = KwV * diag(V0[i]) */ 
	for(i=1;i<=d0;i++)
	{
		gel(KwVlist,i) = zeromatcopy(nKV,nZ);
		for(j=1;j<=nZ-nW;j++)
		{
			P=cW[j];
			gmael(KwVlist,i,P) = FqC_Fq_mul(gel(KwV,P),gcoeff(V0,P,i),T,pe21);
		}
	}

	M = cgetg(2+d0*nW,t_MAT);
	gel(M,1+d0*nW) = mat2col(rho);

	pending = 0;
  worker = strtofunction("PicLift_worker");
  mt_queue_start(&pt, worker);
	for(k=1; k<=d0 || pending; k++)
	{
		mt_queue_submit(&pt,k,k>d0?NULL:mkvecn(10,NW,NZ,NKV,gel(KwVlist,k),VFlist,uv,AinvB,CAinv,T,pe21));
		done = mt_queue_get(&pt, &workid, &pending);
    if(done)
    {
			for(j=1;j<=nW;j++)
			{
				gel(M,(workid-1)*nW+j) = gel(done,j);
			}
    }
	}
	mt_queue_end(&pt);
	KM = matkerpadic(M,T,p,e21);
	n = lg(KM)-1;
	if(n!=d0+1)	printf("WARNING: dim ker M = %ld (expected %ld)\n",n,d0+1);
	
	if(cmpii(pe21,powiu(l,g+1))<=0)
	{
		printf("Lift by mul\n");
		W = PicLift_RandLift(W1,KM,V0,nZ,nW,d0,T,p,e21,pe1,pe2,pe21);
		W = PicMul(J2,W,pe21,0);
		return mkvec2(gerepileupto(av,W),P0_hint);
	}

	Wlifts = cgetg(g+2,t_VEC);
  K = cgetg(g+2,t_MAT);
	todo_c0 = 1;
	worker = strtofunction("PicLiftTors_worker");
	av1 = av;
	while(1)
	{
		if(todo_c0)
		{
			av = av1;
			/* Find coords of 0 */
  		/* TODO pass this as argument */
			c0 = PicChart(J2,W0,itou(P0));
			nc = lg(c0)-1;
			/* TODO could be NULL */
			/* Find index to dehomogenise */
			red = NULL;
			for(k0=1;k0<=nc;k0++)
			{
				red = gel(c0,k0);
				if(!ZX_is0mod(red,p)) break;
			}
			c0 = FqV_Fq_mul(c0,ZpXQ_inv(red,T,p,e2),T,pe2);
			K0 = utoi(k0);
			todo_c0 = 0;
			av2 = av;
		}
		av = av2;
		/* Find g+1 lifts */
		/* TODO avma = av1; */
		pending = 0;
		mt_queue_start(&pt,worker);
		for(i=1;i<=g+1||pending;i++)
		{
			if(i<=g+1)
			{
				randseed = utoi(pari_rand());
				args = mkvecn(19,J2,W1,l,KM,c0,V0,D0,NW,NZ,K0,T,p,E2,pe2,E21,pe21,pe1,randseed,P0);
				mt_queue_submit(&pt,i,args);
			}
			else mt_queue_submit(&pt,i,NULL);
			done = mt_queue_get(&pt,&workid,&pending);
			if(done)
			{
				gel(Wlifts,workid) = gel(done,1);
				gel(K,workid) = gel(done,2);
			}
		}
		mt_queue_end(&pt);
		Ktors = matkerpadic(K,T,p,e21);
    n = lg(Ktors)-1;
		if(n>1)
		{
			printf("Dim ker tors = %ld, changing charts\n",n);
			P0 = addiu(P0,1);
			todo_c0 = 1;
			continue;
		}
		Ktors = gel(Ktors,1);
		red = gel(Ktors,1);
		for(i=2;i<=g+1;i++)
		{
			red = FpX_add(red,gel(Ktors,i),pe2);
		}
		if(ZX_is0mod(red,p))
		{
			printf("Sum of Ktors is zero!\n");
			continue;
		}
		Ktors = FqC_Fq_mul(Ktors,ZpXQ_inv(red,T,p,e2),T,pe2);
		for(i=1;i<=g+1;i++) gel(Wlifts,i) = FqM_Fq_mul(gel(Wlifts,i),gel(Ktors,i),T,pe2);
		W = gel(Wlifts,1);
		for(i=2;i<=g+1;i++) W = FpXM_add(W,gel(Wlifts,i),pe2);
		if(P0_hint == NULL)
		{ /* The chart might not be a diffeo, need to check we got all right */
			if(!PicIsZero(J2,PicMul(J2,W,l,0)))
			{
				printf("Not actually l-torsion!!! Changing charts\n");
      	P0 = addiu(P0,1);
      	todo_c0 = 1;
      	continue;
			}
		}
		return gerepilecopy(av,mkvec2(W,P0));
	}
}

GEN PicLiftTors(GEN J, GEN W, long eini, GEN l)
{
	pari_sp av=avma;
	ulong e,efin,e2;
	GEN Je,p,P0;

	efin = Jgete(J);
	e = eini;
	p = Jgetp(J);
	P0 = NULL;
	while(e<efin)
	{
		e2 = 2*e;
		if(e2>efin) e2 = efin;
		pari_printf("Lifting from prec O(%Ps^%lu) to O(%Ps^%lu)\n",p,e,p,e2);
		Je = e2<efin ? PicRed(J,e2) : J;
		W = PicLiftTors_2(Je,W,e,l,P0);
		P0 = gel(W,2);
		if(P0)
		{
			W = gerepileupto(av,W);
			P0 = gel(W,2);
			W = gel(W,1);
		}
		else
		{
			W = gerepileupto(av,gel(W,1));
		}
		e = e2;
	}
	return gerepilecopy(av,W);
}
