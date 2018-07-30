#include "linalg.h"
#include "pic.h"

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

GEN PicLift_RandLift(GEN W1, GEN KM, GEN V0, ulong nZ, ulong nW, ulong d0, GEN T, GEN p, long e21, GEN pe1, GEN pe2, GEN pe21)
{
	pari_sp av = avma;
	GEN K,red,W;
	ulong n,j,k,P;
	/* Find a random solution to the inhomogeneous system */
  do
  {
    K = RandVec_padic(KM,T,p,pe21);
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
}

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
	KM = matkerpadic_hint(M,T,p,e21,pe21,d0+1);
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
