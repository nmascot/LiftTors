#include "linalg.h"
#include "pic.h"

GEN PicLift_worker(GEN NW, GEN NZ, GEN NKV, GEN KwVk, GEN VFlist, GEN uv, GEN AinvB, GEN CAinv, GEN T, GEN pe21)
{
	GEN M,dK,dKi,ABCD;
	ulong nW,nZ,nKV;
	ulong j,P,i,e;

	nW = itou(NW);
	nZ = itou(NZ);
	nKV = itou(NKV);
	M = cgetg(nW+1,t_MAT);
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
	return M;
}

GEN PicLiftTors_worker(GEN J, GEN W1, GEN l, GEN KM, GEN c0, GEN V0, GEN D0, GEN NW, GEN NZ, GEN K0, GEN T, GEN p, GEN E2, GEN pe2, GEN E21, GEN pe21, GEN pe1, GEN randseed)
{
	pari_sp av = avma;
	GEN K,red,W,c;
	ulong d0,nW,nZ,nc,k0,n,i,j,k,P;
	long e2,e21;
	setrand(randseed);
	d0 = itou(D0);
	nW = itou(NW);
	nZ = itou(NZ);
	nc = lg(c0)-1;
	k0 = itou(K0);
	e2 = itos(E2);
	e21 = itos(E21);
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
  /* Mul by l, get coordinates, and compare them to those of W0 */
  c = PicChart(J,PicMul(J,W,l,0));
  c = FqV_Fq_mul(c,ZpXQ_inv(gel(c,k0),T,p,e2),T,pe2);
  for(i=1;i<=nc;i++)
  {
    gel(c,i) = ZX_Z_divexact(FpX_sub(gel(c,i),gel(c0,i),pe2),pe1);
  }
	return gcopy(mkvec2(W,c));
	return gerepilecopy(av,mkvec2(W,c));
}

GEN PicLiftTors_2(GEN J2, GEN W1, long e1, GEN l)
{
	pari_sp av1,av=avma;
	GEN V,KV,W0,T,p,pe1,pe2,pe21;
	GEN col,K,wV,KwV,Ainv,AinvB,CAinv,rho,ABCD,uv;
	GEN F,VF,VFlist,sW,cW,Vs,V0,KwVlist,M,KM,worker,done;
	GEN c0,Wlifts,W,red;
	GEN NW,NZ,NKV,D0,K0,E2,E21,randseed,args;
	ulong g,d0,nW,nV,nKV,nZ,nc,r;
	long e2,e21,pending,workid;
	ulong e,i,j,k,P,k0;

	struct pari_mt pt;
  static entree ep_worker_1={"PicLift_worker",0,(void*)PicLift_worker,1,"GGGGGGGGGG",""};
  static entree ep_worker_2={"PicLiftTors_worker",0,(void*)PicLiftTors_worker,1,"GGGGGGGGGGGGGGGGGG",""};
	static int worker_1 = 0;
	static int worker_2 = 0;
	if(worker_1==0)
	{
		pari_add_function(&ep_worker_1);
		worker_1 = 1;
	}
	if(worker_2==0)
	{
		pari_add_function(&ep_worker_2);
		worker_2 = 1;
	}

	JgetTpe(J2,&T,&pe2,&p,&e2);
	if(e2<=e1)
	{
		return FpXM_red(W1,pe2);
	}

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
	printf("Ker...\n");
	KM = matkerpadic_hint(M,T,p,e21,pe21,d0+1);
	printf("Dim ker M = %ld\n",lg(KM)-1);
	
	/* Find coords of 0 */
	/* TODO pass this as argument */
	c0 = PicChart(J2,W0);
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

	Wlifts = cgetg(g+2,t_VEC);
	K = cgetg(g+2,t_MAT);
	worker = strtofunction("PicLiftTors_worker");
	av1 = avma;
	do
	{
		/* Find g+1 lifts */
		avma = av1;
		pending = 0;
		mt_queue_start(&pt,worker);
		for(i=1;i<=g+1||pending;i++)
		{
			if(i<=g+1)
			{
				randseed = utoi(pari_rand());
				args = mkvecn(18,J2,W1,l,KM,c0,V0,D0,NW,NZ,K0,T,p,E2,pe2,E21,pe21,pe1,randseed);
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
		for(j=2;j<=g+1;j++)
		{
			for(i=1;i<=nc;i++)
			{
				gcoeff(K,i,j) = FpX_sub(gcoeff(K,i,j),gcoeff(K,i,1),pe21);
			}
		}
		K = matkerpadic(K,T,p,e21);
		i = lg(K)-1;
		printf("Dim ker tors = %ld\n",i);
		/* Find a col in K with 1st entry invertible */
		k = 0;
		for(j=1;j<=i;j++)
		{
			red = gcoeff(K,1,j);
			if(!ZX_is0mod(red,p))
			{
				k = j;
				K = gel(K,k);
				break;
			}
		}
	} while(k==0);
	K = FqC_Fq_mul(K,ZpXQ_inv(red,T,p,e21),T,pe21);
	W = gel(Wlifts,1);
	for(i=2;i<=g+1;i++)
	{
		W = FpXM_add(W,FqM_Fq_mul(FpXM_sub(gel(Wlifts,i),gel(Wlifts,1),pe2),gel(K,i),T,pe2),pe2);
	}
	return gerepileupto(av,W);
}

GEN PicLiftTors(GEN J, GEN W, long eini, GEN l)
{
	pari_sp av=avma;
	ulong e,efin,e2;
	GEN Je,p;

	efin = Jgete(J);
	e = eini;
	p = Jgetp(J);
	while(e<efin)
	{
		e2 = 2*e;
		if(e2>efin) e2 = efin;
		pari_printf("Lifting from prec O(%Ps^%lu) to O(%Ps^%lu)\n",p,e,p,e2);
		Je = e2<efin ? PicRed(J,e2) : J;
		W = gerepileupto(av,PicLiftTors_2(Je,W,e,l));
		e = e2;
	}
	return W;
}
