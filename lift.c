#include "linalg.h"
#include "pic.h"

GEN PicLiftTors_2(GEN J2, GEN W1, ulong e1, GEN l)
{
	pari_sp av1,av2,av=avma;
	GEN V,KV,W0,T,p,pe1,pe2,pe21;
	GEN col,K,wV,KwV,Ainv,AinvB,CAinv,rho,ABCD,uv;
	GEN F,VF,VFlist,sW,cW,Vs,V0,KwVlist,M,KM,dK,dKi;
	GEN c0,c,Wlifts,W,red;
	ulong g,d0,nW,nV,nKV,nZ,nc,e2,e21,r;
	ulong e,i,j,k,n,P,k0;

	e2 = Jgete(J2);
	pe2 = Jgetpe(J2);
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
	T = JgetT(J2);
	p = Jgetp(J2);
	pe1 = powiu(p,e1);
	e21 = e2-e1;
	pe21=powiu(p,e21);

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
	printf("a1");
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
	pari_printf("\nsW=%Ps\ncW=%Ps\nV0=%Ps\n",sW,cW,V0);
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

	n = 0;
	for(j=1;j<=nW;j++)
	{
		for(k=1;k<=d0;k++)
		{
			if(j==1)
			{
				dK = cgetg(nZ+1,t_MAT);
				for(P=1;P<=nZ;P++) gel(dK,P) = cgetg(nKV*(nW-1)+1,t_COL);
				for(i=2;i<=nW;i++)
				{
					dKi = FqM_mul(gel(KwVlist,k),gel(VFlist,i),T,pe21);
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
				ABCD = M2ABCD_1block(gel(KwVlist,k),(j-1)*nKV,0,uv);
			}
			dK = FpXM_sub(FqM_mul(gel(ABCD,1),AinvB,T,pe21),gel(ABCD,2),pe21);
			dK = FqM_mul(CAinv,dK,T,pe21);
			dK = FpXM_add(gel(ABCD,4),dK,pe21);
			dK = FpXM_sub(dK,FqM_mul(gel(ABCD,3),AinvB,T,pe21),pe21);
			gel(M,++n) = mat2col(dK);
		}
	}
	KM = matkerpadic(M,T,p,e21); /* TODO accel */ /* TODO varn */
	printf("Dim ker M=%ld",lg(KM)-1);
	
	/* FInd coords of 0 */
	c0 = PicChart(J2,W0);
	nc = lg(c0)-1;
	/* TODO could be NULL */
	/* Find index to dehomogenise */
	red = NULL;
	for(k0=1;k0<=nc;k0++)
	{
		red = gel(c0,k0);
		if(!FpX_is0modp(red,p)) break;
	}
	c0 = FqV_Fq_mul(c0,ZpXQ_inv(red,T,p,e2),T,pe2);

	Wlifts = cgetg(g+2,t_VEC);
	av1 = avma;
	do
	{
		/* Find g+1 lifts */
		avma = av1;
		for(i=1;i<=g+1;i++)
		{
			av2 = avma;
			do
			{
				avma = av2;
				K = RandVec_padic(KM,T,p,pe21);
				red = gel(K,1+d0*nW);
			}while(FpX_is0modp(red,p));
			red = ZpXQ_inv(red,T,p,e21);
			setlg(K,d0*nW+1);
			K = FqV_Fq_mul(K,red,T,pe21);
			n = 0;
			W = zeromatcopy(nZ,nW);
			for(j=1;j<=nW;j++)
			{
				for(k=1;k<=d0;k++)
				{
					n++;
					for(P=1;P<=nZ;P++)
					{
						gcoeff(W,P,j) = FpX_add(gcoeff(W,P,j),Fq_mul(gel(K,n),gcoeff(V0,P,k),T,pe21),pe21);
					}
				}
			}
			gel(Wlifts,i) = FpXM_add(W1,ZXM_Z_mul(W,pe1),pe2);
		}
		/* Find a combination of them which is l-torsion */
		K = cgetg(g+2,t_MAT);
		for(j=1;j<=g;j++)
		{
			c = PicChart(J2,PicMul(J2,gel(Wlifts,1),l,0));
			c = FqV_Fq_mul(c,ZpXQ_inv(gel(c,k),T,p,e2),T,pe2);
			for(i=1;i<=nc;i++)
			{
				gel(c,i) = ZX_Z_divexact(FpX_sub(gel(c,i),gel(c0,i),pe2),pe1);
			}
			if(j>1)
			{
				for(i=1;i<=nc;i++)
				{
					gel(K,j) = FpX_sub(gel(c,i),gcoeff(K,i,1),pe21);
				}
			}
			gel(K,i) = c;
		}
		K = matkerpadic(K,T,p,e21);
		i = lg(K)-1;
		printf("Dim ker tors = %ld",i);
		/* Find a col in K with 1st entry invertible */
		k = 0;
		for(j=1;j<=i;j++)
		{
			red = gcoeff(K,1,j);
			if(!FpX_is0modp(red,p))
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

GEN PicLiftTors(GEN J, GEN W, GEN l, ulong eini)
{
	pari_sp av=avma;
	ulong e,efin,e2;
	GEN Je;

	efin = Jgete(J);
	e = eini;
	while(e<efin)
	{
		e2 = 2*e;
		if(e2>efin) e2 = efin;
		Je = e2<efin ? PicRed(J,e2) : J;
		W = gerepileupto(av,PicLiftTors_2(Je,W,e,l));
		e = e2;
	}
	return W;
}
