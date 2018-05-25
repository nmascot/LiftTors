#include "linalg.h"
#include "pic.h"

GEN PicRed(GEN J, ulong e1); /* TODO in pic.c ? */

GEN PicLiftTors_2(GEN J2, GEN W1, ulong e1, GEN l)
{
	pari_sp av1,av=avma;
	GEN V,KV,T,p,pe1,pe2,pe21;
	GEN col,K,wV,KwV,Ainv,AinvB,CAinv,rho,ABCD,uv;
	GEN F,VF,VFlist,sW,cW,Vs,V0,KwVlist,M,KM,dK,dKi;
	ulong g,d0,nW,nV,nKV,nZ,e2,e21,r;
	ulong e,i,j,k,n,P;

	e2 = Jgete(J2);
	pe2 = Jgetpe(J2);
	if(e2<=e1)
	{
		return FpXM_red(W1,pe2);
	}

	V = JgetV(J2);
	KV = JgetKV(J2);
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
	pari_printf("%Ps\n",pe21);

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
		gel(K,P) = cgetg(nW*nKV+1,t_COL);
		for(j=1;j<=nW;j++)
		{
			if(j==1)
			{
				for(e=1;e<=nKV;e++)
				{
					gcoeff(K,e,P) = gcoeff(KV,e,P);
				}
			}
			else
			{
				for(e=1;e<=nKV;e++)
				{
					gcoeff(K,(j-1)*nKV+e,P) = Fq_mul(gcoeff(W1,P,j),gcoeff(KwV,e,P),T,pe2);
				}
			}
		}
	}
	r = nZ-(d0+1-g);
	uv = FqM_MinorCompl(K,T,p);
	ABCD = M2ABCD(K,uv);
	Ainv = ZpXQM_inv(gel(ABCD,1),T,p,e2);
	CAinv = FqM_mul(gel(ABCD,3),Ainv,T,pe2);
	AinvB = FqM_mul(Ainv,gel(ABCD,2),T,pe2);
	rho = FqM_mul(CAinv,gel(ABCD,2),T,pe2);
	rho = RgM_sub(gel(ABCD,4),rho); /* size d0+1-g,nW*nKV-r */
	rho = gdiv(rho,pe1);

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
	printf("b");
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
	V0 = FqM_mul(V,matkerpadic(Vs,T,p,e21),T,pe21); /* subspace f V whose rows in sW are 0 */
	V0 = gerepileupto(av1,V0); /* # = nV-nW = d0 */
	printf("c");
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

	printf("d");
	M = cgetg(2+d0*nW,t_MAT);
	gel(M,1+d0*nW) = mat2col(rho);

	n = 0;
	for(j=1;j<=nW;j++)
	{
		for(k=1;k<=d0;k++)
		{
			dK = zeromatcopy(nW*nKV,nZ);
			if(j==1)
			{
				for(i=2;i<=nW;i++)
				{
					dKi = FqM_mul(gel(KwVlist,k),gel(VFlist,i),T,pe21);
					for(P=1;P<=nZ;P++)
					{
						for(e=1;e<=nKV;e++)
						{
							gcoeff(dK,(i-1)*nKV+e,P) = gcoeff(dKi,e,P);
						}
					}
				}
			}
			else
			{
				for(P=1;P<=nZ;P++)
				{
					for(e=1;e<=nKV;e++)
					{
						gcoeff(dK,(j-1)*nKV+e,P) = gcoeff(gel(KwVlist,k),e,P);
					}
				}
			}
			ABCD = M2ABCD(dK,uv);
			dK = FpXM_sub(FqM_mul(gel(ABCD,1),AinvB,T,pe21),gel(ABCD,2),pe21);
			dK = FqM_mul(CAinv,dK,T,pe21);
			dK = FpXM_add(gel(ABCD,4),dK,pe21);
			dK = FpXM_sub(dK,FqM_mul(gel(ABCD,3),AinvB,T,pe21),pe21);
			gel(M,++n) = mat2col(dK);
		}
	}
	printf("f");
	KM = matkerpadic(M,T,p,e21); /* TODO accel */ /* TODO varn */
	printf("Dim ker M=%ld",lg(KM)-1);
	return gerepileupto(av,KM);
}
