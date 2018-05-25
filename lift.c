#include "linalg.h"
#include "pic.h"

GEN PicRed(GEN J, ulong e1); /* TODO in pic.c ? */

GEN PicLiftTors_2(GEN J2, GEN W1, ulong e1, GEN l)
{
	pari_sp av1,av=avma;
	GEN V,KV,T,p,pe1,pe2,pe21;
	GEN col,K,wV,KwV,Ainv,AinvB,CAinv,rho,ABCD,uv;
	ulong g,d0,nW,nV,nKV,nZ,e2,r;
	ulong e,j,P;

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
	pe21=powiu(p,e2-e1);
	printf("a");

	av1 = avma;

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
	KwV = gerepileupto(av1,KwV);
	printf("b");

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
	printf("d");
	ABCD = M2ABCD(K,uv);
	Ainv = ZpXQM_inv(gel(ABCD,1),T,p,e2);
	CAinv = FqM_mul(gel(ABCD,3),Ainv,T,pe2);
	AinvB = FqM_mul(Ainv,gel(ABCD,2),T,pe2);
	rho = FqM_mul(CAinv,gel(ABCD,2),T,pe2);
	printf("e");
	rho = RgM_sub(gel(ABCD,4),rho);

}
