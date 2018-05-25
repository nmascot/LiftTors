#include "linalg.h"
#include "pic.h"

GEN PicRed(GEN J, ulong e1); /* TODO in pic.c ? */

GEN PicLiftTors_2(GEN J2, GEN W1, ulong e1, GEN l)
{
	pari_sp av1,av=avma;
	GEN V,KV,T,p,pe1,pe2,pe21;
	GEN K,wV,KwV,Ainv,AinvB,CAinv,rho,ABCD,uv;
	ulong g,d0,nW,nV,nKV,nZ,e2,r;
	ulong e,j,P;


	/*TODO test e2>e1 */

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
	pe2 = Jgetpe(J2);
	pe21=powiu(p,e2-e1);

	av1 = avma;

	wV = cgetg(nW+1,t_MAT);
	for(j=1;j<=nW;j++)
	{
		gel(wV,j) = cgetg(nZ+1,t_COL);
		for(P=1;P<=nZ;P++)
		{
			gcoeff(wV,P,j) = Fq_mul(gcoeff(W1,P,1),gcoeff(V,P,j),T,pe2);
		}
	}
	KwV = mateqnpadic(wV,T,p,e2);
	KwV = gerepileupto(av1,KwV);

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
	uv = FindMinorCompl(K,T,p); /* TODO */
	ABCD = M2ABCD(K,uv);
	Ainv = matinvzq(gel(ABCD,1),T,pe2);
	CAinv = FqM_mul(gel(ABCD,3),Ainv,T,pe2);
	AinvB = FqM_mul(Ainv,gel(ABCD,2),T,pe2);
	rho = FqM_mul(CAinv,gel(ABCD,2),T,pe2);
	rho = FpXM_sub(gel(ABCD,4),rho,pe2);
	/*for(j=1;j<=nZ-r;j++)
	{
		for(i=1;i<=nW*nKV;i++)
		{
	*/


}
