# include "linalg.h"
#include "pic.h"

GEN FindSuppl(GEN V, GEN W, GEN T, GEN p, GEN pe, int sq)
/* /!\ Shallow */
{
	pari_sp av1,av = avma;
	GEN S,v1,v2,col;
	ulong i,j,nW,nV,nS,nZ;
	nW = lg(W)-1;
	nV = lg(V)-1;
	if(sq) nV = 2*nV+g-1;
	nZ = lg(gel(V,1))-1;
	nS = nV-nW;
	S = cgetg(nV+1,t_MAT);
	av1 = avma;
	do
	{
		avma = av1;
		if(sq)
		{
			for(j=1;j<=nS;j++)
			{
				col = cgetg(nZ+1,t_COL);
				v1 = RandVec_padic(V,T,p,pe);
				v2 = RandVec_padic(V,T,p,pe);
				for(i=1;i<=nZ;i++) gel(col,i) = Fq_mul(gel(v1,i),gel(v2,i),T,pe);
				gel(S,j) = col;
			}
		}
		else
		{
			for(j=1;j<=nS;j++) gel(S,j) = RandVec_padic(V,T,p,pe);
		}
		for(j=1;j<=nW;j++) gel(S,j+nS) = gel(W,j);
	}while(FqM_rank(S,T,p)<nV);
	if(sq) S = gerepilecopy(av,S);
	return S;
}

GEN detratio(GEN K, GEN T, GEN p, ulong e, GEN pe)
{
	pari_sp av = avma;
	GEN K1,K2,col1,col2,M;
	ulong d0,i,j;
	d0 = lg(K)-1;
	K1 = cgetg(d0+1,t_MAT);
  K2 = cgetg(d0+1,t_MAT);
  for(j=1;j<=d0;j++)
  {
    col1 = cgetg(d0+1,t_COL);
    col2 = cgetg(d0+1,t_COL);
    for(i=1;i<=d0;i++)
    {
      gel(col1,i) = gcoeff(K,i,j);
      gel(col2,i) = gcoeff(K,d0+i,j);
    }
    gel(K1,j) = col1;
    gel(K2,j) = col2;
  }
  M = FqM_mul(K2,ZpXQM_inv(K1,T,p,e),T,pe);
  M = ZpXQM_det(M1,T,p,e);
	return gerepileupto(av,M);
}

GEN PicNorm(GEN J, GEN F, GEN E)
{
	pari_sp av = avma;
	ulong g,d0,e,nS1,nV2,nZ;
	ulong i,j;
	GEN V,T,p,pe;
	GEN WEV,V1,V2,col1,col2,M1,M2,K,K1,K2,M;

	g = Jgetg(J);
	d0 = Jgetd0(J);
	V = JgetV(J);
	T = JgetT(J);
	p = Jgetp(J);
	e = Jgete(J);
	pe = Jgetpe(J);


	WEV = DivAdd(V,WE,5*d0+1-g,T,p,e,pe,0);
	V1 = FindSuppl(V,WE,T,p,pe,0);
	V2 = FindSuppl(V2,WEV,T,p,pe,1);

	M = cgetg(nS1+nV2+1,t_MAT);
	for(j=1;j<=nS1;j++) gel(M,j) = gel(V1,J);
	for(j=1;j<=nV2;j++) gel(M,nV2+j) = gel(V2,J);
	K = matkerpadic(M,T,p,e);
	K1 = cgetg(d0+1,t_MAT);
	K2 = cgetg(d0+1,t_MAT);
	for(j=1;j<=d0;j++)
	{
		col1 = cgetg(d0+1,t_COL);
		col2 = cgetg(d0+1,t_COL);
		for(i=1;i<=d0;i++)
		{
			gel(col1,i) = gcoeff(K,i,j);
			gel(col2,i) = gcoeff(K,d0+i,j);
		}
		gel(K1,j) = col1;
		gel(K2,j) = col2;
	}
	M1 = FqM_mul(K2,ZpXQM_inv(K1,T,p,e),T,pe);
	M1 = ZpXQM_det(M1,T,p,e);

	for(j=1;j<=nS1;j++)
	{
		for(i=1;i<=nZ;i++)
		{
			gcoeff(M,i,j) = Fq_mul(gcoeff(M,i,j),gel(F,i),T,pe);
		}
	}
	K = matkerpadic(M,T,p,e);
  K1 = cgetg(d0+1,t_MAT);
  K2 = cgetg(d0+1,t_MAT);
  for(j=1;j<=d0;j++)
  {
    col1 = cgetg(d0+1,t_COL);
    col2 = cgetg(d0+1,t_COL);
    for(i=1;i<=d0;i++)
    {
      gel(col1,i) = gcoeff(K,i,j);
      gel(col2,i) = gcoeff(K,d0+i,j);
    }
    gel(K1,j) = col1;
    gel(K2,j) = col2;
  }
  M2 = FqM_mul(K2,ZPXQM_inv(K1,T,p,e),T,pe);
	M2 = ZpXQM_det(M2,T,p,e);
	
	return gerepileupto(av,ZpXQ_div(M2,M1,T,pe,p,e));
}
