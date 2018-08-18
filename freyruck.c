#include "linalg.h"
#include "exp.h"
#include "pic.h"
#include "hyper.h"

GEN FindSuppl(GEN W, ulong nS, GEN V, GEN Vbis, GEN T, GEN p, GEN pe)
/* /!\ Shallow */
/* Look for suppl of W of dim nS in V*Vbis, or in V if Vbis=NULL */
{
	pari_sp av1,av = avma;
	GEN S,v1,v2,col;
	ulong i,j,nW,nZ;
	nW = lg(W)-1;
	nZ = lg(gel(V,1))-1;
	S = cgetg(nS+nW+1,t_MAT);
	av1 = avma;
	do
	{
		avma = av1;
		if(Vbis)
		{
			for(j=1;j<=nS;j++)
			{
				col = cgetg(nZ+1,t_COL);
				v1 = RandVec_padic(V,T,p,pe);
				v2 = RandVec_padic(Vbis,T,p,pe);
				for(i=1;i<=nZ;i++) gel(col,i) = Fq_mul(gel(v1,i),gel(v2,i),T,pe);
				gel(S,j) = col;
			}
		}
		else
		{
			for(j=1;j<=nS;j++) gel(S,j) = RandVec_padic(V,T,p,pe);
		}
		for(j=1;j<=nW;j++) gel(S,j+nS) = gel(W,j);
	}while(FqM_rank(S,T,p)<nS+nW);
	if(Vbis) S = gerepilecopy(av,S);
	return S;
}

GEN detratio(GEN K, GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN K1,K2,col1,col2;
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
	if(e==1)
		return gerepileupto(av,Fq_div(FqM_det(K2,T,p),FqM_det(K1,T,p),T,p));
	return gerepileupto(av,Fq_mul(ZpXQM_det(K2,T,p,e),ZpXQ_inv(ZpXQM_det(K1,T,p,e),T,p,e),T,pe));
}

GEN PicNorm(GEN J, GEN F, GEN WE)
{
	pari_sp av = avma;
	ulong g,d0,nS,nV5,nZ;
	ulong i,j;
	long e;
	GEN V,V3,T,p,pe;
	GEN WEV3,V1,V2,M1,M2,M;

	g = Jgetg(J);
	d0 = Jgetd0(J);
	V = JgetV(J);
	V3 = JgetV3(J);
	JgetTpe(J,&T,&pe,&p,&e);
	d0 = Jgetd0(J);
	g = Jgetg(J);
	nS = lg(V)-lg(WE); /* codim WE in V = deg E */
	nV5 = 5*d0+1-g;
	nZ = lg(gel(V,1))-1;

	WEV3 = DivAdd(V3,WE,nV5-nS,T,p,e,pe,0);

	V1 = FindSuppl(WE,nS,V,NULL,T,p,pe);
	V2 = FindSuppl(WEV3,nS,V,V3,T,p,pe);

	M = cgetg(nS+nV5+1,t_MAT);
	for(j=1;j<=nS;j++) gel(M,j) = gel(V1,j);
	for(j=1;j<=nV5;j++) gel(M,nS+j) = gel(V2,j);
	M1 = detratio(matkerpadic(M,T,p,e),T,p,e,pe);
	if(ZX_is0mod(M1,p)) pari_err(e_MISC,"D intersects D0");

	for(j=1;j<=nS;j++)
	{
		for(i=1;i<=nZ;i++)
		{
			gcoeff(M,i,j) = Fq_mul(gcoeff(M,i,j),gel(F,i),T,pe);
		}
	}
	M2 = detratio(matkerpadic(M,T,p,e),T,p,e,pe);
	if(ZX_is0mod(M2,p)) pari_err(e_MISC,"F has zeros on D");
	
	return gerepileupto(av,ZpXQ_div(M2,M1,T,pe,p,e));
}

GEN PicFreyRuckMulti(GEN J, GEN Wtors, GEN l, GEN Wtest, GEN W0, GEN C)
/* Pair the l-tors pt Wtors against the pts in Wtest */ 
{
	pari_sp av = avma;
	GEN WtorsM,Fq1,H,col,WA,WB,res,s;
	GEN T,p,pe,KV;
	long e;
	ulong nC,ntest;
	ulong c,d,i,j;
	
	JgetTpe(J,&T,&pe,&p,&e);
	Fq1 = GetFq1(T);
	KV = JgetKV(J);
	nC = lg(C);
	ntest = lg(Wtest);
	WtorsM = cgetg(nC,t_VEC);
	gel(WtorsM,1) = Wtors;
	/*H = cgetg(nC,at_MAT);*/
	H = zeromatcopy(ntest-1,nC-1);
	gel(H,1) = cgetg(ntest,t_COL);
	for(d=1;d<ntest;d++) gcoeff(H,d,1) = Fq1;
	for(c=2;c<nC;c++)
	{
		i = gmael(C,c,2)[1];
		j = gmael(C,c,2)[2];
		WA = i?gel(WtorsM,i):JgetW0(J);
		WB = j?gel(WtorsM,j):JgetW0(J);
		res = PicChord(J,WA,WB,3);
		gel(WtorsM,c) = gel(res,1);
		s = gel(res,2);
		col = cgetg(ntest,t_COL);
		for(d=1;d<ntest;d++)
		{
			gel(col,d) = Fq_mul(PicNorm(J,s,gel(Wtest,d)),gcoeff(H,d,i),T,pe);
			if(j) gel(col,d) = Fq_mul(gel(col,d),gcoeff(H,d,j),T,pe);
			gel(col,d) = ZpXQ_inv(gel(col,d),T,p,e);
		}
		gel(H,c) = FqC_Fq_mul(col,PicNorm(J,s,W0),T,pe);
	}
	s = DivSub(JgetW0(J),gel(WtorsM,nC-1),KV,1,T,p,e,pe,2);
	s = gel(s,1);
	col = gel(H,nC-1);
	for(d=1;d<ntest;d++)
	{
		gel(col,d) = Fq_mul(gel(col,d),PicNorm(J,s,gel(Wtest,d)),T,pe);
	}
	col = FqC_Fq_mul(col,ZpXQ_inv(PicNorm(J,s,W0),T,p,e),T,pe);
	return gerepileupto(av,col);
}

GEN Fq_zeta_l(GEN T, GEN p, GEN l)
{
	pari_sp av = avma;
	GEN q,q1,m,z;
	q = powiu(p,degree(T));
  q1 = subii(q,gen_1);
  if(!gequal0(modii(q1,l))) pari_err(e_MISC,"No l-th roots of 1");
  m = divii(q1,l);
  z = gener_FpXQ(T,p,NULL);
  z = Fq_pow(z,m,T,p);
	return gerepileupto(av,z);
}

GEN Fq_mu_l_log(GEN x, GEN z, GEN T, GEN p, GEN l)
{
	pari_sp av = avma;
	ulong n = 0;
	GEN q,q1,m,zn,y;
  q = powiu(p,degree(T));
  q1 = subii(q,gen_1);
  m = divii(q1,l);
	y = Fq_pow(x,m,T,p);
  zn = gen_1;
  while(!gequal(y,zn))
  {
    zn = Fq_mul(zn,z,T,p);
    n++;
  }
  return gerepileupto(av,utoi(n));
}

GEN PicTorsRels(GEN J, GEN Wtors, GEN l, ulong excess)
{
	pari_sp av = avma;
	ulong ntors,ntest,i;
	GEN T,p,z;
	GEN W0,C,Wtest,R;
	struct pari_mt pt;
	GEN worker,done;
	long pending,j,workid;

	if(Jgete(J)>1) pari_err(e_IMPL,"case e>1");
	T = JgetT(J);
	p = Jgetp(J);
	z = Fq_zeta_l(T,p,l);
	W0 = JgetW0(J);
	W0 = PicChord(J,W0,W0,1);
	C = AddChain(l,0);
	ntors = lg(Wtors)-1;
	ntest = ntors+excess;
	Wtest = cgetg(ntest+1,t_VEC);
	for(i=1;i<=ntest;i++) gel(Wtest,i) = PicChord(J,PicRand0(J),PicRand0(J),1); /* TODO sufficient? */
	R = cgetg(ntors+1,t_MAT);
	pending = 0;
	worker = strtofunction("PicFreyRuckMulti");
	mt_queue_start(&pt,worker);
	for(j=1;j<=ntors||pending;j++)
	{
		mt_queue_submit(&pt,j,j<=ntors?mkvecn(6,J,gel(Wtors,j),l,Wtest,W0,C):NULL);
		done = mt_queue_get(&pt,&workid,&pending);
		if(done) gel(R,workid) = done;
	}
	mt_queue_end(&pt);
	for(j=1;j<=ntors;j++)
	{
		for(i=1;i<=ntest;i++)
		{
			gcoeff(R,i,j) = Fq_mu_l_log(gcoeff(R,i,j),z,T,p,l);
		}
	}
	return gerepileupto(av,FpM_ker(R,l));
}
