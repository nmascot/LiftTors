#include "linalg.h"
#include "exp.h"
#include "pic.h"
#include "hyper.h"

GEN FindSuppl(GEN V, GEN W, GEN T, GEN p, GEN pe, int sq, ulong nV2)
/* /!\ Shallow */
{
	pari_sp av1,av = avma;
	GEN S,v1,v2,col;
	ulong i,j,nW,nV,nS,nZ;
	nW = lg(W)-1;
	nV = sq? nV2 : lg(V)-1;
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

GEN detratio(GEN K, GEN T, GEN p, long e, GEN pe)
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
	M = ZpXQM_det(M,T,p,e);
	return gerepileupto(av,M);
}

GEN PicNorm(GEN J, GEN F, GEN WE)
{
	pari_sp av = avma;
	ulong g,d0,nS1,nV2,nZ;
	ulong i,j;
	long e;
	GEN V,T,p,pe;
	GEN WEV,V1,V2,M1,M2,M;

	g = Jgetg(J);
	d0 = Jgetd0(J);
	V = JgetV(J);
	JgetTpe(J,&T,&pe,&p,&e);
	d0 = Jgetd0(J);
	g = Jgetg(J);
	nS1 = d0;
	nV2 = 6*d0+1-g;
	nZ = lg(gel(V,1))-1;

	WEV = DivAdd(V,WE,5*d0+1-g,T,p,e,pe,0);
	V1 = FindSuppl(V,WE,T,p,pe,0,0);
	V2 = FindSuppl(V,WEV,T,p,pe,1,6*d0+1-g);

	M = cgetg(nS1+nV2+1,t_MAT);
	for(j=1;j<=nS1;j++) gel(M,j) = gel(V1,j);
	for(j=1;j<=nV2;j++) gel(M,nS1+j) = gel(V2,j);
	M1 = detratio(matkerpadic(M,T,p,e),T,p,e,pe);
	if(ZX_is0mod(M1,p)) pari_err(e_MISC,"D intersects D0");

	for(j=1;j<=nS1;j++)
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
	Fq1 = mkpoln(1,gen_1);
	setvarn(Fq1,varn(T));
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

GEN PicTorsRels(GEN J, GEN Wtors, GEN l, ulong excess)
{
	pari_sp av = avma;
	ulong ntors,ntest,n,i;
	GEN T,p,q,q1,z,m;
	GEN W0,C,Wtest,R,x,zn;
	struct pari_mt pt;
	GEN worker,done;
	long pending,j,workid;

	if(Jgete(J)>1) pari_err(e_IMPL,"case e>1");
	T = JgetT(J);
	p = Jgetp(J);
	q = powiu(p,degree(T));
	q1 = subii(q,gen_1);
	if(!gequal0(modii(q1,l))) pari_err(e_MISC,"No l-th roots of 1");
	m = divii(q1,l);
	z = gener_FpXQ(T,p,NULL);
	z = Fq_pow(z,m,T,p);
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
			x = Fq_pow(gcoeff(R,i,j),m,T,p);
			n = 0;
			zn = gen_1;
			while(!gequal(x,zn))
			{
				zn= Fq_mul(zn,z,T,p);
				n++;
			}
			gcoeff(R,i,j) = utoi(n);
		}
	}
	return gerepileupto(av,FpM_ker(R,l));
}
