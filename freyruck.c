#include "linalg.h"
#include "exp.h"
#include "pic.h"

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
				v1 = RandVec_1(V,pe);
				v2 = RandVec_1(Vbis,pe);
				for(i=1;i<=nZ;i++) gel(col,i) = Fq_mul(gel(v1,i),gel(v2,i),T,pe);
				gel(S,j) = col;
			}
		}
		else
		{
			for(j=1;j<=nS;j++) gel(S,j) = RandVec_1(V,pe);
		}
		for(j=1;j<=nW;j++) gel(S,j+nS) = gel(W,j);
	}while(FqM_rank(S,T,p)<nS+nW);
	if(Vbis) S = gerepilecopy(av,S);
	return S;
}

GEN detratio(GEN K, GEN T, GEN p, long e, GEN pe)
{ /* K=(K1|K2) -> det(K2)/det(K1) */
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

GEN PicNorm(GEN J, GEN F1, GEN F2, GEN WE, ulong n)
{ /* F,F1 in V_n, WE = V(-E) -> F1/F2(E) */
	pari_sp av = avma;
	ulong g,d0,nS,nZ,nVn2;
	ulong i,j;
	long e;
	GEN V,Vn,T,p,pe;
	GEN WEVn,S,SS,M1,M2,M;

	printf("u");
	g = Jgetg(J);
	d0 = Jgetd0(J);
	V = JgetV(J);
	JgetTpe(J,&T,&pe,&p,&e);
	d0 = Jgetd0(J);
	g = Jgetg(J);
	nVn2 = (n+2)*d0+1-g; /* dim V_{n+2} = L((n+2)*D0) */
	nS = lg(V)-lg(WE); /* codim WE in V = deg E */
	nZ = lg(gel(V,1))-1;
	printf("\nPicNorm %lu\n",n);
	Vn = NULL;
	switch(n)
	{
		case 1:
			Vn = JgetV1(J);
			break;
		case 2:
			Vn = JgetV(J);
			break;
		case 3:
			Vn = JgetV3(J);
			break;
		default:
			pari_err(e_MISC,"In PicNorm, n must be 1,2, or 3");
	}

	M = cgetg(lg(Vn)+1,t_MAT);
	for(j=1;j<lg(Vn);j++) gel(M,j) = gel(Vn,j);
	gel(M,lg(Vn)) = F1;
	printf("Test F1 in V%lu:%ld\n",n,lg(matkerpadic(M,T,p,e))-1);
	gel(M,lg(V)) = F2;
	printf("Test F2 in V%lu:%ld\n",n,lg(matkerpadic(M,T,p,e))-1);

	WEVn = DivAdd(Vn,WE,nVn2-nS,T,p,e,pe,0); /* Vn*WE = V_{n+2}(-E) */

	S = FindSuppl(WE,nS,V,NULL,T,p,pe); /* V = V(-E) + S, dim S = nS */
	SS = FindSuppl(WEVn,nS,V,Vn,T,p,pe); /* V_{n+2} = V_{n+2}(-E) + S, dim SS = nS */

	M = cgetg(nS+nVn2+1,t_MAT);
	for(j=1;j<=nS;j++)
	{
		gel(M,j) = cgetg(nZ+1,t_COL);
		for(i=1;i<=nZ;i++)
			gcoeff(M,i,j) = Fq_mul(gcoeff(S,i,j),gel(F1,i),T,pe);
	}
	for(j=1;j<=nVn2;j++) gel(M,nS+j) = gel(SS,j);
	M1 = detratio(matkerpadic(M,T,p,e),T,p,e,pe);
	if(ZX_is0mod(M1,p))
	{
		if(DEBUGLEVEL) err_printf("PicNorm: F1 has zeros on D, giving up\n"); /* TODO */
		avma = av;
		return NULL;
	}

	for(j=1;j<=nS;j++)
	{
		for(i=1;i<=nZ;i++)
		{
			gcoeff(M,i,j) = Fq_mul(gcoeff(S,i,j),gel(F2,i),T,pe);
		}
	}
	/* TODO should be useless, use above */
	/*M = cgetg(nS+nV5+1,t_MAT);
  for(j=1;j<=nS;j++)
  {
    gel(M,j) = cgetg(nZ+1,t_COL);
    for(i=1;i<=nZ;i++)
      gcoeff(M,i,j) = Fq_mul(gcoeff(V1,i,j),gel(F,i),T,pe);
  }
  for(j=1;j<=nV5;j++) gel(M,nS+j) = gel(V2,j);*/

	/*pari_printf("M=\n%Ps\n",M);*/
	/*pari_printf("detratio of\n%Ps\n",matkerpadic(M,T,p,e));*/
	M2 = detratio(matkerpadic(M,T,p,e),T,p,e,pe);
	if(ZX_is0mod(M2,p))
	{
		if(DEBUGLEVEL) err_printf("PicNorm: F2 has zeros on D, giving up\n");
		avma = av;
		return NULL;
	}
	
	pari_printf("M2=%Ps  ",M2);
	printf("typ2:%ld  ",typ(M2));
	printf("v");
	return gerepileupto(av,ZpXQ_div(M1,M2,T,pe,p,e)); /* TODO can save divisions by taking detratios together */
}

GEN PicFreyRuckMulti1(GEN J, GEN Wtors, GEN l, GEN Wtest, GEN W0, GEN F1, GEN F2, GEN F3, GEN C)
/* Pair the l-tors pt Wtors against the pts in Wtest */
{
	pari_sp av = avma;
	GEN WtorsM,Fq1,H,col,WA,WB,res,s,N;
	GEN T,p,pe,KV;
	long e;
	ulong nC,ntest;
	ulong c,d,i,j;
	ulong n=0;
	GEN Fn=NULL;
	
	printf("a");
	JgetTpe(J,&T,&pe,&p,&e);
	Fq1 = GetFq1(T);
	KV = JgetKV(J,2);
	nC = lg(C);
	ntest = lg(Wtest);
	WtorsM = cgetg(nC,t_VEC);
	gel(WtorsM,1) = Wtors;
	H = zeromatcopy(ntest-1,nC-1);
	gel(H,1) = cgetg(ntest,t_COL);
	for(d=1;d<ntest;d++) gcoeff(H,d,1) = Fq1;
	for(c=2;c<nC;c++)
	{
		printf("b");
		i = gmael(C,c,2)[1];
		j = gmael(C,c,2)[2];
		if(i)
		{
			WA = gel(WtorsM,i);
			if(j)
			{
				WB = gel(WtorsM,j);
				res = PicChord(J,WA,WB,3);
				n = 3;
				Fn = F3;
			}
			else
			{
				res = PicNeg(J,WA,3);
				n = 2;
				Fn = F2;
			}
		}
		else
		{
			if(j==0) pari_err(e_MISC,"Two zeros in Frey-Ruck AddFlip chain, not supposed to happen!");
			WB = gel(WtorsM,j);
			res = PicNeg(J,WB,3);
			n = 2;
			Fn = F2;
		}
		gel(WtorsM,c) = gel(res,1);
		s = gel(res,2); /* in Vn, n=2 or 3 as above */
		col = cgetg(ntest,t_COL);
		for(d=1;d<ntest;d++)
		{
			printf("c");
			N = PicNorm(J,s,Fn,gel(Wtest,d),n);
			printf("0");
			if(N==NULL)
			{
				av = avma;
				return NULL;
			}
			printf("1");
			gel(col,d) = Fq_mul(N,gcoeff(H,d,i),T,pe);
			printf("2");
			if(j) gel(col,d) = Fq_mul(gel(col,d),gcoeff(H,d,j),T,pe);
			printf("3");
			gel(col,d) = ZpXQ_inv(gel(col,d),T,p,e);
			printf("4");
		}
		printf("d");
		N = PicNorm(J,s,Fn,W0,n);
		if(N==NULL)
    {
      av = avma;
      return NULL;
    }
		gel(H,c) = FqC_Fq_mul(col,N,T,pe);
	}
	printf("e");
	/* WtorsM[nC-1] ~ 0, find section that makes this lin equiv happen */
	s = DivSub(DivMul(F1,JgetV1(J),T,pe),gel(WtorsM,nC-1),KV1,1,T,p,e,pe,2); /* TODO get KV1 */
	s = gel(s,1);
	col = gel(H,nC-1);
	for(d=1;d<ntest;d++)
	{
		printf("f");
		N = PicNorm(J,s,Fn,gel(Wtest,d),n);
		if(N==NULL)
    {
      av = avma;
      return NULL;
    }
		gel(col,d) = Fq_mul(gel(col,d),N,T,pe);
	}
	printf("g");
	N = PicNorm(J,s,Fn,W0,n);
	if(N==NULL)
  {
    av = avma;
    return NULL;
  }
	printf("h");
	col = FqC_Fq_mul(col,ZpXQ_inv(N,T,p,e),T,pe);
	return gerepileupto(av,col);
}

GEN PicFreyRuckMulti(GEN J, GEN Wtors, GEN l, GEN Wtest, GEN W0, GEN C)
/* Pair the l-tors pt Wtors against the pts in Wtest */
{
	pari_sp av = avma,av1;
	GEN res=NULL,Wtors1=Wtors,Wtest1,W01=W0,V1,T,pe,F1,F2,F3;
	ulong ntest,i,nZ;
	ntest = lg(Wtest);
	T = JgetT(J);
	pe = Jgetpe(J);
	V1 = JgetV1(J);
	Wtest1 = cgetg(ntest,t_VEC);
	F1 = gel(V1,1);
	nZ = lg(F1);
	F2 = cgetg(nZ,t_COL);
	F3 = cgetg(nZ,t_COL);
	av1 = avma;
	for(i=1;i<ntest;i++)
		gel(Wtest1,i) = gel(Wtest,i);
	for(i=1;i<nZ;i++)
	{
		gel(F2,i) = Fq_mul(gel(F1,i),gel(F1,i),T,pe);
		gel(F3,i) = Fq_mul(gel(F2,i),gel(F1,i),T,pe);
	}

	for(;;)
	{
		res = PicFreyRuckMulti1(J,Wtors1,l,Wtest1,W01,F1,F2,F3,C);
		if(res) return gerepileupto(av,res);
		if(DEBUGLEVEL) err_printf("Error in Frey-Ruck, retrying\n");
		avma = av1;
		Wtors1 = PicNeg(J,Wtors,1);
		W01 = PicNeg(J,W0,1);
		for(i=1;i<ntest;i++)
			gel(Wtest1,i) = PicNeg(J,gel(Wtest,i),1);
		F1 = RandVec_1(V1,pe);
		for(i=1;i<nZ;i++)
  	{
    	gel(F2,i) = Fq_mul(gel(F1,i),gel(F1,i),T,pe);
    	gel(F3,i) = Fq_mul(gel(F2,i),gel(F1,i),T,pe);
  	}
	}
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
