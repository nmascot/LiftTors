#include "zn.h"
#include "../linalg.h"

GEN l2(GEN EN, GEN P, GEN Q, GEN T, GEN pe, GEN p, long e)
{
	P = GetCoef(EN,P);
	Q = GetCoef(EN,Q);
	return ZpXQ_div(ZX_sub(gel(Q,2),gel(P,2)),ZX_sub(gel(Q,1),gel(P,1)),T,pe,p,e);
}

GEN l1(GEN EN, GEN P, GEN Q, GEN T, GEN pe, GEN p, long e)
{
	pari_sp av = avma;
	ulong N,n;
	GEN S,nP;

	/* TODO methode Kamal addchain */
	N = lg(EN)-1;
	S = l2(EN,P,Q,T,pe,p,e);
	nP = P;
	for(n=1;n<N;n++)
	{
		S = ZX_add(S,l2(EN,P,ZC_add(Q,nP),T,pe,p,e));
		nP = ZC_add(nP,P);
	}
	return gerepileupto(av,S);
}

GEN M2_worker(GEN vw, GEN Ml1, GEN TH, GEN Mpts, GEN T, GEN pe)
{
	pari_sp avs;
	ulong nZ, nTH;
	GEN v,w,C,Cs,M,vM,wM;
	ulong s,h;

	nZ = lg(Mpts);
	nTH = lg(TH);
	v = gel(vw,1);
	w = gel(vw,2);

	C = cgetg(nZ,t_COL);
	for(s=1;s<nZ;s++)
	{
		avs = avma;
		Cs = pol_0(varn(T));
		for(h=1;h<nTH;h++)
		{
			M = ZM_mul(gel(TH,h),gel(Mpts,s));
			vM = ZV_ZM_mul(v,M);
			wM = ZV_ZM_mul(w,M);
			Cs = ZX_add(Cs,ZX_mul(GetCoef(Ml1,vM),GetCoef(Ml1,wM)));
		}
		Cs = Fq_red(Cs,T,pe);
		gel(C,s) = gerepileupto(avs,Cs);
	}
	return C;
}

GEN M2mat(GEN M2gens, GEN Ml1, GEN TH, GEN MPts, GEN T, GEN pe)
{
	pari_sp av = avma;
	GEN M2;
	ulong d,j;
	struct pari_mt pt;
	GEN vFixedParams,worker,done;
	long pending,workid;

	d = lg(M2gens);
	vFixedParams = mkvecn(5,Ml1,TH,MPts,T,pe);
	worker = snm_closure(is_entry("M2_worker"),vFixedParams);
	pending = 0;
	M2 = cgetg(d,t_MAT);
	//for(j=1;j<d;j++) gel(M2,j) = M2_worker(gel(M2gens,j),Ml1,TH,MPts,T,pe);
	mt_queue_start_lim(&pt,worker,d-1);
	for(j=1;j<d||pending;j++)
	{
		mt_queue_submit(&pt,j,j<d?mkvec(gel(M2gens,j)):NULL);
		done = mt_queue_get(&pt,&workid,&pending);
		if(done) gel(M2,workid) = done;
	}
	mt_queue_end(&pt);
	return gerepilecopy(av,M2);
}
