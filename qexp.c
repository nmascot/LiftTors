#include "zn.c"

GEN E1qexp(GEN v, ulong N, GEN zpows, ulong B, GEN T, GEN pe, GEN p, long e)
{ /* v=[c,d] mod N, zpows = powers of primitive Nth root of 1: q-exp of E_1^[c,d] up to O(qN^B) */
	/* TODO use t_SER ? */
	pari_sp av = avma;
	GEN E,Fq0,a0,zd;
	ulong a,b,c,d,m,n;

	if(B==0) return cgetg(t_VEC,1);

	Fq0 = mkpoln(0);
	setvarn(Fq0,varn(T));

	c = umodiu(gel(v,1),N);
	d = umodiu(gel(v,2),N);

	E = cgetg(B+1,t_VEC);
	for(m=1;m<=B;m++) gel(E,m) = Fq0;

	if(c)
	{ /* a0 = 1/2 - c/N */
		a0 = subii(Fp_inv(gen_2,pe),Fp_div(utoi(c),utoi(N),pe));
		a0 = mkpoln(1,a0);
		setvarn(a0,varn(T));
	}
	else
	{ /* a0 = 1/2 * (1+z^d)/(1-z^d) */
		m = ZNnorm(d,N);
		zd = gel(zpows,m);
		a0 = ZpXQ_div(ZX_Z_add(zd,gen_1),ZX_Z_mul(Z_ZX_sub(gen_1,zd),gen_2),T,pe,p,e);
	}
	gel(E,1) = a0;

	/* sum_{a>0,b>0} if(a==c mod N, z^(b*d) * q^(a*b)) - if(a==-c mod N, z^(-b*d) * q^(a*b)) */
	for(a=(c==0?N:c);a<B;a+=N) /* Case a==c mod N */
	{
		for(b=1;(n=a*b+1)<=B;b++)
		{
			m = ZNnorm(b*d,N);
			gel(E,n) = ZX_add(gel(E,n),gel(zpows,m));
		}
	}
	for(a=N-c;a<B;a+=N) /* Case a==-c mod N */
  {
    for(b=1;(n=a*b+1)<=B;b++)
    {
      m = ZNneg(b*d,N);
      gel(E,n) = ZX_sub(gel(E,n),gel(zpows,m));
    }
  }

	return gerepilecopy(av,E);
}

GEN TrE2qexp(GEN vw, ulong N, GEN H, GEN M, ulong w, GEN zpows, ulong B, GEN T, GEN pe, GEN p, long e)
{ /* vw=[v,w] -> qexp of Tr_H(E_1^v * E_1^w) | M in terms of qw up to O(qw^B) */
	pari_sp av = avma;
	ulong Nw,Nwi,BN,nH,h,i,j;
	GEN Fq0,E,hM,fv,fw;

	if(B==0) return cgetg(1,t_VEC);
	Nw = N/w;
	BN = (B-1)*N/w+1;

	Fq0 = mkpoln(0);
  setvarn(Fq0,varn(T));

	E = cgetg(B+1,t_VEC);
	for(i=1;i<=B;i++) gel(E,i) = Fq0;

	nH = lg(H);
	for(h=1;h<nH;h++)
	{
		hM = ZM_mul(gel(H,h),M);
		fv = E1qexp(ZV_ZM_mul(gel(vw,1),hM),N,zpows,BN,T,pe,p,e);
		fw = E1qexp(ZV_ZM_mul(gel(vw,2),hM),N,zpows,BN,T,pe,p,e);
		/* TODO use fast series multiplication */
		/* f[i] = sum_j fv[j]*fw[Nw*i-j] */
		for(i=0;i<B;i++)
		{
			Nwi = Nw*i;
			for(j=0;j<=Nwi;j++)
				gel(E,i+1) = ZX_add(gel(E,i+1),Fq_mul(gel(fv,j+1),gel(fw,Nwi+1-j),T,pe));
		}
	}

	return gerepilecopy(av,E);
}
