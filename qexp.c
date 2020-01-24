#include "zn.c"

GEN E1qexp(GEN v, ulong N, GEN zpows, ulong B, GEN T, GEN pe, GEN p, long e)
{ /* v=[c,d] mod N, zpows = powers of primitive Nth root of 1: q-exp of E_1^[c,d] up to O(qN^B) */
	/* TODO use t_SER ? */
	pari_sp av = avma;
	GEN E,Fq0,a0,zd;
	ulong a,b,c,d,m,n,t;
	
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
	for(t=(c==0?1:0);(a=N*t+c)<B;t++) /* Case a==c mod N */
	{
		for(b=1;(n=a*b+1)<=B;b++)
		{
			m = ZNnorm(b*d,N);
			gel(E,n) = ZX_add(gel(E,n),gel(zpows,m));
		}
	}
	for(t=1;(a=N*t-c)<B;t++) /* Case a==-c mod N */
  {
    for(b=1;(n=a*b+1)<=B;b++)
    {
      m = ZNneg(b*d,N);
      gel(E,n) = ZX_sub(gel(E,n),gel(zpows,m));
    }
  }

	return gerepilecopy(av,E);
}
