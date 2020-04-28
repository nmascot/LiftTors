#include<pari/pari.h>

ulong ZNnorm(long x, ulong N)
{
	x = umodsu(x,N);
	return x?x:N;
}

GEN GetCoef(GEN A, GEN v)
{
	long N;
	long i,j;

	N = lg(A)-1;
	//printf("N=%ld\n",N);
	i = smodis(gel(v,1),N);
	j = smodis(gel(v,2),N);
	return gcoeff(A,i?i:N,j?j:N);
}

ulong ZNneg(long x, ulong N)
{
	ulong y;
	y = umodsu(x,N);
	return y?N-y:N;
}
