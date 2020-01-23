#include<pari/pari.h>

GEN GetCoef(GEN A, GEN v)
{
	ulong N;
	long i,j;

	N = lg(A)-1;
	//printf("N=%ld\n",N);
	i = itos(gel(v,1));
	j = itos(gel(v,2));
	//printf("i=%ld, j=%ld\n",i,j);
	if(i>0) i -= ((i-1)/N)*N;
	else i += N*(1+(-i)/N);
	if(j>0) j -= ((j-1)/N)*N;
	else j += N*(1+(-j)/N);
	//printf("i=%ld, j=%ld\n",i,j);
	return gcoeff(A,i,j);
}

long ZNnorm(long x, ulong N)
{
	long y;
	y = x % N;
	if(y==0) { y = N; }
	return y;
}

long ZNneg(long x, ulong N)
{
	long y;
	y = N - (x%N);
	return y;
}
