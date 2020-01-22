#include<pari/pari.h>

GEN GetCoef(GEN A, GEN v)
{
	ulong N;
	long i,j;

	N = lg(A)-1;
	//printf("N=%ld\n",N);
	/* TODO improve this horror */
	i = itos(gel(v,1));
	j = itos(gel(v,2));
	//printf("i=%ld, j=%ld\n",i,j);
	if(i<0) i = N-((-i)%N);
	else { i=i%N; if(i==0) i=N; }
	if(j<0) j = N-((-j)%N);
  else { j=j%N; if(j==0) j=N; }
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
