#include<pari/pari.h>

ulong ZNnorm(long x, ulong N)
{
  if(x>0) x -= ((x-1)/N)*N;
  else x += N*(1+(-x)/N);
	return x;
}

GEN GetCoef(GEN A, GEN v)
{
	ulong N,i,j;

	N = lg(A)-1;
	//printf("N=%ld\n",N);
	i = itos(gel(v,1));
	if(i>0) i -= ((i-1)/N)*N;
  else i += N*(1+(-i)/N);
	j = itos(gel(v,2));
	if(j>0) j -= ((j-1)/N)*N;
  else j += N*(1+(-j)/N);
	return gcoeff(A,i,j);
}

ulong ZNneg(long x, ulong N)
{
	if(x<0)
	{
		x = -x;
		x -= N*((x-1)/N);
	}
	else x = N*(1+(x/N))-x;
	return x;
}
