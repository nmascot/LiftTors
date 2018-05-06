#include<pari/pari.h>

GEN RandVec(GEN A, GEN t)
{
	pari_sp av = avma;
	long n = lg(A);
	GEN v;
	unsigned long i;
	v = cgetg(n,t_COL);
	for(i=1;i<n;i++)
	{
		gel(v,i) = genrand(t);
	}
	return gerepileupto(av,gmul(A,v));
}

GEN Hsort(GEN A, GEN p)
{
	pari_sp av;
	GEN red;
	unsigned long off=0,i=1,j=1,n,m,all0=1;
	n = lg(A);
	av = avma;
	for(j=1;j<n;j++)
	{
		all0 = 1;							
  	m = lg(gel(A,j));
  	for(i=1;i<m;i++)
		{
			avma = av;
			red = gcoeff(A,i,j);
			if(typ(red) == t_INT)
				red = Fp_red(red,p);
			else
			 red = FpX_red(red,p);
			if(!gequal0(red))
			{
				all0 = 0;
				break;
			}
		}
		if(all0 == 1)
    {
			off++;
		}
		else
		{
			gel(A,j-off) = gel(A,j);
		}
	}
	setlg(A,n-off);
	return A;
}

GEN matkerMod(GEN A, GEN p, GEN T, GEN e)
{
	pari_sp av;
	GEN K;
	if(gequal1(e))
		return FqM_ker(A,T,p);
	av = avma;
	K = ZpXQM_ker(A,T,p,itos(e),NULL);
	return K;
	K = Hsort(K,p);
	return gerepilecopy(av,K);
}
