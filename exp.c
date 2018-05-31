#include<pari/pari.h>

GEN NAF(GEN n)
{
	pari_sp av = avma;
	GEN B,B2,C,N;
	ulong l,i;

	B = binary_zv(n);
	l = lg(B);
	B2 = cgetg(l+2,t_VECSMALL);
	C = cgetg(l+2,t_VECSMALL);
	N = cgetg(l+1,t_VECSMALL);
	for(i=1;i<l;i++)
	{
		B2[i] = B[l-i];
	}
	C[1] = B2[l] = B2[l+1] = 0;
	for(i=1;i<=l;i++)
	{
		C[i+1] = (C[i]+B2[i]+B2[i+1])/2;
		N[i] = B2[i]+C[i]-2*C[i+1];
	}
	if(N[l]==0) setlg(N,l);
	return gerepileupto(av,N);
}

GEN AddChain(GEN n, long signmatters)
{
	pari_sp av = avma;
	GEN N,A;
	ulong l,i,j;
	long sn,m;

	sn = signe(n);
	setsigne(n,1);
	N = NAF(n);
	l = lg(N);
	A = cgetg(2*l,t_VEC);
	gel(A,1) = mkvecsmall3(1,0,-1);
	j = 1;
	m = 1;
	for(i=l-2;i;i--)
	{
		j++;
		m *= -2;
		gel(A,j) = mkvecsmall3(m,j-1,j-1);
		if(N[i])
		{
			j++;
			if(m*N[i]>0)
			{
				m = -(m+1);
				gel(A,j) = mkvecsmall3(m,j-1,1);
			}
			else
			{
				m = 1-m;
				gel(A,j) = mkvecsmall3(m,j-1,-1);
			}
		}
	}
	if(signmatters && m*sn<0)
	{
		j++;
		gel(A,j) = mkvecsmall3(-m,j-1,0);
	}
	setlg(A,j+1);
	avma = av;
	return gerepileupto(av,A);
}
