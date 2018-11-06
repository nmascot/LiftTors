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
	return gerepilecopy(av,N);
}

GEN AddChain(GEN n, long signmatters)
{
	pari_sp av = avma;
	GEN N,A,m;
	ulong l,i,j,jm1;
	long sn;

  if(equalis(n,3))
	{
		A = mkvecn(4,
				mkvec2(gen_1,mkvecsmall2(0,-1)),
				mkvec2(gen_m2,mkvecsmall2(1,1)),
				mkvec2(gen_m1,mkvecsmall2(1,0)),
				mkvec2(stoi(3),mkvecsmall2(3,2))
			);
		return gerepilecopy(av,A);
	}
  if(equalis(n,-3))
  { 
    A = mkvecn(4,
        mkvec2(gen_1,mkvecsmall2(0,-1)),
        mkvec2(gen_m2,mkvecsmall2(1,1)),
        mkvec2(gen_2,mkvecsmall2(2,0)),
        mkvec2(stoi(-3),mkvecsmall2(3,1))
      );
    return gerepilecopy(av,A);
  }

	sn = signe(n);
	setsigne(n,1);
	N = NAF(n);
	l = lg(N);
	A = cgetg(2*l,t_VEC);
	gel(A,1) = mkvec2(gen_1,mkvecsmall2(0,-1));
	j = 1;
	jm1 = 0;
	m = gen_1;
	for(i=l-2;i;i--)
	{
		j++;
		m = mulii(m,gen_m2);
		gel(A,j) = mkvec2(m,mkvecsmall2(j-1,j-1));
		if(N[i])
		{
			j++;
			if(signe(m)*N[i]>0)
			{
				m = negi(addiu(m,1));
				gel(A,j) = mkvec2(m,mkvecsmall2(j-1,1));
			}
			else
			{
				m = subui(1,m);
				if(jm1==0)
				{
					jm1 = j;
					gel(A,jm1) = mkvec2(gen_m1,mkvecsmall2(1,0));
					j++;
					gel(A,j) = mkvec2(m,mkvecsmall2(j-2,jm1));
				}
				else gel(A,j) = mkvec2(m,mkvecsmall2(j-1,jm1));
			}
		}
	}
	if(signmatters && signe(m)*sn<0)
	{
		j++;
		gel(A,j) = mkvec2(negi(m),mkvecsmall2(j-1,0));
	}
	setlg(A,j+1);
	avma = av;
	return gerepilecopy(av,A);
}
