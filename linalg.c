#include<pari/pari.h>

long FpX_is0modp(GEN x,GEN p)
{
	pari_sp av = avma;
	GEN red;
	long res;
	red = (typ(x)==t_POL?FpX_red(x,p):Fp_red(x,p));
	res = gequal0(red);
	avma = av;
	return res;
}

GEN FpXM_red(GEN A, GEN p)
{
	long m,n,i,j;
	GEN B,c;
	RgM_dimensions(A,&m,&n);
	B = cgetg(n+1,t_MAT);
	for(j=1;j<=n;j++)
	{
		gel(B,j) = cgetg(m+1,t_COL);
		for(i=1;i<=m;i++)
		{
			c = gcoeff(A,i,j);
			gcoeff(B,i,j) = (typ(c)==t_POL?FpX_red(c,p):Fp_red(c,p));
		}
	}
	return B;
}

GEN FpXM_add(GEN A, GEN B, GEN p)
{
	long m,n,i,j;
	GEN C;
	RgM_dimensions(A,&m,&n);
	C = cgetg(n+1,t_MAT);
	for(j=1;j<=n;j++)
  {
    gel(C,j) = cgetg(m+1,t_COL);
    for(i=1;i<=m;i++)
    {
      gcoeff(C,i,j) = FpX_add(gcoeff(A,i,j),gcoeff(B,i,j),p);
    }
  }
  return C;
}

GEN FpXM_sub(GEN A, GEN B, GEN p)
{
  long m,n,i,j;
  GEN C;
  RgM_dimensions(A,&m,&n);
	C = cgetg(n+1,t_MAT);
  for(j=1;j<=n;j++)
  {
    gel(C,j) = cgetg(m+1,t_COL);
    for(i=1;i<=m;i++)
    {
      gcoeff(C,i,j) = FpX_sub(gcoeff(A,i,j),gcoeff(B,i,j),p);
    }
  }
  return C;
}

GEN FqV_Fq_mul(GEN v, GEN a, GEN T, GEN p)
{
	ulong n,i;
	GEN av;
	n = lg(v);
	av = cgetg(n,t_VEC);
	for(i=1;i<n;i++)
	{
		gel(av,i) = Fq_mul(a,gel(v,i),T,p);
	}
	return av;
}

GEN FqM_Fq_mul(GEN v, GEN a, GEN T, GEN p)
{
  ulong n,i;
  GEN av;
  n = lg(v);
  av = cgetg(n,t_MAT);
  for(i=1;i<n;i++)
  {
    gel(av,i) = FqC_Fq_mul(gel(v,i),a,T,p);
  }
  return av;
}

GEN ZXM_Z_mul(GEN A, GEN a)
{
	long m,n,i,j;
	GEN B,col;
	RgM_dimensions(A,&m,&n);
	B = cgetg(n+1,t_MAT);
	for(j=1;j<=n;j++)
	{
		col = cgetg(m+1,t_COL);
		for(i=1;i<=m;i++)
		{
			gel(col,i) = ZX_Z_mul(gcoeff(A,i,j),a);
		}
		gel(B,i) = col;
	}
	return B;
}

GEN RandVec_padic(GEN A, GEN T, GEN p, GEN pe)
{
	pari_sp av = avma;
	unsigned long m,n,i,j;
	long dT,vT;
	GEN v,b,c;

	dT = degree(T);
	vT = varn(T);
	n = lg(A);
	m = lg(gel(A,1));
	dT = lg(T);
	v = cgetg(m,t_COL);
	for(j=1;j<n;j++)
	{
		b = random_FpX(dT-1,vT,p);
		for(i=1;i<m;i++)
		{
			c = Fq_mul(b,gcoeff(A,i,j),T,pe);
			if(j==1) gel(v,i) = c;
			else gel(v,i) = Fq_add(gel(v,i),c,T,pe);
		}
		v = gerepilecopy(av,v);
	}
	return v;
}

GEN Hsort(GEN A, GEN p)
{
	ulong off=0,i=1,j=1,n,m,all0=1;
	n = lg(A);
	for(j=1;j<n;j++)
	{
		all0 = 1;							
		m = lg(gel(A,j));
		for(i=1;i<m;i++)
		{
			if(!FpX_is0modp(gcoeff(A,i,j),p))
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

GEN matkerpadic(GEN A, GEN T, GEN p, long e)
{
	pari_sp av;
	GEN K;
	if(e==1) return FqM_ker(A,T,p);
	av = avma;
	K = ZpXQM_ker(A,T,p,e,NULL);
	K = Hsort(K,p);
	return gerepilecopy(av,K);
}

GEN mateqnpadic(GEN A, GEN T, GEN p, long e)
{
	pari_sp av = avma;
	return gerepilecopy(av,shallowtrans(matkerpadic(shallowtrans(A),T,p,e)));
}

GEN matimagepadic(GEN A, GEN T, GEN p, long e)
{
  pari_sp av;
  GEN K;
  if(e==1) return FqM_image(A,T,p);
  av = avma;
  K = ZpXQM_image(A,T,p,e,NULL);
  return K;
  K = Hsort(K,p);
  return gerepilecopy(av,K);
}

GEN matF(GEN A, GEN T, GEN p, long e)
{
	pari_sp av = avma;
	GEN B;
	ulong r,n,j;
	r = lg(A);
	/* reduce A mod p, and supplement */
	B = FpXM_red(A,p);
	B = FqM_suppl(B,T,p);
	/* TODO necessary? Un-reduce the A-part */
	for(j=1;j<r;j++)
	{
		gel(B,j) = gel(A,j);
	}
	/* invert */
	B = ZpXQM_inv(B,T,p,e);
	/* return the first #A rows */
	n = lg(B);
	for(j=1;j<n;j++)
	{
		setlg(gel(B,j),r);
	}
	return gerepilecopy(av,B);
}

GEN mat2col(GEN A)
{
	unsigned long n,m,i,j=1;
	GEN C;
	n = lg(A)-1;
	if(n==0) return cgetg(1,t_COL);
	m = lg(gel(A,1))-1;
	C = cgetg(n*m+1,t_COL);
	for(i=1;i<=m;i++)
	{
		for(j=1;j<=n;j++)
		{
			gel(C,(i-1)*n+j) = gcoeff(A,i,j);
		}
	}
	return C;
}

GEN col2mat(GEN C, unsigned long m, unsigned long n)
{
	GEN A,Aj;
	unsigned long i=1,j=1;
	A = cgetg(n+1,t_MAT);
	for(j=1;j<=n;j++)
	{
		Aj = cgetg(m+1,t_COL);
		for(i=1;i<=m;i++)
		{
			gel(Aj,i) = gel(C,(i-1)*n+j);
		}
		gel(A,j) = Aj;
	}
	return A;
}

GEN M2ABCD(GEN M, GEN uv)
{
	GEN u,v,res,A,col;
	unsigned long m,n,i,j,p,q;
	res = cgetg(5,t_VEC);
	for(p=1;p<=2;p++)
	{
		u = gmael(uv,1,p);
		m = lg(u);
		for(q=1;q<=2;q++)
		{
			v = gmael(uv,2,q);
			n = lg(v);
			A = cgetg(n,t_MAT);
			for(j=1;j<n;j++)
			{
				col = cgetg(m,t_COL);
				for(i=1;i<m;i++)
				{
					gel(col,i) = gcoeff(M,u[i],v[j]);
				}
				gel(A,j) = col;
			}
			gel(res,q+2*(p-1)) = A;
		}
	}
	return res;
}

GEN M2ABCD_1block(GEN M, ulong top, ulong left, GEN uv)
/* Same as above, but all zeros except for block M starting at top+1,left+1 */
/* /!\ Not suitable for gerepile */
{
  GEN u,v,res,A,col;
	long m,n;
  ulong bot,right,i,j,p,q,ui,vj;
  res = cgetg(5,t_VEC);
	RgM_dimensions(M,&m,&n);
	bot = top+m;
	right = left+n;
  for(p=1;p<=2;p++)
  {
    u = gmael(uv,1,p);
    m = lg(u);
    for(q=1;q<=2;q++)
    {
      v = gmael(uv,2,q);
      n = lg(v);
      A = cgetg(n,t_MAT);
      for(j=1;j<n;j++)
      {
				col = cgetg(m,t_COL);
				vj = v[j];
        for(i=1;i<m;i++)
        {
					ui = u[i];
					if(vj>left && vj<=right && ui>top && ui<=bot) gel(col,i) = gcoeff(M,ui-top,vj-left);
					else gel(col,i) = gen_0;
        }
				gel(A,j) = col;
      }
      gel(res,q+2*(p-1)) = A;
    }
  }
  return res;
}

GEN VecSmallCompl(GEN v, ulong n)
{
	ulong iv,ic,m;
	GEN c;
	c = cgetg(n+2-lg(v),t_VECSMALL);
	iv = ic = 1;
	for(m=1;m<=n;m++)
	{
		if(m<v[iv])	c[ic++]=m;
		else iv++;
	}
	return c;
}

GEN FqM_MinorCompl(GEN A, GEN T, GEN p)
{
	pari_sp av=avma;
	GEN IJ,uv;
	long m,n;
	RgM_dimensions(A,&m,&n);
	IJ = FqM_indexrank(A,T,p);
	uv = cgetg(3,t_VEC);
	gel(uv,1) = cgetg(3,t_VEC);
	gel(uv,2) = cgetg(3,t_VEC);
	gmael(uv,1,1) = gel(IJ,1);
	gmael(uv,1,2) = VecSmallCompl(gel(IJ,1),m);
	gmael(uv,2,1) = gel(IJ,2);
	gmael(uv,2,2) = VecSmallCompl(gel(IJ,2),n);
	return gerepilecopy(av,uv);
}
