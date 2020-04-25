#include<linalg.h>

GEN GetFq1(GEN T)
{
	GEN Fq1;
	Fq1 = mkpoln(1,gen_1);
	Fq1[1] = 0;
	setsigne(Fq1,1);
  setvarn(Fq1,varn(T));
	return Fq1;
}

GEN GetFq0(GEN T)
{
  GEN Fq0;
  Fq0 = mkpoln(0);
	Fq0[1] = 0;
	setsigne(Fq0,0);
  setvarn(Fq0,varn(T));
  return Fq0;
}

GEN Z2Fq(GEN x, GEN T)
{
	if(gequal0(x)) return GetFq0(T);
  GEN y = mkpoln(1,x);
	y[1] = 0;
  setsigne(y,1);
  setvarn(y,varn(T));
  return y;
}

long ZX_is0mod(GEN x, GEN p)
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
      gcoeff(C,i,j) = FpX_add(gcoeff(A,i,j),gcoeff(B,i,j),p);
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
      gcoeff(C,i,j) = FpX_sub(gcoeff(A,i,j),gcoeff(B,i,j),p);
  }
  return C;
}

GEN ZXM_add(GEN A, GEN B)
{
  long m,n,i,j;
  GEN C;
  RgM_dimensions(A,&m,&n);
  C = cgetg(n+1,t_MAT);
  for(j=1;j<=n;j++)
  {
    gel(C,j) = cgetg(m+1,t_COL);
    for(i=1;i<=m;i++)
      gcoeff(C,i,j) = ZX_add(gcoeff(A,i,j),gcoeff(B,i,j));
  }
  return C;
}

GEN ZXM_sub(GEN A, GEN B)
{
  long m,n,i,j;
  GEN C;
  RgM_dimensions(A,&m,&n);
  C = cgetg(n+1,t_MAT);
  for(j=1;j<=n;j++)
  {
    gel(C,j) = cgetg(m+1,t_COL);
    for(i=1;i<=m;i++)
      gcoeff(C,i,j) = ZX_sub(gcoeff(A,i,j),gcoeff(B,i,j));
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
		gel(B,j) = col;
	}
	return B;
}

GEN FpXC_add(GEN A, GEN B, GEN p)
{
	ulong n = lg(A),i;
	GEN C;
	C = cgetg(n,t_COL);
	for(i=1;i<n;i++) gel(C,i) = FpX_add(gel(A,i),gel(B,i),p);
	return C;
}

GEN FpXC_sub(GEN A, GEN B, GEN p)
{
	ulong n = lg(A),i;
	GEN C;
	C = cgetg(n,t_COL);
	for(i=1;i<n;i++) gel(C,i) = FpX_sub(gel(A,i),gel(B,i),p);
	return C;
}

GEN ZXC_add(GEN A, GEN B)
{
  ulong n = lg(A),i;
  GEN C;
  C = cgetg(n,t_COL);
  for(i=1;i<n;i++) gel(C,i) = ZX_add(gel(A,i),gel(B,i));
  return C;
}

GEN ZXC_sub(GEN A, GEN B)
{
  ulong n = lg(A),i;
  GEN C;
  C = cgetg(n,t_COL);
  for(i=1;i<n;i++) gel(C,i) = ZX_sub(gel(A,i),gel(B,i));
  return C;
}


GEN RandVec_1(GEN A, GEN pe)
{
  pari_sp av = avma;
  ulong n,j;
  GEN v;
  n = lg(A);
  v = NULL;
	do{
  	for(j=1;j<n;j++)
  	{
			if(random_Fl(2))
			{
				if(v==NULL)
				{
					v = gcopy(gel(A,j));
				}
				else
				{
					if(random_Fl(2)) v = ZXC_sub(v,gel(A,j));
					else v = ZXC_add(v,gel(A,j));
				}
    	}
  	}
	} while(v==NULL);
  return gerepileupto(av,v);
}

GEN RandVec_padic(GEN A, GEN T, GEN p, GEN pe)
{
	pari_sp av = avma;
	unsigned long m,n,i,j;
	long dT,vT;
	GEN v,b,c;

	dT = lg(T);
	vT = varn(T);
	n = lg(A);
	m = lg(gel(A,1));
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
			if(!ZX_is0mod(gcoeff(A,i,j),p))
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

GEN ZpXQMinv(GEN A, GEN T, GEN pe, GEN p, long e)
{
	pari_sp av = avma, avk;
	ulong n,i,j,k;
	GEN B,col,l,col2,C;
	GEN I;
	GEN Fq0,Fq1;

	n = lg(A)-1;
	I = cgetg(n,t_VECSMALL); /* Vector of permutation of rows */
	for(i=1;i<=n;i++) I[i] = i;
	Fq0 = GetFq0(T);
	Fq1 = GetFq1(T);
	B = cgetg(n+1,t_MAT);
	/* Phase 1: to U1 form */
	for(k=n;k;k--)
	{
		avk = avma;
		col = cgetg(2*n+1,t_COL);
		/* Cols n to k+1 processed, now do col k */
		for(i=1;i<=n;i++)
		{
			gel(col,i) = gcoeff(A,i,k);
			gel(col,i+n) = i==k? Fq1 : Fq0;
		}
		/* Col is now vcat Ak,Ik */
		/* The last k cols of B are known */
		/* Upper part of B is UT1, lower part is LT */
		/* Reduce Ck uing Cj for j=n..k+1 */
		/* Actually only need coeffs for i<k and i>=n+k */
		for(j=n;j>k;j--)
		{
			l = Fq_red(ZX_neg(gel(col,I[j])),T,pe);
			/* Ck += l*Cj */
			for(i=1;i<=j;i++)
				gel(col,I[i]) = ZX_add(gel(col,I[i]),ZX_mul(l,gcoeff(B,I[i],j)));
			/* Useless
			for(i=j;i<=n;i++)
				gel(col,I[i]) = Fq0; */
			for(i=n+j;i<=2*n;i++)
				gel(col,i) = ZX_add(gel(col,i),ZX_mul(l,gcoeff(B,i,j)));
		/* TODO gerepile possible ici */
			col = gerepileupto(avk,col);
		}
		for(i=1;i<=k;i++) gel(col,I[i]) = Fq_red(gel(col,I[i]),T,pe);
		for(i=n+k+1;i<=2*n;i++) gel(col,i) = Fq_red(gel(col,i),T,pe);
		/* Now coefs k+1..n of col are 0 */
		/* Find pivot above k (hopefully k) */
		for(i=k;i;i--)
		{
			l = gel(col,I[i]);
			if(ZX_is0mod(l,p)==0) break;
		}
		if(i!=k)
		{
			j = I[k]; I[k] = I[i]; I[i] = j;
		}
		/* Divide by pivot */
		l = ZpXQ_inv(l,T,p,e);
		col2 = cgetg(2*n+1,t_COL);
		for(i=1;i<k;i++) gel(col2,I[i]) = Fq_mul(gel(col,I[i]),l,T,pe);
		gel(col2,I[k]) = Fq1;
		for(i=k+1;i<=n;i++) gel(col2,I[i]) = Fq0;
		for(i=n+1;i<n+k;i++) gel(col2,i) = Fq0;
		for(i=n+k;i<=2*n;i++) gel(col2,i) = Fq_mul(gel(col,i),l,T,pe);
		gel(B,k) = gerepileupto(avk,col2);
	}
	/* Phase 2 : to 1 form */
	/* Upper B is UT1, imagine we reduce it to 1 but do not actually do it */
	/* Then the inverse is lower B with cols permuted by I */
	C = cgetg(n+1,t_MAT);
	gel(C,I[1]) = cgetg(n+1,t_COL);
	for(i=1;i<=n;i++)
		gcoeff(C,i,I[1]) = gcoeff(B,i+n,1);
	for(k=2;k<=n;k++)
	{
		avk = avma;
		col = cgetg(n+1,t_COL);
		for(i=1;i<=n;i++)
			gel(col,i) = gcoeff(B,i+n,k);
		for(j=1;j<k;j++)
		{
			l = FpX_neg(gcoeff(B,I[j],k),pe);
			/* Ck += l*Cj */
			for(i=1;i<=n;i++) gel(col,i) = ZX_add(gel(col,i),ZX_mul(l,gcoeff(C,i,I[j])));
			/* TODO gerepile possible ici */
			col = gerepileupto(avk,col);
		}
		col2 = cgetg(n+1,t_COL);
		for(i=1;i<=n;i++) gel(col2,i) = Fq_red(gel(col,i),T,pe);
		gel(C,I[k]) = gerepileupto(avk,col2);
	}
	/* Ensure C is suitable for gerepile */
	for(i=1;i<=n;i++)
    gcoeff(C,i,I[1]) = gcopy(gcoeff(C,i,I[1]));
	return gerepileupto(av,C);
}

GEN matkerpadic_safe(GEN A, GEN T, GEN p, long e)
{
  pari_sp av = avma;
	GEN K;
  if(e==1) return FqM_ker(A,T,p);
  K = ZpXQM_ker(A,T,p,e,NULL);
  K = Hsort(K,p);
  return gerepilecopy(av,K);
}

GEN matkerpadic(GEN A, GEN T, GEN pe, GEN p, long e)
{ /* Assumes good red, i.e. the rank does not decrease mod p */
	pari_sp av = avma, av1;
	GEN IJ,I,J,J1,P,A1,A2,B,K;
	ulong n,r,i,j;
	GEN Fq0,Fqm1;
	if(e==1) return FqM_ker(A,T,p);
	n = lg(A)-1;
	IJ = FqM_indexrank(A,T,p);
	I = gel(IJ,1); /* Rows spanning the eqns, ignore others (good red) */
	J = gel(IJ,2); /* Cols forming invertible block */
	r = lg(J)-1;
	J1 = VecSmallCompl(J,n); /* Other cols */
	P = cgetg(n+1,t_VECSMALL);
	for(j=1;j<=r;j++) P[j] = J[j];
	for(j=1;j<=n-r;j++) P[r+j] = J1[j];
	av1 = avma;
	A1 = cgetg(r+1,t_MAT); /* Invertible block */
	for(j=1;j<=r;j++)
	{
		gel(A1,j) = cgetg(r+1,t_COL);
		for(i=1;i<=r;i++) gcoeff(A1,i,j) = gcoeff(A,I[i],P[j]);
	}
	B = gerepileupto(av1,ZpXQMinv(A1,T,pe,p,e));
	A2 = cgetg(n-r+1,t_MAT); /* Other block */
	for(j=1;j<=n-r;j++)
  {
    gel(A2,j) = cgetg(r+1,t_COL);
    for(i=1;i<=r;i++) gcoeff(A2,i,j) = gcoeff(A,I[i],P[j+r]);
  }
	/* K = vcat of A1^-1*A2, -Id_{n-r}, with perm P^-1 of rows */
	B = gerepileupto(av1,FqM_mul(B,A2,T,pe));
	Fq0 = GetFq0(T);
	Fqm1 = Z2Fq(gen_m1,T);
	K = cgetg(n-r+1,t_MAT);
	for(j=1;j<=n-r;j++)
	{
		gel(K,j) = cgetg(n+1,t_COL);
		for(i=1;i<=r;i++) gcoeff(K,P[i],j) = gcoeff(B,i,j);
		for(i=r+1;i<=n;i++)
			gcoeff(K,P[i],j) = j+r==i?Fqm1:Fq0;
	}
	return gerepilecopy(av,K);
}

GEN mateqnpadic(GEN A, GEN T, GEN pe, GEN p, long e)
{ /* Assumes good red, i.e. the rank does not decrease mod p */
	pari_sp av = avma, av1;
	GEN IJ,I,J,I1,P,A1,A2,B,E;
	ulong n,r,i,j;
	GEN Fq0,Fqm1;
	if(e==1)
		return gerepilecopy(av,shallowtrans(FqM_ker(shallowtrans(A),T,p)));
	n = lg(gel(A,1))-1;
	IJ = FqM_indexrank(A,T,p);
  I = gel(IJ,1); /* Rows forming invertible block */
  J = gel(IJ,2); /* Columns spanning the space, ignore others (good red) */
  r = lg(I)-1;
  I1 = VecSmallCompl(I,n); /* Other rows */
  P = cgetg(n+1,t_VECSMALL); /* First I then I1 */
  for(i=1;i<=r;i++) P[i] = I[i];
  for(i=1;i<=n-r;i++) P[r+i] = I1[i];
	av1 = avma;
  A1 = cgetg(r+1,t_MAT); /* Invertible block */
	for(j=1;j<=r;j++)
  {
    gel(A1,j) = cgetg(r+1,t_COL);
    for(i=1;i<=r;i++) gcoeff(A1,i,j) = gcoeff(A,P[i],J[j]);
	}
	B = gerepileupto(av1,ZpXQM_inv(A1,T,p,e));
	A2 = cgetg(r+1,t_MAT); /* Other block */
  for(j=1;j<=r;j++)
  {
    gel(A2,j) = cgetg(n-r+1,t_COL);
    for(i=1;i<=n-r;i++) gcoeff(A2,i,j) = gcoeff(A,P[i+r],J[j]);
  }
	/* E = hcat of A2*A1^-1, -Id_{n-r}, with perm P^-1 of cols*/
  B = gerepileupto(av1,FqM_mul(A2,B,T,pe));
  Fq0 = GetFq0(T);
  Fqm1 = Z2Fq(gen_m1,T);
  E = cgetg(n+1,t_MAT);
  for(j=1;j<=r;j++)
  {
    gel(E,P[j]) = cgetg(n-r+1,t_COL);
		for(i=1;i<=n-r;i++)
			gcoeff(E,i,P[j]) = gcoeff(B,i,j);
	}
	for(j=r+1;j<=n;j++)
	{
		gel(E,P[j]) = cgetg(n-r+1,t_COL);
    for(i=1;i<=n-r;i++) gcoeff(E,i,P[j]) = i+r==j?Fqm1:Fq0;
  }
  return gerepilecopy(av,E);
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

GEN RgM_drop_rows(GEN A, GEN I)
{
	GEN B;
	ulong i,j,k,iB,iI;
	long m,n;
	RgM_dimensions(A,&m,&n);
	k = lg(I)-1;
	B = cgetg(n,t_MAT);
	for(j=1;j<=n;j++)
	{
		gel(B,j) = cgetg(m-k+1,t_COL);
		iI=iB=1;
		for(i=1;i<=m;i++)
		{
			if(I[iI]==i) iI++;
			else gcoeff(B,iB++,j)=gcoeff(A,i,j);
		}
	}
	return B;
}

GEN Subspace_normalize(GEN V, GEN I, GEN T, GEN pe, GEN p, long e, long drop)
{ /* V represents a subspace, I list of rows. Change basis so that I-block of V==1. */
	pari_sp av = avma;
	GEN P;
	ulong n,i,j;
	n = lg(V);
	P = cgetg(n,t_MAT);
	for(j=1;j<n;j++)
	{
		gel(P,j) = cgetg(n,t_COL);
		for(i=1;i<n;i++)
			gcoeff(P,i,j) = gcoeff(V,I[i],j);
	}
	P = ZpXQM_inv(P,T,p,e);
	V = FqM_mul(V,P,T,pe);
	if(drop)
	{
		V = RgM_drop_rows(V,I);
		return gerepilecopy(av,V);
	}
	else return gerepileupto(av,V);
}
