#include "linalg.h"
#include "pic.h"

ulong c2i(GEN c, ulong l)
{
	ulong i,d,n;
	d = lg(c)-1;
	i = 0;
	for(n=1;n<=d;n++)
	{
		i *= l;
		i += (c[n]%l);
	}
	return i?i:upowuu(l,d);
}

GEN i2c(ulong i, ulong l, ulong d)
{
	GEN c;
	ulong n;
	c = cgetg(d+1,t_VECSMALL);
	for(n=1;n<=d;n++)
	{
		c[d+1-n] = i%l;
		i /= l;
	}
	return c;
}

ulong Chordi(ulong i1, ulong i2, ulong l, ulong d)
{
	pari_sp av = avma;
	GEN c1,c2,c3;
	ulong n;
	c1 = i2c(i1,l,d);
	c2 = i2c(i2,l,d);
	c3 = cgetg(d+1,t_VECSMALL);
	for(n=1;n<=d;n++)
		c3[n] = (2*l-(c1[n]+c2[n]))%l;
	n = c2i(c3,l);
	avma = av;
	return n;
}

ulong ActOni(GEN m, ulong i, ulong l)
{
	pari_sp av = avma;
	GEN c;
	ulong d;
	d = lg(m)-1;
	c = i2c(i,l,d);
	c = Flm_Flc_mul(m,c,l);
	i = c2i(c,l);
	avma = av;
	return i;
}

GEN TorsSpaceFrob_worker(GEN J, GEN W1, GEN X1, GEN W2, GEN X2)
{
	pari_sp av = avma;
	GEN W;
	ulong x,x1,x2,y;
	printf("Into Tors_worker\n");
	x1 = itou(X1);
	if(W2==gen_0)
	{
		W = PicNeg(J,W1,0);
		for(y=1;y<=x1;y++) W = gerepileupto(av,PicFrob(J,W));
		printf("Exiting Tors_worker (0)\n");
		return W;
	}
	x2 = itou(X2);
	if(x1>=x2)
	{
		x=x2;
		for(y=1;y<=x1-x2;y++) W1 = gerepileupto(av,PicFrob(J,W1));
	}
	else
	{
    x=x1;
    for(y=1;y<=x2-x1;y++) W2 = gerepileupto(av,PicFrob(J,W2));
  }
	W = PicChord(J,W1,W2,0);
	for(y=1;y<=x;y++) W = gerepileupto(av,PicFrob(J,W));
	printf("Exiting Tors_worker\n");
  return W;
}

GEN TorsSpaceFrobEval(GEN J, GEN gens, GEN cgens, ulong l, GEN matFrob)
{
	pari_sp av = avma;
	ulong d,ld,ndone,ngens,nmodF,ntodo,i,j,k,n,ij,ik,xj,xk,x;
	GEN VmodF,ImodF,ZmodF,ImodF2,done,mfrob,c,todo,todon,args,Wj,Wk;
	struct pari_mt pt,pt2;
	GEN worker,worker2,res;
  long pending;

	d = lg(matFrob)-1; // Dim of rep
	ld = upowuu(l,d);
	VmodF = cgetg(ld+1,t_VEC);
	ImodF = cgetg(ld+1,t_VECSMALL);
	done = cgetg(ld+1,t_VEC);
	for(n=1;n<ld;n++)
		gel(done,n) = NULL;
	gel(done,ld) = mkvecsmall2(0,0);
	ndone = 1;
	nmodF = 0;
	ngens = lg(gens)-1;
	worker = strtofunction("TorsSpaceFrob_worker");
	args = cgetg(6,t_VEC);
  gel(args,1) = J;
	// matFrob, version Flm
	mfrob = cgetg(d+1,t_MAT);
  for(j=1;j<=d;j++)
  {
    gel(mfrob,j) = cgetg(d+1,t_VECSMALL);
    for(k=1;k<=d;k++)
      gel(mfrob,j)[k] = itos(liftint(gcoeff(matFrob,k,j)));
  }
	// Throw in generators and their Frob orbits
	c = cgetg(d+1,t_VECSMALL);
	for(n=1;n<=ngens;n++)
	{
		nmodF++;
		gel(VmodF,nmodF) = gel(gens,n);
		for(i=1;i<=d;i++)
			c[i] = itou(gmael(cgens,n,i));
		ImodF[nmodF] = i = c2i(c,l);
		gel(done,i) = mkvecsmall2(nmodF,0);
		ndone++;
		for(x=1;;x++)
		{
			i = ActOni(mfrob,i,l);
			if(gel(done,i)) break;
			gel(done,i) = mkvecsmall2(nmodF,x);
			ndone++;
		}
	}
	// Loop until everything is covered
	todo = cgetg(ld,t_VEC);
	while(ndone<ld)
	{
		printf("Main loop, touched %lu/%lu\n",ndone,ld);
		// TODO gerepile
		ntodo = 0;
		for(j=1;j<ld;j++)
		{
			for(k=j;k<=ld;k++)
			{
				if(gel(done,j)&&gel(done,k)&&gel(done,j)!=gen_m1&&gel(done,k)!=gen_m1) // Loop over pair of known pts
				{
					i = Chordi(j,k,l,d);
					if(gel(done,i)==NULL) // Do we already know their chord?
					{
						// Record the computation that needs to be done
						ij = gel(done,j)[1];
						xj = gel(done,j)[2];
						ik = gel(done,k)[1];
						xk = gel(done,k)[2];
						Wj = ij?gel(VmodF,ij):gen_0;
						Wk = ik?gel(VmodF,ik):gen_0;
						ntodo++;
						gel(todo,ntodo) = mkvecn(5,utoi(i),Wj,utoi(xj),Wk,utoi(xk));
						// Mark that we are about to get this Frob orbit
						while(gel(done,i)==NULL)
						{
							gel(done,i) = gen_m1;
							i = ActOni(mfrob,i,l);
						}
					}
				}
			}
		}
		// Execute operations in J in parallel
		mt_queue_start(&pt,worker);
		printf("Computing %lu new points\n",ntodo);
	  for(n=1;n<=ntodo||pending;n++)
  	{
			if(n<=ntodo)
			{
				todon = gel(todo,n);
				i = itos(gel(todon,1));
				for(j=2;j<=5;j++)
					gel(args,j) = gel(todon,j);
    		mt_queue_submit(&pt,i,args);
			}
			else mt_queue_submit(&pt,0,NULL);
    	res = mt_queue_get(&pt,(long*)&i,&pending);
    	if(res)
			{
				nmodF++;
				// Record new Frob orbit
				gel(VmodF,nmodF) = res;
				ImodF[nmodF] = i;
				// Mark the orbit as known
				for(x=0;gel(done,i)==gen_m1;x++)
				{
					gel(done,i) = mkvecsmall2(nmodF,x);
					ndone++;
					i = ActOni(mfrob,i,l);
				}
			}
  	}
  	mt_queue_end(&pt);
	}

	printf("Evaluating %lu points\n",nmodF);
	args = cgetg(3,t_VEC);
	gel(args,1) = J;
	ZmodF = cgetg(nmodF+1,t_VEC);
	worker2 = strtofunction("PicEval");
	mt_queue_start(&pt2,worker2);
	for(n=1;n<=nmodF||pending;n++)
	{
		if(n<=nmodF)
		{
			gel(args,2) = gel(VmodF,n);
			printf("Submit %lu\n",n);
			mt_queue_submit(&pt2,n,args);
		}
		else mt_queue_submit(&pt2,0,NULL);
		res = mt_queue_get(&pt2,(long*)&i,&pending);
		if(res)
		{
			printf("Getting %lu\n",i);
			gel(ZmodF,i) = res; // TODO gerepile
			printf("%p\n",gel(ZmodF,i));
			pari_printf("%Ps\n",gel(ZmodF,i));
			RgM_dimensions(gel(ZmodF,i),(long*)&j,(long*)&k);
			printf("%lu %lu\n",j,k);
		}
	}
	printf("U\n");
	pari_printf("%Ps\n",gel(ZmodF,1));
	mt_queue_end(&pt2);
	printf("End eval\n");
	pari_printf("%Ps\n",gel(ZmodF,1));

	ImodF2 = cgetg(nmodF+1,t_VECSMALL);
	for(n=1;n<=nmodF;n++) ImodF2[n] = ImodF[n];
	res = mkvec2(ZmodF,ImodF2);
	printf("End Tors\n");
	pari_printf("ImodF: %Ps\n",ImodF2);
	for(i=1;i<=nmodF;i++)
	{
		printf("%lu\n",i);
		printf("%p\n",gel(ZmodF,i));
		pari_printf("%Ps\n",gel(ZmodF,i));
	}
	pari_printf("ZmodF: %Ps\n",ZmodF);
	return gerepilecopy(av,res); // TODO do not copy
}

GEN PolExpID(GEN Z, GEN T, GEN pe) /* bestappr of prod(x-z), z in Z */
{
  pari_sp av = avma;
  GEN f,a;
  f = FqV_roots_to_pol(Z,T,pe,0);
  if(poldegree(f,varn(T))>0) pari_err(e_MISC,"Irrational coefficient: %Ps",f);
  f = simplify_shallow(f);
  f = gmodulo(f,pe);
  a = bestappr(f,NULL);
  return gerepilecopy(av,mkvecn(3,Z,f,a));
}

GEN OnePol(GEN N, GEN D, GEN ImodF, GEN Jfrobmat, ulong l, GEN QqFrobMat, GEN T, GEN pe)
{ /* Actually returns a vector of n1*n2 pols (all elem. symm. fns) */
  pari_sp av = avma, av1;
  GEN R,Z,F,Fi,z;
  ulong d,ld,j,k,n,i1,i2,i;
  long n1,n2,n12;
  n = lg(N);
  RgM_dimensions(gel(N,1),&n2,&n1);
  /* N, D vectors of length n-1 of n2*n1 matrices
   * N[i]/D[i] = value at pt indexed by ImodF[i]
   * Get the others by applying Frob */
  d = lg(Jfrobmat)-1;
  ld = upowuu(l,d);
  n12 = n1*n2;
  R = cgetg(n,t_VEC);
  Z = cgetg(n12+1,t_VEC);
  for(k=1;k<n;k++)
  {
    i=1;
    for(i1=1;i1<=n1;i1++)
    {
      for(i2=1;i2<=n2;i2++)
      {
        gel(Z,i) = Fq_mul(gmael3(N,k,i1,i2),gmael3(D,k,i1,i2),T,pe);
        i++;
      }
    }
    gel(R,k) = FqV_roots_to_pol(Z,T,pe,0);
  }
  R = gerepileupto(av,R);
  F = cgetg(n12+1,t_VEC);
  Z = cgetg(ld,t_VEC);
  for(i=0;i<n12;i++)
  {
    av1 = avma;
    for(j=1;j<ld;j++) gel(Z,j) = NULL; /* Mark roots as unknown */
    for(k=1;k<n;k++)
    {
      z = polcoef(gel(R,k),i,0);
      j = ImodF[k];
      for(;;)
      {
        gel(Z,j) = z;
        j = ActOni(Jfrobmat,j,l);
        if(gel(Z,j)) break;
        z = Frob(z,QqFrobMat,T,pe);
      }
    }
    Fi = PolExpID(Z,T,pe);
    if(n12>1) Fi = gerepileupto(av1,Fi);
    gel(F,i+1) = Fi;
  }
  return gerepileupto(av,F);
}

GEN AllPols(GEN Z, ulong l, GEN JFrobMat, GEN QqFrobMat, GEN T, GEN pe, GEN p, long e)
{
  pari_sp av = avma;
  GEN F,ImodF,Jfrobmat,Ft,F1,f,pols,args;
  ulong d,nF,lF,npols,n,i,j,i1,i2,m,k;
  long n1,n2,n12;
  struct pari_mt pt;
  GEN worker,done;
  long pending,workid;

	printf("Into AllPols\n");
	F = gel(Z,1);
	ImodF = gel(Z,2);
	d = lg(JFrobMat)-1;
  Jfrobmat = cgetg(d+1,t_MAT); /* JFrobMat, version Flm */
  for(j=1;j<=d;j++)
  {
    gel(Jfrobmat,j) = cgetg(d+1,t_VECSMALL);
    for(k=1;k<=d;k++)
      gel(Jfrobmat,j)[k] = itos(liftint(gcoeff(JFrobMat,k,j)));
  }
  nF = lg(F); /* Number of vectors */
  RgM_dimensions(gel(F,1),&n2,&n1);
  n12 = n1*n2;
  lF = lg(gmael3(F,1,1,1))-1; /* Size of each vector */
  /* F = list of nF-1 matrices of size n2*n1 of vectors of size lF */
  Ft = cgetg(lF,t_VEC);
  /* Ft[j,i,i1,i2] = F[i,i1,i2,j] */
  for(j=1;j<lF;j++)
  {
    gel(Ft,j) = cgetg(nF,t_VEC);
    for(i=1;i<nF;i++)
    {
      gmael(Ft,j,i) = cgetg(n1+1,t_MAT);
      for(i1=1;i1<=n1;i1++)
      {
        gmael3(Ft,j,i,i1) = cgetg(n2+1,t_COL);
        for(i2=1;i2<=n2;i2++)
          gmael4(Ft,j,i,i1,i2) = gmael4(F,i,i1,i2,j);
      }
    }
  }
  F1 = cgetg(lF,t_VEC);
  npols = 0;
  for(i=1;i<lF;i++) /* Find the i such that the ith coord of all the vectors in all the matrices are invertible, and store there inverses */
  {
    npols++;
    gel(F1,i) = cgetg(nF,t_VEC);
    for(j=1;j<nF;j++)
    {
      gmael(F1,i,j) = cgetg(n1+1,t_MAT);
      for(i1=1;i1<=n1;i1++)
      {
        gmael3(F1,i,j,i1) = cgetg(n2+1,t_COL);
        for(i2=1;i2<=n2;i2++)
        {
          f = gmael4(F,j,i1,i2,i);
          if(ZX_is0mod(f,p))
          {
            gel(F1,i) = NULL;
            npols--;
            i2=n2+1;i1=n1+1;j=nF;
          }
          else
            gmael4(F1,i,j,i1,i2) = ZpXQ_inv(f,T,p,e); /* TODO do it later in case we give up this i */
        }
      }
    }
  }
  pending = 0;
  worker = strtofunction("OnePol");
  args = cgetg(9,t_VEC);
  gel(args,3) = ImodF;
  gel(args,4) = Jfrobmat;
  gel(args,5) = utoi(l);
  gel(args,6) = QqFrobMat;
  gel(args,7) = T;
  gel(args,8) = pe;
  npols *= (lF-2)*n12;
  pols = cgetg(npols+1,t_VEC);
  mt_queue_start(&pt,worker);
  done = NULL;
  for(i=j=m=n=1;i<lF||pending;n++,j++)
  {
    if(j==lg(F1))
    {
      j=1;
      i++;
    }
    if(gel(F1,j)==NULL || i==j) continue; /* Skip if denom=0 or if denom=num */
    if(i<lF)
    {
      gel(args,1) = gel(Ft,i);
      gel(args,2) = gel(F1,j);
      mt_queue_submit(&pt,n,args);
    }
    else mt_queue_submit(&pt,0,NULL);
    done = mt_queue_get(&pt,&workid,&pending);
    if(done)
    {
      for(k=1;k<=n12;k++)
      {
        gel(pols,m) = gel(done,k);
        m++;
      }
    }
  }
  mt_queue_end(&pt);
  return gerepileupto(av,pols);
}
