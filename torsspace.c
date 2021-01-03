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

GEN TorsSpaceFrob_worker(GEN W1, GEN X1, GEN W2, GEN X2, GEN J)
{
	pari_sp av = avma;
	GEN W;
	ulong x1,x2;
	x1 = itou(X1);
	if(W2==gen_0)
	{
		W = PicNeg(J,W1,0);
		W = PicFrobPow(J,W,x1);
		return gerepileupto(av,W);
	}
	x2 = itou(X2);
	W1 = PicFrobPow(J,W1,x1);
	W2 = PicFrobPow(J,W2,x2);
	W = PicChord(J,W1,W2,0);
  return gerepileupto(av,W);
}

GEN M2Flm(GEN A, ulong l, ulong m, ulong n)
{ /* Reduce m*n matrix A mod l -> ve of vecsmalls */
	GEN a;
	ulong j,k;
	a = cgetg(n+1,t_MAT);
  for(j=1;j<=n;j++)
  {
    gel(a,j) = cgetg(m+1,t_VECSMALL);
    for(k=1;k<=m;k++)
      gel(a,j)[k] = itos(liftint(gcoeff(A,k,j)));
  }
	return a;
}

GEN TorsSpaceFrobEval(GEN J, GEN gens, GEN cgens, ulong l, GEN matFrob, GEN matAuts)
{
	pari_sp av = avma;
	ulong d,ld,ndone,ngens,nmodF,newmodF,new,ntodo,i,j,k,m,n,ij,ik,xj,xk,x,nAuts,a;
	GEN vJ,vW,VmodF,ImodF,ZmodF,ImodF2,done,mfrob,mauts,c,todo,W,Wj,Wk;
	struct pari_mt pt;
	GEN worker,res;
  long pending;

	d = lg(matFrob)-1; // Dim of rep
	ld = upowuu(l,d);
	VmodF = cgetg(ld+1,t_VEC); // list of W's mod Frob
	ImodF = cgetg(ld+1,t_VECSMALL); // list of indices of these W's
	done = cgetg(ld+1,t_VEC); // All space: each entry is either
	// NULL: we don't have this pt yet
	// gen_m1: we will get this pt by applying Auts and Frob
	// Vecmall([n,m]): this pt is Frob^m * VmodF[n]
	// Vecsmall([0,0]): this is the origin
	for(n=1;n<ld;n++) // Initalise
		gel(done,n) = NULL;
	gel(done,ld) = mkvecsmall2(0,0);
	ndone = 1; // Number of pts; we stop when this reaches ld
	nmodF = 0; // Number of pts mod Frob
	ngens = lg(gens)-1;
	// Prepare cloure for parallel code
	vJ = cgetg(2,t_VEC);
  gel(vJ,1) = J;
	worker = snm_closure(is_entry("TorsSpaceFrob_worker"),vJ);
	// matFrob and matAuts, version Flm
	mfrob = M2Flm(matFrob,l,d,d);
	nAuts = lg(matAuts);
	mauts = cgetg(nAuts,t_VEC);
	for(a=1;a<nAuts;a++) gel(mauts,a) = M2Flm(gel(matAuts,a),l,d,d);
	// Now we are ready to begin
	// We will always keep done closed under Frob
	
	// Throw in generators and their Frob orbits
	c = cgetg(d+1,t_VECSMALL);
	for(n=1;n<=ngens;n++)
	{
		nmodF++; // new Frob orbit
		gel(VmodF,nmodF) = gel(gens,n); // store W in VmodF
		for(i=1;i<=d;i++) // Convert coords into index
			c[i] = itou(gmael(cgens,n,i));
		ImodF[nmodF] = i = c2i(c,l);
		gel(done,i) = mkvecsmall2(nmodF,0); // store [nmodF,0] in done[i]
		ndone++; // got new pt
		for(x=1;;x++) // go through Frob orbit
		{
			i = ActOni(mfrob,i,l); // Move index by Frob
			if(gel(done,i)) break; // end of orbit?
			gel(done,i) = mkvecsmall2(nmodF,x); // Record [nmodF,x], meaning x-th translate of this Frob orbit, in done
			ndone++; // got new pt
		}
	}
	// Loop until everything is covered
	todo = cgetg(ld,t_VEC); // list of operations to be executed in parallel; used later
	while(ndone<ld)
	{
		printf("Main loop, touched %lu/%lu\n",ndone,ld);
		// Use Auts as much as possible
		do
		{
			newmodF = 0; // number of new Fob orbs we get at this iteration. If it's still 0 in the end, we are done with Auts.
			for(i=1;i<=nmodF;i++) // Try to apply to all pts mod Frob...
			{ // TODO only new orbits need to be checked
				j = ImodF[i];
				for(a=1;a<nAuts;a++) // ... all auts
				{
					k = ActOni(gel(mauts,a),j,l);
					W = gel(done,k);
					if(W==NULL || W==gen_m1) // does this yield a new pt, and therefore a new Frob orbit?
					{
						W = gel(VmodF,i); // this new pt, before Aut
						W = PicAut(J,W,a); // after Aut
						nmodF++; // new Frob orbit
						newmodF++; // new Frob orbit at this iteration
						gel(VmodF,nmodF) = W; // record representative of Frob orbit
						ImodF[nmodF] = k; // and its index
						gel(done,k) = mkvecsmall2(nmodF,0); // store [nmodF,0] in done
						ndone++; // new pt
						for(x=1;;x++) // go through Frob orb
						{
				      k = ActOni(mfrob,k,l); // Move index by Frob
							W = gel(done,k);
      				if(W && W!=gen_m1) break; // end of orbit?
      				gel(done,k) = mkvecsmall2(nmodF,x); // Record [nmodF,x], meaning x-th translate of this Frob orbit, in done
      				ndone++; // new pt
    				}
					}
				}
			}
			//printf("Applying auts gave %lu new Frob orbits\n",newmodF);
		} while(newmodF);
		printf("Done with auts, now touched %lu/%lu\n",ndone,ld);
		if(ndone==ld) break; // Are we done?
		// TODO gerepile?
		// Now use group law
		// Prepare list of operations, we will run them in parallel later
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
						gel(todo,ntodo) = mkvec2(mkvecn(4,Wj,utoi(xj),Wk,utoi(xk)),utoi(i));
						// Mark that we are about to get this Frob orbit...
						while(gel(done,i)==NULL)
						{
							gel(done,i) = gen_m1;
							i = ActOni(mfrob,i,l);
						}
						// ...and all it translates by auts
						do
				    {
      				new = 0; // number of new pts we get at this iteration. If it's still 0 in the end, we are done with Auts.
      				for(m=1;m<=ld;m++) // Try to apply to all pts...
      				{
								if(gel(done,m)==NULL) continue; // ... that we already have or will have ...
        				for(a=1;a<nAuts;a++) // ... all auts
        				{
          				n = ActOni(gel(mauts,a),m,l);
          				while(gel(done,n)==NULL) // will this yield a new pt, and therefore a new Frob orbit?
          				{
										gel(done,n) = gen_m1;
										new++;
										n = ActOni(mfrob,n,l);
									}
								}
            	}
    				} while(new);
					}
				}
			}
		}
		// Execute operations in J in parallel
		printf("Computing %lu new points\n",ntodo);
		mt_queue_start_lim(&pt,worker,ntodo);
	  for(n=1;n<=ntodo||pending;n++)
  	{
			if(n<=ntodo) mt_queue_submit(&pt,itos(gmael(todo,n,2)),gmael(todo,n,1));
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
	vW = cgetg(2,t_VEC);
	ZmodF = cgetg(nmodF+1,t_VEC);
	worker = snm_closure(is_entry("RREval_worker"),vJ);
	mt_queue_start_lim(&pt,worker,nmodF);
	for(n=1;n<=nmodF||pending;n++)
	{
		if(n<=nmodF)
		{
			gel(vW,1) = gel(VmodF,n);
			mt_queue_submit(&pt,n,vW);
		}
		else mt_queue_submit(&pt,0,NULL);
		res = mt_queue_get(&pt,(long*)&i,&pending);
		if(res)
			gel(ZmodF,i) = res; // TODO gerepile
	}
	mt_queue_end(&pt);

	ImodF2 = cgetg(nmodF+1,t_VECSMALL);
	for(n=1;n<=nmodF;n++) ImodF2[n] = ImodF[n];
	res = mkvec2(ZmodF,ImodF2);
	return gerepilecopy(av,res); // TODO do not copy
}

GEN PolExpID(GEN Z, GEN T, GEN pe) /* bestappr of prod(x-z), z in Z */
{
  pari_sp av;
  GEN res,f;
	res = cgetg(4,t_VEC);
	av = avma;
  f = FqV_roots_to_pol(Z,T,pe,0);
  if(poldegree(f,varn(T))>0) pari_err(e_MISC,"Irrational coefficient: %Ps",f);
  f = simplify_shallow(f);
  f = gmodulo(f,pe);
	f = gerepileupto(av,f);
	gel(res,2) = f;
	gel(res,3) = bestappr(f,NULL);
	gel(res,1) = gmodulo(gmodulo(Z,T),pe);
	return res;
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
  pari_sp av = avma, avj;
  GEN F,ImodF,Jfrobmat,Ft,F1,f,pols,args;
  ulong d,nF,lF,npols,n,i,j,j0,i1,i2,m,k;
	int All0;
  long n1,n2,n12;
  struct pari_mt pt;
  GEN worker,done;
  long pending,workid;

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
  /* Ft[j,i,i1,i2] = F[i,i1,i2,j], keeping only the j such that F[.,.,.,j] not identically 0 */
  for(j=j0=1;j<lF;j++)
  {
		All0 = 1;
		avj = avma;
    gel(Ft,j0) = cgetg(nF,t_VEC);
    for(i=1;i<nF;i++)
    {
      gmael(Ft,j0,i) = cgetg(n1+1,t_MAT);
      for(i1=1;i1<=n1;i1++)
      {
        gmael3(Ft,j0,i,i1) = cgetg(n2+1,t_COL);
        for(i2=1;i2<=n2;i2++)
				{
          f = gmael4(Ft,j0,i,i1,i2) = gmael4(F,i,i1,i2,j);
					if(gequal0(f)==0) All0 = 0;
				}
      }
    }
		if(All0) avma = avj; /* Drop this j */
    else j0++;
  }
	printf("Reducing lF from %lu to %lu\n",lF-1,j0-1);
	lF = j0;
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
          f = gmael4(Ft,i,j,i1,i2);
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
  mt_queue_start_lim(&pt,worker,npols);
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
