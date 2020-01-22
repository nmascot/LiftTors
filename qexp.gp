E1qexp(v,N,z,B,Tpe,var='q)=
/* v=[c,d] mod N, z primitive Nth root of 1: q-exp of E_1^[c,d] up to O(qN^B) */
{
	my([c,d]=v,E,a,b,[T,pe,p,e]=Tpe);
	E=vector(B);
	c=c%N;
  d=d%N;
  E[1] = if(c==0,
    \\1/2*(1+z^d)/(1-z^d),
		ZpXQ_div(liftall(1+z^d),liftall(2*(1-z^d)),T,pe,p,e),
    Mod(1/2-c/N,pe);
    \\1/2-c/N;
  );
	for(n=1,B-1,
		fordiv(n,a,
			b= n/a;
			if(Mod(a,N)==Mod(c,N),
				E[n+1] += z^((b*d)%N)
			);
			if(Mod(a,N)==-Mod(c,N),
        E[n+1] -= z^((-b*d)%N)
      );
		)
	);
	Ser(E,var);
}

E1coef(v,N,z,n,Tpe)=
{
	my([c,d]=v,a,b,[T,pe,p,e]=Tpe,an);
	c=c%N;
  d=d%N;
	if(n==0,
		an = if(c==0,
    	ZpXQ_div(liftall(1+z^d),liftall(2*(1-z^d)),T,pe,p,e),
    	Mod(1/2-c/N,pe)
		)
	,
		an = 0;
		fordiv(n,a,
      b = n/a;
      if(Mod(a,N)==Mod(c,N),an += z^((b*d)%N));
      if(Mod(a,N)==-Mod(c,N),an -= z^((-b*d)%N))
		)
  );
	an;
}

E1qexp_formal(s,N,z=Mod('z,polcyclo(N,'z)),B,var='q)=
/* s=[c,d] mod N, z primitive Nth root of 1: q-exp of E_1^[c,d] up to O(qN^B) */
{
  my([c,d]=s,E,a,b);
  E=vector(B);
  c=c%N;
  d=d%N;
  E[1] = if(c==0,
    1/2*(1+z^d)/(1-z^d),
    \\ZpXQ_div(liftall(1+z^d),liftall(2*(1-z^d)),T,pe,p,e),
    \\Mod(1/2-c/N,pe);
    1/2-c/N;
  );
  for(n=1,B-1,
    fordiv(n,a,
      b= n/a;
      if(Mod(a,N)==Mod(c,N),
        E[n+1] += z^((b*d)%N)
      );
      if(Mod(a,N)==-Mod(c,N),
        E[n+1] -= z^((-b*d)%N)
      );
    )
  );
  Ser(E,var);
}

TrE_qexp(v,N,TH,M,w,z,B,Tpe)=
{ \\ q-exp of Tr_H(prod_i(E1^v[i]))|M in terms of qw up to O(qw^B)
	my(f,qprec,fw);
	if(B==0,return([]));
	qprec = (B-1)*N/w+1;
	f=
		sum(g=1,#TH,
			prod(i=1,#v,
				E1qexp(v[i]*TH[g]*M,N,z,qprec,Tpe,'x)
			)
		);
	/* DEBUG */
	for(n=0,qprec-1,
		if(Mod(n,N/w) && polcoeff(f,n),error("Bad qexp"))
	);
	\\ Contract qN = qw^(N/w) -> qw
	vector(B,n,polcoeff(f,(n-1)*N/w));
}
	

M4RRsubspace(N,Hlist,TH,z,Tpe,M4gens,Cusps,D)=
{ \\ TODO version with forms def / QQ
	my(nCusps,nForms,nTH=#TH,d,A,i,M,w,f,[T,pe,p,e]=Tpe);
	nCusps = #Cusps;
	nForms = #M4gens;
	d = vecsum(D); \\ Degree of divisor
	A = matrix(d,nForms); \\ Matrix of eqn = q-exps
	i = 0;
	for(s=1,nCusps,
		if(D[s]==0,next);
		[M,w] = GammaHCuspData(Cusps[s],N,Hlist); \\ M.oo = s, width w
		for(j=1,nForms,
			f = TrE_qexp(M4gens[j],N,TH,M,w,z,D[s],Tpe); \\ qexp of j-th gen of M4 at cusp s
			for(k=1,D[s],A[i+k,j] = f[k]) \\ Copy its coeffs
		);
		i += D[s];
	);
	matkerpadic(liftall(A),T,p,e);
}

