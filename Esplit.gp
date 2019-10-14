Esplit(p,N,d)= \\ Look for an ell. curve / Fp such tht E[N] splits over Fq, q=p^d.
{
  my(q=p^d,q2=(q-1)/2,T=ffinit(p,d,'t),a,b,E,ap,X,g,M,L,V,l,v,lv,c,D);
  if(Mod(p^d,N)!=1,error("Impossible by Weil pairing"));
	/*m = valuation(N,2);
	M = N/2^m;*/
  while(1,
    print("---- New curve ----");
    a=random(Mod(1,p));
    b=random(Mod(1,p));
    if(4*a^3+27*b^2==0,
			print("Singular");
			next
		);
    E=ellinit([a,b]);
		ap = ellap(E);
		g=gcd(ap^2-4*p,N);
		L=factor(g)[,1]; \\ Primes l s.t. Frob_p not 1/2simple on l-tors
		V=apply(l->valuation(N,l),L); \\ Multiplicities in N
		M=N/prod(i=1,#L,L[i]^V[i]); \\ Part of N where Frob 1/2simple
		if(M>1, \\ N=M*prod(l_i^v_i)
			if(polsym(Mod('x^2-ap*'x+p,M),d)[d+1]!=2, \\ Check if Frob^d=1 on E[M] by testing a^d+b^d==2
        print("Frob^",d," not trivial on E[",M,"]");
        next
      )
		);
		for(i=1,#L, \\ Check if Frob^d unipotent on E[l^v]
			l=L[i];
			v=V[i];
			\\print("Checking ",l,"^",v);
			lv=l^v;
			if(l!=2,
				c=ap/Mod(2,lv); \\ Frob = c*unipotent on E[l^v]
				if(c^d!=1, \\ Check if Frob^d unipotent on E[l^v]
					print("Frob^",d," not unipotent on E[",l,"^",v,"]");
					next(2)
				);
			);
			if(Mod(d,lv)==0,next); \\ d kills unipotents in this case
			D = elldivpol(E,lv); \\ Check manually that E[l^v] def over Fq
			X=polrootsmod(D,[T,p]);
			if(#X<poldegree(D),
				print("E[",l,"^",v,"]/+-1 not split: ",#X," roots out of ",poldegree(D));
				next(2)
			);
			if(l!=2,
				for(j=1,#X,
					if((X[i]^3+a*X[i]+b)^q2!=1,
						print("E[",l,"^",v,"] not split");
						next(3)
					)
				)
			)
		);
		return([a,b])
  );
}

GetPt(x,f)=[x,sqrt(subst(f,'x,x))];

GetLocalBasis(E,t1,f,l,v)=
{
	my(D,X,P,Q,z);
	D = elldivpol(E,l^v);
	X = polrootsmod(D,t1);
	while(1,
		P = GetPt(X[random(#X)+1],f);
		Q = GetPt(X[random(#X)+1],f);
		z = ellweilpairing(E,P,Q,l^v);
		if(z^(l^(v-1))!=1,
			return([P,Q,z])
		)
	);
}

GetBasis(E,t1,f,N)=
{
	my(faN,Pk,Qk,zk,P,Q,z,l,v);
	faN = factor(N);
	P = [0];
	Q = [0];
	z = 1;
	for(k=1,#faN~,
		[l,v] = faN[k,];
		[Pk,Qk,zk] = GetLocalBasis(E,t1,f,l,v);
		P = elladd(E,P,Pk);
		Q = elladd(E,Q,Qk);
		z *= zk^(N/l^v)
	);
	[P,Q,z];
}
