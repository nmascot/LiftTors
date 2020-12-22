CycleDecomp(s)=
{ \\ Decomp of perm [s(i)] into disjoint cycles
	my(n=#s,done=vector(n),j,cyc=List(),c);
	for(i=1,n,
		if(done[i]==0,
			j=i;
			c=List();
			until(j==i,
				listput(c,j);
				done[j]=1;
				j=s[j]
			);
			listput(cyc,Vecsmall(c));
		)
	);
	Vec(cyc);
}

PermRecomp(s)=
{ \\ Prod of disjoint cycles -> Perm [s[i]]
	my(N,P);
	N = vecmax(apply(vecmax,s));
	P = vector(N,i,i);
	for(i=1,#s,
		N = #s[i];
		for(j=1,N-1,
			P[s[i][j]] = s[i][j+1]
		);
		P[s[i][N]] = s[i][1]
	);
	Vecsmall(P);
}

SubPerm(s,M)=
{ \\ Perm [s[i]] acting on 1..N, M<=N -> [t[i]], T
	\\ T subset (possibly reordered) of 1..N stable under s, #T>=M but close, t perm induced on T
	my(N,c,nc,L,l,r,t,T,d,m=0);
	N = #s;
	c = CycleDecomp(s); \\ Cycle decomp
	L = apply(x->#x,c); \\ Length of cycles
	r = vecsort(L,,5); \\ Permutation sorting L
	nc = #c; \\ # of cycles
	c = List(vector(nc,i,c[r[i]])); \\ Sort c and L
	L = List(vector(nc,i,L[r[i]])); \\ so that L decreasing
	T = List();
	t = List();
	while(m<M, \\ m=#T
		for(i=1,nc, \\ Look for the largest cycle we can throw in
			l = L[i];
			if(l+m<=M, \\ Does this cycle fit?
				d = c[i];
				listpop(c,i);
				listpop(L,i);
				nc--;
				for(j=1,l-1,
					listput(T,d[j]);
					listput(t,m+j+1);
				);
				listput(T,d[l]);
				listput(t,m+1);
				m+=l;
				next(2);
			)
		);
		\\ No cycle fits, so we have to exceed M
		d = c[nc];
    for(j=1,l-1,
      listput(T,d[j]);
      listput(t,m+j+1);
    );
    listput(T,d[l]);
    listput(t,m+1);
    m+=l;
	);
	[Vecsmall(t),Vecsmall(T)];
}

SubPerm_multi(S,M)=
{ \\ Perms S=[s[i]] acting on 1..N, M<=N -> [t[i]], T
  \\ T subset (possibly reordered) of 1..N stable under S, #T>=M but close, t perm induced on T
	my(N,Orbs,seen,iseen,P,Orb,SOrb,nOrb,n,m,find,Sub,SubS,l);
	N = #S[1];
	seen = vector(N);
	iseen = 0;
	Orbs=[];
	while(iseen<N,
		iseen++;
		if(seen[iseen],next);
		\\ Start new orbit
		P=iseen;
		Orb=[P];
		SOrb = vector(#S,i,vector(N));
		nOrb=1;
		n=1;
		while(n<=nOrb, \\ Loop over orbit
			for(i=1,#S, \\ for each perm
				if(SOrb[i][n]==0, \\ do we know what happens to this point by this perm?
					m=n; \\ No. Then let us see.
					P=Orb[m]; \\ This point
					while(1,
						P = S[i][P]; \\ Image by this perm.
						seen[P]=1; \\ Mark this point as visited
						find = select(Q->Q==P,Orb,1); \\ Is it in Orb?
						if(#find==1,
							find = find[1]; \\ It is in Orb, in this position.
							SOrb[i][m] = find; \\ Record info in perm
							break
						);
						\\ It is not in Orb yet
						Orb=concat(Orb,[P]); \\ Add it
						nOrb++; \\ Orb size increases
						SOrb[i][m] = nOrb; \\ Record info in perm
						m = nOrb;
					)
				)
			);
			n++
		);
		\\ Found a complete orbit.
		Orb=Orb[1..nOrb]; \\ Shorten vectors
		SOrb = apply(s->s[1..nOrb],SOrb);
		\\ Record new orbit
		Orbs = concat(Orbs,[[Orb,SOrb]])
	);
	\\ So we have split the set into orbs.
	Orbs = vecsort(Orbs,o->#o[1],4); \\ Sort orbs by decreasing size
	\\ Now select orbits to form subset of size >=M
	m=0;
	Sub=[];
	SubS=vector(#S,i,[]);
	for(n=1,#Orbs,
		l=#Orbs[n][1]; \\ Size of orbit
		if(m+l<=M, \\ does it fit?
			Sub=concat(Sub,Orbs[n][1]);
			for(i=1,#S,
				SubS[i] = concat(SubS[i],apply(x->x+m,Orbs[n][2][i]))
			);
			Orbs[n]=0; \\ Mark as used
			m+=l
		)
	);
	\\ Maybe still m<M. So throw in orbits form smaller to larger
	n=1+#Orbs;
	while(m<M,
		n--;
		if(Orbs[n]==0,next);
		l=#Orbs[n][1]; \\ Size of orbit
    Sub=concat(Sub,Orbs[n][1]);
    for(i=1,#S,
    	SubS[i] = concat(SubS[i],apply(x->x+m,Orbs[n][2][i]))
    );
    Orbs[n]=0; \\ Mark as used
    m+=l
	);
	[Vecsmall(Sub),apply(s->Vecsmall(s),SubS)];
}
