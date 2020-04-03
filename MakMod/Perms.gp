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

