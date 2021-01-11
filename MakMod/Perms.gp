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

SubPerm_multi(S,M)=
{ /* Perms S=[s[i]] acting on 1..N, M<=N -> [T,ST]
   T subset (possibly reordered) of 1..N stable under S and with #T>=M but close, 
	ST perms induced by S on T */
	my(nS=#S,N,Orbs,seen,iseen,P,Orb,SOrb,nOrb,n,m,find,Sub,SubS,l);
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
			for(i=1,nS, \\ for each perm
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
	\\ Shuffle them randomly
  for(n=1,#Orbs,
    i=1+random(#Orbs);
    j=1+random(#Orbs);
    Orb = Orbs[i];
    Orbs[i]=Orbs[j];
    Orbs[j]=Orb;
  );
	\\ Sort them by decreasing size
	Orbs = vecsort(Orbs,o->#o[1],4); \\ Sort orbs by decreasing size
	\\ Now select orbits to form subset of size >=M
	m=0;
	Sub=[];
	SubS=vector(nS,i,[]);
	for(n=1,#Orbs,
		l=#Orbs[n][1]; \\ Size of orbit
		if(m+l<=M, \\ does it fit?
			Sub=concat(Sub,Orbs[n][1]);
			for(i=1,nS,
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
    for(i=1,nS,
    	SubS[i] = concat(SubS[i],apply(x->x+m,Orbs[n][2][i]))
    );
    Orbs[n]=0; \\ Mark as used
    m+=l
	);
	[Vecsmall(Sub),apply(s->Vecsmall(s),SubS)];
}

PermConcat(s,t)=
{
  my(n=#s,u=t,m=#t);
  for(i=1,m,u[i]+=n);
  concat([s,u]);
}
