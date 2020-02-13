read("ZN.gp");

SL2lift(M)= \\ Finds M' in SL2(Z), M=M' mod (|M|-1)
{
  my(U,V,D,a,b);
  [U,V,D]=matsnf(M,1); \\ U*M*V = D = diag(a,b), so |M|=ab
  a=D[1,1];
  b=D[2,2];
  U^-1*[1,-1;1-b,b]*[1,0;1-a,1]*[1,b;0,1]*V^-1;
} \\ [1,-1;1-b,b]*[1,0;1-a,1]*[1,b;0,1] = [a, ab-1; 1-ab, b+(1-ab)b]

{
BotToSL2(s,N)= \\ Finds [*,*;c';d'] in SL2(Z) with c~c', d~d' mod N
  my(c,d,u,v,g,h);
  [c,d]=s;
  [u,v,g]=gcdext(c,d);
  h=lift(1/Mod(g,N));
  SL2lift([v*h,-u*h;c,d]);
}

GammaHCuspData(s,N,Hlist)= \\ (c,d) -> [a,b;c,d],width
{
  my(c=s[1],M,a,g,x=1);
  g = gcd(c^2,N);
  g = N/g;
  M = BotToSL2(s,N);
  a = M[1,1];
  while(#select(y->Mod(y,N)==Mod(1+a*c*g*x,N),Hlist,1)==0,x++);
  [M,g*x];
}

ANH(N,Hlist)=
{ \\ Find all (u,v) s.t. gcd(u,v,N)=1 / H
	my(A=List(),tag=matrix(N,N),nH=#Hlist,n=0);
  for(u=0,N-1,
    for(v=0,N-1,
      if(GetCoef(tag,[u,v]),next);
      if(gcd([u,v,N])==1,
				n+=1;
        listput(A,[u,v]);
				for(h=1,nH,tag[ZNnorm(Hlist[h]*u,N),ZNnorm(Hlist[h]*v,N)]=n);
      )
    )
  );
  [Vec(A),tag];
}

GammaHCusps(N,Hlist)= \\ Reps (c,d) of all cusps of GammaH, plus data to find rep in equiv class
{
	/* * Reps (c,d) of all cusps of GammaH
		 * Galois orbits
		 * Vector of indices of cusps such that there is M = [*,*;c,d] in SL(2,Z) such that f|M has rat coefs for all f def/Q
		 * Matrices [*,*;c,d] in SL(2,Z), satifying above condition if possible
		 * Galois orbits
		 * Widths
		 * Tags: (c',d') -> index of equivalent representative
	*/
  my(Cusps=List(),CuspsGal=List(),Qqexp=List(),Mats=List(),Widths=List(),GalOrb,tag=matrix(N,N),nH=#Hlist,n=0,g,g2,h,M,a,Q,w);
	for(c=0,N-1, \\ c in Z/NZ
		g = gcd(c,N);
		g2 = N/gcd(c^2,N);
		GalOrb = List(); \\ Galois orbits: two cups are Galois-conj iff. they have the same c mod H
    for(d=0,g-1, \\ d in (Z/gZ)*
      if(GetCoef(tag,[c,d]),next); \\ Already visited this class
      if(gcd(g,d)==1, \\ Found new class
				listput(Cusps,[c,d]); \\ Record that class
				listput(GalOrb,[c,d]); \\ Add it to its Galois orbit
				n++; \\ Index of that class
				for(x=0,N/g-1, \\ Mark equivalent cusps
					for(i=1,nH,
						h = Hlist[i];
						tag[ZNnorm(h*c,N),ZNnorm(h*d+x*g,N)]=n \\ up to H
					)
				);
				M = BotToSL2([c,d],N); \\ Matrix [a,b;c,d]
				\\ The other choices are [1,t;0,1]*M
				if(Mod(2*c*d,N), \\ If qexps / Q, then necessarily N|2cd
					listput(Mats,M)
				,
					\\ Qqexp iff can choose t so that for all invertible x, ad(x-1)+1 in H
					for(i=1,N,
						Q = 1;
  					a = M[1,1];
						for(x=2,N-1,
							if(gcd(x,N)>1,next);
							if(#select(y->Mod(y,N)==Mod(a*d*(x-1)+1,N),Hlist,1)==0,Q=0;break)
						);
						if(Q,break);
						M=M*[1,1;0,1]
					);
          if(Q,listput(Qqexp,n));
          listput(Mats,M);
				);
				w=1; \\ Compute width: g2 * min w such that 1+acg2w in H
  			while(#select(y->Mod(y,N)==Mod(1+a*c*g2*w,N),Hlist,1)==0,w++);
				listput(Widths,g2*w);
			)
		);
		if(#GalOrb,listput(CuspsGal,Vec(GalOrb)));
	);
	[Vec(Cusps),vecsort(Vec(CuspsGal),x->#x,4),Vec(Qqexp),Vec(Mats),Vec(Widths),tag];
}

GammaHmodN(N,Hlist)= \\ reps of Gamma_H(N) / Gamma(N) in SL2(Z)
{
	my(nH=#Hlist,G=vector(nH*N),i=1,Mh);
	for(h=1,nH,
		Mh = Hlist[h];
		Mh = liftint(matdiagonal([1/Mh,Mh]));
		Mh = SL2lift(Mh); \\ Mh = [h^-1,*;N*,h]
		for(x=0,N-1,
			G[i] = [1,x;0,1]*Mh;
			i++
		)
	);
	G;
}
