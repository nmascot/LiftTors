read("ZN.gp");

SL2lift(M)= \\ Finds M' in SL2(Z), M=M' mod (|M|-1)
{
  my(U,V,D,a,b);
  [U,V,D]=matsnf(M,1); \\ U*M*V = D = diag(a,b), so |M|=ab
  a=D[1,1];
  b=D[2,2];
  U^-1*[1,-1;1-b,b]*[1,0;1-a,1]*[1,b;0,1]*V^-1;
} \\ [1,-1;1-b,b]*[1,0;1-a,1]*[1,b;0,1] = [a, ab-1; 1-ab, b+(1-ab)b]

ANH(N,Hlist)=
{ \\ Find all (u,v) s.t. gcd(u,v,N)=1 / H
	my(A=List(),done=matrix(N,N),nH=#Hlist);
  for(u=0,N-1,
    for(v=0,N-1,
      if(GetCoef(done,[u,v]),next);
      if(gcd([u,v,N])==1,
        listput(A,[u,v]);
				for(h=1,nH,done[ZNnorm(Hlist[h]*u,N),ZNnorm(Hlist[h]*v,N)]=1);
      )
    )
  );
  Vec(A);
}

GammaHCusps(N,Hlist)= \\ Reps of all cusps of GammaH
{
	my(Cusps=List(),done=matrix(N,N),nH=#Hlist);
	  for(u=0,N-1,
    for(v=0,N-1,
      if(GetCoef(done,[u,v]),next);
      if(gcd([u,v,N])==1,
        listput(Cusps,[u,v]);
				for(i=1,nH,
					h = Mod(Hlist[i],N);
					h1 = 1/h;
					for(x=1,N,
						done[ZNnorm(h1*(u+x*v),N),ZNnorm(h*v,N)]=1
					)
				)
      )
    )
  );
	Vec(Cusps);
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

/*{
ToCusp(s,N)= \\ Finds g=[p',*;q';*] in SL2(Z) with p~p', q~q' mod N, i.e. g.oo = (p,q)
  my(p,q,u,v,g,h);
  [p,q]=s;
  [u,v,g]=gcdext(p,q);
  h=lift(1/Mod(g,N));
  SL2lift([p,-v*h;q,u*h]);
}*/

{
BotToSL2(s,N)= \\ Finds [*,*;c';d'] in SL2(Z) with c~c', d~d' mod N
  my(c,d,u,v,g,h);
  [c,d]=s;
  [u,v,g]=gcdext(c,d);
  h=lift(1/Mod(g,N));
  SL2lift([c,-v*h;d,u*h]);
}

/*ActOnCuspi(g,i)=\\ g in SL2 -> the j s.t. g*cusp#i = cusp#j
{
my(c=Cusps[i]);
c=c*g;
c=liftint(Mod(c,N));
select(x->x==c,Cusps,1)[1];
}*/
