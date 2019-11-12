read("Etors.gp");
read("Dimensions.gp");

GetCoef(A,a)=my(i,j,N=#A);[i,j]=liftall(Mod(a,N));if(i==0,i=N);if(j==0,j=N);A[i,j];
ZNnorm(x,N)=my(y=x%N);if(y==0,N,y);
ZNneg(x,N)=my(y=lift(Mod(-x,N)));if(y==0,N,y);

SL2lift(M)= \\ Finds M' in SL2(Z), M=M' mod (|M|-1)
{
  my(U,V,D,a,b);
  [U,V,D]=matsnf(M,1); \\ U*M*V = D = diag(a,b), so |M|=ab
  a=D[1,1];
  b=D[2,2];
  U^-1*[1,-1;1-b,b]*[1,0;1-a,1]*[1,b;0,1]*V^-1;
} \\ [1,-1;1-b,b]*[1,0;1-a,1]*[1,b;0,1] = [a, ab-1; 1-ab, b+(1-ab)b]

XNCusps(N)= \\ List of cusps of X(N)
{	\\ Find all (u,v) s.t. gcd(u,v,N)=1 / +-1
	my(Cusps=List(),done=matrix(N,N));
	for(u=0,N-1,
  	for(v=0,N-1,
    	if(GetCoef(done,[u,v]),next);
    	if(gcd([u,v,N])==1,
      	listput(Cusps,[u,v]);
      	done[ZNnorm(u,N),ZNnorm(v,N)]=1;
      	done[ZNnorm(-u,N),ZNnorm(-v,N)]=1;
    	)
  	)
	);
	Vec(Cusps);
}

{
ToCusp(s,N)= \\ Finds [p',*;q';*] in SL2(Z) with p~p', q~q' mod N
  my(p,q,u,v,g,h);
  [p,q]=s;
  [u,v,g]=gcdext(p,q);
  h=lift(1/Mod(g,N));
  SL2lift([p,-v*h;q,u*h]);
}

{
Cusp2X1(s,N)= \\ Finds [*,*;c';d'] in SL2(Z) with c~c', d~d' mod N
  my(c,d,u,v,g,h);
  [c,d]=s;
  [u,v,g]=gcdext(c,d);
  h=lift(1/Mod(g,N));
  SL2lift([c,-v*h;d,u*h]);
}

ActOnCuspi(g,i)=\\ g in SL2 -> the j s.t. g*cusp#i = cusp#j
{
my(c=Cusps[i]);
c=c*g;
c=liftint(Mod(c,N));
select(x->x==c,Cusps,1)[1];
}

E1atCusp(s,N,z)= \\ Constant coef of Eis series E_1^s of level N
{ \\ epxressed in terms of prim root z
  my([c,d]=s);
  c=c%N;
  d=d%N;
  if(c==0,
    1/2*(z^d+1)/(z^d-1),
    1/2-c/N
  );
}

l2(EN,P,Q,T,pe,p,e)=  \\ Slope of line (PQ)
{
	my([xP,yP]=GetCoef(EN,P),[xQ,yQ]=GetCoef(EN,Q));
	ZpXQ_div(liftall(yQ-yP),liftall(xQ-xP),T,pe,p,e);
}
\\ TODO methode Kamal addchain
l1(EN,P,Q,T,pe,p,e)=Mod(Mod(sum(n=0,#EN-1,l2(EN,P,Q+n*P,T,pe,p,e)),T),pe); \\ sum_n l2(P,Q+nP) ~ sum_n l(P) +l(Q+nP)-l(Q+nP) ~ l*l(P)

ModJacInit(N,H,p,a,e)=
{ \\ J_H(N) over Zq/p^e, q=p^a
	my(E,P0,Q0,zN,T,pe=p^e,d,d1,d2,Cusps,M21,M21gens,M21basis,v,w,MP,B);
	\\ Get a curve E and a basis of E[N]
	[E,P0,Q0,zN] = EBasis(N,p,a,e);
	T = zN.mod;
	\\ Write down all N-torsion: : this is a naive level structure alpha: (Z/NZ)² ~ E[N]
	EN=matrix(N,N); \\ [[ m P0 + n Q0 ]]
  EN[1,N]=P0;
  EN[N,1]=Q0;
  for(x=2,N-1,
    EN[x,N]=elladd_padic(E.a4,EN[x-1,N],P0,T,pe,p,e);
    EN[N,x]=elladd_padic(E.a4,EN[N,x-1],Q0,T,pe,p,e);
  );
  for(x=1,N-1,
    for(y=1,N-1,
      EN[x,y]=elladd_padic(E.a4,EN[x,N],EN[N,y],T,pe,p,e)
    )
  );
	\\ Matrix of l(P) for P in E[N]
	print("Ml1");
	Ml1=matrix(N,N);
	for(x=1,N-1,Ml1[x,N]=l1(EN,[x,0],[0,1],T,pe,p,e)); \\ P=alpha(x,0) -> Q=alpha(0,1)
	for(x=1,N,for(y=1,N-1,Ml1[x,y]=l1(EN,[x,y],[1,0],T,pe,p,e))); \\ P=alpha(x,y), y!=0 -> Q=alpha(1,0)
	print("Gamma1");
	\\ Find a basis for M2(Gamma1(N))
	\\ TODO add variables to my()
	\\ TODO GammaH(N)
	d=DimMk(2,1,N); \\ dim M2(Gamma1(N))
	Pts = XNCusps(N); \\ List of vectors (c,d) mod N,+-1
	\\ P_g = P_g' on X_1(N) <=> g,g' have same bottom row
	d1=min(floor(1.2*d),#Pts); \\ TODO # pts at which we observe lin indep
	d2=d1; \\ # gens
	M21=matrix(#Pts,d2);
	M21gens=vector(d2);
	while(1,
		print("Attempt");
		\\ Take d2 forms in M2(Gamma(N)
  	for(j=1,d2, \\ of the form E_1^v * E_1^w : j = index of gen
			print("Prod");
			v=Pts[1+random(#Pts)];
    	w=Pts[1+random(#Pts)];
    	M21gens[j]=[v,w];
			\\ symmetrise them by sum slashing [1,x;0,1]
    	for(P=1,#Pts,
      	MP=Cusp2X1(Pts[P],N); \\ matrix in SL2(Z) having bottow row P
				\\ sum_x f_v f_w | [1,x;0,1] at P
      	M21[P,j]=sum(x=1,N,GetCoef(Ml1,v*[1,x;0,1]*MP)*GetCoef(Ml1,w*[1,x;0,1]*MP))
    	)
 	  );
  	\\ See if we span all of M2(Gamma1) by checking values at d1 cusps TODO
		print("linalg");
  	B=matindexrank(Mod(M21,p))[2]; \\ working mod p for efficiency
  	if(#B>d,error("Bug M2(Gamma1)")); \\ Not supposed to happen
  	if(#B==d,break); \\ This is what we want
  	print("Retrying: the products of Eis series of wt 1 span a subspace of dim ",#B," out of ",d)
	);
	\\ Extract basis
	M21 = vecextract(M21,B);
	M21basis = vecextract(M21gens,B); 
	[M21,M21basis,Pts];
}
