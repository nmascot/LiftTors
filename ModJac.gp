read("LMod.gp");
read("Etors.gp");
\\read("Dimensions.gp");
read("GammaH.gp");
read("Perms.gp");

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

DivAdd1(A,B,dimres,p,excess)=
{
	my(m,n,C,a,b,r);
	m = #A[,1];
	n = dimres+excess;
	C = matrix(m,n);
	while(1,
		for(j=1,n,
			a = A[,1+random(#A)];
			b = B[,1+random(#B)];
			for(i=1,m,
				C[i,j] = a[i]*b[i]
			)
		);
		r = matindexrank(Mod(C,p))[2];
		if(#r == dimres,return(vecextract(C,r)));
		print1("@",#r,"/",dimres);
	);
}
		

ModJacInit(N,H,p,a,e)=
{ \\ J_H(N) over Zq/p^e, q=p^a
	my(Hlist,Hlist1,Lp,g,Cusps,nCusps,CuspTags,E,P0,Q0,zN,MFrobE,tMFrobE,T,pe=p^e,d,d1,Pts,nPts,PtTags,MPts,M2,M2gens,v,w,B,d0,M4,M6,KV,KV3);
	\\ Get H and H/+-1
  [Hlist,Hlist1] = GetHlist(N,H);
	if(Mod(6*N*#Hlist,p)==0,error("Bad p"));
	Lp = LMod(N,H,p);
	g = poldegree(Lp)/2;
	print("Genus ",g);
	\\ Get a curve E and a basis of E[N]
	[E,P0,Q0,zN,MFrobE] = EBasis(N,p,a,e);
	tMFrobE = mattranspose(MFrobE);
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
	print("M2(GammaH)");
	\\ Find a basis for M2(GammaH(N))
	[Cusps,CuspTags] = GammaHCusps(N,Hlist);
	nCusps = #Cusps;
	print(nCusps," cusps");
	print("Action of Frob on Cusps:");
	d = g+nCusps-1; \\ dim M2(GammaH(N))
	[Pts,PtTags] = ANH(N,Hlist); \\ List of vectors (c,d) mod N,H
	nPts = #Pts;
	print(nPts," points in the fibre of X_H(N) -> X(1)");
  PtsFrob = Vecsmall(apply(P->GetCoef(PtTags,P*tMFrobE),Pts)); \\ Frob([c,d]) = [c,d]*(t^MFrobE)
  print("Action of Frob on fibre:");
	print(PtsFrob); \\ TODO use this
	MPts = apply(s->BotToSL2(s,N),Pts); \\ Matrices having these bottom rows
	\\ P_g = P_g' on X_H(N) <=> g,g' have same bottom row mod H
	d1=min(ceil(1.2*d),#Pts); \\ # gens
	M2=matrix(#Pts,d1);
	M2gens=vector(d1);
	TH = GammaHmodN(N,Hlist1); \\ elts of SL2(Z) representing GammaH mod N,+-1
	while(1,
		print("Attempt");
		\\ Take d1 forms in M2(Gamma(N)
  	for(j=1,d1, \\ of the form E_1^v * E_1^w : j = index of gen
			print("Prod");
			v=Pts[1+random(#Pts)];
    	w=Pts[1+random(#Pts)];
    	M2gens[j]=[v,w];
			\\ symmetrise them by sum slashing the transversal
    	for(P=1,nPts,
				\\ sum_x f_v f_w | T at P for T in TH
      	M2[P,j]=sum(i=1,#TH,GetCoef(Ml1,v*TH[i]*MPts[P])*GetCoef(Ml1,w*TH[i]*MPts[P]))
    	)
 	  );
  	\\ See if we span all of M2(GammaH) by checking full rank
		print("linalg");
  	B=matindexrank(Mod(M2,p))[2]; \\ working mod p for efficiency
  	if(#B>d,error("Bug M2(GammaH)")); \\ Not supposed to happen
  	if(#B==d,break); \\ This is what we want
  	print("Retrying: the products of Eis series of wt 1 span a subspace of dim ",#B," out of ",d)
	);
	\\ Extract basis
	M2 = vecextract(M2,B);
	M2basis = vecextract(M2gens,B);
	\\ Pull down to merom fns by div, TODO dangerous
	f2 = M2[,1]; \\ TODO ensure def/Q TODO necessary?
  for(i=1,#f2,
		u = liftall(f2[i]);
		u = ZpXQ_inv(u,T,p,e);
    M2[i,] *= u;
  );
	\\[M2,Pts,M2basis];
	print("M4(GammaH)");
	d0 = d+g-1;
	M4 = DivAdd1(M2,M2,2*d0+1-g,p,d);
	print("M6(GammaH)");
	M6 = DivAdd1(M4,M2,3*d0+1-g,p,d);
	print("Eqn mats");
	/*f2 = M2[,1]; \\ TODO ensure def/Q TODO necessary?
	for(i=1,#f2, \\ Push to weight 6
		M2[i,] *= f2[i]^2;
		M4[i,] *= f2[i]
	);*/
	M2 = liftall(M2);
	M4 = liftall(M4);
	M6 = liftall(M6);
	KV = matkerzq(M4~,T,p,e)~;
	KV3 = matkerzq(M6~,T,p,e)~;
	\\J = [f,g,d0,L,T,p,e,pe,Frob,V,KV,W0,Z,FrobCyc,V3,KV3,EvData];
	FrobMat = ZpXQ_FrobMat(T,p,e,pe);
	[0,g,d0,[],T,p,e,pe,FrobMat,M4,KV,M2,[],PtsFrob,M6,KV3,[]];
	\\ TODO: EvData
}
