read("LMod.gp");
read("Etors.gp");
\\read("Dimensions.gp");
read("GammaH.gp");
\\read("Perms.gp");
\\read("qexp.gp");

/*l2(EN,P,Q,Tpe)=  \\ Slope of line (PQ)
{
	my([xP,yP]=GetCoef(EN,P),[xQ,yQ]=GetCoef(EN,Q),[T,pe,p,e]=Tpe);
	ZpXQ_div(liftall(yQ-yP),liftall(xQ-xP),T,pe,p,e);
}
\\ TODO methode Kamal addchain
l1(EN,P,Q,Tpe)=my([T,pe,p,e]=Tpe);Mod(Mod(sum(n=0,#EN-1,l2(EN,P,Q+n*P,Tpe)),T),pe); \\ sum_n l2(P,Q+nP) ~ sum_n l(P) +l(Q+nP)-l(Q+nP) ~ l*l(P)
*/
DivAdd1(A,B,dimres,p,excess,flag)=
{ \\ Mult RR spaces A and B.
	\\ Takes dimres + excess products A[u]*B[v] with random u,v
	\\ If flag, also return the list of pairs [u,v]
	my(nA=#A,nB=#B,m,n,C,u,v,uv,a,b,r);
	m = #A[,1];
	n = dimres+excess;
	C = matrix(m,n);
	uv = vector(n);
	while(1,
		for(j=1,n,
			u = 1+random(nA);
			v = 1+random(nB);
			a = A[,u];
			b = B[,v];
			uv[j] = [u,v];
			for(i=1,m,
				C[i,j] = a[i]*b[i]
			)
		);
		r = matindexrank(Mod(C,p))[2];
		if(#r == dimres,
			C = vecextract(C,r);
			uv = vecextract(uv,r);
			return(if(flag,[C,uv],C))
		);
		print1("@",#r,"/",dimres," ");
	);
}

BalancedDiv(d,degs)=
{
	my(n=#degs,s=vecsum(degs),D,q);
	q = d\s;
	d -= q*s;
	D = vector(n,i,q);
	for(i=1,n,
		if(d>=degs[i],
			d -= degs[i];
			D[i] +=1
		)
	);
	D;
}

DivPerturb(D,degs)=
{
	my(n=#degs,D2,i);
	d = sum(j=1,n,D[j]*degs[j]);
	D2 = BalancedDiv(d-1,degs);
	i=n;
	while(degs[i]==1,
		if(D2[i]+1!=D[i],
			D2[i]+=1;
			return(D2)
		);
		i-=1;
	);
	error("I don't know how to perturb this divisor");
}

Divo2Div(Do,Orbs,tags,n)=
{
	my(D,nO,o,no);
	nO = #Orbs;
	D = vector(n);
	for(i=1,nO,
		o = Orbs[i];
		no = #o;
		for(j=1,no,
			D[GetCoef(tags,o[j])] = Do[i];
		)
	);
	D;
}

MRRsubspace(M4qexps,D,T,p,e)=
{
	my(K,i=1,nD=vecsum(D),ncusps=#D);
	K = matrix(nD,#M4qexps[1]);
	for(s=1,ncusps,
		for(j=1,D[s],
			K[i,]=liftall(M4qexps[s][j,]);
			i++;
		)
	);
	matkerpadic_safe(K,T,p,e);
}

ModJacInit(N,H,p,a,e)=
{ \\ J_H(N) over Zq/p^e, q=p^a
	my(Hlist,Hlist1,Lp,g,Cusps,nCusps,CuspsGl,CuspsGalDegs,CuspTags,E,P0,Q0,zN,zNpows,MFrobE,tMFrobE,T,pe=p^e,Tpe,d,d1,Pts,nPts,PtTags,MPts,M2,M2gens,v,w,M,qprec,M2qexps,B,d0,M4,M4gens,M4qexps,E1,E2,U1,U2,M6,KV,f2,W0,J,CuspsQ);
	\\ Get H and H/+-1
  [Hlist,Hlist1] = GetHlist(N,H);
	if(Mod(6*N*#Hlist,p)==0,error("Bad p"));
	Lp = LMod(N,H,p);
	g = poldegree(Lp)/2;
	print("Genus ",g);
	\\ Get a curve E and a basis of E[N]
	[E,P0,Q0,zN,MFrobE] = EBasis(N,p,a,e);
	zNpows = vector(N);
	zNpows[1] = liftall(zN);
	for(i=2,N,
		zNpows[i] = liftall(zNpows[i-1]*zN)
	);
	tMFrobE = mattranspose(MFrobE);
	T = zN.mod;
	Tpe = [T,pe,p,e];
	print("E:y²=x³+",E.a4,"x+",E.a6," ",p,"-adically");
	\\print("Order ", N,": ",mordroot(x^2-ellap(E,p)*x+p,N));
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
	\\EN = liftall(EN);
	\\ Matrix of l(P) for P in E[N]
	print("Ml1");
	Ml1=matrix(N,N);
	for(x=1,N-1,Ml1[x,N]=l1(EN,[x,0],[0,1],T,pe,p,e)); \\ P=alpha(x,0) -> Q=alpha(0,1)
	for(x=1,N,for(y=1,N-1,Ml1[x,y]=l1(EN,[x,y],[1,0],T,pe,p,e))); \\ P=alpha(x,y), y!=0 -> Q=alpha(1,0)
	Ml1 = Mod(Mod(Ml1,T),pe);
	print("M2(GammaH)");
	\\ Find a basis for M2(GammaH(N))
	[Cusps,CuspsGal,CuspTags] = GammaHCusps(N,Hlist);
	nCusps = #Cusps;
	print(nCusps," cusps");
	CuspsGalDegs = apply(o->#o,CuspsGal);
	print("Degrees of Galois orbits: ",CuspsGalDegs);
	print("Action of Frob on Cusps:");
	d = g+nCusps-1; \\ dim M2(GammaH(N))
	[Pts,PtTags] = ANH(N,Hlist); \\ List of vectors (c,d) mod N,H
	nPts = #Pts;
	print(nPts," points in the fibre of X_H(N) -> X(1)");
  PtsFrob = Vecsmall(apply(P->GetCoef(PtTags,liftint(P*tMFrobE)),Pts)); \\ Frob([c,d]) = [c,d]*(t^MFrobE)
  print("Action of Frob on fibre: ",PtsFrob);
	MPts = apply(s->BotToSL2(s,N),Pts); \\ Matrices having these bottom rows
	\\ P_g = P_g' on X_H(N) <=> g,g' have same bottom row mod H
	d1=min(ceil(1.2*d),#Pts); \\ # gens
	M2=matrix(#Pts,d1);
	M2gens=vector(d1);
	TH = GammaHmodN(N,Hlist1); \\ elts of SL2(Z) representing GammaH mod N,+-1
	while(1,
		print("Attempt");
		/* qprec = 100*N;
		M2q = matrix(qprec,d1); DEBUG */
		\\ Take d1 forms in M2(Gamma(N))
  	for(j=1,d1, \\ of the form E_1^v * E_1^w : j = index of gen
			print1("Prod ");
			v=Pts[1+random(#Pts)];
    	w=Pts[1+random(#Pts)];
    	M2gens[j]=[v,w];
			\\ symmetrise them by sum slashing the transversal
    	for(P=1,nPts,
				\\ sum_x f_v f_w | T at P for T in TH
      	M2[P,j]=sum(i=1,#TH,GetCoef(Ml1,v*TH[i]*MPts[P])*GetCoef(Ml1,w*TH[i]*MPts[P]))
    	);
			/* DEBUG
			fvw = sum(i=1,#TH,E1qexp(v*TH[i],N,zN,qprec,Tpe,'x)*E1qexp(w*TH[i],N,zN,qprec,Tpe,'x));
			for(i=1,qprec,M2q[i,j]=polcoef(fvw,i-1)); */
 	  );
    /* DEBUG */
		printf("Checking M2 qexps");
		KM2 = Mod(Mod(matkerpadic(liftall(M2),T,p,e),T),pe);
		qprec = 20; \\ TODO adjust
		for(s=1,nCusps,
    	M2qexps = matrix(qprec,d1);
    	[M,w] = GammaHCuspData(Cusps[s],N,Hlist);
    	for(j=1,d1,
      	M2qexps[,j] = Mod(Mod(TrE2qexp(M2gens[j],N,TH,M,w,zNpows,qprec,T,pe,p,e)~,T),pe)
    	);
    	if(M2qexps*KM2!=0,error("qexp BAD!!!!"),print("Cusp ",s," OK"))
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
	M2gens = vecextract(M2gens,B);
	print("M2 qexps");
	qprec = 10; \\ TODO adjust
	M2qexps = vector(nCusps);
	for(s=1,nCusps,
		print1(s," ");
		M2qexps[s] = matrix(qprec,d);
		[M,w] = GammaHCuspData(Cusps[s],N,Hlist);
		for(j=1,d,
			M2qexps[s][,j] = Mod(Mod(TrE2qexp(M2gens[j],N,TH,M,w,zNpows,qprec,T,pe,p,e)~,T),pe)
		)
	);
	/*\\ Pull down to merom fns by div, TODO dangerous
	f2 = M2[,1]; \\ TODO ensure def/Q TODO necessary?
  for(i=1,#f2,
		u = liftall(f2[i]);
		u = ZpXQ_inv(u,T,p,e);
    M2[i,] *= u;
  );*/
	\\[M2,Pts,M2basis];
	print("M4(GammaH)");
	d0 = d+g-1;
	d = 2*d0+1-g;
	[M4,M4gens] = DivAdd1(M2,M2,2*d0+1-g,p,d,1);
	print("M4 qexps");
	M4qexps = vector(nCusps);
	for(s=1,nCusps,
		M4qexps[s] = matrix(qprec,d);
		for(j=1,d,
			for(n=0,qprec-1,
				M4qexps[s][n+1,j] = sum(k=0,n,M2qexps[s][k+1,M4gens[j][1]]*M2qexps[s][n-k+1,M4gens[j][2]])
			)
		)
	);
	print("Eval data");
	E1o = BalancedDiv(d0-g,CuspsGalDegs);
	E2o = DivPerturb(E1o,CuspsGalDegs);
	E1 = Divo2Div(E1o,CuspsGal,CuspTags,nCusps);
	E2 = Divo2Div(E2o,CuspsGal,CuspTags,nCusps);
	U1 = MRRsubspace(M4qexps,E1,T,p,e);
	U1 = liftall(M4*U1);
	U2 = MRRsubspace(M4qexps,E2,T,p,e);
	U2 = liftall(M4*U2);
	\\breakpoint();
	print("M6(GammaH)");
	M6 = DivAdd1(M4,M2,3*d0+1-g,p,d,0);
	print("Eqn mats");
	/*f2 = M2[,1]; \\ TODO ensure def/Q TODO necessary?
	for(i=1,#f2, \\ Push to weight 6
		M2[i,] *= f2[i]^2;
		M4[i,] *= f2[i]
	);*/
	V = apply(liftall,[M2,M4,M6]);
	KV = apply(x->matkerzq(x~,T,p,e)~,V); \\ TODO parallel
	\\ W0 = f*M2 c M4, f in M2
	\\ TODO need f def / Q ?
	f2 = M2[,1];
	W0 = M2;
	for(i=1,#f2,W0[i,] *= f2[i]);
	W0 = liftall(W0);
	\\J = [f,g,d0,L,T,p,e,pe,FrobMat,[V],[KV],W0,EvData,Z,FrobCyc];
	FrobMat = ZpXQ_FrobMat(T,p,e,pe);
	J=[0,g,d0,[],T,p,e,pe,FrobMat,V,KV,W0,[[U1],[U2]],[],PtsFrob];
	\\CuspsQ = [GetCoef(CuspTags,o[1]) | o<-CuspsGal,#o==1]; \\ Cusps def / Q
	CuspsQ = select(s->gcd(s[1],N)==1,Cusps,1); \\ DEBUG test N=16
	[J,vecextract(M4qexps,CuspsQ),vecextract(Cusps,CuspsQ)];
	\\[J,M4qexps,Cusps]; \\ DEBUG
}

PicEval(J,W)=
{
	my(Z);
	Z = PicEval0(J,W)[1,1];
	Z = concat(apply(M->liftall(M*Z),M4Q)); \\ TODO pass M4Q
	matrix(1,1,i,j,Z);
}

PicEvalDbg(J,W)=
{
	my(p=Jgetp(J),pe=Jgetpe(J),T=JgetT(J),e=Jgete(J),FrobMat=JgetFrobMat(J),Wp,Z,Zp,Y,Yp);
	\\J = PicRed(J,1);
	print("Checking PicEval commutes with Frob, precision ",O(p^e));
	Wp = PicFrob(J,W);
	Z = PicEval0(J,W)[1,1];
	Zp = PicEval0(J,Wp)[1,1];
	for(s=1,#M4Q,
		Y = liftall(M4Q[s]*Z);
		Yp = liftall(M4Q[s]*Zp);
		Y = apply(y->ZpXQ_div(y,Y[#Y],T,pe,p,e),Y);
		Yp = apply(y->ZpXQ_div(y,Yp[#Yp],T,pe,p,e),Yp);
		if(apply(y->Frob(y,FrobMat,T,pe),Y)!=Yp,print("Cusp ",s," BAD!!!"),print("Cusp ",s," OK"))
	);
}
