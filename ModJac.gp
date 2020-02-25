read("LMod.gp");
read("Etors.gp");
\\read("Dimensions.gp");
read("GammaH.gp");
read("Perms.gp");
\\read("qexp.gp");
read("Div.gp");

ModJacInit(N,H,p,a,e)=
{ \\ J_H(N) over Zq/p^e, q=p^a
	my(Hlist,Hlist1,TH);
	my(Lp,g,Cusps,nCusps,CuspsGal,CuspsGalDegs,CuspTags);
	my(E,P0,Q0,zN,zNpows,MFrobE,tMFrobE,T,pe=p^e,EN,todo,done,x,y,Ml1);
	my(d,d1,Pts,nPts,PtTags,MPts,M2,M2gens,v,w,M,qprec,M2qexps,B,d0,C0,U0,V1,V2,V3,V2gens,V1qexps,V2qexps,E1,E2,U1,U2,KV,f2,W0,J,CuspsQ);
	my(IU,MU);
	\\ Get H and H/+-1
  [Hlist,Hlist1] = GetHlist(N,H);
	if(Mod(6*N*#Hlist,p)==0,error("Bad p"));
	Lp = LMod(N,H,p);
	g = poldegree(Lp)/2;
	print("Genus ",g);
	\\ Get a curve E and a basis of E[N]
	[E,P0,Q0,zN,MFrobE] = EBasis(N,p,a,e);
	for(i=1,a,
		if(MFrobE^i==1,
			print("MFrobE has order ",i);
			break
		)
	);
	zNpows = vector(N);
	zNpows[1] = liftall(zN);
	for(i=2,N,
		zNpows[i] = liftall(zNpows[i-1]*zN)
	);
	tMFrobE = mattranspose(MFrobE);
	T = zN.mod;
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
	todo = List();
	for(x=1,N-1,for(y=1,N-1,listput(todo,[x,y])));
	done = parapply(X->elladd_padic(E.a4,EN[X[1],N],EN[N,X[2]],T,pe,p,e),Vec(todo));
	for(i=1,#todo,[x,y]=todo[i];EN[x,y]=done[i]);
	\\ Matrix of l(P) for P in E[N]
	print("Ml1");
	\\ TODO parallelise
	Ml1=matrix(N,N);
	\\ TODO do we use non coprime entries?
	todo = List();
	for(x=1,N-1,listput(todo,[[x,N],[0,1]])); \\ P=alpha(x,0) -> Q=alpha(0,1)
	for(x=1,N,for(y=1,N-1,listput(todo,[[x,y],[1,0]]))); \\ P=alpha(x,y), y!=0 -> Q=alpha(1,0)
	done = parapply(X->l1(EN,X[1],X[2],T,pe,p,e),Vec(todo));
	for(i=1,#todo,[x,y]=todo[i][1];Ml1[x,y]=done[i]); 
	Ml1 = Mod(Mod(Ml1,T),pe);
	print("M2(GammaH)");
	\\ Find a basis for M2(GammaH(N))
	[Cusps,CuspsGal,CuspsQ,CuspsMats,CuspsWidths,CuspTags] = GammaHCusps(N,Hlist);
	nCusps = #Cusps;
	print(nCusps," cusps");
	CuspsGalDegs = apply(o->#o,CuspsGal);
	print("Degrees of Galois orbits: ",CuspsGalDegs);
	d = g+nCusps-1; \\ dim M2(GammaH(N))
	[Pts,PtTags] = ANH(N,Hlist); \\ List of vectors (c,d) mod N,H
	nPts = #Pts;
	print(nPts," points in the fibre of X_H(N) -> X(1)");
  PtsFrob = Vecsmall(apply(P->GetCoef(PtTags,liftint(P*tMFrobE)),Pts)); \\ Frob([c,d]) = [c,d]*(t^MFrobE)
  print("Action of Frob on fibre: ",PtsFrob);
	MPts = apply(s->BotToSL2(s,N),Pts); \\ Matrices having these bottom rows
	\\ P_g = P_g' on X_H(N) <=> g,g' have same bottom row mod H
	d1=min(ceil(1.2*d),#Pts); \\ # gens
	TH = GammaHmodN(N,Hlist1); \\ elts of SL2(Z) representing GammaH mod N,+-1
	while(1,
		print("Attempt");
		\\ Take d1 forms in M2(Gamma(N))
		M2gens=vector(d1,j,[Pts[1+random(#Pts)],Pts[1+random(#Pts)]]);
		\\ of the form E_1^v * E_1^w
		M2=matconcat(
			parapply(
				vw->apply(MP->
					sum(i=1,#TH,
						GetCoef(Ml1,vw[1]*TH[i]*MP)*GetCoef(Ml1,vw[2]*TH[i]*MP)
					)
				,MPts)~
			,M2gens)
		);
    /* DEBUG */
		/*printf("Checking M2 qexps");
		KM2 = Mod(Mod(matkerpadic(liftall(M2),T,p,e),T),pe);
		qprec = 20; \\ TODO adjust
		for(s=1,nCusps,
    	M2qexps = matrix(qprec,d1);
    	[M,w] = GammaHCuspData(Cusps[s],N,Hlist);
    	for(j=1,d1,
      	M2qexps[,j] = Mod(Mod(TrE2qexp(M2gens[j],N,TH,M,w,zNpows,qprec,T,pe,p,e)~,T),pe)
    	);
    	if(M2qexps*KM2!=0,error("qexp BAD!!!!"),print("Cusp ",s," OK"))
		);*/
  	\\ See if we span all of M2(GammaH) by checking full rank
		print("linalg");
  	B=matindexrank(Mod(liftint(M2),p))[2]; \\ working mod p for efficiency
  	if(#B>d,error("Bug M2(GammaH)")); \\ Not supposed to happen
  	if(#B==d,break); \\ This is what we want
  	print("Retrying: the products of Eis series of wt 1 span a subspace of dim ",#B," out of ",d)
	);
	\\ Extract basis
	M2 = vecextract(M2,B);
	M2gens = vecextract(M2gens,B);
	print("M2 qexps");
	qprec = 8; \\ TODO adjust
	M2qexps = vector(nCusps);
	for(s=1,nCusps,
		print1(s," ");
		M = CuspsMats[s];
		w = CuspsWidths[s];
		M2qexps[s] = matconcat(
			parapply(
				vw->Mod(Mod(
					TrE2qexp(vw,N,TH,M,w,zNpows,qprec,T,pe,p,e)~
				,T),pe)
			,M2gens)
		)
	);
	print("Pruning");
	\\ Prune: M2 -> S2(3 cusps) = M2(-C0)
	d0 = 2*g+1;
	C0o = BalancedDiv(nCusps-3,CuspsGalDegs);
	C0 = Divo2Div(C0o,CuspsGal,CuspTags,nCusps);
	U0 = MRRsubspace(M2qexps,C0,T,p,e);
	V1qexps = parapply(page->page*U0,M2qexps);
	\\ Prune more: no need to evaluate at that many points
	[PtsFrob,Pts] = SubPerm(PtsFrob,5*d0+1);
	V1 = vecextract(M2,Pts,vector(#M2,i,i))*U0;
	print("M4(GammaH)");
	d = 2*d0+1-g;
	[V2,V2gens] = DivAdd1(V1,V1,2*d0+1-g,p,d0,1);
	print("M4 qexps");
	V2qexps = parapply(
		page->matrix(qprec,d,n,j,
			sum(k=0,n-1,page[k+1,V2gens[j][1]]*page[n-k,V2gens[j][2]])
		)
	,V1qexps);
	print("Eval data");
	E1o = BalancedDiv(d0-g,CuspsGalDegs); \\ d0-g = g+1
	E2o = DivPerturb(E1o,CuspsGalDegs);
	E1 = Divo2Div(E1o,CuspsGal,CuspTags,nCusps);
	E2 = Divo2Div(E2o,CuspsGal,CuspTags,nCusps);
	export(MRRsubspace);
	[U1,U2] = parapply(Ei->
		liftall(V2*MRRsubspace(V2qexps,2*C0+Ei,T,p,e)), \\ TODO: eqns for 2*C0 already satified. Make MRRsubspace more flexible.
	[E1,E2]);
	print("M6(GammaH)");
	V3 = DivAdd1(V2,V1,3*d0+1-g,p,d0,0);
	print("Eqn mats");
	V = apply(liftall,[V1,V2,V3]);
	IU = matindexrank(Mod(V2,p))[1];
	MU = vecextract(V2,IU,"..");
	MU = ZpXQM_inv(liftall(MU),T,p,e);
	MU = matconcat(vecextract(V2qexps,CuspsQ)~)*MU;
	MU = liftall(MU);
	KV = parapply(x->mateqnpadic(x,T,p,e),V);
	\\ W0 = f*M2 c M4, f in M2
	f2 = V1[,1];
	W0 = V1;
	for(i=1,#f2,W0[i,] *= f2[i]);
	W0 = liftall(W0);
	FrobMat = ZpXQ_FrobMat(T,p,e,pe);
	J=[0,g,d0,[],T,p,e,pe,FrobMat,V,KV,W0,[[U1],[U2],IU,MU],[],PtsFrob];
}

/*PicEval(J,W)=
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
}*/
