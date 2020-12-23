read("LMod.gp");
read("Etors.gp");
read("GammaH.gp");
read("Perms.gp");
read("Div.gp");

ModJacInit(N,H,p,a,e,qprec,Lp,UseTp)=
{ \\ J_H(N) over Zq/p^e, q=p^a
	my(Hlist,Hlist1,TH);
	my(g,Cusps,nCusps,CuspsGal,CuspsGalDegs,CuspTags,CuspQexp,CuspsQexp_list,CuspsGalDiamp,CuspsGalDiampDegs);
	my(E,P0,Q0,zN,zNpows,MFrobE,tMFrobE,T,pe=p^e,EN,todo,done,x,y,Ml1);
	my(d,d1,Pts,nPts,PtTags,MPts,PtsFrob,PtsDiamp,Perms);
	my(E1o,E2o,E1,E2,C0o,C0,Emax);
	my(M2,M2gens,v,w,M,sprec,M2qexps,B,d0,U0,V1,V2,V3,V2gens,V1qexps,V2qexps,U1,U2,KV,f2,W0);
	my(M4Q,IU,MU,Auts);
	\\ Get H and H/+-1
  [Hlist,Hlist1] = GetHlist(N,H);
  if(Mod(6*N*#Hlist,p)==0,error("Bad p"));
	if(qprec==0,qprec=10);
	if(Lp==0,Lp = LMod(N,H,p));
	g = poldegree(Lp)/2;
  print("Genus ",g);
	\\ Get data about cusps
	[Cusps,CuspsGal,CuspsQexp_list,CuspsMats,CuspsWidths,CuspTags] = GammaHCusps(N,Hlist);
	nCusps = #Cusps;
  print(nCusps," cusps");
  CuspsGalDegs = apply(o->#o,CuspsGal);
  print("Degrees of Galois orbits of cusps: ",CuspsGalDegs);
  if(UseTp,
		CuspsGalDiamp = GammaHCusps_GalDiam_orbits(p,Cusps,CuspsGal,CuspTags);
		CuspsGalDiampDegs = apply(o->#o,CuspsGalDiamp);
  	print("Degrees of Galois,<p> orbits of cusps: ",CuspsGalDiampDegs)
	);
  CuspQexp=vector(nCusps); \\ Vec of bits: can we have rational q-exps at this cusp?
  for(i=1,#CuspsQexp_list,CuspQexp[CuspsQexp_list[i]]=1);
  \\ Get data about fibre
  [Pts,PtTags] = ANH(N,Hlist); \\ List of vectors (c,d) mod N,H. Represents fibre XH->X(1).
  nPts = #Pts;
  print(nPts," points on the fibre of X_H(",N,") -> X(1)");
	if(UseTp,
		PtsDiamp = Vecsmall(apply(v->GetCoef(PtTags,p*v),Pts)) \\ Perm of fibre induced by <p>
	);
	MPts = apply(s->BotToSL2(s,N),Pts); \\ Matrices having these bottom rows
	\\ P_g = P_g' on X_H(N) <=> g,g' have same bottom row mod H
	\\ Get a curve E and a basis of E[N]
	print("Finding elliptic curve E such that E[",N,"] is defined over GF(",p,"^",a,")");
	[E,P0,Q0,zN,MFrobE] = EBasis(N,p,a,e);
	T = zN.mod;
  print("E:y²=x³+",E.a4,"x+",E.a6," with accuracy O(",p,"^",e,") and residual degree a=",a);
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
	\\ Get action of Frob on fibre
	tMFrobE = mattranspose(MFrobE);
  PtsFrob = Vecsmall(apply(P->GetCoef(PtTags,liftint(P*tMFrobE)),Pts)); \\ Frob([c,d]) = [c,d]*(t^MFrobE)
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
	Ml1=matrix(N,N);
	\\ TODO do we use non coprime entries?
	todo = List();
	for(x=1,N-1,listput(todo,[[x,N],[0,1]])); \\ P=alpha(x,0) -> Q=alpha(0,1)
	for(x=1,N,for(y=1,N-1,listput(todo,[[x,y],[1,0]]))); \\ P=alpha(x,y), y!=0 -> Q=alpha(1,0)
	done = parapply(X->l1(EN,X[1],X[2],T,pe,p,e),Vec(todo));
	for(i=1,#todo,[x,y]=todo[i][1];Ml1[x,y]=done[i]); 
	\\ Find a basis for M2(GammaH(N))
	d = g+nCusps-1; \\ dim M2(GammaH(N))
	print("M2(GammaH) (dim ",d,")");
	d1=min(ceil(4*d/3),#Pts); \\ # gens
	TH = GammaHmodN(N,Hlist1); \\ elts of SL2(Z) representing GammaH mod N,+-1
	while(1,
		\\ Take d1 forms in M2(Gamma(N))
		M2gens=vector(d1,j,[Pts[1+random(#Pts)],Pts[1+random(#Pts)]]);
		\\ of the form E_1^v * E_1^w
		M2 = M2mat(M2gens,Ml1,TH,MPts,T,pe);
		M2 = Mod(M2,T);
  	\\ See if we span all of M2(GammaH) by checking full rank
		B=matindexrank(Mod(M2,p))[2]; \\ working mod p for efficiency
  	if(#B>d,error("Bug M2(GammaH)")); \\ Not supposed to happen
  	if(#B==d,break); \\ This is what we want
  	print("Retrying: the products of Eis series of wt 1 span a subspace of dim ",#B," out of ",d);
		d1 += d-#B;
	);
	\\ Extract basis
	M2 = Mod(vecextract(M2,B),pe);
	M2gens = vecextract(M2gens,B);
	\\ TODO
	\\ Prepare divisors to know min qprec
	\\ Prune: M2 -> S2(3 cusps) = M2(-C0)
	\\ TODO sort
	if(UseTp,
  	C0o = BalancedDivInf(nCusps-3,CuspsGalDiampDegs);
  	C0 = Divo2Div(C0o,CuspsGalDiamp,CuspTags,nCusps);
		print("Wanted deg C0=",nCusps-3,", got ",vecsum(C0))
	,
		C0o = BalancedDiv(nCusps-3,CuspsGalDegs);
    C0 = Divo2Div(C0o,CuspsGal,CuspTags,nCusps)
	);
	d0 = 2*g-2+nCusps-vecsum(C0);
	\\ Evaluation
	E1o = BalancedDiv(d0-g,CuspsGalDegs);
  E2o = DivPerturb(E1o,CuspsGalDegs);
  E1 = Divo2Div(E1o,CuspsGal,CuspTags,nCusps);
  E2 = Divo2Div(E2o,CuspsGal,CuspTags,nCusps);
	Emax = 2*C0+vector(nCusps,i,max(E1[i],E2[i]));
	M2qexps = vector(nCusps);
	print1("M2 qexps, cusp ");
	for(s=1,nCusps,
		print1(s," ");
		M = CuspsMats[s];
		w = CuspsWidths[s];
		sprec = if(CuspQexp[s],max(Emax[s],qprec),Emax[s]);
		M2qexps[s] = matconcat(
			parapply(
				vw->Mod(Mod(
					TrE2qexp(vw,N,TH,M,w,zNpows,sprec,T,pe,p,e)~
				,T),pe)
			,M2gens)
		)
	);
	print("\nPruning: dim ",d-vecsum(C0),"; eval on >=",5*d0+1," pts");
	\\ Prune: M2 -> S2(>=3 cusps) = M2(-C0)
	U0 = MRRsubspace(M2qexps,C0,0,T,pe,p,e);
	V1qexps = parapply(page->page*U0,M2qexps);
	\\ Prune more: no need to evaluate at that many points
	if(UseTp,
		[Pts,Perms] = SubPerm_multi([PtsFrob,PtsDiamp],5*d0+1);
		[PtsFrob,PtsDiamp] = Perms;
	,
		[Pts,Perms] = SubPerm_multi([PtsFrob],5*d0+1);
    [PtsFrob] = Perms
	);
	print("Wanted nZ=",5*d0+1,", got ",#Pts);
	V1 = vecextract(M2,Pts,vector(#M2,i,i))*U0;
	d = 2*d0+1-g;
	print("M4(GammaH)(-2C0) (dim ",d,")");
	[V2,V2gens] = DivAdd1(V1,V1,2*d0+1-g,p,d0,1);
	print("M4 qexps");
	V2qexps = parapply(
		page->matrix(matsize(page)[1],d,n,j,
			sum(k=0,n-1,page[k+1,V2gens[j][1]]*page[n-k,V2gens[j][2]])
		)
	,V1qexps);
	print("M6(GammaH)(-3C0) (dim ",3*d0+1-g,")");
  V3 = DivAdd1(V2,V1,3*d0+1-g,p,ceil(3*d0/2),0);
  V = apply(liftall,[V1,V2,V3]);
	print("Eval data");
	export(MRRsubspace);
	[U1,U2] = parapply(Ei->
		liftall(V2*MRRsubspace(V2qexps,Ei,2*C0,T,pe,p,e)),
	[E1,E2]);
	IU = matindexrank(Mod(liftint(V2),p))[1]; \\ rows of V2 forming invertible block
	MU = vecextract(V2,IU,".."); \\ this invertible block
	MU = ZpXQMinv(liftall(MU),T,pe,p,e); \\ inverse
	M4Q = apply(i->V2qexps[i][1+2*C0[i]..#V2qexps[i]~,],CuspsQexp_list)~; \\ q-exps of forms in V2 used for eval
	MU = liftall(matconcat(M4Q)*MU); \\ apply this change of basis
	print("Eqn mats");
	KV = parapply(x->mateqnpadic(x,T,pe,p,e),V);
	\\ W0 = f*M2 c M4, f in M2
	f2 = V1[,1];
	W0 = V1;
	for(i=1,#f2,W0[i,] *= f2[i]);
	W0 = liftall(W0);
	FrobMat = ZpXQ_FrobMat(T,p,e,pe);
	Auts=if(UseTp,[PtsDiamp],[]);
	[0,g,d0,[],T,p,e,pe,FrobMat,V,KV,W0,[[U1],[U2],IU,MU],[],PtsFrob,Auts];
}

PicTp(J,W)=
{ /* Action of Tp on ModJac (p = the prime s.t. we work p-adically) */
  my(a,W1,W2); /* TODO optimise */
	if(JgetAutsCyc(J)==[],
		error("This Jacobian does not contain the extra data required to compute the action of Tp. Pass UseTp=1 to ModJacInit().")
	);
  a = poldegree(JgetT(J));
  W1 = PicFrob(J,W);
  W2 = PicMul(J,W,Jgetp(J),2);
  for(i=1,a-1,W2=PicFrob(J,W2)); /* TODO implement FrobInv */
  W2 = PicAut(J,W2,1);
  PicAdd(J,W1,W2);
}
export(PicTp);
