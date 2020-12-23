read("GalRep.gp");
read("ModJac.gp");

mfyt(f,l,coeffs)=
{
	my(N,k,eps,P,Phi,z,tz,Z,T,Y,L);
	[N,k,eps,P,Phi]=mfparams(f);
	[z,tz,k]=rnfequation(Phi,liftpol(P),1);
	\\ tz root of Phi, y+k*tz root of z
	tz=liftall(tz);
	Z=polrootsmod(z,l); \\ Possible values of Z mod l
	T=apply(z->subst(tz,'y,z),Z); \\ Corresponding values of t
	Y=vector(#Z,i,Z[i]-k*T[i]); \\ Corresponding values of y
	L = [1..#Z]; \\ List of possibilities
	for(i=1,#coeffs,
		[n,a] = coeffs[i];
		ayt = liftall(mfcoef(f,n));
		L = select(j->subst(subst(ayt,'t,T[j]),'y,Y[j])==a,L);
	);
	if(#L==0,error("Inconsistent data"));
	if(#L>1,error("Insufficient data"));
	L=L[1];
	[Y[L],T[L]];
}


mfbestp(f,l,coeffs,pmax,UseTp)=
{
	my(N,k,eps,P,Phi,yl,tl,o,ZNX,pmin=5,lN,H,listp,qf,best,p,ap,epsp,chi,Lp,Psi,a1,a2,a,mulTp);
	UseTp=1; \\ TODO
	[N,k,eps,P,Phi]=mfparams(f);
	[yl,tl] = mfyt(f,l,coeffs);
	[ZNX,eps] = znchar(f);
	o = charorder(ZNX,eps);
	if(type(pmax)=="t_VEC",[pmin,pmax]=pmax;pmin=max(5,pmin));
	lN = if(k==2,N,l*N);
	H = select(h->gcd(lN,h)==1,[1..lN-1]);
	H = select(h->Mod(h^(k-2),l)*chareval(ZNX,eps,h,[tl,o])==1,H);
	listp = select(p->Mod(l*N,p),primes([pmin,pmax]));
	print("qexp");
	qf = mfcoefs(f,pmax);
	print("Lp");
	Lp = LMod_multi(lN,H,listp);
	if(Lp[1][1]==1,error("The representation is a power of the cyclotomic character"));
	best = [];
	for(i=1,#listp,
		p = listp[i];
		if(Mod(l*N,p)==0,next);
		ap = subst(subst(liftall(qf[p+1]),'t,tl),'y,yl);
		epsp = chareval(ZNX,eps,p,[tl,o]);
		chi = 'x^2-ap*'x+p^(k-1)*epsp;
		Psi = Mod(Lp[i][1],l)/chi;
		if(poldegree(denominator(Psi)),error("Bug in charpoly"));
		if(poldegree(gcd(chi,Psi)),
			print("p=",p, " has Frob multiplicity");
			if(UseTp,
				mulTp=vecsum(apply(M->#matker(Mod(M,l)-ap),Lp[i][2]));
				print("p=",p, " has Tp multiplicity=",mulTp);
				if(mulTp>1,next)
			,
			next)
		);
		a1 = mordroot(chi,l);
		a2 = znorder(Mod(p,lN));
		a = lcm(a1,a2);
		print("p="p,": Needs deg ",a," (",a1," to split rep, ",a2," for roots of 1), log #J=",round(log(polresultant(Lp[i][1],'x^a-1))),")");
		a = lcm(a1,a2);
		if(best==[] || a < best[2], best=[p,a,Lp[i][1],ap,chi]);
	);
	[H,best];
}

mfgalrep(f,l,coeffs,pmax,D,qprec,UseTp,threadlim)=
{
	my(N,k,H,best,p,a,Lp,ap,chi,chiOO,gcdchi,e,pe,J,CuspsQ,J1,B,matFrob,LinTests,R,TpB,TpR,Tp,KTp,i,M,WB,cWB,TI,Z,AF,def_threads,t0);
	if(threadlim,
		def_threads = default(nbthreads);
		if(#threadlim!=4,error("Give 4 thread limits: Jacobian initialisation, generation of points mod p, p-adic lift, and evaluation"))
	);
	if(threadlim,default(nbthreads,threadlim[1]));
	t0 = [getabstime(),getwalltime()];
	print("--> Finding prime p with small residual degree");
	[H,best] = mfbestp(f,l,coeffs,pmax,UseTp);
	[p,a,Lp,ap,chi] = best;
	gcdchi = gcd(Lp,chi);
	\\ If cst, can disable Tp
	print("Chosen p=",p,", residual degree ",a);
	print("Time choice p: ",timestr(~t0));
	e=2^ceil(log((log(2)+2*D*log(10))/log(p))/log(2)); \\ Power of 2 such that sqrt(1/2*p^e)>10^D
	pe = p^e;
	[N,k] = mfparams(f)[1..2];
	if(k>2,N*=l);
	print("\n--> Initialising modular Jacobian");
	J=ModJacInit(N,H,p,a,e,qprec,Lp,UseTp);
	print("Time ModJacInit: ",timestr(~t0));
	print("Size J: ",mysize(sizebyte(J)));
	if(threadlim,default(nbthreads,threadlim[2]));
	print("\n--> Getting basis of representation space over F_",p);
	J1 = PicRed(J,1);
	if(UseTp,
		print("Using Tp");
		chiOO = chi;
		while(1,
			gcdchi = gcd(Lp/chiOO,chi);
			if(poldegree(gcdchi)==0,break);
			chiOO *= gcdchi
		);
		B = TorsBasis(J1,l,Lp,chiOO,1);
		[B,matFrob,LinTests,R,FRparams] = B;
		print("The matrix of Frob_",p," is");
  	printp(centerlift(matfrobenius(Mod(matFrob,l))));
  	i=1;M=Mod(matFrob,l);
  	while(M!=1,M*=matFrob;i++);
  	print("It has order ",i);
		TpB = parapply(W->PicTp(J1,W),B);
		TpR = matconcat(apply(W->Tors_TestPt(J1,W,l,LinTests,FRparams),TpB)); \\ TODO parapply?
		Tp = Mod(R,l)^(-1)*TpR; \\ Matrix of Tp on B
		print("The matrix of Tp is");
    printp(centerlift(Tp));
		KTp = matker(Tp-ap);
		print("Eigenspace of dim ",#KTp);
		if(#KTp!=2,error("Bug: Tp does not cut the rep"));
		B = vector(2,j,PicLC(J1,centerlift(KTp[,j]),B)); \\ Eigenspace
		matFrob = matker(matconcat([KTp,matFrob*KTp]));
		matFrob = -matFrob[3..4,]^(-1)*matFrob[1..2,];
	,
		[B,matFrob] = TorsBasis(J1,l,Lp,chi,0); \\ Basis of the mod p^1 space and matrix of Frob_p
	);
	print("The matrix of Frob_",p," is");
	printp(centerlift(matfrobenius(Mod(matFrob,l))));
	i=1;M=Mod(matFrob,l);
	while(M!=1,M*=matFrob;i++);
	print("It has order ",i);
	[WB,cWB] = TorsSpaceFrobGen(J1,l,B,matFrob);
	print("Time getting basis over F_",p,": ",timestr(~t0));
	J1 = B = 0;
	if(threadlim,default(nbthreads,threadlim[3]));
	print("\n--> Lifting ",#WB," points ",p,"-adically, target O(",p,"^",e,")");
	WB = apply(W->PicLiftTors(J,W,1,l),WB);
	print("Time lifting: ",timestr(~t0));
	print("Size W: ",mysize(sizebyte(WB[1])));
	if(threadlim,default(nbthreads,threadlim[4]));
	print("\n--> All of representation space");
	Z = TorsSpaceFrobEval(J,WB,cWB,l,matFrob);
	WB = 0;
	print("Time span and eval: ",timestr(~t0));
	print("\n--> Expansion and identification");
	AF = TorsSpaceGetPols(Z,l,matFrob,JgetFrobMat(J),JgetT(J),pe,p,e);
	print("Time polynomials: ",timestr(~t0));
	if(threadlim,default(nbthreads,def_threads));
	AF = AF[1];
	[AF[3],AF[1]];
}
