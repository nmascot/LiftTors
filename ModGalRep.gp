read("install.gp");
read("GalRep.gp");

IDpolcyclo(Phi,N)=
{
	my(var=variable(Phi));
	fordiv(eulerphi(N),d,
		if(eulerphi(d)!=poldegree(Phi),next);
		if(Phi == polcyclo(d,var),return(d))
	);
}

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


mfbestp(f,l,coeffs,pmax)=
{
	my(N,k,eps,P,Phi,yl,tl,o,ZNX,lN,H,best,ap,epsp,chi,Lp,Psi,a1,a2,a);
	[N,k,eps,P,Phi]=mfparams(f);
	[yl,tl] = mfyt(f,l,coeffs);
	o = IDpolcyclo(Phi,N);
	ZNX = znstar(N,1);
	eps = znchar(f);
	lN = if(k==2,N,l*N);
	H = select(h->gcd(lN,h)==1,[1..lN-1]);
	H = select(h->Mod(h^(k-2),l)*chareval(eps[1],eps[2],h,[tl,o])==1,H);
	best = [];
	forprime(p=5,pmax,
		if(Mod(l*N,p)==0,next);
		ap = mfcoef(f,p);
		ap = subst(subst(liftall(ap),'t,tl),'y,yl);
		epsp = chareval(eps[1],eps[2],p,[tl,o]);
		chi = 'x^2-ap*'x+p^(k-1)*epsp;
		Lp = LMod(lN,H,p);
		Psi = Mod(Lp,l)/chi;
		if(poldegree(denominator(Psi)),error("Bug in charpoly"));
		if(poldegree(gcd(chi,Psi)),print("p=",p, " has multiplicity");next);
		a1 = mordroot(chi,l);
		a2 = znorder(Mod(p,lN));
		a = lcm(a1,a2);
		print("p="p,": Needs deg ",a," (",a1," to split rep, ",a2," for roots of 1), log #J=",round(log(polresultant(Lp,'x^a-1))),")");
		a = lcm(a1,a2);
		if(best==[] || a < best[2], best=[p,a,Lp,chi]);
	);
	[H,best];
}

mfgalrep(f,l,coeffs,pmax,D)=
{
	my(N,k,H,best,p,a,Lp,chi,e,J,CuspsQ,J1,B,matFrob,i,M,WB,cWB,TI,Z,AF);
	[H,best] = mfbestp(f,l,coeffs,pmax);
	[p,a,Lp,chi] = best;
	e=2^ceil(log((log(2)+2*D*log(10))/log(p))/log(2)); \\ Power of 2 such that sqrt(1/2*p^e)>10^D
	print("O(",p,"^",e,")");
	[N,k] = mfparams(f)[1..2];
	if(k>2,N*=l);
	J=ModJacInit(N,H,p,a,e);
	J1 = PicRed(J,1);
	[B,matFrob] = TorsBasis(J1,l,Lp,chi); \\ Basis of the mod p^1 space and matrix of Frob_p
	print("The matrix of Frob is");
	printp(centerlift(matfrobenius(Mod(matFrob,l))));
	i=1;M=Mod(matFrob,l);
	while(M!=1,M*=matFrob;i++);
	print("It has order ",i);
	[WB,cWB] = TorsSpaceFrobGen(J1,l,B,matFrob);
	print("\n--> Lifting ",#WB," points ",p,"-adically");
	WB = apply(W->PicLiftTors(J,W,1,l),WB);
	print("\n--> All of T");
	TI = TorsSpaceFrob(J,WB,cWB,l,matFrob);
	print("\n--> Evaluation of ",#TI[2]," points");
	Z = TorsSpaceFrobEval(J,TI,l,2,matFrob);
	AF = TorsSpaceGetPols(J,Z);
	AF[1];
}
