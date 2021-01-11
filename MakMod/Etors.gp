SplitDeg(F,p)=
{
	my(x=variable(F),y=Mod(x,F)^p,i=1);
	while(y!=x,
		y=y^p;
		i++;
	);
	i;
}

EtorsDefInDeg(E,l,v,p,T)=
{
	my(d=poldegree(T),D,X,q2,r);
  D = elldivpol(E,l^v)/elldivpol(E,l^(v-1)); \\ Check manually that E[l^v] def over Fq
	r = SplitDeg(D,p);
	if(Mod(d,r),return(0));
  X = polrootsmod(D,[T,p]);
  if(#X<poldegree(D),
  	if(default(debug),print("E[",l,"^",v,"]/+-1 not split: ",#X," roots out of ",poldegree(D)));
    return(0)
  );
  if(l!=2 || v!=1,
		q2 = (p^d-1)/2;
  	for(j=1,#X,
    	if((X[j]^3+E.a4*X[j]+E.a6)^q2!=1,
      	if(default(debug),print("E[",l,"^",v,"] not split"));
      	return(0)
    	)
		)
  );
	r;
}

Esplit(p,N,d,Avoidj)= \\ Look for an ell. curve / Fp such that E[N] splits over Fq, q=p^d.
{
  my(q=p^d,q2=(q-1)/2,T=ffinit(p,d,'t),a,b,E,j,ap,nu,X,g,fabad,M,fagood,l,v,c,r,M1);
  if(Mod(p^d,N)!=1,error("Impossible by Weil pairing"));
	faN = factor(N);
  while(1,
    if(default(debug),print("---- New curve ----"));
    a=random(Mod(1,p));
    b=random(Mod(1,p));
		if(a==0||b==0,next); \\ Avoid ramification X(N)->X(1)
    if(4*a^3+27*b^2==0,
			if(default(debug),print("Singular"));
			next
		);
    E=ellinit([a,b]);
		j=E.j;
		for(i=1,#Avoidj,
			if(Mod(Avoidj[i],p)==j,next(2))
		);
		ap = ellap(E);
		nu = polsym(Mod('x^2-ap*'x+p,N),d);
		if(Mod(q+1-nu[d+1],N^2),next); \\ Must have N² | #E
		g=gcd(ap^2-4*p,N);
		fabad=factor(g); \\ Primes l s.t. Frob_p not 1/2simple on l-tors
		fabad[,2]=apply(l->valuation(N,l),fabad[,1]); \\ Set multiplicities in N
		M=N/factorback(fabad); \\ Part of N where Frob 1/2simple
		fagood = factor(M);
		if(M>1,
      if(default(debug),print("Checking Frob^",d," trivial on E[",M1,"]"));
			if(Mod(nu[d+1],M)!=2, \\ Check if Frob^d=1 on E[M1] by testing a^d+b^d==2
        if(default(debug),print("Frob^",d," not trivial on E[",M1,"]"));
        next
      )
		);
		for(i=1,#fabad~, \\ Check if Frob^d unipotent on E[l]
			l=fabad[i,1];
			if(l==2,next);
			if(default(debug),print("Checking Frob^",d," unipotent on E[",l,"]"));
			\\c=ap/Mod(2,lv); \\ Frob = c*unipotent on E[l]
			if(Mod(ap,l)^d!=Mod(2,l)^d, \\ Check if Frob^d unipotent on E[l]
				if(default(debug),print("Frob^",d," not unipotent on E[",l,"]"));
				next(2)
			);
		);
		/* So now Frob^d trivial on E[l] but maybe not E[l^v] for l good,
			and unipotent on E[l] but maybe not E[l^v] for l bad.
		 We now check for each l if Frob^d trivial on E[l^v],
		 and also determine the order of Frob on E[l^v]/-1. */
		/* TODO l=2? */
		c=1;
		for(i=1,#fagood~,
			[l,v] = fagood[i,];
			if(v==1,
				\\ Already know Frob^d trivial on E[l^v]=E[l]. Check order on E[l]/-1 using charpoly.
				for(r=1,d, \\ TODO fordiv?
					if(Mod(nu[r+1],l)^2==4 && Mod(p,l)^r==1,
						print("E[",l,"^",v,"]/+-1 splits in deg ",r);
						c = lcm(c,r);
						break
					)
				)
			,
				r =	EtorsDefInDeg(E,l,v,p,T);
				if(r==0,next(2));
				c = lcm(c,r);
			)
		);
		for(i=1,#fabad~,
			[l,v] = fabad[i,];
			\\ if v==1 && Mod(d,l)==0, then Frob^d automatically trivial, but how do we get r?
			r = EtorsDefInDeg(E,l,v,p,T);
      if(r==0,next(2));
      c = lcm(c,r);
		);
		if(c==d,return([a,b]))
  );
}

\\GetPt(x,f)=[x,sqrt(subst(f,'x,x))]; \\ Given x, find a point [x,y] on the curve

GetOrder(E,P,l)= \\ Given l and a point of order l^v, find v
{
	my(n=0);
	while(P!=[0],
		P = ellmul(E,P,l);
		n+=1;
	);
	n;
}

GetzOrder(z,l)= \\ Same as above for roots of 1
{
  my(n=0);
  while(z!=1,
    z = z^l;
    n+=1;
  );
  n;
}

zlog(x,z)=
{
	my(n=0,zn=1);
	while(x!=zn,
		zn *= z;
		n++
	);
	n;
}

ELocalBasis(E,t1,ab,l,v,D)= \\ l,v -> Basis [P,Q] of E[l^v], plus its Weil pairing, mat of Frob, and the l^v-division polynomial 
{
	my(p=t1.p,lv=l^v,X,x,y,P,Q,z,FP,FQ,MFrob);
	if(default(debug),print("---- Getting basis of E[",l,"^",v,"] ----"));
	X = polrootsmod(D,t1);
	while(1,
		x = X[random(#X)+1];
		y = sqrt(x^3+ab[1]*x+ab[2]);
		P = [x,y];
		x = X[random(#X)+1];
    y = sqrt(x^3+ab[1]*x+ab[2]);
    Q = [x,y];
		\\print("P has order ",l,"^",GetOrder(E,P,l));
		\\print("Q has order ",l,"^",GetOrder(E,P,l));
		z = ellweilpairing(E,P,Q,lv);
		if(default(debug),print("z has order ",l,"^",GetzOrder(z,l)));
		if(z^(l^(v-1))!=1,
			FP = [P[1]^p,P[2]^p];
			FQ = [Q[1]^p,Q[2]^p];
			MFrob = matrix(2,2);
			MFrob[1,1] = zlog(ellweilpairing(E,FP,Q,lv),z);
			MFrob[2,1] = -zlog(ellweilpairing(E,FP,P,lv),z);
			MFrob[1,2] = zlog(ellweilpairing(E,FQ,Q,lv),z);
			MFrob[2,2] = -zlog(ellweilpairing(E,FQ,P,lv),z);
			return([P,Q,1/z,Mod(MFrob,lv)]) \\ 1/z because of bug #2187.
		)
	);
}

ffLiftRoot(F,a,e)=
{
	my(p=a.p,T=a.mod,b=a.pol);
	b = padicappr(F,Mod(b,T)+O(p^e))[1];
	Mod(liftint(b),p^e);
}

LiftTorsPt(P,ab,D,e)=
{
	my(x,y,T,p);
	[x,y] = P;
	T = x.mod;
	p = x.p;
	x = ffLiftRoot(D,x,e);
	y = ffLiftRoot('y^2-liftint(x^3+ab[1]*x+ab[2]),y,e);
	[x,y];
}

EBasis(N,p,T,t1,e,Avoidj)=
{ /* T irr mod p, t1 = ffinit(p,T). Get [E/Q,P,Q,eN(P,Q),MFrob] mod p^e */
	my(pe=p^e,Tpe,ab,f,EQ,EFq,faN,Pk,Qk,zk,P,Q,z,l,v,lv,D,MFrob,MFroblv);
	Tpe=[T,pe,p,e];
	ab = liftint(Esplit(p,N,poldegree(T),Avoidj));
	EQ=ellinit(ab);
	EFq=ellinit(ab,t1);
	faN = factor(N);
	P = [0];
	Q = [0];
	z = 1;
	MFrob = vector(#faN~);
	for(k=1,#faN~,
		[l,v] = faN[k,];
		lv = l^v;
		D = elldivpol(EQ,lv)/elldivpol(EQ,lv/l);
		[Pk,Qk,zk,MFroblv] = ELocalBasis(EFq,t1,ab,l,v,D);
		[Pk,Qk] = apply(R->liftall(LiftTorsPt(R,ab,D,e)),[Pk,Qk]);
		zk = ffLiftRoot('x^lv-1,zk^(N/lv),e);
		P = elladd_padic(ab[1],P,Pk,T,pe,p,e);
		Q = elladd_padic(ab[1],Q,Qk,T,pe,p,e);
		z *= zk;
		MFrob[k] = MFroblv;
	);
	[EQ,P,Q,z,chinese(MFrob)];
}

EMl1(E,N,P0,Q0,m,T,pe,p,e)=
{
	my(EN,todo,done,Ml1,x,y);
	\\ Write down all N-torsion: : this is a naive level structure alpha: (Z/NZ)² ~ E[N]
	EN=matrix(N,N); \\ [[ m P0 + n Q0 ]]
  EN[m,N]=P0;
  EN[N,1]=Q0;
  for(x=2,N-1,
    EN[ZNnorm(m*x,N),N]=elladd_padic(E.a4,EN[ZNnorm(m*(x-1),N),N],P0,T,pe,p,e);
    EN[N,x]=elladd_padic(E.a4,EN[N,x-1],Q0,T,pe,p,e);
  );
  todo = List();
  for(x=1,N-1,for(y=1,N-1,listput(todo,[x,y])));
  done = parapply(X->elladd_padic(E.a4,EN[X[1],N],EN[N,X[2]],T,pe,p,e),Vec(todo));
  for(i=1,#todo,[x,y]=todo[i];EN[x,y]=done[i]);
	Ml1=matrix(N,N);
  \\ TODO do we use non coprime entries?
  todo = List();
  for(x=1,N-1,listput(todo,[[x,N],[0,1]])); \\ P=alpha(x,0) -> Q=alpha(0,1)
  for(x=1,N,for(y=1,N-1,listput(todo,[[x,y],[1,0]]))); \\ P=alpha(x,y), y!=0 -> Q=alpha(1,0)
  done = parapply(X->l1(EN,X[1],X[2],T,pe,p,e),Vec(todo));
  for(i=1,#todo,[x,y]=todo[i][1];Ml1[x,y]=done[i]);
	Ml1;
}

GetMl1(N,Pts,PtTags,p,T,t1,e,zNpref,Avoidj)=
{
	my(a=poldegree(T),pe=p^e,E,P0,Q0,zN,m,MFrobE,tMFrobE,PtsFrob,Ml1);
	[E,P0,Q0,zN,MFrobE] = EBasis(N,p,T,t1,e,Avoidj);
	print("E:y²=x³+",E.a4,"x+",E.a6," with accuracy O(",p,"^",e,") and residual degree a=",a);
	\\ Get action of Frob on fibre
	tMFrobE = mattranspose(MFrobE);
	if(zNpref,
		m = zlog(zN,zNpref);
		tMFrobE = matdiagonal([1/m,1])*tMFrobE*matdiagonal([m,1]);
		zN = 0
	,
		m=1
	);
  PtsFrob = Vecsmall(apply(P->GetCoef(PtTags,liftint(P*tMFrobE)),Pts)); \\ Frob([c,d]) = [c,d]*(t^MFrobE)
	if(default(debug),print("Ml1"));
	Ml1 = EMl1(E,N,P0,Q0,m,T,pe,p,e);
	[E,Ml1,PtsFrob,zN];
}
