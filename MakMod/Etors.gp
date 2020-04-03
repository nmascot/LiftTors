\\read("elladd_padic.gp");

Esplit(p,N,d)= \\ Look for an ell. curve / Fp such tht E[N] splits over Fq, q=p^d.
{
  my(q=p^d,q2=(q-1)/2,T=ffinit(p,d,'t),a,b,E,ap,X,g,M,L,V,l,v,lv,c,D,faM,M1);
  if(Mod(p^d,N)!=1,error("Impossible by Weil pairing"));
	/*m = valuation(N,2);
	M = N/2^m;*/
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
		ap = ellap(E);
		g=gcd(ap^2-4*p,N);
		L=factor(g)[,1]; \\ Primes l s.t. Frob_p not 1/2simple on l-tors
		V=apply(l->valuation(N,l),L); \\ Multiplicities in N
		M=N/prod(i=1,#L,L[i]^V[i]); \\ Part of N where Frob 1/2simple
		faM = factor(M);
		M1 = prod(i=1,#faM~,faM[i,1]); \\ Radical of M
		if(M1>1,
      if(default(debug),print("Checking Frob^",d," trivial on E[",M1,"]"));
			if(polsym(Mod('x^2-ap*'x+p,M1),d)[d+1]!=2, \\ Check if Frob^d=1 on E[M] by testing a^d+b^d==2
        if(default(debug),print("Frob^",d," not trivial on E[",M1,"]"));
        next
      )
		);
		for(i=1,#L, \\ Check if Frob^d unipotent on E[l^v]
			l=L[i];
			v=V[i];
			if(default(debug),print("Checking Frob^",d," unipotent on E[",l,"^",v,"]"));
			lv=l^v;
			if(l!=2,
				c=ap/Mod(2,lv); \\ Frob = c*unipotent on E[l^v]
				if(c^d!=1, \\ Check if Frob^d unipotent on E[l^v]
					if(default(debug),print("Frob^",d," not unipotent on E[",l,"^",v,"]"));
					next(2)
				);
			);
		);
		L = matconcat([L,V]); \\ facto prod l_i^v_i
		L = select(t->Mod(d,t[1]^t[2]),L~)~; \\ d kill unipotents if l_i^v_i | d
		faM = select(t->t[2]>1,faM~)~; \\ Only keep repeated factors
		\\ So now Frob^d unipotent on E[l^v] for [l,v] in L,
		\\ and Frob^d trivial on E[l] but maybe not E[l^v] for [l,v] in faM
		L = matconcat([L,faM]~);
		for(i=1,#L~,
			if(default(debug),print("Checking Frob^",d," trivial on E[",L[i,1],"^",L[i,2],"]"));
			[l,v] = L[i,];
			lv = l^v;
			D = elldivpol(E,lv)/elldivpol(E,l^(v-1)); \\ Check manually that E[l^v] def over Fq
			X=polrootsmod(D,[T,p]);
			if(#X<poldegree(D),
				if(default(debug),print("E[",l,"^",v,"]/+-1 not split: ",#X," roots out of ",poldegree(D)));
				next(2)
			);
			if(lv!=2,
				for(j=1,#X,
					if((X[i]^3+a*X[i]+b)^q2!=1,
						if(default(debug),print("E[",l,"^",v,"] not split"));
						next(3)
					)
				)
			)
		);
		return([a,b])
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
	print("---- Getting basis of E[",l,"^",v,"] ----");
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
		print("z has order ",l,"^",GetzOrder(z,l));
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

EBasis(N,p,a,e)=
{
	my(T,t1,pe=p^e,Tpe,ab,f,EQ,EFq,faN,Pk,Qk,zk,P,Q,z,l,v,lv,D,MFrob,MFroblv);
	T=ffinit(p,a,'t);
	t1=ffgen(T,'t);
	T = lift(T);
	Tpe=[T,pe,p,e];
	ab = liftint(Esplit(p,N,a));
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
		print(MFroblv);
		MFrob[k] = MFroblv;
	);
	[EQ,P,Q,z,chinese(MFrob)];
}
