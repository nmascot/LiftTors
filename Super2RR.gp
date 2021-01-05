LSuper(n,d,m)= /* L(n*oo) for y^m=f(x) where f has deg d */
{ /* Assumes f sqfree and (m,d) coprime, then x has deg m and y has deg d */
	my(L=List(),mb=min(m-1,floor(n/d)));
	for(b=0,mb,
		for(a=0,(n-d*b)/m,
			listput(L,'x^a*'y^b)
		)
	);
	Vec(L);
}

Super2RR(f0,m,P)= \\ y^m = f(x), requies d,m coprime and f sqfree
{ \\ P must be a rational point
	my(v,f,d,g,d0,L,LL,L1,L2,E2,Auts,KnownAuts);
	v = variable(f0);
	f = subst(f0,v,'x); \\ make sure var(f) == 'x
	if(!issquarefree(f),error(f0," is not squarefree"));
	d = poldegree(f);
	if(gcd(d,m)>1,error(m," is not coprime with the degree of ",f0));
	if(subst(f,'x,P[1])!=P[2]^m,error("The point ",P," is not on the curve y^",m,"=",f));
	g = (d-1)*(m-1)/2;
	d0 = 2*g+1;
	L = LSuper(d0,d,m);
	LL = LSuper(2*d0,d,m);
	L1 = LSuper(d0+g,d,m);
	L2 = LSuper(d0+g+1,d,m);
	E2 = subst(subst(L2,'x,P[1]),'y,P[2]);
	L2 = L2*matker(Mat(E2));
	if(Mod(m,2),
		Auts=[];
		KnownAuts=0;
	,
		Auts=[['x,-'y,1]];
		KnownAuts=if(m==2,[-1],0);
	);
	[y^m-f,Auts,KnownAuts,g,d0,L,LL,L1,L2];
}
