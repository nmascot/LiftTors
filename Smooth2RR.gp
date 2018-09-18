LSmooth(n,d,x,y)=
{
	my(L=List(),mu=min(d-1,n));
	for(b=0,mu,
		for(a=0,n-b,
			listput(L,x^a*y^b)
		)
	);
	Vec(L);
}

SmoothGeneric(f0,d,P01,P02)=
{
	my(A,f,x,y,u,v,w,P1,P2);
	[x,y] = variables(f0);
	if(polcoeff(f0,d,x)!=0 && polcoeff(f0,d,y)!=0,return([f0,P01,P02]));
	f = f0;
	A = matid(3);
	while(1,
		A = Mat(0);
		while(matdet(A)==0,
			A = matrix(3,3,i,j,random(3)-1);
		);
		F = matrix(d+1,d+1,i,j,polcoeff(polcoeff(f0,i-1,x),j-1,y));
		[u,v,w] = A*[x,y,1]~;
		f = sum(i=0,d,sum(j=0,d,F[i+1,j+1]*u^i*v^j*w^(d-i-j)));
		if(polcoeff(f,d,x)==0 || polcoeff(f,d,y)==0,next);
		A = A^-1;
		P1 = P01;
		for(i=1,#P1,
			[u,v,w] = A*P1[i]~;
			if(w==0,next(2));
			P1[i] = [u,v]/w;
		);
		P2 = P02;
		for(i=1,#P2,
			[u,v,w] = A*P2[i]~;
			if(w==0,next(2));
			P2[i] = [u,v]/w;
		);
		return([f,P1,P2]);
	);
}

TotalDeg(f,x,y)=
{
	my(d=-1,t);
	for(i=0,poldegree(f,x),
		t = polcoeff(f,i,x);
		if(t,
			d = max(d,poldegree(t,y)+i)
		)
	);
	d;
}

Smooth2RR(f0,P01,P02)=
{ \\ P01,P02 should be lists of rat pts (TODO for now distinct)
	\\ d even: of size d/2-1
	\\ d odd : of size d-1
	my(x,y,d,g,d0,L,M,L1,L2);
	[x,y] = variables(f0);
	d = TotalDeg(f0,x,y);
	[f,P1,P2] = SmoothGeneric(f0,d,P01,P02);
	g = (d-1)*(d-2);
	g = g/2;
	d0 = (d-2)*d;
	L = LSmooth(d-2,d,x,y);
	LL = LSmooth(2*(d-2),d,x,y);
	M = LSmooth(if(d%2,(3*d-5)/2,3*(d/2-1)),d,x,y);
	L1 = apply(P->subst(subst(M,x,P[1]),y,P[2]),P1)~;
	L1 = M*matker(matconcat(L1));
	L2 = apply(P->subst(subst(M,x,P[1]),y,P[2]),P2)~;
	L2 = M*matker(matconcat(L2));
	[f,g,d0,L,LL,L1,L2];
}

