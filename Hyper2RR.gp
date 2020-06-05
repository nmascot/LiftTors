read("PtIsOnCurve.gp");

HyperRR(n,g,u,v)=
{
	concat(vector(n+1-poldegree(u,'x),i,'x^(i-1)*u),(y-v)*vector(n-g,i,x^(i-1)));
}

Hyper2RR(f0,P1,P2)= /* y^2=f0(x). P1,P2 rat pts, not conjugate by hyper invol. */
{
	my([x1,y1]=P1,[x2,y2]=P2,f,h,g,d0,L,LL,L1,L2);
	f = subst(f0,variable(f0),'x);
	if(!PtIsOnHyperCurve(f,P1),error("The point ",P1," is not on this hyperelliptic curve."));
	if(!PtIsOnHyperCurve(f,P2),error("The point ",P2," is not on this hyperelliptic curve."));
	if(type(f)=="t_VEC" && #f==2, \\ Long Weierstrass equation
		[f,h]=f;
		y1=2*y1+subst(h,'x,x1);
		y2=2*y2+subst(h,'x,x2);
		f=4*f+h^2
	);
	d0 = poldegree(f);
	if(poldegree(f)%2, \\ Odd degree: change model
		while(polcoef(f,0)==0 || x1==0 || x2==0,
			f=subst(f,x,x+1);
			x1-=1;
			x2-=1
		);
		f='x*polrecip(f);
		d0 += 1;
		x1 = 1/x1;
		x2 = 1/x2;
		y1 /= x1^(d0/2);
		y2 /= x2^(d0/2);
	);
	g=(d0-2)/2;
	L=HyperRR(g+1,g,1,0);
	LL=HyperRR(2*g+2,g,1,0);
	if(g%2,
		L1=HyperRR(3*(g+1)/2,g,'x-x1,y1);
		L2=HyperRR(3*(g+1)/2,g,'x-x2,y2);
	,
		L1=HyperRR(3*g/2+2,g,('x-x1)*('x-x2),(y2-y1)/(x2-x1)*'x+(y1*x2-y2*x1)/(x2-x1));
		L2=HyperRR(3*g/2+1,g,1,0)
	);
	[y^2-f,g,d0,L,LL,L1,L2];
}


RREval(J,L)=
{
  my(W,Z,T,p);
  Z=JgetZ(J);
  T=JgetT(J);
  p=Jgetp(J);
  W=matrix(#Z,#L,P,j,subst(subst(L[j],'x,Z[P][1]),'y,Z[P][2]));
  liftall(Mod(Mod(W,p),T));
}

MakIn2(J,A,B)=
{ \\ [{A,y-B}-OO] in Makdisi form
  my(x='x,y='y,W1,W2);
  \\ V2 = M2 = L(2D0) = L(6OO)
  W1=RREval(J,[A,x*A,x^2*A,y-B,x*(y-B)]); \\ {A,y-B}+2*OO
  W2=RREval(J,[1,x,x^2,x^3,y]); \\ 3*OO
  PicSub(J,W1,W2); \\ {A,y-B}-OO
}

Pts2Cantor(P,Q,T,p)=
{ \\ Cantor form of P+Q
  my(xP,yP,xQ,yQ,A,B);
  [xP,yP]=Mod(Mod(P,p),T);
  [xQ,yQ]=Mod(Mod(Q,p),T);
  if(xP==xQ,error("xP=xQ not implemented in Cantor"));
  A=(x-xP)*(x-xQ);
  B=((yQ-yP)*x+(xQ*yP-xP*yQ))/(xQ-xP);
  [A,B];
}

MakIn1(J,P,Q,P0)=
{ \\ Class of [P-Q] in Makdisi form
  \\ P0 auxiliary point
  \\ P,Q,P0 must have distinct x-coords (for convenience)
  my(Z,T,p,A,B);
  T=JgetT(J);
  p=Jgetp(J);
  if(Mod(Mod(P,p),T)==Mod(Mod(Q,T),p),return(JgetW0(J)));
  if(Mod(Mod(P,p),T)==Mod(Mod(P0,T),p),error("P=P0 not implemented"));
  if(Mod(Mod(Q,p),T)==Mod(Mod(P0,T),p),error("Q=P0 not implemented"));
  [A,B] = Pts2Cantor(P,P0,T,p);
  WP = RREval(J,[A,x*A,x^2*A,y-B,x*(y-B)]); \\ P+P0+2*OO
  [A,B] = Pts2Cantor(Q,P0,T,p);
  WQ = RREval(J,[A,x*A,x^2*A,y-B,x*(y-B)]); \\ Q+P0+2*OO
  PicSub(J,WP,WQ); \\ P-Q
}

HyperRandPt(f,T,p)=
{ \\ Random point on C(Fq), Fq=Fp[t]/T(t)
  my(u,v);
  while(1,
    u = random(p*T);
    u = Mod(Mod(u,p),T);
    v = subst(f,'x,u);
    if(issquare(v,&v),
      return([u,v])
    )
  );
}

RandPtFromC(J,l,P0)=
{ \\ Find a point in J[l] of the from M*[P-P0], P rand pt on curve
  my(f,a,T,p);
  f = 'y^2-J[1];
  f = subst(f,'y,0);
  T = JgetT(J);
  p = Jgetp(J);
  a = poldegree(T);
  my(NJ,v,M);
  NJ = polresultant(hyperellcharpoly(Mod(f,p)),'x^a-1);
  v = valuation(NJ,l);
  M = NJ/l^v;
  my(Q,R,W,lW);
	while(1,
  	Q = HyperRandPt(f,T,p);
  	while(1,
    	R = HyperRandPt(f,T,p);
    	if(Q[1]!=R[1] && P0[1]!=R[1],break)
  	);
  	W = MakIn1(J,Q,P0,R);
  	W = PicMul(J,W,M,2);
		if(PicIsZero_val(J,W)==0,break);
		print("RandPtFromC got 0, retrying");
	);
  lW = PicMul(J,W,l,2);
  while(PicIsZero_val(J,lW)==0,
    W = lW;
    M *= l;
    lW = PicMul(J,W,l,2);
  );
  [Q,M,W];
}

ApplyEndo(J,A,B,P,M)=
{ \\ Apply endo to M*[P-P0]
  my(AP,BP,W);
  AP = subst(subst(A,'u,P[1]),'v,P[2]);
  BP = subst(subst(B,'u,P[1]),'v,P[2]);
  W = MakIn2(J,AP,BP);
  PicMul(J,W,M,2);
}

RandTorsPtEndo(J,A,B,l,P0)=
{
	my(P,M,W,TW);
	[P,M,W] = RandPtFromC(J,l,P0);
	TW = ApplyEndo(J,A,B,P,M);
	[W,TW];
}

