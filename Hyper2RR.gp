read("PtIsOnCurve.gp");

HyperRR(n,g,u,v)=
{
	concat(vector(n+1-poldegree(u,'x),i,'x^(i-1)*u),(y-v)*vector(n-g,i,x^(i-1)));
}

Hyper2RR(f0,P1,P2)= /* y^2=f0(x). P1,P2 rat pts, not conjugate by hyper invol. */
{ /* TODO allow extra auts */
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
	[y^2-f,[[x,-y,1]],g,d0,L,LL,L1,L2];
}
