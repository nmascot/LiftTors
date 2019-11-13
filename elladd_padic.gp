install("ZpXQ_div","GGGGGL");

elladd_padic(a4,P,Q,T,pe,p,e)=
{
	if(P==[0],return(Q));
	if(Q==[0],return(P));
	my([xP,yP]=P,[xQ,yQ]=Q,dx,dy,l,m,xR,yR);
	if(Mod(xP,p)==Mod(xQ,p),
		print("elladd_padic dangerous");
		if(xP!=xQ,error("elladd_padic bad red"));
		if(yP+yQ==0,return([0]));
		dx = 2*yP;
		dy = 3*xP^2+a4;
	,
	 dx = xQ-xP;
	 dy = yQ-yP
	);
	dx = liftall(dx);
	dy = liftall(dy);
	l = ZpXQ_div(dy,dx,T,pe,p,e);
	l = Mod(Mod(l,T),pe);
	m = yP - l*xP;
	xR = l^2-(xP+xQ);
	yR = l*xR+m;
	[xR,-yR];
}


ellmul_padic(a4,P,n,T,pe,p,e)=
{
	my(Q,m);
	if(P==[0],return([0]));
	if(n==0,return([0]));
	if(n==1,return(P));
	if(n<0,
		Q = ellmul_padic(a4,P,-n,T,pe,p,e);
		if(Q==[0],return(Q));
		return([Q[1],-Q[2]])
	);
	if(n%2,
		m = (n-1)/2;
		Q = ellmul_padic(a4,P,m,T,pe,p,e);
		Q = elladd_padic(a4,Q,Q,T,pe,p,e);
		elladd_padic(a4,P,Q,T,pe,p,e)
	,
		m = n/2;
		Q = ellmul_padic(a4,P,m,T,pe,p,e);
    elladd_padic(a4,Q,Q,T,pe,p,e)
	);
}
