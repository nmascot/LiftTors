DimData(c,N)= \\ [Pindex,e2,e3,eOO]
{
	my(fa,n,p,v,a);
	if(N==1,return([1,1,1,1]));
	fa=factor(N);
	n=#fa~;
	if(c==-1,
		a=if(N==2,3,N^2/2*prod(k=1,n,[p,v]=fa[k,];(1-1/p^2)));
		return([N*a,0,0,a])
	);
	if(N==2,return([3,1,0,2])); \\ Gamma1(2) eq Gamma0(2)
	if(c==0,
		return([N*prod(k=1,n,[p,v]=fa[k,];1+1/p),
			if(N%4,prod(k=1,n,[p,v]=fa[k,];1+if(p==2,0,(-1)^((p-1)/2))),0),
			if(N%9,prod(k=1,n,[p,v]=fa[k,];1+centerlift(Mod(p,3))),0),
			prod(k=1,n,[p,v]=fa[k,];if(v%2,2*p^((v-1)/2),(p+1)*(p^(v/2-1))))]);
	);
	if(c==1,
		if(N==3,return([4,0,1,2]));
		if(N==4,return([6,0,0,3]));
		return([N^2/2*prod(k=1,n,[p,v]=fa[k,];1-1/p^2),0,0,N/2*prod(k=1,n,[p,v]=fa[k,];(1-1/p)*(v+1-(v-1)/p))]);
	);
}

Genus(c,N)=my([d,e2,e3,e0]=DimData(c,N));1+d/12-e2/4-e3/3-e0/2;

DimSk(k,c,N)=
{
  my(d,e2,e3,e0,g);
  if(k<=0,return(0));
	[d,e2,e3,e0]=DimData(c,N);
  g=1+d/12-e2/4-e3/3-e0/2;
  if(k==2,return(g));
  if(k%2,
		if(c==0||N<=2,return(0)); \\ -1 in Gamma
		if(k>=3,
    	return((k-1)*(g-1)+floor(k/3)*e3+if(c==1&&N==4,(3*k-5)/2,(k/2-1)*e0))
		);
		\\ Now k=1
		if(e0>2*g-2,return(0));
		error("Unimplemented");
	,
		(k-1)*(g-1)+floor(k/4)*e2+floor(k/3)*e3+(k/2-1)*e0;
  );
}

DimEk(k,c,N)=
{
  my(d,e2,e3,e0,g);
  if(k<0,return(0));
  if(k==0,return(1));
  [d,e2,e3,e0]=DimData(c,N);
  g=1+d/12-e2/4-e3/3-e0/2;
	er=if(c==1 && N==4,2,e0); \\ Regular cusps
  if(k==2,return(e0-1));
  if(k%2,
    if(c==0||N<=2,return(0));
		er*if(k>=3,1,1/2);
  ,
    e0
  );
}

DimMk(k,c,N)=DimSk(k,c,N)+DimEk(k,c,N);
