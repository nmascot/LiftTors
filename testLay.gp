read("install.gp");
read("MakTorsSpace.gp");
f=4*WeiRed(x^5+x^4,x^3+x+1);
p=61;a=6;e=100;l=7;
\\C = x^2-2*x-1;
\\C = x^2+3*x+3;
f=x^6+x+1;p=29;
C=4;a=18;l=2;
J=HyperInit(f,p,a,1);

PicLC(J,C,W)=
{
	/* TODO efficiency */
	my(S);
	if(#C==0,return(JgetW0(J)));
	S = PicMul(J,W[1],C[1],2);
	for(i=2,#C,
		S = PicAdd(J,S,PicMul(J,W[i],C[i],2));
	);
	S;
}

TorsOrd(J,W,l)=
{
	my(v=0,lW);
	v = 0;
	lW = W;
	while(!PicIsZero(J,lW),
		v += 1;
		W = lW;
		lW = PicMul(J,W,l,0)
	);
	[W,v];
}

RandTorsPt(J,f,M,chiC)=
{
	my(W);
	while(1,
		W = HyperPicRand(J,f);
		W = PicMul(J,W,M,0);
		if(chiC,W = PicFrobPoly(J,W,chiC));
		if(!PicIsZero(J,W),return(W));
	);
}

PicTorsTrueRel(J,W,l)=
{
	my(R,ex=0);
	print("Looking for rels btween ",#W," pts");
	while(1,
		R = PicTorsRels(J,W,l,ex);
		print("Found ",#R," rels");
		if(#R==0,return(R));
		if(#R==1,
			if(PicIsZero(J,PicLC(J,R[,1],W)),return(R));
			print("False positive");
		);
		ex += 1
	);
}

TorsBasis(J,f,p,a,l,C)=
{
	my(d,chi,N,M,v,chiC,W,o,T,BW,Bo,BT,R,m,S);
	chi = hyperellcharpoly(Mod(f,p));
	N = polresultant(chi,'x^a-1);
	v = valuation(N,l);
	M = N/l^v;
	d = poldegree(C);
	if(d<=0,
		d = C;
		chiC = 0;
	,
		chiC = lift(Mod(chi,l)/C);
		fa = [C,chiC];
		fa = polhensellift(chi,fa,l,v);
		chiC = fa[2];
		chiC = centerlift(Mod(chiC,l^v))
	);
	BW = vector(d);
	Bo = vector(d);
	BT = vector(d);
	r = 0;
	while(r<d,
		print("Status:",Bo[1..r]);
		print("Getting new point");
		W = RandTorsPt(J,f,M,chiC);
		[T,o] = TorsOrd(J,W,l);
		print("It has order l^",o);
		r += 1;
		BW[r] = W;
		Bo[r] = o;
		BT[r] = T;
		if(r==1,next);
		while(1,
			print("Looking for relations...");
			R = PicTorsTrueRel(J,BT[1..r],l);
			print("Found relation: ",R);
			if(#R == 0,print("Good, no relation");next(2));
			R = R[,1];
			m = vecmin([Bo[i]|i<-[1..r],R[i]]);
			print("Dividing relation ",R," by l^",m); 
			S = vector(r,i,if(R[i],l^(Bo[i]-m)*R[i],0));
			W = PicLC(J,S,BW[1..r]);
			[T,o] = TorsOrd(J,W,l);
			print("gives point of order l^",o);
			if(o==0,
				print("Giving up this point.");
				r -= 1;
				break
			);
			BW[r] = W;
			Bo[r] = o;
			BT[r] = T;
		);
	);
	BT;
}

B = TorsBasis(J,f,p,a,l,C);

/*N = ordJ(f,p,a);
v = valuation(N,l);
M = N/l^v;
W11 = HyperPicRand(J,f);
W11 = PicMul(J,W11,M,0);
W21 = HyperPicRand(J,f);
W21 = PicMul(J,W21,M,0);
chi = hyperellcharpoly(Mod(f,p));
chiC = liftint(chi/Mod(C,l));
fa = [C,chiC];
fa = polhensellift(chi,fa,l,v);
W11 = PicFrobPoly(J,W11,fa[2]);
W21 = PicFrobPoly(J,W21,fa[2]);*/


