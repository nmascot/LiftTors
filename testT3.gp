read("install.gp");
read("MakTorsSpace.gp");

WeiRed(f,h)=
{
 my(F=f+(h/2)^2);
 2^poldegree(F)*subst(F,x,x/2);
}


p=7;a=8;e=100;
f=x^8 + x + 3;
p=79;a=6;e=1024;l=5;d=4;
f=x^6+x+1;

f=WeiRed(x^6-3*x^5+2*x^4+x^3-x,1);
p=7;l=3;a=4;e=100;d=2;

J=HyperInit(f,p,a,e);
J1=PicRed(J,1);

GalRepBasis(J,l,C)=
{
	my(d,B,nB,W,W0);
	d=if(C,poldegree(C),2*Jgetg(J));
	B=vector(d);
	nB=0;
	W0=JgetW0(J);
	while(nB<d,
		W = HyperPicRandTors(J,f,l,C);
		W = PicChord(J,W,W0,1);
		nB+=1;
		B[nB] = W;
		print("Got new point.");
		print(nB);
		if(nB>1,
			R = []; \\ PicTorsRels(J,B[1..nB],l,1);
			if(#R,print("Unfortunately,\n",R);nB-=1);
		);
	);
	B;
}

WB = GalRepBasis(J1,l,'x^2+3*'x+7);
print("Lifting");
my(J=J,l=l);WB = parapply(W->PicLiftTors(J,W,1,l),WB);
print("All the space");
V = TorsSpace(J,WB,l);
print("Evaluation");
my(J=J);Z=parapply(W->HyperPicEval(J,W)[1],V[1..l^d-1]);
F=factorback(apply(u->'x-u,Mod(Z,Jgetpe(J))*Mod(1,JgetT(J))));
F=liftall(F)*(1+O(p^e));
bestappr(F);
