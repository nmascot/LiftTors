read("install.gp");
read("MakTorsSpace.gp");
p=7;a=8;e=512;l=3;C=x^2+3*x+7;d=2;
f=WeiRed(x^6-3*x^5+2*x^4+x^3-x,1);
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
			R = PicTorsRels(J,B[1..nB],l,1);
			if(#R,print("Unfortunately,\n",R);nB-=1);
		);
	);
	B;
}

WB = GalRepBasis(J1,l,C);
print("Lifting");
my(J=J,l=l);WB = parapply(W->PicLiftTors(J,W,1,l),WB);
print("All the space");
V = TorsSpace(J,WB,l);
print("Evaluation");
my(J=J);Z=parapply(W->HyperPicEval(J,W)[1],V[1..l^d-1]);
F=factorback(apply(u->'x-u,Mod(Z,Jgetpe(J))*Mod(1,JgetT(J))));
F=liftall(F)*(1+O(p^e));
bestappr(F);
