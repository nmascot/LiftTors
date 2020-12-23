read("install.gp");

S=mfinit([16,2,Mod(5,16)],0); \\ New cuspidal space of level 16, weight 2, nebentypus Conrey index 5
f=mfeigenbasis(S)[1]; \\ The unique eigenform (dim S = 1)
l=5;
coeffs=[[3,1]];
rangep=[100,105];
qprec=3;

X=mfgalrep(f,l,coeffs,[58,60],10,qprec,1);
print(X[1]);
breakpoint();
print("--> Finding prime p with small residual degree");
[H,best] = mfbestp(f,l,coeffs,rangep);
[p,a,Lp,chi] = best;
print("Chosen p=",p,", residual degree ",a);
e=8;
pe = p^e;
[N,k] = mfparams(f)[1..2];
if(k>2,N*=l);
print("\n--> Initialising modular Jacobian");
J=ModJacInit(N,H,p,a,e,qprec,Lp);
print("Size J: ",mysize(sizebyte(J)));
J1 = PicRed(J,1);

NJ=polresultant(Lp,'x^a-1);
M=NJ/l^valuation(NJ,l);
export(J1);
export(M);
z = Fq_zeta_l(JgetT(J),Jgetp(J),l);
AddC = AddChain(l,0);
W0 = JgetW0(J);
W0 = PicChord(J,W0,W0,1);
FRparams = [AddC,W0,z];

RandTors(J,l,M)=
{
	my(W,lW);
	W = PicRand(J);
	W = PicMul(J,W,M,0);
	lW = PicMul(J,W,l,0);
	while(PicIsZero_val(J,lW)==0,
			print("m");
		W = lW;
		lW = PicMul(J,W,l,0);
	);
	print("M");
	W;
}

R=matrix(4,4);
{
while(matdet(Mod(R,l))==0,
		print("a");
B=vector(4,i,RandTors(J1,l,M));print("b");
Tests=vector(4,i,PicRand(J1));
R=matconcat(apply(W->Tors_TestPt(J1,W,l,Tests,FRparams),B));print("c");
printp(R);
print(matdet(R));
);
}
print("d");
TpB=apply(W->PicTp(J1,W),B);
TpR=matconcat(apply(W->Tors_TestPt(J1,W,l,Tests,FRparams),TpB));

S=liftint(Mod(R^(-1)*TpR,l))

\\[B,matFrob] = TorsBasis(J1,l,Lp,chi); \\ Basis of the mod p^1 space and matrix of Frob_p

\\W=PicRand(J);
