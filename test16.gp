default(parisizemax,2G)

read("install.gp");
read("GalRep.gp");

FieldOfDef(J,W)=
{
	for(n=1,10^10,
		if(PicEq(J,W,PicFrobPoly(J,W,'x^n)),return(n))
	);
}


N=16;
p=43;
a=4;
l=5;
chi=x^2+1;
\\ x^2-1: Gal rep = Borel with eigneval 1,chi, chi of cond 80: 31->1, 21->3, 17->2, Frob_43 -> [1,0;0,-1]
\\ Other choice: x^2+1: Galrep = Borel psi chi5, chi5 cyclo mod 5, psi mod 16: -1->1, 5->2, Frob_43 -> [2,0;0,3]
e=32;

[J,M4Q,CuspsQ]=ModJacInit(N,1,p,a,e);
Lp = LMod(N,1,p);
J1 = PicRed(J,1);
NJ=polresultant(Lp,x^a-1);

/*Cl=AddChain(l,0);
W2=PicRand(J1);
W22=PicMul(J1,W2,2,1);
W23=PicMul(J1,W2,3,1);
W0=PicMul(J1,W2,NJ,0);
T=JgetT(J);
FR=PicFreyRuckMulti(J1,WT,l,[W2,W22,W23],W0,Cl);
r=(p^a-1)/l;
apply(z->Mod(Mod(z,T),p)^r,FR)*/

[B,matFrob] = TorsBasis(J1,l,Lp,chi); \\ Basis of the mod p^1 space and matrix of Frob_p
print("The matrix of Frob is");
printp(centerlift(matfrobenius(Mod(matFrob,l))));
i=1;M=Mod(matFrob,l);
while(M!=1,M*=matFrob;i++);
print("It has order ",i);

[WB,cWB] = TorsSpaceFrobGen(J1,l,B,matFrob);
print("\n--> Lifting ",#WB," points ",p,"-adically");
WB = apply(W->PicLiftTors(J,W,1,l),WB);

print("\n--> All of T");
TI = TorsSpaceFrob(J,WB,cWB,l,matFrob);
print("\n--> Evaluation of ",#TI[2]," points");
export(M4Q);
export(PicEval);
Z = TorsSpaceFrobEval(J,TI,l,2,matFrob);
AF = TorsSpaceGetPols(J,Z);
