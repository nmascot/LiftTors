read("install.gp");
read("GalRep.gp");

FieldOfDef(J,W)=
{
	for(n=1,10^10,
		if(PicEq(J,W,PicFrobPoly(J,W,'x^n)),return(n))
	);
}

time=getwalltime();

N=13;
p=73;
a=4;
l=13;
chi=x^2+7*x+5;
Lp = LMod(N,1,p);
if(poldegree(gcd(Mod(Lp,l)/Mod(chi,l),Mod(chi,l))),error("Chi not coprime with its cofactor"));
e=64;

[J,M4Q,CuspsQ]=ModJacInit(N,1,p,a,e);
J1 = PicRed(J,1);

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
print("\n--> Getting polynomials");
AF = TorsSpaceGetPols(J,Z);
print(strtime(getwalltime()-time));
