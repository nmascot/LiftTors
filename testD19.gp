read("install.gp");
read("../qMak/GalRep.gp");

FieldOfDef(J,W)=
{
	for(n=1,10^10,
		if(PicEq(J,W,PicFrobPoly(J,W,'x^n)),return(n))
	);
}

time = getwalltime();

N=l=19;
p=107;
a=6;
chi=x^2+10*x+8;
Lp = LMod(N,1,p);
if(poldegree(gcd(Mod(Lp,l)/Mod(chi,l),Mod(chi,l))),error("Chi not coprime with its cofactor"));
e=256;

[J,M4Q,CuspsQ]=ModJacInit(N,1,p,a,e);
J1 = PicRed(J,1);

default(debug,1);

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
