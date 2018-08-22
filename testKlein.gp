setrand(29);
read("install.gp");
read("galrep.gp");
f=x^3*y+y^3*z+z^3*x;
f = subst(f,z,x+z);
f = subst(f,z,y+z);
f = subst(f,z,1);
f = subst(f,x,x+1);
\\p=43;a=4;e=512;l=3;
p=41;a=2;e=32;l=2;
C=0;d=6;
J=PlaneInit(f,p,a,e);
J1=PicRed(J,1);

print("\n--> Computing local L factor at p=",p," by point counting...");
Lp = PlaneZeta(f,p);
print("\n--> Getting basis of T mod ",p)
WB = TorsBasis(J1,f,p,a,l,Lp,C);
print("\n--> Lifting ",p,"-adically");
my(J=J,l=l);WB = parapply(W->PicLiftTors(J,W,1,l),WB);
print("\n--> All of T");
V = TorsSpace(J,WB,l);
print("\n--> Evaluation");
U=PlaneEval0_data(J,[1,0,0],[[-1,0]],[1,0,0],[[-1,-1]]);
my(J=J,U=U);Z=parapply(W->PlaneEval0(J,W,U),V[1..l^d-1]);
A = AllPols0(Z,JgetT(J),p,e,Jgetpe(J));
\\print(#A," candidate polynomials");
AI = select(x->x[3]!=[],A);
\\print(#AI," identified polynomials");
AS = select(x->#Set(x[1])==#(x[1]),AI);
\\print(#AS," faithful polynomials");
AS = vecsort(AS,x->sizebyte(x[3]));
F=AS[1][3]
factor(F)
