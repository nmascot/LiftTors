read("install.gp");
read("galrep.gp");
\\f = x^4 - 5*x^3 + 9*x^2 + (-y - 7)*x + (2*y^4 + y^2 + 2);
f = 8*x^4 + (4*y - 32)*x^3 + (6*y^2 - 12*y + 48)*x^2 + (-12*y^2 + 12*y - 104)*x + (y^4 + y^3 + 6*y^2 - 4*y + 80);
p=23;a=6;e=128;l=2;
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
U=PlaneEval0_data(J,[1,0,0],[[1,0]],[1,0,0],[[1,-1]]);
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
