read("install.gp");
read("galrep.gp");
\\ Xns+(13)
f = (-y-z)*x^3+(2*y^2+y*z)*x^2+(-y^3+z*y^2-2*z^2*y+z^3)*x+2*z^2*y^2-3*z^3*y;
f = subst(f,z,x+z);
f = subst(f,z,y+z);
f = subst(f,z,1);
f = -6*x^4 + (-24*y - 37)*x^3 + (-20*y^2 - 113*y - 84)*x^2 + (-7*y^3 - 67*y^2 - 178*y - 83)*x + (-y^4 - 13*y^3 - 57*y^2 - 94*y - 30);
\\ Note : End(J) = O_K, K = Q(zeta7)+
\\ Rat pts and hyper section
P1=[-1,0];
P2=[-1,-1];

/*f = subst(f,z,x+z);
f = subst(f,z,y+z);
f = subst(f,z,1);
f = subst(f,x,x+1);*/
\\p=43;a=4;e=512;l=3;
p=53;a=7;e=256;l=13;
C=x^2+3*x+1;d=2;
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
U=PlaneEval0_data(J,[1,0,0],[P1],[1,0,0],[P2]);
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
