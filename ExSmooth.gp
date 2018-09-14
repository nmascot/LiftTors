\\ Load the package modules
read("install.gp");
read("galrep.gp");
\\ Equation for the curve
f = x^3*y+y^3*z+z^3*x;
\\ The equation must be generic: x^d and y^d must be part of the affine equation
\\ so we change the equation a bit:
f = subst(f,z,x+z);
f = subst(f,z,y+z);
f = subst(f,z,1);
f = subst(f,x,x+1);
\\ To evaluate points on the Jacobian, we need two divisors of degree d-g (=1 in this case)
\\ whose support does not intersect the line x=0
\\ Here we use the rational points (-1,0) and (-1,-1)
P1 = [[-1,0]];
P2 = [[-1,-1]];
l = 2; \\ The representation we want is in the 2-torsion of the Jacobian
C = 0; \\ We wnat all the l-torsion, not a subspace.
d = 6; \\ Dimension of the representation
p = 29; \\ We work 29-adically
e = 32; \\ to accuracy O(29^32)
a = 3; \\ Degree of the unramified extension of Qp over which the torsion points are defined
\\ If instead we want the 3-torsion, we can take: l=3;p=43;e=512;a=4;

J = PlaneInit(f,p,a,e); \\ Jacobian with accuracy O(p^e)
J1 = PicRed(J,1); \\ Reduction mod p

print("\n--> Computing local L factor at p=",p," by point counting...");
Lp = PlaneZeta(f,p);
print("\n--> Getting basis of T mod ",p)
WB = TorsBasis(J1,l,Lp,C);
print("\n--> Lifting ",p,"-adically");
my(J=J,l=l); WB = parapply(W->PicLiftTors(J,W,1,l),WB);
print("\n--> All of T");
T = TorsSpace(J,WB,l);
print("\n--> Evaluation");
U = PlaneEval0_data(J,[1,0,0],P1,[1,0,0],P2); \\ Precomputation for the evaluation of a rational map J -> A1
my(J=J,U=U); Z = parapply(W->PlaneEval0(J,W,U),T[1..l^d-1]); \\ Evaluation of the torsion points in parallel
A = AllPols0(Z,JgetT(J),p,e,Jgetpe(J)); \\ p-adic approximation of a set of polynomials which all define a subfield of the field cut out by the representation, with equality iff. no repeated roots
print(#A," candidate polynomials");
AI = select(x->x[3]!=[],A); \\ Drop the approximations that could not be identified as rationals
print(#AI," identified polynomials");
AS = select(x->#Set(x[1])==#(x[1]),AI); \\ Drop the polynomials having multiple roots
print(#AS," faithful polynomials");
AS = vecsort(AS,x->sizebyte(x[3])); \\ Take the simplest of all polynomials
F=AS[1][3]
factor(F)
