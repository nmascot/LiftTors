\\ Load the package modules
read("install.gp");
read("galrep.gp");
read("Smooth2RR.gp");

f = x^3*y+y^3+x; \\ Equation for the curve (must be smooth)
\\ To evaluate points on the Jacobian, we need two divisors of degree d-g (=1 in this case)
\\ Here we use the rational points (0:0:1) and (1:0:0)
P1 = [0,0,1];
P2 = [1,0,0];
\\ Quick intial precomputation: make the equation ufficiently generic, compute genus g, and some data that will be needed later
[f,g,d0,L,LL,L1,L2] = Smooth2RR(f,[P1],[P2]);

l = 2; \\ The representation we want is in the 2-torsion of the Jacobian
p = 5; \\ We choose to work 5-adically
e = 64; \\ to accuracy O(5^64)
Lp = PlaneZeta(f,p); \\ Local L factor at p
a = 6; \\ Degree of the unramified extension of Qp over which the torsion points are defined
\\ If instead we want the 3-torsion, we can take:
\\l=3;p=43;e=512;a=4;
C = 0; \\ We want all the l-torsion, not a subspace.
d = 6; \\ Dimension of the representation

J=RRInit2(f,g,d0,L,LL,1,p,a,e); \\ Jacobian with accuracy O(p^e)
J1 = PicRed(J,1); \\ Reduction mod p
U=RREvalInit(J,[L1,L2]); \\ Data to evaluate points on the Jacobian

print("\n--> Getting basis of T mod ",p)
[B,matFrob] = TorsBasis(J1,l,Lp,C);
[F,Q] = matfrobenius(Mod(matFrob,l),2); \\ F = rational canonical form of matFrob, Q = transition matrix
Q = lift(Q^(-1));
NF = apply(poldegree,matfrobenius(F,1)); \\ List of degrees of elementary divisors
WB = List(); cWB = List();
n = 1;
{for(i=1,#NF,
	c=Q[,n];
	listput(cWB,c);
	c=centerlift(Mod(c,l));
	listput(WB,PicLC(J,c,B));
	n+=NF[i]
);}
WB = Vec(WB); \\ Generating set of T under Frob
cWB = Vec(cWB); \\ Coordinates of these generators on B
print("\n--> Lifting ",#WB," points ",p,"-adically");
{if(#WB > Jgetg(J),
	my(J=J,l=l); WB = parapply(W->PicLiftTors(J,W,1,l),WB); \\ More efficient in parallel
,
	WB = apply(W->PicLiftTors(J,W,1,l),WB); \\ Less efficient in parallel (TODO tune)
);}
print("\n--> All of T");
TI = TorsSpaceFrob(J,WB,cWB,l,matFrob);
print("\n--> Evaluation of ",#TI[2]," points");
Z = TorsSpaceFrobEval(J,TI,U,l,d);
print("\n--> Expansion and identification");
A = AllPols(Z,JgetT(J),p,e,Jgetpe(J)); \\ p-adic approximation of a set of polynomials which all define a subfield of the field cut out by the representation, with equality iff. no repeated roots
print(#A," candidate polynomials");
AI = select(x->x[3]!=[],A); \\ Drop the approximations that could not be identified as rationals
print(#AI," identified polynomials");
AS = select(x->#Set(x[1])==#(x[1]),AI); \\ Drop the polynomials having multiple roots
print(#AS," faithful polynomials");
AS = vecsort(AS,x->sizebyte(x[3])); \\ Take the simplest of all polynomials
F=AS[1][3]
factor(F)
