\\ Load the package modules
read("install.gp");
read("galrep.gp");

\\ Equation for the curve
f = x^6 + 4*x^5 + 6*x^4 + 2*x^3 + x^2 + 2*x + 1;
print("C : y² = ",f);
l = 7; \\ The representation occurs in the 7-torsion of the Jacobian
p = 17; \\ We choose to get it 17-adically
e = 128; \\ Target p-adic accuracy is O(17^128)
C = x^2-x-2; \\ Char.poly. of the Frobenius at p (another possible choice is x^2-2*x-1)
d = poldegree(C);
a = mordroot(C,l); \\ Degree of the unramified extension of Qp over which the torsion points are defined
print("T = part of J[",l,"] where Frob_",p," acts by ",C);
J = HyperInit(f,p,a,e); \\ Jacobian with accuracy O(p^e)
J1 = PicRed(J,1); \\ Reduction mod p
Lp = hyperellcharpoly(Mod(f,p)); \\ Local L factor of the curve at p, needed to know the number of points on the Jacobian mod p

print("\n--> Getting a basis of T mod ",p);
WB = TorsBasis(J1,f,p,a,l,Lp,C);
print("\n--> Lifting this basis ",p,"-adically");
WB = apply(W->PicLiftTors(J,W,1,l),WB);
print("\n--> All of T");
T = TorsSpace(J,WB,l);
print("--> Evaluation");
U = HyperPicEvalData(J); \\ Precomputation for the evaluation of a rational map J -> A1
my(J=J,U=U); Z = parapply(W->HyperPicEval(J,W,U)[1],T[1..l^d-1]); \\ Evaluation of the torsion points in parallel
F = factorback(apply(u->'x-u,Mod(Z,Jgetpe(J))*Mod(1,JgetT(J)))); \\ Corresponding polynomial mod p^e
F = liftpol(F);
FF = bestappr(F) \\ Identify the coefficients as rationals

\\ Now, for the projective version of the representation:
PZ = A2P1(Z,l,1,JgetT(J),p^e); \\ Gather the roots of F along the fibers of A2 -> P1
G = factorback(apply(u->'x-u,Mod(PZ,Jgetpe(J))*Mod(1,JgetT(J)))); \\ Get corresponding polynomial mod p^e
G = liftpol(G);
GG = bestappr(G) \\ Identify it
