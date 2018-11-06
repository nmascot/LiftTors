\\ Load the package modules
read("install.gp");
read("galrep.gp");
read("Hyper2RR.gp");

\\ Equation for the curve
f = x^6 + 4*x^5 + 6*x^4 + 2*x^3 + x^2 + 2*x + 1;
print("C : y² = ",f);
l = 7; \\ The representation occurs in the 7-torsion of the Jacobian
p = 17; \\ We choose to get it 17-adically
e = 32; \\ Target p-adic accuracy is O(17^32)
C = x^2-x-2; \\ Char.poly. of the Frobenius at p (another possible choice is x^2-2*x-1)
d = poldegree(C);
a = mordroot(C,l); \\ Degree of the unramified extension of Qp over which the torsion points are defined
print("T = part of J[",l,"] where Frob_",p," acts by ",C);
Lp = hyperellcharpoly(Mod(f,p)); \\ Local L factor of the curve at p, needed to know the number of points on the Jacobian mod p
P1=[-1,1];P2=[0,1]; \\ We need TODO
[f,g,d0,L,LL,L1,L2]=Hyper2RR(f,P1,P2); \\ Precomputation of a generic-enough equation of the curve and of some Ri
J=PicInit(f,g,d0,L,LL,1,p,a,e);
U=PicEvalInit(J,[L1,L2]);
J1 = PicRed(J,1); \\ Reduction mod p

print("\n--> Getting a basis of T mod ",p);
[B,matFrob] = TorsBasis(J1,l,Lp,C);
[WB,cWB] = TorsSpaceFrobGen(J1,l,B,matFrob); \\ Generating set of T under Frob and coordinates of these generators on B
print("\n--> Lifting ",#WB," points ",p,"-adically");
{if(#WB > Jgetg(J),
  my(J=J,l=l); WB = parapply(W->PicLiftTors(J,W,1,l),WB); \\ More efficient in parallel
,
  WB = apply(W->PicLiftTors(J,W,1,l),WB); \\ Less efficient in parallel (TODO tune)
);}
print("\n--> All of T");
TI = TorsSpaceFrob(J,WB,cWB,l,matFrob);
print("\n--> Evaluation of ",#TI[2]," points");
Z = TorsSpaceFrobEval(J,TI,U,l,d,matFrob);
print("\n--> Expansion and identification");
AF = TorsSpaceGetPols(J,Z); \\ List of polynomials defining the representation
F=AF[1][3]; \\ Simplest polynomial
print(F);
\\ Now, for the projective version of the representation:
Z=AF[1][1]; \\ Roots of F
PZ = A2P1(Z,l,1,JgetT(J),p^e); \\ Gather the roots of F along the fibers of A2 -> P1
G = factorback(apply(u->'x-u,Mod(PZ,Jgetpe(J))*Mod(1,JgetT(J)))); \\ Get corresponding polynomial mod p^e
G = liftpol(G);
GG = bestappr(G) \\ Identify it
