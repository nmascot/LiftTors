\\ Load the package modules
read("TorsHensel.gp");

\\ Equation for the curve: y²+h(x)*y=f(x)
f = x^5+x^4; h = x^3+x+1; \\ X1(13)
l = 7; \\ The representation occurs in the 7-torsion of the Jacobian
p = 17; \\ We choose to get it 17-adically
e = 32; \\ Target p-adic accuracy is O(17^32)
chi = x^2-x-2; \\ Char.poly. of the Frobenius at p (another possible choice is x^2-2*x-1)
P1 = [-1,1]; P2 = [0,0]; \\ We need two non-conjugate rational points
[F,ZF] = HyperGalRep([f,h],l,p,e,P1,P2,chi);
printf("\nPolynomial defining the linear representation:\n");
print(F);
[G,ZG,VG] = ProjPol(ZF,l,2,0);
printf("\nPolynomial defining the projective representation:\n");
print(G);
printf("\nReduced version:\n");
print(polredabs(G));
