\\ Load the package modules
read("TorsHensel.gp");

\\ Equation for the curve: y²+h(x)*y=f(x)
f = x^5+x^4;
h = x^3+x+1;
l = 7; \\ The representation occurs in the 7-torsion of the Jacobian
p = 19; \\ We choose to get it 19-adically
e = 128; \\ Target p-adic accuracy is O(19^128)
chi = x^2+5*x+4; \\ Char.poly. of the Frobenius at p
P1 = [1,Mod(w-2,w^2-w-4)]; P2 = [0,0]; \\ We need two non-conjugate points
[F,ZF] = HyperGalRep([f,h],l,p,e,P1,P2,chi);
print(F);
[G,ZG,VG] = ProjPol(ZF,l,2,0);
print(G);
print(polredabs(G));
