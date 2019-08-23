\\ Load the package modules
read("TorsHensel.gp");

\\ Equation for the curve: y^m=f(x)
f = x^4+x+1; m=3;
l = 2; \\ The representation occurs in the 7-torsion of the Jacobian
p = 11; \\ We choose to get it 17-adically
e = 128; \\ Target p-adic accuracy is O(17^32)
chi = 0; \\ Char.poly. of the Frobenius at p (another possible choice is x^2-2*x-1)
P = [0,1]; \\ We need a rational point
[F,ZF] = SuperGalRep(f,m,l,p,e,P,chi);
print(F);
