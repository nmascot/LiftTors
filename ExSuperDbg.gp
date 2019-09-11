/* This example triggers several bugs (too many symmetries?) */

\\ Load the package modules
read("TorsHensel.gp");

\\ Equation for the curve: y^m=f(x)
f = x^4+1; m=3;
l = 2; \\ The representation occurs in the 2-torsion of the Jacobian
p = 11; \\ We choose to get it 11-adically
e = 64; \\ Target p-adic accuracy is O(11^64)
chi = 0; \\ Char.poly. of the Frobenius at p (Use to pick only a piece of the whole 2-torsion; here 0 means we want all the 2-torsion.)
P = [0,1]; \\ We need a rational point
[F,ZF] = SuperGalRep(f,m,l,p,e,P,chi);
print(F);
factor(F)
