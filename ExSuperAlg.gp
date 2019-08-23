\\ Load the package modules
read("TorsHensel.gp");

\\ Equation for the curve: y^m=f(x)
f = x^3+x+1; m=4;
\\ Its Jacobian contains the elliptic curve y^2=f(x).
l = 5; \\ The representation occurs in the 2-torsion of the Jacobian
p = 13; \\ We choose to get it 13-adically, because the image of Frob_p has low order for this p
e = 16; \\ Target p-adic accuracy is O(13^16)
chi = x^2+4*x+p; \\ Char.poly. of the Frobenius at p, used here to pick the piece of the 5-torsion coming from the elliptic curve quotient y^2=f(x) of this curve
P = [Mod(w,w^2+1),1]; \\ We can also use algebraic points, although this is much less efficient
[F,ZF] = SuperGalRep(f,m,l,p,e,P,chi);
F
