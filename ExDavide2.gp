\\ Load the package modules
read("TorsHensel.gp");

f = x^4+y^4+6*y^2+4*y+2;  \\ Equation for the curve (must be smooth)
\\ To evaluate points on the Jacobian, we need two divisors of degree d-g (=1 in this case)
\\ Here we use these points, whcih are define over different cubic fields:
P1 = [1,Mod(w,w^4+1),0];
P2 = [Mod(w,w^4+2),0,1];

l = 2; \\ The representation we want is in the 2-torsion of the Jacobian
p = 29; \\ We choose to work 29-adically
e = 2048; \\ to accuracy O(29^2048)
chi = 0; \\ We want all the l-torsion, not a subspace.

[F,ZF] = SmoothGalRep(f,l,p,e,[P1],[P2],chi,4);
print(F);
print("Factoring the polynomial");
fa=factor(F);
printp(fa);
print("Reducing the factors");
fa[,1]=parapply(polredbest,fa[,1]);
fa
