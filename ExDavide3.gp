\\ Load the package modules
read("TorsHensel.gp");

f = 4*x^4+4*x^3+6*x-6*y^4-9;  \\ Equation for the curve (must be smooth)
\\ To evaluate points on the Jacobian, we need two divisors of degree d-g (=1 in this case)
\\ Here we use these points, whcih are define over different cubic fields:
P1 = [Mod(w,w^2+6),0,2];
P2 = [Mod(w-1,w^2-7),0,2];

l = 2; \\ The representation we want is in the 2-torsion of the Jacobian
p = 11; \\ We choose to work 11-adically
e = 1024; \\ to accuracy O(11^1024)
chi = 0; \\ We want all the l-torsion, not a subspace.

[F,ZF] = SmoothGalRep(f,l,p,e,[P1],[P2],chi,2);
print(F);
print("Factoring the polynomial");
fa=factor(F);
printp(fa);
print("Reducing the factors");
fa[,1]=parapply(polredbest,fa[,1]);
fa
