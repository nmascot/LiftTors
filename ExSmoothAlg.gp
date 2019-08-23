\\ Load the package modules
read("TorsHensel.gp");

f = x^3*y+y^3+x; \\ Equation for the curve (must be smooth)
\\ To evaluate points on the Jacobian, we need two divisors of degree d-g (=1 in this case)
\\ Here we use these points, which are defined over different cubic fields:
P1 = [1,Mod(w,w^3+w+1),1];
P2 = [-1,Mod(w,w^3-w-1),1];

l = 2; \\ The representation we want is in the 2-torsion of the Jacobian
p = 5; \\ We choose to work 5-adically
e = 2048; \\ to accuracy O(5^2048)
chi = 0; \\ We want all the l-torsion, not a subspace.

[F,ZF] = SmoothGalRep(f,l,p,e,[P1],[P2],chi);
print(F);
print("Factoring the polynomial");
fa=factor(F);
printp(fa);
print("Reducing the factors");
fa[,1]=parapply(polredbest,fa[,1]);
fa
