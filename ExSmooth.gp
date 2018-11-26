\\ Load the package modules
read("TorsHensel.gp");

f = x^3*y+y^3+x; \\ Equation for the curve (must be smooth)
\\ To evaluate points on the Jacobian, we need two divisors of degree d-g (=1 in this case)
\\ Here we use the rational points (0:0:1) and (1:0:0)
P1 = [0,0,1];
P2 = [1,0,0];

l = 2; \\ The representation we want is in the 2-torsion of the Jacobian
p = 5; \\ We choose to work 5-adically
e = 64; \\ to accuracy O(5^64)
\\ If instead we want the 3-torsion, we can take:
\\l=3;p=43;e=512;
chi = 0; \\ We want all the l-torsion, not a subspace.

[F,ZF] = SmoothGalRep(f,l,p,e,[P1],[P2],chi);
print(F);
factor(F)
