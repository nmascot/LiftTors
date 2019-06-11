\\ Load the package modules
read("TorsHensel.gp");

f = x^3*y+y^3+x; \\ Equation for the curve (must be smooth)
\\ To evaluate points on the Jacobian, we need two divisors of degree d-g (=1 in this case)
\\ Here we use the rational points (0:0:1) and (1:0:0)
P1 = [0,0,1];
P2 = [1,0,0];

l = 3; \\ The representation we want is in the 3-torsion of the Jacobian
p = 43; \\ We choose to work 43-adically
e = 256; \\ to accuracy O(43^256)
\\ Actually this accuracy is too low, but the algorithm will automatically crank it up to e=512 once it realises this.
\\ Setting e=512 in the first place makes the computation faster. Caveat: unnecessarily large values of e slow the computation down!
chi = 0; \\ We want all the l-torsion, not a subspace.
a = 4; \\ The local L factor of the curve indicates that the 3-torsion is defined over the unramified extension of degree 12 of Q_43, so normally the algorithm would automatically pick this degree.
\\ Here we help it by insisting to work in the extension of degree 4 instead, which is faster (and the 3-torsion is indeed actually defined there).
\\ Caveat: Choosing a degree which is too small will trap the algorithm in an infite loop searching for torsion points that do not exist in that degree.
\\ It is usually better NOT to pass this argument (which amounts to passing a=0), in which case the algorithm automatically chooses an appropriate value, which is usually sharp.
[F,ZF] = SmoothGalRep(f,l,p,e,[P1],[P2],chi,a);
faF=factor(F)
