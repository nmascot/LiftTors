read("install.gp");
read("GalRep.gp");

f = x^6 + 47*x^5 + 365*x^4 + 865*x^3 + 400*x^2 + 38*x - 4;
P0 = [0,2*iF];
A = x^2 + (4/5*iF/(u^3 + 6/5*u^2 + u)*v + (3/5*u^3 - 52/5*u^2 + 33/5*u - 8/5)/(u^3 + 6/5*u^2 + u))*x - 4/5*iF/(u^4 + 6/5*u^3 + u^2)*v + (-8/5*u^3 + 46/5*u^2 - 38/5*u + 8/5)/(u^4 + 6/5*u^3 + u^2);
B = ((-278/25*u^3 + 68/5*u^2 + 426/25*u - 8/25)/(u^6 + 12/5*u^5 + 86/25*u^4 + 12/5*u^3 + u^2)*v + (396/25*iF*u^6 + 6496/25*iF*u^5 + 24652/25*iF*u^4 + 19416/25*iF*u^3 + 2976/25*iF*u^2 + 928/25*iF*u - 16/25*iF)/(u^6 + 12/5*u^5 + 86/25*u^4 + 12/5*u^3 + u^2))*x + (318/25*u^3 - 252/25*u^2 - 338/25*u + 48/25)/(u^7 + 12/5*u^6 + 86/25*u^5 + 12/5*u^4 + u^3)*v + (-2*iF*u^7 - 596/25*iF*u^6 - 6384/25*iF*u^5 - 24316/25*iF*u^4 - 3858/5*iF*u^3 - 3176/25*iF*u^2 - 1132/25*iF*u + 96/25*iF)/(u^7 + 12/5*u^6 + 86/25*u^5 + 12/5*u^4 + u^3);

l = 11;
p = 89;
a = 3;
e = 1024;
EndoPol = 3*'x-1;

ip = sqrt(Mod(-1,p));
Ap = subst(A,iF,ip);
Bp = subst(B,iF,ip);
P0p = subst(P0,iF,ip);

P1 = [-2,Mod(10*w,w^2+10)];
P2 = [-3,Mod(10*w,w^2+10)];

X=HyperGalRep_Endo(f,l,p,a,e,[Ap,Bp],P0p,P1,P2,EndoPol);

F=X[1]; \\ Poly of degree l²-1 ncoding linear representation
G=ProjPol(X[2],l,2,1)[1]; \\ Poly of degree l+1 encoding projective representation
G=polredbest(G) \\ Nicer version (would be too slow with F)

