read("install.gp");
N=13;p=29;a=3;e=4;
J=ModJacInit(N,1,p,a,e);
T=JgetT(J);t=variable(T);Red(x)=Mod(Mod(x,T),p);
Lp = LMod(N,1,p);
J1 = PicRed(J,1);
NJ1 = polresultant(Lp,'x^a-1);
W0 = JgetW0(J1);
