read("install.gp");
read("galrep.gp");
f=x^3*y+y^3+x;
f = subst(f,y,x+y);
f = subst(f,x,x+y+1);
p=5;a=6;e=1;l=2;
C = 0; d=6;
J1=PlaneInit(f,p,a,e);

WB = TorsBasis(J1,f,p,a,l,PlaneZeta(f,p),C);
