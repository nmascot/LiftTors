read("install.gp");
read("galrep.gp");
f=x^3*y+y^3+x;
f = subst(f,y,x+y);
f = subst(f,x,x+y+1);
p=5;a=6;e=256;l=3;
C = 0; d=6;
J=PlaneInit(f,p,a,e);
J1=PicRed(J,1);

WB = TorsBasis(J1,f,p,a,l,PlaneZeta(f,p),C);
print("Lifting");
my(J=J,l=l);WB = parapply(W->PicLiftTors(J,W,1,l),WB);
print("All the space");
V = TorsSpace(J,WB,l);
print("Evaluation");
my(J=J);Z=parapply(W->PlaneEval(J,W,[1,2,3],[[-1,0]],[1,-1,-1],[2,1,-1]),V[1..l^d-1]);
F=factorback(apply(u->'x-u,Mod(Z,Jgetpe(J))*Mod(1,JgetT(J))));
F=liftpol(F);
bestappr(F)
