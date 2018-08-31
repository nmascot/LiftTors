read("install.gp");
read("galrep.gp");

f=WeiRed(-2*x^5 + 2*x^4 - 5*x^3 + 2*x^2 - 2*x,-x^3 - 1);
f=4*subst(f,x,2*x)/2^6;
print("C : y² = ",f);
p=17;a=6;e=256;d=4;l=2;
\\C = x^2-2*x-1;
C = 0;
print("T : part of J[",l,"] where Frob_",p," acts by ",C);
J=HyperInit(f,p,a,e);
J1=PicRed(J,1);

print("\n--> Getting a basis of T mod ",p);
WB = TorsBasis(J1,f,p,a,l,hyperellcharpoly(Mod(f,p)),C);
print("\n--> Lifting ",p,"-adically");
my(J=J,l=l);WB = /*par*/apply(W->PicLiftTors(J,W,1,l),WB);
print("\n--> All of T");
V = TorsSpace(J,WB,l);
print("--> Evaluation");
my(J=J,U=HyperPicEvalData(J));Z=parapply(W->HyperPicEval(J,W,U)[2],V[1..l^d-1]);
F=factorback(apply(u->'x-u,Mod(Z,Jgetpe(J))*Mod(1,JgetT(J))));
F=liftpol(F);
FF=bestappr(F)
PZ=A2P1(Z,l,1,JgetT(J),p^e);
G=factorback(apply(u->'x-u,Mod(PZ,Jgetpe(J))*Mod(1,JgetT(J))));
G=liftpol(G);
GG=bestappr(G)
