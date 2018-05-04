read("Frob.gp");
read("LinAlg.gp");
read("MakHyper.gp");
read("MakLift.gp");
read("galrep.gp");
read("Weil.gp");
f=WeiRed(x^6-3*x^5+2*x^4+x^3-x,1);
\\ 85805.a.85805.1
\\ Split: 131.a.1=q-q^3-2q^4-2q^5-q^7 et 665.a.1=q-2q^2-3q^3+2q^4-q^5+6q^6-3q^7
p=101;
a=2;
NJp=ordJ(f,p,a);
X=PicInit(f,p,a,1);
[f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X;
A=vector(6,i,RandPt(f,p,T,1));
B=vector(6,i,RandPt(f,p,T,1));
WA=V*matker(RReval(A,3/2*d0,df));
WB=V*matker(RReval(B,3/2*d0,df));
E=vector(6,i,RandPt(f,p,T,1));
WE=V*matker(RReval(E,3/2*d0,df));
[WC,sABC]=PicChord(X,WA,WB,1);
PicEvalNorm(X,sABC,WE)
WE2=PicMul(X,WE,29);
PicEvalNorm(X,sABC,WE2)

/*
for(j=1,11,for(i=1,6,print(liftall(makfneval(X,WE[,j],E[i])))))
for(i=1,6,print(liftall(makfneval(X,sABC,E[i]))))
*/
