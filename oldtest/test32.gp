read("Frob.gp");
read("LinAlg.gp");
read("MakHyper.gp");
read("MakLift.gp");
read("galrep.gp");
read("Weil.gp");
f=WeiRed(x^6-3*x^5+2*x^4+x^3-x,1);
\\ 85805.a.85805.1
\\ Split: 131.a.1=q-q^3-2q^4-2q^5-q^7 et 665.a.1=q-2q^2-3q^3+2q^4-q^5+6q^6-3q^7
p=29;
a=6;
l=271;
NJp=ordJ(f,p,a);
X=PicInit(f,p,a,1);
[f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X;
W1=PicRandTors(X,l);
W2=PicRandTors(X,l);
/*W=PicRand(X);
WE=PicMul(X,W,29);
[W3,F]=PicChord(X,W1,W2,1);
ZE=makfnratzeros(X,WE[,1]);
for(i=2,#WE,ZE=setintersect(Set(ZE),Set(makfnratzeros(X,WE[,i]))));
print("Found ",#ZE," rat, zeros in E");
print("Verif:");
print(matrix(#ZE,#WE,i,j,liftall(makfneval(X,WE[,j],ZE[i]))));
print("Values of F at these points:");
for(i=1,4,print(liftall(makfneval(X,F,ZE[i]))));*/
