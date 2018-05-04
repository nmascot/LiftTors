read("Frob.gp");
read("LinAlg.gp");
read("MakHyper.gp");
read("MakLift.gp");
read("Weil.gp");
f=WeiRed(x^6-3*x^5+2*x^4+x^3-x,1);
\\ 85805.a.85805.1
\\ Split: 131.a.1=q-q^3-2q^4-2q^5-q^7 et 665.a.1=q-2q^2-3q^3+2q^4-q^5+6q^6-3q^7
p=7;
a=4;
NJp=ordJ(f,p,a);
X=PicInit(f,p,a,100);
X1=PicRed(X,1);
W1=PicRandTors(X1,3,'x^2+3*'x+7);
W2=PicRandTors(X1,3,'x^2+3*'x+7);
