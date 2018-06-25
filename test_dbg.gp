p=7;f=x^6-2*x+3;e=1;a=4;
f=x^8 + x + 3;
J=HyperInit(f,p,a,e);
W=HyperPicRand(J,f);
T=J[3];
t=variable(T);
nZ=#J[11];
Fqred(x)=Mod(x,T)*Mod(1,p);
V=J[8];
KV=J[9];
J1 = PicRed(J,1);
W1=HyperPicRandTors(J1,f,5,0);
W2=PicMul(J1,W1,2,0);

