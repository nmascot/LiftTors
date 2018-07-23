read("install.gp");
read("MakTorsSpace.gp");
p=7;a=8;e=1000;
f=x^8 + x + 3;
p=13;a=3;e=16;l=37;
\\f=x^6+x+1;
J=HyperInit(f,p,a,e);
J1=PicRed(J,1);

W = HyperPicRandTors(J1,f,l,0);
print("Lifting");
W=PicLiftTors(J,W,1,l);
