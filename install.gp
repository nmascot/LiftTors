install("matkerpadic","GGGL",,"./liblinalg.so");
install("mateqnpadic","GGGL",,"./liblinalg.so");
install("VecSmallCompl","GU","VecSmallCompl","./liblinalg.so");
install("FqM_MinorCompl","GGG","FqM_MinorCompl","./liblinalg.so");
install("matF","GGGU","matF","./liblinalg.so");

install("NAF","G","NAF","./libexp.so");
install("AddChain","GL","AddChain","./libexp.so");

install("HyperInit","GGUL","HyperInit","./libhyper.so");
install("HyperRandPt","GGGUG","HyperRandPt","./libhyper.so");
install("ordJ","GGU","ordJ","./libhyper.so");
install("HyperPicRand","GG","HyperPicRand","./libhyper.so");
install("HyperPicRandTors","GGGG","HyperPicRandTors","./libhyper.so");
install("HyperPicRandDbg","GG","HyperPicRandDbg","./libhyper.so");

install("PicRed","GU","PicRed","./libpic.so");
install("PicChord","GGGL","PicChord","./libpic.so");
install("PicAdd","GGG","PicAdd","./libpic.so");
install("PicSub","GGG","PicSub","./libpic.so");
install("PicNeg","GG","PicNeg","./libpic.so");
install("PicMul","GGGL","PicMul","./libpic.so");
install("PicFrob","GG","PicFrob","./libpic.so");
install("PicFrobPoly","GGG","PicFrobPoly","./libpic.so");
install("PicEq","lGGG","PicEq","./libpic.so");
install("PicIsZero","lGG","PicIsZero","./libpic.so");
install("PicChart","GG","PicChart","./libpic.so");
install("PicRand0","G","PicRand0","./libpic.so");

install("PicLiftTors_2","GGLG","PicLiftTors_2","./liblift.so");
install("PicLiftTors","GGLG","PicLiftTors","./liblift.so");

install("PicTorsRels","GGGU","PicTorsRels","./libfreyruck.so");
install("PicNorm","GGG","PicNorm","./libfreyruck.so");

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
