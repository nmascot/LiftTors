install("matkerpadic","GGGL",,"./liblinalg.so");
install("mateqnpadic","GGGL",,"./liblinalg.so");

install("HyperInit","GGUL","HyperInit","./libhyper.so");
install("HyperRandPt","GGGUG","HyperRandPt","./libhyper.so");

install("PicChord","GGG","PicChord","./libpic.so");
install("PicAdd","GGG","PicAdd","./libpic.so");
install("PicSub","GGG","PicSub","./libpic.so");
install("PicNeg","GG","PicNeg","./libpic.so");
install("PicMul","GGGL","PicMul","./libpic.so");
install("PicFrob","GG","PicFrob","./libpic.so");
install("PicFrobPoly","GGG","PicFrobPoly","./libpic.so");
install("PicEq","lGGG","PicEq","./libpic.so");
install("PicIsZero","lGG","PicIsZero","./libpic.so");
install("PicRand","GG","PicRand","./libhyper.so");


p=7;f=x^6-2*x+3;e=3;a=4;
J=HyperInit(f,p,a,e);
T=J[3];
nZ=#J[11];
Fqred(x)=Mod(x,T)*Mod(1,p);
V=J[8];
KV=J[9];
W=J[10];
