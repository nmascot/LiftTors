install("matkerpadic","GGGL",,"./liblinalg.so");
install("mateqnpadic","GGGL",,"./liblinalg.so");

install("HyperInit","GGUL","HyperInit","./libhyper.so");
install("HyperRandPt","GGGUG","HyperRandPt","./libhyper.so");

install("PicChord","GGG","PicChord","./libpic.so");
install("PicAdd","GGG","PicAdd","./libpic.so");
install("PicSub","GGG","PicSub","./libpic.so");
install("PicNeg","GG","PicNeg","./libpic.so");
install("PicMul","GGGL","PicMul","./libpic.so");
install("DivAdd","GGUGGLGU","DivAdd","./libpic.so");
install("RandVec_padic","GGGG","RandVec_padic","./liblinalg.so");

p=7;f=x^6-2*x+3;e=1;a=6;
J=HyperInit(f,p,a,e);
T=J[3];
nZ=#J[11];
Fqred(x)=Mod(x,T)*Mod(1,p);
V=J[8];
KV=J[9];
W=J[10];
