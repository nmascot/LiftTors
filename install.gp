WeiRed(f,h)=
{
 my(F=f+(h/2)^2);
 2^poldegree(F)*subst(F,x,x/2);
}

install("matkerpadic","GGGL",,"./liblinalg.so");
install("mateqnpadic","GGGL",,"./liblinalg.so");
install("VecSmallCompl","GU","VecSmallCompl","./liblinalg.so");
\\install("FqM_MinorCompl","GGG","FqM_MinorCompl","./liblinalg.so");
install("matF","GGGU","matF","./liblinalg.so");

install("NAF","G","NAF","./libexp.so");
install("AddChain","GL","AddChain","./libexp.so");

install("HyperInit","GGUL","HyperInit","./libhyper.so");
install("HyperRandPt","GGGUG","HyperRandPt","./libhyper.so");
install("ordJ","GGU","ordJ","./libhyper.so");
install("HyperPicRand","GG","HyperPicRand","./libhyper.so");
install("HyperPicRandTors","GGGG","HyperPicRandTors","./libhyper.so");
install("HyperPicEvalData","G","HyperPicEvalData","./libhyper.so");
install("HyperPicEval","GGG","HyperPicEval","./libhyper.so");
\\install("HyperPicRandDbg","GG","HyperPicRandDbg","./libhyper.so");

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

install("JgetW0","G","JgetW0","./libpic.so");
install("Jgetg","lG","Jgetg","./libpic.so");
install("JgetT","G","JgetT","./libpic.so");
install("Jgetpe","G","Jgetpe","./libpic.so");

\\install("PicLiftTors_2","GGLG","PicLiftTors_2","./liblift.so");
install("PicLift_worker","UUUGGGGGGG","PicLift_worker","./liblift.so");
install("PicLiftTors_worker","GGGGGGUUUUGGLGLGGG","PicLiftTors_worker","./liblift.so");
install("PicLiftTors","GGLG","PicLiftTors","./liblift.so");

install("PicNorm","GGG","PicNorm","./libfreyruck.so");
install("PicFreyRuckMulti","GGGGGG","PicFreyRuckMulti","./libfreyruck.so");
install("PicTorsRels","GGGU","PicTorsRels","./libfreyruck.so");
