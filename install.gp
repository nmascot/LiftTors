default(sopath,".:..");

WeiRed(f,h)=
{
 my(F=f+(2)^2);
 2^poldegree(F)*subst(F,x,2);
}

install("PicDeflate_U","GGL",,"liblift.so");
install("PicInflate_U","GG",,"liblift.so");
install("PicMember","lGG",,"libpic.so");


install("matkerpadic","GGGGL",,"liblinalg.so");
install("mateqnpadic","GGGGL",,"liblinalg.so");
install("VecSmallCompl","GU","VecSmallCompl","liblinalg.so");
install("matF","GGGU","matF","liblinalg.so");

install("NAF","G","NAF","libexp.so");
install("AddChain","GL","AddChain","libexp.so");

install("PlaneZeta","GU","PlaneZeta","libzeta.so");
install("SuperZeta","GUU","SuperZeta","libzeta.so");

install("RRInit","GUUGGGUL","PicInit","librr.so");
install("RREval","GG","PicEval","librr.so");
install("RREval_worker","GG",,"librr.so");
install("Jlift","GU",,"librr.so");

install("PicRed","GU","PicRed","libpic.so");
install("PicChord","GGGL","PicChord","libpic.so");
install("PicAdd","GGG","PicAdd","libpic.so");
install("PicSub","GGG","PicSub","libpic.so");
install("PicNeg","GG","PicNeg","libpic.so");
install("PicMul","GGGL","PicMul","libpic.so");
install("PicFrob","GG","PicFrob","libpic.so");
install("PicFrobPoly","GGG","PicFrobPoly","libpic.so");
install("PicEq","lGGG","PicEq","libpic.so");
install("PicIsZero","lGG","PicIsZero","libpic.so");
install("PicChart","GGU","PicChart","libpic.so");
install("PicRand0","G","PicRand","libpic.so");

install("JgetW0","G","JgetW0","libpic.so");
install("Jgetg","lG","Jgetg","libpic.so");
install("JgetT","G","JgetT","libpic.so");
install("Jgetpe","G","Jgetpe","libpic.so");
install("Jgetp","G","Jgetp","libpic.so");
install("Jgete","lG","Jgete","libpic.so");
install("JgetFrobMat","G","JgetFrobMat","libpic.so");
install("Frob","GGGG","Frob","libpic.so");

install("PicLift_worker","GUGGGGG",,"liblift.so");
install("PicLiftTors_Chart_worker","GGGGGGGGGLGUG",,"liblift.so");
install("PicLiftTors","GGLG",,"liblift.so");

install("PicNorm","GGGGU","PicNorm","libfreyruck.so");
install("PicFreyRuckMulti","GGGGGG","PicFreyRuckMulti","libfreyruck.so");
install("PicTorsRels","GGGU","PicTorsRels","libfreyruck.so");
install("Fq_zeta_l","GGG","Fq_zeta_l","libfreyruck.so");
install("Fq_mu_l_log","GGGGG","Fq_mu_l_log","libfreyruck.so");

install("TorsSpaceFrob_worker","GGGGG",,"libtorsspace.so");
install("TorsSpaceFrobEval","GGGUG",,"libtorsspace.so");
install("OnePol","GGGGUGGG",,"libtorsspace.so");
install("AllPols","GUGGGGGL",,"libtorsspace.so");
install("c2i","uGU",,"libtorsspace.so");
install("i2c","UUU",,"libtorsspace.so");
install("Chordi","uUUUU",,"libtorsspace.so");
install("ActOni","uGUU",,"libtorsspace.so");
