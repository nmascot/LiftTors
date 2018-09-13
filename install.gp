WeiRed(f,h)=
{
 my(F=f+(h/2)^2);
 2^poldegree(F)*subst(F,x,x/2);
}

install("matkerpadic","GGGL",,"./liblinalg.so");
install("mateqnpadic","GGGL",,"./liblinalg.so");
install("VecSmallCompl","GU","VecSmallCompl","./liblinalg.so");
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
\\install("RReval","GUUGG","RReval","./libhyper.so");

install("PlaneInit","GGUL","PlaneInit","./libplanereg.so");
install("PlaneZeta","GU","PlaneZeta","./libplanereg.so");
install("PlaneEval","GGGGGG","PlaneEval","./libplanereg.so");
install("PlaneEval0_data","GGGGG","PlaneEval0_data","./libplanereg.so");
install("PlaneEval0","GGG","PlaneEval0","./libplanereg.so");
install("AllPols0","GGGLG","AllPols0","./libplanereg.so");


install("RRInit","GUUGGGUL","RRInit","./librr.so");
\\install("FnEvalAt","GGGGGGG","FnEvalAt","./librr.so");
\\install("CurveRandPt","GGGLG","CurveRandPt","./librr.so");
install("RREvalInit","GG","RREvalInit","./librr.so");
install("RREval","GGG","RREval","./librr.so");

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
install("PicChart","GGU","PicChart","./libpic.so");
install("PicRand0","G","PicRand0","./libpic.so");

install("JgetW0","G","JgetW0","./libpic.so");
install("Jgetg","lG","Jgetg","./libpic.so");
install("JgetT","G","JgetT","./libpic.so");
install("Jgetpe","G","Jgetpe","./libpic.so");
install("Jgetp","G","Jgetp","./libpic.so");
install("Jgete","lG","Jgete","./libpic.so");

install("PicLift_worker","UUUGGGGGGG","PicLift_worker","./liblift.so");
install("PicLiftTors_worker","GGGGGGUUUUGGLGLGGGU","PicLiftTors_worker","./liblift.so");
install("PicLiftTors","GGLG","PicLiftTors","./liblift.so");

install("PicNorm","GGG","PicNorm","./libfreyruck.so");
install("PicFreyRuckMulti","GGGGGG","PicFreyRuckMulti","./libfreyruck.so");
install("PicTorsRels","GGGU","PicTorsRels","./libfreyruck.so");
install("Fq_zeta_l","GGG","Fq_zeta_l","./libfreyruck.so");
install("Fq_mu_l_log","GGGGG","Fq_mu_l_log","./libfreyruck.so");

/*f=x^3*y+y^3+x;
f = subst(f,y,x+y);
f = subst(f,x,x+y+1);
p=5;a=6;e=4;l=3;
J=PlaneInit(f,p,a,e);
W = PicRand0(J);
PlaneEval(J,W,[1,2,3],[[-1,0]],[1,-1,-1],[2,1,-1])*/

