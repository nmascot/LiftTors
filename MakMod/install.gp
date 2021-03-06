default(path,".:..:~:~/gp");

read("../install.gp");
install("ZpXQ_FrobMat","GGLG",,"libpic.so");

install("ZNnorm","uLU","ZNnorm","libzn.so");
install("ZNneg","uLU","ZNneg","libzn.so");
install("GetCoef","GG","GetCoef","libzn.so");
install("l1","GGGGGGL","l1","libmodjac.so");
install("M2_worker","GGGGGG",,"libmodjac.so");
install("M2mat","GGGGGG",,"libmodjac.so");
install("PicTp","GG",,"libmodjac.so");
install("elladd_padic","GGGGGGL","elladd_padic","libelladd_padic.so");
install("TrE2qexp","GUGGUGUGGGL","TrE2qexp","libqexp.so");

read("ModGalRep.gp");
