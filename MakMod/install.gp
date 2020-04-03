install("ZpXQ_inv","GGGL");
install("ZpXQ_div","GGGGGL");
install("ZpXQM_inv","GGGL");

WeiRed(f,h)=
{
 my(F=f+(h/2)^2);
 2^poldegree(F)*subst(F,x,x/2);
}

read("../qMak/install.gp");
install("ZpXQ_FrobMat","GGLG",,"./libpic.so");
install("matkerpadic_safe","GGGL",,"./liblinalg.so");

install("ZNnorm","uLU","ZNnorm","./libzn.so");
install("ZNneg","uLU","ZNneg","./libzn.so");
install("GetCoef","GG","GetCoef","./libzn.so");
install("l1","GGGGGGL","l1","./libmodjac.so");
install("elladd_padic","GGGGGGL","elladd_padic","./libelladd_padic.so");
install("E1qexp","GUGUGGGL","CE1qexp","./libqexp.so");
install("TrE2qexp","GUGGUGUGGGL","TrE2qexp","./libqexp.so");

read("ModGalRep.gp");
