read("install.gp");
read("MakTorsSpace.gp");

WeiRed(f,h)=
{
 my(F=f+(h/2)^2);
 2^poldegree(F)*subst(F,x,x/2);
}


f=WeiRed(x^6-3*x^5+2*x^4+x^3-x,1)+2;
p=11;l=3;a=6;e=10;d=2;

J=HyperInit(f,p,a,1);
W1 = HyperPicRandTors(J,f,l,0);
W2 = HyperPicRandTors(J,f,l,0);
PicTorsRels(J,[W1,W2],l,1)
W2 = PicChord(J,W1,W1,1);
PicTorsRels(J,[W1,W2],l,1)

