i=ffgen(Mod('i^2+1,3));
[f,Z,P3,G,C,h]=read("H2/H2final.txt");

FindFrob(p)=
{
 my(a,u,ev,M);
 a=Mod(x,f)*Mod(1,p);
 u=a^p*subst(h,'x,a);
 u=trace(u);
 ev=parapply(g->subst(g,'x,u),G);
 ev=select(y->y==0,ev,1);
 if(#ev==0,error("No possible Frobenius"));
 if(#ev>1,print("Found ",#ev," possible Frobenius");return(0));
 ev=ev[1];
 kronecker(p,3)*C[ev][1];
}
