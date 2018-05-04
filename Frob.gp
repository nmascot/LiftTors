FrobPoly(p,T,e)= \\ returns Mod(Frob(Mod(t,T)),p^e)
{
 my(t=variable(T));
 liftpol(Mod(padicappr(T,liftint(Mod(t^p,T))+O(p^e))[1],p^e));
}

Frobapply(a,T,Frob)= \\ a=Mod(A,T), Frob=Frob(Mod(t,T)), returns Frob(a)
{
 my(t=variable(T));
 Mod(subst(liftpol(a),t,Frob),T);
}

FrobApply(a,T,Frob)= \\ a=Mod(A,T), Frob=Frob(Mod(t,T)), returns Frob(a)
{
 my(t=variable(T));
 apply(u->Mod(subst(liftpol(u),t,Frob),T),a);
}
/*Frob(a,p,e)=
{
 Mod(padicappr(liftint(charpoly(a)),liftint(a^p)+O(p^e))[1],p^e);
}*/

/*FrobMat(p,T,e)=
{
 my(n=poldegree(T),t=variable(T),a,P,M);
 a=Mod(t,T);
 b=Mod(padicappr(liftint(charpoly(a)),liftint(a^p)+O(p^e))[1],p^e);
 P=M=matrix(n,n);
 for(i=1,n-1,M[i+1,i]=1);
 M[1,n]=1;
 P[1,1]=1;
}*/


