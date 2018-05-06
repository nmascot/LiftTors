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
