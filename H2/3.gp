Fa=x*y*(x^2-1)*(y^2-1)*(x^2-y^2+a*x*y);
F=subst(Fa,a,2);
F=subst(F,y,l*x);
F=F/x^4;
A=ellfromeqn(y^2-F);

ellDivpol(A,m)=
{
	my(E,D,a1,a2,a3,a4,a6);
	[a1,a2,a3,a4,a6]=A;
	E=ellinit(A);
	D=elldivpol(E,m);
	polresultant(D,y^2+a1*x*y+a3*y-(x^3+a2*x^2+a4*x+a6),x);
}

F=ellDivpol(A,3);
F=subst(F,l,x);
