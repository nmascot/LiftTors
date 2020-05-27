A2<x,y>:=AffinePlane(QQ);

RRbasis:=function(D)
	L,l:=RiemannRochSpace(D);
	L:=[l(Basis(L)[i]) : i in [1..Dimension(L)]];
	return L;
end function;

RRGP:=procedure(D0,E1,E2,path)
	Write(path,"{f=");
	Write(path,Equation(Curve(D0)));
	Write(path,";g=");
	Write(path,Genus(Curve(D0)));
	Write(path,";d0=");
	Write(path,Degree(D0));
	Write(path,";L=");
	Write(path,RRbasis(D0));
	Write(path,";LL=");
	Write(path,RRbasis(2*D0));
	Write(path,";L1=");
	Write(path,RRbasis(2*D0-E1));
	Write(path,";L2=");
	Write(path,RRbasis(2*D0-E2));
	Write(path,";}");
end procedure;

Lp:=function(f,p)
	return LPolynomial(Curve(AffinePlane(GF(p)),f));
end function;

LpGP:=procedure(f,p,path)
	P:=Lp(f,p);
  _<x>:=Parent(P);
  P:=Parent(P)!Reverse(Coefficients(P));
	Write(path,"{L");
	Write(path,p);
	Write(path,"=");
	Write(path,P);
	Write(path,";}");
end procedure;
