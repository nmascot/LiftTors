install("ZpXQ_div","GGGGGL");

\\ Choose p^e and N
p=29;
e=4;
pe=p^e;
N=12;
PhiN=polcyclo(N);
\\ Choose E
a4=a6=1;
fE='x^3+a4*'x+a6;
E=ellinit([a4,a6]);
\\ Choose q=p^a splitting E[N]
D=elldivpol(E,N);
a=8;
t1=ffgen([p,a],'t);
E1=ellinit(E,t1); \\ E / Fq
T=t1.mod;
\\ Frob_p on Qq
FrobT=padicappr(T,Mod(t^p,T)+O(p^e))[1];
Frob(x)=Mod(subst(liftpol(x),'t,FrobT),T);

\\ The x of E[N]
xD=polrootsmod(D,t1);

GetPt(x,f)=[x,sqrt(subst(f,'x,x))]; \\ x -> [x,y] on E

\\ Lift P in E[N](Fq) to E[N](Qq). f Wei, D divpoly.
LiftTorsPt(P,f,D,T,p,e)=
{
 my(a,b);
 if(P==0,return([0]));
 a=P[1].pol;
 if(subst(f,'x,P[1])==0,
  a=padicappr(f,Mod(a,T)+O(p^e))[1];
  return([liftall(a),0]);
 );
 a=padicappr(D,Mod(a,T)+O(p^e))[1];
 b=sqrt(subst(f,'x,P[1]));
 if(b!=P[2],b=-b;if(b!=P[2],error("b")));
 b=b.pol;
 b=padicappr('y^2-subst(f,'x,a),Mod(b,T)+O(p^e))[1];
 liftall([a,b]);
}

\\ E / Fq=Fp(t) -> Weil paring e_N(P,Q) in Qq \\ TODO pas utile?
/*PtRed(P,t)=subst(Mod(liftpol(P),t.p),variable(t.mod),t);
eWeil(E,P,Q,N,t)=padicappr('x^N-1,Mod(ellweilpairing(ellinit(E,t),PtRed(P,t),PtRed(Q,t),N).pol,t1.mod)+O((t.p)^padicprec(polcoef(liftpol(P[1]),0),t.p)))[1];*/
printf("Gens");
{ \\ Get generators P0,Q0 of E[N](Qq)
while(1,
 P0=GetPt(xD[1+random(#xD)],fE); \\ Random elt P0 of E[N](Qq)
 Q0=GetPt(xD[1+random(#xD)],fE);
 zetaN=ellweilpairing(E1,P0,Q0,N);
 if(subst(PhiN,'x,zetaN)==0,break)
);
}
printf("Tors");
EN=matrix(N,N); \\ [[ m P0 + n Q0 ]]: this is a naive level structure alpha: (Z/NZ)² ~ E[N]
{
EN[1,N]=P0;
EN[N,1]=Q0;
for(x=2,N-1,
	EN[x,N]=elladd(E1,EN[x-1,N],P0);
	EN[N,x]=elladd(E1,EN[N,x-1],Q0);
);
for(x=1,N-1,
	for(y=1,N-1,
		EN[x,y]=elladd(E1,EN[x,N],EN[N,y])
	)
);
}
\\ Lift to Qq
/*export(T);
export(p);
export(e);
export(D);
export(fE);
export(LiftTorsPt);*/
EN=apply(P->liftall(LiftTorsPt(P,fE,D,T,p,e)),EN);
GetCoef(A,a)=my(i,j);[i,j]=liftall(Mod(a,N));if(i==0,i=N);if(j==0,j=N);A[i,j];

l2(P,Q)=my([xP,yP]=GetCoef(EN,P),[xQ,yQ]=GetCoef(EN,Q));ZpXQ_div(yQ-yP,xQ-xP,T,pe,p,e); \\ Slope of line (PQ)
\\ TODO methode Kamal addchain
l1(P,Q)=sum(n=0,N-1,l2(P,Q+n*P)); \\ sum_n l2(P,Q+nP) ~ sum_n l(P) +l(Q+nP)-l(Q+nP) ~ l*l(P)

printf("lambdas");
\\ Matrix of l(P) for P in E[N]
Ml1=matrix(N,N);
for(x=1,N-1,Ml1[x,N]=l1([x,0],[0,1])); \\ P=alpha(x,0) -> Q=alpha(0,1)
for(x=1,N,for(y=1,N-1,Ml1[x,y]=l1([x,y],[1,0]))); \\ P=alpha(x,y), y!=0 -> Q=alpha(1,0)

printf("SL2");
\\ Find all of SL(2,l), TODO that's 12x too many pts, restrict.
SL2=List();
{
for(a=0,N-1,
	for(b=0,N-1,
		for(c=0,N-1,
			for(d=0,N-1,
				if(a*d-b*c==Mod(1,N),listput(SL2,[a,b;c,d]))
			)
		)
	)
);
SL2=Vec(SL2);
}

printf("Cusps");
\\ Find all (u,v) s.t. gcd(u,v,N)=1
Cusps=List();
{
for(u=0,N-1,
	for(v=0,N-1,
		if(gcd([u,v,N])==1,listput(Cusps,[u,v]))
	)
);
Cusps=Vec(Cusps);
}

ActOnCuspi(g,i)=\\ g in SL2 -> the j s.t. g*cusp#i = cusp#j
{
my(c=Cusps[i]);
c=c*g;
c=liftint(Mod(c,N));
select(x->x==c,Cusps,1)[1];
}

{
\\ Matrix representing Eis1(Gamma1(N))
E11=matrix(#Cusps,#Cusps);
for(c=1,#Cusps, \\ this col = l(c)
	for(P=1,#Cusps, \\ this row = values at pt alpha*g, g=SL2[P]
		[a,b]=Cusps[P];
		[u,v,d]=gcdext(a,b);
		x=liftint(1/Mod(d,N));
		E11[P,c]=sum(n=1,N,GetCoef(Ml1,Cusps[c]*[a,b;-v*x+a*n,u*x+b*n]))
	)
);
}

{
\\ Matrix representing Eis1(Gamma(N))
E1=matrix(#SL2,#Cusps);
for(c=1,#Cusps, \\ this col = l(c)
   for(P=1,#SL2, \\ this row = values at pt alpha*g, g=SL2[P]
      E1[P,c]=GetCoef(Ml1,Cusps[c]*SL2[P]);
   )
);
}

\\ Values at OO
z = Mod(Mod(zetaN.pol,p),T);
{
E1_OO = matrix(#SL2,#Cusps);
for(P=1,#SL2,
for(v=1,#Cusps,
	[c,d] = Cusps[v]*SL2[P];
	E1_OO[P,v] = if(c%N==0,
		1/2*(z^d+1)/(z^d-1),
		Mod(1/2-frac(c/N),p))
)
);
}

{
\\ Matrix representing Eis2(Gamma(N))
E2=matrix(#SL2,#Cusps);
for(c=1,#Cusps, \\ this col = l(c)
   for(P=1,#SL2, \\ this row = values at pt alpha*g, g=SL2[P]
      E2[P,c]=GetCoef(EN,Cusps[c]*SL2[P])[1];
   )
);
}

/*M=Mod(Mod(Ml1,p),T);
{
\\ Matrix representing S2(Gamma(N))
S2=matrix(#SL2,N^4);
c=0;
for(x1=0,N-1,
	for(y1=0,N-1,
		print([x1,y1]);
		for(x2=0,N-1,
			for(y2=0,N-1,
				P1=[x1,y1];P2=[x2,y2];
				P3=-P1-P2;
				if(Mod(P1,N)==0,next);
				if(Mod(P2,N)==0,next);
				if(Mod(P3,N)==0,next);
				c++;
   			for(P=1,#SL2, \\ this row = values at pt alpha*g, g=SL2[P]
      			S2[P,c]=GetCoef(M,P1*SL2[P])*GetCoef(M,P2*SL2[P])
						+GetCoef(M,P2*SL2[P])*GetCoef(M,P3*SL2[P])
						+GetCoef(M,P3*SL2[P])*GetCoef(M,P1*SL2[P]);
				)
			)
		)
	)
);
}*/

E1=Mod(Mod(E1,p),T);
E11=Mod(Mod(E11,p),T);
E2=Mod(Mod(E2,p),T);
/*M2=matrix(#SL2,(#Cusps)^2);
{
for(c1=1,#Cusps,
	for(c2=1,#Cusps,
		c=(c1-1)*(#Cusps)+c2;
		for(P=1,#SL2,M2[P,c]=E1[P,c1]*E1[P,c2])
	)
);
}*/

E=matconcat([E1,E1_OO]~);
matrank(E)
