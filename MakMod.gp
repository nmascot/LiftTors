install("ZpXQ_div","GGGGGL");
read("Dimensions.gp");

\\ Choose p^e and N
N=18;
\\ Choose E
a4=a6=1;
fE='x^3+a4*'x+a6;
E=ellinit([a4,a6]);
\\ Choose q=p^a splitting E[N]
D=elldivpol(E,N);
D/=content(D);

BestpSplitE(fE,D,N,NE,pmax)=
{
	my(m=0,D2,p,a,p0=0,a0=0,fa);
	while(1,
		D2=polresultant(D,('y-m*'x)^2-fE,'x);
		D2=D2/gcd(D2,D2');
		if(poldegree(D2)==N^2-1,break);
		m=if(m>0,-m,-m+1);
	);
	D2/=content(D2);
	forprime(p=2,pmax,
		if(Mod(N*NE,p)==0,next);
		if(p0 && znorder(Mod(p,N))>=a0,next);
		fa=factormod(D2,p,1);
		print([p,vecmax(fa[,2]),[p0,a0]]);
		if(vecmax(fa[,2])>1,next);
		a=vecmax(fa[,1]);
		if(p0==0,[p0,a0]=[p,a]);
		if(a<a0,[p0,a0]=[p,a])
	);
	[p0,a0];
}

[p,a]=BestpSplitE(fE,D,N,E.disc,1000);
e=4;
pe=p^e;
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

ENtors(E,fE,N,D,xD)=
{
	my(PhiN,P0,Q0,zetaN,EN);
	PhiN=polcyclo(N);
	\\ Get generators P0,Q0 of E[N](Fq)
	while(1,
		P0=GetPt(xD[1+random(#xD)],fE); \\ Random elt P0 of E[N](Fq)
 		Q0=GetPt(xD[1+random(#xD)],fE);
 		zetaN=ellweilpairing(E,P0,Q0,N);
 		if(subst(PhiN,'x,zetaN)==0,break)
	);
	EN=matrix(N,N); \\ [[ m P0 + n Q0 ]]: this is a naive level structure alpha: (Z/NZ)² ~ E[N]
	EN[1,N]=P0;
	EN[N,1]=Q0;
	for(x=2,N-1,
  	EN[x,N]=elladd(E,EN[x-1,N],P0);
  	EN[N,x]=elladd(E,EN[N,x-1],Q0);
	);
	for(x=1,N-1,
  	for(y=1,N-1,
    	EN[x,y]=elladd(E,EN[x,N],EN[N,y])
  	)
	);
	[EN,zetaN];
}

[EN,zetaN] = ENtors(E1,fE,N,D,xD);

\\ Lift to Qq
/*export(T);
export(p);
export(e);
export(D);
export(fE);
export(LiftTorsPt);*/
EN=apply(P->liftall(LiftTorsPt(P,fE,D,T,p,e)),EN);

GetCoef(A,a)=my(i,j);[i,j]=liftall(Mod(a,N));if(i==0,i=N);if(j==0,j=N);A[i,j];
ZNnorm(x,N)=my(y=x%N);if(y==0,N,y);
ZNneg(x,N)=my(y=lift(Mod(-x,N)));if(y==0,N,y);

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
\\ Find all (u,v) s.t. gcd(u,v,N)=1 / +-1
Cusps=List();
done=matrix(N,N);
{
for(u=0,N-1,
	for(v=0,N-1,
		if(GetCoef(done,[u,v]),next);
		if(gcd([u,v,N])==1,
			listput(Cusps,[u,v]);
			done[ZNnorm(u,N),ZNnorm(v,N)]=1;
			done[ZNnorm(-u,N),ZNnorm(-v,N)]=1;
		)
	)
);
Cusps=Vec(Cusps);
}

SL2lift(M)=
{
	my(U,V,D,a,b);
	[U,V,D]=matsnf(M,1);
	a=D[1,1];
	b=D[2,2];
	U^-1*[1,-1;1-b,b]*[1,0;1-a,1]*[1,b;0,1]*V^-1;
}

{
ToCusp(s,N)= \\ Finds [p',*;q';*] in SL2(Z) with p~p', q~q' mod N
	my(p,q,u,v,g,h);
	[p,q]=s;
	[u,v,g]=gcdext(p,q);
	h=lift(1/Mod(g,N));
	SL2lift([p,-v*h;q,u*h]);
}

{
Cusp2X1(s,N)= \\ Finds [*,*;c';d'] in SL2(Z) with c~c', d~d' mod N
  my(c,d,u,v,g,h);
  [c,d]=s;
  [u,v,g]=gcdext(c,d);
  h=lift(1/Mod(g,N));
  SL2lift([c,-v*h;d,u*h]);
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
\\ Values at cusps
z = Mod(Mod(zetaN.pol,p),T);
E1atCusp(s,N,z)=
{
	my([c,d]=s);
	c=c%N;
	d=d%N;
	if(c==0,
		1/2*(z^d+1)/(z^d-1),
		1/2-c/N
	);
}

\\ TODO
E1 = Mod(Mod(E1,p),T);

{
\\ Generators of M2(Gamma(N))
d=DimMk(2,-1,N); \\ dim
d2=min(2*d,#SL2); \\ # pts at which we observe lin indep
while(1,
	Mtest=matrix(d2,d);
	M2gens=vector(d);
	for(i=1,d,
		u=random(#Cusps)+1;
		v=random(#Cusps)+1;
		M2gens[i]=[u,v];
		for(P=1,d2, \\ this row = values at pt alpha*g, g=SL2[P]
    	Mtest[P,i]=E1[P,u]*E1[P,v];
  	);
	);
	if(matrank(Mtest)==d,break);
	print("M2 fail");
);
}

/*TODO test: On M2(Gamma(N)), sum_c cusp f(c) = 0? OK!
{
for(i=1,#Cusps,
	for(j=1,#Cusps,
		S=0;
		for(s=1,#Cusps,
    	c=Cusps[s];
    	Mc=ToCusp(c,N);
			S+=E1atCusp(Cusps[i]*Mc,N,z)*E1atCusp(Cusps[j]*Mc,N,z);
		);
		print([i,j,liftall(S)]);
	)
);
}
 TODO end test */

M21=matrix(#Cusps,d);
{
	for(P=1,#Cusps,
		MP=Cusp2X1(Cusps[P],N);
		for(i=1,d,
			[v,w]=M2gens[i];
			v=Cusps[v];
			w=Cusps[w];
			M21[P,i]=sum(x=1,N,GetCoef(Ml1,v*[1,x;0,1]*MP)*GetCoef(Ml1,w*[1,x;0,1]*MP))
		)
	);
}
M21=Mod(Mod(M21,p),T);
\\ Extract basis
M21big=M21;
B=matindexrank(M21big)[2];
d=#B;
M21=matrix(#Cusps,d);
M2basis=vector(d,j,M2gens[B[j]]);
{
for(i=1,#Cusps,
	for(j=1,d,
		M21[i,j]=M21big[i,B[j]]
	)
);
}

e01=DimData(1,N)[4];
M21_cusps=matrix(e01,d); \\ Values at cusps / width = residues
{
	n=1;
	done=vector(#Cusps);
	Cusps1=vector(e01);
	for(m=1,#Cusps,
		if(done[m],next);
		s=Cusps[m];
		Cusps1[n]=s;
		h=gcd(N,s[2]); \\ Width of cusp wrt Gamma1(N)
		for(x=1,N,
			sx=s*[1,0;x,1]; \\ ([1,x;0,1]*s~)~
			done[select(t->Mod(t,N)==sx||Mod(-t,N)==sx,Cusps,1)[1]]=1
		);
		Ms=ToCusp(s,N);
		for(i=1,d,
			[v,w]=M2basis[i];
      v=Cusps[v];
      w=Cusps[w];
			M21_cusps[n,i]=sum(x=1,N,E1atCusp(v*[1,x;0,1]*Ms,N,z)*E1atCusp(w*[1,x;0,1]*Ms,N,z))/h
		);
		n++
	);
}

M21=Mod(Mod(M21,p),T);

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

/*M2=matrix(#SL2,(#Cusps)^2);
{
for(c1=1,#Cusps,
	for(c2=1,#Cusps,
		c=(c1-1)*(#Cusps)+c2;
		for(P=1,#SL2,M2[P,c]=E1[P,c1]*E1[P,c2])
	)
);
}*/
