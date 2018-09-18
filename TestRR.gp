read("install.gp");
read("galrep.gp");
read("Smooth2RR.gp");
l=2;
f=x^3*y+y^3+x;
P1=[0,0,1];
P2=[1,0,0];
[f,g,d0,L,L1,L2] = Smooth2RR(f,[P1],[P2]);
p=5;
C=0;
d=6;
a=6;
e=32;
RR_rescale(L,p)=
{
	my(n,A,M);
	n = #L;
	M = L;
	/*until(Mod(matdet(A),p),
		A = matrix(n,n,i,j,random(p))
	);
	print("Change by:");
	print(A);
	M = L*A;*/
	for(i=1,#M,M[i]*=p^-valuation(M[i],p));
	M;
}

Li = [L1,L2];
Lp = PlaneZeta(f,p);

J=RRInit(f,g,d0,L,1,p,a,e);
U=RREvalInit(J,Li);
J1 = PicRed(J,1); \\ Reduction mod p

[B,matFrob] = TorsBasis(J1,l,Lp,C);
[F,Q] = matfrobenius(Mod(matFrob,l),2); \\ F = rational canonical form of matFrob, Q = transition matrix
Q = lift(Q^(-1));
NF = apply(poldegree,matfrobenius(F,1)); \\ List of degrees of elementary divisors
WB = List(); cWB = List();
n = 1;
{for(i=1,#NF,
  c=Q[,n];
  listput(cWB,c);
  c=centerlift(Mod(c,l));
  listput(WB,PicLC(J1,c,B));
  n+=NF[i]
);}
WB = Vec(WB); \\ Generating set of T under Frob
cWB = Vec(cWB); \\ Coordinates of these generators on B
print("\n--> Lifting ",#WB," points ",p,"-adically");
{if(#WB > Jgetg(J),
  my(J=J,l=l); WB = parapply(W->PicLiftTors(J,W,1,l),WB); \\ More efficient in parallel
,
  WB = apply(W->PicLiftTors(J,W,1,l),WB); \\ Less efficient in parallel (TODO tune)
);}
print("\n--> All of T");
[T,ImodF] = TorsSpaceFrob(J,WB,cWB,l,matFrob);
print("\n--> Evaluation of ",#ImodF," points");
Z = vector(l^d-1,i,[]);
my(J=J,T=T,U=U);ZmodF = parapply(i->RREval(J,T[i],U),ImodF);
{for(n=1,#ImodF,
  i = ImodF[n];
  z = ZmodF[n];
  while(1,
    Z[i] = z;
    i = ActOni(matFrob,i,l);
    if(Z[i] != [],break);
    z = apply(x->Frob(x,JgetFrobMat(J),JgetT(J),Jgetpe(J)),z);
  )
);}
print("\n--> Expansion and identification");
A = AllPols0(Z,JgetT(J),p,e,Jgetpe(J)); \\ p-adic approximation of a set of polynomials which all define a subfield of the field cut out by the representation, with equality iff. no repeated roots
print(#A," candidate polynomials");
AI = select(x->x[3]!=[],A); \\ Drop the approximations that could not be identified as rationals
print(#AI," identified polynomials");
AS = select(x->#Set(x[1])==#(x[1]),AI); \\ Drop the polynomials having multiple roots
print(#AS," faithful polynomials");
AS = vecsort(AS,x->sizebyte(x[3])); \\ Take the simplest of all polynomials
print(AS[1][3]);

