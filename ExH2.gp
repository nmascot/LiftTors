read("install.gp");
read("galrep.gp");
l=3;
read("H2/vap.gp");
p=11;
Lp=eval(Str("L_",l,"_",p));
\\C=x^2+7;
iA=4;
if(A[iA][2]!=p,error("Wrong A"));
C=char2(A[iA]);
d=poldegree(C);
a=mordroot(C,l);
if(Mod(a,2),a=2*a);
e=4;
X=read("H2/RR.txt");f=X[1];L=X[2];LL=X[3];Bad=X[4];L1=X[5];L2=X[6];g=7;d0=16;

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

L = RR_rescale(L,p);
LL = RR_rescale(LL,p);
L1 = RR_rescale(L1,p);
L2 = RR_rescale(L2,p);
Li = [L1,L2];
Bad2 = lcm(apply(S->lcm(apply(f->denominator(content(f)),S)),[L,L1,L2]));
Bad *= Bad2;

J=RRInit2(f,g,d0,L,LL,Bad,p,a,e);
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
A = AllPols(Z,JgetT(J),p,e,Jgetpe(J)); \\ p-adic approximation of a set of polynomials which all define a subfield of the field cut out by the representation, with equality iff. no repeated roots
print(#A," candidate polynomials");
AI = select(x->x[3]!=[],A); \\ Drop the approximations that could not be identified as rationals
print(#AI," identified polynomials");
AS = select(x->#Set(x[1])==#(x[1]),AI); \\ Drop the polynomials having multiple roots
print(#AS," faithful polynomials");
AS = vecsort(AS,x->sizebyte(x[3])); \\ Take the simplest of all polynomials
print(AS[1][3]);

