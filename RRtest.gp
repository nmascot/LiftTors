read("install.gp");
read("galrep.gp");
l=3;
read("H2/vap.gp");
p=29;
a=4;
e=1;
/*Lp=eval(Str("L_",l,"_",p));
\\C=x^2+7;
C=char2(A[4]);
d=poldegree(C);
a=mordroot(C,l);
e=8;*/
X=read("H2/RR.txt");f=X[1];L=X[2];Bad=X[3];L1=X[4];L2=X[5];

RR_rescale(L,p)=
{
	my(n,A,M);
	n = #L;
	M = L;
	for(i=1,#M,M[i]*=p^-valuation(M[i],p));
	M;
}

L = RR_rescale(L,p);
L1 = RR_rescale(L1,p);
L2 = RR_rescale(L2,p);
Li = [L1,L2];
Bad2 = lcm(apply(S->lcm(apply(f->denominator(content(f)),S)),[L,L1,L2]));
Bad = lcm(Bad,Bad2);

J=RRInit(f,7,16,L,Bad,p,a,e);
U=RREvalInit(J,Li);
/*J1 = PicRed(J,1); \\ Reduction mod p

print("\n--> Getting basis of T mod ",p)
WB = TorsBasis(J1,f,p,a,l,Lp,C);
print("\n--> Lifting ",p,"-adically");
my(J=J,l=l); WB = parapply(W->PicLiftTors(J,W,1,l),WB);
print("\n--> All of T");
T = TorsSpace(J,WB,l);
print("\n--> Evaluation");
my(J=J,U=U); Z = parapply(W->PlaneEval0(J,W,U),T[1..l^d-1]); \\ Evaluation of the torsion points in parallel
A = AllPols0(Z,JgetT(J),p,e,Jgetpe(J)); \\ p-adic approximation of a set of polynomials which all define a subfield of the field cut out by the representation, with equality iff. no repeated roots
print(#A," candidate polynomials");
AI = select(x->x[3]!=[],A); \\ Drop the approximations that could not be identified as rationals
print(#AI," identified polynomials");
AS = select(x->#Set(x[1])==#(x[1]),AI); \\ Drop the polynomials having multiple roots
print(#AS," faithful polynomials");
AS = vecsort(AS,x->sizebyte(x[3])); \\ Take the simplest of all polynomials
print(AS[1][3]);*/

