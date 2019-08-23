read("install.gp");
read("GalRep.gp");
l=3;
read("H2/vap.gp");
p=11;
\\ Get local L factor and charpoly of Frob
Lp=eval(Str("L_",l,"_",p));
C=char2(select(a->a[2]==p,A)[1]);
d=poldegree(C);
a=mordroot(C,l);
if(Mod(a,2),a=2*a); \\ Ensure Fq contains mu3
E=ellinit([0, -1, 0, -4, 4]); \\ Elliptic factor of the Jacobian
e=1024;
d0=16;
C=subst(C,x,-x); \\ Twist
X=read("H2/RR.txt");f=X[1];L=X[2];LL=X[3];Bad=X[4];L1=X[5];L2=X[6];g=7;

RR_rescale(L,p)=
{
	my(n,A,M);
	n = #L;
	M = L;
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

J=PicInit(f,g,d0,L,LL,Bad,p,a,e);
U=PicEvalInit(J,Li);
J1 = PicRed(J,1); \\ Reduction mod p

[B,matFrob] = TorsBasis(J1,l,Lp,C);
[WB,cWB] = TorsSpaceFrobGen(J1,l,B,matFrob); \\ Generating set of T under Frob and coordinates of these generators on B
print("\n--> Lifting ",#WB," points ",p,"-adically");
{if(#WB > Jgetg(J),
  my(J=J,l=l); WB = parapply(W->PicLiftTors(J,W,1,l),WB); \\ More efficient in parallel
,
  WB = apply(W->PicLiftTors(J,W,1,l),WB); \\ Less efficient in parallel (TODO tune)
);}
print("\n--> All of T");
TI = TorsSpaceFrob(J,WB,cWB,l,matFrob);
print("\n--> Evaluation of ",#TI[2]," points");
Z = TorsSpaceFrobEval(J,TI,U,l,d,matFrob);
print("\n--> Expansion and identification");
AF = TorsSpaceGetPols(J,Z); \\ List of polynomials defining the representation
F = AF[1][3]
