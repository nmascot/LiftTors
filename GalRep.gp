read("TorsSpace.gp");
read("TorsGen.gp");
read("Hyper2RR.gp");
read("Smooth2RR.gp");
read("Super2RR.gp");

mordroot1(f,p)=
\\ Computes the order of x in Fp[x]/(f). Assumes f irreducible mod p.
{
 my(x=variable(f),N,fa,l,v);
 N=p^poldegree(f)-1;
 fa=factor(N);
 for(i=1,#fa~,
  [l,v]=fa[i,];
  while(v,
   if(Mod(x^(N/l)-1,p)%Mod(f,p),break);
   N/=l;
   v-=1
  )
 );
 N;
}

mordroot(f,p)=
\\ Computes an upper bound for the order of x in Fp[x]/f. Gives exact value if f sqfree mod p.
{
	my(fa,N,e);
	fa=factormod(f,p);
	N=lcm(apply(g->mordroot1(g,p),fa[,1]));
	e=vecmax(fa[,2]);
	if(e>1,
	 e=p^ceil(log(e)/log(p));
	 warning("mordroot returning upper bound ",N*e,", but order may be as low as ",N);
	 N*=e
	);
	N;
}

PicIsTorsion(J,W,N)=PicIsZero(J,PicFrobPoly(J,W,N)); \\ For debugging purposes. N can be an integer or a polynomial.

TorsSpaceGetPols(J,Z)=
\\ Given a vector Z of evaluations of the points of a submodule T of J[l],
\\ Computes polynomials defining the Galois representation afforded by T
\\ and returns then ordered by height,
\\ each polynomial being given by a triplet:
\\ [list of p-adic roots, p-adic approximation, and rational identification (non-rigorous)]
{
  my(A,AI,AF);
  A = AllPols(Z,JgetT(J),Jgetp(J),Jgete(J),Jgetpe(J)); \\ p-adic approximation of a set of polynomials which all define a subfield of the field cut out by the representation, with equality iff. no repeated roots
	print(#A," candidate polynomials");
  AI = select(x->x[3]!=[],A); \\ Drop the approximations that could not be identified as rationals
  print(#AI," identified polynomials");
  AF = select(x->#Set(x[1])==#(x[1]),AI); \\ Drop the polynomials having multiple roots
  print(#AF," faithful polynomials");
	if(#AF==0 && #A==#AI,error("None of the evaluation maps gives a squarefree polynomial. Try again with different points."));
  vecsort(AF,x->sizebyte(x[3]));
}

RR_rescale(L,p)=
{
  my(n,A,M);
  n = #L;
  M = L;
  for(i=1,#M,M[i]*=p^-valuation(M[i],p));
  M;
}

GalRep(C,l,p,e,Lp,chi,force_a)=
/* Main function.
	 Given C=[f,g,d0,L,LL,L1,L2,Bad]
	 where f(x,y)=0 defines a curve C of genus g
	 and L=L(D0), LL=L(2D0), L1=L(2D0-E1), L(2D0-E2)
	 Riemann-Roch spaces where D0, E1 and E2 are effective
	 of degrees d0, d0-g, and d0-g,
	 where d0 > 2*g,
	 Computes the Galois representation afforded by
	 the piece of l-torsion of the Jacobian
	 on which Frob_p has charpoly chi
	 (chi=0 means take all the l-torsion)
	 by working at p-adic accuracy O(p^e).
	 Lp must be the(monic) local L-factor of the Jacobian at p,
	 and if chi is nonzero,
	 we must have chi || (Lp mod l).*/
{
	my([f,g,d0,L,LL,L1,L2,Bad]=C,d,J,J1,U,B,matFrob,WB,cWB,TI,Z,AF,F,ZF,M,i,e1=1);
	/* TODO rescale to remove denoms */
	L = RR_rescale(L,p);
  LL = RR_rescale(LL,p);
  L1 = RR_rescale(L1,p);
  L2 = RR_rescale(L2,p);
  Bad *= lcm(apply(S->lcm(apply(f->denominator(content(f)),S)),[L,L1,L2]));
  if(chi,
		print("T = part of J[",l,"] where Frob_",p," acts by ",chi);
		d = poldegree(chi); \\ Dimension of representation
		a = if(force_a,force_a,mordroot(chi,l)) \\ q = p^a
	,
		print("T = all of J[",l,"]");
		d=2*g;
		a = if(force_a,force_a,mordroot(Lp,l))
	);
	print("Working with q=",p,"^",a);
	J=PicInit(f,g,d0,[L,LL,L1,L2],Bad,p,a,e);
	J1 = PicRed(J,1); \\ Reduction mod p
	[B,matFrob] = TorsBasis(J1,l,Lp,chi); \\ Basis of the mod p^1 space and matrix of Frob_p
	print("The matrix of Frob is");
	printp(centerlift(matfrobenius(Mod(matFrob,l))));
	i=1;M=Mod(matFrob,l);
	while(M!=1,M*=matFrob;i++);
	print("It has order ",i);
	if(i<a,warning("Therefore working in degree a=",a," is not optimal. Consider restarting the computation while forcing a=",i,"."));
	[WB,cWB] = TorsSpaceFrobGen(J1,l,B,matFrob); \\ Generating set of T under Frob and coordinates of these generators on B
	while(1,
		print("\n--> Lifting ",#WB," points ",p,"-adically");
		if(#WB > Jgetg(J),
  		my(J=J,e1=e1,l=l); WB = parapply(W->PicLiftTors(J,W,e1,l),WB); \\ More efficient in parallel
		,
  		WB = apply(W->PicLiftTors(J,W,e1,l),WB); \\ Less efficient in parallel (TODO tune)
		);
		print("\n--> All of T");
		TI = TorsSpaceFrob(J,WB,cWB,l,matFrob);
		print("\n--> Evaluation of ",#TI[2]," points");
		Z = TorsSpaceFrobEval(J,TI,l,d,matFrob);
		print("\n--> Expansion and identification");
		AF = TorsSpaceGetPols(J,Z); \\ List of polynomials defining the representation
		if(AF!=[],break);
		e1=e;
		e*=2;
		warning("Could not identify any squarefree polynomial. Increasing p-adic accuracy: ",O(p^e),".");
		J = Jlift(J,e);
	);
	F = AF[1][3];
	if(#variables(F)>1,error("F has more than one variable"));
	ZF = apply(z->Mod(apply(c->c+O(p^e),z),JgetT(J)),AF[1][1]);
	[F,ZF];
}

HyperGalRep(f,l,p,e,P1,P2,chi,force_a)=
/* Computes the Galois representation afforded by
   the piece of l-torsion of the Jacobian
   of the hyperelliptic curve C:y²=f(x)
	 (in case f=[P,Q], the curve C:y²+Q(x)*y=P(x))
	 on which Frob_p has charpoly chi
   (chi=0 means take all the l-torsion)
   by working at p-adic accuracy O(p^e).
	 P1 and P2 must be two points of C(Q)
	 which are not conjugate under the hyperelliptic involution.
   If chi is nonzero,
   we must have chi || (Lp mod l)
	 where Lp is the local L factor at p.*/
{
	my(Lp,C);
	Lp = hyperellcharpoly(Mod(f,p)); \\ Local L factor of the curve at p, needed to know the number of points on the Jacobian mod p
	C = Hyper2RR(f,P1,P2);
	C=concat(C,['y]);
	GalRep(C,l,p,e,Lp,chi,force_a);
}

SmoothGalRep(f,l,p,e,P1,P2,chi,force_a)=
/* Computes the Galois representation afforded by
   the piece of l-torsion of the Jacobian
   of the plane curve f(x,y)=0
   on which Frob_p has charpoly chi
   (chi=0 means take all the l-torsion)
   by working at p-adic accuracy O(p^e).
	 Assumes f(x,y)=0 is smooth in P^2 and has good reduction at p.
   P1 and P2 must be two sets of points of C(Q) (of the form [x,y,z])
	 of size d-g TODO, where d=deg(f) and g=(d-1)(d-2)/2
   If chi is nonzero,
   we must have chi || (Lp mod l)
   where Lp is the local L factor at p.*/
{
  my(Lp,C);
  C = Smooth2RR(f,P1,P2);
	Lp = PlaneZeta(C[1],p); \\ Local L factor at p
  C=concat(C,[1]);
  GalRep(C,l,p,e,Lp,chi,force_a);
}

SuperGalRep(f,m,l,p,e,P,chi,force_a)=
/* Computes the Galois representation afforded by
   the piece of l-torsion of the Jacobian
   of the superelliptic curve y^m=f
   on which Frob_p has charpoly chi
   (chi=0 means take all the l-torsion)
   by working at p-adic accuracy O(p^e).
   Requires f squarefree mod p and m coprime with deg(f).
	 If chi is nonzero,
   we must have chi || (Lp mod l)
   where Lp is the local L factor at p. */
{
	my(Lp,C);
	if(!issquarefree(Mod(f,p)),error(f," is not squarefree mod ",p));
	C = Super2RR(f,m,P);
	Lp = SuperZeta(f,m,p);
	C = concat(C,['y]);
	GalRep(C,l,p,e,Lp,chi,force_a);
}

HyperBestp(f,l,pmax)=
{
	my(D,P,A,a,i);
	if(type(f)=="t_VEC",
		D = poldisc(4*f[1]+f[2]^2)
	,
		D = poldisc(f)
	);
	D *= l;
	P = primes([3,pmax]);
	P = select(p->Mod(D,p),P);
	export(mordroot,mordroot1);
	A = parapply(p->mordroot(hyperellcharpoly(Mod(f,p)),l),P);
	a = vecmin(A,&i);
	[P[i],a];
}

SmoothBestp(f0,D,l,pmax)=
{
	\\ TODO compute D
	my(x,y,d,f,P,A,a,i);
	[x,y] = variables(f0);
  d = TotalDeg(f0,x,y);
  f = SmoothGeneric(f0,d)[1];
	D *= l;
	P = primes([3,pmax]);
	P = select(p->Mod(D,p),P);
	export(mordroot,mordroot1);
	A = parapply(p->mordroot(PlaneZeta(f,p),l),P);
	a = vecmin(A,&i);
	[P[i],a];
}
