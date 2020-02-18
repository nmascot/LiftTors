read("TorsSpace.gp");
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


PicLC(J,C,W)=
\\ Computes sum_i C[i]*W[i] in J
\\ C vector of integer coefficients, W vector of points on J
{
  my(I1,C1,W1,S,n);
	\\ Remove zeros
  I1=select(c->c,C,1);
	C1=vecextract(C,I1);
	W1=vecextract(W,I1);
	n=#I1;
  if(n==0,return(JgetW0(J)));
  S = PicMul(J,W1[1],(-1)^(n-1)*C1[1],2);
  for(i=2,n,
		\\ Now S = (-1)^(n-i+1) sum_{j<i} C1[j]*W1[j] 	
		S = PicChord(J,S,PicMul(J,W1[i],(-1)^(n-i+1)*C1[i],2),0);
		\\ Now S = (-1)^(n-i) sum_{j<=i} C1[j]*W1[j] 	
  );
	S;
}

TorsOrd(J,W,l)=
\\ Given that W is an l-power torsion point of J,
\\ finds v s.t. the order of W is l^v,
\\ and returns [l^(v-1)W, v]
{
  my(v=0,lW);
  v = 0;
  lW = W;
  while(!PicIsZero(J,lW),
    v += 1;
    W = lW;
    lW = PicMul(J,W,l,0)
  );
  [W,v];
}

RandTorsPt(J,l,a,Lp,Chi,Phi,seed)=
{
  my(N,v,M,Psi,fa,W,o,T,lT);
	setrand(seed);
	N = polresultant(Lp,if(Phi,Phi,'x^a-1));
	v = valuation(N,l);
  M = N/l^v;
  while(1,
    W = PicRand(J);
		if(Phi,W = PicFrobPoly(J,W,('x^a-1)/Phi));
    W = PicMul(J,W,M,0);
    if(Chi,
			Psi = lift(Mod(Lp,l)/Chi); \\ cofactor
    	fa = [Chi,Psi];
    	fa = polhensellift(Lp,fa,l,v); \\ lift cofactor l-adically
    	Psi = fa[2];
    	Psi = centerlift(Mod(Psi,l^v)); \\ center mod l^v 
			W = PicFrobPoly(J,W,Psi)
		);
		o = 0;
		T = W;
		lT = T;
		while(!PicIsZero(J,lT),
    	o += 1;
    	T = lT;
    	lT = PicMul(J,T,l,0)
  	);
		if(o,return([W,o,T,if(Phi,poldegree(Phi),a)]));
  );
}

PicIsTorsion(J,W,N)=PicIsZero(J,PicFrobPoly(J,W,N)); \\ For debugging purposes. N can be an integer or a polynomial.

TorsBasis(J,l,chi,C)=
\\ Computes a basis B of the subspace T of J[l] on which Fron acts with charpoly C
\\ Assumes chi = charpoly(Frob|J), so C | chi
\\ If C==0, then we take T=J[l]
\\ Also computes the matrix M of Frob w.r.t B, and returns the vector [B,M]
{
  my(i,a,d,Phi,iPhi,UsedPhi,Batch,nBatch,iBatch,iFrob,W,o,T,BW,Bo,BT,R,KR,Rnew,KRnew,Wtest,Wnew,am,S,AddC,W0,z);
	Phi = 0; \\ List of Cyclo Pols used to accelerated exponentiation if a%l
	iPhi = 0; \\ Index of last used elt of Phi
	nBatch = 0; \\ Size of current batch
	iBatch = 0; \\ Position in current batch
	iFrob = 0; \\ Number of mes Frob has been applied to current point of current batch
	iFrobMax = 0; \\ Upper bound on deg of minpoly of Frob on current point on current batch
	a = poldegree(JgetT(J)); \\ work over Fq=Fp[t]/T, q=p^a
  d = if(C,poldegree(C),poldegree(chi)); \\ dim of rep
	if(Mod(a,l),
		Phi = vecsort(divisors(a),,4); \\ Divisors of a in reverse order
		Phi = apply(polcyclo,Phi); \\ Corresponding cyclotomic pols
  	if(C,Phi = select(phi->poldegree(gcd(Mod(phi,l),Mod(C,l))),Phi)) \\ Keep the ones compatible with charpoly C
	);
  BW = vector(d); \\ list of l-power tors pts
  Bo = vector(d); \\ list of exponents of orders
  BT = vector(d); \\ list of l-tors pts
	Wtest = vector(d,i,PicChord(J,PicRand(J),PicRand(J),1)); \\ list of pts to pair l-tors with
	R = matrix(d,d); \\ matrix of pairings
	AddC = AddChain(l,0);
	W0 = JgetW0(J);
  W0 = PicChord(J,W0,W0,1); \\ Non-trivial origin, needed for pairings
	z = Fq_zeta_l(JgetT(J),Jgetp(J),l); \\ primitive l-th root of 1, to linearize parings
  r = 0; \\ dim of l-tors obtained so far

  while(r<d,
    print("Getting new point");
		iFrob += 1;
		if(iFrob>=iFrobMax, \\ No use applying Frob anymore
			iFrob = 0;
			print(" from batch");
			if(iBatch>=nBatch, \\ Have we used all of the previous batch?
				nBatch = max(ceil((d-r)/a),ceil(default(nbthreads)/2)); \\ If yes, make new batch
				print("  Generating a new batch of ",nBatch," points in parallel");
				my(RandTorsPt=RandTorsPt,seed=vector(nBatch,i,random()));
				if(Phi,
					UsedPhi = vector(nBatch,i,Phi[(iPhi+i-1)%#Phi+1]);
					iPhi += nBatch;
					print("   Using Phi=",UsedPhi)
				,
					UsedPhi = vector(nBatch,i,0);
				);
				Batch = parvector(nBatch,i,RandTorsPt(J,l,a,chi,C,UsedPhi[i],seed[i]));
				print("  Batch of points generated.");
				iBatch=1;
    	,
				iBatch+=1; \\ Else just move to next pt in batch
			);
			[W,o,T,iFrobMax]=Batch[iBatch];
    	print(" It has order l^",o);
		,
			print(" from Frob"); \\ Get new pt by Frob (fast)
			W = PicFrob(J,W);
			T = PicFrob(J,T)
		);
		\\ (H) At this point, the left d*r block of R has full rank r
    r += 1;
    BW[r] = W;
    Bo[r] = o;
    BT[r] = T;
    while(1, \\ Loop until indep new pt is found, or until we cannot exploit current pt any further
    	print(" Computing pairings...");
			\\ New col of R
			Rnew = PicFreyRuckMulti(J,T,l,Wtest,W0,AddC);
			for(i=1,d,
				R[i,r] = Fq_mu_l_log(Rnew[i],z,JgetT(J),Jgetp(J),l)
			);
    	print(" Looking for relations...");
    	KR = centerlift(matker(Mod(R[,1..r],l)));
			if(default(debug)>=2,
				print("R=");
				printp(R);
				print("KR=");
				printp(KR)
			);
			if(#KR>1,error("Bug in TorsSpace, please report")); \\ Not supposed to happen by (H)
			if(#KR==0,
				print(" Good, no relation");
				\\breakpoint();
				next(2)
			);
			KR=KR[,1];
			\\ So we have a pseudo-relation
	    print(" Found pseudo-relation");
			\\ Either this is an actual relation, or our linear forms are not independent
			if(PicIsZero(J,PicLC(J,KR,BT[1..r]))==0,
        print(" Good, it does not actually hold.");
				\\ So our linear forms are not independent.
				\\ Find a new one, and replace one the appropriate old one with it
				print(" Changing linear tests so that we don't get a false positive again.");
				\\ Note: we could exit there, but it is better to fix the linear forms
				\\ so as to compute the matrix of Frobenius later on.
				until(Mod(Rnew*KR,l), \\ loop until new pairing breaks pseudo-relation
          Wnew = PicChord(J,PicRand(J),PicRand(J),1); \\ New rand test pt
          Rnew = parapply(w->PicFreyRuckMulti(J,w,l,[Wnew],W0,AddC)[1],BT[1..r]); \\ Pair it with tors pts
          Rnew = apply(x->Fq_mu_l_log(x,z,JgetT(J),Jgetp(J),l),Rnew); \\ Discrete logs
					if(default(debug),print("New test gives parings ",Rnew));
        );
				\\ So now we have r+1 forms of rank r.
				\\ Find one that can be removed.
        Rnew = matconcat([R[,1..r],Rnew]~); \\ Combine with other pairings
				KRnew = centerlift(matker(Mod(Rnew~,l)))[,1]; \\ Find relation between forms
				i=1;
				while(KRnew[i]==0,i++);
				if(default(debug)>=2,print("Dropping form number ",i));
				\\ Replace the i-th old form with the new one
				Wtest[i] = Wnew;
				for(j=1,r,R[i,j] = Rnew[d+1,j]);
				if(default(debug)>=2,
					print("So now R=");
					printp(R)
				);
				next(2)
			,
				\\ So we have an actual relation.
				\\ Try to use it to make a new point.
      	m = vecmin([Bo[i]|i<-[1..r],KR[i]]);
				if(default(debug)>=2,
					print(" Bo=",Bo);
					print(" m=",m)
				);
				if(m>1,
      		if(default(debug)>=2,print(" Dividing relation ",KR~," by ",l));
					iFrob = 0;
      		S = vector(r,i,if(KR[i],l^(Bo[i]-m)*KR[i],0));
      		W = PicLC(J,S,BW[1..r]);
      		[T,o] = TorsOrd(J,W,l);
					if(default(debug)>=2,
						print(" Making comb ",S);
      			print(" gives point of order l^",o)
					);
					if(PicEq(J,BT[r],T), \\ This can happen if we have generated a Frob-stable suspace
						print("The new point is equal to the old one, giving up this point");
						iFrob = iFrobMax; \\ Reset indices
          	r -= 1; \\ Erase data about this point
          	break
					);
					BW[r] = W;
      		Bo[r] = o;
      		BT[r] = T;
				,
        	print(" Giving up this point");
					\\ Nothing more to do with this point. Start over with a new one.
					iFrob = iFrobMax; \\ Reset indices
        	r -= 1; \\ Erase data about this point
        	break
      	)
			);
    );
  );
	print("Found basis, now computing the matrix of Frobenius");
	\\ Apply Frobenius to BT, and pair
	FR = parapply(T->PicFreyRuckMulti(J,PicFrob(J,T),l,Wtest,W0,AddC),BT);
	FR = apply(x->Fq_mu_l_log(x,z,JgetT(J),Jgetp(J),l),matconcat(FR));
	\\ Find relations between pairings
	K = matker(Mod(matconcat([R,FR]),l));
	K = centerlift(-K[1..r,]);
  [BT,K];
}

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
