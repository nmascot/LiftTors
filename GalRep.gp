read("TorsSpace.gp");
read("Hyper2RR.gp");
read("Smooth2RR.gp");

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
	 print("mordroot warning: returning upper bound ",N*e,", but order may be as low as ",N);
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

RandTorsPt(J,l,M,chiC,seed)=
{
  my(W,o,T,lT);
	setrand(seed);
  while(1,
    W = PicRand(J);
    W = PicMul(J,W,M,0);
    if(chiC,W = PicFrobPoly(J,W,chiC));
		o = 0;
		T = W;
		lT = T;
		while(!PicIsZero(J,lT),
    	o += 1;
    	T = lT;
    	lT = PicMul(J,T,l,0)
  	);
		if(o,return([W,o,T]));
  );
}

TorsBasis(J,l,chi,C)=
\\ Computes a basis B of the subspace T of J[l] on which Fron acts with charpoly C
\\ Assumes chi = charpoly(Frob|J), so C | chi
\\ If C==0, then we take T=J[l]
\\ Also computes the matrix M of Frob w.r.t B, and returns the vector [B,M]
{
  my(a,d,N,M,v,chiC,Batch,nBatch,iBatch,iFrob,W,o,T,BW,Bo,BT,R,KR,Rnew,KRnew,Wtest,Wnew,am,S,AddC,W0,z);
	iBatch = nBatch = 0;
	a = poldegree(JgetT(J));
	iFrob = a-1;
	N = polresultant(chi,'x^a-1);
  v = valuation(N,l);
  M = N/l^v;
  if(C,
		d = poldegree(C);
    chiC = lift(Mod(chi,l)/C); \\ TODO test
    fa = [C,chiC];
    fa = polhensellift(chi,fa,l,v);
    chiC = fa[2];
    chiC = centerlift(Mod(chiC,l^v))
	,
		d = poldegree(chi);
    chiC = 0
  );
  BW = vector(d);
  Bo = vector(d);
  BT = vector(d);
	R = matrix(d,d);
	AddC = AddChain(l,0);
	W0 = JgetW0(J);
  W0 = PicChord(J,W0,W0,1);
	Wtest = vector(d,i,PicChord(J,PicRand(J),PicRand(J),1));
	z = Fq_zeta_l(JgetT(J),Jgetp(J),l);
  r = 0;
  while(r<d,
    /*print("Status:",Bo[1..r]);*/
    print("Getting new point");
		iFrob += 1;
		if(iFrob==a,
			print(" from batch");
			if(iBatch==nBatch,
				nBatch = max(ceil((d-r)/a),ceil(default(nbthreads)/2));
				print("  Generating a new batch of ",nBatch," points in parallel");
				my(RandTorsPt=RandTorsPt,seed=vector(nBatch,i,random()));
				Batch = parvector(nBatch,i,RandTorsPt(J,l,M,chiC,seed[i]));
				print("  Batch of points generated.");
				iBatch=1;
    	,
				iBatch+=1;
			);
			[W,o,T]=Batch[iBatch];
    	print(" It has order l^",o);
		,
			print(" from Frob");
			W = PicFrob(J,W);
			T = PicFrob(J,T)
		);
    r += 1;
    BW[r] = W;
    Bo[r] = o;
    BT[r] = T;
    while(1,		
    	print(" Looking for relations...");
			R[,r] = apply(x->Fq_mu_l_log(x,z,JgetT(J),Jgetp(J),l),PicFreyRuckMulti(J,T,l,Wtest,W0,AddC));
    	KR = centerlift(matker(Mod(R[,1..r],l)));
			if(#KR==0,
				print(" Good, no relation");
				next(2)
			);
			while(#KR>1,
				print("  Adding a linear test");
				Wnew = PicChord(J,PicRand(J),PicRand(J),1);
				Rnew = parapply(w->PicFreyRuckMulti(J,w,l,[Wnew],W0,AddC)[1],BT[1..r]);
				Rnew = apply(x->Fq_mu_l_log(x,z,JgetT(J),Jgetp(J),l),Rnew);
				Rnew = matconcat([R,Rnew]~);
				KRnew = centerlift(matker(Mod(Rnew[,1..r],l)));
				if(#KRnew<#KR,
					R = Rnew;
					KR = KRnew;
					Wtest = concat(Wtest,[Wnew])
				)
			);
      KR = KR[,1];
	    print("Found pseudo-relation ",KR~);
			if(PicIsZero(J,PicLC(J,KR,BT[1..r]))==0,
				print(" Good, it does not actually hold");
				next(2)
			);
      m = vecmin([Bo[i]|i<-[1..r],KR[i]]);
			if(m>1,
      	print(" Dividing relation ",KR," by l");
				iFrob = 0;
      	S = vector(r,i,if(KR[i],l^(Bo[i]-m)*KR[i],0));
      	W = PicLC(J,S,BW[1..r]);
      	[T,o] = TorsOrd(J,W,l);
      	print(" gives point of order l^",o)
			,
        print(" Giving up this point");
				iFrob = a-1;
        r -= 1;
        break
      );
      BW[r] = W;
      Bo[r] = o;
      BT[r] = T;
    );
  );
	print("Found basis, now computing the matrix of Frobenius");
	\\ Now compute matrix of Frobenius
	\\ First of all, make sure we have enough linear tests
	while(#KR,
  	print("  Adding a linear test");
    Wnew = PicChord(J,PicRand(J),PicRand(J),1);
    Rnew = parapply(w->PicFreyRuckMulti(J,w,l,[Wnew],W0,AddC)[1],BT[1..r]);
    Rnew = apply(x->Fq_mu_l_log(x,z,JgetT(J),Jgetp(J),l),Rnew);
    Rnew = matconcat([R,Rnew]~);
    KRnew = centerlift(matker(Mod(Rnew[,1..r],l)));
    if(#KRnew<#KR,
      R = Rnew;
      KR = KRnew;
      Wtest = concat(Wtest,[Wnew])
    )
  );
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

GalRep(C,l,p,e,Lp,chi)=
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
	my([f,g,d0,L,LL,L1,L2,Bad]=C,d,J,J1,U,B,matFrob,WB,cWB,TI,Z,AF,F,ZF);
	/* TODO rescale to remove denoms */
	L = RR_rescale(L,p);
  LL = RR_rescale(LL,p);
  L1 = RR_rescale(L1,p);
  L2 = RR_rescale(L2,p);
  Bad *= lcm(apply(S->lcm(apply(f->denominator(content(f)),S)),[L,L1,L2]));
  if(chi,
		print("T = part of J[",l,"] where Frob_",p," acts by ",chi);
		d = poldegree(chi); \\ Dimension of representation
		a = mordroot(chi,l) \\ q = p^a
	,
		print("T = all of J[",l,"]");
		d=2*g;
		a = mordroot(Lp,l)
	);
	J=PicInit(f,g,d0,L,LL,Bad,p,a,e);
	U=PicEvalInit(J,[L1,L2]); \\ Evaluation data
	J1 = PicRed(J,1); \\ Reduction mod p
	[B,matFrob] = TorsBasis(J1,l,Lp,chi); \\ Basis of the mod p^1 space and matrix of Frob_p
	print("The matrix of Frob is");
	printp(centerlift(matfrobenius(Mod(matFrob,l))));
	[WB,cWB] = TorsSpaceFrobGen(J1,l,B,matFrob); \\ Generating set of T under Frob and coordinates of these generators on B
	print("\n--> Lifting ",#WB," points ",p,"-adically");
	if(#WB > Jgetg(J),
  	my(J=J,l=l); WB = parapply(W->PicLiftTors(J,W,1,l),WB); \\ More efficient in parallel
	,
  	WB = apply(W->PicLiftTors(J,W,1,l),WB); \\ Less efficient in parallel (TODO tune)
	);
	print("\n--> All of T");
	TI = TorsSpaceFrob(J,WB,cWB,l,matFrob);
	print("\n--> Evaluation of ",#TI[2]," points");
	Z = TorsSpaceFrobEval(J,TI,U,l,d,matFrob);
	print("\n--> Expansion and identification");
	AF = TorsSpaceGetPols(J,Z); \\ List of polynomials defining the representation
	F = AF[1][3];
	ZF = apply(z->Mod(apply(c->c+O(p^e),z),JgetT(J)),AF[1][1]);
	[F,ZF];
}

HyperGalRep(f,l,p,e,P1,P2,chi)=
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
	GalRep(C,l,p,e,Lp,chi);
}

SmoothGalRep(f,l,p,e,P1,P2,chi)=
/* Computes the Galois representation afforded by
   the piece of l-torsion of the Jacobian
   of the hyperelliptic curve f(x,y)=0
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
  GalRep(C,l,p,e,Lp,chi);
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
