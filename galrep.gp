read("MakTorsSpace.gp");

mordroot(f,p)=
\\ Computes the order of x in Fp[x]/(f). Fails if disc(f)=0 mod p.
{
 my(x=variable(f),N,fa,l,v);
 if(issquarefree(Mod(f,p))==0,error("Not squarefree!"));
 N=lcm(factormod(f,p,1)[,1]);
 N=p^N-1;
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

PicLC(J,C,W)=
\\ Computes sum_i C[i]*W[i] in J
\\ C vector of coefficients, W vector of points on J
{
  /* TODO efficiency */
  my(S);
  if(#C==0,return(JgetW0(J)));
  S = [];
  for(i=1,#C,
    if(C[i],
			if(S==[],
				S = PicMul(J,W[i],C[i],2);
			,
				S = PicAdd(J,S,PicMul(J,W[i],C[i],2));
			)
		)
  );
  if(S==[],JgetW0(J),S);
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
