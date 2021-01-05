RandTorsPt(J,l,a,Lp,Chi,Phi,seed)=
{
	/* Get a random nonzero pt in J[l] using randseed = seed
		 Actually returns [W,o,T,B(x)],
		 where W has order l^o,
		 T = l^(o-1)*W has order l,
		 B(x) monic multiple of Ann T.
		 Parameters:
		 Lp = charpoly(Frob|J)
		 a = we are working in J(F_p^a)
		 Chi: Actually get a point in J[l,Chi]
		 Phi: work in the submodule J[Phi] of J=J[x^a-1]
	*/
	my(A,B,W,N,v,M,chi,Psi,d,T,o);
	setrand(seed);
	A = B = 'x^a-1;
	W = PicRand(J,0);
	if(Phi,
		W = PicFrobPoly(J,W,B/Phi); \\ Project onto J[Phi]
		A = B = Phi
	);
  \\ Order of the submodule we work in
  N = polresultant(Lp,B);
  v = valuation(N,l); \\ N = l^v*M
  if(v==0,error("No l-torsion!"));
  M = N/l^v; \\ Cofactor, used to project on l-power part
	W = PicMul(J,W,M,0); \\ Project onto l-power part (main bottleneck!)
	B = gcd(Mod(Lp,l),Mod(B,l)); \\ now [l^...]W in J[l,B]
	if(Chi,
		chi = gcd(Chi,B);
		Psi = B/chi;
		d = poldegree(Psi);
		if(d,
			Psi /= polcoef(Psi,d);
			if(v>1,
				\\ lift cofactor l-adically
      	Psi = Mod(polhensellift(A,[liftint(A/Psi),liftint(Psi)],l,v)[2],l^v);
			);
			Psi = centerlift(Psi);
			W = PicFrobPoly(J,W,Psi) \\ Project onto J[chi]
		);
		B = chi;
		d = poldegree(B);
		B = centerlift(B/polcoef(B,d))
	);
	[T,o] = TorsOrd(J,W,l);
  if(o,return([W,o,T,B]));
  if(default(debug),print("RandTorsPt got zero (Phi=",Phi,")"));
  return(0);
}

Tors_UpdateLinTests(J,BT,Tnew,l,LinTests,R,FRparams)=
{
	/* BT contains d=#BT pts of J[l]
		 Tnew pt of J[l]
		 Assume that the matrix R = <BT[j]|LinTests[i]> has rank d. (H)
		 If Tnew is indep of BT, return [1,[LinTests,R']]
		 where R' of rank d+1
		 (may require modifying LinTests).
		 Else, return [0, relation between BT and Tnew].
	*/
	my(d,Rnew,BT2,R2,KR2,rel,NewTest,i);
	d = #BT;
	print("  Computing pairings...");
  \\ New col of R
  Rnew = TorsTestPt(J,Tnew,l,LinTests,FRparams);
  print("  Looking for relations...");
	BT2 = concat([BT,[Tnew]]);
	R2 = matconcat([R,Rnew]);
  KR2 = centerlift(matker(Mod(R2,l)));
	if(#KR2>1,error("Bug in Tors_UpdateLinTests, please report")); \\ Not supposed to happen by (H)
  if(#KR2==0,
    print("  Good, no relation");
		return([1,[LinTests,R2]]);
	);
	rel = KR2[,1];
	print(" Found pseudo-relation: ",rel);
	if(PicIsZero_val(J,PicLC(J,rel,BT2))==Jgete(J),
		print(" It actually holds.");
		return([0,rel])
	);
	print(" Good, it does not actually hold.");
  \\ So our linear forms are not independent.
  \\ Find a new one, and replace one the appropriate old one with it
  print("  Changing linear tests so that we don't get a false positive again.");
  until(Mod(Rnew*rel,l), \\ loop until new pairing breaks pseudo-relation
		NewTest = PicChord(J,PicRand(J,0),PicRand(J,0),1);
		Rnew = parapply(T->TorsTestPt(J,T,l,[NewTest],FRparams)[1],BT2);
		if(default(debug),print("  New test gives parings ",Rnew));
	);
  \\ So now we have r+1 forms of rank r.
  \\ Find one that can be removed.
  KR2 = matker(Mod(R2~,l))[,1]; \\ Find relation between forms
  i=1;
  while(KR2[i]==0,i++);
  if(default(debug)>=2,print("  Dropping form number ",i));
  \\ Replace the i-th old form with the new one
  LinTests[i] = NewTest;
	R2[i,] = Rnew;
  if(default(debug)>=2,
    print("  So now R=");
    printp(R2)
  );
	return([1,[LinTests,R2]]);
}

GuessColFromCharpoly(A,chi)=
{
	my(n=#A,B,chi0,M);
	for(i=1,n,A[i,n] = 0);
	A[n,n] = -(polcoef(chi,n-1)+trace(A));
	chi0 = chi-charpoly(A);
	B = matadjoint(x-A);
	M = matrix(n-1,n);
	for(i=1,n-1,
		for(k=1,n-1,
			M[k,i] = polcoef(B[n,i],k-1)
		)
	);
	for(k=1,n-1,M[k,n] = polcoef(chi0,k-1));
	M = matker(M);
	if(#M!=1,print("Unable to guess last column from charpoly");return(0));
	for(i=1,n-1,A[i,n] = M[i,1]/M[n,1]);
	A;
}

TorsGetMatAuts(J,KnownAuts,B,l,LinTests,R,FRparams)=
{
	my(nAuts,matAuts,Ri);
	nAuts = #JgetAutsCyc(J);
	if(KnownAuts==0,KnownAuts=vector(nAuts));
	matAuts = vector(nAuts);
	for(i=1,nAuts,
		if(KnownAuts[i],
			print(" Aut #",i," claimed to act as ",KnownAuts[i]);
			matAuts[i] = Mod(KnownAuts[i]*matid(#B),l);
		,
			Ri = parapply(W->TorsTestPt(J,PicAut(J,W,i),l,LinTests,FRparams),B);
			matAuts[i] = Mod(R,l)^(-1)*matconcat(Ri)
		)
	);
	matAuts;
}
	
TorsBasis(J,l,Lp,chi,KnownAuts,GetPairings)=
{
	/* Computes a basis B of the subspace T of J[l] on which Frob acts with charpoly chi
		 Assumes Lp = charpoly(Frob|J), so chi | Lp
		 If chi==0, then we take T=J[l]
		 Also computes the matrix M of Frob and list of matrices MA of Auts w.r.t B, and returns the vector [B,M,MA]
		 If GetPairings=1, returns [B,M,WT,R,X], where WT is a list of dim T points on J,
		 R is the matrix of pairings <B_j,WT_i>,
		 which is square and guaranted to be nonsingular,
		 and X are params needed for Frey-Rück.
	*/
	/* TODO use Auts */
	my(a,d,Phi,BW,Bo,BT,LinTests,R,matFrob,matAuts,AddC,W0,z,FRparams,r,iPhi,nBatch,UsedPhi,Batch);
	my(W,o,T,B,iFrob,res,rel,m,S,M);
	a = poldegree(JgetT(J)); \\ work over Fq=Fp[t]/T, q=p^a
  d = if(chi,poldegree(chi),poldegree(Lp)); \\ dim T
	if(Mod(a,l),
    Phi = vecsort(divisors(a),,4); \\ Divisors of a in reverse order
    Phi = apply(polcyclo,Phi); \\ Corresponding cyclotomic pols
    Phi = select(phi->poldegree(gcd(Mod(phi,l),Mod(if(chi,chi,Lp),l))),Phi) \\ Keep the ones compatible with charpoly chi
  );
	BW = vector(d); \\ list of l-power tors pts
  Bo = vector(d); \\ list of exponents of orders
  BT = vector(d); \\ list of l-tors pts
	LinTests = parvector(d,i,PicChord(J,PicRand(J,i+random()),PicRand(J,i+d+random()),1)); \\ list of pts to pair l-tors with
	R = matrix(d,0); \\ matrix of pairings
	matFrob = Mod(matrix(d,d),l);
	AddC = AddChain(l,0);
  W0 = JgetW0(J);
  W0 = PicChord(J,W0,W0,1); \\ Non-trivial origin, needed for pairings
  z = Fq_zeta_l(JgetT(J),Jgetp(J),l); \\ primitive l-th root of 1, to linearize parings
	FRparams = [AddC,W0,z];
  r = 0; \\ dim of l-tors obtained so far
  iPhi = 0; \\ Index of last used elt of Phi
  nBatch = 0; \\ Size of current batch
	while(1,
		\\ Make new batch
		nBatch = max(ceil((d-r)/a),ceil(default(nbthreads)/2));
		print(" Generating a new batch of ",nBatch," points in parallel");
    my(RandTorsPt=RandTorsPt,seed=vector(nBatch,i,random()));
    if(Phi,
      UsedPhi = vector(nBatch,i,Phi[(iPhi+i-1)%#Phi+1]); \\ TODO monitor dim of each Phi-part to choose Phi as needed
      iPhi += nBatch;
      print("   Using Phi=",UsedPhi)
    ,
      UsedPhi = vector(nBatch,i,0);
    );
    Batch = parvector(nBatch,i,RandTorsPt(J,l,a,Lp,chi,UsedPhi[i],seed[i]));
    print(" Batch of points generated.");
		\\ Loop through batch
		for(iBatch=1,nBatch,
			if(Batch[iBatch]==0,next); \\ Did we unluckily genrate the 0 pt?
			[W,o,T,B] = Batch[iBatch];
			print(" Picking point ",iBatch," from batch; its order is ",l,"^",o);
			while(1,
				iFrob = 0;
				while(1,
					r++;
					BW[r] = W; Bo[r] = o; BT[r] = T;
					res = Tors_UpdateLinTests(J,BT[1..r-1],T,l,LinTests,R,FRparams);
					if(res[1]==0, \\ Linear dependency
						rel = res[2];
						if(iFrob, \\ Have we added pts to BT in this loop? If yes, update matFrob.
							for(i=1,iFrob-1,matFrob[r-i,r-1-i]=Mod(1,l));
							for(i=1,r-1,matFrob[i,r-1]=Mod(-rel[i]/rel[r],l));
						);
						break
					);
					[LinTests,R] = res[2]; \\ Update Lintests and R
					iFrob++;
					if(B && iFrob==poldegree(B), \\ We know that next time will be dependent and the relation
						print(" B = ",centerlift(B), " says we'll get linear dependency next time");
						for(i=1,poldegree(B)-1,matFrob[r+1-i,r-i]=Mod(1,l));
            for(i=0,poldegree(B)-1,matFrob[r+1-poldegree(B)+i,r]=Mod(-polcoef(B,i),l));
						if(r==d,
							matAuts = TorsGetMatAuts(J,KnownAuts,BT,l,LinTests,R,FRparams);
							return(if(GetPairings,[BT,matFrob,matAuts,LinTests,R,FRparams],[BT,matFrob,matAuts]))
						);
						break(2);
					);
					if(r==d,
						for(i=1,iFrob-1,matFrob[r+1-i,r-i]=Mod(1,l));
						print("Guessing last column of matFrob from charpoly");
						M = GuessColFromCharpoly(matFrob,Mod(if(chi,chi,Lp),l));
						if(M,
							matFrob = M
						,
							M = TorsTestPt(J,PicFrob(J,BT[d]),l,LinTests,FRparams);
							rel = matker(Mod(matconcat([R,M]),l));
							if(#rel!=1,error("Bug in TorsGen, please report"));
							rel = rel[,1];
							for(i=1,d,matFrob[i,d] = -rel[i]/rel[d+1])
						);
						matAuts = TorsGetMatAuts(J,KnownAuts,BT,l,LinTests,R,FRparams);
						return(if(GetPairings,[BT,matFrob,matAuts,LinTests,R,FRparams],[BT,matFrob,matAuts]));
					);
					\\ Apply Frob and start over
					print(" Applying Frobenius");
					W = PicFrob(J,W);
					T = if(o==1,W,PicFrob(J,T));
				);
				\\ A relation broke the above loop. Try to use it to make a new point.
				m = vecmin([Bo[i]|i<-[1..r],rel[i]]);
        if(default(debug)>=2,
          print(" Bo=",Bo);
          print(" m=",m)
        );
				if(m<=1,
					\\ Cannot divide -> Give up
					print(" Giving up this point");
					r--; \\ Erase data about this point
					break;
				);
        S = vector(r,i,if(rel[i],l^(Bo[i]-m)*rel[i],0));
        W = PicLC(J,S,BW[1..r]);
				B = 0;
        r--; \\ Erase data about previous point, start over with this new one
        [T,o] = TorsOrd(J,W,l);
				print(" Dividing relation ",rel~," gives point of order ",l,"^",o);
				if(o==0,break)
			)
		)
	);
}
