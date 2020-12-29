TpClosure(J,l,BW,matTp,W,LinTests,R,FRparams)=
{
	my(dim,extradim,X,matTp_new);
	dim = #BW;
	extradim=0;
	while(1,
		X = Tors_UpdateLinTests(J,BW,W,l,LinTests,R,FRparams);
		if(X[1]==0, \\ W in span BW
			if(extradim,
				matTp_new = matrix(dim+extradim,dim+extradim,i,j,if(i<=dim && j<=dim,matTp[i,j]));
				for(i=1,extradim-1,matTp_new[dim+i+1,dim+i] = 1);
				for(i=1,dim+extradim,matTp_new[i,dim+extradim]=-X[2][i]/X[2][dim+extradim+1])
			,
				matTp_new = matTp
			);
			return([extradim,BW,matTp_new,LinTests,R]);
		);
		[LinTests,R]=X[2];
		extradim++;
		BW = concat(BW,[W]);
		W = PicTp(J,W);
	);
}	

TpEigen(J,l,Lp,chi,ap)=
{
	\\ TODO use l-power tors
	my(chiOO,dim,maxdim,a,BW,Bo,BT,LinTests,R,matTp,AddC,W0,z,FRparams,r,iPhi,nBatch);
	my(UsedPhi,Batch,W,o,T,B,iFrob,extradim,Eig,FrobR,matFrob);
	chiOO = Mod(chi,l);
  while(1,
    gcdchi = gcd(Lp/chiOO,chi);
    if(poldegree(gcdchi)==0,break);
    chiOO *= gcdchi
  );
	dim = poldegree(chi);
  maxdim = poldegree(chiOO); \\ upper bound on dim
	a = poldegree(JgetT(J)); \\ work over Fq=Fp[t]/T, q=p^a
  if(Mod(a,l),
    Phi = vecsort(divisors(a),,4); \\ Divisors of a in reverse order
    Phi = apply(polcyclo,Phi); \\ Corresponding cyclotomic pols
    Phi = select(phi->poldegree(gcd(Mod(phi,l),Mod(if(chi,chi,Lp),l))),Phi) \\ Keep the ones compatible with charpoly chi
  );
  \\BW = vector(maxdim); \\ list of l-power tors pts
  \\Bo = vector(maxdim); \\ list of exponents of orders
  BT = []; \\ list of l-tors pts
  LinTests = parvector(maxdim,i,PicChord(J,PicRand(J,i+random()),PicRand(J,i+maxdim+random()),1)); \\ list of pts to pair l-tors with
  R = matrix(0,0); \\ matrix of pairings
  matTp = Mod(matrix(0,0),l);
  AddC = AddChain(l,0);
  W0 = JgetW0(J);
  W0 = PicChord(J,W0,W0,1); \\ Non-trivial origin, needed for pairings
  z = Fq_zeta_l(JgetT(J),Jgetp(J),l); \\ primitive l-th root of 1, to linearize parings
  FRparams = [AddC,W0,z];
	r=0; \\ Dim generated so far
  iPhi = 0; \\ Index of last used elt of Phi
  nBatch = 0; \\ Size of current batch
  while(1,
    \\ Make new batch
    nBatch = max(maxdim-r,ceil(default(nbthreads)/2));
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
			W = Batch[iBatch];
      if(W==0,
				print(" Point ",iBatch, " of batch is 0.");
				next); \\ Did we unluckily generate the 0 pt?
      [W,o,T,B] = W;
      print(" Picking point ",iBatch," from batch; its order is ",l,"^",o);
      iFrob = 0;
			while(1,
				[extradim,BT,matTp,LinTests,R] = TpClosure(J,l,BT,matTp,T,LinTests,R,FRparams);
				print("Tp closure gives extra dim ",extradim);
				if(extradim==0,break);
				print("Now the matrix of Tp is");
				printp(matTp);
				Eig = matker(Mod(matTp-ap,l));
				print("Eigenspace of dim ",#Eig);
				if(#Eig==dim,
					BT = parvector(dim,j,PicLC(J,centerlift(Eig[,j]),BT));
					R = R*Eig;
					FrobR = matconcat(apply(W->Tors_TestPt(J,PicFrob(J,W),l,LinTests,FRparams),BT));
					matFrob = Mod(R,l)^(-1)*FrobR;
					return([BT,matFrob]);
				);
				iFrob++;
        if(B && iFrob==poldegree(B), \\ We know that next time will be dependent and the relation
          print(" B = ",centerlift(B), " says we'll get linear dependency next time");
					break;
				);
				print("  Applying Frobenius");
				T = PicFrob(J,T);
			)
		)
	);
}			
