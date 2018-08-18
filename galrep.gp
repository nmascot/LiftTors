read("MakTorsSpace.gp");

Z2pol(Z)=factorback(apply(u->'x-u,Z));

mordroot(f,p)=
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
{
  /* TODO efficiency */
  my(S);
  if(#C==0,return(JgetW0(J)));
  S = PicMul(J,W[1],C[1],2);
  for(i=2,#C,
    S = PicAdd(J,S,PicMul(J,W[i],C[i],2));
  );
  S;
}

TorsOrd(J,W,l)=
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

RandTorsPt(J,f,M,chiC)=
{
  my(W);
  while(1,
    W = PicRand0(J); \\ TODO
    W = PicMul(J,W,M,0);
    if(chiC,W = PicFrobPoly(J,W,chiC));
    if(!PicIsZero(J,W),return(W));
  );
}

PicTorsTrueRel(J,W,l)=
{
  my(R,ex=0);
	\\print("Getting rels between ",#W);
  while(1,
    R = PicTorsRels(J,W,l,ex);
		\\print("Found ",#R);
    if(#R==0,return(R));
    if(#R==1,
      if(PicIsZero(J,PicLC(J,R[,1],W)),return(R))
		/*,
		  return(Mat());*/
    );
		\\print("False positive");
    ex += 1
  );
}

TorsBasis(J,f,p,a,l,chi,C)=
{
  my(d,N,M,v,chiC,W,o,T,BW,Bo,BT,R,KR,Rnew,KRnew,Wtest,Wnew,am,S,AddC,W0,z);
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
	Wtest = vector(d,i,PicChord(J,PicRand0(J),PicRand0(J),1));
	z = Fq_zeta_l(JgetT(J),Jgetp(J),l);
  r = 0;
  while(r<d,
    print("Status:",Bo[1..r]);
    print("Getting new point");
    W = RandTorsPt(J,f,M,chiC);
    [T,o] = TorsOrd(J,W,l);
    print(" It has order l^",o);
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
				Wnew = PicChord(J,PicRand0(J),PicRand0(J),1);
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
      	S = vector(r,i,if(KR[i],l^(Bo[i]-m)*KR[i],0));
      	W = PicLC(J,S,BW[1..r]);
      	[T,o] = TorsOrd(J,W,l);
      	print(" gives point of order l^",o)
			,
        print(" Giving up this point");
        r -= 1;
        break
      );
      BW[r] = W;
      Bo[r] = o;
      BT[r] = T;
    );
  );
  BT;
}
