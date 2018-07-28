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
  my(d,N,M,v,chiC,W,o,T,BW,Bo,BT,R,m,S);
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
  r = 0;
  while(r<d,
    print("Status:",Bo[1..r]);
    print("Getting new point");
    W = RandTorsPt(J,f,M,chiC);
    [T,o] = TorsOrd(J,W,l);
    print("It has order l^",o);
    r += 1;
    BW[r] = W;
    Bo[r] = o;
    BT[r] = T;
    if(r==1,next);
    while(1,
      print("Looking for relations...");
      R = PicTorsTrueRel(J,BT[1..r],l);
      print("Found relation: ",R);
      if(#R == 0,print("Good, no relation");next(2));
      R = R[,1];
      m = vecmin([Bo[i]|i<-[1..r],R[i]]);
      print("Dividing relation ",R," by l^",m);
      S = vector(r,i,if(R[i],l^(Bo[i]-m)*R[i],0));
      W = PicLC(J,S,BW[1..r]);
      [T,o] = TorsOrd(J,W,l);
      print("gives point of order l^",o);
      if(o==0,
        print("Giving up this point.");
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
