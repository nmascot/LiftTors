\\GetCoef(A,a)=my(i,j,N=#A);[i,j]=liftall(Mod(a,N));if(i==0,i=N);if(j==0,j=N);A[i,j];
ZNnorm(x,N)=my(y=liftint(x)%N);if(y==0,N,y);
ZNneg(x,N)=my(y=lift(Mod(-x,N)));if(y==0,N,y);

ZNXspan(S,N)= \\ Span of S in (Z/NZ)*
{
	my(s,o,si,H,G);
	if(#S==0,return([Mod(1,N)]));
	s = Mod(S[1],N);
	o = znorder(s);
	H = ZNXspan(S[2..#S],N);
	G = List();
	si = Mod(1,N);
	for(i=1,o,
		si *= s;
		for(j=1,#H,
			listput(G,si*H[j])
		)
	);
	Vec(Set(G));
}

GetHlist(N,H)= \\ Span of H and -1, and reps of this gp / +-1
{ \\ Special codewords: H=0 -> (Z/NZ)*, H=1 -> +-1
	my(Hlist,Hlist1,Hdone);
	if(H==1,return([[1,-1],[1]]));
	Hlist = if(H==0,Mod(select(x->gcd(x,N)==1,vector(N,i,i)),N),ZNXspan(concat([H,[-1]]),N));
  Hlist1 = List();
  Hdone = vector(#Hlist);
  for(i=1,#Hlist,
    if(Hdone[i]==0,
      listput(Hlist1,Hlist[i]);
      Hdone[i] = 1;
      Hdone[select(h->Mod(h+Hlist[i],N)==0,Hlist,1)[1]] = 1;
    )
	);
	Hlist1 = Vec(Hlist1);
	[Hlist,Hlist1];
}

