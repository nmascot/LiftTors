PicNorm2(X,F,G,WE)= \\ prod of F/G(P) for P in E, or 0 in case F(P) or G(P) is 0 for some P
\\ Assumes F,G in V
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,typ=Mod(T,p^e),
  S1,V1,r,d,v1,v2,S2,V2,K,K1,K2,MF,FS1,MG,mF,mG);
 \\V2=DivAdd(V,V,6*d0+1-g,p,T,e); \\ V*V=H0(6D0)
 WEV=DivAdd(WE,V,5*d0+1-g,p,T,e); \\ V*W_E=H0(6D0-E)
 \\ O_E ~ V/W_E ~ V2/WEV
 \\ Want mat of Id and of mul by F
 \\ Find S s.t. V2 = S oplus W_E
 while(1,
  S1=matrix(#Z,d0);
  \\ Choose d0 elts of V at random
  for(j=1,d0,S1[,j]=RandVec(V,typ));
  V1=matconcat([S1,WE]);
  r=matrank(Mod(V1,p));
  d=3*d0+1-g;
  if(r==d,break,print1("@");print1([r,d]))
 );
 \\ Find S s.t. V2 = S oplus V*W_E
 while(1,
  S2=matrix(#Z,d0);
  \\ Choose d0 elts of V at random
  for(j=1,d0,
   v1=RandVec(V,typ);
   v2=RandVec(V,typ);
   for(P=1,#Z,
    S2[P,j]=v1[P]*v2[P]
   )
  );
  V2=matconcat([S2,WEV]);
  r=matrank(Mod(V2,p));
  d=6*d0+1-g;
  if(r==d,break,print1("@");print1([r,d]))
 );
 \\ Find coords of d0 first vectors of F*V1 on basis of V2 and project on d0 first vectors
 FS1=S1;
 for(P=1,#Z,FS1[P,]*=F[P]);
 K=matkerMod(matconcat([FS1,V2]),p,T,e);
 K1=K[1..d0,];
 K2=K[d0+1..2*d0,];
 K1=matinvMod(K1,p,T,e);
 MF=-K2*K1; \\ Matrix of mul F
 mF=matdetMod(MF,p,T,e);
 if(Mod(mF,p)==0,error("F has zeros on D");return(0));
 \\ Now same for mul by G
 FS1=S1;
 for(P=1,#Z,FS1[P,]*=G[P]);
 K=matkerMod(matconcat([FS1,V2]),p,T,e);
 K1=K[1..d0,];
 K2=K[d0+1..2*d0,];
 K1=matinvMod(K1,p,T,e);
 MG=-K2*K1; \\ Matrix of mul G
 mG=matdetMod(MG,p,T,e);
 if(Mod(mG,p)==0,error("G has zeros on D");return(0));
 mF/mG;
}

/*
PicCertifZero(X,W)= \\ W=W_D, assumes D~D0, finds s in V s.t. [s]=D-D0
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X);
 DivSub(W0,W,KV,1,p,T,e)[,1];
}

AddChain(n)=
{
 if(n==0,return([0,-1,-1]));
 if(n==1,return([[1,-1,-1]]));
 if(n==-1,return([[1,-1,-1],[-1,1,0]]));
 if(n%2,
  C=AddChain((n+1)\2);
  concat([C,[[-(n+1),#C,#C],[n,#C+1,1]]])
 ,
  C=AddChain(-n\2);
  concat([C,[[n,#C,#C]]])
 );
}*/

PicNorm2small(X,F,G,WE)= \\ prod of F/G(P) for P in E, or 0 in case F(P) or G(P) is 0 for some P
\\ Assumes F,G in V
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,typ=Mod(T,p^e),
  S1,V1,r,d,v1,v2,S2,V2,K,K1,K2,MF,FS1,MG,mF,mG);
 S1=RReval(Z,d0/2,df); \\ H0(D0)
 V20=RReval(Z,d0,df); \\ H0(2*D0)
 while(1,
  S2=matrix(#Z,d0);
  \\ Choose d0 elts of V at random
  for(j=1,d0,
   S2[,j]=RandVec(V20,typ)
  );
  V2=matconcat([S2,S1]);
  r=matrank(Mod(V2,p));
  d=2*d0+1-g;
  if(r==d,break,print1("@");print1([r,d]))
 );
 \\ Find coords of d0 first vectors of F*V1 on basis of V2 and project on d0 first vectors
 FS1=S1;
 for(P=1,#Z,FS1[P,]*=F[P]);
 K=matkerMod(matconcat([FS1,V2]),p,T,e);
 K1=K[1..d0,];
 K2=K[d0+1..2*d0,];
 K1=matinvMod(K1,p,T,e);
 MF=-K2*K1; \\ Matrix of mul F
 mF=matdetMod(MF,p,T,e);
 if(Mod(mF,p)==0,error("F has zeros on D");return(0));
 \\ Now same for mul by G
 FS1=S1;
 for(P=1,#Z,FS1[P,]*=G[P]);
 K=matkerMod(matconcat([FS1,V2]),p,T,e);
 K1=K[1..d0,];
 K2=K[d0+1..2*d0,];
 K1=matinvMod(K1,p,T,e);
 MG=-K2*K1; \\ Matrix of mul G
 mG=matdetMod(MG,p,T,e);
 if(Mod(mG,p)==0,error("G has zeros on D");return(0));
 mF/mG;
}

PicWeil2(X,WA,WB,l)= \\ WA, WB l-torsion, returns Weil pairing
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,C,HA,HB,H0,WA_,WB_,a,ha,b,hb,i,j,W,sA,sB);
 /*W00=PicChord(X,W0,W0);
 nu=PicCertifZero(X,W00); \\ [nu]=D'0-D0
 nu=(PicEvalNorm(X,nu,WA)/PicEvalNorm(X,nu,WB))^l;*/
 C=AddChain(l);
 HA=HB=H0=WA_=WB_=vector(#C);
 HA[1]=HB[1]=H0[1]=1;
 WA_[1]=WA;
 WB_[1]=WB;
 for(c=2,#C,
  [i,j]=C[c][2..3];
  [W,sA]=PicChord(X,WA_[i],if(j,WA_[j],W0),3);
  WA_[c]=W;
  HA[c]=1/(HA[i]*if(j,HA[j],1)*PicEvalNorm(X,sA,WB));
  [W,sB]=PicChord(X,WB_[i],if(j,WB_[j],W0),3);
  WB_[c]=W;
  HB[c]=1/(HB[i]*if(j,HB[j],1)*PicEvalNorm(X,sB,WA));
  H0[c]=1/(H0[i]*if(j,H0[j],1)*PicNorm2(X,sA,sB,W0))
 );
 a=PicCertifZero(X,WA_[#C]);
 b=PicCertifZero(X,WB_[#C]);
 (HA[#C]*PicEvalNorm(X,a,WB))/(HB[#C]*PicEvalNorm(X,b,WA))/(H0[#C]*PicNorm2small(X,a,b,W0));
}

dlog(x,a)=my(n=0,y=1);while(y!=x,n+=1;y*=a);n;
