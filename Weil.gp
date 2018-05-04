PicEvalNorm(X,F,WE)= \\ prod of F(P) for P in E, or 0 if one of f(P) is 0 or oo
\\ Assumes F in V
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,typ=Mod(T,p^e),
  WEV,S1,V1,r,d,v1,v2,S2,V2,K,K1,K2,M1,FS1,M2,m1,m2);
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
 \\ Find coords of d0 first vectors of V1 on basis of V2 and project on d0 first vectors
 K=matkerMod(matconcat([S1,V2]),p,T,e);
 K1=K[1..d0,];
 K2=K[d0+1..2*d0,];
 K1=matinvMod(K1,p,T,e);
 M1=-K2*K1; \\ Matrix of Id
 m1=matdetMod(M1,p,T,e);
 if(Mod(m1,p)==0,print("D intersects D0");return(0));
 \\ Now same for mul by F
 FS1=S1;
 for(P=1,#Z,FS1[P,]*=F[P]);
 K=matkerMod(matconcat([FS1,V2]),p,T,e);
 K1=K[1..d0,];
 K2=K[d0+1..2*d0,];
 K1=matinvMod(K1,p,T,e);
 M2=-K2*K1; \\ Matrix of mul F
 m2=matdetMod(M2,p,T,e);
 if(Mod(m2,p)==0,print("F has zeros on D");return(0));
 m2/m1;
}

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
}

PicWeil(X,WA,WB,l)= \\ WA, WB l-torsion, returns Weil pairing
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,W00,nu,C,HA,HB,WA_,WB_,i,j,W,s,a,b,ha,hb);
 W00=W0;
 while(1,
  W00=PicChord(X,W00,W0,1);
  nu=PicCertifZero(X,W00); \\ [nu]=D'0-D0
  a=PicEvalNorm(X,nu,WA);
  if(a==0,next);
  b=PicEvalNorm(X,nu,WB);
  if(b==0,next);
  nu=(a/b)^l;
  break
 );
 C=AddChain(l);
 HA=HB=WA_=WB_=vector(#C);
 HA[1]=HB[1]=1;
 WA_[1]=WA;
 WB_[1]=WB;
 for(c=2,#C,
  \\print(c);
  [i,j]=C[c][2..3];
  \\print("a");
  [W,s]=PicChord(X,WA_[i],if(j,WA_[j],W0),7);
  WA_[c]=W;
  HA[c]=PicEvalNorm(X,s,W00)/(HA[i]*if(j,HA[j],1)*PicEvalNorm(X,s,WB));
  \\print("b");
  [W,s]=PicChord(X,WB_[i],if(j,WB_[j],W0),7);
  WB_[c]=W;
  HB[c]=PicEvalNorm(X,s,W00)/(HB[i]*if(j,HB[j],1)*PicEvalNorm(X,s,WA));
 );
 a=PicCertifZero(X,WA_[#C]);
 ha=HA[#C]*PicEvalNorm(X,a,WB)/PicEvalNorm(X,a,W00);
 b=PicCertifZero(X,WB_[#C]);
 hb=HB[#C]*PicEvalNorm(X,b,WA)/PicEvalNorm(X,b,W00);
 ha/hb*nu;
}

dlog(x,a)=my(n=0,y=1);while(y!=x,n+=1;y*=a);n;
