error!
PicNorm(X,F,WE)= \\ prod of F(P) for P in E, or 0 if one of F(P) is 0 or oo
\\ Assumes F in V
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,typ=Mod(T,p^e),
  WEV,S1,V1,r,d,v1,v2,S2,V2,K,K1,K2,M1,FS1,M2,m1,m2);
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

FqPrimroot(p,T)=Mod(ffprimroot(ffgen(Mod(T,p))).pol,T)*Mod(1,p);
dlog(x,a)=my(n=0,y=1);while(y!=x,n+=1;y*=a);n;
dlogmul(x,a,q,l)=Mod(dlog(x^((q-1)/l),a^((q-1)/l)),l);


PicFreyRuck(X,WT,l,WD,z)= \\ WT l-torsion, returns Frey-Ruck pairing
\\ If z!=0, then z must be a prim root of Fq; uses it for dlog. 
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,W00,C,H,WTm,i,j,s,q=p^poldegree(T));
 if(e>1,error("Not implemented."));
 if(q%l!=1,error("No l-th roots of 1"));
 W00=PicChord(X,W0,W0,1);
 C=AddChain(l);
 H=WTm=vector(#C);
 H[1]=1;
 WTm[1]=WT;
 for(c=2,#C,
  [i,j]=C[c][2..3];
  [WTm[c],s]=PicChord(X,WTm[i],if(j,WTm[j],W0),3);
  H[c]=PicNorm(X,s,W00)/(H[i]*if(j,H[j],1)*PicNorm(X,s,WD));
 );
 s=PicCertifZero(X,WTm[#C]);
 H[#C]*=PicNorm(X,s,WD)/PicNorm(X,s,W00);
 if(z,
  dlogmul(H[#C],z,q,l)
 ,
  H[#C]
 );
}

PicFreyRuckMulti(X,WT,l,W,z)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,W00,C,H,WTm,i,j,s,q=p^poldegree(T));
 if(e>1,error("Not implemented."));
 if(q%l!=1,error("No l-th roots of 1"));
 W00=PicChord(X,W0,W0,1);
 C=AddChain(l);
 WTm=vector(#C);
 H=matrix(#W,#C);
 for(d=1,#W,H[d,1]=1);
 WTm[1]=WT;
 for(c=2,#C,
  [i,j]=C[c][2..3];
  [WTm[c],s]=PicChord(X,WTm[i],if(j,WTm[j],W0),3);
  for(d=1,#W,H[d,c]=1/(H[d,i]*if(j,H[d,j],1)*PicNorm(X,s,W[d])));
  H[,c]*=PicNorm(X,s,W00);
 );
 s=PicCertifZero(X,WTm[#C]);
 for(d=1,#W,H[d,#C]*=PicNorm(X,s,W[d]));
 H[,#C]/=PicNorm(X,s,W00);
 if(z,
  apply(h->dlogmul(h,z,q,l),H[,#C])
 ,
  H[,#C]
 );
}

PicTorsRels(X,WT,l,excess=0)= \\ T vector of l-tors pts, g prim root of Fq
\\ Returns a matrix of relations containing all the relations satisfied by T (and probably no more)
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,q=p^poldegree(T),z,D,E);
 if(e>1,error("Not implemented."));
 if(q%l!=1,error("No l-th roots of 1"));
 z=FqPrimroot(p,T);
 D=vector(#T+excess,i,PicRand2(X)); \\ Find random elts in J
 \\ Pair them with elts of WT
 E=matconcat(vector(#WT,j,PicFreyRuckMulti(X,WT[j],l,D,z)));
 matker(E);
}
