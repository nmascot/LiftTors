PicRed(X,e1)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e2,Frob]=X);
 if(e1>e2,error("Cannot perform this reduction"));
 [f,df,Mod(V,p^e1),Mod(KV,p^e1),Mod(W0,p^e1),Mod(Z,p^e1),FrobCyc,g,d0,p,T,e1,Mod(Frob,p^e1)];
}

MakEqn(X,W)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,wV,KwV,K,Kj);
 wV=V;for(P=1,#Z,wV[P,]*=W[P,1]); \\ w*V, w = 1st col of W
 KwV=matkerMod(wV~,p,T,e)~;
 K=KV;
 for(j=2,#W,
  Kj=KwV;
  for(P=1,#Z,Kj[,P]*=W[P,j]);
  K=matconcat([K;Kj])
 );
 K;
}

MakEqnRand(X,W)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,w=0,wV,KwV,K,Kj);
 while(Mod(w,p)==0,w=RandVec(W,Mod(T,p^e))); \\ Get random elt from W
 wV=V;for(P=1,#Z,wV[P,]*=w[P]); \\ w*V
 KwV=matkerMod(wV~,p,T,e)~;
 K=KV;
 for(j=1,#W,
  Kj=KwV;
  for(P=1,#Z,Kj[,P]*=W[P,j]);
  K=matconcat([K;Kj])
 );
 K;
}

MakEqnRho(X,W,uv1=0)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,K,nZ,r,uv,A,B,C,D);
 nZ=#Z;
 neqns=(#W)*(nZ-#V);
 K=MakEqn(X,W);
 r=nZ-(d0+1-g);
 if(uv1==0,
  print("Warning: generating new uv");
  uv=FindMinorCompl(Mod(K,p))
 ,
  uv=uv1
 );
 [A,B,C,D]=M2ABCD(K,uv);
 Ainv=matinvMod(A,p,T,e);
 D-C*Ainv*B;
}

MakEqnRhoRand(X,W)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e,Frob]=X,K,nZ,r,uv,A,B,C,D);
 nZ=#Z;
 K=MakEqnRand(X,W);
 r=nZ-(d0+1-g);
 uv=FindMinorCompl(Mod(K,p));
 [A,B,C,D]=M2ABCD(K,uv);
 Ainv=matinvMod(A,p,T,e);
 D-C*Ainv*B;
}

PicDefEqn(j,d0,nW,neqnsV,nZ,KwVlist,VFlist,uv,AinvB,CAinv)=
{
 my(res=vector(d0),TK,HA,HB,HC,HD);
 for(k=1,d0,
  \\ We now deform by adding V0[k] to W[j]
  TK=vector(nW,j,matrix(neqnsV,nZ))~;
  if(j==1,
   for(i=2,nW,
    TK[i]=-KwVlist[k]*VFlist[i]
   )
  ,
   TK[j]=KwVlist[k]
  );
  [HA,HB,HC,HD]=M2ABCD(matconcat(TK),uv);
  res[k]=mat2col(HD-HC*AinvB+CAinv*(HA*AinvB-HB))
 );
 res;
}

PicLift_2(X2,W1,e1)=
{
 \\ TODO some steps can be done at lower prec
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,e2,Frob]=X2,nZ=#Z,r=nZ-(d0+1-g),neqnsV=nZ-#V,neqns=(#W1)*neqnsV);
 my(V1,W,wV,KwV,Kj,K);
 my(A,B,C,D,Ainv,CAinv,AinvB,rho,HA,HB,HC,HD,Mj); \\ TODO
 my(Trho,n,VF,VFlist,V0,M,KMq,H);
 my(Vmn,KwVmn,cW,sW);
 V1=Mod(V,p^e1);
 \\ Express W1 as a subspace of V1
 W=matrix(#V,#W1);
 for(j=1,#W1,
  K=matkerMod(matconcat([V1,W1[,j]]),p,T,e1);
  W[,j]=-K[1..#V,1]/K[1+#V,1]
 );
 W=V*Mod(liftint(W),p^e2); \\ Lift W as a supspace of V
 print("Calcul de rho");
 wV=V;for(P=1,nZ,wV[P,]*=W[P,1]); \\ w*V, w = 1st col of W
 KwV=matkerMod(wV~,p,T,e2)~;
 K=KV;
 for(j=2,#W,
  Kj=KwV;
  for(P=1,nZ,Kj[,P]*=W[P,j]);
  K=matconcat([K;Kj])
 );
 \\ scramble so that top left r*r block is invertible
 uv=FindMinorCompl(Mod(K,p));
 [A,B,C,D]=M2ABCD(K,uv);
 Ainv=matinvMod(A,p,T,e2);
 CAinv=C*Ainv;
 AinvB=Ainv*B;
 rho=D-CAinv*B;
 \\return([R,rho,W,S]);
 /* Compute the deformation of EqnW when an entry of W changes */
 print("Preparation pour Trho");
 Vmn=Mod(V,p^(e2-e1));
 VF=Vmn*matF(Mod(wV,p^(e2-e1)),p,T,e2-e1);
 VFlist=vector(#W,i,if(i==1,0,VF));
 for(i=2,#W,
  for(P=1,nZ,
   VFlist[i][,P]*=W[P,i] \\ TODO W1?
  )
 );
 /* Find which deformations stay in V and leave a minor of W fixed (rigidification)*/
 cW=matimagecompl(Mod(W1,p)~); \\ The rows we're going to change
 sW=VecsmallCompl(cW,nZ); \\ These rows of W form a minor
 V0=Vmn*matkerMod(matconcat(vector(#W,i,Vmn[sW[i],])~),p,T,e2-e1); \\ the d0-dim'l subspace of V whose sW rows are 0
 KwVmn=Mod(KwV,p^(e2-e1));
 KwVlist=vector(d0,k,matrix(neqnsV,nZ));
 for(k=1,d0,
  for(h=1,nZ-#W,
   Q=cW[h];
   KwVlist[k][,Q]=V0[Q,k]*KwVmn[,Q]
  )
 );
 Trho=matrix(nZ,#W,i,j,matrix(neqns-r,nZ-r));
 print("Calcul de Trho");
 \\TK=matrix(#W,nZ,Q,j,matrix(neqns,nZ));
 M=matrix((nZ-r)*(neqns-r),1+d0*(#W));
 M[,d0*(#W)+1]=mat2col(Mod(liftint(rho)/p^e1,p^(e2-e1)));
 n=0;
 my(PicDefEqn=PicDefEqn,nW=#W,d0=d0,neqnsV=neqnsV,nZ=nZ,KwVlist=KwVlist,VFlist=VFlist,uv=uv,AinvB=AinvB,CAinv=CAinv);
 Mj=parvector(#W,j,PicDefEqn(j,d0,nW,neqnsV,nZ,KwVlist,VFlist,uv,AinvB,CAinv));
 for(j=1,#W,
  \\print("w",j);
  for(k=1,d0,
   /*\\ We now deform by adding V0[k] to W[j]
   n+=1;
   TK=vector(#W,j,matrix(neqnsV,nZ))~;
   if(j==1,
    for(i=2,#W,
     TK[i]=-KwVlist[k]*VFlist[i]
    )
   ,
    TK[j]=KwVlist[k]
   );
   [HA,HB,HC,HD]=M2ABCD(matconcat(TK),uv);
   M[,n]=mat2col(HD-HC*AinvB+CAinv*(HA*AinvB-HB))*/
   n+=1;
   M[,n]=Mj[j][k];
  )
 );
 print("Noyau");
 KM=if(e1==1,matker(Mod(liftpol(M),T)),matkerMod(M,p,T,e2-e1));
 print("Dim noyau: ",#KM);
 \\ Find col whose last entry is nonzero
 for(k=1,#KM,
  if(Mod(KM[1+d0*(#W),k],p),
   K=KM[1..d0*(#W),k]/KM[1+d0*(#W),k];
   break
  )
 );
 n=0;
 H=matrix(nZ,#W);
 for(j=1,#W,
  for(k=1,d0,
   n+=1;
   H[,j]+=K[n]*V0[,k]
  )
 );
 W+Mod(p^e1*liftint(H),p^e2);
}

PicLiftTors(X,W,eini=1)=
{
 my([f,df,V,KV,W0,Z,FrobCyc,g,d0,p,T,eX,Frob]=X,e=eini,e2,W2=W,Xe);
 while(e<eX,
  e2=min(2*e,eX);
  print("Lift from O(",p,"^",e,") to O(",p,"^",e2,")");
  Xe=PicRed(X,e2);
  W2=PicLift_2(Xe,W2,e);
  print("Mul by ",p,"^",e2-e,"...");
  W2=PicMul(Xe,W2,p^(e2-e));
  e=e2
 );
 W2;
}
  
